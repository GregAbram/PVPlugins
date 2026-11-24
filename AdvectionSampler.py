"""
Advecting sampler -- offline Lagrangian point clouds for a concentration
field carried by a flow.

Given a time series with a velocity field and a concentration field (amount
per m^3, e.g. smoke in mg/m^3) on a fixed mesh, this produces a time series
of points in which every point stands for about the same amount of material,
q -- as MolarPointSampler.py does for one time step -- but the points move
with the flow and keep their identity from frame to frame.

    pvpython AdvectionSampler.py data/fire_plume/fire_plume.pvd \\
        --velocity velocity --concentration smoke --q 8 --diffusivity 2 \\
        --out out/plume

How it works
------------
Points start as a MolarPointSampler sample of the first step.  Then, for each
interval between saved steps:

  1. Advect.  Every point is moved with 4th-order Runge-Kutta through the
     velocity field, interpolated linearly in space (within cells) and in
     time (between the two saved steps).  Points that leave the domain are
     removed (outflow).  An optional random walk adds turbulent diffusion.
  2. Reconcile, at the end of the interval only (where the field is real
     data, not a blend of two steps).  The domain is divided into bins; in
     each bin the amount carried by the points is compared with the field's
     amount.  Small mismatches are absorbed by scaling the points' weights
     (by at most --max-growth up); larger deficits -- new smoke from a
     source, or smoke arriving through an inflow boundary -- produce new
     points; bins with too much lose weight.  Then a weight window keeps
     every weight near q: a point lighter than q/2 survives with probability
     w/q (and is set to q), a point heavier than 2q splits into pieces of
     about q, placed a little apart.
  3. Spread it over the interval.  Corrections are ramped across the
     interval's display frames instead of happening at the saved step:
     weights change geometrically, culled points fade out, split-off points
     fade in.  New points are "emitted" throughout the interval: each is
     placed where the field needs it at the end of the interval and traced
     backwards along the flow by a random fraction of the interval, so near
     a source it appears at the source and travels to its place.

Points carry their history: BirthTime / BirthPosition say when and where the
material was emitted (split-off points inherit their parent's), so Age shows
how old the smoke is and BirthPosition shows which part of a source produced
it.

Output (in --out)
-----------------
particles.pvd + particles_NNNNN.vtp   one point cloud per display frame
                                      (--frames-per-step per saved interval)
  ParticleId     persistent id (int64)
  Weight         amount this point stands for (about q)
  Fade           0..1 display weight: ramps in for new points and out for
                 culled ones; use it for opacity or size
  Age            time since the material was emitted (s)
  BirthTime      time of emission (s)
  BirthPosition  where it was emitted
  <conc>         concentration at the point, interpolated in space and time
  Velocity       flow velocity at the point
residual.pvd + residual_NNNNN.vti     per saved step, on the bin grid:
  FieldAmount    the field's amount in the bin
  ParticleAmount amount the advected points brought into the bin
  Residual       FieldAmount - ParticleAmount: what the flow does not
                 explain (sources, sinks, unresolved mixing)
  BirthAmount, CulledAmount   what reconciliation added / removed
diagnostics.csv                       one row per saved step (see --help)

Tuning
------
diagnostics.csv tells you how well the flow explains the data: "explained"
near 1 and low "churn" are good.  If the data spreads more than advection
alone accounts for (numerical or modelled diffusion), a modest --diffusivity
helps; on the synthetic plume K = 2 m^2/s raised explained from 0.74 to 0.77
and cut culls 3x, while K = 6 made it worse.  --max-growth trades splits
(pieces placed beside their parent) for births (emitted along paths); 1.5
puts births at about the fire's actual emission rate there.  Bins should
hold several points each (the first log line reports points per occupied
bin); use --bin-size to adjust.

Assumptions: fixed mesh; velocity in m/s; times in s; mesh coordinates times
--length-scale give meters.  Image data uses a fast numpy interpolator; any
other mesh uses vtkProbeFilter.  Multiblock input needs merging first.
Requires MolarPointSampler.py in the same directory.
"""
import argparse
import csv
import importlib.util
import math
import os
import sys
import time as _time

import numpy as np
from paraview import simple
from vtkmodules.util import numpy_support as ns
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (vtkCellArray, vtkDataObject, vtkImageData,
                                           vtkPolyData)
from vtkmodules.vtkFiltersCore import vtkCellDataToPointData
from vtkmodules.vtkFiltersCore import vtkProbeFilter
from vtkmodules.vtkIOXML import vtkXMLImageDataWriter, vtkXMLPolyDataWriter

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_sampler():
    """Import MolarPointSampler.py (same directory) for its sampling core."""
    spec = importlib.util.spec_from_file_location(
        "MolarPointSampler", os.path.join(HERE, "MolarPointSampler.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mps = _load_sampler()


# ==========================================================================
# Reading the time series
# ==========================================================================
class Series(object):
    """A ParaView-readable time series on a fixed mesh."""

    def __init__(self, path, velocity, concentration):
        self.reader = simple.OpenDataFile(path)
        if self.reader is None:
            sys.exit("cannot open %s" % path)
        self.times = list(self.reader.TimestepValues or [0.0])
        self.vname, self.cname = velocity, concentration

    def step(self, i):
        """(dataset, velocity at points (npts,3), concentration array, on_points)"""
        self.reader.UpdatePipeline(self.times[i])
        do = self.reader.GetClientSideObject().GetOutputDataObject(0)
        if not do.IsA("vtkDataSet"):
            sys.exit("input is %s; merge blocks first" % do.GetClassName())
        ds = do.NewInstance()
        ds.DeepCopy(do)
        pd, cd = ds.GetPointData(), ds.GetCellData()
        if pd.GetArray(self.vname) is not None:
            vel = ns.vtk_to_numpy(pd.GetArray(self.vname)).astype(np.float64)
        elif cd.GetArray(self.vname) is not None:
            c2p = vtkCellDataToPointData()
            c2p.SetInputData(ds)
            c2p.Update()
            vel = ns.vtk_to_numpy(c2p.GetOutput().GetPointData().GetArray(self.vname)).astype(np.float64)
        else:
            sys.exit("no velocity array '%s'" % self.vname)
        if vel.ndim != 2 or vel.shape[1] != 3:
            sys.exit("velocity '%s' must have 3 components" % self.vname)
        if pd.GetArray(self.cname) is not None:
            conc, on_points = pd.GetArray(self.cname), True
        elif cd.GetArray(self.cname) is not None:
            conc, on_points = cd.GetArray(self.cname), False
        else:
            sys.exit("no concentration array '%s'" % self.cname)
        conc = np.maximum(np.nan_to_num(ns.vtk_to_numpy(conc).astype(np.float64)), 0.0)
        return ds, vel, conc, on_points


# ==========================================================================
# Interpolation: locate points once, interpolate step-a and step-b values
# ==========================================================================
class ImageBackend(object):
    """Trilinear interpolation on vtkImageData, vectorized in numpy."""

    def __init__(self, img):
        self.origin = np.array(img.GetOrigin(), float)
        self.spacing = np.array(img.GetSpacing(), float)
        self.dims = np.array(img.GetDimensions(), int)
        if (self.dims < 2).any():
            raise ValueError("2D image data: use a 3D grid")
        self.lo = self.origin
        self.hi = self.origin + self.spacing * (self.dims - 1)

    def _locate(self, X):
        f = (X - self.origin) / self.spacing
        valid = ((X >= self.lo - 1e-9) & (X <= self.hi + 1e-9)).all(axis=1)
        i0 = np.clip(np.floor(f).astype(np.int64), 0, self.dims - 2)
        t = np.clip(f - i0, 0.0, 1.0)
        return i0, t, valid

    def point_values(self, X, arrays):
        """Interpolate point arrays (each (npts,) or (npts,c)) at X."""
        i0, t, valid = self._locate(X)
        nx, ny = self.dims[0], self.dims[1]
        out = [0.0 for _ in arrays]
        for dx in (0, 1):
            wx = t[:, 0] if dx else 1 - t[:, 0]
            for dy in (0, 1):
                wy = t[:, 1] if dy else 1 - t[:, 1]
                for dz in (0, 1):
                    wz = t[:, 2] if dz else 1 - t[:, 2]
                    idx = (i0[:, 0] + dx) + nx * ((i0[:, 1] + dy) + ny * (i0[:, 2] + dz))
                    wgt = wx * wy * wz
                    for n, a in enumerate(arrays):
                        v = a[idx]
                        out[n] = out[n] + (wgt[:, None] * v if v.ndim == 2 else wgt * v)
        return out, valid

    def cell_values(self, X, arrays):
        i0, _, valid = self._locate(X)
        nx, ny = self.dims[0] - 1, self.dims[1] - 1
        idx = i0[:, 0] + nx * (i0[:, 1] + ny * i0[:, 2])
        return [a[idx] for a in arrays], valid


class ProbeBackend(object):
    """Any mesh: vtkProbeFilter (point arrays interpolated, cell arrays by
    containing cell, vtkValidPointMask for inside/outside)."""

    def __init__(self, ds):
        self.mesh = ds.NewInstance()
        self.mesh.CopyStructure(ds)

    def _probe(self, X, arrays, attrs_getter):
        src = self.mesh.NewInstance()
        src.CopyStructure(self.mesh)
        for n, a in enumerate(arrays):
            va = ns.numpy_to_vtk(np.ascontiguousarray(a, dtype=np.float64), deep=1)
            va.SetName("a%d" % n)
            attrs_getter(src).AddArray(va)
        pts = vtkPoints()
        pts.SetData(ns.numpy_to_vtk(np.ascontiguousarray(X, dtype=np.float64), deep=1))
        pd = vtkPolyData()
        pd.SetPoints(pts)
        pr = vtkProbeFilter()
        pr.SetInputData(pd)
        pr.SetSourceData(src)
        pr.Update()
        o = pr.GetOutput().GetPointData()
        valid = ns.vtk_to_numpy(o.GetArray(pr.GetValidPointMaskArrayName())).astype(bool)
        return [ns.vtk_to_numpy(o.GetArray("a%d" % n)).astype(np.float64)
                for n in range(len(arrays))], valid

    def point_values(self, X, arrays):
        return self._probe(X, arrays, lambda d: d.GetPointData())

    def cell_values(self, X, arrays):
        return self._probe(X, arrays, lambda d: d.GetCellData())


class Flow(object):
    """Velocity and concentration over one interval, linear in time."""

    def __init__(self, backend, va, vb, ca, cb, conc_on_points, length_scale):
        self.b, self.va, self.vb, self.ca, self.cb = backend, va, vb, ca, cb
        self.cpts = conc_on_points
        self.L = length_scale

    def velocity(self, X, alpha):
        """Velocity in mesh units per second, and inside-domain mask."""
        (a, b), valid = self.b.point_values(X, [self.va, self.vb])
        return ((1 - alpha) * a + alpha * b) / self.L, valid

    def sample(self, X, alpha):
        """(velocity in m/s, concentration) at X for display."""
        (a, b), _ = self.b.point_values(X, [self.va, self.vb])
        get = self.b.point_values if self.cpts else self.b.cell_values
        (ca, cb), _ = get(X, [self.ca, self.cb])
        return (1 - alpha) * a + alpha * b, (1 - alpha) * ca + alpha * cb


def rk4(flow, X, t_a, dt_total, t, dt):
    """One RK4 step from time t by dt (negative = backwards).  Times are
    converted to the interval fraction alpha for the velocity blend.
    Returns new positions and a mask of points that stayed inside."""
    al = lambda tt: (tt - t_a) / dt_total
    k1, v1 = flow.velocity(X, al(t))
    k2, v2 = flow.velocity(X + 0.5 * dt * k1, al(t + 0.5 * dt))
    k3, v3 = flow.velocity(X + 0.5 * dt * k2, al(t + 0.5 * dt))
    k4, v4 = flow.velocity(X + dt * k3, al(t + dt))
    Xn = X + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    _, v5 = flow.velocity(Xn, al(t + dt))
    return Xn, v1 & v2 & v3 & v4 & v5


# ==========================================================================
# Field amounts on simplices and bins
# ==========================================================================
class Geometry(object):
    """Simplices of the (fixed) mesh, their volumes, and the bin grid."""

    def __init__(self, ds, length_scale, bin_size):
        sx = mps.simplicize(ds, mps._DIM_VOLUME)
        if sx.tg.GetNumberOfPoints() != ds.GetNumberOfPoints():
            raise ValueError("mesh splitting changed the point numbering")
        self.P, self.conn, self.source_cell = sx.P, sx.conn, sx.source_cell
        self.vol = mps.simplex_measure(self.P, self.conn, length_scale)
        b = np.array(ds.GetBounds(), float)
        self.lo, self.hi = b[0::2], b[1::2]
        ext = self.hi - self.lo
        self.h = (ext.prod() / max(ds.GetNumberOfCells(), 1)) ** (1 / 3.0)   # typical cell size
        self.bin = bin_size if bin_size else 4 * self.h
        self.bdims = np.maximum(np.ceil(ext / self.bin - 1e-9).astype(int), 1)
        self.nbins = int(self.bdims.prod())
        self.sbin = self.bin_of(self.P[self.conn].mean(axis=1))

    def bin_of(self, X):
        i = np.clip(np.floor((X - self.lo) / self.bin).astype(np.int64), 0, self.bdims - 1)
        return i[:, 0] + self.bdims[0] * (i[:, 1] + self.bdims[1] * i[:, 2])

    def amounts(self, conc, on_points):
        """Per-simplex amount and corner values (None for cell data)."""
        if on_points:
            cn = conc[self.conn]
            return self.vol * cn.mean(axis=1), cn
        return self.vol * conc[self.source_cell], None


# ==========================================================================
# Output
# ==========================================================================
def write_points(path, X, arrays):
    pd = vtkPolyData()
    pts = vtkPoints()
    pts.SetData(ns.numpy_to_vtk(np.ascontiguousarray(X, dtype=np.float32), deep=1))
    pd.SetPoints(pts)
    n = len(X)
    verts = vtkCellArray()
    verts.SetData(ns.numpy_to_vtkIdTypeArray(np.arange(n + 1, dtype=np.int64), deep=1),
                  ns.numpy_to_vtkIdTypeArray(np.arange(n, dtype=np.int64), deep=1))
    pd.SetVerts(verts)
    for name, a in arrays:
        va = ns.numpy_to_vtk(np.ascontiguousarray(a), deep=1)
        va.SetName(name)
        pd.GetPointData().AddArray(va)
    w = vtkXMLPolyDataWriter()
    w.SetFileName(path)
    w.SetInputData(pd)
    w.SetDataModeToAppended()
    w.Write()


def write_bins(path, geom, arrays):
    img = vtkImageData()
    img.SetDimensions(*(geom.bdims + 1))
    img.SetSpacing(geom.bin, geom.bin, geom.bin)
    img.SetOrigin(*geom.lo)
    for name, a in arrays:
        va = ns.numpy_to_vtk(np.ascontiguousarray(a, dtype=np.float64), deep=1)
        va.SetName(name)
        img.GetCellData().AddArray(va)
    w = vtkXMLImageDataWriter()
    w.SetFileName(path)
    w.SetInputData(img)
    w.SetDataModeToAppended()
    w.Write()


def write_pvd(path, entries):
    with open(path, "w") as fp:
        fp.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n<Collection>\n')
        for t, name in entries:
            fp.write('  <DataSet timestep="%.9g" file="%s"/>\n' % (t, name))
        fp.write("</Collection>\n</VTKFile>\n")


# ==========================================================================
# The particle system
# ==========================================================================
class Particles(object):
    FIELDS = ("X", "w", "pid", "tbirth", "pbirth", "appear")

    def __init__(self):
        self.X = np.zeros((0, 3))
        self.w = np.zeros(0)
        self.pid = np.zeros(0, np.int64)
        self.tbirth = np.zeros(0)
        self.pbirth = np.zeros((0, 3))
        self.appear = np.zeros(0)     # when it starts fading in (display only)
        self.next_id = 0

    def __len__(self):
        return len(self.w)

    def add(self, X, w, tbirth, pbirth, appear):
        n = len(X)
        self.X = np.vstack([self.X, X])
        self.w = np.concatenate([self.w, w])
        self.pid = np.concatenate([self.pid, np.arange(self.next_id, self.next_id + n)])
        self.next_id += n
        self.tbirth = np.concatenate([self.tbirth, tbirth])
        self.pbirth = np.vstack([self.pbirth, pbirth])
        self.appear = np.concatenate([self.appear, appear])

    def keep(self, mask):
        for f in self.FIELDS:
            setattr(self, f, getattr(self, f)[mask])


def seed_for(seed, i):
    return (seed * 1000003 + 7919 * (i + 1)) & 0x7FFFFFFF


def sample_amounts(geom, amount, cn, q, seed):
    """MolarPointSampler-style sample of per-simplex amounts -> positions."""
    n = mps.stochastic_counts(amount, q, seed)
    tet, k = mps.expand_counts(n)
    X, _ = mps.place_points(seed, geom.P, geom.conn, tet, k, None if cn is None else cn[tet])
    return X


class Run(object):
    def __init__(self, args):
        self.a = args
        self.series = Series(args.input, args.velocity, args.concentration)
        self.times = self.series.times
        last = len(self.times) - 1 if args.last is None else min(args.last, len(self.times) - 1)
        self.steps = list(range(args.first, last + 1))
        ds, vel, conc, cpts = self.series.step(self.steps[0])
        self.geom = Geometry(ds, args.length_scale, args.bin_size)
        self.backend = ImageBackend(ds) if (isinstance(ds, vtkImageData) and not args.probe) \
            else ProbeBackend(ds)
        self.cur = (vel, conc, cpts)
        self.p = Particles()
        self.rng = np.random.default_rng(args.seed)
        self.frames = []
        self.res_entries = []
        self.diag = []
        os.makedirs(args.out, exist_ok=True)

    # ------------------------------------------------------------ helpers
    def log(self, msg):
        if not self.a.quiet:
            print(msg, flush=True)

    def field_bins(self, amount):
        return np.bincount(self.geom.sbin, weights=amount, minlength=self.geom.nbins)

    def particle_bins(self, X, w):
        return np.bincount(self.geom.bin_of(X), weights=w, minlength=self.geom.nbins)

    def write_frame(self, t, X, w, fade, pid, tbirth, pbirth, flow, alpha):
        if flow is not None and len(X):
            vel, c = flow.sample(X, alpha)
        else:
            vel, c = np.zeros((len(X), 3)), np.zeros(len(X))
        name = "particles_%05d.vtp" % len(self.frames)
        f32 = lambda a: np.asarray(a, dtype=np.float32)      # display data: float32 halves the files
        write_points(os.path.join(self.a.out, name), X, [
            ("ParticleId", pid.astype(np.int64)), ("Weight", f32(w)), ("Fade", f32(fade)),
            ("Age", f32(t - tbirth)), ("BirthTime", f32(tbirth)), ("BirthPosition", f32(pbirth)),
            (self.a.concentration, f32(c)), ("Velocity", f32(vel))])
        self.frames.append((t, name))

    # ------------------------------------------------------------ initial
    def initialize(self):
        vel, conc, cpts = self.cur
        t0 = self.times[self.steps[0]]
        amount, cn = self.geom.amounts(conc, cpts)
        X = sample_amounts(self.geom, amount, cn, self.a.q, seed_for(self.a.seed, -1))
        w = np.full(len(X), self.a.q)
        # Match every bin exactly so later residuals measure transport only.
        M, W = self.field_bins(amount), self.particle_bins(X, w)
        f = np.divide(M, W, out=np.ones_like(M), where=W > 0)
        w = w * f[self.geom.bin_of(X)]
        self.p.add(X, w, np.full(len(X), t0), X.copy(), np.full(len(X), -np.inf))
        flow = Flow(self.backend, vel, vel, conc, conc, cpts, self.a.length_scale)
        self.write_frame(t0, self.p.X, self.p.w, np.ones(len(X)), self.p.pid,
                         self.p.tbirth, self.p.pbirth, flow, 0.0)
        self.log("t=%g: %d points, total %.6g (field %.6g), %.1f points per occupied bin"
                 % (t0, len(X), w.sum(), amount.sum(),
                    len(X) / max(1, (self.particle_bins(X, w) > 0).sum())))

    # ------------------------------------------------------------ interval
    def interval(self, i_a, i_b):
        a = self.a
        t_a, t_b = self.times[i_a], self.times[i_b]
        T = t_b - t_a
        va, ca, cpts = self.cur
        ds, vb, cb, cpts_b = self.series.step(i_b)
        if cpts_b != cpts:
            sys.exit("concentration association changes between steps")
        flow = Flow(self.backend, va, vb, ca, cb, cpts, a.length_scale)
        F = a.frames_per_step
        vmax = max(np.linalg.norm(va, axis=1).max(), np.linalg.norm(vb, axis=1).max(), 1e-12)
        nsub = int(math.ceil(T * vmax / a.length_scale / (a.cfl * self.geom.h) / F)) * F
        nsub = max(nsub, F)
        dt = T / nsub
        sub_per_frame = nsub // F
        frame_t = t_a + T * np.arange(1, F + 1) / F
        rng = np.random.default_rng([a.seed, i_b])

        # ---- 1. advect existing points ------------------------------------
        p = self.p
        X = p.X.copy()
        alive = np.ones(len(p), bool)
        frames_X = []
        t = t_a
        for s in range(nsub):
            Xn, inside = rk4(flow, X[alive], t_a, T, t, dt)
            if a.diffusivity > 0:
                Xn += math.sqrt(2 * a.diffusivity * dt) / a.length_scale * rng.standard_normal(Xn.shape)
            if a.ground_z is not None:
                Xn[:, 2] = np.maximum(Xn[:, 2], a.ground_z + 1e-6 * self.geom.h)
            idx = np.nonzero(alive)[0]
            X[idx] = Xn
            alive[idx[~inside]] = False
            t = t_a + (s + 1) * dt
            if (s + 1) % sub_per_frame == 0:
                frames_X.append((X.copy(), alive.copy()))
        outflow = p.w[~alive].sum()
        Xb, wb = X[alive], p.w[alive]

        # ---- 2. reconcile at t_b --------------------------------------------
        amount_b, cn_b = self.geom.amounts(cb, cpts)
        M = self.field_bins(amount_b)
        W = self.particle_bins(Xb, wb)
        resid = M - W
        f = np.minimum(np.divide(M, W, out=np.zeros_like(M), where=W > 0), a.max_growth)
        deficit = np.maximum(M - f * W, 0.0)
        bins_b = self.geom.bin_of(Xb)
        w_new = wb * f[bins_b]

        # births: field amount in deficit bins, sampled from the field there
        share = np.divide(deficit, M, out=np.zeros_like(M), where=M > 0)
        birth_amount = amount_b * share[self.geom.sbin]
        Xn_b = sample_amounts(self.geom, birth_amount, cn_b, a.q, seed_for(a.seed, i_b))

        # weight window on the survivors
        light = w_new < 0.5 * a.q
        survive = rng.random(len(w_new)) < w_new / a.q
        culled = light & ~survive
        w_new[light & survive] = a.q
        culled_amount_bins = np.bincount(bins_b[culled], weights=wb[culled] * f[bins_b[culled]],
                                         minlength=self.geom.nbins)
        heavy = w_new > 2 * a.q
        nsplit = np.where(heavy, np.round(w_new / a.q).astype(np.int64), 1)

        # ---- 3. emit births along their back-traced paths -------------------
        nb = len(Xn_b)
        tau = rng.random(nb) * T                      # how long before t_b it was emitted
        tbirth_b = t_b - tau
        pbirth_b = Xn_b.copy()
        bframes = [None] * F
        bframes[F - 1] = Xn_b.copy()
        Y = Xn_b.copy()
        bin_ok = np.ones(nb, bool)                    # still inside going backwards
        t = t_b
        for s in range(nsub):
            go = bin_ok & (tbirth_b < t - 1e-12)      # not yet back at its birth time
            if not go.any():
                break
            Yn, inside = rk4(flow, Y[go], t_a, T, t, -dt)
            idx = np.nonzero(go)[0]
            # leaving backwards = it entered through the boundary at time t
            tbirth_b[idx[~inside]] = np.maximum(tbirth_b[idx[~inside]], t)
            pbirth_b[idx[~inside]] = Y[idx[~inside]]
            bin_ok[idx[~inside]] = False
            ok = idx[inside]
            Y[ok] = Yn[inside]
            t = t_b - (s + 1) * dt
            reached = ok[tbirth_b[ok] >= t - 1e-12]
            pbirth_b[reached] = Y[reached]
            j = F - 1 - (s + 1) // sub_per_frame      # frame index at time t
            if (s + 1) % sub_per_frame == 0 and j >= 0:
                bframes[j] = Y.copy()
        tbirth_b = np.maximum(tbirth_b, t_a)

        # ---- 4. frames for this interval ------------------------------------
        # Display weight ramps geometrically from w0 to its end value: the
        # scaled weight for culled points (while Fade goes to 0), otherwise
        # the corrected weight, divided among the pieces if it will split
        # (the split-off pieces fade in during the next interval).
        w0 = p.w
        w_end = np.zeros(len(p))
        w_end[alive] = np.where(culled, wb * f[bins_b], w_new / nsplit)
        ratio = np.divide(w_end, w0, out=np.ones_like(w0), where=w0 > 0)
        ratio = np.where(ratio > 0, ratio, 1.0)
        culled_full = np.zeros(len(p), bool)
        culled_full[np.nonzero(alive)[0][culled]] = True
        fade_time = a.fade if a.fade is not None else 0.5 * T
        for j in range(F):
            tj, alpha = frame_t[j], (j + 1) / F
            Xj, aj = frames_X[j]
            show = aj
            wj = w0 * ratio ** alpha
            fade = np.ones(len(p))
            fade = np.minimum(fade, np.clip((tj - p.appear) / fade_time, 0, 1))
            fade[culled_full] = 1 - alpha
            parts = [(Xj[show], wj[show], fade[show], p.pid[show], p.tbirth[show], p.pbirth[show])]
            # births visible once emitted
            bshow = tbirth_b <= tj + 1e-12
            if bframes[j] is not None and bshow.any():
                bfade = np.clip((tj - tbirth_b[bshow]) / fade_time, 0, 1)
                parts.append((bframes[j][bshow], np.full(bshow.sum(), a.q), bfade,
                              p.next_id + np.nonzero(bshow)[0], tbirth_b[bshow], pbirth_b[bshow]))
            cat = [np.concatenate([q[k] for q in parts]) for k in range(6)]
            self.write_frame(tj, *cat, flow=flow, alpha=alpha)

        # ---- 5. new state at t_b ---------------------------------------------
        keep = np.zeros(len(p), bool)
        keep[np.nonzero(alive)[0][~culled]] = True
        new_w = w_new[~culled]
        new_n = nsplit[~culled]
        p.keep(keep)
        p.X = Xb[~culled]
        p.w = new_w / new_n
        nsurv = len(p)
        # births first: the frames above already showed them with these ids
        p.add(Xn_b, np.full(nb, a.q), tbirth_b, pbirth_b, tbirth_b)
        # children of splits: same history, placed a little apart, fade in
        par = np.repeat(np.arange(nsurv), new_n - 1)
        if len(par):
            off = rng.standard_normal((len(par), 3))
            off *= (a.split_jitter or 0.5 * self.geom.h) / np.linalg.norm(off, axis=1, keepdims=True)
            Xc = p.X[par] + off * rng.random((len(par), 1)) ** (1 / 3.0)
            _, ok = self.backend.point_values(Xc, [vb[:, 0]])
            Xc[~ok] = p.X[par][~ok]
            p.add(Xc, p.w[par], p.tbirth[par], p.pbirth[par], np.full(len(par), t_b))
        if len(p) > a.max_particles:
            sys.exit("%d particles exceeds --max-particles; raise --q" % len(p))
        self.cur = (vb, cb, cpts)

        # ---- diagnostics ------------------------------------------------------
        name = "residual_%05d.vti" % i_b
        write_bins(os.path.join(a.out, name), self.geom, [
            ("FieldAmount", M), ("ParticleAmount", W), ("Residual", resid),
            ("BirthAmount", birth_amount_bins(self.geom, birth_amount)),
            ("CulledAmount", culled_amount_bins)])
        self.res_entries.append((t_b, name))
        row = dict(time=t_b, points=len(p), field_total=M.sum(), carried=W.sum(),
                   particle_total=p.w.sum(), outflow=outflow,
                   explained=1 - np.abs(resid).sum() / max(M.sum(), 1e-300),
                   births=nb, birth_amount=birth_amount.sum(), culls=int(culled.sum()),
                   splits=int((new_n - 1).sum()), substeps=nsub,
                   churn=(nb + int(culled.sum())) / max(len(p), 1))
        self.diag.append(row)
        self.log("t=%g: %d pts | explained %.3f | births %d culls %d splits %d | "
                 "outflow %.4g | particles %.6g vs field %.6g"
                 % (t_b, len(p), row["explained"], nb, row["culls"], row["splits"],
                    outflow, row["particle_total"], row["field_total"]))

    def run(self):
        t0 = _time.time()
        self.initialize()
        for i_a, i_b in zip(self.steps[:-1], self.steps[1:]):
            self.interval(i_a, i_b)
        write_pvd(os.path.join(self.a.out, "particles.pvd"), self.frames)
        write_pvd(os.path.join(self.a.out, "residual.pvd"), self.res_entries)
        with open(os.path.join(self.a.out, "diagnostics.csv"), "w", newline="") as fp:
            if self.diag:
                wr = csv.DictWriter(fp, fieldnames=list(self.diag[0].keys()))
                wr.writeheader()
                wr.writerows(self.diag)
        self.log("wrote %d frames to %s in %.1f s" % (len(self.frames), self.a.out, _time.time() - t0))


def birth_amount_bins(geom, birth_amount):
    return np.bincount(geom.sbin, weights=birth_amount, minlength=geom.nbins)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Lagrangian point clouds of a concentration field carried by a flow.",
        epilog="diagnostics.csv columns: time, points, field_total, carried (amount the "
               "advected points brought), particle_total (after reconciliation), outflow, "
               "explained (1 - sum|field - carried| / field, per bin: 1 = the flow explains "
               "everything), births, birth_amount, culls, splits, substeps, churn "
               "((births + culls) / points).")
    ap.add_argument("input", help="time series readable by ParaView (e.g. .pvd)")
    ap.add_argument("--velocity", required=True, help="velocity array (m/s), point or cell")
    ap.add_argument("--concentration", required=True, help="concentration array, point or cell")
    ap.add_argument("--q", type=float, required=True, help="amount per point")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--frames-per-step", type=int, default=4,
                    help="display frames per saved interval (default 4)")
    ap.add_argument("--bin-size", type=float, default=None,
                    help="reconciliation bin size in mesh units (default 4 cells)")
    ap.add_argument("--max-growth", type=float, default=1.5,
                    help="largest weight increase per interval before new points are "
                         "made instead (default 1.5; lower = fewer splits, more births)")
    ap.add_argument("--cfl", type=float, default=0.5,
                    help="max cells a point moves per RK4 sub-step (default 0.5)")
    ap.add_argument("--diffusivity", type=float, default=0.0,
                    help="random-walk diffusivity K in m^2/s (default 0: off)")
    ap.add_argument("--ground-z", type=float, default=None,
                    help="keep points above this z (mesh units)")
    ap.add_argument("--fade", type=float, default=None,
                    help="fade-in time for new points, s (default half an interval)")
    ap.add_argument("--split-jitter", type=float, default=None,
                    help="how far split-off points are placed, mesh units (default half a cell)")
    ap.add_argument("--length-scale", type=float, default=1.0, help="meters per mesh unit")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--first", type=int, default=0, help="first saved step")
    ap.add_argument("--last", type=int, default=None, help="last saved step")
    ap.add_argument("--max-particles", type=int, default=5000000)
    ap.add_argument("--probe", action="store_true",
                    help="use vtkProbeFilter even for image data (for testing)")
    ap.add_argument("--quiet", action="store_true")
    return ap.parse_args(argv)


if __name__ == "__main__":
    Run(parse_args()).run()
