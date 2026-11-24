"""
Synthetic smoke plume from a fire in a breeze, for testing AdvectionSampler.py.

Solves the advection-diffusion equation for a smoke concentration c (mg/m^3)

    dc/dt + div(c u) = S(x, t) + K laplacian(c)

on a regular grid with an analytic, time-varying wind u(x, t):
  * a breeze along +x with a power-law profile (zero at the ground) that
    gusts in strength and slowly veers (a y component that swings),
  * a buoyant updraft above the fire whose strength follows the fire's
    output, tilted downwind and widening with height.
The fire is a pulsing source S near the ground.  Because c is computed from
u with a conservative scheme, the two are consistent the way real simulation
output is -- which is what the particle reconciliation needs to be tested on.

Numerics: node-centered finite volumes, first-order upwind fluxes (positive,
conservative; its numerical diffusion stands in for unresolved turbulence),
explicit diffusion, CFL-limited sub-steps.  Inflow boundaries bring in clean
air; outflow is zero-gradient; no flux through the ground.

Output: <out>/fire_plume.pvd and one .vti per saved step, with point arrays
"velocity" (m/s) and "smoke" (mg/m^3).  Coordinates in meters, times in
seconds from the first saved step (after a spin-up).

    pvpython synthetic/make_fire_plume.py --out data/fire_plume
"""
import argparse
import os

import numpy as np
from vtkmodules.util import numpy_support as ns
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.vtkIOXML import vtkXMLImageDataWriter

FIRE = np.array([60.0, 80.0])       # fire center (x, y), m
Q0 = 2000.0                          # mean smoke output, mg/s


def fire_output(t):
    """Pulsing smoke output Q(t), mg/s."""
    return Q0 * np.maximum(0.1, 1 + 0.6 * np.sin(2 * np.pi * t / 25)
                           + 0.2 * np.sin(2 * np.pi * t / 7))


def wind(t, X, Y, Z):
    """Breeze + fire updraft at time t on node coordinates X, Y, Z."""
    profile = (np.maximum(Z, 0.0) / 10.0) ** 0.2
    u = 3.0 * (1 + 0.3 * np.sin(2 * np.pi * t / 40)) * profile
    v = 0.8 * np.sin(2 * np.pi * t / 55) * profile
    xc = FIRE[0] + 0.6 * Z                       # plume tilts downwind
    sig = 10.0 + 0.15 * Z                        # and widens with height
    r2 = (X - xc) ** 2 + (Y - FIRE[1]) ** 2
    w_peak = 6.0 * (fire_output(t) / Q0) ** (1 / 3)
    w = w_peak * np.exp(-r2 / (2 * sig ** 2)) * (1 - np.exp(-Z / 3.0)) * np.exp(-Z / 120.0)
    return u, v, w


def source_shape(X, Y, Z, dV):
    """Normalized emission distribution near the ground (sums to 1 over dV)."""
    g = np.exp(-((X - FIRE[0]) ** 2 + (Y - FIRE[1]) ** 2) / (2 * 6.0 ** 2)
               - (Z - 2.0) ** 2 / (2 * 3.0 ** 2))
    g[Z > 10] = 0
    return g / (g.sum() * dV)


def flux_divergence(c, vel, axis, h):
    """div of the upwind flux c*vel along one axis, with inflow = clean air
    and outflow = zero-gradient at the two domain faces."""
    cm = np.moveaxis(c, axis, 0)
    vm = np.moveaxis(vel, axis, 0)
    vf = 0.5 * (vm[1:] + vm[:-1])                         # interior faces
    f_int = np.where(vf > 0, vf * cm[:-1], vf * cm[1:])
    f_lo = np.where(vm[0] < 0, vm[0] * cm[0], 0.0)[None]  # outflow only
    f_hi = np.where(vm[-1] > 0, vm[-1] * cm[-1], 0.0)[None]
    f = np.concatenate([f_lo, f_int, f_hi], axis=0)
    return np.moveaxis((f[1:] - f[:-1]) / h, 0, axis)


def laplacian(c, h):
    """Zero-flux (Neumann) Laplacian."""
    p = np.pad(c, 1, mode="edge")
    return ((p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1] + p[1:-1, :-2, 1:-1]
             + p[1:-1, 1:-1, 2:] + p[1:-1, 1:-1, :-2] - 6 * c) / h ** 2)


def write_vti(path, shape, h, u, v, w, c):
    img = vtkImageData()
    img.SetDimensions(*shape)
    img.SetSpacing(h, h, h)
    img.SetOrigin(0, 0, 0)
    vel = np.stack([a.ravel(order="F") for a in (u, v, w)], axis=1).astype(np.float32)
    va = ns.numpy_to_vtk(vel, deep=1)
    va.SetName("velocity")
    ca = ns.numpy_to_vtk(c.ravel(order="F").astype(np.float32), deep=1)
    ca.SetName("smoke")
    img.GetPointData().AddArray(va)
    img.GetPointData().AddArray(ca)
    img.GetPointData().SetActiveScalars("smoke")
    wr = vtkXMLImageDataWriter()
    wr.SetFileName(path)
    wr.SetInputData(img)
    wr.SetDataModeToAppended()
    wr.SetCompressorTypeToZLib()
    wr.Write()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default="data/fire_plume")
    ap.add_argument("--extent", type=float, nargs=3, default=[300.0, 160.0, 80.0],
                    help="domain size x y z (m)")
    ap.add_argument("--spacing", type=float, default=4.0, help="grid spacing (m)")
    ap.add_argument("--steps", type=int, default=41, help="saved steps")
    ap.add_argument("--interval", type=float, default=5.0, help="seconds between saved steps")
    ap.add_argument("--spinup", type=float, default=60.0, help="seconds simulated before saving")
    ap.add_argument("--diffusivity", type=float, default=0.5, help="K (m^2/s)")
    args = ap.parse_args(argv)

    h = args.spacing
    shape = tuple(int(round(e / h)) + 1 for e in args.extent)
    X, Y, Z = np.meshgrid(*[np.arange(n) * h for n in shape], indexing="ij")
    G = source_shape(X, Y, Z, h ** 3)
    c = np.zeros(shape)
    K = args.diffusivity
    os.makedirs(args.out, exist_ok=True)

    t = -args.spinup
    entries = []
    for step in range(args.steps):
        t_save = step * args.interval
        while t < t_save - 1e-9:
            u, v, w = wind(t, X, Y, Z)
            vmax = max(np.abs(u).max(), np.abs(v).max(), np.abs(w).max(), 1e-6)
            dt = min(0.4 * h / vmax, 0.15 * h * h / max(K, 1e-12), t_save - t)
            div = (flux_divergence(c, u, 0, h) + flux_divergence(c, v, 1, h)
                   + flux_divergence(c, w, 2, h))
            c = c + dt * (-div + fire_output(t) * G + K * laplacian(c, h))
            c = np.maximum(c, 0.0)
            t += dt
        u, v, w = wind(t_save, X, Y, Z)
        name = "fire_plume_%03d.vti" % step
        write_vti(os.path.join(args.out, name), shape, h, u, v, w, c)
        entries.append((t_save, name))
        print("t=%6.1f s  total smoke %.4g mg  max %.4g mg/m^3"
              % (t_save, c.sum() * h ** 3, c.max()))

    with open(os.path.join(args.out, "fire_plume.pvd"), "w") as fp:
        fp.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n<Collection>\n')
        for tv, name in entries:
            fp.write('  <DataSet timestep="%g" file="%s"/>\n' % (tv, name))
        fp.write("</Collection>\n</VTKFile>\n")


if __name__ == "__main__":
    main()
