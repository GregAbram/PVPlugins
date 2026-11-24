"""
Tests for AdvectionSampler.py.  Run with ParaView's Python:

    /Applications/ParaView-6.1.1.app/Contents/bin/pvpython tests/test_advection_sampler.py

1. Exact transport: a Gaussian blob carried by solid-body rotation plus a
   slow rise.  The velocity is linear in space and constant in time, so
   trilinear interpolation is exact, and the concentration at every step is
   known analytically.  Points must stay on their circles, rise at the right
   rate, and need almost no correction.  The probe backend (any mesh) must
   give the same paths as the image backend.
2. Fire plume: synthetic/make_fire_plume.py (short run).  Totals must track
   the field, the flow must explain most of the change, new material must
   come from the fire, and ids must persist from frame to frame.

Writes temporary data under the system temp directory and removes it.
Prints PASS/FAIL per check and exits non-zero on failure.
"""
import csv
import glob
import os
import shutil
import sys
import tempfile

import numpy as np
from vtkmodules.util import numpy_support as ns
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.vtkIOXML import vtkXMLImageDataWriter, vtkXMLPolyDataReader

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, ROOT)
import AdvectionSampler  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print("%s  %s  %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)
    if not ok:
        failures.append(name)


def read_frame(path):
    r = vtkXMLPolyDataReader()
    r.SetFileName(path)
    r.Update()
    o = r.GetOutput()
    d = {o.GetPointData().GetArrayName(i): ns.vtk_to_numpy(o.GetPointData().GetArray(i))
         for i in range(o.GetPointData().GetNumberOfArrays())}
    d["X"] = ns.vtk_to_numpy(o.GetPoints().GetData()).astype(np.float64)
    return d


def diagnostics(out):
    rows = list(csv.DictReader(open(os.path.join(out, "diagnostics.csv"))))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def run(*args):
    AdvectionSampler.Run(AdvectionSampler.parse_args(list(args) + ["--quiet"])).run()


tmp = tempfile.mkdtemp(prefix="advect_test_")
try:
    # ------------------------------------------------------------ 1. rotation
    OMEGA, W0, CEN = 2 * np.pi / 60, 0.05, np.array([20.0, 20.0])
    BLOB, SIG = np.array([28.0, 20.0, 6.0]), 3.0
    dims = (41, 41, 21)
    X, Y, Z = np.meshgrid(*[np.arange(n, dtype=float) for n in dims], indexing="ij")
    rot = os.path.join(tmp, "rot")
    os.makedirs(rot)
    entries = []
    for i in range(11):
        t = 2.0 * i
        ang = OMEGA * t                                  # blob center at time t
        c0 = np.array([CEN[0] + (BLOB[0] - CEN[0]) * np.cos(ang) - (BLOB[1] - CEN[1]) * np.sin(ang),
                       CEN[1] + (BLOB[0] - CEN[0]) * np.sin(ang) + (BLOB[1] - CEN[1]) * np.cos(ang),
                       BLOB[2] + W0 * t])
        conc = 100 * np.exp(-((X - c0[0]) ** 2 + (Y - c0[1]) ** 2 + (Z - c0[2]) ** 2) / (2 * SIG ** 2))
        vel = np.stack([-OMEGA * (Y - CEN[1]), OMEGA * (X - CEN[0]), np.full(X.shape, W0)], axis=-1)
        img = vtkImageData()
        img.SetDimensions(*dims)
        for name, a in (("velocity", np.stack([vel[..., k].ravel(order="F") for k in range(3)], axis=1)),
                        ("c", conc.ravel(order="F"))):
            va = ns.numpy_to_vtk(np.ascontiguousarray(a), deep=1)
            va.SetName(name)
            img.GetPointData().AddArray(va)
        w = vtkXMLImageDataWriter()
        w.SetFileName(os.path.join(rot, "rot_%02d.vti" % i))
        w.SetInputData(img)
        w.Write()
        entries.append((t, "rot_%02d.vti" % i))
    AdvectionSampler.write_pvd(os.path.join(rot, "rot.pvd"), entries)

    out = os.path.join(tmp, "rot_out")
    run(os.path.join(rot, "rot.pvd"), "--velocity", "velocity", "--concentration", "c",
        "--q", "0.5", "--out", out)
    dg = diagnostics(out)
    frames = sorted(glob.glob(os.path.join(out, "particles_*.vtp")))
    check("rotation: frame count", len(frames) == 1 + 10 * 4, "%d" % len(frames))
    check("rotation: totals track the field",
          np.abs(dg["particle_total"] / dg["field_total"] - 1).max() < 0.01,
          "max rel err %.4f" % np.abs(dg["particle_total"] / dg["field_total"] - 1).max())
    check("rotation: flow explains the change", dg["explained"].min() > 0.9,
          "min explained %.3f" % dg["explained"].min())
    check("rotation: little correction needed", dg["churn"].max() < 0.03,
          "max churn %.4f, births %d, culls %d" % (dg["churn"].max(), dg["births"].sum(), dg["culls"].sum()))
    f0, f1 = read_frame(frames[0]), read_frame(frames[-1])
    common, i0, i1 = np.intersect1d(f0["ParticleId"], f1["ParticleId"], return_indices=True)
    r0 = np.hypot(f0["X"][i0, 0] - CEN[0], f0["X"][i0, 1] - CEN[1])
    r1 = np.hypot(f1["X"][i1, 0] - CEN[0], f1["X"][i1, 1] - CEN[1])
    a0 = np.arctan2(f0["X"][i0, 1] - CEN[1], f0["X"][i0, 0] - CEN[0])
    a1 = np.arctan2(f1["X"][i1, 1] - CEN[1], f1["X"][i1, 0] - CEN[0])
    dang = np.angle(np.exp(1j * (a1 - a0)))
    check("rotation: most points survive", len(common) > 0.95 * len(i0), "%d of %d" % (len(common), len(f0["X"])))
    check("rotation: points stay on their circles", np.abs(r1 - r0).max() < 1e-3,
          "max radius change %.2e" % np.abs(r1 - r0).max())
    check("rotation: turned by omega*t", np.abs(dang - OMEGA * 20).max() < 1e-3,
          "angle %.5f vs %.5f" % (np.median(dang), OMEGA * 20))
    check("rotation: rose by w*t", np.abs(f1["X"][i1, 2] - f0["X"][i0, 2] - W0 * 20).max() < 1e-6)
    check("rotation: concentration interpolated",
          np.abs(f1["c"] - 100 * np.exp(-np.sum((f1["X"] - [CEN[0] + (BLOB[0] - CEN[0]) * np.cos(OMEGA * 20),
                                                              CEN[1] + (BLOB[0] - CEN[0]) * np.sin(OMEGA * 20),
                                                              BLOB[2] + W0 * 20]) ** 2, axis=1)
                                         / (2 * SIG ** 2))).max() < 5.0,
          "trilinear vs exact Gaussian, within 5% of the peak (100)")

    out_p = os.path.join(tmp, "rot_probe")
    rp = AdvectionSampler.Run(AdvectionSampler.parse_args(
        [os.path.join(rot, "rot.pvd"), "--velocity", "velocity", "--concentration", "c",
         "--q", "0.5", "--out", out_p, "--last", "2", "--probe", "--quiet"]))
    check("probe backend selected", isinstance(rp.backend, AdvectionSampler.ProbeBackend))
    rp.run()
    fi = read_frame(sorted(glob.glob(os.path.join(out, "particles_*.vtp")))[8])
    fp = read_frame(sorted(glob.glob(os.path.join(out_p, "particles_*.vtp")))[8])
    common, ii, ip = np.intersect1d(fi["ParticleId"], fp["ParticleId"], return_indices=True)
    check("probe backend matches image backend",
          len(common) > 0.95 * len(fi["X"]) and np.abs(fi["X"][ii] - fp["X"][ip]).max() < 1e-4,
          "%d common, max diff %.2e" % (len(common), np.abs(fi["X"][ii] - fp["X"][ip]).max()))

    # ---------------------------------------------------------- 2. fire plume
    plume = os.path.join(tmp, "plume")
    sys.path.insert(0, os.path.join(ROOT, "synthetic"))
    import make_fire_plume
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        make_fire_plume.main(["--out", plume, "--steps", "9"])
    out = os.path.join(tmp, "plume_out")
    run(os.path.join(plume, "fire_plume.pvd"), "--velocity", "velocity",
        "--concentration", "smoke", "--q", "8", "--out", out, "--diffusivity", "2")
    dg = diagnostics(out)
    frames = sorted(glob.glob(os.path.join(out, "particles_*.vtp")))
    check("plume: frame and residual counts",
          len(frames) == 1 + 8 * 4 and len(glob.glob(os.path.join(out, "residual_*.vti"))) == 8)
    check("plume: totals track the field",
          np.abs(dg["particle_total"] / dg["field_total"] - 1).max() < 0.02,
          "max rel err %.4f" % np.abs(dg["particle_total"] / dg["field_total"] - 1).max())
    check("plume: flow explains most of the change", dg["explained"].mean() > 0.65,
          "mean explained %.3f" % dg["explained"].mean())
    check("plume: births about the fire's emission",
          0.5 < dg["birth_amount"].mean() / (2000 * 5) < 2.0,
          "%.0f mg per interval vs 10000 emitted" % dg["birth_amount"].mean())
    last = read_frame(frames[-1])
    new = last["BirthTime"] > 0
    bp = last["BirthPosition"][new]
    near = (np.hypot(bp[:, 0] - 60, bp[:, 1] - 80) < 25) & (bp[:, 2] < 20)
    check("plume: new material comes from the fire", near.mean() > 0.6,
          "%.0f%% born within 25 m of the fire" % (100 * near.mean()))
    check("plume: ids unique", len(np.unique(last["ParticleId"])) == len(last["ParticleId"]))
    check("plume: fade in [0,1]", ((last["Fade"] >= 0) & (last["Fade"] <= 1)).all())
    prev = read_frame(frames[-2])
    common, a, b = np.intersect1d(prev["ParticleId"], last["ParticleId"], return_indices=True)
    move = np.linalg.norm(last["X"][b] - prev["X"][a], axis=1)
    vmax = np.linalg.norm(last["Velocity"], axis=1).max()
    check("plume: ids persist, moves bounded by the wind",
          len(common) > 0.8 * len(prev["X"]) and np.median(move) < 1.25 * vmax,
          "%d shared, median move %.2f m (wind up to %.1f m/s, 1.25 s/frame)"
          % (len(common), np.median(move), vmax))
    check("plume: age = time - birth", np.allclose(last["Age"], 40.0 - last["BirthTime"], atol=1e-4))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("\n%d failure(s)" % len(failures) + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
