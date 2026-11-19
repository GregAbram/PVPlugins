"""
Tests for MolarPointSampler.  Run with ParaView's Python:

    /Applications/ParaView-6.1.1.app/Contents/bin/pvpython tests/test_molar_point_sampler.py

Prints one line per check and exits non-zero if any fail.  The expected
error messages printed by the guard tests ("exceeds Max Points", "no 2D
(surface) cells", the "no concentration array selected" warning) are normal,
as are the "vtkExecutive ... returned failure for request" ERR lines that
ParaView prints for those deliberate failures.
"""
import os
import sys

import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from vtkmodules.util import numpy_support as ns
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import VTK_WEDGE, vtkUnstructuredGrid

HERE = os.path.dirname(os.path.abspath(__file__))
LoadPlugin(os.path.join(HERE, "..", "MolarPointSampler.py"), ns=globals())

failures = []


def check(name, ok, detail=""):
    print("%s  %s  %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        failures.append(name)


def fetch(f):
    o = sm.Fetch(f)
    P = ns.vtk_to_numpy(o.GetPoints().GetData()) if o.GetNumberOfPoints() else np.zeros((0, 3))
    fd = o.GetFieldData()
    fdd = {fd.GetArrayName(i): fd.GetArray(i).GetValue(0) for i in range(fd.GetNumberOfArrays())}
    return o, P, fdd


def arr(o, name):
    return ns.vtk_to_numpy(o.GetPointData().GetArray(name))


# ---------------------------------------------------------------- wedge mesh
def tapered_wedges(n=6, taper=0.3):
    """Unit-square triangles extruded to z=1, top shrunk toward the center."""
    xs = np.linspace(0, 1, n + 1)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    bot = np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)]
    top = bot.copy()
    top[:, :2] = 0.5 + (top[:, :2] - 0.5) * (1 - taper)
    top[:, 2] = 1
    P = np.vstack([bot, top])
    off = len(bot)
    idx = lambda i, j: i * (n + 1) + j
    ug = vtkUnstructuredGrid()
    pts = vtkPoints()
    pts.SetData(ns.numpy_to_vtk(P, deep=1))
    ug.SetPoints(pts)
    for i in range(n):
        for j in range(n):
            for t in ((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)),
                      (idx(i, j), idx(i + 1, j + 1), idx(i, j + 1))):
                ug.InsertNextCell(VTK_WEDGE, 6, list(t) + [v + off for v in t])
    return ug, P


def set_cell_array(ug, name, values):
    a = ns.numpy_to_vtk(np.asarray(values, dtype=np.float64), deep=1)
    a.SetName(name)
    ug.GetCellData().AddArray(a)
    ug.Modified()


TAPER = 0.3
ug, Pw = tapered_wedges(taper=TAPER)
nc = ug.GetNumberOfCells()
set_cell_array(ug, "cconst", np.full(nc, 1000.0))
px = ns.numpy_to_vtk(1000.0 * Pw[:, 0], deep=1)
px.SetName("px")
ug.GetPointData().AddArray(px)
src = TrivialProducer()
src.GetClientSideObject().SetOutput(ug)

zz = (np.arange(200000) + 0.5) / 200000
side2 = (1 - TAPER * zz) ** 2
V = side2.mean()                                    # frustum volume

f = MolarPointSampler(Input=src, Concentration=["CELLS", "cconst"], AmountPerPoint=0.01)
o, P, fd = fetch(f)
check("volume: total amount", np.isclose(fd["TotalAmount_mmol"], 1000 * V), "%.3f" % fd["TotalAmount_mmol"])
expected = 1000 * V / 0.01
check("volume: point count", abs(len(P) - expected) < 0.01 * expected, "%d vs %.0f" % (len(P), expected))
lower = side2[zz < 0.5].sum() / side2.sum()
check("volume: uniform in tapered wedges", abs((P[:, 2] < 0.5).mean() - lower) < 0.01,
      "lower-half fraction %.4f vs %.4f" % ((P[:, 2] < 0.5).mean(), lower))
check("volume: SampledDimension", fd["SampledDimension"] == 3)

# Point data, c = 1000 x: sample mean of x must equal the c-weighted mean.
f.Concentration = ["POINTS", "px"]
f.AmountPerPoint = 0.005
o, P, fd = fetch(f)
xs = (np.arange(400) + 0.5) / 400
num = den = 0.0
for zv in (np.arange(400) + 0.5) / 400:
    s = 1 - TAPER * zv
    x = 0.5 + (xs - 0.5) * s
    num += (x * x).sum() * s * s
    den += x.sum() * s * s
check("volume point data: total amount", np.isclose(fd["TotalAmount_mmol"], 500 * V, rtol=1e-6))
check("volume point data: density follows linear field", abs(P[:, 0].mean() - num / den) < 0.005,
      "mean x %.4f vs %.4f" % (P[:, 0].mean(), num / den))
check("volume point data: sampled value interpolated",
      np.abs(arr(o, "px") - 1000 * P[:, 0]).max() < 1e-9)

# Temporal stability: change one cell, others keep identical points.
f.Concentration = ["CELLS", "cconst"]
f.AmountPerPoint = 0.01
o1, P1, _ = fetch(f)
a = np.full(nc, 1000.0)
a[5] = 3000.0
set_cell_array(ug, "cconst", a)
src.MarkModified(src)
o2, P2, _ = fetch(f)
s1, s2 = arr(o1, "SourceCellId"), arr(o2, "SourceCellId")
check("stability: unchanged cells identical", np.array_equal(P1[s1 != 5], P2[s2 != 5]))
old = {tuple(p) for p in P1[s1 == 5]}
new = {tuple(p) for p in P2[s2 == 5]}
check("stability: changed cell keeps its points", old <= new, "%d -> %d" % (len(old), len(new)))

f.Seed = 7
_, P3, _ = fetch(f)
check("seed: different seed, different points", not np.array_equal(P2, P3))
f.Seed = 12345

f.MaxPoints = 100
f.UpdatePipeline()
check("max points: output cleared", f.GetDataInformation().GetNumberOfPoints() == 0)

# ------------------------------------------------------------- wavelet cases
w = Wavelet()                                         # deliberately not updated first
fw = MolarPointSampler(Input=w, AmountPerPoint=2550)
fw.UpdatePipeline()
n = fw.GetDataInformation().GetNumberOfPoints()
check("default array fallback (RTData)", 450 < n < 550, "%d points" % n)

c = Calculator(Input=w, ResultArrayName="C", Function="(RTData/277)^8")
v = Calculator(Input=c, ResultArrayName="Pos", Function="coords")
i = Calculator(Input=v, ResultArrayName="Label", Function="floor(coordsX)", ResultArrayType="Int")
cd = PointDatatoCellData(Input=i, PassPointData=1)
cd.UpdatePipeline()
fp = MolarPointSampler(Input=cd, Concentration=["POINTS", "C"], AmountPerPoint=0.4)
o, P, fd = fetch(fp)
names = {o.GetPointData().GetArrayName(k) for k in range(o.GetPointData().GetNumberOfArrays())}
check("pass arrays: names", {"C", "RTData", "Pos", "Label", "RTData_cell", "Pos_cell",
                             "Label_cell", "SourceCellId"} <= names, str(sorted(names)))
check("pass arrays: linear vector interpolated exactly", np.abs(arr(o, "Pos") - P).max() < 1e-9)
check("pass arrays: integer stays integer", np.issubdtype(arr(o, "Label").dtype, np.integer))
srcd = sm.Fetch(cd)
sid = arr(o, "SourceCellId")
check("pass arrays: cell arrays from source cell",
      np.array_equal(arr(o, "RTData_cell"),
                     ns.vtk_to_numpy(srcd.GetCellData().GetArray("RTData"))[sid]))
fp.PassArrays = 0
o, _, _ = fetch(fp)
check("pass arrays off", o.GetPointData().GetNumberOfArrays() == 2)

# ------------------------------------------------------------------ surfaces
sph = Sphere(Radius=1, ThetaResolution=64, PhiResolution=64)
c1 = Calculator(Input=sph, AttributeType="Cell Data", ResultArrayName="C", Function="100")
c1.UpdatePipeline()
fs = MolarPointSampler(Input=c1, Concentration=["CELLS", "C"], AmountPerPoint=0.1)
o, P, fd = fetch(fs)
area = sm.Fetch(IntegrateVariables(Input=sph)).GetCellData().GetArray("Area").GetValue(0)
expected = 100 * area / 0.1
r = np.linalg.norm(P, axis=1)
check("sphere: point count", abs(len(P) - expected) < 0.02 * expected, "%d vs %.0f" % (len(P), expected))
check("sphere: on the surface", r.min() > 0.99 and r.max() <= 1 + 1e-9)
check("sphere: even split", abs((P[:, 2] > 0).mean() - 0.5) < 0.02)
check("sphere: SampledDimension", fd["SampledDimension"] == 2)

pl = Plane(XResolution=10, YResolution=10)            # quads on [-0.5, 0.5]^2
c2 = Calculator(Input=pl, ResultArrayName="C", Function="coordsX + 0.5")
c2.UpdatePipeline()
fq = MolarPointSampler(Input=c2, Concentration=["POINTS", "C"], AmountPerPoint=0.0001)
o, P, fd = fetch(fq)
check("plane: point count", abs(len(P) - 5000) < 150, "%d" % len(P))
check("plane: density follows linear field", abs(P[:, 0].mean() - (2 / 3 - 0.5)) < 0.01,
      "mean x %.4f" % P[:, 0].mean())
fq.LengthScale = 2
o, P2, _ = fetch(fq)
check("plane: length scale squared", abs(len(P2) / len(P) - 4) < 0.1, "%d -> %d" % (len(P), len(P2)))

wv = Wavelet()
wv.UpdatePipeline()
mix = AppendDatasets(Input=[wv, ExtractSurface(Input=wv)])
mix.UpdatePipeline()
fm = MolarPointSampler(Input=mix, Concentration=["POINTS", "RTData"], AmountPerPoint=2550)
o, P, fd = fetch(fm)
check("mixed: auto picks volume", fd["SampledDimension"] == 3)
fm.CellDimension = "Surface"
o, P, fd = fetch(fm)
check("mixed: surface mode on outer faces",
      fd["SampledDimension"] == 2 and len(P) > 0 and np.isclose(np.abs(P), 10).any(axis=1).all())
fv = MolarPointSampler(Input=wv, Concentration=["POINTS", "RTData"], AmountPerPoint=2550,
                       CellDimension="Surface")
fv.UpdatePipeline()
check("surface mode on pure volume: empty", fv.GetDataInformation().GetNumberOfPoints() == 0)

print("\n%d failure(s)" % len(failures) + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
