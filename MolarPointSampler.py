"""
Molar Point Sampler -- ParaView Python plugin.

Turns a concentration field on a volume mesh (amount per m^3) or a surface
mesh (amount per m^2) into a point cloud in which every point stands for the
same amount of material, q.  Point density is therefore proportional to
concentration: a region with twice the concentration gets twice the points
per unit volume (or area), and a cell twice as large at the same
concentration gets twice the points.  Two materials sampled with the same q
get point counts in proportion to their total moles, so they can be compared
directly.

Tested with ParaView 5.13.3 and 6.1.1 (pvpython).

Loading
-------
GUI: Tools > Manage Plugins > Load New..., pick this file (tick Auto Load to
keep it), then Filters > Sampling > Molar Point Sampler.  After editing this
file, press "Reload Python Module" in the filter's Properties panel.

pvpython:
    from paraview.simple import *
    LoadPlugin("/path/to/MolarPointSampler.py", ns=globals())
    src = OpenDataFile("data.vtu"); src.UpdatePipeline()
    s = MolarPointSampler(Input=src)
    s.Concentration = ["CELLS", "conc"]       # or ["POINTS", ...]
    s.AmountPerPoint = 0.5                     # mmol per point

Properties
----------
Concentration     Point or cell array, single component.  Units are amount
                  per cubic meter (e.g. mmol/m^3) for volumes, per square
                  meter (mmol/m^2) for surfaces; q is then in the same
                  amount unit (mmol).  If nothing is selected -- e.g. a
                  script created the filter before its input had data -- the
                  input's active point scalars are used, else active cell
                  scalars, else the first single-component array, with a
                  warning in Output Messages.
Amount Per Point  q, the amount each point represents.  Expected points =
                  total amount / q.  Must be > 0.
Cell Dimension    Auto (default): 3D cells if the input has any, otherwise
                  2D cells; if it has both (e.g. a volume mesh that also
                  stores boundary faces) the 2D cells are ignored with a
                  warning.  Volume: 3D cells only.  Surface: 2D cells only
                  (triangles, quads, polygons, strips).  1D and 0D cells are
                  never sampled.
Point Data        Only for a point-data Concentration.  Linear (default):
Placement         points follow the interpolated field inside each cell.
                  Uniform Per Cell: each cell's amount is still computed
                  exactly from the point data (same totals), but its points
                  are spread uniformly in the cell, as for cell data.  Points
                  then stay still between time steps -- see Method -- at the
                  cost of detail inside cells, which only shows as density
                  steps at cell faces when cells hold many points each.
                  Refining the mesh (e.g. Tessellate) hides the steps and
                  keeps the stability.
Length Scale      Meters per mesh unit (1 for meters, 0.001 for mm), so that
                  volumes are in m^3 (areas in m^2).  Point counts scale with
                  its cube for volumes, its square for surfaces.
Seed              Changes which random points are drawn, not how many on
                  average.  Same seed + same input = identical output.
Max Points        Safety limit (default 5,000,000).  If the count would
                  exceed it the filter stops, clears its output and reports
                  the smallest q that fits.
Pass Arrays       Carry the input's other arrays onto the points (default
                  on); see Output.

Choosing q
----------
For a target of N points: q = TotalAmount_mmol / N.  Run once with any q,
read TotalAmount_mmol from the output's field data (Information tab), divide.
Counts scale as 1/q.  When comparing materials or animating a time series,
choose q once (e.g. from the largest material or the peak time step) and keep
it fixed; recomputing it per material or per frame hides real differences.

Output (vtkPolyData of vertices)
--------------------------------
Point data:
  <concentration>   the concentration at each point: the cell value, or the
                    linear interpolation of the point values (after clamping,
                    see below) -- interpolated even with Uniform Per Cell
                    placement.
  other arrays      (Pass Arrays on) point arrays interpolated at the point
                    with the same barycentric weights that placed it --
                    floating-point arrays of any component count are
                    interpolated, integer arrays (labels, ids) take the
                    nearest corner's value.  Cell arrays are copied from the
                    containing cell.  If a point and a cell array share a
                    name, the second gets a "_cell" (or "_point") suffix.
                    String/bit arrays and vtkGhostType / vtkValidPointMask
                    are not passed.
  SourceCellId      index of the input cell the point lies in.
Field data:
  TotalAmount_mmol      total amount in the input (after clamping).
  AmountPerPoint_mmol   q as used.
  ClampedNegative_mmol  amount removed by treating negative/NaN values as 0.
  SampledDimension      3 if volume cells were sampled, 2 if surface cells.
  (The _mmol suffix assumes mmol/m^3 input; it is whatever amount unit the
  concentration uses.)

Method
------
  * Every 3D cell is split into tetrahedra (vtkDataSetTriangleFilter); a wedge
    becomes 3 tets.  Placement is uniform in physical space even for tapered
    or twisted cells -- sampling a wedge's parametric coordinates directly
    would bunch points toward its narrow end.  Surfaces work the same way
    with triangles in place of tets: 2D cells are split into triangles and
    everything below applies with 3 corners instead of 4 and area instead
    of volume.  (Below, "tet" means "triangle" for surfaces.)
  * Amount in a tet = volume * concentration.  For point data the field is
    linear inside the tet, so amount = volume * mean of the corner values.
  * Expected points in a tet = amount / q.  The integer count is that value
    rounded down or up at random (stochastic rounding), so every tet gets
    floor or ceil of its expected count and the expected total is exact.
    The actual total varies slightly with the seed (roughly +/- sqrt of the
    number of tets with a fractional count; e.g. 491-509 for a target of 500
    on the Wavelet source).
  * Cell data: points are uniform inside the tet (barycentric weights from a
    Dirichlet(1,1,1,1) draw).  Point data, Linear placement: points follow
    the linear field exactly -- pick a corner in proportion to its value,
    then a barycentric point biased toward it (Dirichlet with that corner's
    parameter 2).  A zero-valued corner draws no extra points.  Point data,
    Uniform Per Cell: each input cell's exact amount is shared among its
    tets by volume, and points are placed as for cell data.
  * All random numbers come from a counter-based hash of (seed, tet, point
    index), not from a stream.  In a time series on a fixed mesh, a cell
    whose count does not change keeps exactly the same points, and a cell
    gaining a point keeps its existing ones, so points appear or vanish only
    where the material changes.  This holds for cell data and for point data
    with Uniform Per Cell placement.  With Linear placement, a point also
    shifts when its tet's corner values change (and can jump toward another
    corner).  If the mesh itself changes between steps, the points are
    re-drawn.

Element types and approximations
--------------------------------
Volumes: any vtkDataSet with 3D cells -- tets, hexes, wedges, pyramids,
voxels, mixed meshes, structured and image data; polyhedra via VTK's own
splitting (untested, non-convex cells may be slightly off).
Surfaces: polydata or unstructured grids of triangles, quads, polygons and
triangle strips, and 2D structured / image data (one-cell-thick extents);
curved surfaces such as spheres are sampled on their flat facets.
Exact for cells with planar faces (planar quads on surfaces).  Non-planar
quads are treated as two triangles, and quadratic/higher-order cells as
straight-edged simplices with a piecewise-linear field -- small differences,
invisible in a point cloud.  Inside a wedge, VTK's own interpolation is not
exactly linear over the 3 tets; the sampler follows the tet-linear version.

Sampling a volume concentration on a slice or boundary surface (rather than
a true per-area field) works for visualization -- density follows
concentration across the surface -- but q is then "per m^2 at 1 m
thickness", not a physical amount.

Data handling
-------------
Negative and NaN concentrations are treated as zero (negative material cannot
be sampled); the amount removed is reported as ClampedNegative_mmol.
Errors and warnings go to the Output Messages window, and on error the output
is cleared rather than left showing a previous result.

Limitations
-----------
  * Multiblock input (e.g. Exodus) needs Merge Blocks first.
  * One concentration array per filter; for several materials use one filter
    each with the same q and color them separately.
  * Applying a transform to the concentration (e.g. a Calculator raising it
    to a power for contrast) makes the points no longer a fixed amount each;
    label such views as enhanced.
"""

import numpy as np

from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain, smhint
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.util import numpy_support as ns
from vtkmodules.vtkCommonCore import vtkPoints, vtkDoubleArray, vtkOutputWindow
from vtkmodules.vtkCommonDataModel import (
    VTK_TETRA,
    VTK_TRIANGLE,
    vtkCellArray,
    vtkDataObject,
    vtkDataSet,
    vtkPolyData,
)
from vtkmodules.vtkFiltersGeneral import vtkDataSetTriangleFilter


# --------------------------------------------------------------------------
# Counter-based random numbers: U(0,1] as a pure function of its keys.
# --------------------------------------------------------------------------
def _mix(x):
    """splitmix64 finalizer on a uint64 array (wraparound is intended)."""
    x = x + np.uint64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def _uniform(seed, tet, k, stream):
    """Deterministic uniform in (0, 1] keyed by (seed, tet id, point index, stream)."""
    with np.errstate(over="ignore"):
        h = _mix(_mix(np.full(tet.shape, seed, dtype=np.uint64)) ^ tet.astype(np.uint64))
        h = _mix(h ^ (k.astype(np.uint64) * np.uint64(8) + np.uint64(stream)))
    return ((h >> np.uint64(11)).astype(np.float64) + 1.0) * (1.0 / 9007199254740992.0)


_COUNT_STREAM = 7   # stream id for the rounding draw; 0-5 are used per point

_DIM_AUTO, _DIM_VOLUME, _DIM_SURFACE = 0, 1, 2   # Cell Dimension property values
_PLACE_LINEAR, _PLACE_CELL = 0, 1                 # Point Data Placement property values


# --------------------------------------------------------------------------
# Sampling core.  Module-level so other tools (AdvectionSampler.py) can import
# this file and reuse exactly the same splitting, counting and placement.
# --------------------------------------------------------------------------
class Simplices(object):
    """Result of simplicize(): the input split into tets or triangles.

    tg           vtkUnstructuredGrid from vtkDataSetTriangleFilter (input point
                 and cell data passed through; cell array SourceCellId added)
    dim          3 (tets) or 2 (triangles)
    sel          indices into tg of the simplices used
    conn         (n, dim+1) point ids of each simplex
    P            (npts, 3) float64 point coordinates
    source_cell  (n,) input cell id of each simplex
    warnings     messages for the caller to report
    """


def simplicize(ds, dimension=_DIM_AUTO):
    """Split ds into simplices of one dimension.  dimension is _DIM_AUTO,
    _DIM_VOLUME or _DIM_SURFACE (see the Cell Dimension property).  Raises
    ValueError if there is nothing of that dimension to sample."""
    tagged = ds.NewInstance()
    tagged.ShallowCopy(ds)
    src_ids = ns.numpy_to_vtk(np.arange(ds.GetNumberOfCells(), dtype=np.int64), deep=1)
    src_ids.SetName("SourceCellId")
    tagged.GetCellData().AddArray(src_ids)

    # 3D cells -> tets, 2D cells -> triangles (lines/vertices are dropped below).
    tri = vtkDataSetTriangleFilter()
    tri.TetrahedraOnlyOff()
    tri.SetInputData(tagged)
    tri.Update()
    tg = tri.GetOutput()

    out = Simplices()
    out.warnings = []
    if tg.GetNumberOfCells() == 0:
        raise ValueError("input has no 2D or 3D cells")
    types = ns.vtk_to_numpy(tg.GetCellTypesArray())
    n3 = int((types == VTK_TETRA).sum())
    n2 = int((types == VTK_TRIANGLE).sum())
    if dimension == _DIM_AUTO:
        dim = 3 if n3 else 2
        if n3 and n2:
            out.warnings.append("input has both 3D and 2D cells; sampling the 3D cells only "
                                "(set Cell Dimension to Surface to sample the 2D cells)")
    else:
        dim = 3 if dimension == _DIM_VOLUME else 2
    if (n3 if dim == 3 else n2) == 0:
        raise ValueError("input has no %s cells" % ("3D (volume)" if dim == 3 else "2D (surface)"))

    s = dim + 1   # corners per simplex: 4 = tet, 3 = triangle
    sel = np.nonzero(types == (VTK_TETRA if dim == 3 else VTK_TRIANGLE))[0]
    offsets = ns.vtk_to_numpy(tg.GetCells().GetOffsetsArray())
    out.conn = ns.vtk_to_numpy(tg.GetCells().GetConnectivityArray())[
        offsets[sel][:, None] + np.arange(s)]
    out.tg, out.dim, out.sel = tg, dim, sel
    out.P = ns.vtk_to_numpy(tg.GetPoints().GetData()).astype(np.float64)
    out.source_cell = ns.vtk_to_numpy(tg.GetCellData().GetArray("SourceCellId"))[sel]
    return out


def simplex_measure(P, conn, length_scale=1.0):
    """Tet volumes (m^3) or triangle areas (m^2), from the number of corners."""
    a, b, c = (P[conn[:, i]] for i in range(3))
    if conn.shape[1] == 4:
        vol = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), P[conn[:, 3]] - a)) / 6.0
    else:
        vol = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    return vol * length_scale ** (conn.shape[1] - 1)


def stochastic_counts(amount, q, seed):
    """Integer points per simplex: amount/q rounded down or up at random,
    keyed by simplex index so the same seed gives the same rounding."""
    ids = np.arange(len(amount), dtype=np.int64)
    u = 1.0 - _uniform(seed, ids, np.zeros(len(amount), np.int64), _COUNT_STREAM)  # [0,1)
    return np.floor(amount / q + u).astype(np.int64)


def expand_counts(n):
    """Per-simplex counts -> (simplex index, index within simplex) per point."""
    tet = np.repeat(np.arange(len(n), dtype=np.int64), n)
    k = np.arange(len(tet), dtype=np.int64) - np.repeat(np.cumsum(n) - n, n)
    return tet, k


def place_points(seed, P, conn, tet, k, corner_values=None):
    """Positions of points (simplex tet[i], index k[i]) and their barycentric
    weights.  corner_values None: uniform in the simplex.  Otherwise an
    (N, s) array of the field at each point's simplex corners, and points
    follow the linear field (corner chosen by value, then a Dirichlet draw
    biased toward it)."""
    N, s = len(tet), conn.shape[1]
    # Barycentric weights: Dirichlet(1,...,1) is uniform in the simplex.
    e = np.stack([-np.log(_uniform(seed, tet, k, st)) for st in range(s)], axis=1)
    if corner_values is not None and N:
        cum = np.cumsum(corner_values, axis=1)
        r = (1.0 - _uniform(seed, tet, k, 4)) * cum[:, s - 1]
        j = np.minimum((cum <= r[:, None]).sum(axis=1), s - 1)
        e[np.arange(N), j] += -np.log(_uniform(seed, tet, k, 5))
    w = e / e.sum(axis=1, keepdims=True)
    return np.einsum("nk,nki->ni", w, P[conn[tet]]), w


# Messages go to ParaView's Output Messages window (print() is easy to miss).
def _fail(msg):
    vtkOutputWindow.GetInstance().DisplayErrorText("Molar Point Sampler: " + msg + "\n")
    return 0


def _warn(msg):
    vtkOutputWindow.GetInstance().DisplayWarningText("Molar Point Sampler: " + msg + "\n")


def _default_array(ds):
    """Active point scalars, else active cell scalars, else first 1-component array."""
    P, C = vtkDataObject.FIELD_ASSOCIATION_POINTS, vtkDataObject.FIELD_ASSOCIATION_CELLS
    for attrs, assoc in ((ds.GetPointData(), P), (ds.GetCellData(), C)):
        s = attrs.GetScalars()
        if s is not None and s.GetNumberOfComponents() == 1 and s.GetName():
            return s.GetName(), assoc
    for attrs, assoc in ((ds.GetPointData(), P), (ds.GetCellData(), C)):
        for i in range(attrs.GetNumberOfArrays()):
            a = attrs.GetArray(i)
            if a is not None and a.GetNumberOfComponents() == 1 and a.GetName():
                return a.GetName(), assoc
    return None, None


@smproxy.filter(label="Molar Point Sampler")
@smhint.xml('<ShowInMenu category="Sampling" />')
@smproperty.xml("""
  <Documentation
      short_help="Point cloud where each point is a fixed amount of material."
      long_help="Samples a concentration field into points whose density is proportional to concentration.">
    Turns a concentration field (amount per cubic meter, e.g. mmol/m^3, on a
    volume mesh; amount per square meter on a surface) into a point cloud in
    which every point stands for the same amount of material, Amount Per
    Point (q). Point density is proportional to concentration, and materials
    sampled with the same q get point counts in proportion to their total
    amount. Expected points = total amount / q; the total is reported in the
    output field data (TotalAmount_mmol), so for N points use
    q = TotalAmount_mmol / N. Works on any 3D cell type (split into
    tetrahedra) and any 2D cell type (split into triangles). See the header
    of MolarPointSampler.py for the method, outputs and limitations.
  </Documentation>""")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class MolarPointSampler(VTKPythonAlgorithmBase):
    """ParaView filter: concentration field -> point cloud, q per point.

    Properties are stored in the attributes set in __init__; RequestData
    does all the work.  See the module docstring for the full description.
    """

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=1, outputType="vtkPolyData")
        self._q = 1.0
        self._length_scale = 1.0
        self._seed = 12345
        self._max_points = 5000000
        self._pass_arrays = True
        self._dimension = _DIM_AUTO
        self._placement = _PLACE_LINEAR

    # ---------------------------------------------------------------- props
    @smproperty.xml("""
      <StringVectorProperty name="Concentration"
                            label="Concentration"
                            command="SetInputArrayToProcess"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list" attribute_type="Scalars"
                         number_of_components="1">
          <RequiredProperties>
            <Property name="Input" function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
        <Documentation>Point or cell array holding concentration, as amount per cubic meter (e.g. mmol/m^3) for volumes or per square meter (mmol/m^2) for surfaces. Amount Per Point is then in the same amount unit. Negative and NaN values are treated as zero.</Documentation>
      </StringVectorProperty>""")
    def SetInputArrayToProcess(self, idx, port, connection, field, name):
        return super().SetInputArrayToProcess(idx, port, connection, field, name)

    @smproperty.doublevector(name="AmountPerPoint", label="Amount Per Point",
                             default_values=1.0)
    @smdomain.xml("<Documentation>q, in mmol: every output point represents this many mmol. "
                  "Expected points = total amount / q, so for N points use "
                  "q = TotalAmount_mmol (output field data) / N. "
                  "Use the same value for every material or time step you want to compare."
                  "</Documentation>")
    def SetAmountPerPoint(self, q):
        if q != self._q:
            self._q = q
            self.Modified()

    @smproperty.intvector(name="CellDimension", label="Cell Dimension", default_values=0)
    @smdomain.xml("""
        <EnumerationDomain name="enum">
          <Entry value="0" text="Auto" />
          <Entry value="1" text="Volume" />
          <Entry value="2" text="Surface" />
        </EnumerationDomain>
        <Documentation>Which cells to sample. Volume: 3D cells, concentration per m^3.
        Surface: 2D cells (triangles, quads, polygons), concentration per m^2.
        Auto: 3D cells if the input has any, otherwise 2D cells.</Documentation>""")
    def SetCellDimension(self, d):
        if d != self._dimension:
            self._dimension = d
            self.Modified()

    @smproperty.intvector(name="PointDataPlacement", label="Point Data Placement",
                          default_values=0)
    @smdomain.xml("""
        <EnumerationDomain name="enum">
          <Entry value="0" text="Linear" />
          <Entry value="1" text="Uniform Per Cell" />
        </EnumerationDomain>
        <Documentation>Only used when Concentration is a point array. Linear: points
        follow the interpolated field inside each cell. Uniform Per Cell: each cell's
        amount is computed exactly from the point data, but its points are spread
        uniformly in the cell, as for cell data. Uniform Per Cell keeps points still
        from one time step to the next (only adding or removing them), at the cost
        of the detail inside each cell; best when cells hold few points each.
        Totals are the same either way.</Documentation>""")
    def SetPointDataPlacement(self, p):
        if p != self._placement:
            self._placement = p
            self.Modified()

    @smproperty.doublevector(name="LengthScale", label="Length Scale",
                             default_values=1.0)
    @smdomain.xml("<Documentation>Meters per mesh unit: coordinates times this give meters. "
                  "1 for meters, 0.001 for millimeters. Point counts scale with its cube "
                  "for volumes and its square for surfaces.</Documentation>")
    def SetLengthScale(self, s):
        if s != self._length_scale:
            self._length_scale = s
            self.Modified()

    @smproperty.intvector(name="Seed", default_values=12345)
    @smdomain.xml("<IntRangeDomain name='range' min='0' max='2147483647'/>"
                  "<Documentation>Selects which random points are drawn; does not change "
                  "how many on average. The same seed and input always give the same "
                  "points, and on a fixed mesh unchanged cells keep their points from one "
                  "time step to the next.</Documentation>")
    def SetSeed(self, seed):
        if seed != self._seed:
            self._seed = seed
            self.Modified()

    @smproperty.intvector(name="MaxPoints", label="Max Points", default_values=5000000)
    @smdomain.xml("<Documentation>Safety limit. If the number of points would exceed "
                  "this, the filter stops, clears its output and reports in Output Messages "
                  "the smallest Amount Per Point that fits.</Documentation>")
    def SetMaxPoints(self, n):
        if n != self._max_points:
            self._max_points = n
            self.Modified()

    @smproperty.intvector(name="PassArrays", label="Pass Arrays", default_values=1)
    @smdomain.xml("<BooleanDomain name='bool'/>"
                  "<Documentation>Carry the input's other arrays onto the points. Point arrays "
                  "are interpolated at each point (integer arrays take the nearest corner's "
                  "value); cell arrays are copied from the cell the point lies in. If a point "
                  "and a cell array share a name, the cell one gets a _cell suffix.</Documentation>")
    def SetPassArrays(self, on):
        on = bool(on)
        if on != self._pass_arrays:
            self._pass_arrays = on
            self.Modified()

    # ----------------------------------------------------------------- data
    @staticmethod
    def _pass_input_arrays(tg, out, conc_name, conc_on_points, corners, w, cells):
        """Interpolate point arrays with the sample's barycentric weights and copy
        cell arrays from the source simplex (cells = its index in tg) onto the
        output points."""
        N = len(cells)
        nearest = w.argmax(axis=1)
        skip = {"SourceCellId", "vtkGhostType", "vtkValidPointMask"}
        used = {conc_name, "SourceCellId"}
        opd = out.GetPointData()

        def add(values, nm, like):
            va = ns.numpy_to_vtk(np.ascontiguousarray(values), deep=1)
            va.SetName(nm)
            va.SetNumberOfComponents(like.GetNumberOfComponents())
            for c in range(like.GetNumberOfComponents()):
                if like.GetComponentName(c):
                    va.SetComponentName(c, like.GetComponentName(c))
            opd.AddArray(va)
            used.add(nm)

        for attrs, is_point in ((tg.GetPointData(), True), (tg.GetCellData(), False)):
            for i in range(attrs.GetNumberOfArrays()):
                a = attrs.GetArray(i)          # None for string / non-numeric arrays
                if a is None or not a.GetName() or a.GetName() in skip:
                    continue
                if a.GetName() == conc_name and is_point == conc_on_points:
                    continue                   # already output (clamped)
                try:
                    v = ns.vtk_to_numpy(a)
                except Exception:
                    continue                   # e.g. bit arrays
                if is_point:
                    vc = v[corners]            # (N, 4[, ncomp])
                    if np.issubdtype(v.dtype, np.floating):
                        r = np.einsum("nk,nk...->n...", w, vc).astype(v.dtype)
                    else:
                        r = vc[np.arange(N), nearest]
                else:
                    r = v[cells]
                nm = a.GetName()
                if nm in used:
                    nm += "_point" if is_point else "_cell"
                add(r, nm, a)

    def RequestData(self, request, inInfo, outInfo):
        """Tetrahedralize, compute amount per tet, round to counts, place
        points, attach arrays and field data.  Returns 0 (and an empty
        output) on error."""
        inp = vtkDataSet.GetData(inInfo[0], 0)
        out = vtkPolyData.GetData(outInfo, 0)

        out.Initialize()   # never leave a stale result on screen after an error

        ainfo = self.GetInputArrayInformation(0)
        name = ainfo.Get(vtkDataObject.FIELD_NAME()) if ainfo.Has(vtkDataObject.FIELD_NAME()) else None
        assoc = ainfo.Get(vtkDataObject.FIELD_ASSOCIATION()) if name else None
        if not name:
            # Nothing selected (e.g. filter created before its input had data):
            # fall back to the active scalars, points first, then cells.
            name, assoc = _default_array(inp)
            if not name:
                return _fail("no concentration array selected and the input has no scalar array")
            _warn("no concentration array selected; using %s array '%s'"
                  % ("point" if assoc == vtkDataObject.FIELD_ASSOCIATION_POINTS else "cell", name))
        if not self._q > 0:
            return _fail("Amount Per Point must be > 0")
        on_points = assoc == vtkDataObject.FIELD_ASSOCIATION_POINTS
        if not on_points and assoc != vtkDataObject.FIELD_ASSOCIATION_CELLS:
            return _fail("concentration must be a point or cell array")

        try:
            sx = simplicize(inp, self._dimension)
        except ValueError as e:
            return _fail(str(e))
        for msg in sx.warnings:
            _warn(msg)
        tg, dim, sel, conn, P, source_cell = sx.tg, sx.dim, sx.sel, sx.conn, sx.P, sx.source_cell
        s = dim + 1
        ntet = len(sel)

        attrs = tg.GetPointData() if on_points else tg.GetCellData()
        arr = attrs.GetArray(name)
        if arr is None or arr.GetNumberOfComponents() != 1:
            return _fail("'%s' missing or not a single-component array" % name)
        raw = np.nan_to_num(ns.vtk_to_numpy(arr).astype(np.float64), nan=0.0)
        if not on_points:
            raw = raw[sel]
        conc = np.maximum(raw, 0.0)

        vol = simplex_measure(P, conn, self._length_scale)

        if on_points:
            cn = conc[conn]                       # (n simplices, s) corner values
            amount = vol * cn.mean(axis=1)
            clamped = (amount - vol * raw[conn].mean(axis=1)).sum()
            if self._placement == _PLACE_CELL:
                # Keep each input cell's exact amount, but spread it over the
                # cell's simplices by size, i.e. uniform density within the cell.
                cell_amount = np.bincount(source_cell, weights=amount)
                cell_vol = np.bincount(source_cell, weights=vol)
                share = np.divide(vol, cell_vol[source_cell],
                                  out=np.zeros_like(vol), where=cell_vol[source_cell] > 0)
                amount = cell_amount[source_cell] * share
        else:
            amount = vol * conc
            clamped = (amount - vol * raw).sum()

        n = stochastic_counts(amount, self._q, self._seed)
        N = int(n.sum())
        if N > self._max_points:
            return _fail("%d points exceeds Max Points (%d); "
                         "raise Amount Per Point to at least %.4g mmol"
                         % (N, self._max_points, amount.sum() / self._max_points))

        tet, k = expand_counts(n)
        cp = cn[tet] if on_points else None
        linear = on_points and self._placement == _PLACE_LINEAR
        pos, w = place_points(self._seed, P, conn, tet, k, cp if linear else None)
        value = (w * cp).sum(axis=1) if (on_points and N) else conc[tet]

        # Output: vertices with the sampled value and the source cell id.
        pts = vtkPoints()
        pts.SetData(ns.numpy_to_vtk(pos, deep=1))
        out.SetPoints(pts)
        verts = vtkCellArray()
        verts.SetData(ns.numpy_to_vtkIdTypeArray(np.arange(N + 1, dtype=np.int64), deep=1),
                      ns.numpy_to_vtkIdTypeArray(np.arange(N, dtype=np.int64), deep=1))
        out.SetVerts(verts)

        va = ns.numpy_to_vtk(value, deep=1)
        va.SetName(name)
        out.GetPointData().AddArray(va)
        out.GetPointData().SetActiveScalars(name)
        if self._pass_arrays and N:
            self._pass_input_arrays(tg, out, name, on_points, conn[tet], w, sel[tet])
        sa = ns.numpy_to_vtk(source_cell[tet].astype(np.int64), deep=1)
        sa.SetName("SourceCellId")
        out.GetPointData().AddArray(sa)

        for fname, fval in (("TotalAmount_mmol", amount.sum()),
                            ("AmountPerPoint_mmol", self._q),
                            ("ClampedNegative_mmol", clamped),
                            ("SampledDimension", dim)):
            fa = vtkDoubleArray()
            fa.SetName(fname)
            fa.InsertNextValue(fval)
            out.GetFieldData().AddArray(fa)
        return 1
