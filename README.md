# ParaView plugins

## Molar Point Sampler

`MolarPointSampler.py`: a ParaView Python filter that turns a concentration
field (amount per m³ on a volume mesh, or per m² on a surface) into a point
cloud in which every point stands for the same amount of material,
**Amount Per Point** (q). Point density follows concentration, and materials
sampled with the same q get point counts in proportion to their total amount.

- Works on any 3D cell type (split into tets) and any 2D cell type (split
  into triangles), with point or cell data.
- Points are stable over time on a fixed mesh: unchanged cells keep their
  points. For point data, set **Point Data Placement** to *Uniform Per Cell*
  to get the same stability (totals are unchanged; detail inside each cell is
  lost).
- Can carry the input's other arrays onto the points.

The header of `MolarPointSampler.py` documents the properties, output, method
and limitations.

**Load:** in ParaView, go to Tools → Manage Plugins → Load New… and pick the
file. The filter is then under Filters → Sampling → Molar Point Sampler.

**For about N points:** run once, read `TotalAmount_mmol` from the output's
field data, and set q = TotalAmount_mmol / N.

## Advecting sampler

`AdvectionSampler.py`: an offline pvpython tool for a time series with a
velocity field and a concentration field, such as smoke from a fire in a
breeze. It writes a point-cloud time series in which each point stands for
about q of material, **moves with the flow**, and keeps its identity from
frame to frame:

- Points are advected with RK4 through velocity interpolated in space and
  time.
- At each saved step, points are reconciled with the field in bins. Small
  mismatches adjust weights. New material (for example from a fire) arrives
  as new points emitted along their paths during the interval. Points in
  clean air fade out.
- Each point carries `ParticleId`, `Weight`, `Fade` (for opacity), `Age`,
  `BirthTime` and `BirthPosition`, the last two showing where and when the
  material was emitted.
- A per-step `Residual` grid shows where the field changed in ways the flow
  doesn't explain (sources, sinks, unresolved mixing), and `diagnostics.csv`
  records how much correction each step needed.

```bash
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython AdvectionSampler.py \
    data/fire_plume/fire_plume.pvd --velocity velocity --concentration smoke \
    --q 8 --diffusivity 2 --out out/plume
```

Then open `out/plume/particles.pvd` in ParaView. The header of
`AdvectionSampler.py` and `--help` describe the method, options and output.
It reuses the Molar Point Sampler's sampling code, so keep both files
together.

**Tuning:** check the `explained` and `churn` columns of `diagnostics.csv`.
If the data diffuses more than the flow alone explains (coarse grids,
modelled turbulence), a modest `--diffusivity` raises `explained` and reduces
culls; too much lowers it again. `--max-growth` trades splits (pieces placed
beside their parent) for births (emitted along paths).

## Synthetic test data

`synthetic/make_fire_plume.py` solves advection–diffusion for smoke from a
pulsing fire in a gusty, veering breeze with a tilted updraft, so the
velocity and concentration are consistent with each other the way real
simulation output is. It takes about 10 s and writes 23 MB:

```bash
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython synthetic/make_fire_plume.py --out data/fire_plume
```

## Tests

Tested with ParaView 5.13.3 and 6.1.1:

```bash
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython tests/test_molar_point_sampler.py
```

```bash
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython tests/test_advection_sampler.py
```

The second test takes about a minute. It checks exact transport (solid-body
rotation, where the answer is known analytically), that the probe backend
for arbitrary meshes matches the image backend, and the fire plume.
