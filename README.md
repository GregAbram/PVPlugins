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
  points.
- Can carry the input's other arrays onto the points.

The header of `MolarPointSampler.py` documents the properties, output, method
and limitations.

**Load:** in ParaView, go to Tools → Manage Plugins → Load New… and pick the
file. The filter is then under Filters → Sampling → Molar Point Sampler.

**For about N points:** run once, read `TotalAmount_mmol` from the output's
field data, and set q = TotalAmount_mmol / N.

**Tests** (tested with ParaView 5.13.3 and 6.1.1):

```bash
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython tests/test_molar_point_sampler.py
```
