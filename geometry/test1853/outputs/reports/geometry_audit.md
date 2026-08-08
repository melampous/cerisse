# Test 1853 geometry build audit

All CAD coordinates are millimetres. `PASS` means that the generated artifact met the checks below; it does not turn a documented public-data reconstruction into the unavailable NASA original CAD.

| Model | OCC solids | STL triangles | Watertight | Volume | Max CAD–STL chord sample |
|---|---:|---:|---|---|---:|
| `test1853_run165_single` | 1 | 353042 | True | True | 0.014864 mm |
| `test1853_run262_263_tri_n3at120` | 1 | 360834 | True | True | 0.014864 mm |
| `test1853_run262_263_tri_n3at240` | 1 | 360872 | True | True | 0.014864 mm |
| `test1853_quad_n3at120` | 1 | 364826 | True | True | 0.014864 mm |
| `test1853_quad_n3at240` | 1 | 364880 | True | True | 0.014864 mm |

## Interpretation boundary

- Run 165 is the center single-nozzle configuration. Runs 247/262/263 use the same tri-nozzle hardware with run-specific flow and roll conditions.
- Runs 307/315 use the quad hardware: center Nozzle 1 plus peripheral Nozzles 2/3/4. The phi=0 CAD/STL masters remain run-neutral; the preparation tool applies the selected run roll.
- The body OML is the documented nominal tangent reconstruction. The public drawings explicitly require missing CAD surfaces.
- Peripheral virtual-exit axial placement is an explicit CFD reconstruction assumption; the physical lip in these solids is the Boolean cone/OML intersection.
- Both Nozzle-3/Nozzle-4 azimuth assignments are exported because the public report does not identify their one-to-one mapping.
- The throat caps are CFD boundary patches. Complete convergent passages, internal plenum, installed sting services, and pressure-port holes are not guessed.

Detailed numeric results are in `build_audit.json`.
