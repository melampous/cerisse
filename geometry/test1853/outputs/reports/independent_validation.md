# Independent geometry validation

| Model | STEP solid | BREP solid | STL watertight | STL volume | GTS self-intersection | Triangles |
|---|---:|---:|---|---|---|---:|
| `test1853_run165_single` | 1 | 1 | True | True | True | 353042 |
| `test1853_run262_263_tri_n3at120` | 1 | 1 | True | True | True | 360834 |
| `test1853_run262_263_tri_n3at240` | 1 | 1 | True | True | True | 360872 |
| `test1853_quad_n3at120` | 1 | 1 | True | True | True | 364826 |
| `test1853_quad_n3at240` | 1 | 1 | True | True | True | 364880 |

All listed models passed CAD re-import, bounds, positive-volume, STL manifold/winding, duplicate/degenerate-face, patch-area, and (unless explicitly skipped) GTS self-intersection gates.

SHA-256 hashes and numeric diagnostics are in `independent_validation.json`.
