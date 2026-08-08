# Independent geometry validation

| Model | STEP solid | BREP solid | STL watertight | STL volume | GTS self-intersection | Triangles |
|---|---:|---:|---|---|---|---:|
| `test1853_single_peripheral_surface_normal` | 1 | 1 | True | True | True | 353080 |

All listed models passed CAD re-import, bounds, positive-volume, STL manifold/winding, duplicate/degenerate-face, patch-area, and (unless explicitly skipped) GTS self-intersection gates.

SHA-256 hashes and numeric diagnostics are in `independent_validation.json`.
