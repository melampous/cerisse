# Legacy isentropic-vortex prototype

This directory contains the original three-dimensional KEEP test.  Its current
initialisation does not satisfy the radial pressure balance of an exact
isentropic vortex, and it must not be used for the thesis verification.

Use the corrected two-dimensional LLF--WENO-Z5 and SSPRK(3,3) case in
[`exm/numerics/isentropic_vortex`](../numerics/isentropic_vortex/README.md).
That directory contains uniform-grid and fixed-AMR inputs, an analytic-error
tool and reproducible run scripts.  The legacy files here are retained only to
avoid silently changing historical examples.
