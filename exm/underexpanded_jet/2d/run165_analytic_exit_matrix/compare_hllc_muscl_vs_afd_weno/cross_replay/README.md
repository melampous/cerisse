# Cross-replay isolation

This directory isolates reconstruction effects from accumulated solution
differences.  The two HLLC operators are cross-restarted at matched physical
time:

- HLLC--MUSCL reads the AFD--HLLC--WENO checkpoint.
- AFD--HLLC--WENO reads the HLLC--MUSCL checkpoint.

Both runs retain the fixed coarse time step `4e-7 s`, the L1 substep
`2e-7 s`, and the original Run165 RZ/no-IBM configuration.
