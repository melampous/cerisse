# Periodic advected contact

This test compares `cns.afd_shock_llf=0` and `1` for the same
AFD HLLC WENO-Z5 executable.  It uses a periodic square density pulse with
constant pressure and velocity.  The two pulse edges are Euler contact
discontinuities.

The test advances one exact transit on uniform grids with 64, 128, 256 and
512 cells.  Both switch values use CFL 0.20, SSPRK(3,3), identical initial
data and the same executable.  The only run-time difference is
`cns.afd_shock_llf`.

The analysis reports density errors, total variation, mixed-cell count,
pressure error, conservation error and an offline reproduction of the AFD
shock gate.  It reports the broad curvature sensor, the relative pressure
jump, the compression test and their final intersection separately.  The
thresholds are 0.08, 0.01 and 0.001, respectively.  These offline counts are
snapshot diagnostics.  They are not time-integrated LLF flux counts.

Run with

```bash
./run_afd_shock_llf_ab.sh
```

Set `CONTACT_AB_RESULTS_ROOT` to select a different new result directory.
The script refuses to overwrite an existing archive.
