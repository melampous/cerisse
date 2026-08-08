# Shu--Osher four-scheme matrix

`run_method_matrix.sh` compares the following Cartesian, single-species Euler
schemes without changing the LLF--WENO-Z5 baseline in `run_study.sh`.

- `llf-wenoz5`
- `llf-teno5`
- `afd-hllc-wenoz5`
- `afd-hllc-teno5`

The script first builds and smoke tests all four executables.  It starts the
full calculations only after every smoke test has passed.  Each method is then
run at `N=200`, `400` and `800` to `t=1.8` with the common `inputs.thesis`
configuration.  The AFD runs explicitly select

```text
cns.afd_correction=1
cns.afd_shock_llf=1
cns.afd_smoothness_threshold=0.08
cns.afd_teno_cutoff=1.0e-4
```

Run directly on a single process with

```bash
CERISSE_BUILD_JOBS=8 ./run_method_matrix.sh
```

If the host requires an MPI launcher, use for example

```bash
CERISSE_BUILD_JOBS=8 \
CERISSE_RUN_PREFIX="mpirun -np 1" \
./run_method_matrix.sh
```

The default archive is `results/y9000x_method_matrix`.  The script refuses to
overwrite it.  Select another archive without editing the script:

```bash
SHU_MATRIX_ROOT=results/y9000x_method_matrix_repeat \
./run_method_matrix.sh
```

Every build and run directory contains the exact command, log and exit status.
The root manifest records the host, source revision, working-tree state and
compiler.  If `numpy`, `yt` and `matplotlib` are available, the script also
runs `analyze_method_matrix.py`.  Otherwise it records the deferred analysis
command so that post-processing can be performed after copying the archive.

The analysis uses each method's `N=800` result only as a resolution-difference
reference for its `N=200` and `N=400` profiles.  These differences must not be
reported as formal errors or as evidence that `N=800` is converged.  This
matrix intentionally does not generate the much larger `N=6400` reference used
by the separate baseline study.

## Skew--JST extension

`run_skew_jst_study.sh` adds the scheme name `skew` without changing the
existing four-scheme archives.  It builds the order-4 conservative
skew-symmetric central flux with JST dissipation.  The default coefficient
pair is `C2=1.5` and `C4=0.016`.  Both values are compile-time selectors.
They are included in the executable name and recorded in the run manifest.

The script performs a smoke test before the `N=200`, `400` and `800` runs.  It
also computes an `N=6400` Skew--JST profile for numerical comparison.  The
analysis reports density-profile differences in the post-shock window
`-1.5 <= x <= 1.5`.  The fine profile is labelled as a numerical reference.
It is not an exact solution and the analysis does not report a formal
convergence order.
