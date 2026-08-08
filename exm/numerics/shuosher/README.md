# Shu--Osher shock--entropy-wave interaction

This directory contains the one-dimensional inviscid benchmark used for the
Chapter 2 shock-capturing check.  `inputs.thesis` uses the canonical domain
`[-5,5]`, places the Mach-3 shock at `x=-4`, and advances to `t=1.8`.

The thesis baseline is deliberately fixed to the same ingredients as the main
Euler calculations:

- single-species calorically perfect gas, with no combustion or turbulence;
- source-free Euler equations (`no_diffusive_t`, `no_source_t`);
- characteristic LLF flux splitting and WENO-Z5 reconstruction;
- SSPRK(3,3) with `CFL=0.3`;
- a prescribed post-shock state at the left boundary and analytic continuation
  of the undisturbed entropy wave at the right boundary.

The last point replaces the historical two-sided `foextrap` setup, which
created a visible left-boundary artefact.  The older `inputs.shu` and
`plotREF` are retained only for legacy comparisons; `plotREF` is an old
N=8912 MUSCL result and must not be used as the thesis reference.

Run the complete uniform-grid study with

```bash
./run_study.sh
```

This produces N=200, 400 and 800 results, plus an N=6400 LLF--WENO-Z5
reference generated with the same boundary treatment.  `analyze.py` reports
density L1, L2 and Linf differences, total variation, and the solver's
time-advance cost.  The N=200-based sequence keeps the initial shock exactly
on a grid face.  The errors are resolution differences for a discontinuous
solution; they must not be presented as a formal order-of-accuracy result.
The script also writes a profile figure and CSV table.  For a
timing table, repeat the runs on an otherwise idle node and report the compiler,
CPU, MPI rank count and the `Run Time advance` value from each `run.log`.

For a quick smoke test:

```bash
mkdir -p smoke
make -j4
./main1d.gnu.MPI.shu_llf-wenoz5.ex inputs.thesis \
  amr.n_cell=128 stop_time=0.02 amr.plot_file=./smoke/plt
```

Alternative schemes remain buildable for method-comparison figures, for
example `make SHU_EULER_SCHEME=afd-hllc-wenoz5`, but those executables are not
the Chapter 2 baseline. The formal compact reconstruction driver is selected
with `llf-wenoz5`, `llf-teno5`, or `llf-teno6`. The last two selectors are
all-fluid verification options and are not production IBM schemes. Explicit
`-old` selectors retain the archived implementation for regression only.
