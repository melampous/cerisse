# Two-dimensional convecting isentropic vortex

This is the smooth, multidimensional Chapter 2 verification case.  It repairs
the older `exm/isentropic_vortex` prototype, which is a three-dimensional KEEP
test and imposes constant pressure rather than the isentropic radial balance.
The legacy directory is left unchanged.

The present case follows the C1.6 convecting-vortex formulation used by Monal
Patel.  The default is the fast vortex (`M=0.5`, `beta=0.2`) in a periodic
`0.1 m x 0.1 m` domain.  The default one-period grid sequence is intentionally
shorter than Patel's 50-period dissipation study.  The temperature perturbation,
density and pressure are constructed consistently, so the vortex is an exact
solution of the source-free Euler equations translated at the uniform stream
velocity.

The numerical configuration is:

- two-dimensional Cartesian grid and a single calorically perfect gas;
- no viscous flux, source term, combustion, turbulence model or IBM;
- characteristic LLF flux splitting with WENO-Z5 reconstruction for the
  thesis baseline;
- SSPRK(3,3), `CFL=0.3`, and periodic boundaries.

## Uniform-grid verification

Run

```bash
./run_study.sh
```

to advance N=40, 80 and 160 grids through one flow-through time.  The final
vortex returns to its initial position.  `analyze.py` evaluates the translated
analytic solution at every active cell centre and writes volume-weighted
density, pressure and velocity errors, mass and energy changes relative to the
discrete initial state, perturbation kinetic-energy retention and
`Run Time advance` to `results/uniform_errors.csv`.
The observed order is meaningful only for this smooth uniform-grid sequence;
with a time step proportional to the mesh spacing, SSPRK(3,3) limits the
asymptotic combined order to three.

For the longer 50-flow-through dissipation test used in the reference thesis,
override the stop time:

```bash
mkdir -p results/long_N120
./main2d.gnu.MPI.iv_llf-wenoz5.ex inputs.uniform \
  stop_time=2.879768262021454e-2 amr.n_cell="120 120" \
  amr.plot_file=./results/long_N120/plt
```

The official slow-vortex variant requires matching changes to the Mach number,
vortex strength, stop time and analysis parameters.  For one flow-through on
an N=80 grid, use

```bash
mkdir -p results/slow_N80
./main2d.gnu.MPI.iv_llf-wenoz5.ex inputs.uniform \
  prob.mach=0.05 prob.beta=0.02 stop_time=5.759536524042908e-3 \
  amr.n_cell="80 80" amr.plot_file=./results/slow_N80/plt
./analyze.py results/slow_N80 --mach 0.05 --beta 0.02 \
  --csv results/slow_N80_errors.csv
```

## Coarse--fine-interface checks

`inputs.amr` creates a fixed level-1 hierarchy over the right half of the
domain, enables refluxing, and starts the vortex in the coarse left half.  Run

```bash
./run_amr_study.sh
```

to measure the additional error and dissipation after the vortex crosses both
coarse--fine interfaces.  The standard input retains the CERISSE default of
limited linear coarse-to-fine interpolation.  Compare the AMR CSV with the
uniform-grid CSV at matched base resolution.

`run_amr_interpolation_matrix.sh` provides a separate smooth-flow diagnostic.
It compares the default linear interpolation with AMReX
`quartic_interp`.  Select the latter with
`cns.amr_state_interp=conservative_quartic`.  This option is restricted to
non-IBM and non-EB verification.  It does not change the production default.
The script also compares global time stepping with level subcycling and can
enable or disable refluxing.

## Relation to the reference implementations

The numerical paths in the reference theses are not identical.  Tin uses
point-state TENO5 as the principal high-order interpolation and also includes
WENO-Z5.  Enson's CERISSE1 source provides both choices.  These methods use
an HLLC Riemann solver and the alternative finite-difference derivative
correction.  This is the AFD--HLLC path in CERISSE.  Its manufactured-solution
result verifies that path.  It does not establish the formal order of the
separate characteristic LLF split-flux path used in this case.

Patel uses global Lax--Friedrichs flux splitting with TENO6 for the reported
AMR vortex study.  The calculation uses fourth-order prolongation, global
time stepping, no coarse--fine flux correction, a CFL number of about 0.75,
and 50 flow-through times.  The base grids are \(60^2\), \(90^2\), and
\(120^2\).  The reported AMR orders are 3.04 and 2.85.  The third-order time
integrator limits that long-time study.

The closest short CERISSE interface check is archived in
`results/amr_reference_config_q1_gts_reflux0_20260729_r1`.  It uses
LLF--WENO-Z5 with the Enson LLF exponent \(q=1\), a face-common scaled
epsilon of \(10^{-6}\), AMReX `quartic_interp`, global time stepping, and no
refluxing.  The grids are \(40^2\), \(80^2\), and \(160^2\).  The vortex is
advanced for one flow-through time at CFL 0.3.  The observed density
\(L_2\) orders are 3.343 and 3.923.  The velocity orders are 3.694 and 4.727.
The pressure orders are 3.785 and 4.566.  The scaled epsilon is a new
CERISSE verification option and is not taken from the reference theses.

This is a short-time assessment of the complete AMR interface treatment under
a Patel-aligned coupling configuration.  The interpolation matrix shows that
limited linear prolongation is an important source of error.  The recovered
rates also depend on the scale-aware LLF--WENO-Z5 regularisation.  The result
is not a reproduction of Patel's 50-flow-through TENO6 calculation.  The two
sets of rates should not be used to rank the schemes.

## Optional Euler-scheme comparison

`IV_EULER_SCHEME` selects one of six compile-time Euler schemes without
changing the baseline scripts:

```text
llf-wenoz5
llf-teno5
afd-hllc-wenoz5
afd-hllc-teno5
skew-central-o4
skew-jst-o4
```

The default is `llf-wenoz5`.  `run_study.sh` and `run_amr_study.sh` also pass
this choice explicitly.  They therefore remain the fixed thesis-baseline
runs. The two Skew selectors use the current order-four coefficient set.
`skew-central-o4` disables both JST terms. `skew-jst-o4` uses
`C2skew=1.5` and `C4skew=0.016`.

Use `run_method_matrix.sh` for the original four-method comparison.  A single invocation
can run one method and one suite, which is convenient for a remote workstation:

```bash
./run_method_matrix.sh llf-wenoz5 uniform
./run_method_matrix.sh llf-teno5 all
./run_method_matrix.sh afd-hllc-wenoz5 all
./run_method_matrix.sh afd-hllc-teno5 all
```

Use `all` as the first argument to run all four methods sequentially.  Set
`CERISSE_MPI_RANKS` to use `mpirun`, for example:

```bash
CERISSE_BUILD_JOBS=8 CERISSE_MPI_RANKS=4 \
  ./run_method_matrix.sh all all
```

Each method and suite has an independent directory below
`results/method_matrix`.  The script refuses to overwrite an existing suite.
It archives the input deck and a case-source snapshot. It also records the
source revision, host, executable checksum, exact replay command, and runtime
parameters. The AFD runs explicitly record and pass

```text
cns.afd_correction           = 1
cns.afd_shock_llf            = 1
cns.afd_smoothness_threshold = 0.08
cns.afd_teno_cutoff          = 1.0e-4
```

The AFD methods use HLLC in smooth regions and retain the configured
characteristic LLF fallback near detected non-smooth regions.  The comparison
should therefore be described as AFD HLLC with an LLF fallback.

## Skew and JST assessment

`run_skew_validation.sh` evaluates the current Skew operator without changing
`src/rhs/Skew.h`. It creates separate executables for the undamped central
operator and the complete Skew--JST method. The script runs the following
tests.

- A 45 degree low-temporal-error sequence on \(40^2\), \(80^2\), \(160^2\),
  and \(320^2\) grids for both executables.
- A one-transit sequence on \(40^2\), \(80^2\), and \(160^2\) grids for
  Skew--JST.
- A five-transit Skew--JST calculation on an \(80^2\) grid.

The dedicated analysis reports \(L_1\), \(L_2\), and \(L_\infty\) errors for
density, both velocity components, and pressure. It also reports perturbation
kinetic energy, peak vorticity, enstrophy, and core circulation. Numerical and
exact vorticity use the same periodic fourth-order diagnostic. The reported
equivalent viscosity is

\[
\nu_{\mathrm{eq},K}
=-\frac{R_v^2}{4t}\log\left(\frac{K^\prime}
{K^\prime_{\mathrm{exact}}}\right).
\]

This expression follows from the energy and enstrophy ratio of the Gaussian
vortex. It is a scale-specific, time-integrated dissipation indicator. It is
not a constant viscosity of the numerical method.

Run the complete assessment with

```bash
CERISSE_BUILD_JOBS=4 CERISSE_MPI_RANKS=4 \
  ./run_skew_validation.sh
```

Use `run_skew_one_transit_comparison.sh` to repeat LLF WENO-Z5, LLF TENO5,
AFD HLLC WENO-Z5, and Skew--JST from the same source snapshot. The script
creates a combined CSV and a three-panel PDF and PNG:

```bash
CERISSE_BUILD_JOBS=4 CERISSE_MPI_RANKS=4 \
  ./run_skew_one_transit_comparison.sh
```

Both drivers refuse to overwrite an existing result root. They archive the
case files, the active flux source, the worktree patch, executable hashes,
exact replay commands, build logs, run logs, and exit status.

## Low-temporal-error spatial convergence

`run_spatial_convergence.sh` provides a separate spatial convergence study for
AFD HLLC WENO-Z5. It uses the uniform grids
\(40^2\), \(80^2\), \(160^2\), and \(320^2\). The vortex is advanced through
0.1 domain lengths. The respective step counts are 64, 256, 1024, and 4096.
The final physical time is identical on every grid and
\(\Delta t\) is proportional to \(\Delta x^2\). The global SSPRK(3,3) error is
therefore \(O(\Delta x^6)\), one order smaller than the nominal fifth-order
spatial error.

Run the study with

```bash
CERISSE_MPI_RANKS=4 \
CERISSE_IV_SPATIAL_RESULTS_ROOT=results/afd_hllc_wenoz5_spatial_convergence \
./run_spatial_convergence.sh
```

The result directory is created once and is never overwritten. The manifest
records the executable hash, grid and step sequences, time-step scaling, and
all AFD runtime options.

The mean-flow angle is controlled by `CERISSE_IV_FLOW_ANGLE_DEG`. The formal
spatial study defaults to 45 degrees. This diagonal case has the same total
Mach number and vortex strength as the axial case:

```bash
CERISSE_MPI_RANKS=4 \
CERISSE_IV_SPATIAL_RESULTS_ROOT=results/afd_hllc_wenoz5_spatial_angle45 \
./run_spatial_convergence.sh
```

The diagonal case keeps the HLLC contact-wave branch away from zero in both
coordinate directions. It is used to assess the nominal order of the smooth
AFD branch. Set `CERISSE_IV_FLOW_ANGLE_DEG=0.0` for the separate axial
contact-wave sensitivity test.

The archived 45 degree study on grids from \(40^2\) to \(320^2\) gave final
density, velocity and pressure \(L_2\) orders of 4.972, 5.146 and 4.988. A
repeat at \(N=160\) with half the time step changed the density error by only
0.013 per cent. A sensor-disabled repeat produced bitwise identical
conservative fields.

The axial companion converged through \(N=240\), but its density \(L_2\) error
increased at \(N=320\). The extra error was local to the vortex core and was
nearly isobaric. It remained after changes to the time step, grid-box layout
and sensor setting. The result is retained as a contact-branch sensitivity.
It is not used to claim uniform fifth-order accuracy.

After all eight suites have produced `errors.csv`, combine the results and
create the uniform-grid and AMR comparison figure with

```bash
python3 analyze_method_matrix.py results/method_matrix
```

The default outputs are `method_matrix_summary.csv` and
`method_matrix_accuracy.png` in the matrix root.  The combined CSV reports
observed density-error rates for both suites and matched AMR-to-uniform error,
cell-count and runtime ratios.  The AMR rates assess the complete AMR
algorithm.  They are not the formal order of the interior reconstruction.

Create a four-method flow-field comparison at the finest requested resolution
with

```bash
python3 plot_flowfield_matrix.py results/method_matrix
```

The upper row shows the uniform-grid density perturbation and perturbation
velocity.  The lower row shows the signed density error after the vortex has
crossed the AMR interfaces.  AMR errors are evaluated at the true leaf-cell
centres and then rendered over each leaf-cell footprint.  No artificial
fine-grid sampling is introduced in the coarse part of the domain.  Both PNG
and PDF outputs are created.  The same command also creates a zoomed
\(N=40\) dimensionless-vorticity comparison.  This coarser field makes the
different levels of numerical dissipation visible.  Vorticity is evaluated
with a periodic centred diagnostic and is not used in the error norms.
The dashed line marks the interior interface at \(x/L=0.5\).  The periodic
domain boundary at \(x/L=0/1\) forms the second coarse-fine interface.

For timing-only repetitions, set `amr.plot_files_output=0` and
`cns.check_state=0`, repeat on the same idle node, and report the median
`Run Time advance`.  Accuracy runs should retain the state check and final
plotfile used above.
