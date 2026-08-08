# AFD HLLC shock fallback A/B test

This case tests a stationary Mach 10 normal shock on a uniform two-dimensional
Cartesian grid.  The transverse direction is periodic.  The initial shock has
an alternating subcell displacement of 0.01 dx.  Exact cell averages are used
in the intersected cells.

One AFD HLLC WENO-Z5 executable is used for every run.  Grid, CFL, Runge Kutta
scheme, final time, initial condition, boundary conditions, reconstruction,
and AFD correction are identical within each pair.  The only changed option is
`cns.afd_shock_llf=0` or `1`.
The enabled branch uses characteristic high-order LLF on sensor-marked
non-smooth stencils that also pass a relative pressure-jump threshold of 0.01
and a normalized compression threshold of 0.001.  It is distinct from the
cell-centred first-order LLF invalid-state safeguard that remains active in
both runs.

An unperturbed 200 by 20 planar control is run for both switch values.  The
seeded comparison uses 200 by 20 and 400 by 40 grids.  The study records
shock-front peak-to-peak distortion, the odd-even Fourier component,
transverse velocity, minimum density and pressure, and complete time
histories.  The offline diagnostic mirrors all three parts of the runtime
sensor.  It reports the broad non-smooth count, pressure and compression gate
counts, and their final intersection for each face direction.  The seeded
initial state must mark all transverse faces at one streamwise column.  This
gives 20 faces on the coarse grid and 40 faces on the fine grid.  The planar
control must mark no transverse face.  It is a focused multidimensional
shock-stability test.  It does not by itself qualify a method for every
supersonic flow.

The frozen 2026-07-29 result and the proposed figure caption are documented
in [RESULTS.md](RESULTS.md).

Build and run:

```bash
make -j8
./run_ab.sh ./main2d.gnu.MPI.afd_gridshock_wenoz5.ex ./inputs \
  ./results/afd_shock_llf_ab
```
