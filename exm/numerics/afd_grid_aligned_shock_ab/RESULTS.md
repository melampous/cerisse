# Mach 10 grid-aligned shock A/B result

> Historical baseline for the broad pre-gate sensor. The authoritative
> joint-sensor rerun is
> `results/afd_joint_sensor_20260729_r2` and is summarised in
> `../AFD_JOINT_SHOCK_SENSOR_20260729.md`.

## Frozen configuration

All six calculations use the same
`main2d.gnu.MPI.afd_gridshock_wenoz5.ex` executable.  Its SHA256 is
`d2be03db662d3e60d5d353456592d625117cdd22b8cbaca17324b547bd8031b3`.
The runs use a Cartesian Level 0 mesh, CFL 0.25, SSPRK3, AFD correction,
WENO-Z5 point reconstruction, and a final time of 0.05.  The sole
method change within each pair is `cns.afd_shock_llf=0` or `1`.
When enabled, this option selects the characteristic high-order LLF branch
on stencils marked non-smooth by the existing AFD sensor.  The cell-centred
first-order LLF branch remains an independent invalid-state safeguard in both
runs.

The frozen source hashes are recorded in
`results/afd_shock_llf_ab_20260729_r1/FROZEN_CORE_SOURCE.sha256`.  They match
the shared A/B source snapshot.

## Result

The unperturbed 200 by 20 controls preserve a planar front.  Their measured
odd-even amplitude and transverse velocity remain zero for both switch
values.  The seeded runs start from the same alternating front displacement.
The initial measured peak-to-peak displacement is 0.01385 cells.

At 200 by 20 cells, the final peak-to-peak shock distortion is 0.1791 cells
with the fallback disabled and 0.02465 cells with the fallback enabled.  The
fallback therefore reduces the final cell-normalized distortion by a factor
of 7.26.

At 400 by 40 cells, the corresponding values are 0.9840 and 0.03932 cells.
The reduction factor is 25.0.  In physical coordinates, the no-fallback
distortion grows from 8.95e-4 on the coarse grid to 2.46e-3 on the fine grid.
The fallback result decreases from 1.23e-4 to 9.83e-5.  This opposite
refinement trend is the main stability result.

All runs finish normally.  No NaN, abort, or non-admissible state is found.
The final minimum pressure remains above 0.9986 and the final minimum density
remains above 0.9998.  The seeded mode is an alternating row-wise displacement.
Without the fallback, the maximum row-to-row pressure range increases from
10.22 on the coarse grid to 51.15 on the fine grid.  The corresponding
streamwise-velocity range increases from 0.988 to 7.560.  With the fallback,
the pressure ranges are 2.409 and 2.168, and the streamwise-velocity ranges
are 0.536 and 0.199.  The transverse velocity remains below 4e-17 in every
seeded run.  The result should therefore be described as odd-even row
decoupling, not as transverse-velocity growth.

## Figure caption

**AFD HLLC response to an odd-even perturbation of a stationary Mach 10
normal shock.**  The grid, CFL, time integrator, initial condition, and
executable are identical within each pair.  Only the shock-region LLF
fallback is changed.  Panel (a) shows the odd-even component of the shock
location.  Panel (b) shows the final front on the 400 by 40 grid.  The
unperturbed 200 by 20 controls retain zero odd-even amplitude for both
settings.  Without the fallback, the seeded displacement grows with time and
becomes stronger under refinement.  The LLF fallback confines the
distortion to less than 0.04 cells at the final time.

The figure is archived as `shock_stability_ab.pdf` and
`shock_stability_ab.png` under the result directory.
