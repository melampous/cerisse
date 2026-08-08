# AFD shock-LLF switch: advected-contact results

> Historical baseline for the broad pre-gate sensor. The authoritative
> joint-sensor rerun is
> `results/afd_joint_sensor_20260729_r2` and is summarised in
> `../AFD_JOINT_SHOCK_SENSOR_20260729.md`.

The authoritative archive is
`results/afd_shock_llf_ab_20260729_r4`.  All eight calculations and the
analysis completed successfully.  The earlier successful `r3` calculation
gave an identical `comparison.csv`.  The `r4` archive is retained because its
analysis source snapshot also includes the initial-sensor diagnostic.

Both switch values use the same AFD HLLC WENO-Z5 executable, uniform grid,
CFL 0.20, SSPRK(3,3), one-transit stop time and all other run-time options.
The only paired-run difference is `cns.afd_shock_llf`.

An automated command audit checked all four pairs.  After removing the
required output-path token and the switch token, every paired command is
identical.  All eight run status files report exit status zero.  The common
executable SHA256 is
`337c2d09d73363edbebc3127a66e2c49958159db66f3b9239d6503801b69c393`.

| Cells | Density L1, switch 0 | Density L1, switch 1 | Ratio 1/0 | Mixed cells, both | Initial sensor faces | Final sensor faces |
|---:|---:|---:|---:|---:|---:|---:|
| 64  | 4.28981e-2 | 4.27733e-2 | 0.99709 | 21 | 10 | 0 |
| 128 | 2.48282e-2 | 2.47931e-2 | 0.99858 | 26 | 10 | 0 |
| 256 | 1.45372e-2 | 1.45098e-2 | 0.99811 | 31 | 10 | 0 |
| 512 | 9.17460e-3 | 9.15442e-3 | 0.99780 | 35 | 10 | 0 |

The corresponding L1 rates for switch 0 are 0.789, 0.772 and 0.664.  The
rates for switch 1 are 0.787, 0.773 and 0.665.  These rates describe a
discontinuous solution and are not a formal-order measurement.

The current sensor flags ten faces in the exact initial state.  Every trigger
comes from density.  Pressure and velocity are constant.  The contact has
become smooth enough that the offline final snapshot has no trigger.  These
are snapshot counts rather than time-integrated flux counts.

The two switch values give the same one-percent mixed-cell count on every
grid.  Their density L1 errors differ by less than 0.3 per cent.  The pressure
L-infinity errors are also similar.  At 512 cells they are 2.366e-3 for
switch 0 and 2.429e-3 for switch 1.  Density remains positive.  Mass,
momentum and total-energy conservation errors remain below 4.5e-13.

This test does not show a material contact-resolution penalty from the switch
over one transit.  It does show that the present density-based criterion
classifies a pure contact as non-smooth even though no compression shock is
present.

Core source SHA256 values:

```text
Afd.h      5d461b134bed5e81f833b07b1b9fe6fb2c4267f79b37728801979bea8e817af9
AfdIBM.h   b1b3318561f19485847fa229432f25912166e2b9315b7ef00e2c71b5d09d1337
Riemann.h  d8a453fb9a886943453b6aae60368a5b6ea43ce7a7da172d987ae69593e3119d
```
