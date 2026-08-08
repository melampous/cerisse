# Test-1853 repair endurance and H100 grid results (2026-08-08)

## Test D: support-hook repair endurance

The repaired jet-on case advanced from chk03400 to Step 4000 at CFL 0.01 with
normal dynamic regridding and finished `complete_pass`:

- exit status 0, 600 coarse steps, final time 0.7056776 ms;
- zero state-check, inadmissible-state, soft-positivity, NaN, or Inf markers;
- final L2 GP and surface targets both zero;
- all recorded coarse-fine missing-target/support counts zero;
- pure shared-GP full-Cartesian manifest, cut-control disabled;
- 15,869,440 final valid cells and 98,217,984 SSPRK43 cell-stages per coarse step.

The final 200-step timing was 6.27 s for ordinary steps and 9.33 s for steps
that rebuilt AMR/IBM data.  With `regrid_int=4 4`, an `lbase=1` rebuild occurred
every two coarse steps.  This validates the support-hook repair under many
regrids, but CFL 0.01 is not a flux/RK robustness result and this field is not a
production solution.

The full AMR leaf-fluid audit gives:

| Step | Tmin (K) | T0min (K) | Total energy (J) |
| ---: | ---: | ---: | ---: |
| 3400 | 21.938 | 64.730 | 1016.919 |
| 3600 | 21.575 | 64.730 | 1021.422 |
| 3800 | 21.903 | 64.730 | 1025.969 |
| 4000 | 20.187 | 64.730 | 1030.551 |

No leaf-fluid cell was below 20 K.  The lowest final static temperature is an
accelerated level-1 fluid cell at `(17.81,-56.18,-95.45)` mm, not the earlier
failure point.  Total mass/energy are not closed-domain conservation residuals
because Run165 uses open boundaries.

The final image also shows that the field remains a staged startup state: the
downstream part of the long domain is still quiescent.  This is acceptable for
the repair endurance objective, not for physical SRP statistics.

## New H100 no-jet L0 field

The successful cold-start run is:

```text
/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/
run165_h100_largebox_nojet_20260808/stage0c_external_l0
```

It used `256x192x192` isotropic 2.25 mm cells over
`[-0.192,0.384] x [-0.216,0.216]^2` m, eight H100 ranks, SFC, blocking factor
32, maximum box size 128, no jet, and L0 only.  AMReX produced exactly eight
boxes, from `128x64x64` through `128x128x128`.

The run reached 0.5 ms at Step 1339 and finished `complete_pass`, with no state,
thermodynamic, or support errors.  HBM usage was approximately 61.4-62.5 GiB
per GPU.  The all-step median was 0.448 s and mean was 0.485 s; final ordinary
steps were about 0.54-0.56 s.  The previous fragmented `256x176x176` no-jet L0
case was about 1.06 s/step, so the large-box layout is roughly twice as fast
despite 19% more base cells.

The final Tmin is 64.634 K and no cell is below 20 K.  The measured axis bow
shock is at `x=-24.375 mm`, a 26.816 mm standoff from the STL nose.  At
`y/z=+/-30 mm`, all four detected locations are `x=-17.625 mm`; the y-pair,
z-pair, and mean-y/mean-z differences are zero at this 2.25 mm resolution.

The near-body bow shock is established and symmetric, but the final image shows
that the downstream startup front has not left the domain.  The checkpoint is
therefore a good initial condition for a no-jet L1 equilibration stage, not a
globally steady external-flow solution.

Two prior directories stopped intentionally before Step 1 and contain no flow
result.  `stage0_external_l0` caught a jet-on executable used with jet-off input;
`stage0b_external_l0` caught the invalid historical `stage_positivity=1` option
without the 3-D shared-face limiter.  Both fail-closed checks behaved correctly.

## Next gate

Restart `stage0c_external_l0/output/run165/chk01339` with jet still disabled and
`max_level=1`.  First perform one dynamic hierarchy build with strict GP/surface
support audits and measure boxes, cells, and HBM against the remaining roughly
18 GiB per GPU.  Continue long enough for the downstream startup front to leave
the region of interest, then export and freeze the accepted hierarchy.  Only
after that should the jet-on executable add the pure-fluid L2 nozzle core and
run the CFL 0.05/0.10/0.20/0.30 matrix.

No L1 or jet-on continuation was submitted as part of this result.
