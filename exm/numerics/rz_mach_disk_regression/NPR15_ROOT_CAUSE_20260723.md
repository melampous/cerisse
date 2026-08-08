# NPR15 Mach-disk root-cause audit, 2026-07-23

## Conclusion

The apparent 20-cell near-axis Mach-disk cliff is an analysis error, not a
20-cell deformation of the primary Mach disk.

The old tracker selected the largest positive density gradient independently
in each radial row. In the first three rows it selected a downstream
recirculation/stagnation compression. From the fourth row onward it selected
the actual Mach disk. Connecting those different features produced the
vertical cliff. Extending the comparison to `r/D=0.3--0.5` additionally mixed
the Mach disk with the barrel shock.

The field still contains a strong post-disk recirculation structure. That
feature is real in the computed R-Z solution and explains the visually rounded
cap, but it is not Mach-disk deformation.

## Exact identity evidence

The H-off restart from `chk05000` reproduces historical `plt05098` at
`t=2.250413726331 ms`. Relative L2 differences are `1.97e-13` for density,
`2.63e-13` for pressure, `4.58e-13` for radial velocity, and `8.86e-14` for
axial velocity.

At the axis:

| feature | z/D | Mach before/after | pressure before/after |
|---|---:|---:|---:|
| primary Mach disk | 2.39648 | 3.984 / 0.847 | 7.86 / 85.9 kPa |
| post-disk compression | 2.47557 | 0.448 / 0.244 | 100.4 / 236.7 kPa |

The downstream density-gradient peak is 1.111 times the primary peak in the
first ring. Across it, axial velocity changes from about `+145 m/s` to
`-89 m/s`, confirming that it bounds a subsonic reversal/recirculation region.

## Corrected geometry

The physical tracker selects the first compression with upstream Mach at
least 1.5 and downstream Mach at most 1.5.

| primary-disk metric | baseline | H=0.125 | wb1 | wb2 |
|---|---:|---:|---:|---:|
| core range, `r/D<=0.08` [L4 cells] | 0.7676 | 0.7282 | 0.7618 | 0.7324 |
| first-8 range [L4 cells] | 0.2556 | 0.1999 | 0.2500 | 0.2226 |
| first-16 range [L4 cells] | 0.5551 | 0.5147 | 0.5494 | 0.5201 |
| secondary/primary axis gradient | 1.1109 | 1.1274 | 1.1107 | 1.1158 |

The H and legacy well-balanced probes do not materially change the primary
disk. H slightly increases the downstream-to-primary gradient ratio.

Across the 30-frame statistics window:

- median instantaneous core range: `1.45` L4 cells;
- time-averaged instantaneous-locus range: `1.95` cells;
- online mean-density-locus range: `2.83` cells.

These are smooth small curvatures. They do not support the former 32.8-cell
mean-axis defect claim.

## Solver decision

The case uses characteristic LLF-WENO-Z5, not HLLC.

- Keep `cns.rz_euler_h_correction=0`.
- Keep `cns.rz_wb_axis_ncell=0`.
- Do not port `wb1`, `wb2`, or pressure quadrature into the production case as
  a response to this visual feature.
- Do not use row-wise maximum density gradient or a center-to-`r/D=0.3--0.5`
  line as a Mach-disk flatness metric.

No IBM production path is involved in this non-IBM underexpanded-jet audit.
No cut-cell, cut-control, aggregate, or EB flow mechanism was introduced.

## Remaining question

The post-disk reversal can be amplified by the enforced 2-D axisymmetric
`m=0` dynamics: coherent vortex rings cannot break through azimuthal modes.
The current tests rule out the proposed H and simple first-ring source fixes,
but do not by themselves prove whether the recirculation amplitude is a
physical laminar R-Z mode, a grid-dependent R-Z artifact, or mainly a
suppressed-3-D-mode limitation.

The discriminating sequence is:

1. repeat the same physical diagnostic at `D/128`, `D/256`, and `D/512`;
2. compare reversal length, minimum axial velocity, secondary/primary
   gradient ratio, and primary-disk locus in physical units;
3. run a full-azimuth 3-D case with small `m=1` and `m=2` perturbations;
4. compare time means and spectra after the transient, without using a narrow
   periodic sector that excludes low azimuthal modes.

## Artifacts

- Peak identity:
  [`output/npr15_peak_identity_r1/npr15_peak_identity.png`](output/npr15_peak_identity_r1/npr15_peak_identity.png)
- Exact-time physical A/B:
  [`output/npr15_physical_disk_ab_r1/npr15_physical_disk_tracks.png`](output/npr15_physical_disk_ab_r1/npr15_physical_disk_tracks.png)
- Corrected transient and mean history:
  [`output/npr15_primary_history_r1/npr15_primary_history.png`](output/npr15_primary_history_r1/npr15_primary_history.png)
- Well-balanced causal probe:
  [`output/npr15_wb_probe_r1/analysis/npr15_wb_probe_fronts.png`](output/npr15_wb_probe_r1/analysis/npr15_wb_probe_fronts.png)

The legacy max-gradient plots under `output/npr15_defect_ab_r2` are retained
only to document how the misidentification occurred.
