# Section 3.5 validation summary (accumulated statistics)

Conventions: sim uncertainty = block SE (3-8 tau-weighted blocks); exp
position interval = +-0.04 D (sampling resolution, not a confidence
interval); spacing interval = +-0.08 D (conservative propagation); Mj-mismatch floor
~1.4% (beta scaling); verdict rule: |dev| < interval -> within;
~interval -> marginal; > interval -> systematic.

## A. Mean shock-cell positions (3D L5 unless noted; exp = phase-mean M1.42)
| Quantity | Sim | Exp | Dev | Interval | Verdict |
|---|---|---|---|---|---|
| x_fc | 0.883 | 0.888 | -0.005 | 0.04 | within |
| x_1 | 1.178 +- 0.007 | 1.199 | -0.021 | 0.04 | within |
| x_1 (2D L6) | 1.105 +- 0.018 | 1.199 | -0.094 | 0.04 | systematic (upstream) |
| S_1 | 1.301 | 1.359 | -0.058 | 0.08 | within |
| S_2 | 1.105 | 1.186 | -0.081 | 0.08 | marginal (short) |
| S_3 | 1.004 | 1.130 | -0.126 | 0.08 | systematic (short) |

## B. Centreline density extrema (fig4 M1.43; see table3p8_truestat.md)
3D L5 amplitudes within 1.2-2.4% of experiment at all six extrema;
positions within 0.08 D. 2D L6: first cell comparable, second/third
maxima low by 22-25% (downstream collapse).

## C. Radial profiles (fig5a M1.44), RMS difference in rho/rho_j
| x/D | 3D L5 | 2D L6 |
|---|---|---|
| 0.60 | 0.031 | 0.032 |
| 0.75 | 0.034 | 0.026 |
| 0.90 | 0.048 | 0.049 |
| 1.05 | 0.035 | 0.098 |
| 1.25 | 0.057 | 0.072 |
| 1.55 | 0.130 | 0.151 |
| 2.00 | 0.046 | 0.107 |
| 3.05 | 0.094 | 0.072 |
First cell (<= 1.25): both formulations 0.03-0.10; post-closure station
1.55 worst for both (axial displacement of the recovery); downstream: 3D
retains 0.05-0.09, 2D 0.07-0.15.

## D. Lip line r/D = 0.63 (fig10bc time-mean), MAE in rho/rho_j over x 0.6-6.9
3D L5: 0.075; 2D L6: 0.137.

## E. Multi-Mach first-cell metrics (2D sweep 99 um vs fig4 five traverses)
| Mj | x_1 sim/exp | S_1 sim/exp |
|---|---|---|
| 1.19 | 0.71 / 0.73 | 0.73 / 0.74 |
| 1.42/1.43 | 1.20 / 1.21 | 1.29 / 1.38 |
| 1.59/1.60 | 1.55 / 1.79 | 1.93 / 1.54 |
| 1.78/1.80 | 1.82 / 1.79 | 2.95 / 2.10 |
x_1 tracks experiment below Mach-disk formation and at 1.8; S_1 follows
Pack and experiment for Mj <= 1.5 and degrades once a Mach disk forms
(post-disk axisymmetric limitation, cf. dimensionality section).

## F. Convective Mach number (fig10d): definitional mismatch
Broadband cross-correlation Uc = 0.52-0.88 c_inf lies along the lower
envelope of the screech-locked phase-speed data (mean 1.06); qualitative
comparison only (no screech feedback in the computation).

## G. Not comparable (excluded with reason)
Screech-band rms density (fig10bc rms; narrowband, no screech in sim);
near-field acoustics (fig10a SPL, m142prs); phase-averaged mode shapes.
