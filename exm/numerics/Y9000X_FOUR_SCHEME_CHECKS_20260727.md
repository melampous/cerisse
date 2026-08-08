# Y9000X four-scheme numerical checks

## Scope and configuration

The checks compare the following four Euler discretisations.

- LLF with WENO-Z5 reconstruction
- LLF with TENO5 reconstruction
- AFD-HLLC with WENO-Z5 reconstruction
- AFD-HLLC with TENO5 reconstruction

All calculations use double precision, a CFL number of 0.3, and SSPRK(3,3).
They contain no viscous flux, turbulence model, combustion model, IBM, or
embedded-boundary flow update. The AFD calculations use
`cns.afd_correction=1`, `cns.afd_shock_llf=1`,
`cns.afd_smoothness_threshold=0.08`, and
`cns.afd_teno_cutoff=1.0e-4`. They are therefore AFD-HLLC calculations with
the configured LLF fallback near detected non-smooth regions.

The calculations were performed on the Y9000X host `Melampous`, which has an
Intel Core i9-13900HX processor. The small one-dimensional Shu-Osher cases use
one process. The two-dimensional vortex cases use four MPI processes.

The LLF-TENO5 implementation uses its code default cutoff of `1.0e-3`.
Consequently, these results compare the four implemented methods. They do not
isolate the effect of the Riemann solver while holding every nonlinear
selection parameter fixed.

## Shu-Osher shock and entropy-wave interaction

All four short smoke tests completed successfully. The formal calculations
use the canonical domain from -5 to 5 and advance to \(t=1.8\). The table
reports differences relative to the \(N=800\) result from the same method.
These values measure resolution sensitivity. They are not errors relative to
an exact or converged solution.

| Method | \(L_1\) difference at \(N=200\) | \(L_1\) difference at \(N=400\) | Density TV at \(N=800\) |
|---|---:|---:|---:|
| LLF WENO-Z5 | \(5.63781\times10^{-2}\) | \(1.08601\times10^{-2}\) | 20.1686 |
| LLF TENO5 | \(5.31165\times10^{-2}\) | \(1.10423\times10^{-2}\) | 20.0208 |
| AFD-HLLC WENO-Z5 | \(5.11502\times10^{-2}\) | \(1.03139\times10^{-2}\) | 20.2086 |
| AFD-HLLC TENO5 | \(4.52613\times10^{-2}\) | \(1.12470\times10^{-2}\) | 20.2644 |

The four \(N=800\) profiles are closely aligned. All methods preserve the
physical minimum density of approximately 0.8. No calculation reports a
state-check failure, NaN, or abort. The AFD-HLLC TENO5 result retains more of
the short-wave response on the \(N=200\) grid. This distinction is much
smaller on the two finer grids.

## Uniform-grid isentropic vortex

The fast vortex uses \(M=0.5\) and \(\beta=0.2\) in a periodic square domain.
It is advanced through one flow-through time. The density error is evaluated
against the translated analytic solution.

| Method | \(L_2(\rho)\), \(40^2\) | \(L_2(\rho)\), \(80^2\) | \(L_2(\rho)\), \(160^2\) | Rates | \(K^\prime/K^\prime_{\mathrm{exact}}\), \(160^2\) |
|---|---:|---:|---:|---:|---:|
| LLF WENO-Z5 | \(1.49186\times10^{-4}\) | \(1.49616\times10^{-5}\) | \(1.77416\times10^{-6}\) | 3.32, 3.08 | 0.999164 |
| LLF TENO5 | \(1.33442\times10^{-4}\) | \(1.32384\times10^{-5}\) | \(1.30327\times10^{-6}\) | 3.33, 3.34 | 0.999235 |
| AFD-HLLC WENO-Z5 | \(1.97341\times10^{-4}\) | \(1.09428\times10^{-5}\) | \(1.14049\times10^{-7}\) | 4.17, 6.58 | 0.999840 |
| AFD-HLLC TENO5 | \(9.31184\times10^{-5}\) | \(3.60448\times10^{-6}\) | \(1.23437\times10^{-7}\) | 4.69, 4.87 | 0.999842 |

Every sequence shows a large and monotonic reduction in analytic error. The
LLF sequences approach the third-order limit expected from SSPRK(3,3) when
the time step is proportional to the grid spacing. The larger apparent AFD
rates over these three grids are pre-asymptotic combined rates. They must not
be reported as the formal spatial order. A separate time-step study would be
required to isolate the spatial order.

## Isentropic vortex crossing AMR interfaces

The AMR hierarchy has one fixed refined level over the right half of the
domain. The refinement ratio is two. Refluxing and level subcycling are
enabled. The vortex starts in the coarse half and crosses both coarse-fine
interfaces during one flow-through.

| Method | \(L_2(\rho)\), \(40^2\) base | \(L_2(\rho)\), \(80^2\) base | \(L_2(\rho)\), \(160^2\) base | Rates | \(K^\prime/K^\prime_{\mathrm{exact}}\), \(160^2\) base |
|---|---:|---:|---:|---:|---:|
| LLF WENO-Z5 | \(1.85373\times10^{-4}\) | \(5.63040\times10^{-5}\) | \(1.70574\times10^{-5}\) | 1.72, 1.72 | 0.985096 |
| LLF TENO5 | \(1.75630\times10^{-4}\) | \(5.11475\times10^{-5}\) | \(1.45772\times10^{-5}\) | 1.78, 1.81 | 0.988282 |
| AFD-HLLC WENO-Z5 | \(2.40352\times10^{-4}\) | \(7.86879\times10^{-5}\) | \(2.34116\times10^{-5}\) | 1.61, 1.75 | 0.986378 |
| AFD-HLLC TENO5 | \(1.97135\times10^{-4}\) | \(6.89966\times10^{-5}\) | \(2.10853\times10^{-5}\) | 1.51, 1.71 | 0.990318 |

All AMR calculations remain stable and conservative to round-off. Across the
complete matrix, the maximum absolute relative mass and energy changes are
\(1.382\times10^{-13}\) and \(1.578\times10^{-13}\), respectively. The AMR
density error at a base resolution of 160 is 9.61 to 205 times the
same-method uniform-grid error. The measured AMR rates are between 1.51 and
1.81. This behaviour is consistent with the current linear coarse-fine
interpolation dominating the error after interface crossings. The AMR case
therefore qualifies stability, conservation, and interface sensitivity. It
does not verify the formal order of the interior WENO or TENO reconstruction.

## Thesis-use decision

The uniform isentropic-vortex sequence is suitable as an analytic accuracy
check. The Shu-Osher sequence is suitable as a qualitative shock-interaction
and resolution-sensitivity comparison. The AMR sequence should be presented
as an interface test and as a documented limitation of the current
coarse-fine treatment. It should not support a claim of fifth-order accuracy
across AMR interfaces.

The recorded advance times are single-run diagnostics. They are not
performance benchmarks. A timing comparison would require repeated runs on
an otherwise idle machine and a statistic such as the median.

## Result archives

- Shu-Osher:
  `shuosher/results/y9000x_four_scheme_20260727_r2`
- Isentropic vortex:
  `isentropic_vortex/results/y9000x_four_scheme_20260727_r2`

The clean result matrix contains 16 successful Shu-Osher run-status records,
12 Shu-Osher formal plotfiles, 24 finalized vortex logs, 24 vortex plotfiles,
and 24 rows in the combined vortex summary.
