# Section 3.3 revision record: Panda-condition cases

Date recorded: 2026-07-17
Purpose: preserve the supplied conversation and turn it into a traceable evidence ledger for the later revision of thesis Section 3.3. This file is not replacement thesis prose and does not by itself promote a reported result to a validated conclusion.

## 1. Preserved source and evidence hierarchy

The complete supplied conversation is preserved byte-for-byte at:

`analysis/panda_benchmark/conversation_archive/2026-07-17_case_configurations_review_and_cd_analysis.txt`

Archive integrity:

- 676 lines, 3,526 words, 41,557 bytes;
- SHA-256: `a7e853cb5a1a39cdf459b8c56bfd504b7d949f349687777b6dcc9014d62625cf`;
- the archived copy was verified with `cmp` against the original Codex attachment.

The principal configuration source is:

`analysis/panda_benchmark/case_configurations.md`

It is an as-built configuration record derived from the four cases' `inputs` and `prob.h`; it is not a results record. Its SHA-256 at the time of this note was `6f54eba34f45fc4a1d129d0c51f91592a7125b2a3482d44e8e72d1dc22669458`.

Status labels used below:

- **CONFIRMED CONFIGURATION**: present in the current as-built configuration record and/or explicitly checked against input definitions in the supplied conversation;
- **CONFIRMED CORRECTION**: a counting or interpretation error explicitly checked in the conversation;
- **REPORTED RESULT**: reported in the conversation and retained for later checking against the underlying arrays, plots or logs;
- **HYPOTHESIS**: a proposed physical or numerical explanation, not a conclusion;
- **PENDING**: requires another calculation, source check or controlled comparison.

If a thesis value conflicts with this ledger, the final as-run input, executable banner, probe header, plotfile header and analysis script take precedence, in that order as applicable. Chat summaries never override those primary records.

## 2. Configuration facts to carry into the Section 3.3 rewrite

### 2.1 Model designation and experimental differences

**CONFIRMED CONFIGURATION**: use the designation

> Panda-condition, prescribed-exit, simplified flush-baffle model

Do not call the calculations `Panda-exact`. The documented differences are:

1. the internal converging-nozzle flow is not resolved;
2. pressure, temperature and velocity are prescribed directly on a circular sonic exit plane;
3. the physical nozzle lip is not represented;
4. the imposed tanh layer, with nominal `delta/D = 0.01`, is a modelling assumption rather than a Panda measurement;
5. the effective discharge coefficient is a property of the imposed profile/discretisation, not an experimental value;
6. the wall is an ideal coplanar adiabatic slip baffle, whereas Panda used a 305 mm metal nozzle block/flange that was important to the screech feedback;
7. the approximately 8 m/s flow in the Rayleigh-scattering density facility was transverse apparatus flow, not axial coflow; the baseline calculation omits it and uses a quiescent ambient.

### 2.2 Shared gas, ambient and exit definition

**CONFIRMED CONFIGURATION**:

- calorically perfect air, `gamma = 1.4`, with Sutherland viscosity and conductivity;
- nominal `R = 287.06 J/(kg K)` in the configuration record;
- `p_amb = 99,780 Pa`, `T_amb = T0 = 297.15 K`, and zero ambient velocity;
- circular choked exit, `M_e = 1`, `D = 25.4 mm`, `R_exit = 12.7 mm`;
- `T_e = 247.625 K` and `U_e approximately 315.5 m/s`;
- inside the aperture, `p = p_e` over the full radius, `u = U_e f(r)`, `T = T0 - u^2/(2 c_p)`, and `rho = p_e/(R T)`;
- `f(r) = 0.5[1 - tanh((r-r0)/delta)]`, with `delta = 254 micrometres` and `r0 = R_exit - 3 delta`, giving nearly zero velocity at the physical lip;
- from the physical lip to radius `6D`, the boundary is an adiabatic slip wall with mirrored ghost parity; outside `6D`, the x-low/z-low corner is prescribed quiescent ambient.

The three pressure-ratio cases are:

| case | NPR | `p_e/p_amb` | `p_e` [Pa] | `M_j` | `T_j` [K] | `U_j` [m/s] | `rho_j` [kg/m3] | nominal `Re_j` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| m119 | 2.393613 | 1.2645 | 126,173 | 1.19 | 231.6 | 363.0 | 1.5011 | 0.92e6 |
| m142 | 3.273446 | 1.7293 | 172,550 | 1.42 | 211.75 | 414.2 | 1.6413 | 1.24e6 |
| m180 | 5.745796 | 3.0354 | 302,871 | 1.80 | 180.3 | 484.5 | 1.9276 | 1.95e6 |

The Reynolds numbers are nominal fully expanded reference values. They must not be described as proof that the imposed exit boundary layer, mass flow or momentum flux matches the experiment.

### 2.3 Shared numerical method

**CONFIRMED CONFIGURATION**:

- compressible Navier--Stokes equations;
- characteristic-wise WENO-Z5 flux-vector splitting with face-local Lax--Friedrichs dissipation and an internal first-order positivity fallback;
- second-order centred viscous fluxes;
- SSP-RK3 with `CFL = 0.30`;
- AMReX refinement ratio 2 and `do_reflux = 1`;
- dt-weighted online first and second moments.

**CONFIRMED CORRECTION**: online statistics comprise 11 components in 2D and 15 components in 3D. Consequently, the output decomposition is 20 components in 2D (`5 conserved + 11 statistics + 4 derived`) and 25 in 3D (`5 conserved + 15 statistics + 5 derived`).

### 2.4 Axisymmetric cases

**CONFIRMED CONFIGURATION**:

- m119, m142 and m180 use the same 2D RZ `inputs` byte-for-byte;
- their `prob.h` files differ only in NPR, except that m180 gates the density-gradient sensor with `nt > 0` during initial hierarchy construction;
- domain: `r in [0,12D]`, `z in [0,32D]`;
- base grid: `192 x 512 = D/16`;
- `max_level = 6`, giving nominal finest spacing `D/1024`;
- density-gradient thresholds 0.05/0.03, level-dependent spatial caps, forced L4 core, L5 shear-layer corridor and L6 lip corridor;
- 32 active pressure probes, all sampled at level 0 every coarse step;
- statistics windows: m119 10.0 ms, m142 9.25 ms and m180 3.5 ms.

The three cases form a tightly controlled NPR series in their governing configuration, but their statistics windows and the m180 initial sensor gate are not identical. The thesis should state both facts.

**CONFIRMED CORRECTION**: the 2D probe total is 32, not 31.

**Grid-history qualification**: the final `D/1024` hierarchy was not active from `t = 0`; the runs developed from coarser seeds and were subsequently uplifted. Section 3.3 must distinguish nominal final resolution from resolution history.

### 2.5 Three-dimensional m142 case

**CONFIRMED CONFIGURATION**:

- Cartesian domain with jet axis x: `x in [0,24D]`, `y,z in [-8D,8D]`;
- base grid `192 x 128 x 128 = D/8`;
- static production BoxArray with levels L1--L5 and about 50.9 million cells;
- nominal finest spacing `D/256`, confined to the lip/shear-layer sleeve;
- the near-axis shock-cell column has no L5 coverage and is principally at L4, i.e. `D/128`;
- the fourth shock peak lies near an L4-to-L3 interface;
- the final static L5 hierarchy was introduced at approximately 5.0 ms and ran to 10.16 ms, so it was active for about 5.16 ms rather than from `t = 0`;
- a one-off broadband non-axisymmetric seed with `m = 1...5`, total azimuthal RMS amplitude `3e-4 U_e`, was applied over 5.0--5.2 ms;
- statistics were accumulated over 6.5--10.16 ms, a 3.66 ms window;
- 216 probes are defined, but only 210 are active: 12 rings x 16 points = 192 ring probes plus 18 single probes;
- q1--q5 are near-lip-line axial stations, not centreline probes; a1--a3 are the centreline group;
- six legacy probes are defined but inactive.

**Boundary qualification**: the open `foextrap` faces are not mathematically non-reflecting. The supplied estimate gives a transverse acoustic round trip of about 1.18 ms across the `+/-8D` domain, so a disturbance could return several times during the statistics window. This does not invalidate mean-flow comparisons, but it precludes unqualified quantitative screech amplitude or phase claims.

## 3. Corrections established in the supplied review

The following corrections were explicitly checked and should be treated as settled unless a primary as-run file later disproves them:

1. 2D/3D online-statistics counts are 11/15, not a common 11.
2. The 3D plotfile split is `5 + 15 + 5 = 25`, not `5 + 14 + derived`.
3. The 2D probe count is 32.
4. The 3D active probe count is 210, although 216 definitions exist.
5. There are 12 active azimuthal rings, not 13.
6. q1--q5 are near-lip-line probes; a1--a3 are on the axis.
7. Agreement of nominal fully expanded `Re_j` with Panda is not agreement of the actual imposed exit flux or boundary layer.
8. The final fine hierarchy was introduced after a coarser development phase in both 2D and 3D.
9. The model must be described as prescribed-exit and simplified-baffle, not Panda-exact.
10. `foextrap` must not be described as a non-reflecting acoustic boundary.

The conversation also records that a matched-window 2D/3D comparison had already been made using approximately 3.55 ms of 2D statistics and 3.66 ms of 3D statistics. This is a **REPORTED RESULT** and should be checked against the relevant plotfile and analysis metadata before citation.

## 4. Result interpretations preserved for later results-section review

These statements belong primarily in the validation/results sections, not in the Section 3.3 method description.

### m119 axisymmetric calculation

**REPORTED RESULT**:

- the 10 ms statistics/record window spans roughly 84 periods at 8.4 kHz;
- no strong, statistically persistent Panda-like 8.4 kHz axisymmetric tone was found;
- reported near- and outer-field peaks around 3.5 and 2.9 kHz do not establish a common propagating tone.

Permitted wording after data recheck:

> The present prescribed-exit, simplified-baffle RZ model did not reproduce the 8.4 kHz axisymmetric screech reported by Panda.

Do not infer that the physical Panda operating condition itself lacks screech.

### m142 axisymmetric calculation

**CONFIRMED INTERPRETIVE LIMIT**: RZ excludes the experimental helical `m = +/-1` mode by construction. The case can support mean shock-cell validation but cannot provide a negative acoustic validation when a 5.375/5.4 kHz helical tone is absent.

### m142 three-dimensional calculation

**REPORTED RESULT**:

- the first cell and first strong recompression, and the broad experimental mean-density shape, compare favourably;
- the 3D mean retains a longer shock train than 2D in the reported matched-window comparison;
- downstream peak locations progressively lead the experiment by about 0.25--0.30D;
- the 3.66 ms window contains about 20 cycles at 5.375 kHz and shows no strong, saturated, persistent `m = +/-1` screech.

The last statement does not exclude a weak or slowly growing global mode. Do not claim that three-dimensionality alone corrected the 2D result, that ILES dissipation is the established cause of shortened cells, or that experimental centreline density alone uniquely identifies a Mach disk.

### m180 axisymmetric calculation

**REPORTED RESULT**: a strongly underexpanded barrel-shock/Mach-disk-type structure and downstream shock train are present. The 3.5 ms statistics window and current comparison data do not yet support quantitative Mach-disk position/diameter validation, complete Zapryagaev validation or a screech-presence/absence conclusion.

## 5. Effective-diameter/discharge-coefficient hypothesis

### 5.1 Motivation

**HYPOTHESIS**: because `r0 = R - 3 delta = 0.47D`, the velocity-profile diameter is approximately `0.94D`, with area ratio `0.94^2 = 0.8836`, close to the reported effective discharge coefficient. If shock-cell length scaled primarily with effective jet diameter, normalisation by nominal `D` could account for roughly 6% of the reported approximately 7% shortened spacing.

This correspondence is suggestive, not proof. In particular, the first cell reportedly agrees well while the phase error accumulates downstream.

### 5.2 Legacy-profile screen

The conversation records the following profile-integrated values:

| legacy profile | `C_d` | `sqrt(C_d)` | reported first counted shock-cell peak `x/D` |
|---|---:|---:|---:|
| top-hat | 1.0000 | 1.0000 | 1.34 |
| tanh 200 micrometres at R | 0.9666 | 0.9832 | 1.24 |
| Pohlhausen/Blasius | 0.9444 | 0.9718 | 1.14 |
| tanh 100 micrometres at R | 0.9831 | 0.9915 | 1.18 |
| current v2 profile, 254 micrometres at `R-3 delta` | 0.8789 | 0.9375 | not part of this legacy family |

The peak detector also returned a near-exit feature at `x/D approximately 0.34`; the conversation's final comparison treated the next peak as the first shock-cell peak. That selection rule must be documented if these numbers are reused.

**REPORTED VERDICT**: the legacy family does not isolate a `sqrt(C_d)` scaling. Its observed 7--15% changes are much larger than the 1--3% change predicted from `sqrt(C_d)`, are not monotonic in `C_d`, and are confounded by a different 10 mm nozzle, a 6D x 24D unflanged domain, different profile centring, coarser `D/256` resolution and short startup-tail averaging windows. The useful conclusion is only that inlet-profile details are a first-order sensitivity.

### 5.3 Controlled discriminator

**REPORTED/PENDING**: an AWS case named `m142_cd0966` was prepared but not submitted. It reportedly changes only the tanh centre from `R-3 delta` to `R`, holding `delta = 254 micrometres`, to raise `C_d` from approximately 0.879 to approximately 0.966. This state was not independently rechecked in this recording turn.

Interpretation planned in the conversation:

- a roughly 4.7% systematic cell-length increase would support effective-diameter scaling;
- if the already well-placed first cell overshoots the experiment, a pure global `C_d` explanation is weakened;
- if the downstream phase-decay pattern is unchanged after only an overall rescaling, `C_d` is likely a leading contributor;
- otherwise inlet momentum thickness/profile shape, grid-interface placement and downstream numerical/acoustic effects remain entangled.

The hypothesis and pending run do not belong in Section 3.3 as established findings. Section 3.3 should define the inlet and its uncertainty; the test and outcome belong in a sensitivity-results subsection.

## 6. Direct map for rewriting Section 3.3

Section 3.3 should contain, in this order:

1. **Scope and model name**: Panda-condition, prescribed-exit, simplified flush-baffle; state that the nozzle interior is not resolved.
2. **Reference conditions and case table**: `D`, ambient state, sonic exit state, NPR, `M_j`, fully expanded reference quantities and the distinction between nominal `Re_j` and actual imposed flux.
3. **Inlet and x-low boundary**: give the tanh formula, `delta`, `r0`, constant exit pressure, constant-total-temperature relation, physical-lip transition to the 6D slip baffle and ambient corner.
4. **Numerical method**: WENO-Z5/local LLF, second-order viscosity, SSP-RK3, CFL, reflux and statistics definition.
5. **2D matrix and mesh**: controlled NPR series, RZ limitation, final nominal resolution, caps/forced zones and actual resolution history.
6. **3D mesh and perturbation**: static BoxArray, where each resolution actually occurs, lack of L5 on the axis, seed definition and final-grid history.
7. **Diagnostics and statistics windows**: correct active probe counts, ring count, probe roles, output components and exact uninterrupted windows.
8. **Known differences and claim boundary**: prescribed exit, unmeasured layer, inferred/discrete `C_d`, simplified flange, omitted transverse facility flow, RZ modal restriction and non-NRBC open boundaries.

Do not place job IDs, GPU throughput, queue history, speculative explanations for result discrepancies or unfinished tests in the main Section 3.3 narrative. Those belong in provenance records, sensitivity results or an appendix.

## 7. Conflicts that must be resolved before thesis insertion

1. **3D lateral extent wording**: the as-built table gives `y,z = +/-8D`, but its note begins “lateral +/-6D equals flange radius”. The likely intended distinction is domain boundary `+/-8D` versus baffle radius `6D`; final `inputs` must be checked and the sentence repaired before reuse.
2. **Discharge coefficient**: `case_configurations.md` reports a common `C_d,eff = 0.8789`, whereas `U_jet/analysis/CASE_SPEC_panda2d_v2_2026-07-13.md` reports banner values about 0.874 for m142 and 0.877 for m119. The thesis should use per-case discrete as-run values from the final banners/flux audit and label 0.8789 as an analytic/profile estimate if that is what it represents.
3. **Gas constant**: the configuration record uses `R = 287.06`, while earlier case notes mention a compiled-code/banner value near 287.10. Use the value from the final executable/source and recompute dependent rounded values consistently.
4. **Actual probe coordinates**: named coordinates are planning labels; thesis tables should use box/index-centre coordinates parsed from final headers.
5. **Conservation wording**: `do_reflux = 1` is a configuration fact, not proof that the composite-grid mass/energy ledger closes.
6. **Acoustic boundary statement**: the 1.18 ms estimate is consistent with a `+/-8D` transverse boundary. It must not be paired with an erroneous claim that the computational side boundary is at `+/-6D`.
7. **Results provenance**: every numerical discrepancy, spectral conclusion and matched-window statement listed in Section 4 must be tied to the exact plotfile/probe range and analysis script before entering the thesis.

## 8. Preserved pending work

The conversation's pending queue is:

1. complete the controlled `C_d`/inlet-profile discriminator;
2. overlay AMR-level boundaries on 3D figures, especially near `x approximately 0.75D`, `4D` and the fourth peak;
3. compute ring-based `m-f` spectra and `m = +/-1` envelopes;
4. compute `k_x-f` diagnostics for downstream KH waves and upstream-propagating waves;
5. output counts and spatial locations of the first-order positivity fallback;
6. create a dedicated `VALIDATION_RESULTS.md` linking every claim to statistics windows, raw unshifted errors, optionally aligned shape errors, spectral stability and figure/data paths;
7. retain TENO conclusions as pending: the available 2D TENO window is short and the 3D TENO branch had not been run in the supplied conversation.

## 9. File index for the next Section 3.3 pass

- complete conversation: `analysis/panda_benchmark/conversation_archive/2026-07-17_case_configurations_review_and_cd_analysis.txt`;
- as-built configuration: `analysis/panda_benchmark/case_configurations.md`;
- earlier 2D specification/evolution log: `U_jet/analysis/CASE_SPEC_panda2d_v2_2026-07-13.md`;
- current extracted fields and comparisons: `analysis/panda_benchmark/`;
- current data audit: `analysis/panda_benchmark/current_data_analysis/CURRENT_DATA_ANALYSIS.md`;
- screech requirements: `analysis/panda_benchmark/SCREECH_REPRODUCTION_REQUIREMENTS.md`;
- project handoff: `U_jet/analysis/HANDOFF_UJET_2026-07-10.md`;
- authoritative thesis worktree, when revision is explicitly requested: `/home/qiaoj/Documents/thesis-overleaf`.

No thesis file was modified while creating this record.
