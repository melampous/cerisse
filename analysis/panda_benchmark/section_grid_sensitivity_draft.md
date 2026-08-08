# §3.X Grid sensitivity, dimensionality effects, and validation scope for the M_j = 1.42 underexpanded jet

<!-- DRAFT v2 — after 3-way audit (numeric vs data files / claim strength / terminology-completeness;
     33 confirmed findings applied). All numbers traceable to
     analysis/panda_benchmark/shockcell/shock_cell_metrics.csv, AXT_*.npz, AX3T_*.npz.
     Figure placeholders:
       [Fig. GRID-2D]  = shockcell/paperfig_grid2d_L1to6.png
       [Fig. GRID-3D]  = m142_3d_grid_levels.png
       [Fig. AXIS-FAM] = m142_2d_gridfam_axis.png (centreline mean-profile overlay across the family)
       [Fig. RMS-FAM]  = shockcell/fig37_axis_mean_rms_bands.png (mean ± rms bands per level)
       [Fig. EXP-CMP]  = comparison with the EPAPS centreline data (e.g. paperfig_expsim_rho_m142.png)
     Cross-reference placeholders: §3.Y = shock-cell spacing/upstream-shift analysis; §2.Z = numerical method. -->

## 3.X.1 Methodology and metric definitions

Two grid families are considered for the M_j = 1.42 case: an axisymmetric (RZ) family of
eight configurations with finest cell sizes from Δx = D_e/16 to D_e/1024, and a
three-dimensional family of five configurations from D_e/16 to D_e/256
([Fig. GRID-2D], [Fig. GRID-3D]). All lengths are normalized by the nozzle exit diameter
D_e = 25.4 mm (the convention used in all figures); pressures by the ambient static
pressure p∞; and densities by the fully expanded (isentropic) jet density ρ_j, the
normalization of the Panda & Seasholtz (1999) datasets. All computations are implicit
LES: the upwind-biased shock-capturing dissipation of the WENO-based scheme provides the
subgrid closure, and no explicit subgrid-scale model is used (§2.Z). Two axisymmetric
members share Δx = D_e/32 but differ in construction — a uniformly refined base grid
(L0hi) versus a single AMR level (L1) — and are retained to expose the sensitivity to the
refinement footprint itself, i.e. to the spatial region covered by the finest AMR level
(their comparison is reported in §3.X.2 and §3.X.4).

All statistics are accumulated in the solver at every time step (time-integrated first
and second moments); no snapshot averaging is used. The statistical sampling windows are
not identical across all members: the axisymmetric L0–L3 runs sample to t = 9.5 ms,
L4/L5 to 8.0 ms, and L6 to 17.25 ms; the three-dimensional L1–L3 runs to t ≈ 10.4 ms, L4
to 9.6 ms, and L5 to 10.2 ms. Consequences of this are noted where relevant, and the
results of this section are presented as a grid-sensitivity assessment rather than a
formal grid-convergence study.

The following metrics are extracted from the time-mean centreline profiles (all ⟨·⟩
quantities in this section are centreline values):

- **x_s** — the closure location of the first shock cell, defined as the axial position of
  maximum d⟨p⟩/dx on the centreline, i.e. where the reflected oblique-shock system
  reaches the axis;
- **x_n** — the axial position of the n-th principal compression maximum of ⟨p⟩ (x_1 the
  first peak exceeding the ambient pressure, excluding sub-ambient local maxima of the
  oscillation at the base of the pressure recovery immediately upstream of the closure);
- **S_n = x_{n+1} − x_n** — the spacings between successive principal compression maxima;
- **upstream and downstream window extrema** — the maximum of ⟨M⟩ over
  [x_s − 0.5 D_e, x_s − 0.02 D_e] and the minimum of ⟨M⟩ over
  [x_s + 0.02 D_e, x_s + 0.6 D_e] (a ±0.02 D_e buffer excludes the discrete jump itself),
  denoted ⟨M⟩_uw and ⟨M⟩_dw, together with ⟨p⟩/p∞ and ⟨ρ⟩/ρ_j at those locations. These
  are non-collocated extrema of time-mean profiles extracted from finite axial windows;
  they are *not* interpreted as one-sided shock-jump states, and no Rankine–Hugoniot
  consistency is implied.

The positional uncertainty is taken as u_x = max(0.01 D_e, Δx/2); the 0.01 D_e floor
reflects the sampling-window sensitivity of the peak locations, which exceeds Δx/2 on the
finest grids. Because the refinement footprints do not cover the entire measurement
domain, the *local* cell size at a measurement station can differ from the nominal finest
resolution: in the three-dimensional L5 configuration the finest level forms an annular
sleeve around the lip shear layer extending to x ≈ 0.75 D_e, so that both
x_s ≈ 0.88 D_e and x_1 ≈ 1.18 D_e lie in D_e/128-resolved cells. The axisymmetric members
do not require this caveat: in the D_e/512 and D_e/1024 grids the first-cell stations are
contained in the finest-level footprint (verified from the stored level boxes). Local
cell sizes at the stations are reported in Table 3.X-1.

## 3.X.2 Mean shock-cell geometry: bounded positional variation

Table 3.X-1 summarises the first-cell metrics for all thirteen configurations; the
corresponding centreline mean-profile overlay is shown in [Fig. AXIS-FAM].

**Table 3.X-1.** First-cell metrics for the thirteen configurations (online statistics;
lengths in units of D_e; pressures normalized by p∞). Δx_loc is the local cell size at
the first-cell stations x_s, x_1.

| Configuration | Δx (nominal) | Δx_loc | x_s | x_1 | S_1 | ⟨M⟩_dw | ⟨p⟩(x_1)/p∞ | p_rms(x_1)/p∞ |
|---|---|---|---|---|---|---|---|---|
| 2D-L0   | D_e/16   | D_e/16   | 0.906 | 1.140 | 1.318 | 1.02 | 1.58 | 0.025 |
| 2D-L0hi | D_e/32   | D_e/32   | 0.922 | 1.134 | 1.351 | 1.01 | 1.60 | 0.027 |
| 2D-L1   | D_e/32   | D_e/32   | 0.891 | 1.144 | 1.325 | 1.01 | 1.60 | 0.044 |
| 2D-L2   | D_e/64   | D_e/64   | 0.883 | 1.175 | 1.222 | 0.97 | 1.61 | 0.145 |
| 2D-L3   | D_e/128  | D_e/128  | 0.895 | 1.101 | 1.386 | 0.64 | 1.77 | 0.214 |
| 2D-L4   | D_e/256  | D_e/256  | 0.873 | 1.131 | 1.153 | 0.77 | 1.73 | 0.172 |
| 2D-L5   | D_e/512  | D_e/512  | 0.886 | 1.103 | 1.253 | 0.48 | 1.79 | 0.219 |
| 2D-L6   | D_e/1024 | D_e/1024 | 0.860 | 1.105 | 1.167 | 0.44 | 1.72 | 0.260 |
| 3D-L1   | D_e/16   | D_e/16   | 0.906 | 1.149 | 1.325 | 1.02 | 1.58 | 0.011 |
| 3D-L2   | D_e/32   | D_e/32   | 0.922 | 1.123 | 1.345 | 1.00 | 1.61 | 0.015 |
| 3D-L3   | D_e/64   | D_e/64   | 0.898 | 1.158 | 1.320 | 0.95 | 1.63 | 0.027 |
| 3D-L4   | D_e/128  | D_e/128  | 0.887 | 1.140 | 1.326 | 0.85 | 1.73 | 0.046 |
| 3D-L5   | D_e/256  | D_e/128  | 0.877 | 1.178 | 1.301 | 0.94 | 1.67 | 0.035 |

The closure location remains within x_s/D_e = 0.860–0.922 and the first compression
maximum within x_1/D_e = 1.101–1.178 across the entire matrix — total spans of 0.062 D_e
and 0.076 D_e respectively (deviations from the 13-configuration ensemble means bounded
by ±0.033 D_e for x_s and ±0.041 D_e for x_1), spanning a 64-fold change in resolution
and a change of dimensionality. Neither sequence is monotonic in Δx, so these bounded
ranges demonstrate *limited sensitivity* of the gross axial shock-cell geometry
(bounded variation, not convergence; cf. §3.X.1). Part of the residual variation is
attributable to the refinement footprint alone: the two D_e/32 members differ by
Δx_s = 0.031 D_e — half the full-matrix span — and by ≈40% in median near-field p_rms.
The first spacing S_1 varies over 1.15–1.39 D_e (a wider relative spread than x_s or
x_1), and the downstream spacings behave differently in the two families: the
three-dimensional sequences decay monotonically (S_1 ≈ 1.30 D_e to S_5 ≈ 0.85 D_e at
L5), whereas in the fine axisymmetric solutions the peaks beyond x_4 are progressively
smeared beyond reliable identification by the resolved unsteadiness, and the extracted
spacings there are not considered reliable.

The weak grid dependence of the leading cells is physically plausible: the quasi-steady
near-field shock-cell structure is controlled by integral parameters — the nozzle
pressure ratio, the fully expanded jet Mach number, and the integral inflow fluxes (a
discrete audit of the inflow mass flux gives ṁ/ṁ_analytic = 1.0000 on every level;
§3.Y) — none of which requires fine resolution. The same argument explains why the
coarsest solutions reproduce the leading time-mean structure at all: they obtain
approximately the right mean geometry while the resolved dynamics are wrong (their
fluctuation fields are essentially empty, §3.X.4), so the agreement of coarse-grid mean
fields with experiment must not be read as model fidelity.

## 3.X.3 Closure-region state: strong sensitivity to resolution and dimensionality

In contrast to the positions, the local state in the closure region is strongly grid- and
dimensionality-dependent. The coarse members of both families (axisymmetric D_e/16–D_e/64,
three-dimensional D_e/16–D_e/64) return a near-sonic downstream window extremum,
⟨M⟩_dw ≈ 0.95–1.02, with first-peak pressures ⟨p⟩(x_1)/p∞ ≈ 1.58–1.63. The fine
axisymmetric solutions (D_e/128–D_e/1024) instead develop a distinctly subsonic minimum,
⟨M⟩_dw = 0.44–0.77, with ⟨p⟩(x_1)/p∞ = 1.72–1.79; the fine three-dimensional solutions
(D_e/128–D_e/256) show the same tendency in weaker form, ⟨M⟩_dw = 0.85–0.94. Grid
refinement therefore substantially strengthens the resolved centreline compression at the
first-cell closure.

Because these quantities are window extrema of time-mean profiles, they cannot by
themselves classify the shock-reflection topology (regular versus Mach reflection). At
M_j = 1.42 the case lies near the boundary at which a Mach disk first forms — Mach-disk
formation in underexpanded jets is typically reported for M_j ≳ 1.5 (cf. Addy 1981) — and
Panda & Seasholtz (1999) describe the first-cell closure at this condition as a conical
compression zone rather than a developed Mach disk. Establishing the instantaneous
topology would require the identification, in instantaneous meridional fields, of a
finite-radius near-normal shock segment with a triple-point ring and slip line, with the
segment radius resolution-independent — this is left outside the claims of the present
section.

## 3.X.4 Resolved unsteadiness

The unsteady content shows the opposite sensitivity hierarchy to the positions
([Fig. RMS-FAM]). Along the axisymmetric family the maximum near-field centreline
density-fluctuation rms grows from approximately 0.09 ρ_j at D_e/16 to 0.25 ρ_j at
D_e/32 — where its location also moves upstream, from x/D_e ≈ 4.6 at D_e/16 to
x/D_e ≈ 3.4 — and reaches 0.34–0.38 ρ_j for D_e/64 and finer grids, with the peak locked
near the first-cell closure (x/D_e ≈ 0.9). The variation within the fine branch (≈12%,
non-monotonic) is within the influence of the differing sampling windows, so the
fine-branch plateau is reported as strongly reduced sensitivity, not statistical
convergence.

Dimensionality affects the fluctuation amplitudes far more than the mean geometry. At the
first compression maximum the axisymmetric fine grids give p_rms/p∞ = 0.17–0.26, whereas
the three-dimensional solutions at comparable local resolution give only 0.03–0.05. This
persistent disparity is consistent with the stronger coherence enforced by the
axisymmetric (m = 0) constraint, which excludes non-axisymmetric vortex-ring instability
and three-dimensional breakdown. Part of the near-shock rms in both families reflects
shock-position oscillation rather than turbulence, since p′ ≈ −(d⟨p⟩/dx)·x_s′ near a
moving compression, where x_s′ denotes the instantaneous fluctuation of the shock
position about its mean; attributing the axisymmetric excess to coherent vortex-ring
impingement on the shock would require conditional averaging or
shock-motion/shear-layer coherence diagnostics, which are not attempted here. The
three-dimensional family itself shows a non-monotonic fluctuation trend (median
near-field centreline p_rms/p∞ of 0.023, 0.081, 0.089, 0.074, 0.068 for L1–L5), peaking
at intermediate resolution.

## 3.X.5 Time-mean centreline supersonic extent

A complementary integral indicator is the time-mean centreline supersonic extent,
X_M = sup{x ≤ x_w : ⟨M⟩(x) ≥ 1}, evaluated on a fixed analysis window x_w = 12 D_e, with
⟨M⟩ constructed from the mean fields as ⟨M⟩ = |⟨u⟩|/√(γR⟨T⟩) (a Mach number of the
means, not the mean of the instantaneous Mach number). The far field is excluded because
the stored mean profiles there are affected by the outflow boundary; weak ⟨M⟩ ≥ 1
patches beyond the window (present on the coarser and finest axisymmetric grids) do not
enter this metric. The definition deliberately avoids the term potential-core length: in
an underexpanded jet ⟨M⟩ crosses unity repeatedly through the shock-cell train, and X_M
responds jointly to entrainment and mixing, shock-cell amplitude and phase, numerical
shock dissipation, and the sampling window. Values that remain supersonic at the window
edge are lower bounds (right-censored, marked ≥); this affects only axisymmetric
members.

**Table 3.X-2.** Time-mean centreline supersonic extent X_M/D_e, keyed by nominal
resolution (both D_e/32 axisymmetric members shown). † the 2D D_e/128 profile is
window-limited: weak supersonic intervals recur beyond x_w (to x/D_e ≈ 15.9).

| Δx | D_e/16 | D_e/32 | D_e/64 | D_e/128 | D_e/256 | D_e/512 | D_e/1024 |
|---|---|---|---|---|---|---|---|
| 3D | 9.91 | 7.42 | 7.84 | 8.75 | 9.94 | — | — |
| 2D | ≥12 | ≥12 (L0hi); ≥12 (L1) | 10.7 | ≥12 † | 7.9 | 4.8 | 4.2 |

The three-dimensional family is non-monotonic: a pronounced shortening from D_e/16
(where the sonic crossing at x/D_e = 9.9 is a genuine interior crossing) to D_e/32,
followed by a monotonic recovery to D_e/256. The shortest extent and the largest
centreline fluctuation levels both occur in the intermediate-resolution range. These
correlated trends are consistent with a resolution-dependent balance between numerical
dissipation, coherent shear-layer dynamics, and three-dimensional breakdown: the coarsest
grid strongly attenuates the resolved instability of the lip shear layer (near-field rms
a factor ≈3 below the finer-grid levels and ≈4 below the intermediate-resolution
maximum), the intermediate grids resolve the dominant coherent structures while their
secondary three-dimensional breakdown may remain under-resolved, and the finer grids
reduce the fluctuation level again while the supersonic extent recovers. The present
diagnostics do not establish a maximum of entrainment efficiency or an energy-cascade
blockage at intermediate resolution, and no such claim is made.

The axisymmetric branch behaves differently: beyond D_e/128 the supersonic extent
shortens monotonically with refinement (7.9, 4.8, 4.2). The axisymmetric constraint
suppresses non-axisymmetric secondary instabilities, so the fine-branch shortening is a
trend consistent with the resolved m = 0 structures (vortex rings) becoming progressively
stronger and more persistent in the absence of the three-dimensional breakdown path; the
axisymmetric sequence should therefore be read as an m = 0-constrained family, not as an
approximation to the three-dimensional refinement path. (Axisymmetric flow retains the
geometric vorticity-stretching term u_r ω_θ / r, so this behaviour is attributed to the
missing azimuthal degrees of freedom rather than to planar-two-dimensional
inverse-cascade phenomenology.)

## 3.X.6 Comparison with the Panda–Seasholtz experiment

The experimental reference values are: first centreline compression maximum
x_1/D_e ≈ 1.20 for M_j = 1.42, from the 37-phase-averaged centreline density dataset
accompanying Panda & Seasholtz (1999) (EPAPS file M142DEN); an independent time-mean
centreline dataset at M_j = 1.43 from the same archive gives 1.213. For the other two
conditions of the validation matrix — quoted here because the same extraction pipeline is
applied across M_j = 1.19–1.80 (§3.Y) — the first compression maximum is
x_1/D_e ≈ 0.77 at M_j = 1.19 and the Mach-disk location x_MD/D_e ≈ 1.39 at M_j = 1.80
(maximum centreline density gradient).

For M_j = 1.42 the finest three-dimensional solution gives x_1/D_e = 1.178
([Fig. EXP-CMP]), a deviation of −1.9% from the experimental 1.20, with
x_s/D_e = 0.877 against an experimental closure estimated from the steepest centreline
density rise of the phase-averaged dataset (0.88–0.90, depending on smoothing and
differentiation stencil) — i.e. just below that range, a deviation of 0.3–2.6% depending
on the reference point — and a first spacing S_1 = 1.30 D_e within 2% of the experimental
mean spacing of 1.28 D_e. Farther downstream the agreement deteriorates in the
systematic way quantified in §3.Y: the computed shock-cell spacing is approximately
7–10% shorter than measured (an accumulated upstream shift of ≈0.3 D_e by
x/D_e ≈ 4–6), and the computed density modulation persists farther downstream than in
the experiment.

The validation claim is therefore restricted as follows:

- **validated (quantitative)** — near-field time-mean shock-cell geometry: x_1 and S_1
  within 2% of experiment at the finest three-dimensional resolution, and the first-cell
  closure location within 0.3–2.6% of the estimated experimental range;
- **semi-quantitative** — time-mean density profiles over the first three to four cells,
  reported together with the 7–10% downstream spacing deficit;
- **not validated** — downstream mixing and shock-cell decay, fluctuation amplitudes
  (the experimental fluctuation data are phase-locked coherent components at the screech
  frequency, not total rms, and are not directly comparable), screech tones, and acoustic
  radiation;
- **axisymmetric (RZ) results** — used only for positional metrics and for the
  sensitivity scans of this section; they are not used to validate three-dimensional
  dynamics.

## 3.X.7 Summary

The thirteen-configuration matrix establishes a clear hierarchy of grid demands. The
axial positions of the leading shock-cell structures are the most robust quantities,
bounded within ±0.033 D_e (x_s) and ±0.041 D_e (x_1) of their ensemble means across all
resolutions and both dimensionalities, and agreeing with the experimental first
compression maximum and leading spacing to within 2% at the finest three-dimensional
resolution, with the closure location just below the estimated experimental range. The
local closure-region state is far more sensitive: only grids of D_e/128 and finer
strengthen the centreline compression toward a subsonic downstream minimum, and the fine
axisymmetric and three-dimensional branches still differ (downstream-window Mach minima
of 0.44–0.77 versus 0.85–0.94). Second-order statistics are the most demanding: they
show strongly reduced sensitivity only for D_e/64 and finer in the axisymmetric family
(a fine-branch variation of ≈12%, within the influence of the differing sampling
windows), remain non-monotonic in the three-dimensional family up to D_e/256, and differ
between dimensionalities by a factor of ≈5 at matched nominal resolution
(p_rms/p∞ = 0.172 versus 0.035 at Δx = D_e/256). Finally, the shock-reflection topology
of the first cell cannot be classified from centreline time-mean data at any of the
present resolutions. Positional metrics may therefore be extracted with bounded
(≤0.041 D_e) uncertainty from comparatively coarse grids — whose solutions reproduce the
mean geometry only, not the resolved dynamics — whereas fluctuation-based conclusions
require the fine three-dimensional configurations, and the axisymmetric model, however
fine, is not a substitute for them.
