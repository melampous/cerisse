# Grid-sensitivity figure set — descriptions, definitions, extraction, captions

> **Superseded-method notice (2026-07-23).** The historical extraction notes
> below retain the former pressure-gradient definition of \(x_{fc}\), block
> statistics, and the former \(0.05D_e/0.071D_e\) experimental intervals.
> They must not be used to regenerate the current thesis figures. The
> canonical method is the raw mean-centreline-density maximum positive
> gradient between the first expansion minimum and subsequent principal
> compression maximum, with local quadratic refinement. Current values and
> assets are recorded in §19 of
> `U_jet/analysis/HANDOFF_UJET_2026-07-10.md` and in the compiled thesis
> Section 3.4.2.

All quantities from the accumulated (per-time-step, on-line) statistics over each
case's frozen-grid statistics window. Normalisation: D_e = 25.4 mm, p_inf =
99780 Pa, rho_j = 1.6413 kg/m3 (fully expanded, Panda convention), U_j = 416.2
m/s (fully expanded, M_j = 1.42, T0 = 300 K). The horizontal resolution axis
is D_e/dx_min and denotes the nominal maximum hierarchy resolution, not the
local resolution of every extracted feature. The 2D sequence spans L0--L6 and
the 3D sequence spans L1--L5. Level naming follows Table 3.5 of the thesis: a
run with maximum refinement level n is labelled "Ln" in all figure panels,
e.g. the L6
axisymmetric grid and the L5 three-dimensional grid.

## Fig. A — paperfig_shockpos_vs_res (mean shock-cell metrics)

Panel (a): x_fc and x_1 vs D_e/dx_min.
- x_fc (first shock-cell closure): argmax of d<p>/dx on the centreline,
  searched in x/D_e in [0.4, 1.6]. This topology-neutral gradient criterion is
  robust to shock-motion smearing.
- x_1 (first principal mean-pressure maximum): searched DOWNSTREAM of x_fc.
  Canonical definition sentence: "Downstream of x_fc, the successive
  principal mean-pressure maxima are denoted by x_n. Each accepted maximum
  must satisfy <p>/p_inf > 1, have a prominence greater than 0.03 in
  normalised pressure, and be separated from the adjacent accepted maximum
  by at least 0.5 D_e. The peak-to-peak shock-cell spacings are then defined
  as S_n = x_{n+1} - x_n."
- Physical chain: expansion -> steepest recompression (x_fc) -> first
  pressure maximum (x_1) -> second pressure maximum (x_2); S_1 = x_2 - x_1
  is a PEAK-TO-PEAK spacing, not a nozzle-to-shock length and not x_fc.
- Profiles sampled at each case's finest AMR level (finest-wins composite),
  box-filtered 0.05 D_e; the discrete feature locations are estimated to
  sub-cell precision by a three-point parabolic fit (this does not increase
  the true spatial resolution of the field).
- Error bars: standard error over 3-8 statistics sub-blocks obtained by
  tau-weighted differencing of the accumulated moments (blocks < 0.2 ms merged).
- Shaded-band centres: positions extracted from the mean of the 37
  phase-resolved M_j = 1.42 centreline density profiles of Panda & Seasholtz
  (1999), giving x_fc = 0.888 and x_1 = 1.199. The conservative half-width
  +-0.05 D_e is half the 0.10 D_e sampling step of the independently acquired
  M_j = 1.43 time-mean traverse.
  Estimated extraction-resolution intervals, NOT author-quoted confidence
  intervals; near-nozzle stray-light systematics (0.5 < x/D_e < 1.5,
  no-flow correction region) are not included.
  Sim rho-based vs p-based positions differ by < 0.01 D_e.

Panel (b): S_n/D_e vs D_e/dx_min, where S_n = x_{n+1} - x_n. The phase-mean
experimental references are S_1..3 = 1.359, 1.186, 1.130. Each coloured band
has the quadrature-propagated half-width sqrt(2)*0.05 = 0.071 D_e, assuming
independent local peak-location errors. The 2D spacing sequences are shown as
light-grey dashed lines because they are retained only for illustration. Open
circle, square and triangle markers distinguish S_1, S_2 and S_3,
respectively.

Readings: all x_fc values lie inside the estimated experimental band. The
finest 3D x_1 differs from the phase-mean experimental extraction by about
0.02 D_e; the finest 2D value is shifted upstream by about 0.09 D_e. The 3D
S_1 values remain within or close to the corresponding experimental band,
whereas S_2 and S_3 shorten on the finer hierarchies. The uncertain 2D S_2 and
S_3 estimates are not used to support a trend.

Caption (EN):
"Hierarchy-level mesh sensitivity of the mean shock-cell metrics. (a) First
shock-cell closure x_fc (maximum of d<p>/dx on the centreline, x/D_e in
[0.4, 1.6]) and
first main compression peak x_1 (first local maximum of <p>/p_inf above the
ambient level; prominence 0.03 p_inf, minimum separation 0.5 D_e) versus the
nominal maximum hierarchy resolution D_e/dx_min. (b) Individual shock-cell
spacings S_1/D_e, S_2/D_e and S_3/D_e. Open circles/dashed: 2D; filled
squares/solid: 3D. The reliability of the 2D S_2 and S_3 estimates is discussed
in the accompanying text rather than encoded using separate markers. Error
bars: standard error over 3-8 statistics sub-blocks
(tau-weighted differencing). Shaded intervals are estimated spatial and
feature-extraction resolution bounds, not confidence intervals reported by
Panda & Seasholtz (1999): their centres are extracted from the mean of the 37
phase-resolved M_j = 1.42 profiles. Position intervals are +-0.05 D_e and
spacing intervals follow from quadrature propagation of the two bounding
peak locations (+-0.071 D_e), assuming independent peak errors."

## Fig. B — paperfig_rmsheat_2x2 (r.m.s. fluctuation heatmaps)

2x2 heatmaps: rows = 2D family (7 levels) and 3D family (5 levels); columns =
p_rms/p_inf (magma) and u_x,rms/U_j (viridis); x in [0, 8] D_e; y = discrete
level bands labelled L0--L6 and L1--L5, with no inter-level interpolation.
- q_rms = sqrt(<q^2> - <q>^2) from accumulated first/second moments:
  statistically converged r.m.s. fluctuation (not snapshot std).
- Colours are absolute non-dimensional amplitudes (no per-row normalisation).
  The 2D and 3D panels share a common limit for each variable, set from the
  99.5th percentile of the combined values.
- Display smoothing 0.05 D_e along x only.
- Centreline r.m.s. does not assess shear-layer resolution (the layer lies at
  r/D_e ~ 0.5): that role belongs to the instantaneous vorticity/schlieren
  figures.

Readings: (i) coarse-grid suppression (D_e/dx_min = 16 rows nearly dark; 2D L0
quasi-steady); (ii) the p_rms stripes (shock stations) align across rows for
D_e/dx_min >= 64: the spatial organisation of the unsteady wave system is similar while
amplitudes grow; (iii) fluctuation onset moves upstream with refinement
(earlier roll-up); (iv) 2D amplitudes are 3-5x the 3D ones at matched nominal
maximum hierarchy resolutions
(two-dimensional over-coherence).

Caption (EN):
"Centreline r.m.s. fluctuations q_rms = (<q^2> - <q>^2)^{1/2} from the
accumulated moments over the full statistics window: (a, c) p_rms/p_inf and
(b, d) u_x,rms/U_j for the 2D (top) and 3D (bottom) grid families. Each row
is one simulation labelled by its AMR hierarchy; rows are discrete bands
without inter-level interpolation. Colours give absolute non-dimensional
amplitudes, and each variable uses the same colour limit in the two coordinate
formulations. Profiles are box-filtered over 0.05 D_e. Centreline r.m.s.
levels do not assess shear-layer resolution, since the shear layer lies at
r/D_e ~ 0.5."

## Fig. C — paperfig_rmsint_vs_res (integral fluctuation metrics)

Panel (a): A_p^{[a,b]} = [ (b-a)^{-1} D_e^{-1} int_{aD_e}^{bD_e}
(p_rms/p_inf)^2 dx ]^{1/2} for [0.05, 2] (the near-nozzle region containing
the first shock cell and the upstream part of the second)
and [2, 6] (downstream shock train and mixing region).
Panel (b): fluctuation-intensity-weighted centroids x_c = int x q_rms^2 dx /
int q_rms^2 dx over [0.05, 8] D_e, for q = p and q = u_x. These are not
source locations or energy centroids.
Error bars: block SE as in Fig. A. Resolvability criterion: adjacent-grid
differences below ~2 SE are within the time-sampling uncertainty.

Readings: x_c(p) moves upstream from 6.3-6.5 (coarse) and remains near 4.6-4.8
on the finer 2D grids. The 2D A_p^{[2,6]} values remain between 0.28 and 0.30
from L3 to L6, while A_p^{[0.05,2]} still grows from L5 to L6
(0.187+-0.006 -> 0.213+-0.009, > 2 SE):
refinement increases the resolved near-field fluctuation intensity without
substantially changing the spatial extent or axial organisation of the
unsteady region — the ILES-consistent statement; "grid convergence of the
r.m.s. field" must not be claimed. One resolvable 3D grid effect: x_c(p)
rises to 5.25+-0.08 on L5 (vs 4.25-4.43 on L3--L4, >> 2 SE): the
finest 3D grid sustains the downstream shock-train fluctuations (x ~ 4-6)
better, consistent with the recovery of the supersonic-core length.

Caption (EN):
"(a) Band-integrated fluctuation intensity A_p over the near-nozzle region
[0.05, 2] D_e and the downstream shock-train region [2, 6] D_e. (b)
Fluctuation-intensity-weighted centroids x_c over [0.05, 8] D_e for p_rms and
u_x,rms.
Symbols as in the shock-position figure; error bars denote the standard
error over 3-8 statistics sub-blocks. Beyond D_e/dx_min ~ 128 the resolved
near-field intensity still increases while the spatial organisation of the
unsteady region is unchanged."

## Figs. D, E — paperfig_axis2d_4var / paperfig_axis3d_4var (appendix)

Four stacked centreline profiles, all levels of one family per figure
(2D: L0--L6; 3D: L1--L5), turbo colour ramp ordered by resolution,
finest level heavier:
(a) M_qbar = |<u>| / sqrt(gamma R <T>): Mach number reconstructed from the
    time-averaged primitive variables (not the mean of instantaneous M);
    M = 1 dotted line.
(b) <rho>/rho_j with the Panda & Seasholtz (1999) centreline density
    (open circles; Rayleigh scattering measures density only).
(c) <u_x>/U_j.
(d) rho_rms/rho_j on a square-scaled ordinate (emphasises peaks; the maxima
    at x ~ 0.88 mark the oscillation of the first shock closure).
Display smoothing 0.05 D_e (means) / 0.15 D_e (rms). Statistics windows:
2D 1.5 ms per case (absolute windows differ per level and are stated in the
text); 3D 1.56-3.7 ms.

Role: complete archival record (appendix); the main text carries Figs. A-C.
Readings: 2D supersonic core shortens monotonically with refinement while
the 3D core length is non-monotonic (minimum on L2, recovery by L5);
the coarse-grid mean field tracks the experimental density surprisingly well
(inviscid wave skeleton + grid-insensitive integral quantities); the maximum
rho_rms locates at the first shock closure once resolved (2D L2 and finer,
3D L4 and finer).

Caption (EN):
"Centreline profiles for all grid levels of the [2D/3D] family: (a) Mach
number reconstructed from the time-averaged primitive variables, (b) mean
density with the Rayleigh-scattering measurement of Panda & Seasholtz
(1999), (c) mean axial velocity, (d) density r.m.s. fluctuation
(square-scaled ordinate). Curves are box-filtered over 0.05 D_e (means) and
0.15 D_e (r.m.s.); colours order the levels by resolution."

## Fig. A0 — paperfig_metricdefs_schlieren (metric definitions, illustration)

Two stacked panels for the 3D L5 solution.
(a) Centreline <p>/p_inf with the operational definitions marked: the
expansion branch, the steepest recompression x_fc (blue dashed), and the
principal mean-pressure maxima x_1..x_4 (red markers), searched downstream
of x_fc with the acceptance criteria above.
(b) Time-mean numerical schlieren exp(-8 |grad<rho>|/|grad<rho>|_p99.5) on
the z = 0 plane, gradient evaluated per AMR level before compositing
(finest-wins); x_fc and x_n carried down as vertical lines, and the
peak-to-peak shock-cell spacings S_1..S_3 annotated as double arrows.
Values (3D L5, full window): x_fc = 0.876, x_n = 1.178, 2.479, 3.584,
4.588; S_n = 1.301, 1.105, 1.004.

Caption (EN):
"Operational definitions of the mean shock-cell metrics, illustrated with
the 3D L5 solution. (a) Centreline mean pressure: the expansion branch ends
at the steepest recompression x_fc = argmax d<p>/dx; downstream of x_fc the
successive principal mean-pressure maxima x_n must satisfy <p>/p_inf > 1,
have a prominence greater than 0.03 in normalised pressure, and be separated
by at least 0.5 D_e. Discrete feature locations are estimated to sub-cell
precision by a three-point parabolic fit. (b) Time-mean numerical schlieren
of the same solution (z = 0 plane; the density gradient is evaluated on each
AMR level before compositing); the peak-to-peak shock-cell spacings
S_n = x_{n+1} - x_n are indicated."
