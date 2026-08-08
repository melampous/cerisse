# Pure GP-IBM Gate 1: source-free smooth curved-wall accuracy

Status date: 2026-07-22

## Decision

\[
\boxed{\text{Gate 1: PASS for the restricted current-source pure GP-IBM}}
\]

The current pure full-Cartesian GP-IBM recovers approximately second-order
fixed-time accuracy for a source-free smooth Euler solution around a circular
slip wall. The conclusion is supported by four asymptotic grids and two
normalized grid phases. It applies to L1/L2 solution norms, the first
wall-adjacent fluid-cell layer, a fixed physical near-wall band, BI pressure,
and reconstructed wall-normal mass flux. It is not a claim that every Linf
sequence is uniformly second order.

No solver, IBM, WENO, BIC, limiter, or case source was changed during this
gate. One clean CUDA executable was used for all runs.

## Method under test

The verification solution is the stationary, source-free two-dimensional
isentropic Euler vortex centred on the immersed circle. With
\(\tilde r^2=((x-x_c)^2+(y-y_c)^2)/R^2\), \(\beta=5\), and
\(\gamma=1.4\),

\[
u=-\frac{\beta}{2\pi}\frac{y-y_c}{R}
   \exp\!\left(\frac{1-\tilde r^2}{2}\right),
\qquad
v= \frac{\beta}{2\pi}\frac{x-x_c}{R}
   \exp\!\left(\frac{1-\tilde r^2}{2}\right),
\]

\[
T=1-\frac{(\gamma-1)\beta^2}{8\gamma\pi^2}
       \exp(1-\tilde r^2),
\qquad
\rho=T^{1/(\gamma-1)},
\qquad
p=\rho T=\rho^\gamma.
\]

The circle has \((x_c,y_c)=(0.5,0.5)\) and \(R=0.2\). It is an exact
streamline, so the wall satisfies

\[
u_n=0,
\qquad
\partial_n u_t=0,
\qquad
\partial_n p=\frac{\rho u_t^2}{R}.
\]

The numerical configuration is the frozen production candidate:

- shared solid-side ghost-cell state;
- cell-average-aware quadratic BI-CWLS;
- expanded visible support, with no degree-one fallback in this campaign;
- numerical-trace Euler-slip curvature pressure closure;
- one-sided entropy extension and ideal-gas EOS recovery;
- one solid-side ghost-cell layer;
- full Cartesian control volumes and face areas;
- LLF-WENO-Z5 and SSPRK(4,3);
- conservative full-cell positivity limiting;
- second-order BI surface recovery.

The runtime manifest reports
`family=pure_shared_gp_full_cartesian`, `pure_gp_production=1`, and
`cut_control_compiled=0`. No cut-control or aggregate update symbol is defined
in the executable.

## Protocol

All runs use the same fixed 4096-segment circular geometry and terminate at
\(t=0.01\). The complete sequence is \(N=64,96,128,192,288\); \(N=64\) is
retained as a pre-asymptotic diagnostic, while formal fits below use
\(N=96,128,192,288\).

Two normalized geometry phases are used:

- `axis`: domain origin \((0,0)\);
- `phase037_023`: domain origin shifted by \((-0.37h,-0.23h)\), while the
  physical circle remains fixed.

Errors are computed against exact conservative cell averages using four-point
Gauss--Legendre quadrature in each coordinate. The reported primitive fields
are \(\rho,u,v,p,T\). The first-wall-layer mask consists of active fluid cells
sharing a Cartesian face with a solid cell. The fixed physical band is
\(0\le r-R\le0.05\).

## Fixed-time solution convergence

The table gives least-squares fitted orders over \(N=96,128,192,288\).

### L2 orders

| Phase and region | rho | u | v | p | T |
|---|---:|---:|---:|---:|---:|
| axis, all fluid | 1.960 | 2.052 | 2.052 | 2.003 | 2.003 |
| shifted, all fluid | 2.092 | 2.076 | 2.074 | 2.018 | 2.006 |
| axis, first wall layer | 1.903 | 1.952 | 1.952 | 1.952 | 1.968 |
| shifted, first wall layer | 2.023 | 2.043 | 2.038 | 2.018 | 1.983 |
| axis, fixed physical band | 1.962 | 2.161 | 2.161 | 2.010 | 2.008 |
| shifted, fixed physical band | 2.091 | 2.226 | 2.223 | 2.055 | 2.015 |

The corresponding L1 fits are 2.00--2.26 globally, 1.90--2.08 in the first
wall layer, and 2.01--2.29 in the fixed physical band. These results satisfy
the complete-state approximate-second-order gate in both phases.

### Linf qualification

The first-wall-layer shifted-phase Linf fits are 1.979--2.047. The axis-phase
fits are 1.837 for density, 1.897 for velocity, 1.891 for pressure, and 1.934
for temperature. Global density and pressure Linf fits are 1.726/1.814 for
the axis phase and 1.777/1.851 for the shifted phase. Their maxima move among
different near-wall cells as the grid is refined. On the finest pair,
\(N=192\to288\), every axis-phase first-layer variable recovers an observed
order above 2.13.

Accordingly, this gate establishes approximately second-order L1/L2
convergence and asymptotic Linf recovery; it does not establish a uniform
strict-second-order Linf bound independent of grid phase.

## BI pressure and impermeability

The independent BI audit reconstructs a fluid-side trace on all 4096 surface
elements. Every trace is valid and quadratic; the support count is 26--37 and
the maximum reported condition number is 43.31.

Fitted orders over \(N=96,128,192,288\) are:

| Phase | rho un L1 | rho un L2 | rho un Linf | signed BI integral | absolute BI flux | BI pressure L2 |
|---|---:|---:|---:|---:|---:|---:|
| axis | 2.564 | 2.618 | 2.680 | 2.231 | 2.564 | 2.248 |
| shifted | 2.630 | 2.626 | 2.595 | 2.368 | 2.630 | 2.268 |

Here `signed BI integral` denotes

\[
\frac{\left|\int_\Gamma \rho u_n\,dS\right|}{\dot m_{\rm ref}},
\]

and `absolute BI flux` denotes

\[
\frac{\int_\Gamma |\rho u_n|\,dS}{\dot m_{\rm ref}}.
\]

At \(N=288\), the independent reconstructed BI values are:

| Phase | rho un L1 | rho un L2 | rho un Linf | signed BI/ref | absolute BI/ref | BI pressure rel. L2 |
|---|---:|---:|---:|---:|---:|---:|
| axis | 4.984e-6 | 5.934e-6 | 1.717e-5 | 6.381e-6 | 7.939e-6 | 2.027e-5 |
| shifted | 4.764e-6 | 6.130e-6 | 2.025e-5 | 6.276e-6 | 7.588e-6 | 2.003e-5 |

The reconstructed physical-wall impermeability error therefore decreases
faster than second order and is not stalled at the 0.08--0.10% level observed
in the one-grid inverse-MOC Gate 0 diagnostic.

The signed Cartesian fluid/solid numerical exchange also converges, with
orders 2.189 and 2.149 over the full \(N=64\ldots288\) sequence. Its absolute
crossing-face exchange remains \(O(1)\), as expected for oblique Cartesian
faces around a curved wall, and is not interpreted as pointwise physical-wall
leakage. The active-fluid Cartesian mass budget closes to at most
\(2.74\times10^{-13}\) after normalization.

## Limiter and reconstruction path

Across both phases there are 246 complete SSPRK steps and 984 forward-Euler
brackets. Every bracket reports `HIGH_ORDER theta_min=1`; every stage
admissibility audit passes. There is no global-theta fallback, all-low final
fallback, or rejected step.

The two phases contain 864 and 859 GP targets, respectively. All 1723 targets
use quadratic BI-CWLS and none uses the linear fallback. Thus the observed
order is not produced by a changing fallback population, and the positivity
limiter is transparent in smooth flow.

## Acceptance and limitations

Gate 1 is closed because:

1. global, fixed-band, and first-wall-layer complete-state L1/L2 errors are
   approximately second order in both phases;
2. BI pressure converges above second order;
3. \(\rho u_n\) L1/L2/Linf and both signed and absolute reconstructed BI mass
   fluxes converge above second order;
4. the two phases have the same asymptotic trend and similar finest-grid
   errors;
5. the limiter is exactly transparent and no BI-CWLS downgrade occurs.

The PASS is restricted to a two-dimensional Cartesian, uniform level-0,
fixed smooth analytic wall, single-component ideal-gas Euler-slip problem
with LLF-WENO-Z5 and SSPRK(4,3). It does not certify feature edges, narrow
gaps, general faceted STL curvature, AMR, R-Z, Navier--Stokes wall heat flux
or shear, three-dimensional geometry, or moving solids. It also does not
replace the Gate 2 smooth-nozzle qualification.

## Evidence

- Raw and derived data:
  `IBM/cases/validation/mms_ibm/diagnostics/gate1_pure_20260722`
- Axis L2 convergence:
  `formal/axis/CONVERGENCE_ASYMPTOTIC.md`
- Shifted L2 convergence:
  `formal/phase037_023/CONVERGENCE_ASYMPTOTIC.md`
- Fixed physical band:
  `formal/FIXED_PHYSICAL_BAND_R0P05.md`
- Fixed-band reproduction tool: `postprocess_fixed_physical_band.py`
- BI impermeability and load audit:
  `formal/IMPERMEABILITY_LOAD_64_288.md`
- Reproduction metadata: `PROVENANCE.md`
- Recursive evidence hashes: `SHA256SUMS`
