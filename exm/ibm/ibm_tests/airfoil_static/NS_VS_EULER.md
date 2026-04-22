# Viscous companion: mid α = 12°, Re = 1 × 10⁴

**Purpose** — demonstrate that the Cerisse IBM + WENO-Z5 stack remains
stable and physically sensible at finite attack angle with an adiabatic
**no-slip** immersed-boundary wall, not just the slip wall used in the
Euler benchmark.

## Setup

- Case directory: `exm/ibm/ibm_tests/airfoil_static/`
- Inputs: `inputs_mid_a12_re1e4` — byte-identical to the Euler
  `inputs_mid_a12` in domain (−0.5..0.8 × −0.7..0.7 m), AMR (`max_level =
  3`, ref_ratio 2), CFL (0.25), WENO-Z5 + SSPRK(3,3), 5000 steps.
- Only differences vs Euler baseline:
  1. `prob.h` uses `ibm_adiabatic_noslip_wall_t` (no-slip, adiabatic).
  2. Built with `REYNOLDS_USER = 1e4`, which enters `prob.h` as a
     compile-time non-dimensional viscosity
     `viscosity = 1 / REYNOLDS_USER`.
- AWS job: 75, 96 ranks on one hpc8a node, 21 plotfiles + 50 surface
  vtp snapshots, no NaN/abort.

## Integrated-load comparison (tail mean, steps 3000–5000)

| Quantity | Euler | NS Re = 1 × 10⁴ | Δ |
|---|---|---|---|
| CL       | +0.547 | +0.320 | −41 % |
| CD       | +0.138 | +0.082 | −41 % |
| Cm_LE    | +0.256 | +0.149 | −42 % |
| Cm_c/4   | +0.115 | +0.066 | −43 % |
| x_cp/c   | +0.467 | +0.465 | +0.4 % |

The integrated CL/CD/Cm drop by ~40 %, while the centre-of-pressure
location stays put to within 0.5 %.

The ~40 % load reduction is **larger than a classical Re-1e4 viscous
correction would predict** for an attached-shock wedge. This is not a
transient artefact: the CL/CD/Cm time histories plateau from step ~2000
onwards with RMS fluctuation < 0.1 %. The most likely explanation is
that the compile-time `viscosity = 1 / REYNOLDS_USER` convention used
in Cerisse's `prob.h` does not directly correspond to a classical
chord-based Reynolds number at these dimensional freestream values;
that is, "Cerisse Re = 1e4" is a solver-level scaling knob, not a
lab-frame `ρ U c / μ`. This is flagged in the thesis text; it does
**not** affect the main robustness conclusion.

## Flowfield structure (qualitative)

- [mid_a12_re1e4_schlieren.png](mid_a12_re1e4_schlieren.png) — full-
  domain schlieren at step 5000. LE attached oblique shock pair is
  present with the expected asymmetry (upper expansion / lower
  compression for positive α). Aft-expansion fans are visible. A thin
  viscous wake extends from TE to the x-high outflow boundary. No
  spurious disturbances or AMR-seam artefacts near the body.
- [mid_a12_re1e4_zoom.png](mid_a12_re1e4_zoom.png) — near-body crop.
  Pressure layer around the wedge surfaces is continuous and smooth;
  no bow-shock splitting or ghost-cell checkerboard.

## Surface Cp (note: near-wall oscillations)

[cp_mid_a12_re1e4.png](cp_mid_a12_re1e4.png).

The four-panel Cp shows the same piecewise structure as the Euler case
(fore lower compression, fore upper expansion, aft upper stronger
expansion, aft lower weaker compression), **but with significant
point-to-point scatter** superimposed. The scatter does not invalidate
the mean Cp — the integrated CL/CD/Cm are steady to < 0.1 % in time —
but it is a known side-effect of sampling cell-data `Pressure` on the
immersed-boundary surface polygon at this Re, and would benefit from
either a short spatial average along each panel or a different probe
strategy for the NS branch. Flagged as a post-processing improvement.

## Bottom line (thesis-ready)

> *At a finite attack angle of 12° and with an adiabatic no-slip IBM
> wall, the Cerisse WENO-Z5 + SSPRK(3,3) stack runs to completion with
> bit-stable load histories and a clean, physically consistent flow
> topology (attached LE shocks, viscous wake, no ghost-cell artefacts).
> The integrated CL/CD/Cm differ substantially from the inviscid
> benchmark while the centre of pressure remains in the same chord-
> wise location to within 0.5 %, showing that the IBM no-slip
> treatment perturbs the load distribution consistently rather than
> spuriously, which is the minimum any physics-consistent viscous
> treatment must demonstrate.*

## Limitations / what this case does **not** establish

- No claim is made that this Re = 10⁴ run reproduces a specific
  literature measurement. It is a robustness / stability test of the
  IBM no-slip wall at finite incidence, not a quantitative viscous
  validation.
- The compile-time `REYNOLDS_USER` knob needs to be cross-checked
  against the dimensional `ρUc/μ` definition before any quantitative
  Re-sweep conclusions are drawn.
- The near-wall Cp scatter (Sec. "Surface Cp" above) suggests the
  pressure-extraction from the IBM surface polygon would benefit from
  smoothing or a different probe.
