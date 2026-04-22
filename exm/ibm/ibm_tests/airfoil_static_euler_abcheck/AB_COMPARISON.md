# A/B equivalence check — rotate geometry vs rotate freestream

**Test case:** mid geometry at a relative attack angle of 12°, M∞ = 2.

| Variant | Geometry input | Freestream |
|---|---|---|
| **A** (baseline)   | `diamond_wedge_mid_a12.dat` (pre-rotated 12° CW) | u∞ = U cos 0°, v∞ = 0 |
| **B** (companion)  | `diamond_wedge_mid_a00.dat` (un-rotated)           | u∞ = U cos 12°, v∞ = U sin 12° |

All other settings (domain, AMR levels, CFL, schemes, run length) are
byte-identical.

## Integrated aerodynamic coefficients (tail-mean, steps 3000–5000)

The A/B case's raw force components are lab-frame (x_hat, y_hat); they
are rotated back by 12° to express lift/drag in the freestream-aligned
frame (the same frame the A case lives in):

      CD_aero = + cos(α) CD_lab + sin(α) CL_lab
      CL_aero = - sin(α) CD_lab + cos(α) CL_lab

| Quantity | A (rotate-geom) | B (rotate-FS) | |Δ|/A |
|---|---|---|---|
| CL       | +0.547311       | +0.543996     | 0.61 % |
| CD       | +0.138099       | +0.137534     | 0.41 % |
| Cm_LE    | +0.255836       | +0.255519     | 0.12 % |
| Cm_c/4   | +0.114820       | +0.115343     | 0.46 % |
| x_cp/c   | +0.467442       | +0.469708     | 0.48 % |

All five quantities match to within 0.6 % — well below the 1 %
threshold set a priori for engineering equivalence.

## Surface Cp comparison

See [cp_abcheck.png](cp_abcheck.png) (A/B) alongside
[../cp_panels/cp_mid_a12.png](../cp_panels/cp_mid_a12.png) (A). Each of
the four wedge panels reproduces to within visual precision:

| Panel      | A Cp    | B Cp    |
|---|---|---|
| fore upper | −0.08   | −0.08   |
| fore lower | +0.65   | +0.64   |
| aft  upper | −0.20   | −0.20   |
| aft  lower | +0.25   | +0.25   |

## Thesis-ready one-liner

> *At the same relative attack angle, the two implementations —
> rotating the immersed geometry vs tilting the freestream — produce
> CL, CD, Cm and surface Cp in agreement to <1 %. Since the rotate-
> geometry convention additionally exercises the IBM ghost-point
> reconstruction on obliquely oriented embedded boundaries, the
> agreement is a direct sanity check on that reconstruction.*

## Files

- Case dir: `exm/ibm/ibm_tests/airfoil_static_euler_abcheck/`
  - `prob.h` (tilted-freestream BCs, ALPHA_FS=12°)
  - `inputs_abcheck`, `run_aws_abcheck.sbatch`
  - `diamond_wedge_mid_a00.dat`
- AWS job: 73, completed 2026-04-22, 21 plotfiles, 50 vtp surfaces.
