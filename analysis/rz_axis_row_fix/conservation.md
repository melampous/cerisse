# Discrete conservation of the LLF-WENO-Z5 operator — measured

**Date:** 2026-08-08 · **Solvers:** `cerisse_local_0808a` (unpatched) and
`cerisse_axfix` (axis-row fix) · **Case:** `/shared/conserv_test`

## Test

Closed-box axisymmetric blast. Domain r,z in [0,1]. Every boundary reflecting
(r=0 Symmetry code 3, all others SlipWall code 4 — pure AMReX mirror conditions,
no `bcnormal` path). Euler fluxes (`no_diffusive_t`), no source term. A
high-pressure blob (rho 10, p 100 vs 1, 1) centred **on the axis** drives strong
shocks repeatedly through the axis-adjacent cell row, so the first-order fallback
under test does fire.

With no flux through any boundary, the volume-weighted sums of rho and rho*E must
be constant for any conservative flux-difference scheme, whatever flux function
each individual face uses. `CNS::printTotal()` supplies them, using
`geom.GetVolume()` (cylindrical cell volumes in RZ).

Radial momentum is excluded: in RZ it has a genuine physical source (the
geometric term), so it is not a conserved quantity. Axial momentum is also not
usable here — the initial condition is symmetric about z = 0.5, so it is zero by
symmetry at all times regardless of the scheme.

## Result — matched physical time T = 0.15

| geometry | N | d(rho)/rho_0 | d(rho E)/(rho E)_0 | steps |
|---|---|---|---|---|
| RZ | 64 | 2.8877e-04 | 1.8522e-03 | 173 |
| RZ | 128 | 7.9321e-05 | 5.1365e-04 | 362 |
| RZ | 256 | 1.7263e-05 | 1.0775e-04 | 737 |
| Cartesian | 64 | 5.4157e-15 | 7.2142e-15 | 192 |
| Cartesian | 128 | 9.7668e-15 | 1.4669e-14 | 401 |
| Cartesian | 256 | 1.8029e-14 | 2.8924e-14 | 818 |

Observed orders: RZ density **1.86, 2.20**; RZ energy **1.85, 2.25**.
Cartesian: −0.85 to −1.02, i.e. the drift grows in proportion to the step count
(192 → 401 → 818 while the drift doubles twice) at roughly 2e-17 per step — this
is floating-point round-off accumulation, not a discretisation error.

## Conclusions

1. **The Cartesian operator is discretely conservative to machine precision**,
   including when the first-order fallback fires. This is expected from the
   structure: the fallback writes one value to `face_flux(face_index, ·)` and
   returns, so the telescoping identity behind Eq. 2.22-2.23 is untouched —
   changing which flux function a face uses cannot break conservation.

2. **The RZ operator is not discretely conservative.** The defect is a consistent
   **O(dr^2)** discretisation error that converges under refinement: halving dr
   reduces it fourfold. At N = 128 over T = 0.15 the mass defect is 0.008 % and
   the energy defect 0.05 %.

3. **Mechanism** — the thesis states the design choice but does not quantify its
   consequence. §2.2.3: "Although the physical radius-weighted flux H_r = r F_r
   vanishes at the axis, the auxiliary half-index flux H_hat_{r,-1/2} **need not
   vanish** ... It is constructed so that the numerical flux difference provides
   a high-order approximation to the radial divergence at the adjacent cell
   centre." The cylindrical volume-weighted sum telescopes to
   `H_hat_{N+1/2} - H_hat_{-1/2}`; with a closed outer wall the first term is
   zero, so a nonzero auxiliary axis flux leaks mass and energy every step.
   §2.3.5 records the same choice for the skew operator: "This axis-face
   auxiliary flux is not forced to zero."

4. **The axis-row fix is not the cause.** Unpatched vs patched, same case:
   density 8.374e-05 vs 8.021e-05, energy 5.245e-04 vs 5.161e-04.

## Usage boundary

- Mass/energy budget analyses in RZ (mixing efficiency, total-enthalpy checks)
  must quantify this defect first; it cannot be assumed zero. Refinement reduces
  it as dr^2 but never to round-off.
- Forcing `H_hat_{-1/2} = 0` would restore exact conservation at the cost of the
  high-order divergence approximation at the axis-adjacent cell. The thesis
  already reports the first-ring radial-momentum residual converging at order
  0.998, so that accuracy has little margin. **This is a conservation-versus-
  near-axis-accuracy trade, not a defect to patch blindly.**
- RZ combined with a shock had no verification case anywhere in the suite
  (Shu-Osher is 1D Cartesian, the RZ MMS is smooth). This closed-box blast is the
  first, and it yields the scaling law above; it is worth keeping as a regression
  test.
