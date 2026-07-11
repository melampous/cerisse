# 3d/npr3_lev4_ns

Underexpanded jet case — auto-generated from `data_organization/manifest.py` on 2026-04-22.

## Identity
- **host**: `local`
- **path**: `/home/qiaoj/testcerisse/cerisse/exm/underexpanded_jet/3d/npr3_lev4_ns`
- **dim**: 3D
- **coord_sys**: 0 (Cartesian)
- **NPR_0** (stagnation/ambient): 3.0
- **amr.max_level**: 4
- **physics**: NS  (compressible Navier-Stokes)

## Flow setup (canonical, may be overridden by `prob.h`)
- Jet exit: `D_e = 10 mm`, `M_jet = 1.0` (sonic exit), `p_jet = NPR_0 / 1.89293 * p_amb` (isentropic from stagnation).
- Ambient: `p_amb = 101325 Pa`, `T_amb = 300 K`, air (γ=1.4).
- `pe/pa = NPR_0 / (1+0.2*M^2)^3.5 = 1.5848` (underexpanded).

## Grid setup
- Domain: `[0, 240] mm × [-60, 60] mm × [-60, 60] mm`  (24De × 12De × 12De)
- L0 = 384×192×192 (`dx=0.625 mm`)
- Finest level L4: `dx ≈ 39 µm`

## Algorithm
- Time stepping: `RK3` (`cns.order_rk = 3`)
- Riemann + reconstruction: as in `prob.h` / `Make.CNS` (typically WENO5-Z + HLLC).
- Viscous terms: ON (compressible NS).
- BCs: `lo = 1 2 2` (inflow on x=0, outflow lateral), `hi = 2 2 2` (outflow).

## Outputs
- `plot/pltXXXXX/` plotfiles, `plot/chkXXXXX/` checkpoints.
- See `inputs` for `amr.plot_int`, `amr.check_int` cadence.

## Notes
- Local source-only mirror of canonical Powell-Tam (Case 2) production case on cx3.

## Companion runs (same NPR, different host/level/physics)
- `cx3:3d/npr3_lev2_euler` (euler, lev=2)
- `cx3:3d/npr3_lev4_ns` (ns, lev=4)
- `aws:3d/npr3_lev4_euler` (euler, lev=4)
- `aws:3d/npr3_lev3_ns` (ns, lev=3)
