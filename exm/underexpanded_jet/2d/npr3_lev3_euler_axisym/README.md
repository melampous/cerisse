# 2d/npr3_lev3_euler_axisym

Underexpanded jet case — auto-generated from `data_organization/manifest.py` on 2026-04-22.

## Identity
- **host**: `local`
- **path**: `/home/qiaoj/testcerisse/cerisse/exm/underexpanded_jet/2d/npr3_lev3_euler_axisym`
- **dim**: 2D
- **coord_sys**: 1 (RZ axisymmetric)
- **NPR_0** (stagnation/ambient): 3.0
- **amr.max_level**: 3
- **physics**: EULER  (inviscid Euler)

## Flow setup (canonical, may be overridden by `prob.h`)
- Jet exit: `D_e = 10 mm`, `M_jet = 1.0` (sonic exit), `p_jet = NPR_0 / 1.89293 * p_amb` (isentropic from stagnation).
- Ambient: `p_amb = 101325 Pa`, `T_amb = 300 K`, air (γ=1.4).
- `pe/pa = NPR_0 / (1+0.2*M^2)^3.5 = 1.5848` (underexpanded).

## Grid setup
- Domain: `[0, 60] mm × [0, 240] mm` (r × z, 6De × 24De)
- L0 = 96×384 (`dx=0.625 mm`); levels: L1=0.313, L2=0.156, L3=0.078 mm

## Algorithm
- Time stepping: `RK3` (`cns.order_rk = 3`)
- Riemann + reconstruction: as in `prob.h` / `Make.CNS` (typically WENO5-Z + HLLC).
- Viscous terms: OFF (Euler).
- BCs: `lo = 3 1` (axis on r=0, inflow on z=0), `hi = 2 2` (outflow).

## Outputs
- `plot/pltXXXXX/` plotfiles, `plot/chkXXXXX/` checkpoints.
- See `inputs` for `amr.plot_int`, `amr.check_int` cadence.

## Notes
- Local source-only mirror of cx3/aws case (prob.h, inputs, GNUmakefile).

## Companion runs (same NPR, different host/level/physics)
- `cx3:2d/npr3_lev3_euler_axisym` (euler, lev=3)
- `cx3:2d/npr3_lev3_euler_axisym_gridconv_coarse` (euler, lev=3)
- `cx3:2d/npr3_lev3_euler_axisym_gridconv_medium` (euler, lev=3)
- `cx3:2d/npr3_lev3_euler_axisym_gridconv_fine` (euler, lev=3)
- `aws:2d/npr3_lev3_euler_axisym` (euler, lev=3)
