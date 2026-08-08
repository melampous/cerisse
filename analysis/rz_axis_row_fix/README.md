# RZ axis-row first-order fallback — fix and acceptance record

**Date:** 2026-08-08 · **Base:** `cerisse_local_0808a` = local working tree,
`src/rhs/Weno.h` md5 `f0ecbaf3267010cb1ebfd95f69dbf04c`
**Patched tree:** `/shared/cerisse_axfix`, md5 `d1444031bf1c7414e5ba60ad5941f6bc`
**Patch:** `rz_axis_row_fallback.patch` (112 lines, +39 net; `git apply --check` clean)

## Symptom

Uniform RZ grid, D_e/256 (3072 x 8192 = 25 165 824 cells), NPR = 15, from t = 0:

```
STATE CHECK FAILED
  level 0  step 15742  time 4.688312418e-4
  min(rho)      =  0.2398013086     (positive)
  min(rho*e_int) = -210218.2154     (negative)
  argmin cell   = (0, 663)   r = dr/2 (axis-adjacent),  z = 2.591 D_e
```

Deterministic and bit-reproducible across restart, rank count (8 -> 4) and I/O
settings. Uniform D_e/128 runs the same problem 0 -> 4 ms clean (71 222 steps).
AMR D_e/256 is clean **only because every such run restarted from a developed
checkpoint at 3.5 ms and never went through the startup transient**.

## Mechanism

1. The Mach disc drifts downstream at ~4.3 D_e/ms while the startup transient
   settles. The failure happens exactly when the disc sweeps through the
   axis-adjacent cell: in 42 steps that cell goes from post-shock subsonic
   (M = 0.05, rho e/E = 0.999) to pre-shock supersonic (M = 2.05, rho e/E = 0.46)
   with u_r reaching -64.7 m/s (elsewhere +-1..7 m/s).
2. The discrete RZ radial divergence (thesis Eq. 2.27) divides by r_j dr. At
   j = 0, r_0 = dr/2, so the divisor is dr^2/2: the outer radial flux is
   amplified by **2/dr relative to Cartesian**. Halving dr **doubles** it.
3. The existing near-axis collar covers only the first interior **radial** face
   (`face_index[0] == radial_axis_face_index + 1`). The disc sweeps in the
   **axial** direction, so the axial faces of that same cell fall back only via
   the generic normalized-density-curvature test, whose threshold (8.0e-2) is a
   grid-independent constant. At D_e/256 the sweeping disc never reaches it.

## Bisection evidence (runtime parameters, no rebuild)

| change | result |
|---|---|
| `cns.llf_fallback_density_jump` 3.0e-3 -> 1.0e-4 (30x wider) | still fails at step 15742; min(rho e) only -210218 -> -61469 |
| `cns.llf_fallback_density_curvature` 8.0e-2 -> 1.0e-2 | **clears the failure point** |
| `cfl` 0.3 -> 0.1 | fails **earlier** (t = 2.88e-4, cell (0,634)) — not a timestep problem |

## Fix

Tighten the curvature threshold **only on the axial faces of the axis-adjacent
cell row**, where the 2/dr amplification lives. Every other face keeps its
current threshold byte-for-byte, so interior accuracy is untouched. A global
threshold reduction was rejected: the axial smooth density jump is O(dz) with no
axis parity protection, so it would degrade the MMS order (the same trap as the
earlier static wide collar, which dropped MMS to second order).

New runtime control, default 1.0e-2:

```
cns.llf_fallback_axis_row_curvature = 1.0e-2
```

Manifest line now reads:

```
[Numerics] local_first_order_llf_fallback=enabled density_jump=0.003
  density_curvature=0.08 axis_row_curvature=0.01
  scope=cartesian_rz_all_fluid_shared_gp first_interior_rz_gate=density_jump_only
  rz_axis_row_axial_faces=tightened_curvature axis_face=excluded
```

## Acceptance

**A — the failure is cleared.** D_e/256 restarted from `chk14000` with default
thresholds:

```
last step 16972  (failure was at 15742)   t = 4.99976e-4 s (requested stop_time)
STATE CHECK FAILED = 0   amrex::Abort/Error = 0   AMReX finalized normally
```

**B — accuracy is provably untouched.** RZ MMS (WENO-Z5, N = 16/32/64/128),
patched vs unpatched:

```
mms_wenoz5  2.605336e-05  6.582529e-07  2.532432e-09  5.703811e-11  order 5.307 8.022 5.472
mms_axfix   2.605336e-05  6.582529e-07  2.532432e-09  5.703811e-11  order 5.307 8.022 5.472

L2 errors bit-identical to the unpatched solver: YES
```

The tightened axis-row threshold never fires in a smooth field, so the patched
solver reproduces the unpatched result bit for bit.

## Not covered

- The thesis methods chapter documents neither the first-order LLF fallback nor
  any positivity statement; both should be added.
- LLF-TENO5 in RZ uses the complete radius-weighted flux plus the geometric
  source rather than the paired-pressure operator, and has no RZ verification at
  all. Whether it shares this failure mode is untested.
- There is no RZ test with a shock anywhere in the verification suite; the
  Shu-Osher case is 1D Cartesian and the RZ MMS is smooth.
