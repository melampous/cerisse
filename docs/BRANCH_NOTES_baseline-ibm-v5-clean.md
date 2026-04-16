# Branch: `baseline/ibm-v5-clean`

**Purpose**: moving-body IBM debugging / FSI development baseline.

**NOT** a general-purpose IBM baseline. See scope note below.

---

## What this branch contains

### Algorithmic fixes (kept, active)
1. `src/ibm/ibm_solver.h` — `fixExposedCells` Pass 1/2/3 uses **single-donor
   selection with directional flow score**; Pass 1 additionally rejects
   donors with `p ≤ 0` or non-finite values.
2. `src/rhs/Weno.h` — `eflux_ibm` falls back to **1st-order LLF** at any
   face whose WENO stencil touches an IBM solid/ghost marker
   (`gl < ng || gr < ng`).

### Debug infrastructure (kept, disabled by default)
- `src/ibm/ibm_tripwire.h` — 4-stage pipeline scanner (T1 step entry,
  T2 post-rebuildIBM, T3 post-fixExposed, T4 post-regrid).
  Reports `min_rho`, `min_p`, `max_sig`, `n_bad` per level.
- Tripwire call sites in `advance.cpp` / `CNS.cpp`, window set to
  `(-1, -1)` so the scan never fires.
- `fixExposedCells` per-cell DIAG capture (`FixExposedDiagRec`),
  same `(-1, -1)` default window.

**To re-enable**: edit the `TRIPWIRE_LO`/`TRIPWIRE_HI` in
`src/tim/advance.cpp` and the inline `-1, -1` in `src/CNS.cpp`, plus
`DIAG_LO`/`DIAG_HI` inside `fixExposedCells`. Rebuild; diagnostic
output appears in stdout.

### Test case
- `exm/fsi/airfoil_pitchup_euler/` — 2D diamond-wedge airfoil, M=2,
  constant pitch after 720 μs settle. Includes `inputs_gpu` (level=1 AMR)
  and `test_level0/inputs_l0` (level=0 isolation).

---

## Scope & limitations

> **The `Weno.h` near_ib fallback affects ALL IBM/EB cases on this branch,
> not only the airfoil pitch-up case.** Static IBM cases will see slightly
> more near-wall dissipation vs `ibm_fsi` HEAD.
>
> This branch is therefore **intentionally dedicated to moving-FSI debugging**.
> Do not merge into `ibm_fsi` without a regression sweep on the static-IBM
> tests or without gating the near_ib fallback behind a runtime switch.

### Verified behaviour
- Static IBM: unchanged (no regressions observed on isolated tests).
- Moving body, `max_level=0`: stable for 5000 steps, reaching ~153° pitch.
- Moving body, `max_level=1` AMR: first-failure at step 2107 (original)
  is resolved; secondary instability now at step ~2397 (~24° AoA);
  **long-time AMR solution is not fully reliable**.

### Not verified
- `max_level ≥ 2`
- Viscous (Navier–Stokes)
- Other FSI geometries (sphere, etc.)

---

## Related commits / tags

| Ref | Commit | Purpose |
|---|---|---|
| `v5-frozen-debug` (tag) | `dc70f8a40` | Full debug snapshot (tripwires + DIAG windows set to real ranges) |
| `baseline/ibm-v5-clean` (this branch) | `d0034f43f` | Debug disabled; algorithmic fixes retained |

---

## Documentation
- `docs/debug_moving_ibm_v5_summary.md` — full debug narrative, timeline,
  limitations, and future-work priorities.
