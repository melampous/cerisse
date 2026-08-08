# RZ near-axis positivity fix for the WENO-Z5 metric flux (Cerisse)

Author: automated debugging session, 2026-08-07.
Status: **verified** — robust (full 1.8 ms run, 0 aborts), 5th order preserved
(MMS bit-identical to unfixed), flow morphology consistent with the 2nd-order
baseline. Verified on AWS p5 / H100, CUDA 12.6, `weno_t<ReconScheme::WenoZ5>`.

This document is written to be handed to another engineer/AI. It gives the
symptom, the root cause, the exact code change (file + line anchors + snippets),
the two subtleties that make or break the fix, the verification, and the
dead-ends that were ruled out with data.

---

## 0. Scope — which scheme this touches

The fix is **only** in `src/rhs/Weno.h`, in the characteristic LLF WENO-Z5
**R-Z metric-flux path** (`weno_t<ReconScheme::WenoZ5, ...>`, the
`ComputeRzPairedPressure && reconstruct_radial_metric_flux` branch).

`grep -lE "reconstruct_radial_metric|rz_paired_pressure|radial_metric_scale"`
returns **nothing** for `src/rhs/Afd.h`, `src/rhs/Riemann.h` (HLLC),
`src/rhs/Rusanov.h`. Those flux schemes do **not** do the high-order r-weighted
metric reconstruction — they handle R-Z through the ordinary split + geometric
pressure source in `compute_rhs.cpp` (2nd order, already robust near the axis).
So **only the WENO-Z5 metric path has this defect and needs this fix.** If AFD
or HLLC are later given a high-order metric R-Z reconstruction, the same
principle (and the same two subtleties in §5) applies.

---

## 1. Symptom

2-D axisymmetric (`geometry.coord_sys = 1`) jet with the new 5th-order metric
R-Z reconstruction aborts on a non-positive state:

```
==================== STATE CHECK FAILED ====================
  level = 3  step = 2703  time = 7.97e-05
  min(rho) = 0.394   min(rho*e_int) = -87728   contains_nan = no
  argmin cell = (0,134)  x = 3.9e-05  y = 0.0105      <- r = dr/2 (axis-adjacent), z = 10.5 mm
amrex::Abort: non-finite / non-positive state (non-IBM check) !!!
```

- Always at the **axis-adjacent first radial cell** `i = 0` (`r = dr/2`).
- Triggered when a **contact/shock discontinuity crosses the axis** (here the
  starting-jet front: fast cold expanded flow meeting slow hot flow on the
  axis). The near-axis feature is a **contact** — density/temperature jump at
  ~constant pressure.
- The old 2nd-order scheme (ordinary split + geometric source) runs the same
  case cleanly — this is purely a defect of the metric reconstruction.

---

## 2. Root cause

The metric path reconstructs the **metric flux `r·F`** (r-weighted) with WENO,
rather than reconstructing `F` and multiplying by `r`. Near the axis:

1. The R-Z radial-momentum update of cell `i=0` is dominated by its single
   inner radial face and carries a geometric amplification factor
   `2/r_c = 2/dr`.
2. The WENO nonlinear weights are computed on **r-modulated samples**
   `r·(U ± F/α)`. At a discontinuity these weights misbehave, producing a
   spurious reconstruction.
3. The spurious reconstruction error, amplified by `2/dr`, is injected almost
   entirely into **radial momentum** → a spurious axis-ward velocity
   (`u_r ≈ -720 m/s` at the failing cell, vs `≈ -1.3 m/s` for the old scheme at
   the same cell/time) → kinetic energy overwhelms internal energy → `rho*e < 0`.

Confirmed by controlled experiments (see §7 dead-ends): it is a **positivity /
reconstruction** failure (not CFL — lower CFL only shrinks the overshoot, it
still crosses zero; not GPU; not a build issue).

---

## 3. The fix (one sentence)

On the innermost interior radial face(s), **when a density-jump sensor detects a
discontinuity**, replace the r-modulated WENO reconstruction with a **monotone
first-order metric LLF flux** built from the two adjacent cells — keeping the
r-weighting (metric-consistent) but using an **unweighted physical-jump
dissipation** so it is positivity-friendly and exact for uniform flow.

- **Dynamic (sensor-gated):** in smooth flow the sensor never fires → the fix is
  inert → the scheme is byte-for-byte the original 5th-order metric scheme. This
  is why the MMS order is untouched (§6).
- **Local:** only the K=1 innermost interior radial face; the rest of the domain
  is unchanged.

---

## 4. The math

For the face between the left cell `L = i-1` and right cell `R = i` at radius
`r_face`, with LLF wave speed `α = max(|u|+c)` over the stencil:

```
              0.5
  f_face[c] = ----- ( r_L · Fadv_L[c]  +  r_R · Fadv_R[c] )   <- metric-weighted flux average
             r_face
            - 0.5 · α · ( U_R[c] - U_L[c] )                   <- UNWEIGHTED physical-jump LLF dissipation
```

- `f_face` is stored in the per-area convention (metric-flux / r_face) that
  `compute_rhs.cpp` expects; it multiplies back by `r_face` for the metric
  divergence.
- `Fadv` is the physical flux with **pressure removed from the radial-momentum
  component** (`Fadv[UMX] = F[UMX] - p`); the pressure is carried separately.
- Pressure face value: `p_face = 0.5·(p_L + p_R)`, written to both
  `face_flux(UMX) = f_face[UMX] + p_face` and `radial_pressure_face_flux`,
  matching the paired-pressure convention in `compute_rhs.cpp` (which recovers
  the advective part as `radial_flux(UMX) - p` and adds the ordinary pressure
  gradient separately).

**Sanity check (uniform state, `u_r=0`):** `Fadv[URHO]=ρ·u_r=0`,
`U_R=U_L` ⇒ `f_face = 0`. The metric flux `r_face·f_face = 0`. Correct — no
spurious flux. (This is the property the naive version broke; see §5.2.)

---

## 5. Exact code change — `src/rhs/Weno.h`

Four insertions inside the face-flux `amrex::ParallelFor` device lambda of the
WENO reconstruction. Line numbers are for the pre-fix working tree; the anchor
line above each block is the unambiguous locator. Net change ≈ 35 lines of
logic. The full unified diff is in `analysis/rz_axis_fix/rz_axis_fix.patch`.

### 5.1 Block A — constants (after line ~1568)

Anchor (insert immediately **after** this line):
```cpp
            amrex::ignore_unused(ibm_markers, face_flux);
```
Insert:
```cpp
            // Near-axis metric-flux robustification constants (lambda-local
            // to avoid a constexpr-if device capture): first-order metric LLF
            // fallback on the innermost interior radial face. threshold < 0 =
            // always on that face; a positive value gates by density jump.
            constexpr int rz_axis_collar_faces = 1;
            constexpr Real rz_axis_shock_sensor_threshold = Real(0.003);
```
> **Must be declared inside the lambda.** Declaring these `constexpr` at the
> enclosing (loop) scope and referencing them inside the `if constexpr`
> paired-pressure branch triggers nvcc:
> *"An extended `__device__` lambda cannot first-capture variable in
> constexpr-if context."*

### 5.2 Block B — storage (after line ~1710)

Anchor (insert after the `sample_metric_scales` initialization loop, i.e. after
its closing `}`):
```cpp
              sample_metric_scales[sample_index] = Real(1.0);
            }
```
Insert:
```cpp
            // First-order metric fallback storage (immediate L/R cells).
            Real fo_cons[2][Closure::NCONS];
            Real fo_flux[2][Closure::NCONS];
            Real fo_density[2] = {Real(0.0), Real(0.0)};
            Real fo_pressure[2] = {Real(0.0), Real(0.0)};
            Real fo_radius[2] = {Real(0.0), Real(0.0)};
```

### 5.3 Block C — save the two adjacent cells (after line ~1742)

Anchor (inside the driver-sample loop, insert after this line):
```cpp
              sample_metric_scales[sample_index] = radial_metric_scale;
```
Insert:
```cpp
              if (reconstruct_radial_metric_flux &&
                  (sample_index == required_ghost_cells - 1 ||
                   sample_index == required_ghost_cells)) {
                const int fo_side = sample_index - (required_ghost_cells - 1);
                fo_density[fo_side] =
                    primitive_states(stencil_cell_index, Closure::QRHO);
                fo_pressure[fo_side] =
                    primitive_states(stencil_cell_index, Closure::QPRES);
                fo_radius[fo_side] = radial_origin +
                    (Real(stencil_cell_index[0]) + Real(0.5)) * radial_spacing;
                for (int fo_c = 0; fo_c < Closure::NCONS; ++fo_c) {
                  fo_cons[fo_side][fo_c] = conservative_state[fo_c];
                  fo_flux[fo_side][fo_c] = physical_flux[fo_c];
                }
              }
```
> `sample_index == required_ghost_cells - 1` is the left cell `i-1`,
> `== required_ghost_cells` is the right cell `i`. `conservative_state` and
> `physical_flux` are the per-sample arrays already computed a few lines above
> by `closure->prims2cons` / `closure->prims2flux`.

### 5.4 Block D — override with the first-order metric flux (after line ~1982)

Anchor — the end of the high-order metric UMX write, immediately **before** the
`return;` that closes the paired-pressure branch:
```cpp
                face_flux(face_index, Closure::UMX) =
                    radial_face > Real(0.0)
                        ? metric_advective_flux_per_radius + pressure_face
                        : Real(0.0);
                return;                                  // <- insert BEFORE this return
```
Insert (between the UMX write above and `return;`):
```cpp
                if (radial_face > Real(0.0) && face_index[0] >= 1 &&
                    face_index[0] <= rz_axis_collar_faces) {
                  const Real fo_dmin = amrex::min(fo_density[0], fo_density[1]);
                  const Real fo_djump = (fo_dmin > Real(0.0))
                      ? std::abs(fo_density[0] - fo_density[1]) / fo_dmin
                      : Real(1.0);
                  if (fo_djump > rz_axis_shock_sensor_threshold) {
                    const Real fo_inv_rface = Real(1.0) / radial_face;
                    const Real fo_alpha = llf_max_wave_speed;
                    Real fo_full[Closure::NCONS];
                    for (int fo_c = 0; fo_c < Closure::NCONS; ++fo_c) {
                      Real fo_adv_l = fo_flux[0][fo_c];
                      Real fo_adv_r = fo_flux[1][fo_c];
                      if (fo_c == Closure::UMX) {
                        fo_adv_l -= fo_pressure[0];
                        fo_adv_r -= fo_pressure[1];
                      }
                      // Metric-weighted flux average with an UN-weighted
                      // physical-jump LLF dissipation (vanishes for a uniform
                      // state; a metric-weighted dissipation would not).
                      fo_full[fo_c] = Real(0.5) * fo_inv_rface *
                          (fo_radius[0] * fo_adv_l + fo_radius[1] * fo_adv_r) -
                          Real(0.5) * fo_alpha *
                          (fo_cons[1][fo_c] - fo_cons[0][fo_c]);
                    }
                    const Real fo_pressure_face =
                        Real(0.5) * (fo_pressure[0] + fo_pressure[1]);
                    for (int fo_c = 0; fo_c < Closure::NCONS; ++fo_c) {
                      face_flux(face_index, fo_c) = fo_full[fo_c];
                    }
                    face_flux(face_index, Closure::UMX) =
                        fo_full[Closure::UMX] + fo_pressure_face;
                    radial_pressure_face_flux(face_index, 0) = fo_pressure_face;
                  }
                }
```
> Gate is `face_index[0] >= 1 && <= rz_axis_collar_faces` and `radial_face > 0`
> — i.e. the first `K` interior radial faces, **excluding** the axis face
> (`radial_face == 0`, which keeps its special metric slot). K=1 = the single
> face bounding cell `i=0`'s outer side, which is where the `2/dr` amplification
> bites.

### 5.5 Parameters (both `constexpr`, in Block A)

| name | value | meaning |
|---|---|---|
| `rz_axis_collar_faces` | `1` | number of innermost interior radial faces the fallback can act on |
| `rz_axis_shock_sensor_threshold` | `0.003` | relative density jump across the face that triggers first-order; `< 0` = always on |

---

## 6. The two subtleties that make or break it

Both were found the hard way (each initial version made the crash **worse** or
merely relocated it):

**6.1 The fallback must be metric-consistent (r-weighted), NOT a plain split.**
Simply switching the innermost face to the ordinary (non-metric) split makes the
face-flux convention inconsistent with its metric neighbours → a **new** crash
at the nozzle at step 8. First-order **metric** LLF (same r-weighting convention)
avoids this.

**6.2 The LLF dissipation must use the UNWEIGHTED physical jump `U_R - U_L`.**
The natural-looking form `α·(r_R·U_R - r_L·U_L)/r_face` does **not** vanish for a
uniform state (because `r_L ≠ r_R`), so it manufactures a large spurious mass
flux `≈ -0.5·α·ρ` even in uniform flow — draining density (`ρ → 0.12`) and
crashing **earlier** the wider the collar. The correct dissipation
`0.5·α·(U_R - U_L)` vanishes for uniform flow (§4 sanity check). This single sign
of the fix — where the r-weighting is and is not applied — is the crux.

---

## 7. Verification (all passed)

| check | method | result |
|---|---|---|
| **robustness** | full run `panda_mj1p42_L3`, `stop_time=1.8e-3` | step 9935 / 1.8 ms, **0 aborts** (was: abort at step 2703) |
| **order preserved** | RZ MMS (`exm/mms/eulerrz`, smooth), fixed vs unfixed, N=64 & N=128 | `max|Δρ| = 0.000` — **byte-identical** ⇒ sensor never fires in smooth flow ⇒ 5th order intact (unfixed metric MMS is 3.96e-11 @ N=128) |
| **morphology** | FIXED vs OLD baseline density field & centreline @ 1.8 ms | same shock-cell train / spacing / far-field; FIXED sharper (5th-order); first valley `z/D=0.36` identical |

**Verification logic for order:** because the fix only overwrites the face flux
inside `if (fo_djump > threshold)`, a byte-identical MMS result proves the branch
was never taken → the reconstruction is untouched in smooth flow → order is
exactly that of the base scheme.

---

## 8. Dead-ends ruled out (with data) — do NOT retry these

| approach | outcome |
|---|---|
| plain-split collar at inner face (static) | **worse** — inconsistent flux convention, crash at nozzle step 8 |
| density-sensor **plain-split** collar | crash step 323 (same class) |
| lower CFL (0.3→0.1) | still crashes at the same cell/time; overshoot shrinks `-87728 → -8434` but still `< 0` ⇒ **positivity, not CFL** |
| boost LLF `α` (×5) on the high-order flux | **worse**, overshoot `-87728 → -177010` (amplifies a bad reconstruction's dissipation) |
| first-order metric with **metric-weighted** dissipation | drains mass (`ρ→0.12`), crashes earlier the wider K — this was the §6.2 bug |
| static (always-on) first-order metric collar | robust but **drops MMS order to ~2nd** (memory: 2.35/1.98/1.93) — never ship static; use the sensor |
| density threshold `0.01` | fixed the first crash but a later contact re-crashed (`ρe=-4256`) via a temporal "sensor-off" accumulation gap; lowering to `0.003` fires consistently and still keeps MMS byte-identical |

---

## 9. Build & test recipe (H100 / CUDA 12.6)

- Build: `DIM=2 COMP=gnu CXXSTD=c++20 USE_MPI=TRUE USE_CUDA=TRUE CUDA_ARCH=90`,
  against a CCCL-guarded AMReX; **CUDA 12.6** (CUDA 13 breaks the AMReX build).
- Robustness case: `exm/underexpanded_jet/2d/panda_mj1p42_L3` (RZ L3),
  `stop_time=1.8e-3`; PASS = runs to end with no `STATE CHECK FAILED`.
- Order case: `exm/mms/eulerrz`, `input_dir/inputs{16,32,64,128}`; PASS = fixed
  run is byte-identical to the unfixed run (density) at N=64 and N=128.

---

## 10. Tuning notes for the next engineer

- If a stronger case still leaks a tiny residual: **lower the threshold** first
  (it stays MMS-identical as long as the smooth per-cell density jump at the
  finest MMS resolution is below it), or set `rz_axis_collar_faces = 2`.
- A cleaner, threshold-free sensor is the **WENO smoothness indicator τ**
  (already computed in the reconstruction) or a Jameson-type normalized second
  difference — both are resolution-aware, unlike a fixed relative-jump
  threshold. Not needed at 0.003 for the tested cases, but it removes the one
  remaining tunable.
