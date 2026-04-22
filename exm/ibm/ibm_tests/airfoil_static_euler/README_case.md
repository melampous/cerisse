# airfoil_static_euler — canonical static double-wedge benchmark

Main Euler validation case for the IBM module. Intended to be written up
as a self-contained thesis subsection:

> The IBM can reproduce canonical oblique-shock / expansion-fan physics
> and integrated aerodynamic loads on a diamond wedge, across both
> thickness distributions and across the attached-to-detached-shock
> transition.

---

## 1. Convention

* Freestream fixed along `+x`, uniform M = 2 everywhere at t = 0.
* `ib.move = 0`. Each angle of attack is an independent static run.
* Positive angle of attack = geometry rotated **clockwise** about the
  origin; leading edge at (-0.05, 0) lifts upward.
* Rotation formula (with `(x_c, y_c) = (0, 0)`):

        x' =  cos(a) * x + sin(a) * y
        y' = -sin(a) * x + cos(a) * y

  Applied by `gen_rotated_airfoils.py`.

## 2. Geometries

Two thickness distributions, same chord (0.10 m) and same maximum
half-thickness (0.00437 m, ~8.7 %).

| name  | vertex order                                              | front θ | rear θ |
|---    |---                                                        |---      |---     |
| mid   | (-0.05, 0) (0, 0.00437) (0.05, 0) (0, -0.00437)           | 5.00 deg| 5.00 deg |
| x018  | (-0.05, 0) (-0.032, 0.00437) (0.05, 0) (-0.032, -0.00437) | 13.64 deg | 3.05 deg |

Both are written in the same vertex order (LE → Upper → TE → Lower → LE),
so the same panel-classification rule works after rotation.

## 3. AoA sweep

6 angles x 2 geometries = 12 .dat files / 12 inputs files:

```
0 deg, 4 deg, 8 deg, 12 deg, 16 deg, 18 deg
```

For the mid geometry the upper-front deflection reaches the Mach-2
detachment threshold at a = 18 deg, so that value sits right on the
transition. For the x018 geometry detachment occurs between a = 8 deg
and a = 12 deg, which is why the same sweep also carries the
"attached -> detached" narrative cleanly.

## 4. Outputs and observables

The four observables the write-up reports are:

* **beta** — leading-edge shock angle (upper and lower surfaces separately)
* **four-panel Cp** — fore-upper, aft-upper, aft-lower, fore-lower
* **CL, CD**
* **Cm and x_cp**

### 4.1 Panel classification

Panels are named by vertex order in the .dat file, not by global y sign.
After rotation the original "lower" vertex in x018 can land in
world-frame y > 0, which would break any y > 0 / y < 0 split.

The four panels are the four consecutive line segments in the
polygon traversal:

| segment        | panel name    |
|---             |---            |
| LE -> Upper    | `fore_upper`  |
| Upper -> TE    | `aft_upper`   |
| TE -> Lower    | `aft_lower`   |
| Lower -> LE    | `fore_lower`  |

Each surface sample is assigned to the panel whose line segment it lies
on; panel identity is decided from the .dat vertex order, never from
the sign of world-frame y. (After rotation, what used to be the
"Lower" vertex can land at world-frame y > 0; a y-sign split would
misclassify it.)

### 4.2 Moment / x_cp convention

* `Cm_LE`  — moment about the body-local leading edge.
* `Cm_c/4` — moment about the body-local quarter-chord
              (x_body = -0.025 m when chord = 0.1 m).
* `x_cp/c = -Cm_LE / CL`, reported only when `|CL| > eps` (eps = 1e-4
  of the expected post-shock `q_inf`).
* At a = 0 the geometry is symmetric and CL is at noise level.
  `x_cp` is NOT reported in that case; only `Cm_LE` and `Cm_c/4`
  (both at noise level) go into the table. This is recorded in the
  results README per-case, not silently hidden.

### 4.3 Time-window averaging

Loads are reported as time averages over a user-picked window:

1. First run with a fine plot interval (every ~100 steps).
2. Inspect `CL(t)`, `CD(t)`, `beta(t)` time histories.
3. Pick a window that starts after the leading-edge shock has
   established (`t > ~2 chord / u_inf ~ 3e-4 s`) and ends before any
   reflection from y = +/- 0.7 returns to the body region.
4. Report mean and standard deviation over that window.

## 5. Files

* `inputs_base`               — shared solver / AMR / BC template
* `inputs_{mid,x018}_a{00..18}` — 12 per-case inputs (only ib.filename
                                   and plot paths differ from inputs_base)
* `diamond_wedge_{mid,x018}.dat`           — two reference geometries
* `diamond_wedge_{mid,x018}_a{00..18}.dat` — 12 pre-rotated geometries
* `gen_rotated_airfoils.py`   — regenerates all rotated .dat files;
                                 embedded sanity checks on chord, signed
                                 area, vertex order
* `prob.h`                    — Euler + WENO-Z5 + slip wall + level-
                                 aware tagging; uniform freestream IC
* `GNUmakefile`               — CPU build; CLIP_MINTEMP = TRUE
