# airfoil_static — viscous no-slip method-oriented benchmark (NS secondary)

Secondary, method-oriented companion to `airfoil_static_euler`.

The main IBM validation against canonical double-wedge shock-expansion
physics is done in `airfoil_static_euler/` (Euler + slip wall + uniform
M = 2). This case adds a viscous / no-slip layer on **exactly the same
geometry, same freestream, same domain**, and sweeps the Reynolds number
as a method stress parameter.

## Scope

> The Re sweep is a **robustness test of the IBM wall treatment**, not a
> quantitative aerodynamic reference. It exists to show that the
> ghost-point adiabatic no-slip wall stays stable, produces a smooth
> surface pressure plateau, and does not introduce spurious oscillation
> near the leading edge as the boundary layer is made progressively
> thinner.

## Configuration

Identical to `airfoil_static_euler` except:

| item                  | value                                                 |
|---                    |---                                                    |
| inviscid flux         | WENO-Z5 (unchanged)                                   |
| diffusive flux        | `viscous_t<methodparm_t, ProbClosures>` (enabled)     |
| transport model       | `transport_const_t<methodparm_t>` (constant μ, λ)     |
| wall model            | `ibm_adiabatic_noslip_wall_t`                         |
| μ                     | 1 / Re (per-case override in `methodparm_t`)          |
| Pr                    | 0.7 (for λ = μ c_p / Pr)                              |
| geometry              | pre-rotated .dat from airfoil_static_euler directory  |

## Reynolds sweep

| Re    | Role                                                             |
|---    |---                                                               |
| 1e4   | Method stress case. Boundary layer is thick relative to chord;    |
|       | the IBM wall sits on about dx/delta ~ O(1).                      |
| 1e5   | Intermediate, still fully laminar in the solver.                 |
| 1e6   | High-Re robustness check (see disclaimer below).                 |

## High-Reynolds disclaimer

> The `Re = 1e6` cases are included as a high-Reynolds-number robustness
> check of the immersed-boundary wall treatment. Given the present grid
> and the absence of a transition or turbulence model, these solutions
> are interpreted qualitatively and are not used as quantitative
> skin-friction benchmarks. At Re = 1e6 the laminar Blasius boundary
> layer thickness δ / c ~ 5 / √Re ~ 5e-3. With the present finest cell
> size (about chord / 315) the boundary layer is represented on fewer
> than 2 cells normal to the wall over most of the chord, which is
> sufficient to detect stability or breakdown of the IBM wall
> treatment but not sufficient to resolve skin friction accurately.

## Why not Sutherland

A fair question is whether `transport_suth_t` (μ(T) Sutherland law)
should replace the constant-μ model here. It was considered and
**declined** for this Re sweep. Two reasons:

1. The purpose of this case is a **Re sweep** as a method-stress axis.
   Sutherland locks μ to the physical-air reference (~1.8e-5 kg/m/s at
   300 K), which removes Re as a user-controllable parameter. To sweep
   Re one would then have to distort ρ∞ or u∞, which would break the
   same-freestream equivalence with the Euler companion case.
2. The expected temperature excursion across a Mach-2 weak oblique
   shock is small (T/T∞ ~ 1.15–1.3), so Sutherland would change μ by
   only 20–40% relative to its freestream value. That magnitude of
   variation is below the method-stress signal this case is trying to
   expose.

Sutherland is therefore kept as an option for a future, separately
scoped **physical flight-regime case** (for example a direct
comparison against Bertram and McCauley double-wedge data), not for
this sweep.

## Files

* `inputs_re1e4`, `inputs_re1e5`, `inputs_re1e6` — per-Re inputs;
  geometry is picked via `ib.filename = ../airfoil_static_euler/diamond_wedge_<name>.dat`
  so that both cases share the same single source of rotated .dat files.
* `prob.h` — unchanged from the historical airfoil_static, except that
  the scheme stack uses `viscous_t` + `transport_const_t` and the wall
  is `ibm_adiabatic_noslip_wall_t`.

## First-run recommendation

Pick two sanity points before running the whole sweep:

* `mid_a00` + Re=1e4
* `mid_a00` + Re=1e6

Compare:

* Surface pressure along the four panels (stable plateau, no ringing).
* `CL(t) / CD(t)` time history — look for overshoot-and-settle pattern
  only, no unbounded drift.
* Wall temperature field near the leading edge — bounded above by the
  adiabatic stagnation estimate.

Only after those two are clean does it pay off to run the remaining
a × Re matrix.
