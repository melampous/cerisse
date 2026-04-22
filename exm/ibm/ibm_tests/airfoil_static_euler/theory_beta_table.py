#!/usr/bin/env python3
"""Generate the theoretical leading-edge shock-angle table for the
mid / x018 geometries across the 6 AoA values used in the sweep.

For each (geometry, alpha) pair we compute the effective flow-deflection
angle on the upper and lower front face and then, from the oblique
theta-beta-M relation

    tan(theta) = 2 cot(beta) *
                 (M^2 sin^2(beta) - 1) / (M^2 (gamma + cos(2 beta)) + 2)

we solve numerically for the weak-solution beta. Detachment happens
at theta > theta_max(M).

Classification of the front face:

    theta_eff > 0  -> compression (weak oblique shock if attached)
    theta_eff < 0  -> expansion (Prandtl-Meyer fan, no shock)

Output is written both to stdout and to `theory_beta.csv`:

    geometry, alpha_deg, side, theta_eff_deg, beta_weak_deg, state

where `state` is one of:
    attached    (weak oblique shock, beta reported)
    marginal    (theta_eff within 0.5 deg of theta_max)
    detached    (no real weak-solution beta)
    expansion   (theta_eff < 0, Prandtl-Meyer)
    zero        (theta_eff ~ 0, no wave)
"""

from __future__ import annotations

import csv
import math
import os
import sys

M_INF = 2.0
GAMMA = 1.4

# Geometry-level front half-angles (in the *body frame*, before rotation).
GEOM_FRONT_HALF_ANGLE = {
    "mid":  math.degrees(math.atan2(0.00437, 0.050)),   # = 5.00 deg (front half = rear half)
    "x018": math.degrees(math.atan2(0.00437, 0.018)),   # = 13.64 deg (LE to x_maxT = 0.018 m)
}

ALPHAS = [0.0, 4.0, 8.0, 12.0, 16.0, 18.0]


def theta_of_beta(beta_deg, M=M_INF, g=GAMMA):
    """Oblique shock: deflection angle theta as a function of wave angle beta."""
    b = math.radians(beta_deg)
    if math.sin(b) <= 1.0 / M:
        return 0.0  # below Mach wave
    num = M * M * math.sin(b) ** 2 - 1.0
    den = M * M * (g + math.cos(2 * b)) + 2.0
    return math.degrees(math.atan2(2.0 / math.tan(b) * num, den))


def theta_max(M=M_INF, g=GAMMA):
    """Maximum deflection angle (beta at ~64 deg for M=2)."""
    best = 0.0
    # sweep beta in 0.05 deg increments between Mach-angle and 90
    mach_ang = math.degrees(math.asin(1.0 / M))
    b = mach_ang
    while b < 90.0:
        t = theta_of_beta(b, M, g)
        if t > best:
            best = t
        b += 0.05
    return best


def beta_weak(theta_deg, M=M_INF, g=GAMMA, tol=1e-4):
    """Invert theta(beta); return the WEAK-solution beta in degrees
    (the smaller root; between Mach angle and theta_max's beta).
    Returns None if no real solution (detached)."""
    if theta_deg <= 0:
        return None
    tmax = theta_max(M, g)
    if theta_deg > tmax:
        return None
    # weak root lies between Mach-angle and beta-at-theta_max (~64 deg for M=2)
    mach_ang = math.degrees(math.asin(1.0 / M))
    lo, hi = mach_ang + 1e-3, 64.06  # strict upper bound for M=2
    # monotonic region; bisection on theta(beta)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        t_mid = theta_of_beta(mid, M, g)
        if t_mid < theta_deg:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def classify(theta_eff, tmax):
    """Return state label."""
    if abs(theta_eff) < 1e-3:
        return "zero"
    if theta_eff < 0:
        return "expansion"
    if theta_eff > tmax:
        return "detached"
    if tmax - theta_eff < 0.5:
        return "marginal"
    return "attached"


def main():
    tmax = theta_max()
    rows = []
    for geom, th_w in GEOM_FRONT_HALF_ANGLE.items():
        for a in ALPHAS:
            for side, sign in (("upper", +1.0), ("lower", -1.0)):
                # effective front-face deflection when flow meets the
                # leading edge: upper face turns flow by (th_w + alpha),
                # lower face by (th_w - alpha). Signs track whether the
                # rotated surface is still compression or has become
                # expansion.
                theta_eff = th_w + sign * a
                # for 'lower', sign=-1: theta_eff = th_w - a
                # for 'upper', sign=+1: theta_eff = th_w + a
                # (the "sign" name is about which face, not about theta itself)
                if side == "lower":
                    theta_eff = th_w - a
                else:
                    theta_eff = th_w + a
                b_weak = beta_weak(theta_eff) if theta_eff > 0 else None
                state = classify(theta_eff, tmax)
                rows.append({
                    "geometry": geom,
                    "alpha_deg": a,
                    "side": side,
                    "theta_eff_deg": theta_eff,
                    "beta_weak_deg": b_weak if b_weak is not None else "",
                    "state": state,
                })

    # stdout table
    print(f"# Mach = {M_INF}, gamma = {GAMMA}")
    print(f"# theta_max = {tmax:.3f} deg (detachment threshold for M={M_INF})")
    print()
    print(f"{'geom':<6} {'alpha':>6} {'side':<6} "
          f"{'theta_eff':>11} {'beta_weak':>11} {'state':<12}")
    for r in rows:
        beta = (f"{r['beta_weak_deg']:.2f}"
                if isinstance(r['beta_weak_deg'], float) else "-")
        print(f"{r['geometry']:<6} {r['alpha_deg']:>6.0f} {r['side']:<6} "
              f"{r['theta_eff_deg']:>10.2f}  {beta:>10}  {r['state']:<12}")

    out_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "theory_beta.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["geometry", "alpha_deg", "side",
                    "theta_eff_deg", "beta_weak_deg", "state"])
        for r in rows:
            w.writerow([r["geometry"], r["alpha_deg"], r["side"],
                        f"{r['theta_eff_deg']:.3f}",
                        (f"{r['beta_weak_deg']:.3f}"
                         if isinstance(r['beta_weak_deg'], float) else ""),
                        r["state"]])
    print(f"\n# wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
