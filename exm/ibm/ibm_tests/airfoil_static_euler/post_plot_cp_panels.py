#!/usr/bin/env python3
"""Plot the four-panel surface Cp distribution for a single case.

Reads the per-snapshot CSVs produced by post_extract_cp.py and plots
    Cp  vs  normalized chordwise body-frame x,   one curve per panel:
        fore_upper (solid, red)
        aft_upper  (solid, orange)
        aft_lower  (dashed, blue)
        fore_lower (dashed, cyan)

Default picks the last snapshot CSV in the `post/` directory. Override
with --step <NNNNN>. Optional --theory flag overlays the two-segment
shock-expansion analytic Cp on the upper and lower front faces, using
theory_beta.csv (needs to exist in the same directory).

Usage
    python3 post_plot_cp_panels.py <case_dir> --alpha 12 --geom mid
    python3 post_plot_cp_panels.py <case_dir> --alpha 12 --geom mid \
            --step 5000 --theory -o cp_mid_a12.png
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

P_INF = 101325.0
RHO_INF = 1.1768
U_INF = 694.43
Q_INF = 0.5 * RHO_INF * U_INF * U_INF
GAMMA = 1.4
M_INF = 2.0

PANEL_STYLE = {
    "fore_upper": dict(color="C3", ls="-",  marker="o", label="fore upper"),
    "aft_upper":  dict(color="C1", ls="-",  marker="s", label="aft upper"),
    "aft_lower":  dict(color="C0", ls="--", marker="v", label="aft lower"),
    "fore_lower": dict(color="C9", ls="--", marker="^", label="fore lower"),
}


def read_cp_csv(path):
    """Return dict: panel -> list of (x_body, Cp)."""
    panels = {k: [] for k in PANEL_STYLE}
    with open(path) as f:
        for row in csv.DictReader(f):
            p = row["panel"]
            panels[p].append((float(row["x_body"]), float(row["Cp"])))
    for p in panels:
        panels[p].sort()
    return panels


def cp_after_oblique_shock(theta_deg, M1=M_INF, g=GAMMA):
    """Return Cp immediately after a weak oblique shock deflecting
    flow by theta_deg. None if detached."""
    if theta_deg <= 0:
        return None
    # bisect beta
    from math import radians, sin, cos, tan, atan2, degrees
    def theta_of_beta(b_deg):
        b = radians(b_deg)
        if sin(b) <= 1.0 / M1: return 0.0
        return degrees(atan2(
            2.0 / tan(b) * (M1 * M1 * sin(b) ** 2 - 1.0),
            M1 * M1 * (g + cos(2 * b)) + 2.0))
    mach_ang = degrees(math.asin(1.0 / M1))
    lo, hi = mach_ang + 1e-3, 64.06
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if theta_of_beta(mid) < theta_deg:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-5: break
    beta_deg = 0.5 * (lo + hi)
    b = math.radians(beta_deg)
    p2_p1 = 1.0 + 2.0 * g / (g + 1.0) * (M1 * M1 * math.sin(b) ** 2 - 1.0)
    # Cp = (p2 - p1) / q_inf = (p2/p1 - 1) * p1 / (0.5 gamma p1 M1^2)
    return (p2_p1 - 1.0) * 2.0 / (g * M1 * M1)


def cp_after_expansion(theta_deg, M1=M_INF, g=GAMMA):
    """Cp after a Prandtl-Meyer expansion through |theta| degrees."""
    if theta_deg >= 0:
        return None
    theta = abs(theta_deg)
    # Prandtl-Meyer function
    def nu(M):
        if M < 1.0:
            return 0.0
        k = math.sqrt((g + 1.0) / (g - 1.0))
        return math.degrees(
            k * math.atan(math.sqrt((g - 1.0) / (g + 1.0) * (M * M - 1.0)))
            - math.atan(math.sqrt(M * M - 1.0))
        )
    nu1 = nu(M1)
    nu2 = nu1 + theta
    # solve for M2 via bisection
    lo, hi = M1, 30.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if nu(mid) < nu2:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-5: break
    M2 = 0.5 * (lo + hi)
    # isentropic: p2/p1 = ((1 + (g-1)/2 M1^2)/(1 + (g-1)/2 M2^2))^(g/(g-1))
    p2_p1 = ((1 + 0.5 * (g - 1) * M1 * M1) /
             (1 + 0.5 * (g - 1) * M2 * M2)) ** (g / (g - 1))
    return (p2_p1 - 1.0) * 2.0 / (g * M1 * M1)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case_dir")
    p.add_argument("--alpha", type=float, required=True)
    p.add_argument("--geom", choices=["mid", "x018"], required=True)
    p.add_argument("--step", type=int, default=None,
                   help="step number; default = latest CSV")
    p.add_argument("--post-dir", default=None)
    p.add_argument("--theory", action="store_true",
                   help="overlay analytic Cp from shock-expansion")
    p.add_argument("-o", "--out", default=None)
    args = p.parse_args(argv)

    post_dir = (args.post_dir or
                os.path.join(os.path.abspath(args.case_dir), "post"))
    csvs = sorted(glob.glob(os.path.join(post_dir, "cp_step*.csv")))
    if not csvs:
        print("no cp_step*.csv under", post_dir); return 1

    if args.step is None:
        chosen = csvs[-1]
    else:
        chosen = [c for c in csvs
                  if int(os.path.basename(c).split("step")[1].split(".")[0])
                     == args.step]
        if not chosen:
            print(f"no CSV at step {args.step}"); return 1
        chosen = chosen[0]

    panels = read_cp_csv(chosen)

    fig, ax = plt.subplots(figsize=(7, 5))
    for name, pts in panels.items():
        if not pts: continue
        xs = np.array([p[0] for p in pts])
        ys = np.array([p[1] for p in pts])
        ax.plot(xs, ys, **PANEL_STYLE[name])

    # theory overlay
    #
    # Convention (positive alpha = nose-up): the upper front face sees a
    # reduced compression or an expansion, the lower front face sees an
    # increased compression. So:
    #   upper θ_eff = θ_w − α   (Prandtl-Meyer when α > θ_w)
    #   lower θ_eff = θ_w + α   (always compressive for α ≥ 0)
    if args.theory:
        th_w = (math.degrees(math.atan2(0.00437, 0.050))
                if args.geom == "mid"
                else math.degrees(math.atan2(0.00437, 0.018)))
        # UPPER front face
        theta_upper = th_w - args.alpha
        if theta_upper > 0:
            cp_fu = cp_after_oblique_shock(theta_upper)
            lbl = f"theory fore_upper (θ={theta_upper:.1f}°)"
        else:
            cp_fu = cp_after_expansion(theta_upper)
            lbl = f"theory fore_upper (exp {abs(theta_upper):.1f}°)"
        if cp_fu is not None:
            ax.axhline(cp_fu, color="C3", ls=":", lw=1.2, label=lbl)
        # LOWER front face
        theta_lower = th_w + args.alpha
        if theta_lower > 0:
            cp_fl = cp_after_oblique_shock(theta_lower)
            lbl = f"theory fore_lower (θ={theta_lower:.1f}°)"
        else:
            cp_fl = cp_after_expansion(theta_lower)
            lbl = f"theory fore_lower (exp {abs(theta_lower):.1f}°)"
        if cp_fl is not None:
            ax.axhline(cp_fl, color="C9", ls=":", lw=1.2, label=lbl)

    ax.set_xlabel("x_body (chord-wise)")
    ax.set_ylabel("Cp")
    ax.set_title(f"Four-panel Cp  {args.geom}  α={args.alpha}°  "
                 f"(from {os.path.basename(chosen)})")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(loc="best", fontsize=8)
    out = args.out or os.path.join(post_dir, f"cp_panels_{args.geom}_a{int(args.alpha):02d}.png")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print("Saved", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
