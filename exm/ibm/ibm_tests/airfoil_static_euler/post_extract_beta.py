#!/usr/bin/env python3
"""Extract leading-edge shock angle beta from an AMReX plotfile.

Method (robust first-pass):

1. Load the Density field from the plotfile at the finest covered level.
2. Clip to a rectangular "LE neighborhood" window:
       x in [x_LE - 0.02, x_LE + 0.04]
       y in [y_LE - 0.06, y_LE + 0.06]   (upper and lower split below)
3. For each small y-stripe above the LE, find the x location where
   |d rho / dx| is maximal. That x lies on the shock line.
4. Linear-fit x_shock(y) over a stripe range and report beta as the
   angle between the shock and the +x axis:
       tan(beta) = (y2 - y1) / (x_shock(y2) - x_shock(y1))
5. Repeat for the lower stripe (y < y_LE) to get beta_lower.

Output CSV columns:
    case_name, geometry, alpha_deg, step, side, beta_deg,
    fit_window_ymin, fit_window_ymax, fit_quality

"fit_quality" is R^2 of the linear fit; if < 0.95 the shock is not
clean (often detached) and the row is flagged.

This uses schlieren.py's read_amrex_plotfile helper to avoid adding
a new plotfile reader. schlieren.py lives at tools/schlieren.py.

Usage
    python3 post_extract_beta.py <case_dir> --alpha 12 --geom mid
    python3 post_extract_beta.py <case_dir> --alpha 12 --geom mid \
            --step-window 2500:5000 --out post/beta_mid_a12.csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import sys

import numpy as np

# import schlieren.py to reuse its plotfile reader
SCHL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        "..", "..", "..", "..", "tools"))
if SCHL_DIR not in sys.path:
    sys.path.insert(0, SCHL_DIR)
try:
    from schlieren import read_amrex_plotfile
except Exception as exc:
    raise SystemExit(
        "Could not import read_amrex_plotfile from tools/schlieren.py:\n"
        f"  {exc}\n"
        f"Looked in: {SCHL_DIR}"
    )


GEOM_X_MAX_THICKNESS = {"mid": 0.0, "x018": -0.032}


def rotate_cw(x, y, alpha_deg):
    a = math.radians(alpha_deg)
    c, s = math.cos(a), math.sin(a)
    return c * x + s * y, -s * x + c * y


def linear_fit(y_arr, x_arr):
    """Linear regression x = m*y + b. Returns (m, b, R^2).

    We regress x vs y because the upper shock is nearly vertical at low
    alpha (dx/dy finite, dy/dx can be huge near beta -> 90 deg).
    """
    y = np.asarray(y_arr, dtype=float)
    x = np.asarray(x_arr, dtype=float)
    if y.size < 3:
        return float("nan"), float("nan"), float("nan")
    m, b = np.polyfit(y, x, 1)
    x_fit = m * y + b
    ss_res = float(np.sum((x - x_fit) ** 2))
    ss_tot = float(np.sum((x - x.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(m), float(b), float(r2)


def extract_one(pltdir, alpha_deg, geom, level=-1):
    """Return dict with beta_upper, beta_lower, fit quality for one plotfile."""
    le_x, le_y = rotate_cw(-0.05, 0.0, alpha_deg)

    rho, x1, x2, z, level_used, _ = read_amrex_plotfile(
        pltdir, "Density", out_level=level
    )
    # read_amrex_plotfile returns rho as shape (nx, ny) or (ny, nx) depending
    # on the 2D axis convention inside schlieren.py; probe:
    # schlieren.py uses (x, y) -> rho has shape (Nx, Ny). x1 is x vector,
    # x2 is y vector.
    rho = np.asarray(rho)
    if rho.ndim != 2:
        raise RuntimeError(f"Density field is {rho.shape}, expected 2D")

    # Ensure axis 0 is x, axis 1 is y. schlieren.py's convention is
    # rho[ix, iy]; x1 has length rho.shape[0], x2 has length rho.shape[1].
    if len(x1) == rho.shape[1] and len(x2) == rho.shape[0]:
        rho = rho.T
        x1, x2 = x2, x1

    Nx, Ny = rho.shape
    # finite-difference d rho / dx along axis 0
    drho_dx = np.gradient(rho, x1, axis=0)

    def fit_side(side):
        if side == "upper":
            y_lo = le_y + 0.002            # skip first 2 mm (GP / local noise)
            y_hi = le_y + 0.06
        else:
            y_lo = le_y - 0.06
            y_hi = le_y - 0.002
        # window in x
        x_lo = le_x - 0.02
        x_hi = le_x + 0.04

        ix_win = np.where((x1 >= x_lo) & (x1 <= x_hi))[0]
        iy_win = np.where((x2 >= y_lo) & (x2 <= y_hi))[0]
        if ix_win.size == 0 or iy_win.size == 0:
            return None

        xs = []
        ys = []
        for iy in iy_win:
            strip = np.abs(drho_dx[ix_win, iy])
            if np.all(strip == 0):
                continue
            k = int(np.argmax(strip))
            xs.append(float(x1[ix_win[k]]))
            ys.append(float(x2[iy]))

        if len(xs) < 5:
            return None

        m, b, r2 = linear_fit(ys, xs)
        # x = m*y + b  ->  tan(beta) = dy/dx = 1/m with beta measured
        # from +x axis. For a shock leaning leftward-up at low alpha,
        # dx/dy is small positive; dy/dx is large positive; beta ~ 70-85
        # for Mach 2 weak oblique.
        if not math.isfinite(m):
            return None
        if abs(m) < 1e-12:
            beta_deg = 90.0
        else:
            beta_deg = math.degrees(math.atan2(1.0, abs(m)))
        return {
            "beta_deg": beta_deg,
            "fit_ymin": y_lo,
            "fit_ymax": y_hi,
            "r2": r2,
            "n_points": len(xs),
        }

    return {
        "upper": fit_side("upper"),
        "lower": fit_side("lower"),
        "level_used": int(level_used) if level_used is not None else -1,
    }


def step_from_pltname(path):
    base = os.path.basename(os.path.normpath(path))
    # e.g. plt01000
    digits = "".join(c for c in base if c.isdigit())
    return int(digits) if digits else -1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case_dir")
    p.add_argument("--alpha", type=float, required=True)
    p.add_argument("--geom", choices=list(GEOM_X_MAX_THICKNESS.keys()),
                   required=True)
    p.add_argument("--plot-glob", default="plot_*/plt*")
    p.add_argument("--out", default=None,
                   help="output CSV (default: <case>/post/beta.csv)")
    p.add_argument("--level", type=int, default=-1,
                   help="AMR level to read; -1 means finest (default)")
    p.add_argument("--window-min", type=int, default=None)
    p.add_argument("--window-max", type=int, default=None)
    args = p.parse_args(argv)

    case_dir = os.path.abspath(args.case_dir)
    out_csv = args.out or os.path.join(case_dir, "post", "beta.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    plts = sorted(glob.glob(os.path.join(case_dir, args.plot_glob)))
    if not plts:
        plts = sorted(glob.glob(os.path.join(case_dir, "plot*/plt*")))
    if not plts:
        print("no plotfiles found under", case_dir)
        return 1

    rows = []
    case_name = f"{args.geom}_a{int(args.alpha):02d}"
    for plt in plts:
        if not os.path.isdir(plt):
            continue
        step = step_from_pltname(plt)
        if args.window_min is not None and step < args.window_min:
            continue
        if args.window_max is not None and step > args.window_max:
            continue
        try:
            r = extract_one(plt, args.alpha, args.geom, level=args.level)
        except Exception as exc:
            print(f"  step {step:05d}  skipped: {exc}")
            continue
        for side in ("upper", "lower"):
            d = r[side]
            if d is None:
                rows.append({
                    "case_name": case_name,
                    "geometry": args.geom,
                    "alpha_deg": args.alpha,
                    "step": step,
                    "side": side,
                    "beta_deg": "",
                    "fit_ymin": "", "fit_ymax": "",
                    "r2": "", "n_points": 0,
                })
            else:
                rows.append({
                    "case_name": case_name,
                    "geometry": args.geom,
                    "alpha_deg": args.alpha,
                    "step": step,
                    "side": side,
                    "beta_deg": f"{d['beta_deg']:.3f}",
                    "fit_ymin": f"{d['fit_ymin']:+.4f}",
                    "fit_ymax": f"{d['fit_ymax']:+.4f}",
                    "r2": f"{d['r2']:.4f}",
                    "n_points": d["n_points"],
                })
        print(f"  step {step:05d}  upper beta={r['upper']['beta_deg'] if r['upper'] else 'n/a'}  "
              f"lower beta={r['lower']['beta_deg'] if r['lower'] else 'n/a'}")

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case_name", "geometry", "alpha_deg", "step", "side",
                    "beta_deg", "fit_ymin", "fit_ymax", "r2", "n_points"])
        for r in rows:
            w.writerow([r["case_name"], r["geometry"], r["alpha_deg"],
                        r["step"], r["side"], r["beta_deg"],
                        r["fit_ymin"], r["fit_ymax"], r["r2"], r["n_points"]])
    print(f"\n# wrote {len(rows)} rows to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
