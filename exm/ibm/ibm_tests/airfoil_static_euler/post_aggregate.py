#!/usr/bin/env python3
"""Aggregate a case matrix into alpha-sweep CL/CD/Cm/beta curves.

Walks through every per-case `post/integrals.csv` (written by
post_integrals.py) and `post/beta.csv` (written by post_extract_beta.py)
in a directory layout like

    <case_dir>/post/integrals.csv
    <case_dir>/post/beta.csv

and produces one summary CSV and a set of comparison plots:

    post_summary/summary_integrals.csv
        case, geom, alpha_deg, step, CL_mean, CD_mean, CmLE_mean, Cmc4_mean,
        xcp_mean, n_steps_averaged, window_min, window_max
    post_summary/summary_beta.csv
        case, geom, alpha_deg, step, side, beta_deg, theory_beta_deg,
        theory_state
    post_summary/CL_vs_alpha.png
    post_summary/CD_vs_alpha.png
    post_summary/Cm_vs_alpha.png
    post_summary/beta_vs_alpha.png

The time-averaging window is picked as the last `--tail-fraction` of
the available snapshots (default 0.4 = last 40%).

Usage
    # run from the airfoil_static_euler directory
    python3 post_aggregate.py --root .  --tail-fraction 0.4
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_case_tag(name):
    """case directory like 'plot_mid_a12' or 'surf_x018_a04' -> ('mid', 12) etc."""
    m = re.match(r"(?:plot|surf|post)_([a-zA-Z0-9]+)_a(\d+)", name)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def read_integrals(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append({
                "step": int(r["step"]),
                "CL": float(r["CL"]),
                "CD": float(r["CD"]),
                "Cm_LE": float(r["Cm_LE"]),
                "Cm_c4": float(r["Cm_c4"]),
                "xcp": float(r["x_cp_over_c"]) if r["x_cp_over_c"] else None,
            })
    rows.sort(key=lambda x: x["step"])
    return rows


def read_beta(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                b = float(r["beta_deg"]) if r["beta_deg"] else None
            except ValueError:
                b = None
            rows.append({
                "step": int(r["step"]),
                "side": r["side"],
                "beta_deg": b,
            })
    rows.sort(key=lambda x: (x["side"], x["step"]))
    return rows


def read_theory(root):
    """Return dict (geom, alpha, side) -> (beta or None, state)."""
    path = os.path.join(root, "theory_beta.csv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for r in csv.DictReader(f):
            beta = float(r["beta_weak_deg"]) if r["beta_weak_deg"] else None
            out[(r["geometry"], int(float(r["alpha_deg"])),
                 r["side"])] = (beta, r["state"])
    return out


def average_tail(values, fraction):
    n = len(values)
    if n == 0:
        return float("nan"), 0, -1, -1
    start = int((1.0 - fraction) * n)
    vs = values[start:]
    return (float(np.mean([v for v in vs if v is not None])),
            len(vs),
            values[start]["step"] if hasattr(values[start], "get") else start,
            -1)


def aggregate_integrals(rows, fraction):
    n = len(rows)
    if n == 0:
        return {}
    start = int((1.0 - fraction) * n)
    tail = rows[start:]
    keys = ["CL", "CD", "Cm_LE", "Cm_c4"]
    means = {}
    for k in keys:
        vals = [r[k] for r in tail]
        means[k] = float(np.mean(vals)) if vals else float("nan")
    xcp_vals = [r["xcp"] for r in tail if r["xcp"] is not None]
    means["xcp"] = (float(np.mean(xcp_vals)) if xcp_vals else None)
    return {
        "CL_mean":   means["CL"],
        "CD_mean":   means["CD"],
        "CmLE_mean": means["Cm_LE"],
        "Cmc4_mean": means["Cm_c4"],
        "xcp_mean":  means["xcp"],
        "n_steps_averaged": len(tail),
        "window_min": tail[0]["step"] if tail else -1,
        "window_max": tail[-1]["step"] if tail else -1,
    }


def aggregate_beta(rows, fraction):
    """Return dict side -> (beta_mean or None, n_averaged)."""
    out = {}
    for side in ("upper", "lower"):
        side_rows = [r for r in rows if r["side"] == side]
        n = len(side_rows)
        start = int((1.0 - fraction) * n)
        tail = side_rows[start:]
        betas = [r["beta_deg"] for r in tail if r["beta_deg"] is not None]
        if betas:
            out[side] = (float(np.mean(betas)), len(betas))
        else:
            out[side] = (None, 0)
    return out


def discover_cases(root):
    """Find case directories; return list of (case_dir, geom, alpha)."""
    out = []
    for d in sorted(os.listdir(root)):
        full = os.path.join(root, d)
        if not os.path.isdir(full):
            continue
        if not (d.startswith("plot_") or d.startswith("post_")):
            continue
        # case is identified by the 'plot_<geom>_a<NN>' form
        geom, alpha = parse_case_tag(d)
        if geom is None:
            continue
        out.append((os.path.dirname(full), geom, alpha, d))
    # dedupe by (geom, alpha)
    seen = set()
    unique = []
    for root_dir, geom, alpha, plot_dir_name in out:
        key = (geom, alpha)
        if key in seen:
            continue
        seen.add(key)
        unique.append((root_dir, geom, alpha))
    return sorted(unique, key=lambda x: (x[1], x[2]))


def plot_vs_alpha(all_rows, key, ylabel, out_path, title):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for geom, marker, color in (("mid", "o", "C0"), ("x018", "s", "C3")):
        rows = [r for r in all_rows if r["geom"] == geom]
        rows.sort(key=lambda r: r["alpha_deg"])
        xs = [r["alpha_deg"] for r in rows]
        ys = [r[key] for r in rows]
        ax.plot(xs, ys, marker=marker, color=color, linestyle="-",
                label=f"{geom}")
    ax.set_xlabel("angle of attack [deg]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print("Saved", out_path)


def plot_beta_vs_alpha(all_rows, theory, out_path):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    colors = {"mid": "C0", "x018": "C3"}
    for geom in ("mid", "x018"):
        rows = [r for r in all_rows if r["geom"] == geom]
        rows.sort(key=lambda r: r["alpha_deg"])
        # two lines: upper/lower
        for side, ls, mk in (("upper", "-", "o"), ("lower", "--", "v")):
            xs, ys, ytheory = [], [], []
            for r in rows:
                b = r["beta_" + side]
                if b is None:
                    continue
                xs.append(r["alpha_deg"])
                ys.append(b)
                key = (geom, int(r["alpha_deg"]), side)
                t = theory.get(key, (None, None))
                ytheory.append(t[0])
            if xs:
                ax.plot(xs, ys, linestyle=ls, marker=mk, color=colors[geom],
                        label=f"{geom} {side} (num)")
                th_x = [x for x, y in zip(xs, ytheory) if y is not None]
                th_y = [y for y in ytheory if y is not None]
                if th_y:
                    ax.plot(th_x, th_y, linestyle=":", color=colors[geom],
                            marker="x",
                            label=f"{geom} {side} (theory)")
    # Reference Mach wave at M=2: μ = asin(1/M) = 30°. On faces where
    # θ_eff < 0 (expansion, state=="expansion"), the extractor can only
    # see the leading Mach wave of the Prandtl-Meyer fan, which sits at
    # μ from the upstream flow direction. Drawing this as a guide line
    # makes the expansion regime visually distinguishable from attached
    # shock regimes.
    mu_M2 = math.degrees(math.asin(1.0 / 2.0))
    ax.axhline(mu_M2, color="gray", linestyle="-.", lw=1.0,
               label=f"Mach wave μ={mu_M2:.0f}° (M=2)")
    ax.set_xlabel("angle of attack [deg]")
    ax.set_ylabel("leading-edge shock angle β [deg]")
    ax.set_title("Shock angle: numerical vs weak-shock theory (M=2)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print("Saved", out_path)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", default=".",
                   help="airfoil_static_euler/ root containing case subdirs")
    p.add_argument("--tail-fraction", type=float, default=0.4,
                   help="fraction of tail snapshots to time-average (default 0.4)")
    p.add_argument("--out-dir", default=None)
    args = p.parse_args(argv)

    root = os.path.abspath(args.root)
    out_dir = args.out_dir or os.path.join(root, "post_summary")
    os.makedirs(out_dir, exist_ok=True)

    cases = discover_cases(root)
    if not cases:
        print("no case directories found under", root); return 1

    theory = read_theory(root)
    all_integrals = []
    all_beta = []
    for (root_dir, geom, alpha) in cases:
        # Conventional dir layout: plot_<geom>_a<NN>/ + post/<same tag>/
        # We look for: <root>/plot_<geom>_a<NN>/post/integrals.csv etc.
        case_name = f"{geom}_a{alpha:02d}"
        case_post = os.path.join(root, f"plot_{case_name}", "post")
        if not os.path.isdir(case_post):
            # fallback: post_<case_name>/
            case_post = os.path.join(root, "post", case_name)
        int_csv = os.path.join(case_post, "integrals.csv")
        beta_csv = os.path.join(case_post, "beta.csv")

        entry = {
            "case": case_name, "geom": geom, "alpha_deg": alpha,
            "CL_mean": None, "CD_mean": None,
            "CmLE_mean": None, "Cmc4_mean": None, "xcp_mean": None,
            "beta_upper": None, "beta_lower": None,
        }

        if os.path.exists(int_csv):
            rows = read_integrals(int_csv)
            agg = aggregate_integrals(rows, args.tail_fraction)
            entry.update(agg)
        if os.path.exists(beta_csv):
            rows = read_beta(beta_csv)
            b = aggregate_beta(rows, args.tail_fraction)
            entry["beta_upper"] = b["upper"][0]
            entry["beta_lower"] = b["lower"][0]

        all_integrals.append(entry)

    # summary CSV
    with open(os.path.join(out_dir, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        cols = ["case", "geom", "alpha_deg", "CL_mean", "CD_mean",
                "CmLE_mean", "Cmc4_mean", "xcp_mean",
                "beta_upper", "beta_lower"]
        w.writerow(cols)
        for e in all_integrals:
            w.writerow([e[c] if e[c] is not None else "" for c in cols])
    print("Saved", os.path.join(out_dir, "summary.csv"))

    # plots
    if any(e["CL_mean"] is not None for e in all_integrals):
        plot_vs_alpha(all_integrals, "CL_mean", "CL",
                      os.path.join(out_dir, "CL_vs_alpha.png"),
                      "Lift coefficient vs angle of attack")
        plot_vs_alpha(all_integrals, "CD_mean", "CD",
                      os.path.join(out_dir, "CD_vs_alpha.png"),
                      "Drag coefficient vs angle of attack")
        plot_vs_alpha(all_integrals, "Cmc4_mean", "Cm (c/4)",
                      os.path.join(out_dir, "Cm_vs_alpha.png"),
                      "Pitching moment about c/4 vs alpha")
    if any(e["beta_upper"] is not None or e["beta_lower"] is not None
           for e in all_integrals):
        plot_beta_vs_alpha(all_integrals, theory,
                           os.path.join(out_dir, "beta_vs_alpha.png"))

    return 0


if __name__ == "__main__":
    sys.exit(main())
