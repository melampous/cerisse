#!/usr/bin/env python3
"""Extract surface Cp on the 4 panels from a Cerisse IBM surf VTP.

Panel classification is done in body-frame coordinates, NEVER from
the sign of world-frame y. The body-frame is obtained by inverse-
rotating the world-frame surface point about the origin by -alpha:

    x_body =  cos(alpha) * x_world - sin(alpha) * y_world
    y_body =  sin(alpha) * x_world + cos(alpha) * y_world

Then with the body-frame 4-vertex diamond:
    LE    = (-0.05, 0)
    Upper = (x_maxT,  +0.00437)     x_maxT = 0 for mid, -0.032 for x018
    TE    = (+0.05, 0)
    Lower = (x_maxT,  -0.00437)

the panel of each sample is decided by

    y_body > 0 and x_body < x_maxT   -> fore_upper
    y_body > 0 and x_body > x_maxT   -> aft_upper
    y_body < 0 and x_body > x_maxT   -> aft_lower
    y_body < 0 and x_body < x_maxT   -> fore_lower

Output: one CSV per surface snapshot (step) with columns
    s_panel, x, y, x_body, y_body, panel_name, Cp
plus a summary CSV that takes the time-averaged Cp over a user-
picked step window.

Usage
    python3 post_extract_cp.py <case_dir> --alpha 12 --geom mid
    python3 post_extract_cp.py <case_dir> --alpha 12 --geom mid \
            --step-window 2500:5000 --out-dir post/

Currently reads the surf VTP files written via ib.plot_surf=1.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
import xml.etree.ElementTree as ET


P_INF  = 101325.0
T_INF  = 300.0
RHO_INF = 1.1768
U_INF  = 694.43
Q_INF  = 0.5 * RHO_INF * U_INF * U_INF   # dynamic pressure, ~2.838e5 Pa

GEOM_X_MAX_THICKNESS = {"mid": 0.0, "x018": -0.032}


def inverse_rotate(x, y, alpha_deg):
    """Map world-frame (x, y) back to body-frame by rotating by -alpha.

    (alpha is the clockwise rotation applied to the geometry.)
    """
    a = math.radians(alpha_deg)
    ca, sa = math.cos(a), math.sin(a)
    return ca * x - sa * y, sa * x + ca * y


def classify_panel(x_body, y_body, x_maxT):
    """Return one of fore_upper / aft_upper / aft_lower / fore_lower."""
    if y_body >= 0.0:
        return "fore_upper" if x_body < x_maxT else "aft_upper"
    else:
        return "fore_lower" if x_body < x_maxT else "aft_lower"


def read_vtp(path):
    """Minimal VTP polydata reader. Returns (points, fields) where
    points is a list of (x, y) tuples (z dropped) and fields is a
    dict {name: list[float]} with per-point scalar data arrays.
    """
    ns = ""  # no namespace
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find(".//Piece")

    # points
    pts_data = piece.find("Points").find("DataArray").text.split()
    ncomp = int(piece.find("Points").find("DataArray").attrib.get(
                "NumberOfComponents", "3"))
    pts = []
    for i in range(0, len(pts_data), ncomp):
        x = float(pts_data[i])
        y = float(pts_data[i + 1])
        pts.append((x, y))

    # point data + cell data (Cerisse writes most scalars to CellData)
    fields = {}
    for tag in ("PointData", "CellData"):
        d = piece.find(tag)
        if d is None:
            continue
        for da in d.findall("DataArray"):
            name = da.attrib["Name"]
            vals = [float(v) for v in da.text.split()]
            fields[name] = vals
    return pts, fields


def step_from_filename(path):
    m = re.search(r"(\d+)\.vtp$", os.path.basename(path))
    return int(m.group(1)) if m else -1


def process_one(vtp_path, alpha_deg, geom, out_csv):
    if geom not in GEOM_X_MAX_THICKNESS:
        raise SystemExit(f"unknown geometry {geom!r}")
    x_maxT = GEOM_X_MAX_THICKNESS[geom]

    pts, fields = read_vtp(vtp_path)
    if "Pressure" not in fields:
        raise SystemExit("Pressure field not found in " + vtp_path)
    P = fields["Pressure"]

    step = step_from_filename(vtp_path)

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "i", "x", "y", "x_body", "y_body",
                    "panel", "p", "Cp"])
        for i, ((x, y), p) in enumerate(zip(pts, P)):
            xb, yb = inverse_rotate(x, y, alpha_deg)
            panel = classify_panel(xb, yb, x_maxT)
            cp = (p - P_INF) / Q_INF
            w.writerow([step, i, x, y, xb, yb, panel, p, cp])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case_dir")
    p.add_argument("--alpha", type=float, required=True,
                   help="AoA in degrees (the clockwise rotation applied "
                        "to the geometry)")
    p.add_argument("--geom", choices=list(GEOM_X_MAX_THICKNESS.keys()),
                   required=True)
    p.add_argument("--surf-glob",
                   default="surf_*/surf/*.vtp",
                   help="glob for VTP files relative to case_dir")
    p.add_argument("--out-dir", default=None,
                   help="where to write per-snapshot CSV "
                        "(default: <case_dir>/post)")
    args = p.parse_args(argv)

    case_dir = os.path.abspath(args.case_dir)
    out_dir = args.out_dir or os.path.join(case_dir, "post")
    os.makedirs(out_dir, exist_ok=True)

    vtps = sorted(glob.glob(os.path.join(case_dir, args.surf_glob)))
    if not vtps:
        # also try looser path (surf has subdir per body)
        vtps = sorted(glob.glob(os.path.join(case_dir, "surf*/*/*.vtp")))
    if not vtps:
        print("no VTP files found under", case_dir)
        return 1

    print(f"# geometry={args.geom}  alpha={args.alpha} deg  "
          f"x_maxT={GEOM_X_MAX_THICKNESS[args.geom]:+.3f}")
    print(f"# found {len(vtps)} snapshots")
    for vtp in vtps:
        step = step_from_filename(vtp)
        out = os.path.join(out_dir, f"cp_step{step:05d}.csv")
        process_one(vtp, args.alpha, args.geom, out)
        print(f"  step {step:05d}  ->  {os.path.relpath(out, case_dir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
