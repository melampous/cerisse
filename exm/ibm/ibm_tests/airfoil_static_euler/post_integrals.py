#!/usr/bin/env python3
"""Integrate surface pressure on a Cerisse IBM surf VTP to CL, CD, Cm,
and x_cp.

Conventions (must match README_case.md):

* Freestream along +x; freestream dynamic pressure q_inf = 0.5 rho u^2.
* Reference chord c = 0.1 m.
* Positive CL = upward force (+y world).
* Cm_LE reference point: body-frame (-0.05, 0), world-frame = rotated.
* Cm_c4 reference point: body-frame (-0.025, 0), world-frame = rotated.
* x_cp / c = -Cm_LE / CL, reported only when |CL| > CL_EPS. For
  |CL| <= CL_EPS we emit x_cp = '' to CSV and only report Cm.

The Euler / slip-wall case has zero shear, so the surface force is
purely the pressure integral:

    F = - sum_i (p_i - p_inf) * n_i * ds_i

where n_i is the outward unit normal and ds_i is the element arc
length. In the VTP the polys are already line-segment faces; each
"face" between two consecutive surface points gives a segment. The
exported fields are per-POINT, so the per-segment pressure is the
average of the two endpoint pressures.

Usage
    python3 post_integrals.py <case_dir> --alpha 12 --geom mid
    python3 post_integrals.py <case_dir> --alpha 12 --geom mid \
            --out post/integrals_mid_a12.csv
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


P_INF   = 101325.0
RHO_INF = 1.1768
U_INF   = 694.43
Q_INF   = 0.5 * RHO_INF * U_INF * U_INF
C_REF   = 0.1
CL_EPS  = 1.0e-4  # fraction of q_inf*c below which x_cp is suppressed


def rotate_cw(x, y, alpha_deg):
    a = math.radians(alpha_deg)
    c, s = math.cos(a), math.sin(a)
    return c * x + s * y, -s * x + c * y


def read_vtp(path):
    """Return (points[(x,y)], fields{name: list}, conn (list of tuples))."""
    tree = ET.parse(path)
    piece = tree.getroot().find(".//Piece")

    da_pts = piece.find("Points").find("DataArray")
    ncomp = int(da_pts.attrib.get("NumberOfComponents", "3"))
    raw = da_pts.text.split()
    pts = [(float(raw[i]), float(raw[i + 1]))
           for i in range(0, len(raw), ncomp)]

    fields = {}
    for tag in ("PointData", "CellData"):
        d = piece.find(tag)
        if d is None:
            continue
        for da in d.findall("DataArray"):
            fields[da.attrib["Name"]] = [float(v) for v in da.text.split()]

    conn = []
    polys = piece.find("Polys")
    if polys is not None:
        c_arr = None
        o_arr = None
        for da in polys.findall("DataArray"):
            if da.attrib["Name"] == "connectivity":
                c_arr = [int(v) for v in da.text.split()]
            elif da.attrib["Name"] == "offsets":
                o_arr = [int(v) for v in da.text.split()]
        if c_arr is not None and o_arr is not None:
            prev = 0
            for off in o_arr:
                conn.append(tuple(c_arr[prev:off]))
                prev = off
    return pts, fields, conn


def step_from_filename(path):
    m = re.search(r"(\d+)\.vtp$", os.path.basename(path))
    return int(m.group(1)) if m else -1


def integrate_one(vtp_path, alpha_deg, geom, case_name):
    pts, fields, conn = read_vtp(vtp_path)
    if "Pressure" not in fields:
        raise SystemExit("no Pressure field in " + vtp_path)
    P = fields["Pressure"]

    # Reference points in world frame.
    le_world = rotate_cw(-0.05, 0.0, alpha_deg)
    qc_world = rotate_cw(-0.025, 0.0, alpha_deg)

    CL_world = 0.0  # +y
    CD_world = 0.0  # +x
    Mz_LE = 0.0     # moment about LE (z-axis into page is positive)
    Mz_QC = 0.0     # moment about c/4

    # Treat each polygon as a directed 2-vertex segment (edge).
    # Segment orientation follows .dat vertex order -> outward normal
    # is to the right of the segment (hand-rule in 2D with z out of
    # page means the body lies to the LEFT of the segment).
    #
    # Cerisse's exporter writes per-body surface polygons consistent
    # with that convention, so the force on the body from pressure is
    #     F = sum (p - p_inf) * n_hat * ds
    # where n_hat points OUTWARD (away from the body).
    # Cerisse exports each wall panel as a 2-vertex Polys entry;
    # Pressure is per-cell. For >2-vertex polys, treat as closed and
    # sum edges (but average P over the two point indices is not
    # meaningful for cell data, so use the cell's single P).
    for cell_idx, poly in enumerate(conn):
        if len(poly) < 2:
            continue
        if len(poly) == 2:
            edges = [(poly[0], poly[1])]
        else:
            edges = [(poly[k], poly[(k + 1) % len(poly)])
                     for k in range(len(poly))]
        p_seg = P[cell_idx]
        dp = p_seg - P_INF
        for i0, i1 in edges:
            if i0 == i1:
                continue
            x0, y0 = pts[i0]
            x1, y1 = pts[i1]
            dx = x1 - x0
            dy = y1 - y0
            ds = math.hypot(dx, dy)
            if ds == 0:
                continue
            # outward normal (body to the left, normal to the right)
            nx = dy / ds
            ny = -dx / ds
            # force contribution per segment
            fx = dp * nx * ds
            fy = dp * ny * ds
            CD_world += fx
            CL_world += fy
            # moment arms
            xm = 0.5 * (x0 + x1)
            ym = 0.5 * (y0 + y1)
            Mz_LE += (xm - le_world[0]) * fy - (ym - le_world[1]) * fx
            Mz_QC += (xm - qc_world[0]) * fy - (ym - qc_world[1]) * fx

    # Non-dimensionalize (per unit span, 2D)
    # Force on body = -sum (p - p_inf) n_hat ds  (n outward), so
    # flip sign of the accumulator above.
    norm = Q_INF * C_REF
    CL = -CL_world / norm
    CD = -CD_world / norm
    Cm_LE = -Mz_LE / (norm * C_REF)
    Cm_c4 = -Mz_QC / (norm * C_REF)

    # x_cp (distance from LE along chord, nondimensionalized by chord).
    # With our sign convention (force ON body, CL>0 for nose-up α),
    # positive Cm_LE about LE is nose-up, which for a lift force at
    # x_cp downstream of the LE gives Mz_LE = +x_cp * CL * c, so
    # x_cp/c = Cm_LE / CL.
    if abs(CL) > CL_EPS:
        x_cp_over_c = Cm_LE / CL
    else:
        x_cp_over_c = None

    return {
        "case": case_name,
        "geom": geom,
        "alpha_deg": alpha_deg,
        "step": step_from_filename(vtp_path),
        "CL": CL,
        "CD": CD,
        "Cm_LE": Cm_LE,
        "Cm_c4": Cm_c4,
        "x_cp_over_c": x_cp_over_c,
        "q_inf": Q_INF,
        "c_ref": C_REF,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case_dir")
    p.add_argument("--alpha", type=float, required=True)
    p.add_argument("--geom", choices=["mid", "x018"], required=True)
    p.add_argument("--surf-glob", default="surf_*/surf/*.vtp")
    p.add_argument("--out", default=None,
                   help="output CSV path (default: <case>/post/integrals.csv)")
    p.add_argument("--window-min", type=int, default=None)
    p.add_argument("--window-max", type=int, default=None)
    args = p.parse_args(argv)

    case_dir = os.path.abspath(args.case_dir)
    out_csv = args.out or os.path.join(case_dir, "post", "integrals.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    vtps = sorted(glob.glob(os.path.join(case_dir, args.surf_glob)))
    if not vtps:
        vtps = sorted(glob.glob(os.path.join(case_dir, "surf*/*/*.vtp")))
    if not vtps:
        print("no VTP files under", case_dir)
        return 1

    case_name = f"{args.geom}_a{int(args.alpha):02d}"
    rows = []
    for vtp in vtps:
        step = step_from_filename(vtp)
        if args.window_min is not None and step < args.window_min:
            continue
        if args.window_max is not None and step > args.window_max:
            continue
        r = integrate_one(vtp, args.alpha, args.geom, case_name)
        rows.append(r)

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "geom", "alpha_deg", "step",
                    "CL", "CD", "Cm_LE", "Cm_c4", "x_cp_over_c"])
        for r in rows:
            w.writerow([r["case"], r["geom"], r["alpha_deg"], r["step"],
                        f"{r['CL']:.6e}", f"{r['CD']:.6e}",
                        f"{r['Cm_LE']:.6e}", f"{r['Cm_c4']:.6e}",
                        "" if r["x_cp_over_c"] is None
                           else f"{r['x_cp_over_c']:.6e}"])

    print(f"# wrote {len(rows)} rows to {out_csv}")
    # quick sanity print
    if rows:
        last = rows[-1]
        print(f"# final step {last['step']}: CL={last['CL']:+.4e} "
              f"CD={last['CD']:+.4e} Cm_LE={last['Cm_LE']:+.4e} "
              f"Cm_c/4={last['Cm_c4']:+.4e} "
              + ("x_cp/c=n/a (|CL|<eps)"
                 if last['x_cp_over_c'] is None
                 else f"x_cp/c={last['x_cp_over_c']:+.4f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
