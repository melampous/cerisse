#!/usr/bin/env python3
"""Generate rotated diamond-wedge .dat files for the static Euler benchmark.

Convention (matches the thesis-level spec):
  * Freestream is fixed along +x.
  * Positive AoA = geometry rotates CLOCKWISE about origin (0,0);
    leading edge at (-0.05, 0) moves UP.
  * Rotation formula (clockwise by angle alpha):
        x' =  cos(a) * x + sin(a) * y
        y' = -sin(a) * x + cos(a) * y
  * All rotations are about the origin; quarter-chord is NOT used.

Two base geometries are supported:
  * 'mid'  : max thickness at midchord (x=0); front and rear half-angle 5.0 deg
  * 'x018' : max thickness at 0.18c from LE (x=-0.032);
             front half-angle 13.64 deg, rear half-angle 3.05 deg

The output files contain one vertex per line ('x y\n'), with the first
vertex repeated at the end to close the polygon (matching the existing
airfoil_static/diamond_wedge.dat format).

Usage:
    python3 gen_rotated_airfoils.py
    python3 gen_rotated_airfoils.py --alphas 0,4,8,12,16,18 --outdir .
"""

from __future__ import annotations

import argparse
import math
import os
import sys


# ---- geometry definitions --------------------------------------------------

# Each base geometry is a list of vertices in consistent order
# (counter-clockwise when y is up). Closure point (first vertex repeated)
# is added when writing.

BASE_GEOMS = {
    "mid": [
        (-0.05, 0.0),
        (0.0, 0.00437),
        (0.05, 0.0),
        (0.0, -0.00437),
    ],
    "x018": [
        (-0.05, 0.0),
        (-0.032, 0.00437),
        (0.05, 0.0),
        (-0.032, -0.00437),
    ],
}


def rotate_cw(verts, alpha_deg, xc=0.0, yc=0.0):
    """Rotate vertices clockwise by alpha_deg about (xc, yc)."""
    a = math.radians(alpha_deg)
    ca, sa = math.cos(a), math.sin(a)
    out = []
    for x, y in verts:
        xs, ys = x - xc, y - yc
        xr = ca * xs + sa * ys
        yr = -sa * xs + ca * ys
        out.append((xc + xr, yc + yr))
    return out


def signed_area(verts):
    """Shoelace; positive for CCW order (y up)."""
    n = len(verts)
    s = 0.0
    for i in range(n):
        x0, y0 = verts[i]
        x1, y1 = verts[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return 0.5 * s


def bbox(verts):
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    return (min(xs), max(xs), min(ys), max(ys))


def write_dat(path, verts):
    """Write one vertex per line, repeat the first at the end to close."""
    with open(path, "w") as f:
        for x, y in verts:
            f.write(f"{x:.10f} {y:.10f}\n")
        f.write(f"{verts[0][0]:.10f} {verts[0][1]:.10f}\n")


def sanity_check(name, base_verts, rotated_verts, alpha_deg):
    """Return (ok, reason_str, metadata_dict)."""
    area_base = signed_area(base_verts)
    area_rot = signed_area(rotated_verts)

    # sign must match (no normal flip)
    if (area_base > 0) != (area_rot > 0):
        return False, "signed-area sign changed", {}

    # magnitude must be preserved under rigid rotation
    if abs(abs(area_rot) - abs(area_base)) > 1e-9:
        return False, f"|area| drift {abs(area_rot) - abs(area_base):.3e}", {}

    # chord length invariant (distance LE -> TE)
    le_b, te_b = base_verts[0], base_verts[2]
    le_r, te_r = rotated_verts[0], rotated_verts[2]
    chord_b = math.hypot(te_b[0] - le_b[0], te_b[1] - le_b[1])
    chord_r = math.hypot(te_r[0] - le_r[0], te_r[1] - le_r[1])
    if abs(chord_r - chord_b) > 1e-9:
        return False, f"chord drift {chord_r - chord_b:.3e}", {}

    meta = {
        "geometry_name": name,
        "alpha_deg": alpha_deg,
        "thickness_location": "midchord" if name == "mid" else "0.18c",
        "bounding_box": bbox(rotated_verts),
        "signed_area": area_rot,
        "chord": chord_r,
        "vertex_order_ok": True,
    }
    return True, "OK", meta


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--alphas", default="0,4,8,12,16,18",
                   help="comma-separated AoA list in degrees")
    p.add_argument("--outdir", default=".",
                   help="output directory for .dat files")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)

    alphas = [float(a) for a in args.alphas.split(",") if a.strip()]
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    print(f"# gen_rotated_airfoils.py")
    print(f"# freestream: fixed along +x")
    print(f"# rotation  : clockwise about (0, 0), positive alpha -> LE up")
    print(f"# outdir    : {outdir}")
    print()

    all_ok = True
    for name, verts in BASE_GEOMS.items():
        # first write the un-rotated base
        base_path = os.path.join(outdir, f"diamond_wedge_{name}.dat")
        write_dat(base_path, verts)
        ok, why, meta = sanity_check(name, verts, verts, 0.0)
        print(f"{os.path.basename(base_path):40s}  {'ok' if ok else why}"
              f"  area={meta.get('signed_area', float('nan')):+.6f}"
              f"  chord={meta.get('chord', float('nan')):.4f}")

        for alpha in alphas:
            rot = rotate_cw(verts, alpha, 0.0, 0.0)
            out = os.path.join(outdir,
                               f"diamond_wedge_{name}_a{int(alpha):02d}.dat")
            write_dat(out, rot)
            ok, why, meta = sanity_check(name, verts, rot, alpha)
            if not ok:
                all_ok = False
            line = (f"{os.path.basename(out):40s}  "
                    f"{'ok' if ok else why}  "
                    f"area={meta.get('signed_area', float('nan')):+.6f}  "
                    f"chord={meta.get('chord', float('nan')):.4f}")
            if args.verbose:
                line += f"  bbox={meta.get('bounding_box')}"
            print(line)

    print()
    if not all_ok:
        print("FAILURE: at least one geometry failed sanity checks")
        sys.exit(1)
    print("All geometries pass: chord length, signed area, orientation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
