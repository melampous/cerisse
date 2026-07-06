#!/usr/bin/env python3
"""
audit_surf_vtp.py -- hard validity audit of IBM surface VTP exports.

For every VTP given, parses ALL point/cell DataArrays (ascii XML VTP as
written by cerisse ibm_solver_io) and reports per-field:
  - NaN / Inf counts
  - Pressure: count P < abs_floor (default 1000 Pa) and P < 0.01*p_inf
    (p_inf given per case via --pinf; --pinf rel uses 0.01*median(P) of
    each file as the relative floor, for nondimensional cases)
  - Temperature: count T <= 0 or NaN
  - Density (rho/Density) if present: rho <= 0 or NaN
  - min / max / median per field
  - bad-cell counts grouped by Level and by any owner id array present
    (Rank / fab / owner / grid)

Output: one summary row per VTP (CSV to --out and/or stdout table) plus a
family rollup (aggregated over all files given in one invocation).

Usage:
  audit_surf_vtp.py --label cyl_e2_g4096 --pinf rel  dir/*.vtp
  audit_surf_vtp.py --label airfoil_re1e4 --pinf 101325 --out audit.csv f1.vtp f2.vtp
  audit_surf_vtp.py --list-arrays f.vtp      # just show available arrays

Exit code 1 if any file has bad cells (NaN / floor violations), else 0.
"""
import argparse
import glob
import math
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

OWNER_ID_HINTS = ("rank", "fab", "owner", "grid", "proc", "box")
ABS_P_FLOOR_DEFAULT = 1000.0  # Pa, meaningful for dimensional cases only


def parse_vtp(path):
    """Return dict {array_name: np.ndarray} for all Point/Cell DataArrays,
    plus '_npoints'/'_ncells'. ascii format only."""
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find(".//Piece")
    out = {
        "_npoints": int(piece.get("NumberOfPoints", 0)),
        "_ncells": int(piece.get("NumberOfPolys", piece.get("NumberOfCells", 0))),
    }
    for section, tag in (("point", "PointData"), ("cell", "CellData")):
        sec = piece.find(tag)
        if sec is None:
            continue
        for da in sec.findall("DataArray"):
            name = da.get("Name")
            if da.get("format", "ascii") != "ascii":
                out[name] = None  # unsupported format; flagged upstream
                continue
            dtype = float if da.get("type", "").startswith("Float") else np.int64
            txt = da.text or ""
            try:
                arr = np.fromstring(txt, sep=" ") if dtype is float else \
                    np.array(txt.split(), dtype=np.int64)
            except Exception:
                arr = np.fromiter((float(t) for t in txt.split()), dtype=float)
            out[name] = arr
            out.setdefault("_sections", {})[name] = section
    return out


def field_stats(a):
    if a is None or a.size == 0:
        return dict(n=0, nan=0, inf=0, mn=math.nan, mx=math.nan, med=math.nan)
    fa = a.astype(float)
    nan = int(np.isnan(fa).sum())
    inf = int(np.isinf(fa).sum())
    fin = fa[np.isfinite(fa)]
    if fin.size == 0:
        return dict(n=a.size, nan=nan, inf=inf, mn=math.nan, mx=math.nan, med=math.nan)
    return dict(n=a.size, nan=nan, inf=inf,
                mn=float(fin.min()), mx=float(fin.max()), med=float(np.median(fin)))


def group_counts(mask, key_arr):
    """counts of True in mask grouped by integer key array."""
    if key_arr is None or key_arr.size != mask.size:
        return {}
    out = {}
    for k in np.unique(key_arr):
        out[int(k)] = int(mask[key_arr == k].sum())
    return out


def audit_file(path, pinf_mode, abs_floor):
    d = parse_vtp(path)
    arrays = {k: v for k, v in d.items() if not k.startswith("_")}
    ncells = d["_ncells"]

    P = arrays.get("Pressure")
    T = arrays.get("Temperature")
    rho = arrays.get("rho", arrays.get("Density"))
    level = arrays.get("Level")
    owner_keys = [k for k in arrays if any(h in k.lower() for h in OWNER_ID_HINTS)
                  and k not in ("Level",)]

    rec = {"file": path, "ncells": ncells,
           "arrays": ",".join(sorted(arrays.keys()))}

    # per-field stats
    for name, a in arrays.items():
        st = field_stats(a)
        rec[f"{name}.nan"] = st["nan"]
        rec[f"{name}.min"] = st["mn"]
        rec[f"{name}.max"] = st["mx"]
        rec[f"{name}.med"] = st["med"]

    # pressure floors
    bad_total = np.zeros(ncells, dtype=bool)
    if P is not None and P.size:
        pnan = ~np.isfinite(P)
        if pinf_mode == "rel":
            med = float(np.median(P[np.isfinite(P)])) if np.isfinite(P).any() else math.nan
            rel_floor = 0.01 * med
            rec["pinf_used"] = f"rel(0.01*med={rel_floor:.4g})"
        else:
            rel_floor = 0.01 * float(pinf_mode)
            rec["pinf_used"] = f"{float(pinf_mode):g}"
        p_lt_abs = np.isfinite(P) & (P < abs_floor)
        p_lt_rel = np.isfinite(P) & (P < rel_floor)
        rec["P.nan"] = int(pnan.sum())
        rec["P.lt_1000Pa"] = int(p_lt_abs.sum())
        rec["P.lt_0.01pinf"] = int(p_lt_rel.sum())
        rec["P.le_0"] = int((np.isfinite(P) & (P <= 0)).sum())
        bad_total |= pnan | p_lt_rel
    if T is not None and T.size:
        t_bad = ~np.isfinite(T) | (T <= 0)
        rec["T.le0_or_nan"] = int(t_bad.sum())
        bad_total |= t_bad
    if rho is not None and rho.size:
        r_bad = ~np.isfinite(rho) | (rho <= 0)
        rec["rho.le0_or_nan"] = int(r_bad.sum())
        bad_total |= r_bad
    else:
        rec["rho.le0_or_nan"] = "n/a"
    # tau fields: NaN only (sign is physical)
    for name in arrays:
        if name.lower().startswith("tau") or name == "dTdn":
            a = arrays[name]
            if a is not None and a.size:
                bad_total |= ~np.isfinite(a)

    rec["bad_cells"] = int(bad_total.sum())
    rec["bad_pct"] = 100.0 * rec["bad_cells"] / max(ncells, 1)
    rec["bad_by_level"] = group_counts(bad_total, level)
    rec["cells_by_level"] = group_counts(np.ones(ncells, dtype=bool), level)
    for ok in owner_keys:
        rec[f"bad_by_{ok}"] = group_counts(bad_total, arrays[ok])
    return rec


def merge_dicts(a, b):
    for k, v in b.items():
        a[k] = a.get(k, 0) + v
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="VTP files (or globs)")
    ap.add_argument("--pinf", default="rel",
                    help="freestream pressure in file units, or 'rel' to use "
                         "0.01*median(P) per file as the relative floor "
                         "(for nondimensional cases). Default: rel")
    ap.add_argument("--abs-floor", type=float, default=ABS_P_FLOOR_DEFAULT,
                    help="absolute P floor count threshold (default 1000; "
                         "only meaningful for dimensional Pa data)")
    ap.add_argument("--label", default="family", help="family label for rollup")
    ap.add_argument("--out", default=None, help="write per-file CSV here")
    ap.add_argument("--list-arrays", action="store_true",
                    help="only list available arrays of first file and exit")
    args = ap.parse_args()

    files = []
    for f in args.files:
        files.extend(sorted(glob.glob(f)) if any(c in f for c in "*?[") else [f])
    if not files:
        print("no files matched", file=sys.stderr)
        return 2

    if args.list_arrays:
        d = parse_vtp(files[0])
        for k, v in d.items():
            if not k.startswith("_"):
                sec = d.get("_sections", {}).get(k, "?")
                print(f"{k:20s} [{sec} data] n={0 if v is None else v.size}")
        return 0

    recs = []
    for f in files:
        try:
            recs.append(audit_file(f, args.pinf, args.abs_floor))
        except Exception as e:
            recs.append({"file": f, "error": str(e)})

    # ---- per-file table ----
    cols = ["file", "ncells", "P.nan", "P.le_0", "P.lt_1000Pa", "P.lt_0.01pinf",
            "T.le0_or_nan", "rho.le0_or_nan", "bad_cells", "bad_pct",
            "Pressure.min", "Pressure.med", "Pressure.max",
            "Temperature.min", "Temperature.med", "Temperature.max",
            "Tau1.nan", "Tau1.min", "Tau1.med", "Tau1.max",
            "Tau2.min", "Tau2.med", "Tau2.max",
            "dTdn.min", "dTdn.med", "dTdn.max",
            "bad_by_level", "bad_by_Rank", "pinf_used"]
    print(f"\n=== per-file audit: {args.label} (pinf={args.pinf}) ===")
    hdr = ["file", "ncells", "Pnan", "P<=0", "P<1kPa", "P<.01pinf", "Tbad",
           "bad", "bad%", "Pmin", "Pmed", "Pmax", "Tmin", "Tmax", "byLevel"]
    print(" ".join(f"{h:>10s}" for h in hdr))
    for r in recs:
        if "error" in r:
            print(f"{os.path.basename(r['file'])}: ERROR {r['error']}")
            continue
        print(" ".join([
            f"{os.path.basename(r['file'])[-28:]:>28s}",
            f"{r['ncells']:>7d}",
            f"{r.get('P.nan', 0):>5}",
            f"{r.get('P.le_0', 0):>5}",
            f"{r.get('P.lt_1000Pa', 0):>7}",
            f"{r.get('P.lt_0.01pinf', 0):>9}",
            f"{r.get('T.le0_or_nan', 0):>5}",
            f"{r.get('bad_cells', 0):>5}",
            f"{r.get('bad_pct', 0):>6.2f}",
            f"{r.get('Pressure.min', math.nan):>11.4g}",
            f"{r.get('Pressure.med', math.nan):>11.4g}",
            f"{r.get('Pressure.max', math.nan):>11.4g}",
            f"{r.get('Temperature.min', math.nan):>10.4g}",
            f"{r.get('Temperature.max', math.nan):>10.4g}",
            f"{r.get('bad_by_level', {})}",
        ]))

    # ---- family rollup ----
    ok = [r for r in recs if "error" not in r]
    tot_cells = sum(r["ncells"] for r in ok)
    tot_bad = sum(r["bad_cells"] for r in ok)
    lvl_bad, lvl_all = {}, {}
    rank_bad = {}
    for r in ok:
        merge_dicts(lvl_bad, r.get("bad_by_level", {}))
        merge_dicts(lvl_all, r.get("cells_by_level", {}))
        merge_dicts(rank_bad, r.get("bad_by_Rank", {}))
    glob_min = {}
    for fld in ("Pressure", "Temperature", "Tau1", "Tau2", "dTdn"):
        mns = [r.get(f"{fld}.min") for r in ok if f"{fld}.min" in r]
        mxs = [r.get(f"{fld}.max") for r in ok if f"{fld}.max" in r]
        meds = [r.get(f"{fld}.med") for r in ok if f"{fld}.med" in r]
        if mns:
            glob_min[fld] = (np.nanmin(mns), np.nanmedian(meds), np.nanmax(mxs))
    print(f"\n=== rollup: {args.label} ===")
    print(f"files={len(ok)} (+{len(recs)-len(ok)} errors)  cells={tot_cells}  "
          f"bad={tot_bad} ({100.0*tot_bad/max(tot_cells,1):.3f}%)")
    for lv in sorted(lvl_all):
        b = lvl_bad.get(lv, 0)
        print(f"  Level {lv}: {b}/{lvl_all[lv]} bad "
              f"({100.0*b/max(lvl_all[lv],1):.3f}%)")
    if any(v for v in rank_bad.values()):
        print(f"  bad by Rank: { {k: v for k, v in rank_bad.items() if v} }")
    for fld, (mn, md, mx) in glob_min.items():
        print(f"  {fld:12s} min={mn:.5g} med~{md:.5g} max={mx:.5g}")
    verdict = "CLEAN" if tot_bad == 0 else "CONTAMINATED"
    print(f"  VERDICT: {verdict}")

    if args.out:
        import csv
        allkeys = []
        for r in recs:
            for k in r:
                if k not in allkeys:
                    allkeys.append(k)
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=allkeys)
            w.writeheader()
            for r in recs:
                w.writerow({k: r.get(k, "") for k in allkeys})
        print(f"CSV written: {args.out}")

    return 1 if tot_bad else 0


if __name__ == "__main__":
    sys.exit(main())
