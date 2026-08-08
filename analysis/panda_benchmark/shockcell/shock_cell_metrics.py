#!/usr/bin/env python3
"""Comprehensive shock-cell metrics from centreline TRUE statistics.
2D (L0, L0hi, L1-L6) and 3D (L1-L5), m142.
Per configuration:
  x1 (first compression-peak = first cell closure), pre/post extrema of
  <M>, <rho>/rho_j, <p>/p_inf around x1, all peak positions x_n and
  spacings S_n, and rms levels (p_rms/p_inf, rho_rms/rho_j) at x1 and
  median over x/D in [1,6].
Peaks: box-smoothed <p> (width ~0.05D), scipy find_peaks (prominence 0.02
p_inf), parabolic sub-cell refinement. M from mean fields.
Outputs: printed table + shock_cell_metrics.csv + .npz"""
import numpy as np
from scipy.signal import find_peaks
from pathlib import Path
import csv

HERE = Path(__file__).parent
D = 0.0254
P_AMB = 99780.0
RHO_J = 1.6413
GAM, RGAS = 1.4, 287.0

CASES_2D = [
    ("2D-L0",   "AXT_L0lo.npz",         1584),
    ("2D-L0hi", "AXT_m142_2d_L0hi.npz",  792),
    ("2D-L1",   "AXT_L1lo.npz",          794),
    ("2D-L2",   "AXT_L2lo.npz",          397),
    ("2D-L3",   "AXT_L3.npz",            198),
    ("2D-L4",   "AXT_L4.npz",             99),
    ("2D-L5",   "AXT_L5.npz",             50),
    ("2D-L6",   "AXT_L6.npz",             25),
]
CASES_3D = [
    ("3D-L1", "AX3T_L1.npz", 1588),
    ("3D-L2", "AX3T_L2.npz",  794),
    ("3D-L3", "AX3T_L3.npz",  397),
    ("3D-L4", "AX3T_L4.npz",  198),
    ("3D-L5", "AX3T_L5.npz",   99),
]

def profiles_2d(fn):
    d = np.load(HERE / fn, allow_pickle=True)
    x = d["z"] / D
    p = d["stat_pressureMEAN"]
    ps = d["stat_pressureSQR"]
    rho = d["stat_DensityMEAN"]
    rs = d["stat_DensitySQR"]
    T = d["stat_temperatureMEAN"]
    u = np.abs(d["stat_y_velocityMEAN"])
    M = u / np.sqrt(GAM * RGAS * np.maximum(T, 1.0))
    prms = np.sqrt(np.maximum(ps - p**2, 0))
    rrms = np.sqrt(np.maximum(rs - rho**2, 0))
    return x, p / P_AMB, rho / RHO_J, M, prms / P_AMB, rrms / RHO_J

def profiles_3d(fn):
    d = np.load(HERE / fn, allow_pickle=True)
    x = d["x"] / D
    p = d["ax_pressureMEAN"]
    ps = d["ax_pressureSQR"]
    rho = d["ax_DensityMEAN"]
    rs = d["ax_DensitySQR"]
    T = d["ax_temperatureMEAN"]
    u = np.abs(d["ax_x_velocityMEAN"])
    M = u / np.sqrt(GAM * RGAS * np.maximum(T, 1.0))
    prms = np.sqrt(np.maximum(ps - p**2, 0))
    rrms = np.sqrt(np.maximum(rs - rho**2, 0))
    return x, p / P_AMB, rho / RHO_J, M, prms / P_AMB, rrms / RHO_J

def boxsmooth(a, w):
    if w < 2:
        return a.copy()
    k = np.ones(w) / w
    return np.convolve(a, k, mode="same")

def interp(x, f, xq):
    return float(np.interp(xq, x, f))

def analyze(tag, x, pn, rn, M, prms, rrms, dx_um):
    dx = dx_um * 1e-6 / D
    m = (x >= 0.45) & (x <= 10.0)
    xs, ps = x[m], pn[m]
    w = max(3, int(round(0.05 / max(dx, 1e-9))))
    w = min(w, 51)
    psm = boxsmooth(ps, w)
    dxg = xs[1] - xs[0]
    dist = max(1, int(round(0.5 / dxg)))
    psm2 = psm.copy()
    psm2[:w] = -np.inf; psm2[-w:] = -np.inf
    pk, props = find_peaks(psm2, prominence=0.03, distance=dist, height=1.0)
    xpk = []
    for j in pk:
        if 1 <= j < len(xs) - 1:
            a, b, c = psm[j-1], psm[j], psm[j+1]
            den = a - 2*b + c
            off = 0.5*(a - c)/den if abs(den) > 1e-12 else 0.0
            xpk.append(xs[j] + off*(xs[1]-xs[0]))
    xpk = np.array(xpk[:6])
    Sn = np.diff(xpk)
    res = {"tag": tag, "dx_um": dx_um, "xpk": xpk, "Sn": Sn}
    if len(xpk) >= 1:
        x1 = xpk[0]
        res["x1"] = x1
        # shock-closure position x_s: steepest mean-pressure rise before x1
        ws = (xs >= 0.5) & (xs <= x1)
        if ws.sum() > 3:
            grad = np.gradient(psm, xs)
            js = np.argmax(grad[ws])
            res["x_s"] = float(xs[ws][js])
        else:
            res["x_s"] = np.nan
        xc = res["x_s"] if np.isfinite(res.get("x_s", np.nan)) else x1
        wpre = (x >= xc - 0.5) & (x <= xc - 0.02)
        wpost = (x >= xc + 0.02) & (x <= xc + 0.6)
        if wpre.any():
            jp = np.argmax(M[wpre])
            res["x_pre"] = float(x[wpre][jp])
            res["M_pre"] = float(M[wpre][jp])
            res["p_pre"] = interp(x, pn, res["x_pre"])
            res["rho_pre"] = interp(x, rn, res["x_pre"])
        if wpost.any():
            jq = np.argmin(M[wpost])
            res["x_post"] = float(x[wpost][jq])
            res["M_post"] = float(M[wpost][jq])
            res["p_post"] = interp(x, pn, res["x_post"])
            res["rho_post"] = interp(x, rn, res["x_post"])
        res["M_x1"] = interp(x, M, x1)
        res["p_x1"] = interp(x, pn, x1)
        res["rho_x1"] = interp(x, rn, x1)
        res["prms_x1"] = interp(x, prms, x1)
        res["rrms_x1"] = interp(x, rrms, x1)
    m16 = (x >= 1) & (x <= 6)
    res["prms_med16"] = float(np.nanmedian(prms[m16]))
    res["rrms_med16"] = float(np.nanmedian(rrms[m16]))
    res["u_x"] = max(0.01, dx / 2)
    return res

ALL = []
for tag, fn, dx in CASES_2D:
    ALL.append(analyze(tag, *profiles_2d(fn), dx))
for tag, fn, dx in CASES_3D:
    ALL.append(analyze(tag, *profiles_3d(fn), dx))

hdr1 = (f"{'case':8s} {'dx':>5s} {'x_s':>6s} {'x1':>6s} {'u_x':>5s} "
        f"{'M_pre':>6s} {'M_x1':>5s} {'M_post':>6s} "
        f"{'p_pre':>6s} {'p_x1':>5s} {'p_post':>6s} "
        f"{'r_pre':>6s} {'r_x1':>5s} {'r_post':>6s} "
        f"{'prms@x1':>7s} {'rrms@x1':>7s} {'prms[1,6]':>9s}")
print(hdr1)
for r in ALL:
    print(f"{r['tag']:8s} {r['dx_um']:5.0f} {r.get('x_s', np.nan):6.3f} {r.get('x1', np.nan):6.3f} {r['u_x']:5.3f} "
          f"{r.get('M_pre', np.nan):6.3f} {r.get('M_x1', np.nan):5.3f} {r.get('M_post', np.nan):6.3f} "
          f"{r.get('p_pre', np.nan):6.3f} {r.get('p_x1', np.nan):5.3f} {r.get('p_post', np.nan):6.3f} "
          f"{r.get('rho_pre', np.nan):6.3f} {r.get('rho_x1', np.nan):5.3f} {r.get('rho_post', np.nan):6.3f} "
          f"{r.get('prms_x1', np.nan):7.3f} {r.get('rrms_x1', np.nan):7.3f} {r['prms_med16']:9.3f}")

print("\npeak positions x_n (x/D) and spacings S_n:")
for r in ALL:
    xs = " ".join(f"{v:6.3f}" for v in r["xpk"])
    ss = " ".join(f"{v:6.3f}" for v in r["Sn"])
    print(f"{r['tag']:8s} x_n: {xs}")
    print(f"{'':8s} S_n: {ss}")

with open(HERE / "shock_cell_metrics.csv", "w", newline="") as f:
    wcsv = csv.writer(f)
    wcsv.writerow(["case", "dx_um", "x_s", "x1", "u_x", "M_pre", "M_x1", "M_post",
                   "p_pre", "p_x1", "p_post", "rho_pre", "rho_x1", "rho_post",
                   "prms_x1", "rrms_x1", "prms_med16", "rrms_med16",
                   "x2", "x3", "x4", "x5", "x6", "S1", "S2", "S3", "S4", "S5"])
    for r in ALL:
        xp = list(r["xpk"]) + [np.nan]*6
        sn = list(r["Sn"]) + [np.nan]*5
        wcsv.writerow([r["tag"], r["dx_um"], r.get("x_s", np.nan), r.get("x1", np.nan), r["u_x"],
                       r.get("M_pre", np.nan), r.get("M_x1", np.nan), r.get("M_post", np.nan),
                       r.get("p_pre", np.nan), r.get("p_x1", np.nan), r.get("p_post", np.nan),
                       r.get("rho_pre", np.nan), r.get("rho_x1", np.nan), r.get("rho_post", np.nan),
                       r.get("prms_x1", np.nan), r.get("rrms_x1", np.nan),
                       r["prms_med16"], r["rrms_med16"]] + xp[1:6] + sn[:5])
np.savez(HERE / "shock_cell_metrics.npz",
         **{r["tag"].replace("-", "_") + "_xpk": r["xpk"] for r in ALL},
         **{r["tag"].replace("-", "_") + "_Sn": r["Sn"] for r in ALL})
print("\nwrote shock_cell_metrics.csv / .npz")
