#!/usr/bin/env python3
"""wall-tanh (a=200 um) three-block acceptance + shock-cell metrics.
Block means recovered from cumulative per-step statistics with
tau-weighting (tau_k = t_k - t_s, t_s = 8.0 ms):
  q_blk(k) = [tau_k q_cum(t_k) - tau_{k-1} q_cum(t_{k-1})] / (tau_k - tau_{k-1})
Peaks: same method as shock_cell_metrics (smoothed <p>, prominence 0.03,
height > p_inf, distance 0.5D, parabolic refinement).
Acceptance: range_k x_n <= u_x = max(0.01D, dx/2); no monotone drift.
Also compares full-window x_n vs tophat-L5 and baseline-L5."""
import numpy as np
from scipy.signal import find_peaks
from pathlib import Path

HERE = Path(__file__).parent
D = 0.0254
P_AMB = 99780.0
T_S = 8.0e-3
FR = [("AXW_b1.npz", None), ("AXW_b2.npz", None), ("AXW_b3.npz", None)]

def peaks_of(x, pn, dxg):
    m = (x >= 0.45) & (x <= 10.0)
    xs, ps = x[m], pn[m]
    w = max(3, int(round(0.05 / dxg)))
    w = min(w, 51)
    k = np.ones(w) / w
    psm = np.convolve(ps, k, mode="same")
    dist = max(1, int(round(0.5 / (xs[1] - xs[0]))))
    psm2 = psm.copy(); psm2[:w] = -np.inf; psm2[-w:] = -np.inf
    pk, _ = find_peaks(psm2, prominence=0.03, distance=dist, height=1.0)
    out = []
    for j in pk:
        if 1 <= j < len(xs) - 1:
            a, b, c = psm[j-1], psm[j], psm[j+1]
            den = a - 2*b + c
            off = 0.5*(a - c)/den if abs(den) > 1e-12 else 0.0
            out.append(xs[j] + off*(xs[1]-xs[0]))
    return np.array(out[:5])

# cumulative fields and taus
cum = []
taus = []
for fn, _ in FR:
    d = np.load(HERE / fn, allow_pickle=True)
    x = d["x"] / D
    cum.append(d["ax_pressureMEAN"])
    taus.append(float(d["time"]) - T_S)
taus = np.array(taus)
dxg = x[1] - x[0]

# block means (tau-weighted differencing)
blocks = [cum[0]]
for k in (1, 2):
    blocks.append((taus[k]*cum[k] - taus[k-1]*cum[k-1]) / (taus[k] - taus[k-1]))

u_x = max(0.01, dxg / 2)
print(f"tau_k [ms]: {taus*1e3}")
print(f"u_x = {u_x:.4f} D\n")
XB = []
for k, blk in enumerate(blocks):
    xp = peaks_of(x, blk / P_AMB, dxg)
    XB.append(xp)
    print(f"block {k+1} [{(T_S+([0]+list(taus))[k])*1e3:.2f},{(T_S+taus[k])*1e3:.2f}]ms  "
          f"x_n: " + " ".join(f"{v:6.3f}" for v in xp))
n_ok = min(len(v) for v in XB)
print(f"\nacceptance on x_1..x_{n_ok}:")
ok_all = True
for n in range(n_ok):
    v = [XB[k][n] for k in range(3)]
    rng = max(v) - min(v)
    mono = (v[0] < v[1] < v[2]) or (v[0] > v[1] > v[2])
    verdict = "PASS" if (rng <= u_x and not (mono and rng > 0.5*u_x)) else "FAIL"
    if verdict == "FAIL":
        ok_all = False
    print(f"  x_{n+1}: {v[0]:.3f} {v[1]:.3f} {v[2]:.3f}  range={rng:.4f}D "
          f"monotone={mono}  {verdict}")
print("\nOVERALL:", "PASS" if ok_all else "FAIL")

# full-window peaks vs tophat / baseline
d3 = np.load(HERE / "AXW_b3.npz", allow_pickle=True)
xw = peaks_of(x, d3["ax_pressureMEAN"] / P_AMB, dxg)
mm = np.load(HERE / "shock_cell_metrics.npz", allow_pickle=True)
print("\nfull-window x_n comparison (x/D):")
print("wall-tanh :", " ".join(f"{v:6.3f}" for v in xw))
tt = np.load(HERE / "STATS3D_tophat.npz" if (HERE/"STATS3D_tophat.npz").exists()
             else HERE.parent / "STATS3D_tophat.npz", allow_pickle=True)
print("(baseline/tophat L5 from earlier metrics: 3D-L5 x_n:",
      " ".join(f"{v:6.3f}" for v in mm["3D_L5_xpk"]), ")")
np.savez(HERE / "WALLTANH_metrics.npz", x=x, taus=taus,
         xpk_blocks=np.array([np.pad(v.astype(float), (0, 5-len(v)),
                                     constant_values=np.nan) for v in XB]),
         xpk_full=xw)
print("\nwrote WALLTANH_metrics.npz")
