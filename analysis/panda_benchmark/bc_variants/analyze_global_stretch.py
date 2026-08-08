#!/usr/bin/env python3
"""Test whether the shock-cell trains of the four nozzle profiles differ
by a single global axial stretch factor alpha, i.e. whether every cell
spacing scales by the same alpha relative to the baseline.

For each profile: x_fc (density closure) and the compression maxima
x_1..x_6 (density landmarks). Define the distance from closure
d_n = x_n - x_fc. If the trains are self-similar under a uniform stretch,
d_n^profile = alpha * d_n^baseline for all n (linear through the origin).
The slope alpha and its R^2 quantify how well a single factor describes
the whole train; the per-cell spacing ratio S_n^profile/S_n^baseline
shows whether alpha is really n-independent.

Figure: (a) d_n^profile vs d_n^baseline with the fitted alpha lines;
(b) per-cell spacing ratio vs cell index."""
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RHO_J = 1.6413
CASES = [("baseline", "Baseline", "#7a7f85"),
         ("tophat", "Top hat", "#0b6ef5"),
         ("walltanh", "Wall tanh", "#2ca02c"),
         ("shift100", "Thin shifted tanh", "#ff6a00")]


def boxfilt(a, x, wD=0.05):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w // 2] = a[:w // 2]; out[-(w // 2):] = a[-(w // 2):]
    return out


def refine(x, y, i):
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    d = y[i - 1] - 2 * y[i] + y[i + 1]
    return x[i] if d == 0 else x[i] + 0.5 * (y[i - 1] - y[i + 1]) / d * (x[1] - x[0])


def peaks(v, nmax=6):
    d = np.load(HERE / f"AXIS8_{v}.npz")
    x = d["x"]; y = boxfilt(d["DensityMEAN"] / RHO_J, x)
    e = np.where((x >= 0.35) & (x <= 1.0))[0]; imin = e[np.argmin(y[e])]
    c = np.where((np.arange(x.size) > imin) & (x <= 1.65))[0]
    imax = c[np.argmax(y[c])]
    g = np.gradient(y, x)
    xfc = refine(x, g, np.arange(imin, imax + 1)[np.argmax(g[imin:imax + 1])])
    dist = max(1, int(round(0.5 / (x[1] - x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, y, -np.inf), prominence=0.02,
                       distance=dist)
    xn = np.array([refine(x, y, i) for i in pk[:nmax]])
    return xfc, xn


P = {v: peaks(v) for v, _, _ in CASES}
xfc_b, xn_b = P["baseline"]
db = xn_b - xfc_b                                    # baseline distances

print("%-18s %6s  alpha   R^2   per-cell S_n/S_n^base"
      % ("case", "x_fc"))
ALPHA = {}
for v, lab, c in CASES:
    xfc, xn = P[v]
    d = xn - xfc
    m = min(len(d), len(db))
    # alpha via least squares through the origin
    alpha = float(np.dot(d[:m], db[:m]) / np.dot(db[:m], db[:m]))
    pred = alpha * db[:m]
    ss = 1 - np.sum((d[:m] - pred)**2) / np.sum((d[:m] - d[:m].mean())**2)
    ALPHA[v] = (alpha, xfc, xn)
    Sb = np.diff(np.concatenate([[xfc_b], xn_b]))
    Sv = np.diff(np.concatenate([[xfc], xn]))
    ratios = (Sv[:m] / Sb[:m])
    print("%-18s %6.3f  %.3f  %.3f  %s"
          % (lab, xfc, alpha, ss,
             " ".join("%.2f" % r for r in ratios)))

# ---- figure ----
plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.6, 3.2))
for ax in (a1, a2):
    ax.set_axisbelow(True); ax.grid(True, ls="--", lw=0.5, color="0.88")
    ax.tick_params(length=2.5)
dd = np.linspace(0, db.max() * 1.05, 50)
a1.plot(dd, dd, color="0.6", lw=0.8, ls=":", zorder=1)
for v, lab, c in CASES:
    alpha, xfc, xn = ALPHA[v]
    d = xn - xfc
    m = min(len(d), len(db))
    a1.plot(db[:m], d[:m], "o", color=c, ms=5, mec="w", mew=0.5, zorder=4,
            label=r"%s ($\alpha$=%.3f)" % (lab, alpha))
    a1.plot(dd, alpha * dd, color=c, lw=1.0, zorder=2)
    Sb = np.diff(np.concatenate([[xfc_b], xn_b]))
    Sv = np.diff(np.concatenate([[xfc], xn]))
    a2.plot(np.arange(1, m + 1), (Sv[:m] / Sb[:m]), "o-", color=c, ms=5,
            mec="w", mew=0.5, lw=1.0, label=lab)
    a2.axhline(alpha, color=c, lw=0.7, ls="--", alpha=0.6)
a1.set_xlabel(r"$d_n^{\mathrm{baseline}} = (x_n-x_{fc})_{\mathrm{base}}$")
a1.set_ylabel(r"$d_n^{\mathrm{profile}}$")
a1.legend(fontsize=6.8, frameon=False, loc="upper left", handlelength=1.4)
a1.text(0.03, 0.03, "(a)", transform=a1.transAxes, fontsize=9)
a2.set_xlabel("cell index $n$")
a2.set_ylabel(r"$S_n^{\mathrm{profile}}/S_n^{\mathrm{baseline}}$")
a2.axhline(1.0, color="0.6", lw=0.8, ls=":")
a2.set_ylim(0.8, 1.35)
a2.legend(fontsize=6.8, frameon=False, loc="upper right", handlelength=1.4)
a2.text(0.03, 0.03, "(b)", transform=a2.transAxes, fontsize=9)
fig.tight_layout(w_pad=1.4)
fig.savefig(HERE / "analyze_global_stretch.png", dpi=250)
fig.savefig(HERE / "analyze_global_stretch.pdf")
print("\nwrote analyze_global_stretch.png/pdf")
