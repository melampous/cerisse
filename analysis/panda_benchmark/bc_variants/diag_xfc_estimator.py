#!/usr/bin/env python3
"""Is the x_fc scatter physical unsteadiness of the first cell, or an unstable
estimator?

The thesis definition takes the RAW density gradient and picks its argmax
between the first expansion minimum and the first compression maximum. If the
gradient has competing local maxima, that argmax can jump by many cells while
the underlying profile barely moves. This compares the raw estimator against
progressively smoothed ones on the same disjoint windows.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
D, RHO_J = 0.0254, 1.6413
T0 = 0.01574728874


def boxfilt(a, x, wD):
    if wD <= 0:
        return a.copy()
    w = max(3, int(round(wD / (x[1]-x[0]))))
    if w % 2 == 0:
        w += 1
    o = np.convolve(a, np.ones(w)/w, mode="same")
    o[:w//2] = a[:w//2]; o[-(w//2):] = a[-(w//2):]
    return o


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def xfc_of(x, rho, sm=0.0):
    r = boxfilt(rho, x, sm)
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > imin) & (x <= 1.7))[0]
    imax = c[np.nanargmax(r[c])]
    g = np.gradient(r, x)
    seg = np.arange(imin, imax+1)
    return refine(x, g, seg[np.nanargmax(g[seg])]), r, g, imin, imax


d = np.load(HERE / "XFC_window_L6.npz")
keys = sorted([k[:-2] for k in d.files if k.endswith("_t")],
              key=lambda s: int("".join(c for c in s if c.isdigit())))
dx = float(d[keys[0]+"_t"][1])
t = np.array([float(d[k+"_t"][0]) for k in keys])
M = np.array([d[k+"_rho"] for k in keys]) / RHO_J
x = (np.arange(M.shape[1])+0.5)*dx/D
W = t - T0

segs, ctr = [], []
for i in range(len(t)-1):
    if t[i+1]-t[i] < 0.15e-3:
        continue
    segs.append((W[i+1]*M[i+1] - W[i]*M[i])/(t[i+1]-t[i]))
    ctr.append(0.5*(t[i]+t[i+1])*1e3)
segs = np.array(segs)
full = M[-1]

print("estimator stability on %d disjoint ~0.243 ms windows" % len(segs))
print("%14s %9s %9s %9s   %s" % ("smoothing", "mean", "std", "spread",
                                 "per-window values"))
best = None
for sm in (0.0, 0.02, 0.04, 0.06, 0.08, 0.12):
    v = np.array([xfc_of(x, s, sm)[0] for s in segs])
    print("%10.2f D_e %9.4f %9.5f %9.5f   %s"
          % (sm, v.mean(), v.std(ddof=1), np.ptp(v),
             " ".join("%.4f" % q for q in v)))
    if best is None or v.std(ddof=1) < best[1]:
        best = (sm, v.std(ddof=1), v)
print("\nfull 1.503 ms window, x_fc vs smoothing:")
for sm in (0.0, 0.02, 0.04, 0.06, 0.08, 0.12):
    print("   %.2f D_e -> %.5f" % (sm, xfc_of(x, full, sm)[0]))

sm_b, sd_b, _ = best
print("\nbest estimator: %.2f D_e smoothing, sigma = %.5f D_e per 0.243 ms"
      % (sm_b, sd_b))
for tgt in (0.005, 0.002):
    print("   window needed for sigma <= %.3f D_e : %.2f ms"
          % (tgt, 0.243*(sd_b/tgt)**2))

# --------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.3))
cm = plt.cm.plasma(np.linspace(0.05, 0.85, len(segs)))
m = (x > 0.3) & (x < 1.6)
for s, c, cc in zip(segs, cm, ctr):
    a1.plot(x[m], s[m], color=c, lw=1.0, label="%.2f ms" % cc)
    g = np.gradient(s, x)
    a2.plot(x[m], g[m], color=c, lw=1.0)
    v = xfc_of(x, s, 0.0)[0]
    a2.axvline(v, color=c, lw=0.7, ls="--")
a1.plot(x[m], full[m], color="k", lw=1.6, label="full 1.50 ms")
a2.plot(x[m], np.gradient(full, x)[m], color="k", lw=1.6)
a1.set_xlabel("$x/D_e$"); a1.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
a2.set_xlabel("$x/D_e$")
a2.set_ylabel(r"$\mathrm{d}\langle\rho\rangle/\mathrm{d}(x/D_e)\,/\rho_j$")
a1.legend(fontsize=6.5, frameon=False, loc="upper left", handlelength=1.3,
          labelspacing=0.25)
for a, lab in ((a1, "(a) mean density"), (a2, "(b) gradient, dashed = raw $x_{fc}$")):
    a.grid(True, ls="--", lw=0.5, color="0.88"); a.set_axisbelow(True)
    a.text(0.97, 0.05, lab, transform=a.transAxes, ha="right", fontsize=7.5)
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "diag_xfc_estimator.png", dpi=260)
fig.savefig(HERE / "diag_xfc_estimator.pdf")
print("\nwrote diag_xfc_estimator.png/pdf")
