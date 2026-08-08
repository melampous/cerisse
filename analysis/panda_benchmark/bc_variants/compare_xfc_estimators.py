#!/usr/bin/env python3
"""Robust replacements for the argmax-of-raw-gradient x_fc estimator.

The centreline is a finest-wins composite across AMR levels, so it carries
level-jump steps that a finite-difference gradient turns into spikes. Taking
the argmax of that gradient is the single most fragile way to locate the
closure. These alternatives never differentiate, or differentiate only after
the location has been fixed by an integral quantity:

  E1 argmax of the raw gradient                (current thesis definition)
  E2 argmax of the gradient, 0.08 D_e smoothed
  E3 midpoint crossing: x where rho = (rho_min+rho_max)/2 on the rising branch
  E4 gradient centroid over the rise, int x g dx / int g dx  (g clipped >= 0)
  E5 least-squares tanh fit to the rising branch, centre of the fit

Each is applied identically to the same disjoint windows, so the spread is a
direct measure of the estimator's sampling noise.
"""
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
D, RHO_J = 0.0254, 1.6413
T0 = 0.01574728874


def boxfilt(a, x, wD):
    if wD <= 0:
        return a.copy()
    w = max(3, int(round(wD/(x[1]-x[0]))))
    if w % 2 == 0:
        w += 1
    o = np.convolve(a, np.ones(w)/w, mode="same")
    o[:w//2] = a[:w//2]; o[-(w//2):] = a[-(w//2):]
    return o


def bracket(x, r):
    """First expansion minimum and first compression maximum."""
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = e[np.nanargmin(r[e])]
    c = np.where((np.arange(x.size) > imin) & (x <= 1.7))[0]
    imax = c[np.nanargmax(r[c])]
    return imin, imax


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def E1(x, r):
    i0, i1 = bracket(x, r); g = np.gradient(r, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def E2(x, r):
    rs = boxfilt(r, x, 0.08)
    i0, i1 = bracket(x, rs); g = np.gradient(rs, x)
    s = np.arange(i0, i1+1)
    return refine(x, g, s[np.nanargmax(g[s])])


def E3(x, r):
    i0, i1 = bracket(x, r)
    half = 0.5*(r[i0] + r[i1])
    seg = r[i0:i1+1]
    k = np.where(seg >= half)[0]
    if not len(k):
        return np.nan
    k = k[0] + i0
    if k == i0:
        return x[i0]
    x0, x1, y0, y1 = x[k-1], x[k], r[k-1], r[k]
    return x0 + (half-y0)*(x1-x0)/(y1-y0)


def E4(x, r):
    i0, i1 = bracket(x, r)
    g = np.clip(np.gradient(r, x)[i0:i1+1], 0, None)
    xs = x[i0:i1+1]
    return float(np.trapezoid(xs*g, xs)/np.trapezoid(g, xs))


def _tanh(x, xc, w, a, b):
    return a + b*np.tanh((x-xc)/w)


def E5(x, r):
    i0, i1 = bracket(x, r)
    xs, ys = x[i0:i1+1], r[i0:i1+1]
    a0 = 0.5*(ys[0]+ys[-1]); b0 = 0.5*(ys[-1]-ys[0])
    try:
        p, _ = curve_fit(_tanh, xs, ys, p0=[0.5*(xs[0]+xs[-1]), 0.05, a0, b0],
                         maxfev=20000)
        return float(p[0])
    except Exception:
        return np.nan


EST = [("E1 raw gradient argmax", E1), ("E2 smoothed argmax", E2),
       ("E3 midpoint crossing", E3), ("E4 gradient centroid", E4),
       ("E5 tanh fit centre", E5)]

d = np.load(HERE / "XFC_window_L6.npz")
keys = sorted([k[:-2] for k in d.files if k.endswith("_t")],
              key=lambda s: int("".join(c for c in s if c.isdigit())))
dx = float(d[keys[0]+"_t"][1])
t = np.array([float(d[k+"_t"][0]) for k in keys])
M = np.array([d[k+"_rho"] for k in keys])/RHO_J
x = (np.arange(M.shape[1])+0.5)*dx/D
W = t - T0

segs, ctr = [], []
for i in range(len(t)-1):
    if t[i+1]-t[i] < 0.15e-3:
        continue
    segs.append((W[i+1]*M[i+1]-W[i]*M[i])/(t[i+1]-t[i]))
    ctr.append(0.5*(t[i]+t[i+1])*1e3)
segs = np.array(segs)
full = M[-1]
TW = 0.243

print("%-24s %8s %9s %9s %9s   %s"
      % ("estimator", "full", "mean", "sigma", "spread", "per-window"))
res = {}
for nm, f in EST:
    v = np.array([f(x, s) for s in segs])
    res[nm] = v
    print("%-24s %8.4f %9.4f %9.5f %9.5f   %s"
          % (nm, f(x, full), np.nanmean(v), np.nanstd(v, ddof=1),
             np.ptp(v[np.isfinite(v)]),
             " ".join("%.4f" % q for q in v)))

print("\nwindow length needed, per estimator (sigma ~ 1/sqrt(T))")
print("%-24s %12s %12s %12s"
      % ("estimator", "s=0.005 De", "s=0.002 De", "0.8 ms gives"))
for nm, f in EST:
    sd = np.nanstd(res[nm], ddof=1)
    print("%-24s %9.2f ms %9.2f ms %9.5f De"
          % (nm, TW*(sd/0.005)**2, TW*(sd/0.002)**2, sd*np.sqrt(TW/0.8)))

print("\nsignals to resolve:  s-ladder 0.052 De ;  Arm B pairs 0.020 De")

# --------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.3),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})
m = (x > 0.55) & (x < 1.35)
cm = plt.cm.plasma(np.linspace(0.05, 0.85, len(segs)))
for s, c in zip(segs, cm):
    a1.plot(x[m], s[m], color=c, lw=0.9)
a1.plot(x[m], full[m], color="k", lw=1.7, label="full 1.50 ms")
cols = ["#7a7f85", "#0b6ef5", "#b3251e", "#2ca02c", "#ff6a00"]
for (nm, f), cc in zip(EST, cols):
    a1.axvline(f(x, full), color=cc, lw=1.1, ls="--")
a1.set_xlabel("$x/D_e$"); a1.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
a1.grid(True, ls="--", lw=0.5, color="0.88"); a1.set_axisbelow(True)
a1.text(0.03, 0.94, "(a) disjoint windows + estimator locations",
        transform=a1.transAxes, va="top", fontsize=7.5)

pos = np.arange(len(EST))
for j, ((nm, f), cc) in enumerate(zip(EST, cols)):
    v = res[nm]
    a2.plot(np.full(len(v), j), v, "o", color=cc, ms=5, mec="w", mew=0.6,
            alpha=0.85)
    sd = np.nanstd(v, ddof=1)
    a2.errorbar(j, np.nanmean(v), yerr=sd, color=cc, capsize=4, lw=1.4,
                marker="_", ms=16)
    a2.text(j, 0.955, "%.4f" % sd, ha="center", fontsize=6.6, color=cc)
a2.set_xticks(pos)
a2.set_xticklabels([n.split()[0] for n in [e[0] for e in EST]], fontsize=8)
a2.set_ylabel("$x_{fc}/D_e$ per 0.24 ms window")
a2.set_ylim(0.78, 0.97)
a2.grid(True, axis="y", ls="--", lw=0.5, color="0.88")
a2.set_axisbelow(True)
a2.text(0.03, 0.05, "(b) $\\sigma$ printed above each", transform=a2.transAxes,
        fontsize=7.5)
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "compare_xfc_estimators.png", dpi=260)
fig.savefig(HERE / "compare_xfc_estimators.pdf")
print("\nwrote compare_xfc_estimators.png/pdf")
