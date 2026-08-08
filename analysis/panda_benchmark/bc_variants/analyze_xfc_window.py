#!/usr/bin/env python3
"""How long an averaging window does x_fc actually need?

The statistics dumps carry RUNNING means accumulated since t0, so

  cumulative window [t0, t_i]  ->  M_i directly
  disjoint window  [t_i, t_j]  ->  [ (t_j-t0) M_j - (t_i-t0) M_i ] / (t_j-t_i)

x_fc is extracted with the thesis definition: the steepest rise of the mean
centreline density between the first expansion minimum and the first
compression maximum, refined parabolically. Identical code path for every
window, so the scatter is pure sampling noise.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
D = 0.0254
RHO_J = 1.6413
T0 = 0.01574728874          # stats restart time, from truestat.sbatch/chk32080


def refine(x, y, i):
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def xfc_of(x, rho):
    """Thesis definition: raw density-gradient steepest rise between the first
    expansion minimum and the first compression maximum."""
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > imin) & (x <= 1.7))[0]
    imax = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)
    seg = np.arange(imin, imax + 1)
    return refine(x, g, seg[np.nanargmax(g[seg])])


d = np.load(HERE / "XFC_window_L6.npz")
keys = sorted([k[:-2] for k in d.files if k.endswith("_t")],
              key=lambda s: int("".join(c for c in s if c.isdigit())))
dx = float(d[keys[0] + "_t"][1])
t = np.array([float(d[k + "_t"][0]) for k in keys])
M = np.array([d[k + "_rho"] for k in keys]) / RHO_J
x = (np.arange(M.shape[1]) + 0.5) * dx / D
W = t - T0                                     # accumulated window length

print("dx = %.2f um   dx/D_e = %.5f" % (dx*1e6, dx/D))
print("\nCUMULATIVE windows  [t0, t_i]")
print("%10s %10s %12s" % ("t [ms]", "window[ms]", "x_fc/D_e"))
xc = []
for k, ti, wi, m in zip(keys, t, W, M):
    v = xfc_of(x, m); xc.append(v)
    print("%10.4f %10.4f %12.5f" % (ti*1e3, wi*1e3, v))
xc = np.array(xc)
print("  spread over the last 4 (window 0.93-1.50 ms): %.5f D_e"
      % np.ptp(xc[-4:]))
print("  drift from window 0.20 ms to 1.50 ms        : %+.5f D_e"
      % (xc[-1] - xc[0]))

print("\nDISJOINT windows  [t_i, t_j]  (independent samples)")
print("%22s %10s %12s" % ("window [ms]", "length[ms]", "x_fc/D_e"))
xd, wd = [], []
for i in range(len(t) - 1):
    num = W[i+1]*M[i+1] - W[i]*M[i]
    seg = num / (t[i+1] - t[i])
    v = xfc_of(x, seg); xd.append(v); wd.append((t[i+1]-t[i])*1e3)
    print("%9.4f - %-9.4f %10.4f %12.5f"
          % (t[i]*1e3, t[i+1]*1e3, (t[i+1]-t[i])*1e3, v))
xd = np.array(xd); wd = np.array(wd)
sel = wd > 0.15                                # drop the short last window
print("  %d windows of ~%.3f ms:  mean %.5f  std %.5f  spread %.5f D_e"
      % (sel.sum(), wd[sel].mean(), xd[sel].mean(), xd[sel].std(ddof=1),
         np.ptp(xd[sel])))
print("\n  s-ladder signal 0.052 D_e ;  Arm B discrimination 0.020 D_e")
sd = xd[sel].std(ddof=1)
for wtar in (0.2, 0.3, 0.4):
    print("  a %.1f ms window -> sigma(x_fc) ~ %.5f D_e  (%.1f%% of the "
          "0.052 signal)" % (wtar, sd*np.sqrt(wd[sel].mean()/wtar),
                             100*sd*np.sqrt(wd[sel].mean()/wtar)/0.052))

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.2))
a1.plot(W*1e3, xc, "o-", color="#1f4e9c", ms=5, lw=1.2, mec="w", mew=0.6)
a1.axhline(xc[-1], color="0.55", lw=0.8, ls=":")
a1.fill_between([0, W[-1]*1e3], xc[-1]-0.001, xc[-1]+0.001, color="#1f4e9c",
                alpha=0.12, lw=0)
a1.set_xlabel("cumulative averaging window [ms]")
a1.set_ylabel("$x_{fc}/D_e$")
a1.set_xlim(0, W[-1]*1e3*1.03)
a1.grid(True, ls="--", lw=0.5, color="0.88"); a1.set_axisbelow(True)
a1.text(0.97, 0.06, "band = $\\pm$0.001 $D_e$\n(extraction resolution)",
        transform=a1.transAxes, ha="right", fontsize=6.8, color="#1f4e9c")
a1.text(0.03, 0.94, "(a) cumulative", transform=a1.transAxes, va="top",
        fontsize=8.5)

ctr = 0.5*(t[:-1] + t[1:])*1e3
a2.plot(ctr, xd, "s-", color="#b3251e", ms=5, lw=1.2, mec="w", mew=0.6)
a2.axhline(xd[sel].mean(), color="0.4", lw=0.9)
a2.fill_between([ctr[0], ctr[-1]], xd[sel].mean()-xd[sel].std(ddof=1),
                xd[sel].mean()+xd[sel].std(ddof=1), color="#b3251e",
                alpha=0.13, lw=0)
a2.set_xlabel("window centre [ms]")
a2.set_ylabel("$x_{fc}/D_e$")
a2.grid(True, ls="--", lw=0.5, color="0.88"); a2.set_axisbelow(True)
a2.text(0.03, 0.94, "(b) disjoint ~0.24 ms windows", transform=a2.transAxes,
        va="top", fontsize=8.5)
a2.text(0.97, 0.06, "$\\sigma$ = %.4f $D_e$" % xd[sel].std(ddof=1),
        transform=a2.transAxes, ha="right", fontsize=7.5, color="#b3251e")
fig.tight_layout(w_pad=1.3)
fig.savefig(HERE / "xfc_window.png", dpi=260)
fig.savefig(HERE / "xfc_window.pdf")
print("\nwrote xfc_window.png/pdf")
