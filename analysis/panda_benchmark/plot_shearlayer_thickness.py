#!/usr/bin/env python3
"""Where each formulation spends its dissipation: shear layer or compressions.

The envelope figure ranks the five formulations by how fast the shock-cell
train decays. It does not say through which channel the dissipation acts. Two
channels are possible and they are distinguishable:

  * indirectly, by thickening the mean shear layer, which weakens the
    reflection that closes each cell and so shortens the train;
  * directly, by damping the compressions themselves.

The vorticity thickness of the mean axial velocity separates them. If a
formulation decays the cell train because its shear layer is thick, the two
quantities move together and the cases fall on one line. A formulation that
damps the compressions directly leaves that line: it decays the train while
keeping a thin shear layer.

    delta_omega(x) = [ max_y <u> - u_amb ] / max_y |d<u>/dy|

evaluated separately on the two sides of the axis and averaged, then smoothed
over x +/- 0.5 D_e. Without that smoothing the shock cells modulate max_y <u>
and the curve oscillates by tens of per cent, which is the shock-cell signal
and not the shear layer.

The y-derivative uses the run-length operator, not np.gradient. The slice is a
finest-wins AMR composite, so away from the finest patches one coarse cell is
replicated across 2, 4 or 8 output rows; a plain difference returns zero inside
each replicated block and the whole jump at its edge, which destroys the
maximum this metric is built on. Each line is run-length encoded, the
derivative is taken between the centres of neighbouring runs and broadcast
back, so the operator always uses the cell that owns the data. The same
operator is used for all five, including the two cases whose files carry an
rfmap, so no case benefits from information the others lack.

Panda measured density by Rayleigh scattering, not velocity, so panel (a) has
no experimental anchor. This figure compares formulations with each other.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

HERE = Path(__file__).resolve().parent
FL, SW = HERE/"bc_variants"/"fields", HERE/"bc_variants"/"sweep_out"
RHO_J, D = 1.6413, 0.0254
C_LLF, C_HLLC, C_SKEW = "#0072B2", "#D55E00", "#009E73"

CASES = [
    dict(k="llfweno",  lab="LLF + WENO-Z5",      c=C_LLF,  ls="-",  lw=1.7,
         flux=0, rec=0, u=FL/"SLM_llfweno_x_velocityMEAN.npz"),
    dict(k="llfteno",  lab="LLF + TENO5",        c=C_LLF,  ls="--", lw=1.7,
         flux=0, rec=1, u=FL/"SLM_llfteno_x_velocityMEAN.npz"),
    dict(k="hllcweno", lab="AFD–HLLC + WENO-Z5", c=C_HLLC, ls="-",  lw=1.7,
         flux=1, rec=0, u=FL/"SLM_hllcweno_x_velocityMEAN.npz"),
    dict(k="hllcteno", lab="AFD–HLLC + TENO5",   c=C_HLLC, ls="--", lw=1.7,
         flux=1, rec=1, u=FL/"SLM_hllcteno_x_velocityMEAN.npz"),
    dict(k="skewjst",  lab="Skew-symmetric + JST", c=C_SKEW, ls="-", lw=2.3,
         flux=None, rec=None, u=SW/"SL_L5_s3d_skewjst_x_velocityMEAN.npz"),
]
# x_half from plot_dissipation_envelope.py, recomputed there from the axis data
XHALF = dict(llfweno=5.989, llfteno=5.200, hllcweno=6.855, hllcteno=8.035,
             skewjst=4.977)


def grad_rf(a, h, axis):
    """Derivative over the cell that owns the data, on an AMR composite."""
    a = a if axis == 1 else a.T
    out = np.zeros_like(a)
    n = a.shape[1]
    for k in range(a.shape[0]):
        row = a[k]
        e = np.r_[0, np.flatnonzero(np.diff(row) != 0)+1, n]
        if e.size < 3:
            continue
        c = 0.5*(e[:-1]+e[1:]-1)*h
        v = row[e[:-1]]
        g = np.gradient(v, c) if v.size > 2 else np.diff(v)/np.diff(c)
        out[k] = np.repeat(g if v.size > 2 else np.r_[g, g[-1]], np.diff(e))
    return out if axis == 1 else out.T


def delta_omega(path):
    z = np.load(path)
    u = z["field"].astype(np.float64)          # (x, y)
    x, y = z["x"], z["y"]
    h = float(z["dx_um"])*1e-6/D               # cell size in D_e
    dudy = grad_rf(u, h, 1)
    amb = np.median(u[np.argmin(np.abs(x-3.0)), np.abs(y) > 1.4])
    du = u.max(axis=1) - amb
    out = np.empty_like(x)
    for side in (1, -1):
        m = (side*y > 0.02) & (np.abs(y) < 1.30)
        g = np.abs(dudy[:, m]).max(axis=1)
        out = du/g if side == 1 else 0.5*(out + du/g)
    # smooth over +/- 0.5 D_e so the shock cells do not modulate the metric
    n = int(round(1.0/(x[1]-x[0])))
    n = n if n % 2 else n+1
    sm = np.convolve(out, np.ones(n)/n, mode="same")
    sm[:n//2] = out[:n//2]
    sm[-(n//2):] = out[-(n//2):]
    return x, sm, amb, h


print("%-22s %8s %8s %8s %8s %8s %8s %9s"
      % ("scheme", "u_amb", "dw(0.5)", "dw(1)", "dw(3)", "dw(4.5)", "dw(6)",
         "cells@6D"))
for c in CASES:
    c["x"], c["dw"], amb, h = delta_omega(c["u"])
    g = lambda q: float(np.interp(q, c["x"], c["dw"]))
    c["dw6"] = g(6.0)
    print("%-22s %8.2f %8.4f %8.4f %8.4f %8.4f %8.4f %9.1f"
          % (c["lab"], amb, g(0.5), g(1.0), g(3.0), g(4.5), g(6.0),
             g(6.0)/(4*h)))          # local cell at 6 D_e is 4x the finest

G = {(c["flux"], c["rec"]): c["dw6"] for c in CASES if c["flux"] is not None}
fx = 0.5*((G[(1, 0)]+G[(1, 1)]) - (G[(0, 0)]+G[(0, 1)]))
rc = 0.5*((G[(0, 1)]+G[(1, 1)]) - (G[(0, 0)]+G[(1, 0)]))
ix = 0.5*((G[(1, 1)]-G[(1, 0)]) - (G[(0, 1)]-G[(0, 0)]))
print("\ntwo-by-two decomposition of delta_omega(6 D_e) [D_e]")
print("  flux path   LLF -> AFD-HLLC     %+.4f" % fx)
print("  reconstruction WENO-Z5 -> TENO5 %+.4f" % rc)
print("  interaction                     %+.4f" % ix)
print("  -> flux path is %.0f times the reconstruction main effect"
      % abs(fx/rc))

# ---- does the 2x2 set fall on a line, and does skew sit off it?
q = [c for c in CASES if c["flux"] is not None]
X = np.array([c["dw6"] for c in q])
Y = np.array([XHALF[c["k"]] for c in q])
p = np.polyfit(X, Y, 1)
res = Y - np.polyval(p, X)
sk = [c for c in CASES if c["flux"] is None][0]
sk_res = XHALF[sk["k"]] - np.polyval(p, sk["dw6"])
print("\nx_half against delta_omega(6 D_e)")
print("  fit through the four two-by-two cases: slope %.2f D_e per unit"
      " delta_omega, rms residual %.3f D_e" % (p[0], np.sqrt((res**2).mean())))
for c, r in zip(q, res):
    print("    %-22s dw6 %.3f  x_half %.3f  residual %+.3f"
          % (c["lab"], c["dw6"], XHALF[c["k"]], r))
print("    %-22s dw6 %.3f  x_half %.3f  residual %+.3f  (%.1f x the fit rms)"
      % (sk["lab"], sk["dw6"], XHALF[sk["k"]], sk_res,
         abs(sk_res)/np.sqrt((res**2).mean())))
print("\n  two inequalities, no regression needed:")
thin = min(c["dw6"] for c in q if c["flux"] == 0)
print("    skew shear layer at 6 D_e is %.1f %% thinner than the thinner LLF"
      " case (%.3f vs %.3f)" % (100*(1-sk["dw6"]/thin), sk["dw6"], thin))
print("    yet its x_half is the shortest of the five (%.2f vs %.2f)"
      % (XHALF[sk["k"]], min(XHALF[c["k"]] for c in q)))

plt.rcParams.update({"font.size": 9.8, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.8, 4.5),
                             gridspec_kw=dict(width_ratios=[1.55, 1.0]))
fig.subplots_adjust(left=0.062, right=0.985, top=0.912, bottom=0.128,
                    wspace=0.215)

for c in CASES:
    m = (c["x"] >= 0.3) & (c["x"] <= 6.5)
    a1.plot(c["x"][m], c["dw"][m], color=c["c"], ls=c["ls"], lw=c["lw"],
            zorder=5, label="%s   %.3f" % (c["lab"], c["dw6"]))
a1.axvline(6.0, color="0.6", lw=0.9, ls=":", zorder=3)
a1.set_xlim(0.3, 6.5)
a1.set_ylim(0, 0.95)
a1.set_xlabel("$x/D_e$")
a1.set_ylabel(r"$\delta_\omega/D_e$")
a1.set_title(r"(a)  vorticity thickness of the mean axial velocity"
             r"   [value at $6D_e$]", fontsize=9.8, pad=5)
a1.legend(loc="upper left", fontsize=8.4, handlelength=2.2)

xx = np.linspace(min(X.min(), sk["dw6"])-0.03, X.max()+0.03, 50)
a2.plot(xx, np.polyval(p, xx), "-", color="0.55", lw=1.2, zorder=3,
        label="fit through the four (slope %.1f)" % p[0])
for c in q:
    a2.plot(c["dw6"], XHALF[c["k"]], "o", ms=10, mfc=c["c"], mec="k", mew=0.8,
            zorder=6, alpha=1.0 if c["rec"] == 0 else 0.5)
a2.plot(sk["dw6"], XHALF[sk["k"]], "s", ms=11, mfc=C_SKEW, mec="k", mew=0.9,
        zorder=7)
a2.annotate("", xy=(sk["dw6"], XHALF[sk["k"]]),
            xytext=(sk["dw6"], np.polyval(p, sk["dw6"])),
            arrowprops=dict(arrowstyle="<->", color="#009E73", lw=1.2),
            zorder=5)
a2.text(sk["dw6"]+0.012, 0.5*(XHALF[sk["k"]]+np.polyval(p, sk["dw6"])),
        "%.2f $D_e$ below\nthe line" % abs(sk_res), fontsize=8.4,
        color="#00805c", va="center")
for c in CASES:
    a2.annotate(c["lab"].replace(" + ", "+").replace("Skew-symmetric", "skew"),
                (c["dw6"], XHALF[c["k"]]), textcoords="offset points",
                xytext=(0, 13), fontsize=8.0, ha="center")
a2.set_xlim(0.53, 0.90)
a2.set_ylim(4.4, 8.9)
a2.set_xlabel(r"$\delta_\omega(6D_e)/D_e$")
a2.set_ylabel(r"$x_{1/2}/D_e$")
a2.set_title("(b)  decay against shear-layer thickness", fontsize=9.8, pad=5)
a2.legend(loc="upper right", fontsize=8.4, handlelength=1.8)
for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"shearlayer_thickness.png", dpi=250)
fig.savefig(HERE/"shearlayer_thickness.pdf")
print("\nwrote shearlayer_thickness.png/pdf")
