#!/usr/bin/env python3
"""What the shift constant s actually does, and why it is the right knob
for an x_fc study.

The production inflow is  f = 0.5[1 - tanh((r - r0)/a)],  r0 = R - s a,
truncated at the physical lip r = R (flange slip wall beyond). Then

    delta_w = 2a                     (exactly independent of s)
    theta   = (a/4)(1 + tanh s)  ->  a/2   frozen to 0.3% for s >= 3
    delta*  = (a/2) ln(1 + e^{2s}) -> s a  grows linearly

so (a, s) set the shear-layer shape and the effective aperture
independently. Panel (a) shows the family at fixed a; panel (b) shows that
the aperture law and the shear-layer laws predict qualitatively different
behaviour along the s ladder - a slope versus a flat line.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

from design_2d_sweep import scalars, fvel, R, D, DX_LIP, EXIST

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
A = 170e-6
SS = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

y = np.array([e[3] for e in EXIST])
S4 = [scalars(e[1] * 1e-6, e[2]) for e in EXIST]
LAWS = []
for key, nm, col in [("dstar_c", r"aperture law  $x_{fc}(\delta^*)$", "#b3251e"),
                     ("dw", r"$\delta_\omega$ law", "#1f4e9c"),
                     ("theta_c", r"$\theta$ law", "#2ca02c")]:
    xv = np.array([s[key] for s in S4]) / D
    m, c = np.polyfit(xv, y, 1)
    LAWS.append((nm, key, col, m, c))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.3))

# ------------------------------------------------------------- panel (a)
RLO, RHI = 0.418, 0.5085
for k, v in enumerate(np.arange(0.0, R + DX_LIP, DX_LIP) / D):
    if RLO <= v <= RHI:
        ax.axhline(v, color="0.88", lw=0.35, zorder=0)
        if k % 2 == 0:
            ax.axhspan(v, min(v + DX_LIP / D, RHI), color="0.968", zorder=-1)
r = np.linspace(0.35 * D, R, 6000)
cm = plt.cm.viridis(np.linspace(0.03, 0.88, len(SS)))
for s, c in zip(SS, cm):
    f = fvel(r, A, "shifted", s)
    ax.plot(np.r_[f, 0.0], np.r_[r / D, R / D], color=c, lw=1.3, zorder=4,
            label="$s$=%g" % s)
fw = fvel(r, 340e-6, "wall")
ax.plot(np.r_[fw, 0.0], np.r_[r / D, R / D], color="k", lw=1.2, ls="--",
        zorder=5, label="wall tanh\n$a$=340 (same $\\delta_\\omega$)")
ax.axhline(R / D, color="k", lw=1.1, zorder=6)
ax.text(0.5, 0.972, "physical lip  $r=R$", transform=ax.transAxes,
        fontsize=7, ha="center")
ax.annotate("", xy=(0.50, 0.5 - 7 * A / D), xytext=(0.50, 0.5),
            arrowprops=dict(arrowstyle="<->", lw=0.8, color="0.3"))
ax.text(0.53, 0.5 - 3.5 * A / D, "$s\\,a$", fontsize=7.5, color="0.3")
ax.set_xlim(-0.03, 1.06); ax.set_ylim(RLO, RHI)
ax.set_xlabel("$u/U_e$"); ax.set_ylabel("$r/D_e$")
leg = ax.legend(fontsize=6.2, frameon=False, loc="lower left",
                handlelength=1.2, labelspacing=0.22, borderpad=0.2,
                title="$a$=170 $\\mu$m fixed  $\\Rightarrow$ $\\delta_\\omega$=340 $\\mu$m fixed")
leg.get_title().set_fontsize(6.2)
ax.text(0.965, 0.035, "(a)", transform=ax.transAxes, ha="right", fontsize=9)

# ------------------------------------------------------------- panel (b)
sl = np.array(SS[3:])
for (nm, key, col, m, c) in LAWS:
    v = np.array([c + m * scalars(A, "shifted", s)[key] / D for s in sl])
    ax2.plot(sl, v, "o-", color=col, ms=5, lw=1.3, mec="w", mew=0.6,
             label=nm, zorder=4)
ax2.axvline(3.0, color="0.55", lw=0.8, ls=":", zorder=1)
ax2.text(3.08, 0.9455, "production\nvalue $s$=3", fontsize=6.6, color="0.35",
         va="top")
ax2.set_xlabel("shift constant $s$   ($r_0 = R - s\\,a$,  $a$=170 $\\mu$m)")
ax2.set_ylabel("$x_{fc}/D_e$ predicted")
ax2.set_xlim(2.7, 7.3); ax2.set_ylim(0.842, 0.950)
ax2.grid(True, ls="--", lw=0.5, color="0.88"); ax2.set_axisbelow(True)
ax2.annotate("", xy=(6.85, 0.8510), xytext=(6.85, 0.9034),
             arrowprops=dict(arrowstyle="<->", lw=0.9, color="#b3251e"))
ax2.text(6.72, 0.877, "0.052 $D_e$", fontsize=7, color="#b3251e",
         rotation=90, va="center", ha="right")
ax2.text(5.0, 0.9105, "$\\delta_\\omega$ and $\\theta$ frozen $\\Rightarrow$ flat",
         fontsize=6.8, color="#1f4e9c", ha="center")
ax2.legend(fontsize=6.6, frameon=False, loc="lower left", handlelength=1.5,
           labelspacing=0.25)
ax2.text(0.965, 0.035, "(b)", transform=ax2.transAxes, ha="right", fontsize=9)

fig.subplots_adjust(left=0.085, right=0.985, top=0.975, bottom=0.155,
                    wspace=0.30)
fig.savefig(HERE / "shift_family.png", dpi=260)
fig.savefig(HERE / "shift_family.pdf")
print("wrote shift_family.png/pdf")
print("%3s %8s %8s %9s %7s %7s %10s" % ("s", "dw[um]", "th_c[um]", "d*_c[um]",
                                        "H", "C_d", "u(R)/Ue"))
for s in SS:
    x = scalars(A, "shifted", s)
    print("%3g %8.1f %8.1f %9.1f %7.2f %7.4f %10.2e"
          % (s, x["dw"] * 1e6, x["theta_c"] * 1e6, x["dstar_c"] * 1e6,
             x["H"], x["Cd"], x["leak"]))
