#!/usr/bin/env python3
"""A family of inlet profiles at fixed displacement thickness.

The shifted tanh has two parameters, and they enter the integral scales
differently:

    delta_omega = 2a                 set by a alone
    theta       ~ a/2                set by a alone
    delta*      ~ s*a + O(a)         set by both

so the locus delta* = const is a curve in the (a, s) plane along which a - and
with it the peak shear and the momentum deficit - is free to move. Walking that
curve gives a set of profiles that share a mass deficit and a discharge
coefficient while their shear differs by a factor of nearly three, with the
functional form never changing. That is the complement of the s-ladder, which
holds a fixed and lets delta* run.

Two things bound how far the curve can be walked, and both are drawn:

    lip leak     f(R-) = (1 - tanh s)/2 grows as s falls; past about 1 % the
                 profile no longer closes on the lip and the boundary
                 condition stops representing a solid edge
    resolution   the production grid resolves the lip at 24.8 um, so
                 delta_omega below ~200 um is carried by fewer than eight cells

Panel (c) shows the locus for several targets, so a different anchor can be
read straight off it.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
GAM, RGAS = 1.4, 8.31446261815324/28.96e-3
CP = GAM*RGAS/(GAM-1.0)
T_AMB, R = 297.15, 0.0127
T_E = T_AMB/1.2
U_E = np.sqrt(GAM*RGAS*T_E)
DX_LIP = 1.5875e-3/64
RHO_E = 1.0/(RGAS*T_E)
r = np.linspace(0.0, R, 200001)
rmm, Rmm = r*1e3, R*1e3


def prof(a, s):
    return 0.5*(1.0 - np.tanh((r - (R - s*a))/a))


def sc(a, s):
    f = prof(a, s)
    u = U_E*f
    rho = 1.0/(RGAS*(T_AMB - 0.5*u*u/CP))
    g = rho*u/(RHO_E*U_E)
    return dict(dw=1.0/np.max(np.abs(np.gradient(f, r))),
                thc=np.trapezoid(g*(1.0-f), r),
                dsc=np.trapezoid(1.0-g, r),
                Cd=2.0*np.trapezoid(g*r, r)/R**2,
                leak=0.5*(1.0 - np.tanh(s)))


def solve_s(a, target):
    lo, hi = 0.15, 60.0
    for _ in range(60):
        m = 0.5*(lo + hi)
        if sc(a, m)["dsc"] < target:
            lo = m
        else:
            hi = m
    return 0.5*(lo + hi)


TARGET = sc(170e-6, 3.0)["dsc"]
MEMBERS = [80, 100, 130, 170, 220]          # um; 80 is the marginal one
FAM = []
print("proposed family at d*_c = %.2f um\n" % (TARGET*1e6))
print("%8s %8s %9s %10s %10s %8s %10s %7s"
      % ("a[um]", "s", "dw[um]", "th_c[um]", "d*_c[um]", "C_d", "f(R-)", "cells"))
for a_um in MEMBERS:
    a = a_um*1e-6
    s = solve_s(a, TARGET)
    q = sc(a, s)
    FAM.append((a_um, s, q))
    print("%8.0f %8.3f %9.1f %10.1f %10.2f %8.4f %10.2e %7.1f"
          % (a_um, s, q["dw"]*1e6, q["thc"]*1e6, q["dsc"]*1e6, q["Cd"],
             q["leak"], q["dw"]/DX_LIP))
print("\n   delta_w  x%.2f      theta_c  x%.2f      d*_c  x%.5f      C_d spread %.4f"
      % (FAM[-1][2]["dw"]/FAM[0][2]["dw"], FAM[-1][2]["thc"]/FAM[0][2]["thc"],
         FAM[-1][2]["dsc"]/FAM[0][2]["dsc"],
         max(q["Cd"] for _, _, q in FAM) - min(q["Cd"] for _, _, q in FAM)))

COL = plt.get_cmap("viridis")(np.linspace(0.08, 0.86, len(FAM)))

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13.4, 5.0))
fig.subplots_adjust(left=0.055, right=0.988, top=0.925, bottom=0.115,
                    wspace=0.255)

for (a_um, s, q), col in zip(FAM, COL):
    f = prof(a_um*1e-6, s)
    g = np.abs(np.gradient(f, rmm))
    lab = "$a$=%d, $s$=%.2f" % (a_um, s)
    a1.plot(f, rmm, color=col, lw=1.8, zorder=5, label=lab)
    a2.plot(rmm, g, color=col, lw=1.8, zorder=5)

a1.axhspan(Rmm, 20.0, color="0.915", zorder=0)
a1.axhline(Rmm, color="0.4", lw=0.9, zorder=2)
a1.set_xlim(-0.02, 1.06)
a1.set_ylim(10.9, 12.95)
a1.set_xlabel("profile function  $f = u/U_e$")
a1.set_ylabel("radius  $r$  [mm]")
a1.set_title(r"(a)  the family,  all at $\delta^*_c$ = %.0f $\mu$m" % (TARGET*1e6),
             fontsize=10.5, pad=6)
a1.legend(loc="lower left", fontsize=9.0, handlelength=1.8, labelspacing=0.35)

a2.axvspan(Rmm, 20.0, color="0.915", zorder=0)
a2.axvline(Rmm, color="0.4", lw=0.9, zorder=2)
a2.set_xlim(10.9, 12.95)
a2.set_xlabel("radius  $r$  [mm]")
a2.set_ylabel(r"$|\mathrm{d}f/\mathrm{d}r|$   [1/mm]")
a2.set_title(r"(b)  peak shear spans $\times$%.1f"
             % (FAM[-1][2]["dw"]/FAM[0][2]["dw"]), fontsize=10.5, pad=6)

# ------------------------------------------------- (c) the (a, s) plane
aa = np.linspace(50e-6, 420e-6, 60)
for tg, lw in ((300e-6, 1.0), (400e-6, 1.0), (TARGET, 2.0),
               (700e-6, 1.0), (900e-6, 1.0)):
    ss = [solve_s(a, tg) for a in aa]
    a3.plot(aa*1e6, ss, color="0.25" if lw > 1.5 else "0.6", lw=lw, zorder=4)
    j = 6
    a3.text(aa[j]*1e6, ss[j], "%.0f" % (tg*1e6), fontsize=8.6,
            color="0.25" if lw > 1.5 else "0.55", ha="left", va="bottom")

s_leak = np.arctanh(1 - 2*0.01)
a3.axhspan(0, s_leak, color="#D55E00", alpha=0.11, zorder=0)
a3.axhline(s_leak, color="#D55E00", lw=1.0, ls="--", zorder=3)
a3.text(415, s_leak-0.12, "lip leak > 1 %", ha="right", va="top",
        fontsize=9, color="#B04A00")
a_res = 8*DX_LIP/2*1e6
a3.axvspan(50, a_res, color="0.55", alpha=0.15, zorder=0)
a3.axvline(a_res, color="0.35", lw=1.0, ls="--", zorder=3)
a3.text(a_res-4, 9.6, "$\\delta_\\omega$ < 8 cells", ha="right", va="top",
        fontsize=9, color="0.3", rotation=90)
for (a_um, s, q), col in zip(FAM, COL):
    a3.plot([a_um], [s], marker="o", color=col, ms=8, mec="k", mew=0.7, zorder=6)
a3.set_xlim(50, 420)
a3.set_ylim(0, 10)
a3.set_xlabel("$a$   [$\\mu$m]")
a3.set_ylabel("$s$")
a3.set_title(r"(c)  loci of constant $\delta^*_c$  [$\mu$m]", fontsize=10.5, pad=6)

for ax in (a1, a2, a3):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"iso_dstar_family.png", dpi=250)
fig.savefig(HERE/"iso_dstar_family.pdf")
print("\nwrote iso_dstar_family.png/pdf")
