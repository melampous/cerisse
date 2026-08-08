#!/usr/bin/env python3
"""OLD-vs-NEW RZ near-axis instability diagnosis (panda_mj1p42_L3 solver A/B).

NEW aborts at L0 step ~348 / L3 step 2703, t=7.97e-5 s, cell (0,134):
r=Delta_r/2=3.9e-5 m (axis-adjacent), z=10.5 mm, min(rho*e)=-87728 (negative).

Three panels:
  (1) radial near-axis u_r at the matched pre-crash step 200 -> OLD == NEW.
  (2) on-axis axial profile at ~crash time (OLD, step 400) -> a strong
      jet-front shock/contact sits on the axis at z~10.4 mm.
  (3) same-cell state at the crash location: OLD (bounded) vs NEW (blow-up).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
g = 1.4
def load(f):
    return np.loadtxt(HERE / f)   # cols: coord, rho, u_r, u_z, p, T

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in"})
C_OLD, C_NEW = "#1f4e9c", "#b3251e"

fig, ax = plt.subplots(1, 3, figsize=(11.2, 3.4))

# ---- Panel 1: radial u_r near axis at step 200 (matched, pre-crash) ----
rn = load("rad_NEW_s200.txt"); ro = load("rad_OLD_s200.txt")
m = rn[:, 0] <= 1.5e-3
ax[0].plot(ro[m, 0]*1e3, ro[m, 2], "-o", color=C_OLD, ms=3, lw=1.3, label="OLD")
ax[0].plot(rn[m, 0]*1e3, rn[m, 2], "--s", color=C_NEW, ms=3, lw=1.3, label="NEW", mfc="none")
ax[0].set_xlabel("r [mm]"); ax[0].set_ylabel(r"$u_r$ [m/s]")
ax[0].set_title("(1) step 200 (t=0.049 ms), z=10.5 mm\nOLD $\\equiv$ NEW — identical start", fontsize=8)
ax[0].legend(fontsize=8, frameon=False); ax[0].grid(True, ls=":", lw=0.5, color="0.85")

# ---- Panel 2: on-axis axial profile at ~crash time (OLD s400) ----
a = load("ax_OLD_s400.txt"); z = a[:, 0]*1e3
mm = (z >= 6) & (z <= 15)
ax2 = ax[1]; ax2b = ax2.twinx()
ax2.plot(z[mm], a[mm, 1], "-", color="#1f4e9c", lw=1.5, label=r"$\rho$")
ax2.plot(z[mm], a[mm, 3]/1000, "-", color="#6b6b6b", lw=1.2, label=r"$u_z$ /1000")
ax2b.plot(z[mm], a[mm, 5], "-", color="#c8781e", lw=1.2, label="T")
ax2.axvline(10.5, color="k", ls="--", lw=0.8)
ax2.annotate("crash z=10.5mm\n(shock/contact\non axis)", xy=(10.5, 1.0), xytext=(11.4, 0.7),
             fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.6))
ax2.set_xlabel("z [mm]"); ax2.set_ylabel(r"$\rho$ [kg/m³],  $u_z$/1000")
ax2b.set_ylabel("T [K]", color="#c8781e")
ax2.set_title("(2) on-axis at t=0.091 ms (OLD)\njet-front discontinuity crosses axis", fontsize=8)
ax2.grid(True, ls=":", lw=0.5, color="0.85")
h1, l1 = ax2.get_legend_handles_labels(); h2, l2 = ax2b.get_legend_handles_labels()
ax2.legend(h1+h2, l1+l2, fontsize=7, frameon=False, loc="upper right")

# ---- Panel 3: same-cell state at crash location (r=Dr/2, z=10.5mm) ----
# OLD from s400 (nearest crash time), first radial cell; NEW from abort diagnostic
old_ur, old_re = load("rad_OLD_s400.txt")[0, 2], load("rad_OLD_s400.txt")[0, 4]/(g-1)
new_ur, new_re = -720.0, -87728.0
ax3 = ax[2]
xs = [0, 1]
axb = ax3.twinx()
ax3.bar([-0.18, 0.82], [old_ur, new_ur], width=0.32, color=[C_OLD, C_NEW], label="u_r")
axb.bar([0.18, 1.18], [old_re, new_re], width=0.32, color=[C_OLD, C_NEW], alpha=0.45)
ax3.axhline(0, color="k", lw=0.6)
ax3.set_xticks([0, 1]); ax3.set_xticklabels(["OLD\n(s400, bounded)", "NEW\n(crash)"], fontsize=8)
ax3.set_ylabel(r"$u_r$ [m/s]  (solid)")
axb.set_ylabel(r"$\rho e$ [J/m³]  (faded)")
ax3.set_title("(3) same cell r=$\\Delta r$/2, z=10.5mm\nOLD $u_r$=-1.3, $\\rho e$>0  |  NEW $u_r$=-720, $\\rho e$<0", fontsize=8)
ax3.annotate(f"-720 m/s", xy=(0.82, -720), xytext=(0.4, -500), fontsize=7,
             arrowprops=dict(arrowstyle="->", lw=0.6))
axb.annotate(r"$\rho e$=-87728", xy=(1.18, -87728), xytext=(0.55, -60000), fontsize=7, color=C_NEW,
             arrowprops=dict(arrowstyle="->", lw=0.6, color=C_NEW))

fig.suptitle("RZ near-axis instability: NEW (metric-weighted radial split) vs OLD (ordinary split + geom. source) — $M_j$=1.42 axisymmetric, L3",
             fontsize=9.5, y=1.02)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"rz_crash_old_vs_new.{ext}", dpi=300, bbox_inches="tight")
print("wrote rz_crash_old_vs_new.png/.pdf")
