"""F1: Mach-disk formation criteria across NPR (Muraoka-style, Figs. 6-9).
F2: Mach-disk location and diameter vs NPR with literature overlays.
All 2D RZ cases: plateau-window online statistics at D_e/256.
m142 3D point: AX3T_L5 online statistics at D_e/256 (t-avg window as in paper).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lit_common import (D, GAM, RGAS, PINF, BASE, AUX, CASES_2D, load2d,
                        axis_state, pocket, disk_radius, mach_field,
                        load_m180_slice, m180_axis_inst)

plt.rcParams.update({"font.size": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True,
                     "ytick.right": True})

# ---------------- gather per-case axis metrics ----------------
rows = []          # npr, x_in, x_out, Mmin, d_disk, Ptmin_ratio, Tmax_ratio
axprofiles = {}
for kind, npr, fn in CASES_2D:
    c = load2d(fn)
    z, p, T, rho, uz, M = axis_state(c)
    axprofiles[npr] = (z, M)
    pk = pocket(z, M)
    P0 = npr * PINF
    Pt = np.where(M <= 1, p * (1 + 0.2 * M**2)**3.5,
                  p * (1 + 0.2 * M**2)**3.5 *
                  (((GAM + 1) * M**2 / ((GAM - 1) * M**2 + 2))**(GAM / (GAM - 1)) *
                   ((GAM + 1) / (2 * GAM * M**2 - (GAM - 1)))**(1 / (GAM - 1))))
    # Pt here = flow total pressure (isentropic below M=1, normal-shock loss above)
    # T0 from the exit state
    i0 = np.argmin(np.abs(z - 0.05))
    T0 = T[i0] * (1 + 0.2 * M[i0]**2)
    if pk:
        x_in, x_out, Mmin, cens = pk
        win = (z > x_in - 0.2) & (z < min(x_out + 0.5, 5.9))
        dd = 2.0 * disk_radius(c, x_in)
    else:
        x_in = x_out = Mmin = np.nan
        cens = False
        win = z < 4.0
        dd = 0.0
    rows.append((npr, x_in, x_out, Mmin, dd,
                 float((Pt[win] / P0).min()), float((T[win] / T0).max()),
                 float(cens)))
    print(f"NPR={npr:7.3f} pocket=({x_in:.3f},{x_out:.3f}) cens={cens} "
          f"Mmin={Mmin:.3f} d/D={dd:.3f} Ptmin/P0={(Pt[win]/P0).min():.3f} "
          f"Tmax/T0={(T[win]/T0).max():.3f}")
R = np.array([[r[i] for r in rows] for i in range(8)])
NPR, XIN, XOUT, MMIN, DDISK, PTMIN, TMAX, CENS = R
CENS = CENS > 0.5

# m142 3D point (online stats, D/256)
a = np.load(AUX + "/shockcell/AX3T_L5.npz")
x3 = a["x"] / D
M3 = a["ax_x_velocityMEAN"] / np.sqrt(GAM * RGAS * a["ax_temperatureMEAN"])
pk3 = pocket(x3, M3)
print("m142-3D pocket:", pk3)

# m180 3D instantaneous point
s = load_m180_slice()
xi, pi, Ti, ui, Mi = m180_axis_inst(s)
pk180 = pocket(xi, Mi)
print("m180-3D inst pocket:", pk180)

THR = (3.03, 3.12)   # Muraoka onset band (NPR)

# ================= F1: criteria figure =================
fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6))
(axA, axB), (axC, axD) = axes

# (a) axis Mach profiles around the threshold
sel_npr = [2.394, 2.800, 3.273446, 3.671, 5.746]
blues = [plt.cm.Blues(v) for v in np.linspace(0.45, 1.0, len(sel_npr))]
for c_, nn in zip(blues, sel_npr):
    z, M = axprofiles[nn]
    axA.plot(z, M, color=c_, lw=1.5, label=f"NPR={nn:.2f} (RZ)")
axA.plot(x3, M3, "--", color="#d95f02", lw=1.6,
         label="NPR=3.27, 3D ($D_e/256$)")
axA.axhline(1.0, color="0.3", lw=0.8, ls=":")
axA.set_xlim(0, 4); axA.set_ylim(0, 4.2)
axA.set_xlabel(r"$x/D_e$"); axA.set_ylabel(r"axis mean Mach number")
axA.legend(fontsize=8, loc="upper right")
axA.text(0.02, 0.93, "(a)", transform=axA.transAxes, fontsize=12)

# (b) subsonic pocket extent vs NPR
has = ~np.isnan(XIN)
axB.axvspan(*THR, color="#fee0b6", alpha=0.8,
            label="Muraoka & Hiejima onset band")
axB.plot(NPR[has], XIN[has], "o-", color="#08306b", ms=5,
         label=r"$x_{\rm in}$ (axis $M{=}1$ downcross)")
ok = has & ~CENS
axB.plot(NPR[ok], XOUT[ok], "s--", color="#2171b5", ms=5, mfc="none",
         label=r"$x_{\rm out}$ (re-acceleration)")
axB.plot(NPR[has & CENS], np.full((has & CENS).sum(), 6.0), "^",
         color="#2171b5", ms=7, mfc="none",
         label=r"$x_{\rm out}>6$ (RZ subsonic core persists)")
if pk3:
    axB.plot(3.273446, pk3[0], "^", color="#d95f02", ms=9, mec="k", mew=0.6,
             label="3D m142 (stats)")
if pk180:
    axB.plot(5.746, pk180[0], "v", color="#7b3294", ms=9, mec="k", mew=0.6,
             label="3D m180 (instantaneous)")
axB.set_xscale("log")
axB.set_xticks([2, 3, 5, 10, 20]); axB.set_xticklabels(["2", "3", "5", "10", "20"])
axB.set_xlabel("NPR"); axB.set_ylabel(r"$x/D_e$")
axB.set_ylim(0.5, 6.5)
axB.legend(fontsize=8, loc="center left")
axB.text(0.02, 0.93, "(b)", transform=axB.transAxes, fontsize=12)

# (c) minimum axis total-pressure ratio
axC.axvspan(*THR, color="#fee0b6", alpha=0.8)
axC.axhline(0.6, color="0.4", lw=0.8, ls=":")
axC.text(1.95, 0.605, "~40% loss at onset (Muraoka)", fontsize=8, color="0.3",
         va="bottom")
axC.plot(NPR, PTMIN, "o-", color="#08306b", ms=5)
axC.set_xscale("log")
axC.set_xticks([2, 3, 5, 10, 20]); axC.set_xticklabels(["2", "3", "5", "10", "20"])
axC.set_xlabel("NPR")
axC.set_ylabel(r"$\min\,(p_t/p_0)$ on axis, first-cell window")
axC.text(0.02, 0.93, "(c)", transform=axC.transAxes, fontsize=12)

# (d) peak axis static temperature
axD.axvspan(*THR, color="#fee0b6", alpha=0.8)
axD.plot(NPR, TMAX, "o-", color="#08306b", ms=5)
axD.set_xscale("log")
axD.set_xticks([2, 3, 5, 10, 20]); axD.set_xticklabels(["2", "3", "5", "10", "20"])
axD.set_xlabel("NPR")
axD.set_ylabel(r"$\max\,(T/T_0)$ on axis, first-cell window")
axD.text(0.02, 0.93, "(d)", transform=axD.transAxes, fontsize=12)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{BASE}/litfig_criteria.{ext}", dpi=250)
plt.close(fig)

# ================= F2: x_MD and d_MD vs NPR =================
# plateau-aggregated x1 from the sweep detector + nprfix metrics
d = np.load(AUX + "/npr_sweep_case/xs_x1_vs_npr_final.npz")
plats = [1.89293, 2.394, 2.800, 3.273446, 3.671, 4.200, 4.900, 5.746]
x1m, x1s = [], []
for pv in plats:
    m = (~d["ramp"]) & (np.abs(d["npr"] - pv) < 1e-6) & np.isfinite(d["x1"])
    tt = d["t"][m]
    if len(tt) > 4:                       # keep the settled last 60 %
        m &= d["t"] > tt[0] + 0.4 * (tt[-1] - tt[0])
    x1m.append(np.nanmean(d["x1"][m])); x1s.append(np.nanstd(d["x1"][m]))
nf = np.load(AUX + "/npr_sweep_case/nprfix_metrics.npz")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.6))
nn = np.linspace(1.9, 22, 300)
ax1.plot(nn, 0.67 * np.sqrt(nn), "-", color="0.45", lw=1.4,
         label=r"$x_{\rm MD}/D_e=0.67\sqrt{\rm NPR}$ (Crist; Muraoka)")
# NPR=1.893 (perfectly expanded): no shock-cell train, x1 detector invalid
ax1.errorbar(plats[1:], x1m[1:], yerr=x1s[1:], fmt="o", color="#2171b5",
             ms=5, capsize=2,
             label=r"$x_1$ first compression peak (RZ sweep)")
ax1.errorbar(nf["npr"], nf["x1"], yerr=nf["x1std"], fmt="o", color="#2171b5",
             ms=5, mfc="none", capsize=2)
ax1.plot(NPR[has], XIN[has], "d-", color="#08306b", ms=6, lw=1.2,
         label=r"$x_{\rm MD}$: axis $M{=}1$ downcross (RZ mean)")
if pk3:
    ax1.plot(3.273446, pk3[0], "^", color="#d95f02", ms=9, mec="k", mew=0.6,
             label="3D m142 (stats)")
if pk180:
    ax1.plot(5.746, pk180[0], "v", color="#7b3294", ms=9, mec="k", mew=0.6,
             label="3D m180 (inst.)")
ax1.plot([5.60], [1.43], "*", color="#238b45", ms=13, mec="k", mew=0.5,
         label=r"Li et al. LES: $H_{\rm MD}/D=1.43$ (NPR 5.60)")
ax1.plot([5.0], [1.5], "P", color="#e7298a", ms=9, mec="k", mew=0.5,
         label=r"Zapryagaev exp.: first cell $\approx 1.5D$ (NPR 5.0)")
ax1.set_xscale("log")
ax1.set_xticks([2, 3, 5, 10, 20]); ax1.set_xticklabels(["2", "3", "5", "10", "20"])
ax1.set_xlabel("NPR"); ax1.set_ylabel(r"$x/D_e$")
ax1.set_ylim(0, 3.6)
ax1.legend(fontsize=7.5, loc="upper left")
ax1.text(0.02, 0.03, "(a)", transform=ax1.transAxes, fontsize=12)

sel = DDISK > 0
ax2.plot(nn, 2.5 * np.log10(nn / 1.89293) - 0.75, "-", color="0.45", lw=1.4,
         label=r"$d_{\rm MD}/D_e=\frac{5}{2}\log_{10}(p_e/p_\infty)-\frac{3}{4}$ (Muraoka)")
ax2.plot(NPR[sel], DDISK[sel], "d-", color="#08306b", ms=6, lw=1.2,
         label=r"subsonic-pocket diameter (RZ mean, $M{<}1$)")
ax2.plot([5.60], [0.38], "*", color="#238b45", ms=13, mec="k", mew=0.5,
         label=r"Li et al. LES: $W_{\rm MD}/D=0.38$ (NPR 5.60)")
ax2.axhline(0, color="0.6", lw=0.6)
ax2.set_xscale("log")
ax2.set_xticks([2, 3, 5, 10, 20]); ax2.set_xticklabels(["2", "3", "5", "10", "20"])
ax2.set_xlabel("NPR"); ax2.set_ylabel(r"$d_{\rm MD}/D_e$")
ax2.legend(fontsize=8, loc="upper left")
ax2.text(0.02, 0.03, "(b)", transform=ax2.transAxes, fontsize=12)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{BASE}/litfig_xmd_dmd.{ext}", dpi=250)
plt.close(fig)
print("saved litfig_criteria + litfig_xmd_dmd")
