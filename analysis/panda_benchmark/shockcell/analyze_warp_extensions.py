#!/usr/bin/env python3
"""Extensions of the 5-node PCHIP axial-alignment diagnosis.

(1) Grid convergence of the phase error: the same 5-peak monotone map is
    built independently for 3D L4 and 3D L5; if the two displacement
    fields coincide, the phase error is grid-converged (physical), not a
    resolution artefact.
(2) Post-alignment amplitude decomposition: per-cell peak-trough
    amplitudes and envelope decay rates of the aligned L5 solution vs the
    experiment (the residual the alignment cannot remove).
(3) Closure with the radial-profile registration: the axis displacement
    field delta(x) = x_exp - x_sim evaluated at the radial-survey columns
    vs the independently fitted radial offsets -Dx*_142.
(4) The same alignment applied to the 2D L6 solution (as many peak pairs
    as are detectable) -- quantifies that the 2D error is amplitude
    death, which no axial alignment can repair.

Same conventions as plot_axis_rho_pchip_warp.py: 0.05 D_e display
smoothing, peak prominence 0.02, x > 0.5, three-point parabolic
refinement, PCHIP through origin + peak pairs, linear extension with the
end slope.  Experiment = m142den phase-mean centreline (already rho/rho_j).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks, hilbert

HERE = Path(__file__).resolve().parent
D_E = 0.0254
RHO_J = 1.6413
N_PEAKS = 5


def moving_average(values, coordinates, width):
    count = max(3, int(round(width / (coordinates[1] - coordinates[0]))))
    if count % 2 == 0:
        count += 1
    result = np.convolve(values, np.ones(count) / count, mode="same")
    result[:count // 2] = values[:count // 2]
    result[-(count // 2):] = values[-(count // 2):]
    return result


def refined_extrema(x, y, n, prominence=0.02, xmin=0.5, troughs=False):
    yy = -y if troughs else y
    idx, _ = find_peaks(yy, prominence=prominence)
    idx = idx[x[idx] > xmin][:n]
    out = []
    for i in idx:
        d = y[i - 1] - 2.0 * y[i] + y[i + 1]
        shift = 0.5 * (y[i - 1] - y[i + 1]) / d if abs(d) > 0 else 0.0
        h = 0.5 * (x[i + 1] - x[i - 1])
        out.append(x[i] + np.clip(shift, -1, 1) * h)
    return np.array(out)


def build_warp(x_sim, rho_sim, exp_peaks, n):
    sim_peaks = refined_extrema(x_sim, rho_sim, n)
    m = min(len(sim_peaks), len(exp_peaks), n)
    nodes_sim = np.concatenate(([0.0], sim_peaks[:m]))
    nodes_exp = np.concatenate(([0.0], exp_peaks[:m]))
    warp = PchipInterpolator(nodes_sim, nodes_exp)
    end_slope = float(warp.derivative()(nodes_sim[-1]))
    inside = x_sim <= nodes_sim[-1]
    x_aligned = np.where(
        inside,
        warp(np.minimum(x_sim, nodes_sim[-1])),
        nodes_exp[-1] + end_slope * (x_sim - nodes_sim[-1]),
    )
    return sim_peaks[:m], x_aligned, end_slope


experiment = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
x_exp = experiment["x"]
rho_exp = experiment["rho"]
exp_peaks = refined_extrema(x_exp, rho_exp, N_PEAKS)

profiles = {}
for tag, fn, xk, rk in (("3D L5", "AX3T_L5.npz", "x", "ax_DensityMEAN"),
                        ("3D L4", "AX3T_L4.npz", "x", "ax_DensityMEAN"),
                        ("2D L6", "AXT_L6.npz", "z", "stat_DensityMEAN")):
    d = np.load(HERE / fn, allow_pickle=True)
    x = d[xk] / D_E
    r = moving_average(d[rk] / RHO_J, x, 0.05)
    profiles[tag] = (x, r)

# ---------------- (1) grid convergence of the displacement field ------
print("=" * 72)
print("(1) Phase-error grid convergence: 3D L4 vs 3D L5, independent warps")
fields = {}
for tag in ("3D L5", "3D L4"):
    x, r = profiles[tag]
    pk, x_al, slope = build_warp(x, r, exp_peaks, N_PEAKS)
    fields[tag] = (x, r, x_al, pk, slope)
    disp = pk - exp_peaks[:len(pk)]
    print("%s: peaks %s" % (tag, "  ".join("%.4f" % v for v in pk)))
    print("      displacement exp-sim %s   end slope %.4f"
          % ("  ".join("%+.4f" % (-v) for v in disp), slope))
x5, r5, xa5, pk5, s5 = fields["3D L5"]
x4, r4, xa4, pk4, s4 = fields["3D L4"]
# displacement fields as functions of the EXP coordinate
grid = np.linspace(0.6, 6.92, 400)
d5 = np.interp(grid, xa5, xa5 - x5)
d4 = np.interp(grid, xa4, xa4 - x4)
print("per-peak L4-L5 pin difference (sim side): %s"
      % "  ".join("%+.4f" % v for v in (pk4 - pk5)))
print("displacement-field difference |d_L4 - d_L5| over 0.6-6.92: "
      "mean %.4f  max %.4f  (field magnitude mean %.4f)"
      % (np.mean(np.abs(d4 - d5)), np.max(np.abs(d4 - d5)),
         np.mean(np.abs(d5))))

# ---------------- (2) post-alignment amplitude decomposition ----------
print("=" * 72)
print("(2) Amplitude decomposition after alignment (3D L5)")
rho_al_exp = np.interp(x_exp, xa5, r5)          # aligned sim on exp grid
exp_tr = refined_extrema(x_exp, rho_exp, N_PEAKS + 1, troughs=True)
sim_tr_al = refined_extrema(x_exp, rho_al_exp, N_PEAKS + 1, troughs=True)
exp_pk_v = np.interp(exp_peaks, x_exp, rho_exp)
sim_pk_v = np.interp(exp_peaks, x_exp, rho_al_exp)   # same pinned locations


def cell_amplitudes(pk_x, pk_v, tr_x, x, y):
    amps = []
    for n, xp in enumerate(pk_x):
        prev = tr_x[(tr_x < xp)]
        if len(prev) == 0:
            amps.append(np.nan); continue
        amps.append(pk_v[n] - np.interp(prev[-1], x, y))
    return np.array(amps)


A_exp = cell_amplitudes(exp_peaks, exp_pk_v, exp_tr, x_exp, rho_exp)
A_sim = cell_amplitudes(exp_peaks, sim_pk_v, sim_tr_al, x_exp, rho_al_exp)
print("cell   x_peak   A_exp    A_sim(aligned)   ratio sim/exp")
for n in range(len(exp_peaks)):
    print("%2d    %6.3f   %6.3f   %6.3f          %6.3f"
          % (n + 1, exp_peaks[n], A_exp[n], A_sim[n], A_sim[n] / A_exp[n]))

# envelope decay: detrend with a one-cell (1.2 D_e) moving average
FITLO, FITHI = 1.0, 6.0
tr_exp = moving_average(rho_exp, x_exp, 1.2)
tr_sim = moving_average(rho_al_exp, x_exp, 1.2)
env_exp = np.abs(hilbert(rho_exp - tr_exp))
env_sim = np.abs(hilbert(rho_al_exp - tr_sim))
m = (x_exp >= FITLO) & (x_exp <= FITHI)
c_exp = np.polyfit(x_exp[m], np.log(env_exp[m]), 1)
c_sim = np.polyfit(x_exp[m], np.log(env_sim[m]), 1)
print("envelope decay over %.1f-%.1f D_e:  exp L_d = %.2f D_e, "
      "aligned sim L_d = %.2f D_e" % (FITLO, FITHI, -1 / c_exp[0],
                                      -1 / c_sim[0]))

# ---------------- (3) closure with the radial registration ------------
print("=" * 72)
print("(3) Axis displacement field vs radial-profile offsets (-Dx*_142)")
RADIAL = [(0.76, +0.020), (0.92, +0.035), (1.08, +0.055),
          (1.56, -0.004), (1.96, +0.052), (3.08, +0.047)]
print("station col   -Dx*_142   axis delta(x)   difference")
for xc, mdx in RADIAL:
    dax = float(np.interp(xc, grid, d5))
    print("  %5.2f       %+6.3f     %+6.3f        %+6.3f"
          % (xc, mdx, dax, mdx - dax))

# ---------------- (4) the same alignment applied to 2D L6 -------------
print("=" * 72)
print("(4) 2D L6: alignment with all detectable peak pairs")
x2, r2 = profiles["2D L6"]
pk2_all = refined_extrema(x2, r2, N_PEAKS)
print("detectable 2D peaks (prom 0.02): %s"
      % "  ".join("%.3f" % v for v in pk2_all))
pk2, xa2, s2 = build_warp(x2, r2, exp_peaks, len(pk2_all))
r2_orig = np.interp(x_exp, x2, r2)
r2_al = np.interp(x_exp, xa2, r2)
rm_o = np.sqrt(np.mean((r2_orig - rho_exp)**2))
rm_a = np.sqrt(np.mean((r2_al - rho_exp)**2))
r5_orig = np.interp(x_exp, x5, r5)
rm5_o = np.sqrt(np.mean((r5_orig - rho_exp)**2))
rm5_a = np.sqrt(np.mean((rho_al_exp - rho_exp)**2))
print("2D L6 r.m.s.: %.3f -> %.3f  (%.0f%% reduction, %d pairs)"
      % (rm_o, rm_a, 100 * (1 - rm_a / rm_o), len(pk2)))
print("3D L5 r.m.s.: %.3f -> %.3f  (%.0f%% reduction, 5 pairs)"
      % (rm5_o, rm5_a, 100 * (1 - rm5_a / rm5_o)))

# ---------------- figure ----------------------------------------------
plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.7,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
fig, (axd, axe) = plt.subplots(1, 2, figsize=(7.2, 3.0))
for a in (axd, axe):
    a.set_axisbelow(True)
    a.grid(True, color="0.87", ls="--", lw=0.5)
    a.tick_params(length=3)

axd.plot(grid, d5, color="#b3251e", lw=1.3,
         label="axis alignment, 3D L5")
axd.plot(grid, d4, color="#e39c94", lw=1.0,
         label="axis alignment, 3D L4")
axd.plot(exp_peaks, exp_peaks - pk5, "o", color="#b3251e", ms=4, mew=0,
         label="peak pairs (L5)")
axd.plot([p[0] for p in RADIAL], [p[1] for p in RADIAL], "s", color="k",
         ms=4.5, mfc="none", mew=0.9,
         label="radial registration $-\\Delta x^\\ast$")
axd.axhline(0, color="0.35", lw=0.7)
axd.set_xlabel("$x/D_e$")
axd.set_ylabel("axial displacement $\\delta/D_e$ (exp $-$ sim)")
axd.set_xlim(0.5, 7)
axd.legend(frameon=False, fontsize=7, loc="upper left", handlelength=1.5)
axd.text(0.02, 0.03, "(a)", transform=axd.transAxes, fontsize=9)

axe.semilogy(x_exp, env_exp, "o", color="k", ms=2.8, mew=0,
             label="experiment")
axe.semilogy(x_exp, env_sim, color="#167c72", lw=1.2,
             label="3D L5, aligned")
for c, col in ((c_exp, "k"), (c_sim, "#167c72")):
    xx = np.linspace(FITLO, FITHI, 50)
    axe.semilogy(xx, np.exp(np.polyval(c, xx)), color=col, lw=0.8, ls=":")
axe.set_xlabel("$x/D_e$")
axe.set_ylabel("oscillation envelope $/\\rho_j$")
axe.set_xlim(0.5, 7)
axe.set_ylim(2e-2, 6e-1)
axe.legend(frameon=False, fontsize=7, loc="upper right", handlelength=1.5)
axe.text(0.02, 0.03, "(b)", transform=axe.transAxes, fontsize=9)
axe.text(0.35, 0.06, "$L_d^{\\rm exp}=%.1f\\,D_e$\n$L_d^{\\rm sim}=%.1f\\,D_e$"
         % (-1 / c_exp[0], -1 / c_sim[0]),
         transform=axe.transAxes, fontsize=7.5, va="bottom")

fig.tight_layout(w_pad=1.2)
out = HERE / "diagnostic_warp_extensions"
fig.savefig(out.with_suffix(".png"), dpi=300)
fig.savefig(out.with_suffix(".pdf"))
print("wrote", out.with_suffix(".png").name)
