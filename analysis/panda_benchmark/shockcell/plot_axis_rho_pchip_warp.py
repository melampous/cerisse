#!/usr/bin/env python3
"""Diagnostic comparison of the original and axially aligned 3D L5 density.

Five-node empirical axial coordinate alignment: a monotone piecewise-cubic
(PCHIP) map through the origin and the first five principal compression
locations, so every peak is pinned exactly,

    xi_aligned = P(xi),  P(0)=0,  P(x_i^sim)=x_i^exp,  i=1..5,

extended linearly beyond the last node with the end slope.  Peak locations
on both sides are refined by three-point parabolic interpolation.  This is
an empirical axial coordinate alignment for phase-error diagnosis only; the
density values are not rescaled.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks


HERE = Path(__file__).resolve().parent
D_E = 0.0254
RHO_J = 1.6413
N_PEAKS = 5


def moving_average(values, coordinates, width):
    """Apply the same 0.05 D_e box filter used in the comparison figure."""
    count = max(3, int(round(width / (coordinates[1] - coordinates[0]))))
    if count % 2 == 0:
        count += 1
    result = np.convolve(values, np.ones(count) / count, mode="same")
    result[:count // 2] = values[:count // 2]
    result[-(count // 2):] = values[-(count // 2):]
    return result


def refined_peaks(x, y, n, prominence=0.02, xmin=0.5):
    """First n peak locations, refined by 3-point parabolic interpolation."""
    idx, _ = find_peaks(y, prominence=prominence)
    idx = idx[x[idx] > xmin][:n]
    out = []
    for i in idx:
        denominator = y[i - 1] - 2.0 * y[i] + y[i + 1]
        shift = 0.0
        if abs(denominator) > 0:
            shift = 0.5 * (y[i - 1] - y[i + 1]) / denominator
        h = 0.5 * (x[i + 1] - x[i - 1])
        out.append(x[i] + np.clip(shift, -1, 1) * h)
    return np.array(out)


def interval_rmse(x, numerical, experimental, lower, upper):
    selected = (x >= lower) & (x <= upper)
    return np.sqrt(np.mean((numerical[selected] - experimental[selected])**2))


simulation = np.load(HERE / "AX3T_L5.npz", allow_pickle=True)
experiment = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")

x_sim = simulation["x"] / D_E
rho_sim = simulation["ax_DensityMEAN"] / RHO_J
rho_sim = moving_average(rho_sim, x_sim, 0.05)
x_exp = experiment["x"]
rho_exp = experiment["rho"]

sim_peaks = refined_peaks(x_sim, rho_sim, N_PEAKS)
exp_peaks = refined_peaks(x_exp, rho_exp, N_PEAKS)

# Monotone map through the origin and the five peak pairs; linear beyond.
nodes_sim = np.concatenate(([0.0], sim_peaks))
nodes_exp = np.concatenate(([0.0], exp_peaks))
warp = PchipInterpolator(nodes_sim, nodes_exp)
end_slope = float(warp.derivative()(nodes_sim[-1]))
inside = x_sim <= nodes_sim[-1]
x_aligned = np.where(
    inside,
    warp(np.minimum(x_sim, nodes_sim[-1])),
    nodes_exp[-1] + end_slope * (x_sim - nodes_sim[-1]),
)

rho_original_at_exp = np.interp(x_exp, x_sim, rho_sim)
rho_aligned_at_exp = np.interp(x_exp, x_aligned, rho_sim)
residual_original = rho_original_at_exp - rho_exp
residual_aligned = rho_aligned_at_exp - rho_exp

rmse_original = np.sqrt(np.mean(residual_original**2))
rmse_aligned = np.sqrt(np.mean(residual_aligned**2))

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

fig, (ax_profile, ax_residual) = plt.subplots(
    2, 1, figsize=(7.2, 4.8), sharex=True,
    gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.08},
)

for axis in (ax_profile, ax_residual):
    axis.set_axisbelow(True)
    axis.grid(True, color="0.87", linestyle="--", linewidth=0.5)
    axis.tick_params(length=3)

domain_original = (x_sim >= x_exp.min()) & (x_sim <= x_exp.max())
domain_aligned = (x_aligned >= x_exp.min()) & (x_aligned <= x_exp.max())

for xp in exp_peaks:
    ax_profile.axvline(xp, color="0.82", linewidth=0.6, zorder=0)

ax_profile.plot(
    x_sim[domain_original], rho_sim[domain_original],
    color="#b3251e", linewidth=1.35, label="3D L5, original",
)
ax_profile.plot(
    x_aligned[domain_aligned], rho_sim[domain_aligned],
    color="#167c72", linewidth=1.35, linestyle="--",
    label="3D L5, aligned (5-node monotone map)",
)
ax_profile.plot(
    x_exp, rho_exp, "o", color="black", markersize=3.6,
    label="Experiment ($M_j=1.42$)", zorder=5,
)
ax_profile.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
ax_profile.set_ylim(0.3, 1.5)
ax_profile.legend(frameon=False, ncol=1, loc="upper right", fontsize=8)

ax_residual.axhline(0.0, color="0.35", linewidth=0.7)
ax_residual.plot(
    x_exp, residual_original, color="#b3251e", linewidth=1.0,
    marker="o", markersize=2.6,
    label=rf"original, r.m.s.={rmse_original:.3f}",
)
ax_residual.plot(
    x_exp, residual_aligned, color="#167c72", linewidth=1.0,
    linestyle="--", marker="o", markersize=2.6,
    label=rf"aligned, r.m.s.={rmse_aligned:.3f}",
)
ax_residual.set_ylabel(r"$(\langle\rho\rangle_{\rm num}-\langle\rho\rangle_{\rm exp})/\rho_j$")
ax_residual.set_xlabel(r"$x/D_e$")
ax_residual.text(
    0.018, 0.93,
    rf"r.m.s.: ${rmse_original:.3f}\;\rightarrow\;{rmse_aligned:.3f}$",
    transform=ax_residual.transAxes, ha="left", va="top", fontsize=8.2,
    bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=1.2),
)

ax_profile.set_xlim(x_exp.min(), x_exp.max())
fig.align_ylabels((ax_profile, ax_residual))
fig.subplots_adjust(left=0.105, right=0.985, top=0.975, bottom=0.105)

output = HERE / "diagnostic_axis_rho_pchip_warp"
fig.savefig(output.with_suffix(".png"), dpi=300)
fig.savefig(output.with_suffix(".pdf"))

mapped_peaks = warp(sim_peaks)
print(f"map slope at origin: {float(warp.derivative()(0.0)):.4f}, "
      f"end slope: {end_slope:.4f}")
print("peak   simulation  experiment  aligned  aligned-exp")
for index, (sim, exp, aligned) in enumerate(
    zip(sim_peaks, exp_peaks, mapped_peaks), start=1
):
    print(f"x{index}     {sim:8.6f}  {exp:8.6f}  {aligned:8.6f}  {aligned-exp:+9.6f}")
print(f"full-range r.m.s.: {rmse_original:.6f} -> {rmse_aligned:.6f}")
for lower, upper in ((0.6, 2.0), (2.0, 4.0), (4.0, 6.92)):
    before = interval_rmse(x_exp, rho_original_at_exp, rho_exp, lower, upper)
    after = interval_rmse(x_exp, rho_aligned_at_exp, rho_exp, lower, upper)
    print(f"{lower:.2f} <= x/D_e <= {upper:.2f}: {before:.6f} -> {after:.6f}")
mono = np.all(np.diff(x_aligned) > 0)
print(f"map monotone: {mono}")
print(f"wrote {output.with_suffix('.png').name} and {output.with_suffix('.pdf').name}")
