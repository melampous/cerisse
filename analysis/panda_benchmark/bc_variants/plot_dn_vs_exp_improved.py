#!/usr/bin/env python3
"""Plot closure-referenced shock-cell phase measures for the four exit profiles."""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTDIR = Path(__file__).resolve().parent
MARKS = np.load(OUTDIR / "shockcell_marks.npz")
N_FIT = 4

PROFILES = {
    "Baseline": {
        "xfc": float(MARKS["baseline_xfc"]),
        "xn": MARKS["baseline_xmax"],
        "color": "#aab0b6",
        "linestyle": "-",
    },
    "Top hat": {
        "xfc": float(MARKS["tophat_xfc"]),
        "xn": MARKS["tophat_xmax"],
        "color": "#0b6ef5",
        "linestyle": "--",
    },
    "Wall tanh": {
        "xfc": float(MARKS["walltanh_xfc"]),
        "xn": MARKS["walltanh_xmax"],
        "color": "#2ca02c",
        "linestyle": "-.",
    },
    "Thin shifted tanh": {
        "xfc": float(MARKS["shift100_xfc"]),
        "xn": MARKS["shift100_xmax"],
        "color": "#ff6a00",
        "linestyle": ":",
    },
}


plt.rcParams.update(
    {
        "font.size": 9.5,
        "axes.labelsize": 9.5,
        "axes.linewidth": 0.7,
        "legend.fontsize": 7.4,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

d_exp = (MARKS["experiment_xmax"] - float(MARKS["experiment_xfc"]))[:N_FIT]
x_line = np.linspace(0.0, 1.08 * d_exp.max(), 100)

fig, (ax_a, ax_b) = plt.subplots(
    1,
    2,
    figsize=(8.4, 2.95),
    gridspec_kw={"width_ratios": [1.12, 0.88], "wspace": 0.24},
)
for ax in (ax_a, ax_b):
    ax.set_axisbelow(True)
    ax.grid(True, color="0.89", linewidth=0.45, linestyle="--")

ax_a.plot(x_line, x_line, color="0.52", linewidth=1.0, linestyle=":", label=r"$y=x$")
ax_b.axhline(0.0, color="0.52", linewidth=1.0, linestyle=":")

print(f"{'Profile':20s} {'alpha':>8s} {'RMS':>8s} {'R2':>10s}")
for name, data in PROFILES.items():
    d_calc = (data["xn"] - data["xfc"])[:N_FIT]
    alpha = float(np.dot(d_calc, d_exp) / np.dot(d_exp, d_exp))
    fit_residual = d_calc - alpha * d_exp
    rms = float(np.sqrt(np.mean(fit_residual**2)))
    r_squared = float(
        1.0 - np.sum(fit_residual**2) / np.sum((d_calc - d_calc.mean()) ** 2)
    )
    print(f"{name:20s} {alpha:8.4f} {rms:8.4f} {r_squared:10.6f}")

    style = {
        "color": data["color"],
        "linestyle": data["linestyle"],
        "linewidth": 1.3,
    }
    ax_a.plot(x_line, alpha * x_line, **style)
    ax_a.plot(
        d_exp,
        d_calc,
        linestyle="none",
        marker="o",
        markersize=5.7,
        markerfacecolor=data["color"],
        markeredgecolor="white",
        markeredgewidth=0.55,
        label=rf"{name} ($\alpha_q={alpha:.3f}$)",
        zorder=5,
    )
    ax_b.plot(x_line, (alpha - 1.0) * x_line, **style)
    ax_b.plot(
        d_exp,
        d_calc - d_exp,
        linestyle="none",
        marker="o",
        markersize=5.7,
        markerfacecolor=data["color"],
        markeredgecolor="white",
        markeredgewidth=0.55,
        zorder=5,
    )

ax_a.set_xlim(0.0, 1.08 * d_exp.max())
ax_a.set_ylim(0.0, 1.15 * d_exp.max())
ax_a.set_xlabel(r"Experimental $d_n^{\mathrm{exp}}$")
ax_a.set_ylabel(r"Calculated $d_n^{(q)}$")
ax_a.legend(
    loc="upper left",
    frameon=False,
    handlelength=1.45,
    handletextpad=0.6,
    labelspacing=0.45,
)

ax_b.set_xlim(-0.05, 1.08 * d_exp.max())
ax_b.set_ylim(-0.30, 0.06)
ax_b.set_xlabel(r"Experimental $d_n^{\mathrm{exp}}$")
ax_b.set_ylabel(r"Cumulative difference $\Delta d_n^{(q)}$", labelpad=3.0)

for ax, panel in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.text(
        0.97,
        0.95,
        panel,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
    )

fig.subplots_adjust(left=0.085, right=0.985, bottom=0.20, top=0.97)
for suffix in ("png", "pdf"):
    fig.savefig(OUTDIR / f"dn_vs_exp_improved.{suffix}", dpi=300)
plt.close(fig)
