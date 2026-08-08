#!/usr/bin/env python3
"""Plot first-cell metrics against nozzle-exit displacement thickness."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parent
DE_UM = 25_400.0

profiles = (
    ("Top hat", 0.0, 0.945, 1.244, "#2474b7", "o"),
    ("Wall tanh", 136.547347, 0.930, 1.415, "#2f9e44", "^"),
    ("Thin shifted tanh", 300.101187, 0.922, 1.508, "#e67e22", "D"),
    ("Baseline", 762.306583, 0.883, 1.300, "#73777b", "s"),
)

delta = np.array([row[1] for row in profiles]) / DE_UM
x_fc = np.array([row[2] for row in profiles])
s_1 = np.array([row[3] for row in profiles])


def linear_fit(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    fitted = intercept + slope * x
    residual = y - fitted
    r_squared = 1.0 - np.sum(residual**2) / np.sum((y - np.mean(y)) ** 2)
    return slope, intercept, r_squared


fit_fc = linear_fit(delta, x_fc)
fit_s1 = linear_fit(delta, s_1)

plt.rcParams.update(
    {
        "font.size": 9.5,
        "axes.linewidth": 0.75,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "mathtext.fontset": "stix",
    }
)

fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.05))
x_fit = np.linspace(0.0, 0.0315, 160)

panels = (
    (
        axes[0],
        x_fc,
        fit_fc,
        r"$x_{fc}/D_e$",
        0.888,
        r"$x_{fc}/D_e=%.4f%+.4f\,\delta^*/D_e$" % (fit_fc[1], fit_fc[0]),
    ),
    (
        axes[1],
        s_1,
        fit_s1,
        r"$S_1/D_e$",
        1.359,
        r"$S_1/D_e=%.4f%+.4f\,\delta^*/D_e$" % (fit_s1[1], fit_s1[0]),
    ),
)

for panel_index, (ax, values, fit, ylabel, experimental, equation) in enumerate(panels):
    slope, intercept, r_squared = fit
    ax.set_axisbelow(True)
    ax.grid(True, color="0.89", linestyle="--", linewidth=0.5)
    ax.plot(
        x_fit,
        intercept + slope * x_fit,
        color="0.18",
        linewidth=1.2,
        label="Linear fit",
        zorder=2,
    )
    ax.axhline(
        experimental,
        color="0.35",
        linestyle=":",
        linewidth=1.0,
        label="Experiment",
        zorder=1,
    )
    for label, _, _, _, color, marker in profiles:
        idx = [row[0] for row in profiles].index(label)
        ax.plot(
            delta[idx],
            values[idx],
            marker=marker,
            markersize=6.8,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.6,
            linestyle="none",
            color=color,
            label=label,
            zorder=4,
        )
    ax.set_xlim(-0.0012, 0.032)
    ax.set_xlabel(r"Displacement thickness, $\delta^*/D_e$")
    ax.set_ylabel(ylabel)
    text_position = (0.43, 0.81) if panel_index == 0 else (0.04, 0.06)
    ax.text(
        *text_position,
        equation + "\n" + rf"$R^2={r_squared:.3f}$",
        transform=ax.transAxes,
        fontsize=8.2,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
    )
    ax.text(
        0.02,
        0.96,
        f"({chr(ord('a') + panel_index)})",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9.5,
    )

handles, labels = axes[0].get_legend_handles_labels()
order = [2, 3, 4, 5, 0, 1]
fig.legend(
    [handles[i] for i in order],
    [labels[i] for i in order],
    loc="upper center",
    bbox_to_anchor=(0.5, 1.015),
    ncol=3,
    handlelength=1.6,
    columnspacing=1.25,
    fontsize=8.2,
)
fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.87), w_pad=1.5)

for suffix in ("pdf", "png"):
    fig.savefig(
        OUT / f"profile_metrics_vs_displacement_thickness.{suffix}",
        dpi=300,
        bbox_inches="tight",
    )

print(
    "x_fc: intercept=%.8f slope=%.8f R2=%.8f"
    % (fit_fc[1], fit_fc[0], fit_fc[2])
)
print(
    "S1:   intercept=%.8f slope=%.8f R2=%.8f"
    % (fit_s1[1], fit_s1[0], fit_s1[2])
)
