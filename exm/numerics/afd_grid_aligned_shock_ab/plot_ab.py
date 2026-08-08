#!/usr/bin/env python3
"""Plot the archived Mach-10 shock-fallback A/B histories."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_history(root: Path, name: str) -> list[dict[str, object]]:
    path = root / "runs" / name / "metrics_history.json"
    return json.loads(path.read_text())


def curve(history: list[dict[str, object]], key: str) -> np.ndarray:
    return np.asarray([float(row[key]) for row in history])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.result_root
    output = args.output or root / "shock_stability_ab"

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
        }
    )

    colors = {0: "#0072B2", 1: "#D55E00"}
    labels = {
        0: r"AFD--HLLC, fallback off",
        1: r"shock LLF fallback on",
    }
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)

    histories: dict[tuple[int, int], list[dict[str, object]]] = {}
    markers = {200: "o", 400: "s"}
    lines = {200: "--", 400: "-"}
    for n in (200, 400):
        history = {
            llf: load_history(root, f"n{n}_seeded_llf{llf}")
            for llf in (0, 1)
        }
        for llf in (0, 1):
            histories[(n, llf)] = history[llf]
            axes[0].semilogy(
                curve(history[llf], "time"),
                curve(history[llf], "shock_odd_even_abs_cells"),
                marker=markers[n],
                markersize=3.2,
                linewidth=1.2,
                linestyle=lines[n],
                color=colors[llf],
                label=rf"$N_x={n}$, {labels[llf]}",
            )
    axes[0].set_xlabel(r"$t$")
    axes[0].set_ylabel(r"odd--even amplitude, $|a_{\rm oe}|/\Delta x$")
    axes[0].grid(True, which="both", linewidth=0.35, alpha=0.35)
    axes[0].legend(loc="best", fontsize=6.8)

    # The finer grid gives the clearest final-front comparison.  Plotting the
    # locus itself avoids reducing the carbuncle mode to a single scalar.
    n = 400
    for llf in (0, 1):
        final = histories[(n, llf)][-1]
        locus = np.asarray(final["shock_locus"], dtype=float)
        dx = float(final["dx"])
        offset = (locus - np.mean(locus)) / dx
        y = (np.arange(locus.size) + 0.5) / locus.size
        axes[1].plot(
            offset,
            y,
            marker=markers[n],
            markersize=2.4,
            linewidth=0.9,
            color=colors[llf],
            label=labels[llf],
        )
    axes[1].axvline(0.0, color="0.55", linewidth=0.7)
    axes[1].set_xlabel(
        r"final front displacement, $(x_s-\overline{x}_s)/\Delta x$"
    )
    axes[1].set_ylabel(r"$y/L_y$")
    axes[1].set_ylim(0.0, 1.0)
    axes[1].grid(True, linewidth=0.35, alpha=0.35)

    for label, ax, xloc, align in zip(
        ("a", "b"), axes, (0.98, 0.02), ("right", "left")
    ):
        ax.text(
            xloc,
            0.96,
            rf"({label})",
            transform=ax.transAxes,
            ha=align,
            va="top",
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
