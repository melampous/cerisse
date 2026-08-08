#!/usr/bin/env python3
"""Plot a same-source one-transit comparison that includes Skew--JST."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


METHODS = (
    ("llf-wenoz5", "LLF WENO-Z5", "#0072B2", "o"),
    ("llf-teno5", "LLF TENO5", "#009E73", "s"),
    ("afd-hllc-wenoz5", "AFD HLLC WENO-Z5", "#D55E00", "^"),
    ("skew-jst-o4", "Skew JST", "#CC79A7", "D"),
)


def read_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    converted: list[dict[str, float | str]] = []
    for row in rows:
        converted.append(
            {
                key: value if key == "plotfile" else float(value)
                for key, value in row.items()
            }
        )
    return sorted(converted, key=lambda row: float(row["base_n"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output-stem", default="one_transit_four_method_comparison")
    args = parser.parse_args()

    all_rows: dict[str, list[dict[str, float | str]]] = {}
    for method, _, _, _ in METHODS:
        all_rows[method] = read_rows(args.root / method / "errors.csv")

    combined_csv = args.root / f"{args.output_stem}.csv"
    first = True
    with combined_csv.open("w", newline="") as stream:
        writer = None
        for method, label, _, _ in METHODS:
            for row in all_rows[method]:
                output = {"method": method, "method_label": label, **row}
                if first:
                    writer = csv.DictWriter(stream, fieldnames=list(output))
                    writer.writeheader()
                    first = False
                assert writer is not None
                writer.writerow(output)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 8.5,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.25,
            "lines.markersize": 4.5,
        }
    )
    figure, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)

    for method, label, colour, marker in METHODS:
        rows = all_rows[method]
        n = [float(row["base_n"]) for row in rows]
        density_l2 = [float(row["density_L2"]) for row in rows]
        ke = [float(row["perturbation_ke_ratio"]) for row in rows]
        enstrophy = [float(row["enstrophy_ratio"]) for row in rows]
        style = {
            "color": colour,
            "marker": marker,
            "markerfacecolor": colour,
            "markeredgecolor": "black",
            "markeredgewidth": 0.55,
            "label": label,
        }
        axes[0].loglog(n, density_l2, **style)
        axes[1].plot(n, ke, **style)
        axes[2].plot(n, enstrophy, **style)

    titles = (
        r"(a) Density error",
        r"(b) Perturbation kinetic energy",
        r"(c) Enstrophy",
    )
    ylabels = (
        r"$L_2(\rho-\rho_{\mathrm{exact}})$",
        r"$K^\prime/K^\prime_{\mathrm{exact}}$",
        r"$Z/Z_{\mathrm{exact}}$",
    )
    for axis, title, ylabel in zip(axes, titles, ylabels):
        axis.set_title(title)
        axis.set_xlabel(r"Cells per direction, $N$")
        axis.set_ylabel(ylabel)
        axis.set_xscale("log", base=2)
        axis.set_xticks([40, 80, 160])
        axis.xaxis.set_major_formatter(ScalarFormatter())
        axis.grid(True, which="major", color="0.85", linestyle=":", linewidth=0.65)
        axis.tick_params(direction="in", top=True, right=True)
    axes[0].set_yscale("log")
    axes[1].axhline(1.0, color="0.35", linestyle="--", linewidth=0.8, zorder=0)
    axes[2].axhline(1.0, color="0.35", linestyle="--", linewidth=0.8, zorder=0)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=4,
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.45,
    )

    png = args.root / f"{args.output_stem}.png"
    pdf = args.root / f"{args.output_stem}.pdf"
    figure.savefig(png, dpi=300, bbox_inches="tight")
    figure.savefig(
        pdf,
        bbox_inches="tight",
        metadata={
            "Title": "Same-source isentropic-vortex method comparison",
            "Author": "CERISSE verification workflow",
        },
    )
    print(f"wrote {combined_csv}")
    print(f"wrote {png}")
    print(f"wrote {pdf}")


if __name__ == "__main__":
    main()
