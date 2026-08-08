#!/usr/bin/env python3
"""Plot the cell-by-cell AFD-HLLC-WENO failure mechanism against HLLC-MUSCL."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "compare_hllc_muscl_vs_afd_weno"


def target_history() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    grouped: dict[str, list[tuple[float, float]]] = {}
    with (DATA / "point_target_4_171_timeline.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            grouped.setdefault(row["scheme"], []).append(
                (float(row["time_us"]), float(row["rhoe_J_m3"]))
            )
    return {
        name: (np.array([v[0] for v in values]), np.array([v[1] for v in values]))
        for name, values in grouped.items()
    }


def stage1_stencil() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    log = DATA / "afd_fallback_off" / "cell_trace_radius3_fail.log"
    pattern = re.compile(
        r"label=rk2_stage1_fe .* cell=\(4,(\d+)\).*"
        r"values=\(rho=([^,]+),mx=([^,]+),my=([^,]+),mz=([^,]+),E=([^\)]+)\)"
        r".*internal_energy_density=([^ ]+)"
    )
    rows: dict[int, tuple[float, float]] = {}
    for line in log.read_text().splitlines():
        match = pattern.search(line)
        if not match:
            continue
        j = int(match.group(1))
        rho = float(match.group(2))
        un = float(match.group(4)) / rho
        pressure = 0.4 * float(match.group(7))
        rows[j] = (un, pressure)
    indices = np.arange(168, 175)
    return (
        indices,
        np.array([rows[j][0] for j in indices]),
        np.array([rows[j][1] for j in indices]),
    )


def main() -> None:
    history = target_history()
    indices, cell_un, cell_pressure = stage1_stencil()

    afd_un = -8860.998302447642
    muscl_un = -734.4107825250052
    afd_pressure = 56.140161965060543
    muscl_pressure = 10041.868146426694
    afd_energy_flux = 107012302603.95641
    llf_energy_flux = 134334938.15265739
    muscl_energy_flux = 106248727.29886031
    dx = 1.06129666667e-3
    dt = 2.0e-7
    afd_face_speed = 8876.985223066904
    muscl_face_speed = 929.7963621413825

    colors = {"muscl": "#2155a6", "afd_off": "#c5283d", "afd_on": "#16835d"}
    labels = {
        "muscl": "HLLC–MUSCL O2",
        "afd_off": "AFD–HLLC–WENO (LLF off)",
        "afd_on": "AFD–HLLC–WENO (LLF on)",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.5), constrained_layout=True)

    ax = axes[0, 0]
    for name in ("muscl", "afd_off", "afd_on"):
        time_us, rhoe = history[name]
        mask = time_us >= 60.0
        ax.plot(time_us[mask], rhoe[mask] / 1000.0, lw=2.0,
                color=colors[name], label=labels[name])
    ax.axvline(73.8, color="black", ls="--", lw=1.2)
    ax.scatter([73.8], [-36341962.12 / 1000.0], marker="x", s=70,
               color=colors["afd_off"], clip_on=False, zorder=5)
    ax.annotate("RK2 stage 2 fails\n$\\rho e=-3.63\\times10^7$ J/m³",
                xy=(73.8, 0.0), xytext=(69.6, 45.0),
                arrowprops={"arrowstyle": "->", "color": "black"}, fontsize=9)
    ax.set_xlim(60.0, 74.0)
    ax.set_ylim(0.0, 100.0)
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("target-cell $\\rho e$ [kJ/m³]")
    ax.set_title("(a) Local precursor at cell (4,171)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(indices, cell_un, "o-", color="#555555", lw=1.8,
            label="AFD stage-1 cell centers")
    ax.axvline(170.5, color="black", ls=":", lw=1.2)
    ax.scatter([170.5], [afd_un], color=colors["afd_off"], marker="v", s=70,
               label="AFD reconstructed right")
    ax.scatter([170.5], [muscl_un], color=colors["muscl"], marker="s", s=55,
               label="MUSCL reconstructed right")
    ax.set_yscale("symlog", linthresh=1000.0)
    ax.set_xlabel("axial cell index $j$")
    ax.set_ylabel("normal velocity [m/s]")
    ax.set_title("(b) WENO creates a 12× velocity overshoot")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[1, 0]
    names = ["AFD–HLLC\nLLF off", "HLLC–MUSCL"]
    pressure = [afd_pressure, muscl_pressure]
    speed = [abs(afd_un), abs(muscl_un)]
    x = np.arange(2)
    width = 0.36
    b1 = ax.bar(x - width / 2, pressure, width, color=[colors["afd_off"], colors["muscl"]],
                alpha=0.75, label="pressure [Pa]")
    ax2 = ax.twinx()
    b2 = ax2.bar(x + width / 2, speed, width, color=[colors["afd_off"], colors["muscl"]],
                 hatch="//", alpha=0.55, label="$|u_n|$ [m/s]")
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax.set_xticks(x, names)
    ax.set_ylabel("reconstructed pressure [Pa]")
    ax2.set_ylabel("reconstructed $|u_n|$ [m/s]")
    ax.set_title("(c) Positive state check passes, but state is extreme")
    ax.grid(axis="y", alpha=0.25)
    ax.legend([b1, b2], ["pressure", "$|u_n|$"], fontsize=8, loc="upper center")

    ax = axes[1, 1]
    flux_names = ["AFD–HLLC\nLLF off", "AFD shock-LLF", "HLLC–MUSCL"]
    flux = [afd_energy_flux, llf_energy_flux, muscl_energy_flux]
    bars = ax.bar(flux_names, flux, color=[colors["afd_off"], colors["afd_on"], colors["muscl"]])
    ax.set_yscale("log")
    ax.set_ylabel("$|\\hat F_E|$ [W/m²]")
    ax.set_title("(d) Shared-face energy-flux amplification")
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, flux):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.18,
                f"{value:.2e}", ha="center", va="bottom", fontsize=8)
    afd_cfl = dt * afd_face_speed / dx
    muscl_cfl = dt * muscl_face_speed / dx
    ax.text(0.03, 0.95,
            f"face CFL: AFD={afd_cfl:.3f}, MUSCL={muscl_cfl:.3f}\n"
            f"AFD/MUSCL flux={afd_energy_flux / muscl_energy_flux:.0f}×",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox={"boxstyle": "round", "fc": "white", "alpha": 0.9})

    fig.suptitle(
        "Run165 all-fluid RZ: cell-by-cell collapse mechanism at 73.8 µs\n"
        "uniform L1, dx=dr=1.0613 mm, fine dt=0.2 µs, no IBM",
        fontsize=13,
    )
    output = DATA / "hllc_muscl_vs_afd_collapse_mechanism.png"
    fig.savefig(output, dpi=220)
    print(output)


if __name__ == "__main__":
    main()
