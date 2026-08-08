#!/usr/bin/env python3
"""Track the flow peak and the dangerous face-171 WENO candidate in time."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_hllc_muscl_afd_weno_pointwise as pointwise
import analyze_hllc_same_state_reconstruction as reconstruction


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_hllc_muscl_vs_afd_weno"
CSV_PATH = RESULTS / "face171_peak_evolution.csv"
FIGURE_PATH = RESULTS / "face171_peak_evolution.png"
GAMMA = 1.4
RADIAL_I = 4
FACE_J = 171
PROFILE_JS = np.arange(168, 175)
PLOT_STEPS = {176, 178, 180, 182, 183, 184}


def make_cells(fields: dict[str, np.ndarray]) -> dict[int, dict[str, float]]:
    cells: dict[int, dict[str, float]] = {}
    for j in PROFILE_JS:
        rho = float(fields["rho"][RADIAL_I, j])
        pressure = float(fields["pressure"][RADIAL_I, j])
        cells[int(j)] = {
            "rho": rho,
            "p": pressure,
            "un": float(fields["uz"][RADIAL_I, j]),
            "ur": float(fields["ur"][RADIAL_I, j]),
            "ut2": 0.0,
            "c": float(np.sqrt(GAMMA * pressure / rho)),
        }
    return cells


def diagnose(
    label: str,
    time_us: float,
    cells: dict[int, dict[str, float]],
    x_centres: np.ndarray,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    pressure = np.asarray([cells[int(j)]["p"] for j in PROFILE_JS])
    mach = np.asarray(
        [
            np.hypot(cells[int(j)]["ur"], cells[int(j)]["un"])
            / cells[int(j)]["c"]
            for j in PROFILE_JS
        ]
    )
    p_index = int(np.argmin(pressure))
    mach_index = int(np.argmax(mach))
    p_j = int(PROFILE_JS[p_index])
    mach_j = int(PROFILE_JS[mach_index])

    _, right, _, right_details = reconstruction.weno_states(cells)
    w0_detail = right_details[0]
    wplus_detail = right_details[1]
    wminus_detail = right_details[2]
    candidate0_rho = float(w0_detail["candidate"][0]) / (1.0 - 1.0 / GAMMA)
    candidate0_wplus = float(wplus_detail["candidate"][0])
    candidate0_wminus = float(wminus_detail["candidate"][0])
    candidate0_pressure = candidate0_wplus + candidate0_wminus
    candidate0_velocity = np.nan
    if candidate0_rho > 0.0 and candidate0_pressure > 0.0:
        candidate0_velocity = (
            candidate0_wplus - candidate0_wminus
        ) / np.sqrt(GAMMA * candidate0_rho * candidate0_pressure)

    target = cells[FACE_J]
    row: dict[str, object] = {
        "label": label,
        "time_us": f"{time_us:.12g}",
        "local_pressure_min_j": p_j,
        "local_pressure_min_x_mm": f"{x_centres[p_j]:.12g}",
        "local_pressure_min_Pa": f"{pressure[p_index]:.17g}",
        "local_mach_max_j": mach_j,
        "local_mach_max_x_mm": f"{x_centres[mach_j]:.12g}",
        "local_mach_max": f"{mach[mach_index]:.17g}",
        "j171_rho": f"{target['rho']:.17g}",
        "j171_pressure_Pa": f"{target['p']:.17g}",
        "j171_un_m_s": f"{target['un']:.17g}",
        "j171_mach": f"{mach[np.where(PROFILE_JS == FACE_J)[0][0]]:.17g}",
        "right_wminus_candidate0_weight": f"{float(wminus_detail['weight'][0]):.17g}",
        "right_candidate0_pressure_Pa": f"{candidate0_pressure:.17g}",
        "right_candidate0_un_m_s": f"{candidate0_velocity:.17g}",
        "right_reconstructed_pressure_Pa": f"{right['p']:.17g}",
        "right_reconstructed_un_m_s": f"{right['un']:.17g}",
    }
    return row, pressure, mach


def main() -> None:
    pointwise.patch_yt_cylindrical_readonly_edge()
    pointwise.yt.funcs.mylog.setLevel(50)

    rows: list[dict[str, object]] = []
    profiles: dict[str, tuple[float, np.ndarray, np.ndarray]] = {}
    x_centres: np.ndarray | None = None

    for step in range(170, 185):
        fields, time_us, _, _, _, x_centres = pointwise.load_fields(
            pointwise.plotfile_for("afd_off", step)
        )
        cells = make_cells(fields)
        row, pressure, mach = diagnose(str(step), time_us, cells, x_centres)
        rows.append(row)
        if step in PLOT_STEPS:
            profiles[str(step)] = (time_us, pressure, mach)

    if x_centres is None:
        raise RuntimeError("No plotfiles were loaded")

    stage_cells = reconstruction.load_stencil()
    stage_row, stage_pressure, stage_mach = diagnose(
        "184_stage1_fe", 73.8, stage_cells, x_centres
    )
    rows.append(stage_row)
    profiles["184_stage1_fe"] = (73.8, stage_pressure, stage_mach)

    with CSV_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    times = np.asarray([float(row["time_us"]) for row in rows])
    x_profile = x_centres[PROFILE_JS]
    colours = plt.cm.viridis(np.linspace(0.05, 0.95, len(profiles)))
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), constrained_layout=True)

    for colour, (label, (time_us, pressure, mach)) in zip(colours, profiles.items()):
        style = "--" if label == "184_stage1_fe" else "-"
        axes[0, 0].plot(
            x_profile, pressure / 1.0e3, marker="o", ms=3.0,
            color=colour, linestyle=style, label=f"{time_us:.1f} us"
        )
        axes[0, 1].plot(
            x_profile, mach, marker="o", ms=3.0,
            color=colour, linestyle=style, label=f"{time_us:.1f} us"
        )

    axes[0, 0].set(
        xlabel="x [mm]", ylabel="pressure [kPa]",
        title="i=4 pressure profile, j=168...174"
    )
    axes[0, 1].set(
        xlabel="x [mm]", ylabel="Mach",
        title="i=4 Mach profile, j=168...174"
    )
    axes[0, 0].legend(ncol=2, fontsize=8)

    pmin_j = np.asarray([int(row["local_pressure_min_j"]) for row in rows])
    machmax_j = np.asarray([int(row["local_mach_max_j"]) for row in rows])
    axes[1, 0].step(times, pmin_j, where="post", marker="o", label="pressure minimum j")
    axes[1, 0].step(times, machmax_j, where="post", marker="s", label="Mach maximum j")
    axes[1, 0].set(
        xlabel="time [microseconds]", ylabel="cell index j",
        title="Movement of the local structure"
    )
    axes[1, 0].set_yticks(PROFILE_JS)
    axes[1, 0].legend(fontsize=8)

    reconstructed_pressure = np.asarray(
        [float(row["right_reconstructed_pressure_Pa"]) for row in rows]
    )
    candidate0_pressure = np.asarray(
        [float(row["right_candidate0_pressure_Pa"]) for row in rows]
    )
    candidate0_weight = np.asarray(
        [float(row["right_wminus_candidate0_weight"]) for row in rows]
    )
    axes[1, 1].semilogy(
        times, reconstructed_pressure, marker="o", label="final reconstructed p_R"
    )
    axes[1, 1].semilogy(
        times, candidate0_pressure, marker="s", label="paired candidate-0 p_R"
    )
    axes[1, 1].set(
        xlabel="time [microseconds]", ylabel="pressure [Pa]",
        title="Face 171 candidate deterioration"
    )
    weight_axis = axes[1, 1].twinx()
    weight_axis.plot(
        times, candidate0_weight, color="tab:red", marker="^",
        linestyle=":", label="w_minus omega_0"
    )
    weight_axis.set_ylabel("candidate-0 weight", color="tab:red")
    weight_axis.tick_params(axis="y", colors="tab:red")
    lines, labels = axes[1, 1].get_legend_handles_labels()
    lines2, labels2 = weight_axis.get_legend_handles_labels()
    axes[1, 1].legend(lines + lines2, labels + labels2, fontsize=8)

    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.25)
    axes[1, 0].axvline(73.6, color="0.5", linestyle=":", linewidth=0.8)
    axes[1, 1].axvline(73.6, color="0.5", linestyle=":", linewidth=0.8)

    fig.suptitle(
        "AFD fallback-off: physical peak motion versus face-reconstruction failure",
        fontsize=13,
    )
    fig.savefig(FIGURE_PATH, dpi=220)
    print(CSV_PATH)
    print(FIGURE_PATH)


if __name__ == "__main__":
    main()
