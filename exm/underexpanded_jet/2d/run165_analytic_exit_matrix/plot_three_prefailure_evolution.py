#!/usr/bin/env python3
"""Trace the AFD-HLLC-WENO precursor against LLF-WENO and HLLC-MUSCL."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_hllc_muscl_afd_weno_pointwise as standard_plot
import plot_three_schemes_73p8us as stage_plot


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_three_pre_failure"
FIGURES = RESULTS / "figures"
STAGES = RESULTS / "stage_fields"
I_RING = 4
J_CELLS = (170, 171, 172)

SCHEMES = {
    "llfweno": (r"LLF-WENO-Z5 ($q=2$)", "#0072B2"),
    "hllcweno": (
        r"AFD-HLLC-WENO-Z5 ($q=2$, fallback off)",
        "#D55E00",
    ),
    "hllcmuscl": ("HLLC-MUSCL O2", "#009E73"),
}


def coarse_path(scheme: str, step: int) -> Path:
    onset_names = {
        "llfweno": "llf",
        "hllcweno": "afd",
        "hllcmuscl": "muscl",
    }
    if 126 <= step <= 150:
        return (
            RESULTS
            / f"onset_replay/{onset_names[scheme]}/plt{step:05d}"
        )
    if scheme == "llfweno":
        return RESULTS / f"llfweno_replay/plt{step:05d}"
    if scheme == "hllcweno":
        return (
            ROOT
            / f"compare_hllc_muscl_vs_afd_weno/afd_fallback_off/plt{step:05d}"
        )
    return ROOT / f"compare_hllc_muscl_vs_afd_weno/muscl/plt{step:05d}"


def final_path(scheme: str) -> Path:
    names = {
        "llfweno": "llfweno_q2_73p8us",
        "hllcweno": "afd_hllcweno_q2_73p8us",
        "hllcmuscl": "hllc_muscl_73p8us",
    }
    return ROOT / "compare_three_at_73p8us/stage_fields" / names[scheme]


def load_standard(path: Path) -> dict[str, np.ndarray | float]:
    fields, time_us, r_edges, z_edges, r_centres, z_centres = (
        standard_plot.load_fields(path)
    )
    fields.update(
        {
            "time_us": time_us,
            "r_edges": r_edges,
            "z_edges": z_edges,
            "r_centres": r_centres,
            "z_centres": z_centres,
            "admissible": (fields["rho"] > 0.0) & (fields["rhoe"] > 0.0),
        }
    )
    return fields


def load_all():
    coarse: dict[str, dict[int, dict[str, np.ndarray | float]]] = {
        scheme: {} for scheme in SCHEMES
    }
    for scheme in SCHEMES:
        for step in range(170, 185):
            coarse[scheme][step] = load_standard(coarse_path(scheme, step))

    fine: dict[str, dict[float, dict[str, np.ndarray | float]]] = {
        scheme: {} for scheme in SCHEMES
    }
    for scheme in SCHEMES:
        for time_tag, time_us in (
            ("72p2", 72.2),
            ("72p6", 72.6),
            ("73p0", 73.0),
            ("73p4", 73.4),
        ):
            fine[scheme][time_us] = stage_plot.load_stage(
                STAGES / f"{scheme}_{time_tag}us"
            )
        fine[scheme][73.8] = stage_plot.load_stage(
            STAGES / f"{scheme}_73p8us_fe"
        )

    final = {
        scheme: stage_plot.load_stage(final_path(scheme))
        for scheme in SCHEMES
    }
    return coarse, fine, final


def chronological_states(coarse, fine, scheme: str):
    states = [coarse[scheme][step] for step in range(170, 185)]
    states.extend(fine[scheme].values())
    return sorted(states, key=lambda state: float(state["time_us"]))


def state_at(coarse, fine, scheme: str, time_us: float):
    coarse_step = round(time_us / 0.4)
    if abs(0.4 * coarse_step - time_us) < 1.0e-9:
        return coarse[scheme][coarse_step]
    return fine[scheme][time_us]


def plot_time_history(coarse, fine, final) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.8), constrained_layout=True)
    for axis, j in zip(axes.flat[:3], J_CELLS):
        for scheme, (label, color) in SCHEMES.items():
            states = chronological_states(coarse, fine, scheme)
            times = [state["time_us"] for state in states]
            pressure = [state["pressure"][I_RING, j] / 1.0e3 for state in states]
            axis.plot(times, pressure, "o-", ms=3.1, lw=1.35,
                      color=color, label=label)
        axis.set_title(f"cell (i,j)=({I_RING},{j})")
        axis.set_ylabel("pressure [kPa]")
        axis.set_xlim(68.0, 73.9)
        axis.grid(alpha=0.22)
    velocity_axis = axes[1, 1]
    for scheme, (label, color) in SCHEMES.items():
        states = chronological_states(coarse, fine, scheme)
        velocity_axis.plot(
            [state["time_us"] for state in states],
            [state["uz"][I_RING, 171] for state in states],
            "o-", ms=3.1, lw=1.35, color=color, label=label,
        )
    velocity_axis.set_title(f"cell (i,j)=({I_RING},171)")
    velocity_axis.set_ylabel("axial velocity [m/s]")
    velocity_axis.set_xlim(68.0, 73.9)
    velocity_axis.grid(alpha=0.22)
    for axis in axes.flat:
        axis.set_xlabel("time [μs]")
        axis.axvspan(72.0, 73.8, color="#D55E00", alpha=0.055)
        axis.axvline(73.8, color="black", ls="--", lw=0.8)
    axes[0, 0].legend(fontsize=8, loc="best")
    axes[0, 1].annotate(
        "HLLC-WENO Heun final:\n−14.54 MPa (off scale)",
        xy=(73.8, fine["hllcweno"][73.8]["pressure"][I_RING, 171] / 1.0e3),
        xytext=(72.15, 15.0), arrowprops={"arrowstyle": "->", "color": "#D55E00"},
        color="#D55E00", fontsize=9,
    )
    velocity_axis.annotate(
        "HLLC-WENO Heun final:\n+26.7 km/s (off scale)",
        xy=(73.8, fine["hllcweno"][73.8]["uz"][I_RING, 171]),
        xytext=(71.9, -620.0), arrowprops={"arrowstyle": "->", "color": "#D55E00"},
        color="#D55E00", fontsize=9,
    )
    fig.suptitle(
        "Pre-failure evolution on radial ring i=4 (r=4.776 mm)\n"
        "73.8 μs points on the curves are admissible first-stage FE predictors",
        fontsize=14,
    )
    output = FIGURES / "run165_hllcweno_prefailure_cell_history.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def plot_profile_sequence(coarse, fine) -> Path:
    selected = (68.0, 70.0, 71.6, 72.6, 73.4, 73.8)
    fig, axes = plt.subplots(2, 3, figsize=(15.2, 8.2), constrained_layout=True)
    for axis, time_us in zip(axes.flat, selected):
        for scheme, (label, color) in SCHEMES.items():
            state = state_at(coarse, fine, scheme, time_us)
            axis.plot(
                state["z_centres"], state["pressure"][I_RING, :] / 1.0e3,
                color=color, lw=1.55, label=label,
            )
        suffix = " (FE predictor)" if time_us == 73.8 else ""
        axis.set_title(f"t = {time_us:.1f} μs{suffix}")
        axis.set_xlim(-30.0, -12.0)
        axis.set_ylim(0.0, 45.0)
        axis.set_xlabel("axial x [mm]")
        axis.set_ylabel("pressure [kPa]")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7.8, loc="best")
    fig.suptitle(
        "Pressure-profile evolution at r=4.776 mm — moving peak/trough precursor",
        fontsize=14,
    )
    output = FIGURES / "run165_hllcweno_prefailure_pressure_profiles.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def plot_rk_jump(coarse, fine, final) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(15.8, 10.8), constrained_layout=True)
    fields = (
        ("pressure", "pressure [Pa]", 1.0e4, (-2.0e7, 5.0e6)),
        ("uz", "axial velocity [m/s]", 1.0e3, (-4.0e4, 4.0e4)),
        ("rhoe", r"internal-energy density $\rho e$ [J/m³]", 2.0e4,
         (-5.0e7, 2.0e7)),
    )
    for row, (scheme, (label, _)) in enumerate(SCHEMES.items()):
        states = (
            (coarse[scheme][184], "73.6 μs accepted", "#555555", "o-"),
            (fine[scheme][73.8], "73.8 μs FE predictor", "#0072B2", "s-"),
            (final[scheme], "73.8 μs Heun final", "#D55E00", "x-"),
        )
        for col, (field, ylabel, linthresh, limits) in enumerate(fields):
            axis = axes[row, col]
            for state, stage_label, color, style in states:
                axis.plot(
                    state["z_centres"][168:175],
                    state[field][I_RING, 168:175],
                    style, color=color, lw=1.45, ms=4.2, label=stage_label,
                )
            axis.axvline(
                final[scheme]["z_centres"][171], color="black", ls="--", lw=0.8
            )
            axis.set_yscale("symlog", linthresh=linthresh)
            axis.set_ylim(*limits)
            axis.set_xlabel("axial x [mm]")
            axis.set_ylabel(ylabel)
            axis.grid(alpha=0.2)
            axis.set_title(f"{label}\n{field}")
    axes[0, 0].legend(fontsize=7.5, loc="best")
    fig.suptitle(
        "The catastrophic spike is created by the second RHS / Heun corrector\n"
        "dashed line: cell j=171; all three FE predictors are still admissible",
        fontsize=14,
    )
    output = FIGURES / "run165_three_schemes_73p8us_rk_stage_jump.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def plot_moving_trough(fine) -> Path:
    records = []
    for step in range(140, 185):
        states = {
            scheme: load_standard(coarse_path(scheme, step))
            for scheme in SCHEMES
        }
        records.append((float(states["hllcweno"]["time_us"]), states))
    for time_us in (72.2, 72.6, 73.0, 73.4, 73.8):
        records.append(
            (time_us, {scheme: fine[scheme][time_us] for scheme in SCHEMES})
        )
    records.sort(key=lambda item: item[0])

    times = []
    locations = []
    indices = []
    pressures = {scheme: [] for scheme in SCHEMES}
    deficits = []
    for time_us, states in records:
        afd_pressure = states["hllcweno"]["pressure"][I_RING, :]
        candidates = [
            j
            for j in range(168, 185)
            if afd_pressure[j] > 1.0e3
            and afd_pressure[j] < afd_pressure[j - 1]
            and afd_pressure[j] < afd_pressure[j + 1]
        ]
        j = min(candidates, key=lambda index: afd_pressure[index])
        times.append(time_us)
        indices.append(j)
        locations.append(states["hllcweno"]["z_centres"][j])
        for scheme in SCHEMES:
            pressures[scheme].append(states[scheme]["pressure"][I_RING, j] / 1.0e3)
        baseline = 0.5 * (
            states["llfweno"]["pressure"][I_RING, j]
            + states["hllcmuscl"]["pressure"][I_RING, j]
        )
        deficits.append(100.0 * (afd_pressure[j] - baseline) / baseline)

    fig, axes = plt.subplots(3, 1, figsize=(11.2, 10.0), constrained_layout=True,
                             sharex=True)
    for scheme, (label, color) in SCHEMES.items():
        axes[0].plot(times, pressures[scheme], "o-", ms=3.0, lw=1.35,
                     color=color, label=label)
    axes[0].set_ylabel("pressure at moving trough [kPa]")
    axes[0].legend(fontsize=8)
    axes[1].plot(times, deficits, "o-", ms=3.0, lw=1.4, color="#D55E00")
    axes[1].axhline(-5.0, color="black", ls="--", lw=0.8)
    axes[1].set_ylabel("HLLC-WENO deficit\nvs mean of references [%]")
    axes[2].step(times, locations, where="mid", color="#7A3E9D", lw=1.5)
    axes[2].scatter(times, locations, c=indices, cmap="viridis", s=17)
    axes[2].set_ylabel("trough axial x [mm]")
    axes[2].set_xlabel("time [μs]")
    for axis in axes:
        axis.grid(alpha=0.22)
        axis.axvline(56.8, color="#D55E00", ls=":", lw=1.0)
        axis.axvline(73.8, color="black", ls="--", lw=0.8)
    axes[1].annotate(
        "first systematic undershoot\n≈56.8 μs",
        xy=(56.8, deficits[int(np.argmin(np.abs(np.asarray(times) - 56.8)))]),
        xytext=(58.0, -12.0),
        arrowprops={"arrowstyle": "->", "color": "#D55E00"},
        color="#D55E00", fontsize=9,
    )
    fig.suptitle(
        "Onset and transport of the HLLC-WENO pressure undershoot\n"
        "the tracked cell changes as the rear jet wave travels upstream",
        fontsize=14,
    )
    output = FIGURES / "run165_hllcweno_moving_trough_onset.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def write_csv(coarse, fine, final) -> Path:
    rows = []
    for scheme in SCHEMES:
        for state in chronological_states(coarse, fine, scheme):
            row = {
                "scheme": scheme,
                "state": "accepted_or_fe_predictor",
                "time_us": float(state["time_us"]),
            }
            for j in J_CELLS:
                row[f"p_j{j}_Pa"] = float(state["pressure"][I_RING, j])
                row[f"uz_j{j}_m_s"] = float(state["uz"][I_RING, j])
                row[f"rhoe_j{j}_J_m3"] = float(state["rhoe"][I_RING, j])
            rows.append(row)
        state = final[scheme]
        row = {"scheme": scheme, "state": "heun_final", "time_us": 73.8}
        for j in J_CELLS:
            row[f"p_j{j}_Pa"] = float(state["pressure"][I_RING, j])
            row[f"uz_j{j}_m_s"] = float(state["uz"][I_RING, j])
            row[f"rhoe_j{j}_J_m3"] = float(state["rhoe"][I_RING, j])
        rows.append(row)
    output = RESULTS / "run165_three_schemes_prefailure_evolution.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    stage_plot.patch_yt_cylindrical_readonly_edge()
    stage_plot.yt.funcs.mylog.setLevel(50)
    coarse, fine, final = load_all()
    outputs = (
        plot_time_history(coarse, fine, final),
        plot_profile_sequence(coarse, fine),
        plot_moving_trough(fine),
        plot_rk_jump(coarse, fine, final),
        write_csv(coarse, fine, final),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
