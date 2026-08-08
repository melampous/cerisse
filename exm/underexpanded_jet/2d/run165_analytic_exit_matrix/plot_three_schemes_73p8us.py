#!/usr/bin/env python3
"""Compare the three Run165 RZ schemes at the exact 73.8 microsecond stage."""

from __future__ import annotations

import csv
import inspect
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, TwoSlopeNorm
import numpy as np
import yt
from yt.frontends.boxlib import data_structures as boxlib_data


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_three_at_73p8us"
FIELDS = RESULTS / "stage_fields"
FIGURES = RESULTS / "figures"
GAMMA = 1.4
P_INF = 534.5806702256042
FAILURE_CELL = (4, 171)

CASES = {
    "llfweno": (
        FIELDS / "llfweno_q2_73p8us",
        r"LLF-WENO-Z5 ($q=2$)",
        "#0072B2",
    ),
    "hllcweno": (
        FIELDS / "afd_hllcweno_q2_73p8us",
        r"AFD-HLLC-WENO-Z5 ($q=2$, fallback off)",
        "#D55E00",
    ),
    "hllcmuscl": (
        FIELDS / "hllc_muscl_73p8us",
        "HLLC-MUSCL O2",
        "#009E73",
    ),
}


def patch_yt_cylindrical_readonly_edge() -> None:
    """Work around the yt/NumPy read-only cylindrical-domain edge issue."""
    source = textwrap.dedent(
        inspect.getsource(boxlib_data.BoxlibDataset._parse_header_file)
    )
    source = source.replace(
        "dre = self.domain_right_edge\n",
        "dre = self.domain_right_edge.copy()\n",
    )
    namespace = dict(boxlib_data.__dict__)
    exec(source, namespace)
    boxlib_data.BoxlibDataset._parse_header_file = namespace[
        "_parse_header_file"
    ]


def load_stage(plotfile: Path) -> dict[str, np.ndarray | float]:
    dataset = yt.load(str(plotfile))
    if int(dataset.index.max_level) != 0:
        raise RuntimeError(f"Expected a standalone stage level: {plotfile}")
    dims = np.asarray(dataset.domain_dimensions, dtype=int).copy()
    dims[2] = 1
    grid = dataset.covering_grid(0, dataset.domain_left_edge, dims)

    def field(name: str) -> np.ndarray:
        return np.asarray(
            grid[("boxlib", name)].to_ndarray()[:, :, 0], dtype=np.float64
        )

    rho = field("Density")
    mx = field("Xmom")
    my = field("Ymom")
    mz = field("Zmom")
    energy = field("Energy")
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        kinetic = 0.5 * (mx * mx + my * my + mz * mz) / rho
        rhoe = energy - kinetic
        pressure = (GAMMA - 1.0) * rhoe
        ur = mx / rho
        uz = my / rho
        sound_speed = np.sqrt(GAMMA * pressure / rho)
        mach = np.hypot(ur, uz) / sound_speed
    admissible = (
        np.isfinite(rho)
        & np.isfinite(rhoe)
        & (rho > 0.0)
        & (rhoe > 0.0)
    )
    mach = np.where(admissible, mach, np.nan)

    nr, nz = rho.shape
    r_edges = np.linspace(
        float(dataset.domain_left_edge[0]),
        float(dataset.domain_right_edge[0]),
        nr + 1,
    ) * 1.0e3
    z_edges = np.linspace(
        float(dataset.domain_left_edge[1]),
        float(dataset.domain_right_edge[1]),
        nz + 1,
    ) * 1.0e3
    r_centres = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    return {
        "rho": rho,
        "mx": mx,
        "my": my,
        "mz": mz,
        "energy": energy,
        "rhoe": rhoe,
        "pressure": pressure,
        "ur": ur,
        "uz": uz,
        "mach": mach,
        "admissible": admissible,
        "r_edges": r_edges,
        "z_edges": z_edges,
        "r_centres": r_centres,
        "z_centres": z_centres,
        "time_us": float(dataset.current_time) * 1.0e6,
    }


def mark_failure(axis, data: dict[str, np.ndarray | float]) -> None:
    i, j = FAILURE_CELL
    axis.scatter(
        [data["z_centres"][j]],
        [data["r_centres"][i]],
        marker="x",
        s=55,
        linewidths=1.7,
        color="#00FFFF",
        zorder=8,
    )


def add_invalid_markers(axis, data: dict[str, np.ndarray | float]) -> None:
    ii, jj = np.nonzero(~data["admissible"])
    if ii.size:
        axis.scatter(
            data["z_centres"][jj],
            data["r_centres"][ii],
            marker="s",
            s=38,
            facecolors="none",
            edgecolors="#FF00FF",
            linewidths=1.3,
            zorder=7,
        )


def plot_full(states: dict[str, dict[str, np.ndarray | float]]) -> Path:
    fig, axes = plt.subplots(3, 2, figsize=(14.2, 11.0), constrained_layout=True)
    mach_image = pressure_image = None
    for row, (key, (_, label, _)) in enumerate(CASES.items()):
        data = states[key]
        mach_image = axes[row, 0].pcolormesh(
            data["z_edges"],
            data["r_edges"],
            data["mach"],
            shading="auto",
            cmap="turbo",
            vmin=0.0,
            vmax=6.0,
        )
        pressure_ratio = np.ma.masked_where(
            ~data["admissible"], data["pressure"] / P_INF
        )
        pressure_image = axes[row, 1].pcolormesh(
            data["z_edges"],
            data["r_edges"],
            pressure_ratio,
            shading="auto",
            cmap="magma",
            norm=LogNorm(vmin=0.5, vmax=250.0),
        )
        for col in range(2):
            axes[row, col].set_xlim(-201.3277, 2.4413)
            axes[row, col].set_ylim(0.0, 67.9230)
            axes[row, col].set_ylabel("r [mm]")
            mark_failure(axes[row, col], data)
            add_invalid_markers(axes[row, col], data)
        axes[row, 0].set_title(f"{label}\nMach number")
        axes[row, 1].set_title(f"{label}\nStatic-pressure ratio $p/p_\\infty$")
    axes[-1, 0].set_xlabel("axial x [mm]")
    axes[-1, 1].set_xlabel("axial x [mm]")
    fig.colorbar(mach_image, ax=axes[:, 0], label="Mach", shrink=0.92)
    fig.colorbar(
        pressure_image, ax=axes[:, 1], label=r"$p/p_\infty$", shrink=0.92
    )
    fig.suptitle(
        "Run165 RZ analytic-exit jet — exact fine-stage comparison at t = 73.8 μs\n"
        "cyan ×: common reference cell (4,171); magenta square: inadmissible cell",
        fontsize=14,
    )
    output = FIGURES / "run165_three_schemes_73p8us_full_mach_pressure.png"
    fig.savefig(output, dpi=210)
    plt.close(fig)
    return output


def plot_zoom(states: dict[str, dict[str, np.ndarray | float]]) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(16.2, 10.7), constrained_layout=True)
    pressure_image = velocity_image = rhoe_image = None
    for row, (key, (_, label, _)) in enumerate(CASES.items()):
        data = states[key]
        pressure_ratio = np.ma.masked_where(
            ~data["admissible"], data["pressure"] / P_INF
        )
        pressure_image = axes[row, 0].pcolormesh(
            data["z_edges"], data["r_edges"], pressure_ratio,
            shading="auto", cmap="magma", norm=LogNorm(0.5, 250.0)
        )
        velocity_image = axes[row, 1].pcolormesh(
            data["z_edges"], data["r_edges"], data["uz"],
            shading="auto", cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-1500.0, vcenter=0.0, vmax=1500.0),
        )
        rhoe_image = axes[row, 2].pcolormesh(
            data["z_edges"], data["r_edges"], data["rhoe"],
            shading="auto", cmap="coolwarm",
            norm=SymLogNorm(
                linthresh=1.0e3, linscale=0.8, vmin=-4.0e7, vmax=4.0e7,
                base=10.0,
            ),
        )
        axes[row, 0].set_title(f"{label}\n$p/p_\\infty$ (inadmissible masked)")
        axes[row, 1].set_title(f"{label}\naxial velocity $u_x$ [m/s]")
        axes[row, 2].set_title(f"{label}\ninternal-energy density $\\rho e$ [J/m³]")
        for col in range(3):
            axes[row, col].set_xlim(-32.0, -8.0)
            axes[row, col].set_ylim(0.0, 15.0)
            axes[row, col].set_ylabel("r [mm]")
            mark_failure(axes[row, col], data)
            add_invalid_markers(axes[row, col], data)
    for axis in axes[-1, :]:
        axis.set_xlabel("axial x [mm]")
    fig.colorbar(pressure_image, ax=axes[:, 0], label=r"$p/p_\infty$", shrink=0.92)
    fig.colorbar(velocity_image, ax=axes[:, 1], label=r"$u_x$ [m/s]", shrink=0.92)
    fig.colorbar(rhoe_image, ax=axes[:, 2], label=r"$\rho e$ [J/m³]", shrink=0.92)
    fig.suptitle(
        "Failure-region comparison at t = 73.8 μs — same grid, same timestep and jet ramp",
        fontsize=14,
    )
    output = FIGURES / "run165_three_schemes_73p8us_failure_zoom.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def plot_profiles(states: dict[str, dict[str, np.ndarray | float]]) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.8), constrained_layout=True)
    i, j_failure = FAILURE_CELL
    specs = (
        ("rho", r"density $\rho$ [kg/m³]", "linear"),
        ("pressure", r"pressure $p$ [Pa]", "symlog"),
        ("uz", r"axial velocity $u_x$ [m/s]", "symlog"),
        ("rhoe", r"internal-energy density $\rho e$ [J/m³]", "symlog"),
    )
    for axis, (field, ylabel, scale) in zip(axes.flat, specs):
        for key, (_, label, color) in CASES.items():
            data = states[key]
            axis.plot(
                data["z_centres"], data[field][i, :], label=label,
                color=color, linewidth=1.45,
            )
        reference = states["hllcweno"]
        axis.axvline(
            reference["z_centres"][j_failure], color="black", linestyle="--",
            linewidth=0.9,
        )
        axis.set_xlim(-32.0, -8.0)
        axis.set_xlabel("axial x [mm]")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.22)
        if scale == "symlog":
            axis.set_yscale("symlog", linthresh=1.0e3)
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle(
        "Line profiles through radial ring i=4 (r=4.776 mm), t = 73.8 μs\n"
        "dashed line: AFD-HLLC-WENO first inadmissible cell j=171",
        fontsize=14,
    )
    output = FIGURES / "run165_three_schemes_73p8us_i4_profiles.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def write_metrics(states: dict[str, dict[str, np.ndarray | float]]) -> Path:
    rows: list[dict[str, float | int | str]] = []
    i, j = FAILURE_CELL
    for key, (_, label, _) in CASES.items():
        data = states[key]
        valid_mach = data["mach"][data["admissible"]]
        rows.append(
            {
                "scheme": key,
                "label": label,
                "time_us": data["time_us"],
                "invalid_cell_count": int(np.count_nonzero(~data["admissible"])),
                "min_density_kg_m3": float(np.nanmin(data["rho"])),
                "min_pressure_Pa": float(np.nanmin(data["pressure"])),
                "min_rhoe_J_m3": float(np.nanmin(data["rhoe"])),
                "max_admissible_Mach": float(np.nanmax(valid_mach)),
                "cell_4_171_density_kg_m3": float(data["rho"][i, j]),
                "cell_4_171_pressure_Pa": float(data["pressure"][i, j]),
                "cell_4_171_uz_m_s": float(data["uz"][i, j]),
                "cell_4_171_rhoe_J_m3": float(data["rhoe"][i, j]),
            }
        )
    output = RESULTS / "run165_three_schemes_73p8us_metrics.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)
    return output


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    patch_yt_cylindrical_readonly_edge()
    yt.funcs.mylog.setLevel(50)
    states = {key: load_stage(entry[0]) for key, entry in CASES.items()}
    times = [float(data["time_us"]) for data in states.values()]
    if max(times) - min(times) > 1.0e-10 or abs(times[0] - 73.8) > 1.0e-10:
        raise RuntimeError(f"Stage times are not exactly matched: {times}")
    outputs = [
        plot_full(states),
        plot_zoom(states),
        plot_profiles(states),
        write_metrics(states),
    ]
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
