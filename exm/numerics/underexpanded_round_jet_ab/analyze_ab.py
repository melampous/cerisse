#!/usr/bin/env python3
"""Analyse the reduced 3D underexpanded round-jet shock-sensor A/B test."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import yt
except ImportError as exc:
    raise SystemExit("analyze_ab.py requires yt") from exc


GAMMA = 1.4
P_AMBIENT = 101325.0
D_EXIT = 0.010
PRESSURE_THRESHOLD = 0.01
COMPRESSION_THRESHOLD = 0.001
EXPECTED_DIMS = np.array([96, 64, 64])
EXPECTED_LO = np.array([0.0, -0.02, -0.02])
EXPECTED_HI = np.array([0.06, 0.02, 0.02])


def latest_plotfile(path: Path) -> Path:
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        item
        for item in path.glob("plot/plt*")
        if pattern.fullmatch(item.name) and (item / "Header").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"no plotfile below {path}")
    return max(candidates, key=lambda item: int(pattern.fullmatch(item.name)[1]))


def field_key(dataset, name: str):
    wanted = name.lower()
    for key in dataset.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"{name!r} is absent from {dataset.field_list}")


def load_case(case_dir: Path) -> dict[str, object]:
    plotfile = latest_plotfile(case_dir)
    dataset = yt.load(str(plotfile))
    dims = np.asarray(dataset.domain_dimensions, dtype=int)
    lo = np.asarray(dataset.domain_left_edge.d, dtype=float)
    hi = np.asarray(dataset.domain_right_edge.d, dtype=float)
    if int(dataset.index.max_level) != 0:
        raise ValueError(f"{plotfile} is not Level-0-only")
    if not np.array_equal(dims, EXPECTED_DIMS):
        raise ValueError(f"{plotfile} has dimensions {dims}, expected {EXPECTED_DIMS}")
    if not (
        np.allclose(lo, EXPECTED_LO, rtol=0.0, atol=1.0e-13)
        and np.allclose(hi, EXPECTED_HI, rtol=0.0, atol=1.0e-13)
    ):
        raise ValueError(f"{plotfile} has unexpected domain {lo} to {hi}")

    grid = dataset.covering_grid(0, dataset.domain_left_edge, dims)

    def values(name: str) -> np.ndarray:
        return np.asarray(grid[field_key(dataset, name)].d, dtype=float)

    rho = values("Density")
    xmom = values("Xmom")
    ymom = values("Ymom")
    zmom = values("Zmom")
    energy = values("Energy")
    ux = xmom / rho
    uy = ymom / rho
    uz = zmom / rho
    pressure = (GAMMA - 1.0) * (
        energy - 0.5 * (xmom * xmom + ymom * ymom + zmom * zmom) / rho
    )
    sound = np.sqrt(GAMMA * pressure / rho)
    mach = np.sqrt(ux * ux + uy * uy + uz * uz) / sound
    dx = (hi - lo) / dims
    x = lo[0] + (np.arange(dims[0]) + 0.5) * dx[0]
    y = lo[1] + (np.arange(dims[1]) + 0.5) * dx[1]
    z = lo[2] + (np.arange(dims[2]) + 0.5) * dx[2]
    return {
        "plotfile": str(plotfile),
        "time": float(dataset.current_time),
        "x": x,
        "y": y,
        "z": z,
        "dx": dx,
        "rho": rho,
        "pressure": pressure,
        "ux": ux,
        "uy": uy,
        "uz": uz,
        "mach": mach,
        "sound": sound,
    }


def central_profile(field: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    iy = np.argsort(np.abs(y))[:2]
    iz = np.argsort(np.abs(z))[:2]
    return np.mean(field[:, iy[:, None], iz[None, :]], axis=(1, 2))


def strongest_compression(case: dict[str, object]) -> dict[str, float]:
    x = np.asarray(case["x"])
    p = central_profile(
        np.asarray(case["pressure"]), np.asarray(case["y"]), np.asarray(case["z"])
    )
    mach = central_profile(
        np.asarray(case["mach"]), np.asarray(case["y"]), np.asarray(case["z"])
    )
    dpdx = np.gradient(p, x)
    window = (x >= 0.4 * D_EXIT) & (x <= 2.5 * D_EXIT)
    candidates = np.flatnonzero(window)
    index = int(candidates[np.argmax(dpdx[candidates])])
    peak = float(dpdx[index])
    half = 0.5 * peak
    left = index
    right = index
    while left > candidates[0] and dpdx[left - 1] >= half:
        left -= 1
    while right < candidates[-1] and dpdx[right + 1] >= half:
        right += 1
    upstream = slice(max(index - 2, 0), index)
    downstream = slice(index + 1, min(index + 3, x.size))
    return {
        "compression_index": index,
        "compression_x_over_De": float(x[index] / D_EXIT),
        "compression_dpdx_Pa_per_m": peak,
        "compression_halfmax_width_over_De": float(
            (right - left + 1) * np.asarray(case["dx"])[0] / D_EXIT
        ),
        "centreline_pressure_upstream_over_pa": float(np.mean(p[upstream]) / P_AMBIENT),
        "centreline_pressure_downstream_over_pa": float(
            np.mean(p[downstream]) / P_AMBIENT
        ),
        "centreline_mach_upstream": float(np.mean(mach[upstream])),
        "centreline_mach_downstream": float(np.mean(mach[downstream])),
    }


def symmetry_metrics(case: dict[str, object]) -> dict[str, float]:
    pressure = np.asarray(case["pressure"])
    y = np.asarray(case["y"])
    z = np.asarray(case["z"])
    iy = np.argsort(np.abs(y))[:2]
    iz = np.argsort(np.abs(z))[:2]
    xy = np.mean(pressure[:, :, iz], axis=2)
    xz = np.mean(pressure[:, iy, :], axis=1)
    difference = (xy - xz) / P_AMBIENT
    return {
        "yz_symmetry_pressure_L2_over_pa": float(np.sqrt(np.mean(difference**2))),
        "yz_symmetry_pressure_Linf_over_pa": float(np.max(np.abs(difference))),
    }


def odd_even_metric(case: dict[str, object]) -> dict[str, float]:
    pressure = np.asarray(case["pressure"])
    x = np.asarray(case["x"])
    y = np.asarray(case["y"])
    z = np.asarray(case["z"])
    yy, zz = np.meshgrid(y, z, indexing="ij")
    core = yy * yy + zz * zz <= (0.75 * D_EXIT) ** 2
    parity_y = (-1.0) ** np.arange(y.size)[:, None]
    parity_z = (-1.0) ** np.arange(z.size)[None, :]
    axial = (x >= 0.4 * D_EXIT) & (x <= 2.5 * D_EXIT)
    amp_y = []
    amp_z = []
    for plane in pressure[axial]:
        centred = plane - np.mean(plane[core])
        amp_y.append(abs(np.mean((centred * parity_y)[core])) / P_AMBIENT)
        amp_z.append(abs(np.mean((centred * parity_z)[core])) / P_AMBIENT)
    return {
        "nearfield_odd_even_y_max_over_pa": float(max(amp_y)),
        "nearfield_odd_even_z_max_over_pa": float(max(amp_z)),
    }


def pressure_compression_gate_proxy(
    case: dict[str, object],
) -> dict[str, float | int]:
    pressure = np.asarray(case["pressure"])
    ux = np.asarray(case["ux"])
    uy = np.asarray(case["uy"])
    uz = np.asarray(case["uz"])
    sound = np.asarray(case["sound"])
    x = np.asarray(case["x"])
    y = np.asarray(case["y"])
    z = np.asarray(case["z"])
    cell_size = np.asarray(case["dx"])
    tiny = np.finfo(float).tiny

    pressure_jump = np.zeros_like(pressure)
    divergence = np.zeros_like(pressure)
    sound_scale = sound.copy()
    for axis, velocity in enumerate((ux, uy, uz)):
        p_minus = np.roll(pressure, 1, axis=axis)
        p_plus = np.roll(pressure, -1, axis=axis)
        p_scale = np.maximum(np.maximum(np.abs(p_minus), np.abs(p_plus)), tiny)
        pressure_jump = np.maximum(pressure_jump, np.abs(p_plus - p_minus) / p_scale)
        u_minus = np.roll(velocity, 1, axis=axis)
        u_plus = np.roll(velocity, -1, axis=axis)
        c_minus = np.roll(sound, 1, axis=axis)
        c_plus = np.roll(sound, -1, axis=axis)
        divergence += (u_plus - u_minus) / (2.0 * cell_size[axis])
        sound_scale = np.maximum(sound_scale, np.maximum(c_minus, c_plus))
    compression = (
        np.maximum(-divergence, 0.0)
        * np.min(cell_size)
        / np.maximum(sound_scale, tiny)
    )

    # Remove wrapped neighbours introduced by np.roll.
    interior = np.zeros_like(pressure, dtype=bool)
    interior[1:-1, 1:-1, 1:-1] = True
    p_gate = pressure_jump > PRESSURE_THRESHOLD
    c_gate = compression > COMPRESSION_THRESHOLD
    pressure_compression = interior & p_gate & c_gate

    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    radius = np.sqrt(yy * yy + zz * zz)
    shock_core = (
        interior
        & (xx >= 0.4 * D_EXIT)
        & (xx <= 2.5 * D_EXIT)
        & (radius <= 0.35 * D_EXIT)
    )
    shear_annulus = (
        interior
        & (xx <= 2.0 * D_EXIT)
        & (np.abs(radius - 0.5 * D_EXIT) <= 2.0 * np.asarray(case["dx"])[1])
    )

    def fraction(mask: np.ndarray, region: np.ndarray) -> float:
        count = int(np.count_nonzero(region))
        return float(np.count_nonzero(mask & region) / count) if count else math.nan

    return {
        "offline_pressure_compression_gate_proxy_cells":
            int(np.count_nonzero(pressure_compression)),
        "offline_pressure_compression_gate_proxy_fraction":
            fraction(pressure_compression, interior),
        "offline_pressure_gate_shear_fraction": fraction(p_gate, shear_annulus),
        "offline_compression_gate_shear_fraction": fraction(c_gate, shear_annulus),
        "offline_pressure_compression_gate_shear_fraction":
            fraction(pressure_compression, shear_annulus),
        "offline_pressure_gate_shock_core_fraction": fraction(p_gate, shock_core),
        "offline_compression_gate_shock_core_fraction": fraction(c_gate, shock_core),
        "offline_pressure_compression_gate_shock_core_fraction":
            fraction(pressure_compression, shock_core),
    }


def case_metrics(label: str, case: dict[str, object]) -> dict[str, object]:
    rho = np.asarray(case["rho"])
    pressure = np.asarray(case["pressure"])
    mach = np.asarray(case["mach"])
    values: dict[str, object] = {
        "variant": label,
        "afd_shock_llf": int(label[-1]),
        "plotfile": str(case["plotfile"]),
        "time_s": float(case["time"]),
        "density_min": float(np.min(rho)),
        "density_max": float(np.max(rho)),
        "pressure_min_Pa": float(np.min(pressure)),
        "pressure_max_Pa": float(np.max(pressure)),
        "mach_min": float(np.min(mach)),
        "mach_max": float(np.max(mach)),
        "finite_state": bool(
            np.all(np.isfinite(rho))
            and np.all(np.isfinite(pressure))
            and np.all(np.isfinite(mach))
        ),
        "positive_density_and_pressure": bool(np.min(rho) > 0.0 and np.min(pressure) > 0.0),
    }
    values.update(strongest_compression(case))
    values.update(symmetry_metrics(case))
    values.update(odd_even_metric(case))
    values.update(pressure_compression_gate_proxy(case))
    return values


def pair_metrics(off: dict[str, object], on: dict[str, object]) -> dict[str, float]:
    if not math.isclose(float(off["time"]), float(on["time"]), rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("paired plotfiles have different final times")

    rho_off = np.asarray(off["rho"])
    rho_on = np.asarray(on["rho"])
    p_off = np.asarray(off["pressure"])
    p_on = np.asarray(on["pressure"])
    m_off = np.asarray(off["mach"])
    m_on = np.asarray(on["mach"])
    velocity_scale = math.sqrt(GAMMA * P_AMBIENT / np.mean(rho_off[:, -4:, -4:]))
    velocity_difference = np.sqrt(
        (np.asarray(on["ux"]) - np.asarray(off["ux"])) ** 2
        + (np.asarray(on["uy"]) - np.asarray(off["uy"])) ** 2
        + (np.asarray(on["uz"]) - np.asarray(off["uz"])) ** 2
    )
    centreline_p_difference = central_profile(
        p_on - p_off, np.asarray(off["y"]), np.asarray(off["z"])
    )
    return {
        "final_time_difference_s": float(float(on["time"]) - float(off["time"])),
        "density_difference_L1_over_mean": float(
            np.mean(np.abs(rho_on - rho_off)) / np.mean(rho_off)
        ),
        "density_difference_Linf_over_mean": float(
            np.max(np.abs(rho_on - rho_off)) / np.mean(rho_off)
        ),
        "pressure_difference_L1_over_pa": float(np.mean(np.abs(p_on - p_off)) / P_AMBIENT),
        "pressure_difference_Linf_over_pa": float(np.max(np.abs(p_on - p_off)) / P_AMBIENT),
        "velocity_difference_L1_over_camb": float(np.mean(velocity_difference) / velocity_scale),
        "velocity_difference_Linf_over_camb": float(np.max(velocity_difference) / velocity_scale),
        "mach_difference_L1": float(np.mean(np.abs(m_on - m_off))),
        "mach_difference_Linf": float(np.max(np.abs(m_on - m_off))),
        "centreline_pressure_difference_Linf_over_pa": float(
            np.max(np.abs(centreline_p_difference)) / P_AMBIENT
        ),
    }


def make_figure(root: Path, cases: dict[str, dict[str, object]]) -> None:
    off = cases["shock_llf_0"]
    on = cases["shock_llf_1"]
    x = np.asarray(off["x"]) / D_EXIT
    y = np.asarray(off["y"]) / D_EXIT
    z = np.asarray(off["z"])
    iz = int(np.argmin(np.abs(z)))
    p_off = np.asarray(off["pressure"])[:, :, iz] / P_AMBIENT
    p_on = np.asarray(on["pressure"])[:, :, iz] / P_AMBIENT
    difference = p_on - p_off
    limit = max(float(np.max(np.abs(difference))), 1.0e-12)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.8), constrained_layout=True)
    common_min = min(float(np.min(p_off)), float(np.min(p_on)))
    common_max = max(float(np.max(p_off)), float(np.max(p_on)))
    for axis, field, title in (
        (axes[0, 0], p_off, r"$\mathrm{afd\_shock\_llf}=0$"),
        (axes[0, 1], p_on, r"$\mathrm{afd\_shock\_llf}=1$"),
    ):
        image = axis.pcolormesh(
            x, y, field.T, shading="auto", cmap="viridis",
            vmin=common_min, vmax=common_max
        )
        axis.set_title(title)
        axis.set_xlabel(r"$x/D_e$")
        axis.set_ylabel(r"$y/D_e$")
        axis.set_aspect("equal")
        fig.colorbar(image, ax=axis, label=r"$p/p_a$")

    difference_image = axes[1, 0].pcolormesh(
        x, y, difference.T, shading="auto", cmap="coolwarm",
        vmin=-limit, vmax=limit
    )
    axes[1, 0].set_title("LLF on minus LLF off")
    axes[1, 0].set_xlabel(r"$x/D_e$")
    axes[1, 0].set_ylabel(r"$y/D_e$")
    axes[1, 0].set_aspect("equal")
    fig.colorbar(difference_image, ax=axes[1, 0], label=r"$\Delta p/p_a$")

    for label, case, style in (
        ("LLF off", off, "-"),
        ("LLF on", on, "--"),
    ):
        pressure_profile = central_profile(
            np.asarray(case["pressure"]), np.asarray(case["y"]), np.asarray(case["z"])
        )
        mach_profile = central_profile(
            np.asarray(case["mach"]), np.asarray(case["y"]), np.asarray(case["z"])
        )
        axes[1, 1].plot(x, pressure_profile / P_AMBIENT, style, label=f"{label}, $p/p_a$")
        axes[1, 1].plot(x, mach_profile, style, alpha=0.65, label=f"{label}, $M$")
    axes[1, 1].set_xlabel(r"$x/D_e$")
    axes[1, 1].set_ylabel(r"$p/p_a$ or $M$")
    axes[1, 1].set_xlim(0.0, 3.0)
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(frameon=False, fontsize=8, ncol=2)

    fig.savefig(root / "jet_sensor_ab_diagnostics.png", dpi=220)
    plt.close(fig)


def write_outputs(
    root: Path,
    rows: list[dict[str, object]],
    pair: dict[str, float],
) -> None:
    (root / "metrics.json").write_text(
        json.dumps({"cases": rows, "pair": pair}, indent=2) + "\n"
    )
    with (root / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (root / "pair_metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value"))
        writer.writerows(pair.items())

    lines = [
        "# Reduced underexpanded round-jet A/B results",
        "",
        "| Metric | LLF off | LLF on |",
        "|---|---:|---:|",
    ]
    selected = (
        "time_s",
        "density_min",
        "pressure_min_Pa",
        "mach_max",
        "compression_x_over_De",
        "compression_halfmax_width_over_De",
        "centreline_mach_upstream",
        "centreline_mach_downstream",
        "nearfield_odd_even_y_max_over_pa",
        "nearfield_odd_even_z_max_over_pa",
        "offline_pressure_compression_gate_shear_fraction",
        "offline_pressure_compression_gate_shock_core_fraction",
    )
    for key in selected:
        lines.append(f"| `{key}` | {rows[0][key]:.8g} | {rows[1][key]:.8g} |")
    lines.extend(("", "## Paired field differences", ""))
    for key, value in pair.items():
        lines.append(f"- `{key}`: {value:.8g}")
    lines.extend(
        (
            "",
            "The compression location is the strongest positive centreline pressure "
            "gradient between 0.4 and 2.5 exit diameters. It is not automatically "
            "classified as a developed Mach disk.",
            "",
            "The pressure-compression values are offline cell-centred proxies. "
            "They do not include the broad smoothness test and are not "
            "instrumented face counts from the solver.",
            "",
        )
    )
    (root / "RESULTS.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    args = parser.parse_args()
    root = args.result_root.resolve()
    cases = {
        label: load_case(root / "runs" / label)
        for label in ("shock_llf_0", "shock_llf_1")
    }
    rows = [case_metrics(label, cases[label]) for label in cases]
    if not all(bool(row["finite_state"]) for row in rows):
        raise SystemExit("non-finite state detected")
    if not all(bool(row["positive_density_and_pressure"]) for row in rows):
        raise SystemExit("non-positive density or pressure detected")
    pair = pair_metrics(cases["shock_llf_0"], cases["shock_llf_1"])
    make_figure(root, cases)
    write_outputs(root, rows, pair)
    print(json.dumps({"cases": rows, "pair": pair}, indent=2))


if __name__ == "__main__":
    main()
