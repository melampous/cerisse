#!/usr/bin/env python3
"""Cell-by-cell matched-time comparison for the Run165 RZ A/B matrix."""

from __future__ import annotations

import csv
import inspect
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import yt
from yt.frontends.boxlib import data_structures as boxlib_data


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_hllc_muscl_vs_afd_weno"
P_INF = 534.5806702256042
GAMMA = 1.4
RHOE_INF = P_INF / (GAMMA - 1.0)
TARGET = (4, 171)

SCHEMES = {
    "muscl": (RESULTS / "early_muscl", RESULTS / "muscl"),
    "afd_off": (
        RESULTS / "early_afd_fallback_off",
        RESULTS / "afd_fallback_off",
    ),
    "afd_on": (
        RESULTS / "early_afd_fallback_on",
        RESULTS / "afd_fallback_on",
    ),
}
LABELS = {
    "muscl": "HLLC-MUSCL O2",
    "afd_off": "AFD-HLLC-WENO-Z5, fallback off",
    "afd_on": "AFD-HLLC-WENO-Z5, fallback on",
}
PAIRS = {
    "afd_off_minus_muscl": "afd_off",
    "afd_on_minus_muscl": "afd_on",
}


def patch_yt_cylindrical_readonly_edge() -> None:
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


def plotfile_for(scheme: str, step: int) -> Path:
    early, main = SCHEMES[scheme]
    directory = early if step <= 20 else main
    result = directory / f"plt{step:05d}"
    if not result.is_dir():
        raise FileNotFoundError(result)
    return result


def load_fields(plotfile: Path):
    dataset = yt.load(str(plotfile))
    level = int(dataset.index.max_level)
    dims = np.asarray(dataset.domain_dimensions, dtype=int).copy()
    dims[:2] *= 2**level
    dims[2] = 1
    grid = dataset.covering_grid(level, dataset.domain_left_edge, dims)

    def field(name: str) -> np.ndarray:
        return np.asarray(
            grid[("boxlib", name)].to_ndarray()[:, :, 0], dtype=np.float64
        )

    rho = field("Density")
    mx = field("Xmom")
    my = field("Ymom")
    mz = field("Zmom")
    energy = field("Energy")
    pressure = field("pressure")
    kinetic = 0.5 * (mx * mx + my * my + mz * mz) / rho
    rhoe = energy - kinetic

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
        "pressure": pressure,
        "rhoe": rhoe,
        "ur": mx / rho,
        "uz": my / rho,
    }, float(dataset.current_time) * 1.0e6, r_edges, z_edges, r_centres, z_centres


def location(array: np.ndarray, r_centres: np.ndarray, z_centres: np.ndarray):
    i, j = np.unravel_index(np.nanargmax(array), array.shape)
    return int(i), int(j), float(r_centres[i]), float(z_centres[j])


def minimum_location(array: np.ndarray, r_centres: np.ndarray, z_centres: np.ndarray):
    i, j = np.unravel_index(np.nanargmin(array), array.shape)
    return int(i), int(j), float(r_centres[i]), float(z_centres[j])


def bbox(mask: np.ndarray, r_centres: np.ndarray, z_centres: np.ndarray):
    indices = np.argwhere(mask)
    if indices.size == 0:
        return (0, np.nan, np.nan, np.nan, np.nan)
    imin, jmin = indices.min(axis=0)
    imax, jmax = indices.max(axis=0)
    return (
        int(indices.shape[0]),
        float(r_centres[imin]),
        float(r_centres[imax]),
        float(z_centres[jmin]),
        float(z_centres[jmax]),
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    patch_yt_cylindrical_readonly_edge()
    yt.funcs.mylog.setLevel(50)

    steps = list(range(0, 21)) + [25, 50, 75, 100, 125] + list(range(150, 185))
    thresholds = [1.0e-12, 1.0e-6, 1.0e-3, 1.0e-2, 1.0e-1]
    state_rows: list[dict] = []
    difference_rows: list[dict] = []
    target_rows: list[dict] = []
    first_crossing: dict[tuple[str, float], dict] = {}
    pointwise: dict[str, list[np.ndarray]] = {
        f"{pair}_{field}": []
        for pair in PAIRS
        for field in ("pressure", "rhoe", "rho", "uz")
    }
    times: list[float] = []
    cached_for_maps: dict[int, dict[str, dict[str, np.ndarray]]] = {}
    map_steps = {1, 25, 150, 175, 184}

    r_edges = z_edges = r_centres = z_centres = None
    for step in steps:
        states: dict[str, dict[str, np.ndarray]] = {}
        scheme_times = []
        for scheme in SCHEMES:
            fields, time_us, re, ze, rc, zc = load_fields(
                plotfile_for(scheme, step)
            )
            states[scheme] = fields
            scheme_times.append(time_us)
            if r_edges is None:
                r_edges, z_edges, r_centres, z_centres = re, ze, rc, zc
            minimum = minimum_location(fields["rhoe"], rc, zc)
            ti, tj = TARGET
            state_rows.append(
                {
                    "step": step,
                    "time_us": f"{time_us:.12g}",
                    "scheme": scheme,
                    "min_rhoe_J_m3": f"{fields['rhoe'].min():.17g}",
                    "min_rhoe_i": minimum[0],
                    "min_rhoe_j": minimum[1],
                    "min_rhoe_r_mm": f"{minimum[2]:.12g}",
                    "min_rhoe_x_mm": f"{minimum[3]:.12g}",
                    "min_rho": f"{fields['rho'].min():.17g}",
                    "max_abs_uz_m_s": f"{np.abs(fields['uz']).max():.17g}",
                }
            )
            target_rows.append(
                {
                    "step": step,
                    "time_us": f"{time_us:.12g}",
                    "scheme": scheme,
                    "i": ti,
                    "j": tj,
                    "r_mm": f"{rc[ti]:.12g}",
                    "x_mm": f"{zc[tj]:.12g}",
                    "rho": f"{fields['rho'][ti,tj]:.17g}",
                    "mx": f"{fields['mx'][ti,tj]:.17g}",
                    "my": f"{fields['my'][ti,tj]:.17g}",
                    "energy_J_m3": f"{fields['energy'][ti,tj]:.17g}",
                    "rhoe_J_m3": f"{fields['rhoe'][ti,tj]:.17g}",
                    "pressure_Pa": f"{fields['pressure'][ti,tj]:.17g}",
                    "ur_m_s": f"{fields['ur'][ti,tj]:.17g}",
                    "uz_m_s": f"{fields['uz'][ti,tj]:.17g}",
                }
            )
        if max(scheme_times) - min(scheme_times) > 1.0e-9:
            raise RuntimeError(f"non-synchronous step {step}: {scheme_times}")
        time_us = scheme_times[0]
        times.append(time_us)

        if step in map_steps:
            cached_for_maps[step] = states

        base = states["muscl"]
        for pair, scheme in PAIRS.items():
            candidate = states[scheme]
            differences = {
                name: candidate[name] - base[name]
                for name in ("pressure", "rhoe", "rho", "uz")
            }
            for name, values in differences.items():
                pointwise[f"{pair}_{name}"].append(values.copy())

            dp = differences["pressure"]
            de = differences["rhoe"]
            dp_abs = np.abs(dp)
            de_abs = np.abs(de)
            p_argmax = location(dp_abs, r_centres, z_centres)
            e_argmax = location(de_abs, r_centres, z_centres)
            row = {
                "step": step,
                "time_us": f"{time_us:.12g}",
                "pair": pair,
                "pressure_l1_over_pinf": f"{np.mean(dp_abs)/P_INF:.17g}",
                "pressure_l2_over_pinf": f"{np.sqrt(np.mean(dp*dp))/P_INF:.17g}",
                "pressure_linf_over_pinf": f"{dp_abs.max()/P_INF:.17g}",
                "pressure_argmax_i": p_argmax[0],
                "pressure_argmax_j": p_argmax[1],
                "pressure_argmax_r_mm": f"{p_argmax[2]:.12g}",
                "pressure_argmax_x_mm": f"{p_argmax[3]:.12g}",
                "rhoe_l1_over_freestream": f"{np.mean(de_abs)/RHOE_INF:.17g}",
                "rhoe_l2_over_freestream": f"{np.sqrt(np.mean(de*de))/RHOE_INF:.17g}",
                "rhoe_linf_over_freestream": f"{de_abs.max()/RHOE_INF:.17g}",
                "rhoe_argmax_i": e_argmax[0],
                "rhoe_argmax_j": e_argmax[1],
                "rhoe_argmax_r_mm": f"{e_argmax[2]:.12g}",
                "rhoe_argmax_x_mm": f"{e_argmax[3]:.12g}",
            }
            normalized_de = de_abs / RHOE_INF
            for threshold in thresholds:
                count, rmin, rmax, zmin, zmax = bbox(
                    normalized_de > threshold, r_centres, z_centres
                )
                tag = f"{threshold:.0e}"
                row[f"cells_rhoe_gt_{tag}"] = count
                row[f"bbox_rmin_mm_{tag}"] = f"{rmin:.12g}"
                row[f"bbox_rmax_mm_{tag}"] = f"{rmax:.12g}"
                row[f"bbox_xmin_mm_{tag}"] = f"{zmin:.12g}"
                row[f"bbox_xmax_mm_{tag}"] = f"{zmax:.12g}"
                key = (pair, threshold)
                if count and key not in first_crossing:
                    first_crossing[key] = {
                        "pair": pair,
                        "threshold": threshold,
                        "step": step,
                        "time_us": time_us,
                        "cells": count,
                        "rmin": rmin,
                        "rmax": rmax,
                        "xmin": zmin,
                        "xmax": zmax,
                    }
            difference_rows.append(row)

    write_csv(RESULTS / "state_minima_timeline.csv", state_rows)
    write_csv(RESULTS / "point_target_4_171_timeline.csv", target_rows)
    write_csv(RESULTS / "field_difference_timeline.csv", difference_rows)

    npz_values = {
        "steps": np.asarray(steps, dtype=np.int32),
        "times_us": np.asarray(times),
        "r_centres_mm": r_centres,
        "x_centres_mm": z_centres,
    }
    npz_values.update(
        {
            name: np.stack(values, axis=0)
            for name, values in pointwise.items()
        }
    )
    np.savez_compressed(RESULTS / "pointwise_difference_fields.npz", **npz_values)

    with (RESULTS / "first_difference_thresholds.txt").open("w") as stream:
        for pair in PAIRS:
            for threshold in thresholds:
                item = first_crossing.get((pair, threshold))
                if item is None:
                    stream.write(
                        f"{pair} |drhoe|/rhoe_inf>{threshold:g}: never\n"
                    )
                else:
                    stream.write(
                        f"{pair} |drhoe|/rhoe_inf>{threshold:g}: "
                        f"step={item['step']} time_us={item['time_us']:.6g} "
                        f"cells={item['cells']} r_mm=[{item['rmin']:.6g},"
                        f"{item['rmax']:.6g}] x_mm=[{item['xmin']:.6g},"
                        f"{item['xmax']:.6g}]\n"
                    )

    # Timeline figure.
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), constrained_layout=True)
    diff_lookup = {(int(row["step"]), row["pair"]): row for row in difference_rows}
    state_lookup = {(int(row["step"]), row["scheme"]): row for row in state_rows}
    target_lookup = {(int(row["step"]), row["scheme"]): row for row in target_rows}
    for pair, color in zip(PAIRS, ("tab:red", "tab:blue")):
        axes[0, 0].semilogy(
            times,
            [float(diff_lookup[(s, pair)]["pressure_l2_over_pinf"]) for s in steps],
            color=color,
            label=pair.replace("_", " "),
        )
        axes[0, 1].semilogy(
            times,
            [float(diff_lookup[(s, pair)]["rhoe_linf_over_freestream"]) for s in steps],
            color=color,
            label=pair.replace("_", " "),
        )
    for scheme, color in zip(SCHEMES, ("black", "tab:red", "tab:blue")):
        axes[1, 0].semilogy(
            times,
            [float(state_lookup[(s, scheme)]["min_rhoe_J_m3"]) for s in steps],
            color=color,
            label=LABELS[scheme],
        )
        axes[1, 1].semilogy(
            times,
            [float(target_lookup[(s, scheme)]["rhoe_J_m3"]) for s in steps],
            color=color,
            label=LABELS[scheme],
        )
    axes[0, 0].set(ylabel=r"$L_2(\Delta p)/p_\infty$", title="Full-field pressure difference")
    axes[0, 1].set(ylabel=r"$L_\infty(\Delta \rho e)/(\rho e)_\infty$", title="Worst cell-wise internal-energy difference")
    axes[1, 0].set(xlabel="time [microseconds]", ylabel=r"global min $\rho e$ [J/m$^3$]", title="Admissibility margin")
    axes[1, 1].set(xlabel="time [microseconds]", ylabel=r"target $\rho e$ [J/m$^3$]", title="Cell (4,171), pre-collapse")
    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.25)
        axis.axvline(73.8, color="0.4", linestyle="--", linewidth=0.8)
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].legend(fontsize=8)
    fig.savefig(RESULTS / "difference_growth_timeline.png", dpi=220)
    plt.close(fig)

    # Spatial maps of the point-wise internal-energy difference.
    fig, axes = plt.subplots(2, len(map_steps), figsize=(17.0, 6.3), constrained_layout=True)
    for column, step in enumerate(sorted(map_steps)):
        states = cached_for_maps[step]
        time_us = step * 0.4
        for row, scheme in enumerate(("afd_off", "afd_on")):
            normalized = np.abs(states[scheme]["rhoe"] - states["muscl"]["rhoe"]) / RHOE_INF
            image = axes[row, column].pcolormesh(
                z_edges,
                r_edges,
                np.log10(normalized + 1.0e-16),
                shading="auto",
                cmap="magma",
                vmin=-8,
                vmax=4,
            )
            axes[row, column].scatter(
                [z_centres[TARGET[1]]], [r_centres[TARGET[0]]],
                marker="x", color="cyan", s=25, linewidths=1.0,
            )
            axes[row, column].set(xlim=(-45, 2.5), ylim=(0, 22))
            axes[row, column].set_title(f"{LABELS[scheme]}\n{time_us:.1f} us", fontsize=8)
            if column == 0:
                axes[row, column].set_ylabel("r [mm]")
            if row == 1:
                axes[row, column].set_xlabel("x [mm]")
    fig.colorbar(image, ax=axes, label=r"$\log_{10}(|\Delta \rho e|/(\rho e)_\infty)$", shrink=0.82)
    fig.savefig(RESULTS / "pointwise_internal_energy_difference_maps.png", dpi=220)
    plt.close(fig)

    # Last accepted-state neighborhood before the failed 73.8-us RK stage.
    states = cached_for_maps[184]
    radius_i = 8
    radius_j = 10
    i0, j0 = TARGET
    islice = slice(max(0, i0 - radius_i), min(len(r_centres), i0 + radius_i + 1))
    jslice = slice(max(0, j0 - radius_j), min(len(z_centres), j0 + radius_j + 1))
    panels = [
        ("muscl", states["muscl"]["rhoe"] / RHOE_INF, "HLLC-MUSCL"),
        ("afd_off", states["afd_off"]["rhoe"] / RHOE_INF, "AFD fallback off"),
        ("afd_on", states["afd_on"]["rhoe"] / RHOE_INF, "AFD fallback on"),
        ("difference", (states["afd_off"]["rhoe"] - states["muscl"]["rhoe"]) / RHOE_INF, "AFD off - MUSCL"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.8), constrained_layout=True)
    for axis, (_, values, title) in zip(axes, panels):
        local = values[islice, jslice]
        limit = np.max(np.abs(local))
        if title == "AFD off - MUSCL":
            shown = axis.imshow(
                local,
                origin="lower",
                aspect="auto",
                cmap="coolwarm",
                vmin=-limit,
                vmax=limit,
                extent=(z_edges[jslice.start], z_edges[jslice.stop], r_edges[islice.start], r_edges[islice.stop]),
            )
        else:
            shown = axis.imshow(
                local,
                origin="lower",
                aspect="auto",
                cmap="viridis",
                extent=(z_edges[jslice.start], z_edges[jslice.stop], r_edges[islice.start], r_edges[islice.stop]),
            )
        axis.scatter([z_centres[j0]], [r_centres[i0]], marker="x", color="red", s=35)
        axis.set(title=title, xlabel="x [mm]")
        fig.colorbar(shown, ax=axis, shrink=0.8)
    axes[0].set_ylabel("r [mm]")
    fig.suptitle(r"Last accepted coarse state, t=73.6 us: $\rho e/(\rho e)_\infty$")
    fig.savefig(RESULTS / "precollapse_target_neighborhood.png", dpi=220)
    plt.close(fig)

    print(RESULTS / "field_difference_timeline.csv")
    print(RESULTS / "pointwise_difference_fields.npz")
    print(RESULTS / "difference_growth_timeline.png")
    print(RESULTS / "pointwise_internal_energy_difference_maps.png")
    print(RESULTS / "precollapse_target_neighborhood.png")


if __name__ == "__main__":
    main()
