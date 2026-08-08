#!/usr/bin/env python3
"""Build audited Run-165 RZ side-view geometries with and without N1.

``with_nozzle`` is the public-data reconstruction of the Run-165 centre N1
divergent passage and throat face.  ``without_nozzle`` is a deliberately
derived numerical control: the nominal R25.4 sphere is continued to its axis
tip, so there is no opening, internal passage, injection patch, or jet.

Both Cerisse geometry files are a single CCW polygon in ``(r_m, x_TSP_m)``
order with implicit last-to-first closure.  This script changes case geometry
and visualization only; it does not implement cut-cell, embedded-boundary,
aggregate, or cut-control flow updates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from dataclasses import asdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853-rz-geometry-ab")

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "test1853-rz-nozzle-geometry-ab-v1"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from build_run165_axisymmetric_section import (
    BUILD_AUDIT,
    PARAMETER_FILE,
    Geometry,
    assemble_profile,
    build_segments,
    derive_geometry,
    exact_volume_mm3,
    load_json,
    polygon_area_centroid_r,
    polygon_is_simple,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "outputs" / "nozzle_geometry_ab"
REFERENCE_WITH_NOZZLE = (
    HERE
    / "outputs"
    / "jet_plane_variants"
    / "run165_srp_throat_choked_profile.dat"
)

SOURCE_CLASS = {
    "axis_closure": "RZ_NUMERICAL_CLOSURE",
    "throat_jet_patch": "NASA_AS_BUILT_PLUS_DERIVED_X",
    "nozzle_divergent_wall": "NASA_AS_BUILT_PLUS_PLACEMENT_ASSUMPTION",
    "closed_sphere_oml": "NOMINAL_PUBLIC_DATA_RECONSTRUCTION_CONTROL",
    "sphere_oml": "NOMINAL_PUBLIC_DATA_RECONSTRUCTION",
    "cone_oml": "NASA_NOMINAL_PLUS_DERIVED_TANGENCY",
    "shoulder_oml": "NASA_RADIUS_PLUS_DERIVED_TANGENCY",
    "cylinder_oml": "NASA_DIMENSIONED",
    "aft_cap": "ARTIFICIAL_IBM_CLOSURE",
}

COLORS = {
    "axis_closure": "#616161",
    "throat_jet_patch": "#f28e2b",
    "nozzle_divergent_wall": "#d62728",
    "closed_sphere_oml": "#0072b2",
    "sphere_oml": "#1769aa",
    "cone_oml": "#1769aa",
    "shoulder_oml": "#2ca02c",
    "cylinder_oml": "#1769aa",
    "aft_cap": "#9467bd",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_without_nozzle_segments(g: Geometry) -> dict[str, np.ndarray]:
    """Return the sealed nominal sphere-cone in CCW ``(r_mm,x_mm)`` order."""
    x_sphere = np.linspace(
        g.nominal_nose_tip_x_mm, g.sphere_cone_x_mm, 180
    )
    r_sphere = np.sqrt(
        np.maximum(
            0.0,
            g.nose_radius_mm**2
            - (x_sphere - g.sphere_center_x_mm) ** 2,
        )
    )
    x_shoulder = np.linspace(
        g.cone_shoulder_x_mm, g.shoulder_center_x_mm, 150
    )
    r_shoulder = g.shoulder_center_r_mm + np.sqrt(
        np.maximum(
            0.0,
            g.shoulder_radius_mm**2
            - (x_shoulder - g.shoulder_center_x_mm) ** 2,
        )
    )
    return {
        "axis_closure": np.array(
            [
                [0.0, g.aft_x_mm],
                [0.0, g.nominal_nose_tip_x_mm],
            ],
            dtype=float,
        ),
        "closed_sphere_oml": np.column_stack((r_sphere, x_sphere)),
        "cone_oml": np.array(
            [
                [g.sphere_cone_r_mm, g.sphere_cone_x_mm],
                [g.cone_shoulder_r_mm, g.cone_shoulder_x_mm],
            ],
            dtype=float,
        ),
        "shoulder_oml": np.column_stack((r_shoulder, x_shoulder)),
        "cylinder_oml": np.array(
            [
                [g.body_radius_mm, g.shoulder_center_x_mm],
                [g.body_radius_mm, g.aft_x_mm],
            ],
            dtype=float,
        ),
        "aft_cap": np.array(
            [
                [g.body_radius_mm, g.aft_x_mm],
                [0.0, g.aft_x_mm],
            ],
            dtype=float,
        ),
    }


def assert_segment_chain(segments: dict[str, np.ndarray]) -> None:
    items = list(segments.items())
    for index, (name, values) in enumerate(items):
        if len(values) < 2 or not np.all(np.isfinite(values)):
            raise RuntimeError(f"invalid segment {name}")
        next_name, next_values = items[(index + 1) % len(items)]
        error = float(np.linalg.norm(values[-1] - next_values[0]))
        if error > 1.0e-11:
            raise RuntimeError(
                f"open segment chain {name}->{next_name}: {error} mm"
            )


def validate_profile(
    profile: np.ndarray, segments: dict[str, np.ndarray]
) -> dict[str, float | int | bool]:
    assert_segment_chain(segments)
    area_mm2, centroid_r_mm = polygon_area_centroid_r(profile)
    edges = np.roll(profile, -1, axis=0) - profile
    edge_lengths = np.linalg.norm(edges, axis=1)
    if area_mm2 <= 0.0:
        raise RuntimeError(f"profile is not CCW: area={area_mm2}")
    if not polygon_is_simple(profile):
        raise RuntimeError("profile is self-intersecting")
    if np.any(~np.isfinite(profile)) or float(np.min(profile[:, 0])) < 0.0:
        raise RuntimeError("profile has non-finite or negative-r coordinates")
    if float(np.min(edge_lengths)) <= 1.0e-14:
        raise RuntimeError("profile has a zero-length edge")
    return {
        "vertices": int(len(profile)),
        "ccw": True,
        "simple_polygon": True,
        "signed_area_mm2": float(area_mm2),
        "centroid_r_mm": float(centroid_r_mm),
        "polyline_revolution_volume_mm3": float(
            2.0 * math.pi * area_mm2 * centroid_r_mm
        ),
        "minimum_edge_mm": float(np.min(edge_lengths)),
        "maximum_edge_mm": float(np.max(edge_lengths)),
        "minimum_r_mm": float(np.min(profile[:, 0])),
        "maximum_r_mm": float(np.max(profile[:, 0])),
        "minimum_x_mm": float(np.min(profile[:, 1])),
        "maximum_x_mm": float(np.max(profile[:, 1])),
    }


def sphere_r2_integral(g: Geometry, x0: float, x1: float) -> float:
    t0 = x0 - g.sphere_center_x_mm
    t1 = x1 - g.sphere_center_x_mm
    return g.nose_radius_mm**2 * (x1 - x0) - (t1**3 - t0**3) / 3.0


def exact_without_nozzle_volume_mm3(g: Geometry) -> float:
    beta = math.radians(g.nozzle_divergent_half_angle_deg)
    divergent_cavity = (
        math.pi
        * (g.virtual_exit_r_mm**3 - g.throat_r_mm**3)
        / (3.0 * math.tan(beta))
    )
    closed_nose_cap = math.pi * sphere_r2_integral(
        g, g.nominal_nose_tip_x_mm, g.virtual_exit_x_mm
    )
    return exact_volume_mm3(g) + divergent_cavity + closed_nose_cap


def write_profile(
    output: Path,
    stem: str,
    profile: np.ndarray,
    outgoing: list[str],
) -> tuple[Path, Path]:
    dat = output / f"{stem}.dat"
    csv_path = output / f"{stem}.csv"
    np.savetxt(dat, profile * 1.0e-3, fmt="%.15e")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "point_id",
                "outgoing_segment",
                "source_class",
                "r_m",
                "x_m",
                "r_mm",
                "x_mm",
            ]
        )
        for index, (point, segment) in enumerate(zip(profile, outgoing)):
            writer.writerow(
                [
                    index,
                    segment,
                    SOURCE_CLASS[segment],
                    f"{point[0] * 1.0e-3:.15e}",
                    f"{point[1] * 1.0e-3:.15e}",
                    f"{point[0]:.12f}",
                    f"{point[1]:.12f}",
                ]
            )
    return dat, csv_path


def draw_profile(
    axis: plt.Axes,
    profile: np.ndarray,
    segments: dict[str, np.ndarray],
    *,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    title: str,
) -> None:
    x = profile[:, 1]
    r = profile[:, 0]
    axis.fill(x, r, color="#dfe5eb", zorder=0)
    axis.fill(x, -r, color="#dfe5eb", zorder=0)
    for name, values in segments.items():
        sx = values[:, 1]
        sr = values[:, 0]
        style = "--" if name == "axis_closure" else "-"
        width = 3.0 if name == "throat_jet_patch" else 2.0
        axis.plot(sx, sr, style, color=COLORS[name], lw=width, zorder=3)
        axis.plot(sx, -sr, style, color=COLORS[name], lw=width, zorder=3)
    axis.axhline(0.0, color="0.45", lw=0.8)
    axis.axvline(0.0, color="0.45", lw=0.8, ls=":")
    axis.set_xlim(*xlim)
    axis.set_ylim(*ylim)
    axis.set_aspect("equal", adjustable="box")
    axis.grid(alpha=0.18)
    axis.set_xlabel("x from theoretical sharp point [mm]")
    axis.set_ylabel("signed radius [mm]")
    axis.set_title(title)


def render_individual(
    output: Path,
    stem: str,
    title: str,
    subtitle: str,
    g: Geometry,
    profile: np.ndarray,
    segments: dict[str, np.ndarray],
) -> tuple[Path, Path]:
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, 6.5),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [2.1, 1.0]},
    )
    draw_profile(
        axes[0],
        profile,
        segments,
        xlim=(-8.0, g.aft_x_mm + 7.0),
        ylim=(-70.0, 70.0),
        title="Complete meridional side view",
    )
    draw_profile(
        axes[1],
        profile,
        segments,
        xlim=(-3.0, 31.0),
        ylim=(-25.0, 25.0),
        title="Forebody/nozzle detail",
    )
    if "throat_jet_patch" in segments:
        axes[1].annotate(
            "N1 throat",
            (g.throat_x_mm, 0.55 * g.throat_r_mm),
            xytext=(17.0, 7.0),
            arrowprops={"arrowstyle": "->", "lw": 0.9},
        )
        axes[1].annotate(
            "N1 physical lip",
            (g.virtual_exit_x_mm, g.virtual_exit_r_mm),
            xytext=(7.0, 13.0),
            arrowprops={"arrowstyle": "->", "lw": 0.9},
        )
    else:
        axes[1].annotate(
            "sealed nominal sphere tip\n(no opening / no jet patch)",
            (g.nominal_nose_tip_x_mm, 0.0),
            xytext=(8.0, 10.0),
            arrowprops={"arrowstyle": "->", "lw": 0.9},
        )
    handles = [
        Patch(facecolor="#dfe5eb", label="IBM solid"),
        Line2D([0], [0], color="#1769aa", lw=2, label="physical OML"),
        Line2D([0], [0], color="#2ca02c", lw=2, label="R2.54 shoulder"),
        Line2D([0], [0], color="#9467bd", lw=2, label="artificial aft cap"),
        Line2D([0], [0], color="#616161", lw=2, ls="--", label="RZ axis closure"),
    ]
    if "throat_jet_patch" in segments:
        handles.extend(
            [
                Line2D([0], [0], color="#d62728", lw=2, label="N1 divergent wall"),
                Line2D([0], [0], color="#f28e2b", lw=3, label="N1 throat patch"),
            ]
        )
    axes[0].legend(handles=handles, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle(f"{title}\n{subtitle}", fontsize=14)
    png = output / f"{stem}.png"
    svg = output / f"{stem}.svg"
    fig.savefig(png, dpi=240)
    fig.savefig(svg, metadata={"Date": None})
    plt.close(fig)
    return png, svg


def render_comparison(
    output: Path,
    g: Geometry,
    with_profile: np.ndarray,
    with_segments: dict[str, np.ndarray],
    without_profile: np.ndarray,
    without_segments: dict[str, np.ndarray],
) -> tuple[Path, Path]:
    fig, axes = plt.subplots(2, 2, figsize=(17, 9), constrained_layout=True)
    rows = (
        (with_profile, with_segments, "A — with resolved centre N1 nozzle"),
        (
            without_profile,
            without_segments,
            "B — sealed no-nozzle numerical control",
        ),
    )
    for row, (profile, segments, title) in enumerate(rows):
        draw_profile(
            axes[row, 0],
            profile,
            segments,
            xlim=(-8.0, g.aft_x_mm + 7.0),
            ylim=(-70.0, 70.0),
            title=title,
        )
        draw_profile(
            axes[row, 1],
            profile,
            segments,
            xlim=(-3.0, 31.0),
            ylim=(-25.0, 25.0),
            title="same-scale forebody detail",
        )
    fig.suptitle(
        "NASA Test 1853 axisymmetric geometry A/B — actual Cerisse polygons\n"
        "same TSP origin, shoulder, cylinder and aft closure; alpha=0 only",
        fontsize=14,
    )
    png = output / "run165_rz_nozzle_geometry_side_view_comparison.png"
    svg = output / "run165_rz_nozzle_geometry_side_view_comparison.svg"
    fig.savefig(png, dpi=240)
    fig.savefig(svg, metadata={"Date": None})
    plt.close(fig)
    return png, svg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    g = derive_geometry(load_json(PARAMETER_FILE))
    with_segments = build_segments(g)
    without_segments = build_without_nozzle_segments(g)
    with_profile, with_outgoing = assemble_profile(with_segments)
    without_profile, without_outgoing = assemble_profile(without_segments)
    with_qa = validate_profile(with_profile, with_segments)
    without_qa = validate_profile(without_profile, without_segments)

    reference = np.loadtxt(REFERENCE_WITH_NOZZLE) * 1.0e3
    if reference.shape != with_profile.shape:
        raise RuntimeError("with-nozzle profile shape differs from audited Case T")
    reference_error_mm = float(np.max(np.abs(reference - with_profile)))
    if reference_error_mm > 1.0e-11:
        raise RuntimeError(
            f"with-nozzle profile differs from audited Case T by {reference_error_mm} mm"
        )

    with_exact = exact_volume_mm3(g)
    without_exact = exact_without_nozzle_volume_mm3(g)
    with_qa["analytic_revolution_volume_mm3"] = with_exact
    without_qa["analytic_revolution_volume_mm3"] = without_exact
    for qa, exact in ((with_qa, with_exact), (without_qa, without_exact)):
        qa["polyline_volume_relative_error"] = abs(
            float(qa["polyline_revolution_volume_mm3"]) - exact
        ) / exact
        if float(qa["polyline_volume_relative_error"]) > 2.0e-7:
            raise RuntimeError(f"polyline volume QA failed: {qa}")

    with_dat, with_csv = write_profile(
        output, "run165_rz_with_nozzle_profile", with_profile, with_outgoing
    )
    without_dat, without_csv = write_profile(
        output,
        "run165_rz_without_nozzle_profile",
        without_profile,
        without_outgoing,
    )
    if sha256(with_dat) != sha256(REFERENCE_WITH_NOZZLE):
        raise RuntimeError("with-nozzle DAT is not byte-identical to audited Case T")

    with_png, with_svg = render_individual(
        output,
        "run165_rz_with_nozzle_side_view",
        "With nozzle: reconstructed Run-165 centre N1",
        "resolved divergent passage; throat disk is the choked injection patch",
        g,
        with_profile,
        with_segments,
    )
    without_png, without_svg = render_individual(
        output,
        "run165_rz_without_nozzle_side_view",
        "Without nozzle: sealed sphere-cone numerical control",
        "not a NASA Run-165 hardware claim; no opening, passage, or jet patch",
        g,
        without_profile,
        without_segments,
    )
    comparison_png, comparison_svg = render_comparison(
        output,
        g,
        with_profile,
        with_segments,
        without_profile,
        without_segments,
    )

    generated = [
        with_dat,
        with_csv,
        without_dat,
        without_csv,
        with_png,
        with_svg,
        without_png,
        without_svg,
        comparison_png,
        comparison_svg,
    ]
    metadata = {
        "status": "GEOMETRY_AB_QA_PASS",
        "purpose": "Run165 RZ resolved-N1 versus sealed no-nozzle geometry control",
        "method_scope": {
            "family": "pure_shared_gp_full_cartesian",
            "impact": "case geometry and visualization only",
            "production_rhs_changed": False,
            "cut_cell_or_eb_enabled": False,
            "aggregate_or_cut_control_enabled": False,
        },
        "coordinate_convention": {
            "dat_columns": ["r_m", "x_from_TSP_m"],
            "axis": "r=0",
            "winding": "CCW",
            "closure": "implicit last vertex to first vertex",
            "freestream": "+x",
            "retropropulsive_jet_when_present": "-x",
        },
        "geometry_parameters": asdict(g),
        "with_nozzle": {
            "hardware_semantics": "Run165 public-data reconstruction",
            "jet_semantics": "M=1 choked state at N1 throat when enabled by the case",
            "audited_case_T_coordinate_max_error_mm": reference_error_mm,
            "qa": with_qa,
        },
        "without_nozzle": {
            "hardware_semantics": "derived sealed numerical control, not a NASA run claim",
            "jet_semantics": "no injection patch and jet must be disabled",
            "sealed_tip_x_mm": g.nominal_nose_tip_x_mm,
            "qa": without_qa,
        },
        "semantic_guards": [
            "The no-nozzle geometry is not Case E virtual-exit injection.",
            "Jet-OFF with the N1 geometry retains the nozzle cavity and is not the sealed control.",
            "The with-nozzle geometry resolves only the divergent passage; no convergent passage or plenum is invented.",
            "Artificial axis and aft-cap closures are excluded from physical loads.",
            "RZ force integration requires 2*pi*r weighting.",
        ],
        "source_sha256": {
            "geometry_parameters.yaml": sha256(PARAMETER_FILE),
            "build_audit.json": sha256(BUILD_AUDIT),
            "audited_case_T.dat": sha256(REFERENCE_WITH_NOZZLE),
            "generator": sha256(Path(__file__).resolve()),
        },
        "outputs_sha256": {path.name: sha256(path) for path in generated},
    }
    metadata_path = output / "run165_rz_nozzle_geometry_ab_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {output}")
    print(
        "[RUN165-RZ-GEOMETRY-AB-PASS] "
        f"with_vertices={len(with_profile)} "
        f"without_vertices={len(without_profile)} "
        f"with_sha={sha256(with_dat)} "
        f"without_sha={sha256(without_dat)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
