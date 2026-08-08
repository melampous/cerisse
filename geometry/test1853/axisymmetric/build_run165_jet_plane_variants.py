#!/usr/bin/env python3
"""Build and draw the two Run-165 RZ jet-plane representations.

Case T resolves the reconstructed N1 divergent passage and applies a choked
state on the throat disk.  Case E removes that internal passage from the fluid
domain and applies the corresponding supersonic state directly on the virtual
exit disk.  Both profiles use the same reconstructed external mold line.

The Cerisse polygon files are written in ``(r_m, x_TSP_m)`` order, counter-
clockwise, with implicit last-to-first closure.  This generator changes only
case geometry and visualization; it does not implement or enable any cut-cell,
embedded-boundary, aggregate, or cut-control-volume flow method.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from dataclasses import asdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853-jet-planes")

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "nasa-test1853-run165-jet-plane-v1"
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
DEFAULT_OUTPUT = HERE / "outputs" / "jet_plane_variants"
REFERENCE_THROAT_DAT = (
    HERE / "outputs" / "run165_axisymmetric_section_preview.dat"
)

SOURCE_CLASS = {
    "axis_closure": "RZ_NUMERICAL_CLOSURE",
    "throat_jet_patch": "NASA_AS_BUILT_PLUS_DERIVED_X",
    "virtual_exit_jet_patch": (
        "NASA_AS_BUILT_EXIT_AREA_PLUS_PLACEMENT_ASSUMPTION"
    ),
    "nozzle_divergent_wall": (
        "NASA_AS_BUILT_PLUS_PLACEMENT_ASSUMPTION"
    ),
    "sphere_oml": "NOMINAL_PUBLIC_DATA_RECONSTRUCTION",
    "cone_oml": "NASA_NOMINAL_PLUS_DERIVED_TANGENCY",
    "shoulder_oml": "NASA_RADIUS_PLUS_DERIVED_TANGENCY",
    "cylinder_oml": "NASA_DIMENSIONED",
    "aft_cap": "ARTIFICIAL_IBM_CLOSURE",
}

COLORS = {
    "axis_closure": "#616161",
    "throat_jet_patch": "#f28e2b",
    "virtual_exit_jet_patch": "#00a6a6",
    "nozzle_divergent_wall": "#d62728",
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


def build_virtual_exit_segments(g: Geometry) -> dict[str, np.ndarray]:
    """Return Case-E segments in CCW ``(r_mm, x_mm)`` traversal order."""
    x_sphere = np.linspace(
        g.virtual_exit_x_mm, g.sphere_cone_x_mm, 150
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
                [0.0, g.virtual_exit_x_mm],
            ],
            dtype=float,
        ),
        "virtual_exit_jet_patch": np.array(
            [
                [0.0, g.virtual_exit_x_mm],
                [g.virtual_exit_r_mm, g.virtual_exit_x_mm],
            ],
            dtype=float,
        ),
        "sphere_oml": np.column_stack((r_sphere, x_sphere)),
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
        if len(values) < 2:
            raise RuntimeError(f"{name} has fewer than two points")
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"{name} contains a non-finite coordinate")
        next_name, next_values = items[(index + 1) % len(items)]
        error = float(np.linalg.norm(values[-1] - next_values[0]))
        if error > 1.0e-11:
            raise RuntimeError(
                f"segment chain is open: {name}->{next_name}, "
                f"endpoint error={error} mm"
            )


def validate_profile(
    profile: np.ndarray, segments: dict[str, np.ndarray]
) -> dict[str, float | int | bool]:
    assert_segment_chain(segments)
    area_mm2, centroid_r_mm = polygon_area_centroid_r(profile)
    if area_mm2 <= 0.0:
        raise RuntimeError(f"profile must be CCW, signed area={area_mm2}")
    if not polygon_is_simple(profile):
        raise RuntimeError("profile is self-intersecting")
    if float(np.min(profile[:, 0])) < 0.0:
        raise RuntimeError("profile contains negative radius")
    edge_lengths = np.linalg.norm(
        np.roll(profile, -1, axis=0) - profile, axis=1
    )
    if float(np.min(edge_lengths)) <= 1.0e-14:
        raise RuntimeError("profile contains a zero-length edge")
    return {
        "vertices": int(len(profile)),
        "signed_area_mm2": float(area_mm2),
        "centroid_r_mm": float(centroid_r_mm),
        "polyline_revolution_volume_mm3": float(
            2.0 * math.pi * area_mm2 * centroid_r_mm
        ),
        "minimum_edge_mm": float(np.min(edge_lengths)),
        "minimum_r_mm": float(np.min(profile[:, 0])),
        "maximum_r_mm": float(np.max(profile[:, 0])),
        "minimum_x_mm": float(np.min(profile[:, 1])),
        "maximum_x_mm": float(np.max(profile[:, 1])),
        "ccw": True,
        "simple_polygon": True,
    }


def write_profile(
    output: Path,
    stem: str,
    profile: np.ndarray,
    outgoing: list[str],
) -> tuple[Path, Path]:
    dat_path = output / f"{stem}.dat"
    csv_path = output / f"{stem}.csv"
    np.savetxt(dat_path, profile * 1.0e-3, fmt="%.15e")
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
        for index, (point, segment) in enumerate(
            zip(profile, outgoing)
        ):
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
    return dat_path, csv_path


def write_rz_axis_safe_profile(
    output: Path,
    stem: str,
    half_profile: np.ndarray,
) -> Path:
    """Write a closed full meridional polygon for an r>=0 RZ domain.

    The usual half-profile is closed by an artificial edge on r=0.  That edge
    is harmless for inside/outside testing away from the axis, but a solid
    ghost cell in the first radial ring can select it as its closest IBM
    surface instead of the nozzle-exit disk.  Mirroring the physical outline
    through r=0 removes that artificial axis edge while preserving exactly the
    same solid subset in the computational half-plane r>=0.
    """
    if len(half_profile) < 4:
        raise RuntimeError("RZ axis-safe profile needs at least four vertices")
    if not (
        abs(float(half_profile[0, 0])) <= 1.0e-14
        and abs(float(half_profile[1, 0])) <= 1.0e-14
    ):
        raise RuntimeError("half-profile does not start with its axis closure")

    positive_outline = half_profile[1:].copy()
    negative_outline = positive_outline[1:][::-1].copy()
    negative_outline[:, 0] *= -1.0
    full_profile = np.vstack((positive_outline, negative_outline))
    if np.any(
        np.linalg.norm(
            np.roll(full_profile, -1, axis=0) - full_profile, axis=1
        )
        <= 1.0e-14
    ):
        raise RuntimeError("RZ axis-safe profile contains a zero-length edge")

    dat_path = output / f"{stem}.dat"
    np.savetxt(dat_path, full_profile * 1.0e-3, fmt="%.15e")
    return dat_path


def draw_profile(
    axis: plt.Axes,
    profile: np.ndarray,
    segments: dict[str, np.ndarray],
    *,
    mirrored: bool,
    fill_color: str,
) -> None:
    x_profile = profile[:, 1]
    r_profile = profile[:, 0]
    axis.fill(x_profile, r_profile, color=fill_color, zorder=0)
    if mirrored:
        axis.fill(x_profile, -r_profile, color=fill_color, zorder=0)
    for name, values in segments.items():
        x = values[:, 1]
        r = values[:, 0]
        style = "--" if name == "axis_closure" else "-"
        width = 3.2 if "jet_patch" in name else 2.0
        axis.plot(
            x,
            r,
            style,
            color=COLORS[name],
            lw=width,
            zorder=3,
        )
        if mirrored:
            axis.plot(
                x,
                -r,
                style,
                color=COLORS[name],
                lw=width,
                zorder=3,
            )
    axis.axvline(0.0, color="0.35", lw=0.9, ls=":")
    axis.grid(alpha=0.18)
    axis.set_aspect("equal", adjustable="box")


def render_single(
    output: Path,
    stem: str,
    title: str,
    subtitle: str,
    g: Geometry,
    profile: np.ndarray,
    segments: dict[str, np.ndarray],
) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(16, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=(2.05, 1.0))
    full = fig.add_subplot(grid[0, 0])
    zoom = fig.add_subplot(grid[0, 1])
    draw_profile(
        full,
        profile,
        segments,
        mirrored=True,
        fill_color="#dfe5eb",
    )
    draw_profile(
        zoom,
        profile,
        segments,
        mirrored=False,
        fill_color="#dfe5eb",
    )
    full.axhline(0.0, color="0.35", lw=0.8)
    full.set_xlim(-10.0, g.aft_x_mm + 7.0)
    full.set_ylim(-70.0, 70.0)
    full.set_xlabel("x from theoretical sharp point [mm]")
    full.set_ylabel("signed radius [mm]")
    full.set_title(title)
    zoom.set_xlim(-1.0, 31.0)
    zoom.set_ylim(-1.0, 67.0)
    zoom.set_xlabel("x from theoretical sharp point [mm]")
    zoom.set_ylabel("r [mm] (Cerisse RZ half-plane)")
    zoom.set_title("Forebody and jet-patch detail")

    if "throat_jet_patch" in segments:
        zoom.annotate(
            "choked throat patch\n"
            f"x={g.throat_x_mm:.6f} mm\n"
            f"Rt={g.throat_r_mm:.5f} mm",
            (g.throat_x_mm, 0.55 * g.throat_r_mm),
            xytext=(16.5, 1.2),
            fontsize=8,
            arrowprops={"arrowstyle": "->", "lw": 0.8},
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.78,
                "pad": 1.5,
            },
        )
        zoom.annotate(
            "resolved 14.981° divergent wall",
            (
                0.5 * (g.virtual_exit_x_mm + g.throat_x_mm),
                0.5 * (g.virtual_exit_r_mm + g.throat_r_mm),
            ),
            xytext=(5.6, 10.0),
            fontsize=8,
            arrowprops={"arrowstyle": "->", "lw": 0.8},
        )
    else:
        zoom.annotate(
            "direct supersonic virtual-exit patch\n"
            f"x={g.virtual_exit_x_mm:.6f} mm\n"
            f"Re={g.virtual_exit_r_mm:.5f} mm",
            (g.virtual_exit_x_mm, 0.55 * g.virtual_exit_r_mm),
            xytext=(7.5, 3.0),
            fontsize=8,
            arrowprops={"arrowstyle": "->", "lw": 0.8},
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.78,
                "pad": 1.5,
            },
        )
        zoom.text(
            7.3,
            8.5,
            "internal divergent passage is not resolved",
            fontsize=8,
            color="#006b6b",
        )

    zoom.annotate(
        "common reconstructed external OML",
        (g.sphere_cone_x_mm, g.sphere_cone_r_mm),
        xytext=(8.2, 13.3),
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )
    full.annotate(
        "artificial aft closure",
        (g.aft_x_mm, 0.70 * g.body_radius_mm),
        xytext=(-122, -4),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )

    handles = [
        Patch(facecolor="#dfe5eb", edgecolor="none", label="IBM solid"),
        Line2D(
            [0], [0], color="#1769aa", lw=2, label="external mold line"
        ),
        Line2D(
            [0],
            [0],
            color="#2ca02c",
            lw=2,
            label="R2.54 shoulder",
        ),
    ]
    if "throat_jet_patch" in segments:
        handles.extend(
            [
                Line2D(
                    [0],
                    [0],
                    color="#d62728",
                    lw=2,
                    label="resolved divergent wall",
                ),
                Line2D(
                    [0],
                    [0],
                    color="#f28e2b",
                    lw=3,
                    label="choked throat patch",
                ),
            ]
        )
    else:
        handles.append(
            Line2D(
                [0],
                [0],
                color="#00a6a6",
                lw=3,
                label="supersonic virtual-exit patch",
            )
        )
    handles.extend(
        [
            Line2D(
                [0],
                [0],
                color="#616161",
                lw=2,
                ls="--",
                label="RZ axis closure",
            ),
            Line2D(
                [0],
                [0],
                color="#9467bd",
                lw=2,
                label="artificial aft cap",
            ),
        ]
    )
    full.legend(handles=handles, loc="lower center", ncol=2, fontsize=8)
    fig.suptitle(
        "NASA Test 1853 Run 165 — 2-D RZ SRP JET REPRESENTATION\n"
        f"{subtitle}; geometry preview, pure shared-GP/full-Cartesian IBM",
        fontsize=13,
    )
    png = output / f"{stem}.png"
    svg = output / f"{stem}.svg"
    fig.savefig(png, dpi=240)
    fig.savefig(svg, metadata={"Date": None})
    plt.close(fig)
    return png, svg


def render_comparison(
    output: Path,
    g: Geometry,
    throat_profile: np.ndarray,
    throat_segments: dict[str, np.ndarray],
    exit_profile: np.ndarray,
    exit_segments: dict[str, np.ndarray],
) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(19, 7.8), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=(1.2, 1.2, 0.85))
    axes = [fig.add_subplot(grid[0, index]) for index in range(3)]

    draw_profile(
        axes[0],
        throat_profile,
        throat_segments,
        mirrored=True,
        fill_color="#dfe5eb",
    )
    draw_profile(
        axes[1],
        exit_profile,
        exit_segments,
        mirrored=True,
        fill_color="#dfe5eb",
    )
    for axis in axes[:2]:
        axis.axhline(0.0, color="0.35", lw=0.8)
        axis.set_xlim(-7.0, g.aft_x_mm + 6.0)
        axis.set_ylim(-70.0, 70.0)
        axis.set_xlabel("x [mm from TSP]")
    axes[0].set_ylabel("signed radius [mm]")
    axes[0].set_title("Case T — resolved nozzle, M=1 at throat")
    axes[1].set_title("Case E — direct state at virtual exit")

    overlay = axes[2]
    overlay.fill(
        exit_profile[:, 1],
        exit_profile[:, 0],
        color="#e8edf2",
        alpha=0.65,
    )
    for name in ("sphere_oml", "cone_oml", "shoulder_oml"):
        values = exit_segments[name]
        overlay.plot(
            values[:, 1],
            values[:, 0],
            color="#1769aa",
            lw=2.2,
        )
    for name in ("nozzle_divergent_wall", "throat_jet_patch"):
        values = throat_segments[name]
        overlay.plot(
            values[:, 1],
            values[:, 0],
            color=COLORS[name],
            lw=3.0 if "patch" in name else 2.4,
        )
    values = exit_segments["virtual_exit_jet_patch"]
    overlay.plot(
        values[:, 1],
        values[:, 0],
        color=COLORS["virtual_exit_jet_patch"],
        lw=3.0,
    )
    overlay.axhline(0.0, color="0.35", lw=0.8)
    overlay.axvline(0.0, color="0.35", lw=0.9, ls=":")
    overlay.grid(alpha=0.18)
    overlay.set_aspect("equal", adjustable="box")
    overlay.set_xlim(-1.0, 25.5)
    overlay.set_ylim(-0.5, 65.0)
    overlay.set_xlabel("x [mm from TSP]")
    overlay.set_ylabel("r [mm]")
    overlay.set_title("Forebody overlay")
    overlay.annotate(
        "Case E exit patch",
        (g.virtual_exit_x_mm, 0.55 * g.virtual_exit_r_mm),
        xytext=(7.0, 3.2),
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )
    overlay.annotate(
        "Case T throat patch",
        (g.throat_x_mm, 0.55 * g.throat_r_mm),
        xytext=(15.7, 6.8),
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )
    overlay.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="#d62728",
                lw=2.4,
                label="Case T divergent wall",
            ),
            Line2D(
                [0],
                [0],
                color="#f28e2b",
                lw=3,
                label="Case T M=1 throat",
            ),
            Line2D(
                [0],
                [0],
                color="#00a6a6",
                lw=3,
                label="Case E supersonic exit",
            ),
            Line2D(
                [0],
                [0],
                color="#1769aa",
                lw=2.2,
                label="common external OML",
            ),
        ],
        loc="upper left",
        fontsize=8,
    )
    fig.suptitle(
        "Run 165 SRP RZ geometry A/B — resolved throat injection versus "
        "direct virtual-exit injection\n"
        "Only the internal jet representation changes; external OML and "
        "pure shared-GP/full-Cartesian method are identical",
        fontsize=13,
    )
    png = output / "run165_srp_jet_plane_comparison.png"
    svg = output / "run165_srp_jet_plane_comparison.svg"
    fig.savefig(png, dpi=240)
    fig.savefig(svg, metadata={"Date": None})
    plt.close(fig)
    return png, svg


def main() -> int:
    output = DEFAULT_OUTPUT.resolve()
    output.mkdir(parents=True, exist_ok=True)
    parameters = load_json(PARAMETER_FILE)
    g = derive_geometry(parameters)

    throat_segments = build_segments(g)
    exit_segments = build_virtual_exit_segments(g)
    throat_profile, throat_outgoing = assemble_profile(throat_segments)
    exit_profile, exit_outgoing = assemble_profile(exit_segments)
    throat_qa = validate_profile(throat_profile, throat_segments)
    exit_qa = validate_profile(exit_profile, exit_segments)

    shared_oml_names = (
        "sphere_oml",
        "cone_oml",
        "shoulder_oml",
        "cylinder_oml",
        "aft_cap",
    )
    shared_oml_max_error_mm = 0.0
    for name in shared_oml_names:
        if throat_segments[name].shape != exit_segments[name].shape:
            raise RuntimeError(f"shared OML shape differs for {name}")
        error_mm = float(
            np.max(
                np.abs(throat_segments[name] - exit_segments[name])
            )
        )
        shared_oml_max_error_mm = max(
            shared_oml_max_error_mm, error_mm
        )
    if shared_oml_max_error_mm != 0.0:
        raise RuntimeError(
            "throat/exit external OML arrays are not identical; "
            f"maximum error={shared_oml_max_error_mm} mm"
        )

    reference_throat = np.loadtxt(REFERENCE_THROAT_DAT) * 1.0e3
    if reference_throat.shape != throat_profile.shape:
        raise RuntimeError(
            "throat profile shape no longer matches the audited Run165 DAT"
        )
    throat_reference_error_mm = float(
        np.max(np.abs(reference_throat - throat_profile))
    )
    if throat_reference_error_mm > 1.0e-11:
        raise RuntimeError(
            "throat profile differs from the audited Run165 DAT by "
            f"{throat_reference_error_mm} mm"
        )

    throat_dat, throat_csv = write_profile(
        output,
        "run165_srp_throat_choked_profile",
        throat_profile,
        throat_outgoing,
    )
    if sha256(throat_dat) != sha256(REFERENCE_THROAT_DAT):
        raise RuntimeError(
            "regenerated throat DAT is not byte-identical to the audited file"
        )
    exit_dat, exit_csv = write_profile(
        output,
        "run165_srp_virtual_exit_profile",
        exit_profile,
        exit_outgoing,
    )
    exit_axis_safe_dat = write_rz_axis_safe_profile(
        output,
        "run165_srp_virtual_exit_profile_rz_axis_safe",
        exit_profile,
    )

    throat_png, throat_svg = render_single(
        output,
        "run165_srp_throat_choked",
        "Case T: resolved divergent passage",
        "choked M=1 state imposed at the reconstructed N1 throat",
        g,
        throat_profile,
        throat_segments,
    )
    exit_png, exit_svg = render_single(
        output,
        "run165_srp_virtual_exit",
        "Case E: unresolved internal nozzle",
        "supersonic isentropic state imposed directly at the virtual exit",
        g,
        exit_profile,
        exit_segments,
    )
    comparison_png, comparison_svg = render_comparison(
        output,
        g,
        throat_profile,
        throat_segments,
        exit_profile,
        exit_segments,
    )

    beta = math.radians(g.nozzle_divergent_half_angle_deg)
    divergent_cavity_mm3 = (
        math.pi
        * (g.virtual_exit_r_mm**3 - g.throat_r_mm**3)
        / (3.0 * math.tan(beta))
    )
    throat_exact_volume_mm3 = exact_volume_mm3(g)
    exit_exact_volume_mm3 = (
        throat_exact_volume_mm3 + divergent_cavity_mm3
    )
    throat_qa["analytic_revolution_volume_mm3"] = (
        throat_exact_volume_mm3
    )
    throat_qa["polyline_volume_relative_error"] = abs(
        float(throat_qa["polyline_revolution_volume_mm3"])
        - throat_exact_volume_mm3
    ) / throat_exact_volume_mm3
    exit_qa["analytic_revolution_volume_mm3"] = exit_exact_volume_mm3
    exit_qa["polyline_volume_relative_error"] = abs(
        float(exit_qa["polyline_revolution_volume_mm3"])
        - exit_exact_volume_mm3
    ) / exit_exact_volume_mm3

    output_files = [
        throat_dat,
        throat_csv,
        exit_dat,
        exit_csv,
        exit_axis_safe_dat,
        throat_png,
        throat_svg,
        exit_png,
        exit_svg,
        comparison_png,
        comparison_svg,
    ]
    metadata = {
        "status": "GEOMETRY_VISUALIZATION_QA_PASS",
        "purpose": (
            "Run165 RZ comparison of a resolved divergent nozzle with "
            "choked-throat injection and an unresolved nozzle with direct "
            "virtual-exit injection"
        ),
        "method_scope": {
            "family": "pure_shared_gp_full_cartesian",
            "impact": (
                "case geometry and jet-patch thermodynamic target only"
            ),
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
            "retropropulsive_jet": "-x",
        },
        "geometry_parameters": asdict(g),
        "shared_external_oml": {
            "segment_names": list(shared_oml_names),
            "coordinate_max_error_mm": shared_oml_max_error_mm,
            "identical": True,
        },
        "case_T_throat_choked": {
            "description": (
                "Audited Run165 axisymmetric geometry with the divergent "
                "passage resolved and a jet patch on the throat disk"
            ),
            "jet_patch_x_mm": g.throat_x_mm,
            "jet_patch_radius_mm": g.throat_r_mm,
            "divergent_wall_half_angle_deg": (
                g.nozzle_divergent_half_angle_deg
            ),
            "reference_dat_sha256": sha256(REFERENCE_THROAT_DAT),
            "reference_coordinate_max_error_mm": (
                throat_reference_error_mm
            ),
            "qa": throat_qa,
        },
        "case_E_virtual_exit": {
            "description": (
                "Same external OML, with the internal divergent passage "
                "removed from the fluid domain and a jet patch on the "
                "virtual-exit disk"
            ),
            "jet_patch_x_mm": g.virtual_exit_x_mm,
            "jet_patch_radius_mm": g.virtual_exit_r_mm,
            "filled_divergent_cavity_volume_mm3": (
                divergent_cavity_mm3
            ),
            "qa": exit_qa,
            "rz_axis_safe_dat": exit_axis_safe_dat.name,
            "rz_axis_safe_semantics": (
                "full mirrored meridional polygon; identical r>=0 solid "
                "region with no artificial r=0 IBM surface edge"
            ),
        },
        "semantic_guards": [
            (
                "Case E is a virtual-exit boundary representation of Run165, "
                "not a claim that the NASA model lacked a nozzle."
            ),
            (
                "Jet-OFF is different from Case E: jet-OFF retains whichever "
                "geometry is selected and changes the jet patch to a wall."
            ),
            (
                "The throat profile resolves only the reconstructed "
                "divergent passage and throat face; no convergent passage or "
                "plenum is invented."
            ),
            (
                "Artificial axis and aft-cap closures are excluded from "
                "physical surface loads."
            ),
        ],
        "source_sha256": {
            "geometry_parameters.yaml": sha256(PARAMETER_FILE),
            "build_audit.json": sha256(BUILD_AUDIT),
            "audited_throat_profile.dat": sha256(
                REFERENCE_THROAT_DAT
            ),
            "generator": sha256(Path(__file__).resolve()),
        },
        "outputs_sha256": {
            path.name: sha256(path) for path in output_files
        },
    }
    metadata_path = output / "run165_srp_jet_plane_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    print(f"wrote {output}")
    print(
        "QA "
        f"throat_sha={sha256(throat_dat)} "
        f"exit_vertices={len(exit_profile)} "
        f"exit_volume_rel={exit_qa['polyline_volume_relative_error']:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
