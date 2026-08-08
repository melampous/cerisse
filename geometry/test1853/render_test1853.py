#!/usr/bin/env python3
"""Create dimensioned diagnostic views for the Test 1853 reconstruction."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from build_test1853 import OUTPUT_ROOT, configurations, derived_body, load_parameters


def body_profile(body: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    x1 = np.linspace(body["nose_tip_x"], body["sphere_cone_tangent_x"], 150)
    r1 = np.sqrt(
        np.maximum(0.0, body["nose_radius"] ** 2 - (x1 - body["sphere_center_x"]) ** 2)
    )
    x2 = np.linspace(body["sphere_cone_tangent_x"], body["cone_shoulder_tangent_x"], 300)
    r2 = x2 * math.tan(body["alpha_rad"])
    x3 = np.linspace(body["cone_shoulder_tangent_x"], body["shoulder_center_x"], 150)
    r3 = body["shoulder_center_r"] + np.sqrt(
        np.maximum(
            0.0,
            body["shoulder_radius"] ** 2 - (x3 - body["shoulder_center_x"]) ** 2,
        )
    )
    x4 = np.array([body["shoulder_center_x"], body["aft_end_x"]])
    r4 = np.full(2, body["body_radius"])
    return np.concatenate([x1, x2, x3, x4]), np.concatenate([r1, r2, r3, r4])


def plot_meridional(p: dict, out: Path) -> None:
    body = derived_body(p)
    x, r = body_profile(body)
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
    axes[0].plot(x, r, color="#173f5f", lw=1.8)
    axes[0].plot(x, -r, color="#173f5f", lw=1.8)
    axes[0].axhline(0.0, color="0.65", lw=0.7)
    axes[0].set_aspect("equal", adjustable="box")
    axes[1].plot(x, r, color="#173f5f", lw=2.0)
    for ax in axes:
        ax.set_ylabel("r [mm]")
        ax.grid(alpha=0.22)
    axes[0].set_xlim(body["nose_tip_x"] - 3, body["aft_end_x"] + 3)
    axes[0].set_title("NASA Test 1853 nominal public-data OML reconstruction")
    axes[0].set_xlabel("x from theoretical sharp point [mm]")
    axes[1].set_xlim(body["nose_tip_x"] - 1.5, 32.0)
    axes[1].set_ylim(-1.0, 67.0)
    axes[1].set_xlabel("x from theoretical sharp point [mm] — forebody zoom")
    markers = [
        ("sphere/cone", body["sphere_cone_tangent_x"], body["sphere_cone_tangent_r"], (22, 12)),
        ("cone/R2.54", body["cone_shoulder_tangent_x"], body["cone_shoulder_tangent_r"], (-82, -33)),
        ("R2.54/cylinder", body["shoulder_center_x"], body["body_radius"], (-88, 26)),
        ("aft-cover start", p["body"]["aft_cover_forward_x_mm"], body["body_radius"], (22, -15)),
    ]
    for label, xm, rm, offset in markers:
        axes[1].scatter([xm], [rm], s=18, zorder=3)
        axes[1].annotate(
            f"{label}\nx={xm:.4f}",
            (xm, rm),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
            arrowprops={"arrowstyle": "-", "lw": 0.6},
        )
    fig.savefig(out / "body_meridional_section.png", dpi=220)
    plt.close(fig)


def plot_layout(p: dict, out: Path) -> None:
    radius = p["body"]["maximum_radius_mm"]
    axis_r = p["nozzle_layout"]["peripheral_axis_radius_mm"]
    variants = p["nozzle_layout"]["mapping_variants"]
    fig, axes = plt.subplots(1, 5, figsize=(19, 4.3), constrained_layout=True)
    titles = [
        "Run 165 single",
        "Tri mapping A",
        "Tri mapping B",
        "Quad mapping A",
        "Quad mapping B",
    ]
    mappings = [
        {"1": None},
        variants["n3_at_120"],
        variants["n3_at_240"],
        {"1": None, **variants["n3_at_120"]},
        {"1": None, **variants["n3_at_240"]},
    ]
    for ax, title, mapping in zip(axes, titles, mappings):
        ax.add_patch(plt.Circle((0, 0), radius, fc="#e9eef2", ec="#173f5f", lw=1.4))
        for number, theta_deg in mapping.items():
            if number == "1":
                y, z = 0.0, 0.0
            else:
                theta = math.radians(float(theta_deg))
                y = axis_r * math.sin(theta)
                z = axis_r * math.cos(theta)
            ax.scatter([y], [z], s=100, c="#ed553b", edgecolor="k", zorder=3)
            ax.text(y, z + 5, f"N{number}", ha="center", va="bottom", fontsize=9)
        ax.arrow(0, 0, 0, radius * 0.83, width=0.35, head_width=2.4, color="0.3")
        ax.text(3, radius * 0.75, r"$\theta=0^\circ$ (+z)", fontsize=8)
        ax.set_title(title)
        ax.set_xlabel("y [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_aspect("equal")
        ax.set_xlim(-70, 70)
        ax.set_ylim(-70, 70)
        # Match NASA DWG 1168296: in the downstream-looking front view,
        # theta=90 (+y) is drawn on screen left.
        ax.invert_xaxis()
        ax.grid(alpha=0.18)
    fig.suptitle("View looking downstream (+x); all nozzle axes are parallel to x", fontsize=11)
    fig.savefig(out / "nozzle_layout_front_view.png", dpi=220)
    plt.close(fig)


def plot_lips(p: dict, root: Path, out: Path) -> None:
    name = "test1853_run262_263_tri_n3at120"
    path = root / "reference" / name / "physical_lip_reconstruction.csv"
    rows: dict[int, list[tuple[float, float, float]]] = {2: [], 3: [], 4: []}
    with path.open("r", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            number = int(row["nozzle"])
            rows[number].append(
                (float(row["psi_deg"]), float(row["x_mm"]), float(row["local_rho_mm"]))
            )
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    for number, values in rows.items():
        data = np.asarray(values)
        axes[0].plot(data[:, 0], data[:, 1], label=f"Nozzle {number}")
        axes[1].plot(data[:, 0], data[:, 2], label=f"Nozzle {number}")
        virtual_r = 0.5 * p["nozzles"][str(number)]["virtual_exit_diameter_mm"]
        axes[1].axhline(virtual_r, color="0.65", lw=0.5, ls="--")
    axes[0].set_ylabel("physical lip x [mm]")
    axes[1].set_ylabel("local lip radius rho [mm]")
    axes[1].set_xlabel("local azimuth psi [deg]; 0 = radially outward")
    axes[0].set_title("Peripheral physical lip is a skew Boolean intersection, not a circular plane")
    for ax in axes:
        ax.grid(alpha=0.22)
        ax.legend(ncol=3, fontsize=8)
    fig.savefig(out / "peripheral_physical_lip_variation.png", dpi=220)
    plt.close(fig)


def render_solids(root: Path, out: Path) -> None:
    import pyvista as pv

    models = {
        "run165_single_isometric.png": "test1853_run165_single.stl",
        "run262_263_tri_isometric.png": "test1853_run262_263_tri_n3at120.stl",
        "quad_isometric.png": "test1853_quad_n3at120.stl",
    }
    for image_name, stl_name in models.items():
        mesh = pv.read(root / "stl" / stl_name)
        plotter = pv.Plotter(off_screen=True, window_size=(1400, 850))
        plotter.set_background("white")
        plotter.add_mesh(
            mesh,
            color="#c9dce9",
            smooth_shading=True,
            show_edges=False,
            specular=0.22,
            specular_power=18,
        )
        plotter.camera_position = [(-180.0, 180.0, 140.0), (35.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
        plotter.add_axes()
        plotter.screenshot(out / image_name)
        plotter.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    p = load_parameters()
    # Force schema/config evaluation before plotting.
    configurations(p)
    out = args.output / "figures"
    out.mkdir(parents=True, exist_ok=True)
    plot_meridional(p, out)
    plot_layout(p, out)
    plot_lips(p, args.output, out)
    render_solids(args.output, out)
    print(f"PASS: wrote diagnostic figures to {out}")


if __name__ == "__main__":
    main()
