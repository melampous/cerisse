#!/usr/bin/env python3
"""Render three orthographic views and one oblique view of the tilted nozzle."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853-surface-normal")

import numpy as np
import pyvista as pv
from PIL import Image


HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE / "outputs" / "single_peripheral_surface_normal_build"
MODEL = "test1853_single_peripheral_surface_normal"


def load_inputs(root: Path) -> tuple[pv.PolyData, pv.PolyData, pv.PolyData, dict]:
    stl = root / "stl" / f"{MODEL}.stl"
    throat = root / "patches" / MODEL / "nozzle_1_throat_cap.stl"
    virtual_exit = (
        root
        / "reference"
        / MODEL
        / "nozzle_1_virtual_exit_REFERENCE_ONLY.stl"
    )
    audit_path = root / "reports" / "build_audit.json"
    missing = [path for path in (stl, throat, virtual_exit, audit_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing render inputs: " + ", ".join(map(str, missing)))

    audit = json.loads(audit_path.read_text())[MODEL]
    nozzle = audit["nozzles"][0]
    return pv.read(stl), pv.read(throat), pv.read(virtual_exit), nozzle


def add_geometry(
    plotter: pv.Plotter,
    body: pv.PolyData,
    throat: pv.PolyData,
    virtual_exit: pv.PolyData,
    nozzle: dict,
    title: str,
) -> None:
    forebody = body.clip(
        normal=(1.0, 0.0, 0.0), origin=(82.0, 0.0, 0.0), invert=True
    )
    plotter.add_mesh(
        forebody,
        color="#cbd7de",
        smooth_shading=True,
        show_edges=False,
        specular=0.24,
        specular_power=20,
    )
    plotter.add_mesh(
        throat,
        color="#d64b32",
        smooth_shading=False,
        show_edges=True,
        edge_color="#6f1d15",
        line_width=1.0,
    )
    virtual_exit_outline = virtual_exit.extract_feature_edges(
        boundary_edges=True,
        non_manifold_edges=False,
        feature_edges=False,
        manifold_edges=False,
        clear_data=True,
    )
    plotter.add_mesh(
        virtual_exit_outline,
        color="#23818d",
        line_width=5.0,
        render_lines_as_tubes=True,
        lighting=False,
    )

    throat_center = np.asarray(nozzle["throat_center_mm"], dtype=float)
    exit_center = np.asarray(nozzle["virtual_exit_center_mm"], dtype=float)
    jet = np.asarray(nozzle["jet_direction_unit"], dtype=float)
    jet /= np.linalg.norm(jet)
    passage_axis = pv.Line(throat_center, exit_center, resolution=1)
    arrow = pv.Arrow(
        start=exit_center + 0.8 * jet,
        direction=jet,
        scale=25.0,
        shaft_radius=0.025,
        tip_radius=0.08,
        tip_length=0.25,
    )
    plotter.add_mesh(
        passage_axis,
        color="#27343b",
        line_width=5.0,
        render_lines_as_tubes=True,
    )
    plotter.add_mesh(arrow, color="#16825d", smooth_shading=True)
    plotter.add_text(title, position="upper_left", font_size=14, color="#17242c")
    plotter.add_text(
        "teal: virtual-exit reference\ngreen: jet direction",
        position="lower_left",
        font_size=9,
        color="#29434a",
    )
    plotter.set_background("#f7f8f8")


def render(root: Path, output: Path) -> None:
    body, throat, virtual_exit, nozzle = load_inputs(root)
    views = [
        (
            "side_xz",
            "Side: x-z",
            [(35.0, -260.0, 0.0), (35.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
        ),
        (
            "top_xy",
            "Top: x-y",
            [(35.0, 0.0, 260.0), (35.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        ),
        (
            "front_yz",
            "Front: y-z",
            [(-255.0, 0.0, 0.0), (25.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
        ),
        (
            "oblique",
            "Oblique",
            [(-185.0, -175.0, 135.0), (30.0, 0.0, 5.0), (0.0, 0.0, 1.0)],
        ),
    ]

    output.parent.mkdir(parents=True, exist_ok=True)
    panels: list[Path] = []
    for slug, title, camera in views:
        panel = output.parent / f"{output.stem}_{slug}.png"
        plotter = pv.Plotter(off_screen=True, window_size=(1100, 800))
        add_geometry(plotter, body, throat, virtual_exit, nozzle, title)
        plotter.camera_position = camera
        plotter.enable_parallel_projection()
        plotter.camera.zoom(1.08)
        plotter.enable_anti_aliasing("ssaa")
        plotter.screenshot(panel)
        plotter.close()
        panels.append(panel)

    images: list[Image.Image] = []
    for panel in panels:
        with Image.open(panel) as source:
            images.append(source.convert("RGB"))
    width, height = images[0].size
    gutter = 8
    canvas = Image.new(
        "RGB", (2 * width + gutter, 2 * height + gutter), color=(174, 185, 191)
    )
    for index, source in enumerate(images):
        x = (index % 2) * (width + gutter)
        y = (index // 2) * (height + gutter)
        canvas.paste(source, (x, y))
    canvas.save(output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_ROOT
        / "figures"
        / "test1853_single_peripheral_surface_normal_three_views_oblique.png",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    render(args.root.resolve(), args.output.resolve())
    print(f"PASS: rendered {args.output.resolve()}")


if __name__ == "__main__":
    main()
