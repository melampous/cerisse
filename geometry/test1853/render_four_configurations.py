#!/usr/bin/env python3
"""Render the four conceptual 3-D Test 1853 nozzle configurations."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-test1853-four-configs")

import pyvista as pv
from PIL import Image


HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"


@dataclass(frozen=True)
class Model:
    stem: str
    title: str
    root: Path
    output_stem: str

    @property
    def stl(self) -> Path:
        return self.root / "stl" / f"{self.stem}.stl"

    @property
    def patch_dir(self) -> Path:
        return self.root / "patches" / self.stem


def models(ring_root: Path) -> list[Model]:
    return [
        Model(
            "test1853_run165_single",
            "1 center nozzle",
            OUTPUTS,
            "configuration_1_center",
        ),
        Model(
            "test1853_run262_263_tri_n3at120",
            "3 peripheral nozzles",
            OUTPUTS,
            "configuration_3_peripheral",
        ),
        Model(
            "test1853_quad_n3at120",
            "1 center + 3 peripheral",
            OUTPUTS,
            "configuration_center_plus_3",
        ),
        Model(
            "test1853_four_peripheral_ring90",
            "4 peripheral nozzles (90 deg ring)",
            ring_root,
            "configuration_4_peripheral_ring",
        ),
    ]


def require_inputs(items: list[Model]) -> None:
    missing: list[Path] = []
    for item in items:
        if not item.stl.is_file():
            missing.append(item.stl)
        patches = sorted(item.patch_dir.glob("nozzle_*_throat_cap.stl"))
        if not patches:
            missing.append(item.patch_dir)
    if missing:
        raise FileNotFoundError("Missing render inputs: " + ", ".join(map(str, missing)))


def add_model(plotter: pv.Plotter, item: Model, title: bool = True) -> None:
    body = pv.read(item.stl)
    forebody = body.clip(normal=(1.0, 0.0, 0.0), origin=(82.0, 0.0, 0.0), invert=True)
    plotter.add_mesh(
        forebody,
        color="#cbd7de",
        smooth_shading=True,
        show_edges=False,
        specular=0.28,
        specular_power=22,
    )
    for patch in sorted(item.patch_dir.glob("nozzle_*_throat_cap.stl")):
        plotter.add_mesh(
            pv.read(patch),
            color="#d64b32",
            smooth_shading=False,
            show_edges=True,
            edge_color="#7c2117",
            line_width=0.6,
        )
    if title:
        plotter.add_text(item.title, position="upper_left", font_size=15, color="#17242c")
        plotter.add_text("red: throat boundary", position="lower_left", font_size=9, color="#6b281f")
    plotter.set_background("#f7f8f8")
    plotter.camera_position = [
        (-250.0, 95.0, 75.0),
        (20.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
    ]
    plotter.camera.zoom(1.05)
    plotter.enable_anti_aliasing("ssaa")


def render_individual(item: Model, output: Path) -> None:
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1100))
    add_model(plotter, item)
    plotter.add_axes(line_width=2, color="#27343b")
    plotter.screenshot(output / f"{item.output_stem}.png")
    plotter.close()


def render_montage(items: list[Model], output: Path) -> None:
    images: list[Image.Image] = []
    for item in items:
        with Image.open(output / f"{item.output_stem}.png") as source:
            images.append(source.convert("RGB"))
    width, height = images[0].size
    gutter = 12
    canvas = Image.new(
        "RGB",
        (2 * width + gutter, 2 * height + gutter),
        color=(174, 185, 191),
    )
    for index, source in enumerate(images):
        x = (index % 2) * (width + gutter)
        y = (index // 2) * (height + gutter)
        canvas.paste(source, (x, y))
    canvas.save(output / "four_3d_nozzle_configurations.png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ring-root",
        type=Path,
        default=OUTPUTS / "four_peripheral_ring90_build",
    )
    parser.add_argument("--output", type=Path, default=OUTPUTS / "figures")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    items = models(args.ring_root.resolve())
    require_inputs(items)
    args.output.mkdir(parents=True, exist_ok=True)
    for item in items:
        render_individual(item, args.output)
    render_montage(items, args.output)
    print(f"PASS: rendered {len(items)} configurations under {args.output}")


if __name__ == "__main__":
    main()
