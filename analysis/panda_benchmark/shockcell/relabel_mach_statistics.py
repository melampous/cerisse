#!/usr/bin/env python3
"""Relabel the mesh-study Mach fluctuation diagnostic without changing data.

The archived panels use total resolved velocity variance divided by the local
mean sound speed.  They therefore show M_{u,rms}, not the r.m.s. of the
instantaneous Mach number.  This script changes only the two publication
labels in each existing high-resolution figure and writes matching PDF files.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
FILES = (
    "paperfig_mach2d_main",
    "paperfig_mach3d_main",
    "paperfig_mach2d_L0to6",
    "paperfig_mach3d_L1to5",
)


def relabel(stem):
    path = HERE / f"{stem}.png"
    image = Image.open(path).convert("RGBA")
    width, height = image.size

    # The column and colour-bar labels lie on the white figure margin. Erase
    # only those margin boxes; no plotted field, contour, tick, or colour
    # scale is altered.
    draw = ImageDraw.Draw(image)
    draw.rectangle(
        (int(0.455 * width), 0, int(0.555 * width), int(0.031 * height)),
        fill="white",
    )
    draw.rectangle(
        (int(0.760 * width), 0, int(0.910 * width), int(0.031 * height)),
        fill="white",
    )
    draw.rectangle(
        (
            int(0.335 * width),
            int(0.979 * height),
            int(0.535 * width),
            height,
        ),
        fill="white",
    )
    draw.rectangle(
        (
            int(0.770 * width),
            int(0.979 * height),
            int(0.900 * width),
            height,
        ),
        fill="white",
    )

    overlay = plt.figure(
        figsize=(width / 100.0, height / 100.0),
        dpi=100,
        facecolor=(1, 1, 1, 0),
    )
    overlay.text(
        0.505,
        0.995,
        r"$M_{\overline{\mathbf{q}}}$",
        ha="center",
        va="top",
        fontsize=18,
    )
    overlay.text(
        0.835,
        0.995,
        r"$M_{u,\mathrm{rms}}$",
        ha="center",
        va="top",
        fontsize=18,
    )
    overlay.text(
        0.435,
        0.004,
        r"$M,\;M_{\overline{\mathbf{q}}}$",
        ha="center",
        va="bottom",
        fontsize=18,
    )
    overlay.text(
        0.835,
        0.004,
        r"$M_{u,\mathrm{rms}}$",
        ha="center",
        va="bottom",
        fontsize=18,
    )
    overlay.canvas.draw()
    rgba = Image.frombuffer(
        "RGBA",
        overlay.canvas.get_width_height(),
        overlay.canvas.buffer_rgba(),
        "raw",
        "RGBA",
        0,
        1,
    )
    plt.close(overlay)

    image = Image.alpha_composite(image, rgba).convert("RGB")
    image.save(path, dpi=(300, 300))
    image.save(HERE / f"{stem}.pdf", "PDF", resolution=300.0)
    print(f"updated {path.name} and {stem}.pdf")


for name in FILES:
    relabel(name)
