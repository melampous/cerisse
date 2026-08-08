#!/usr/bin/env python3
"""Quantify whether the Skew local-LLF collar removes or moves R-Z filaments."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[3]
HELPERS = REPO / "analysis/skew_axis_diagnosis/field_audit/compare_first8_rings.py"
CASES = {
    "fallback_off": HERE / "fallback_off_500us/plt01250",
    "fallback_on_K0": HERE / "fallback_on_k0_500us/plt01250",
    "fallback_on_K1": HERE / "fallback_on_500us/plt01250",
    "fallback_on_K2": HERE / "collar_k2_500us/plt01250",
    "fallback_on_K4": HERE / "collar_k4_500us/plt01250",
}
FIELDS = ("pressure", "density", "temperature", "radial_velocity", "mach")
FLOORS = {
    "pressure": 1.0,
    "density": 1.0e-12,
    "temperature": 1.0e-12,
    "radial_velocity": 1.0,
    "mach": 1.0e-12,
}


def load_helpers():
    spec = importlib.util.spec_from_file_location("first8_helpers", HELPERS)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {HELPERS}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    helper = load_helpers()
    helper.OUTPUT = HERE / "collar_sweep_audit"
    helper.OUTPUT.mkdir(parents=True, exist_ok=True)
    axis = helper.load_axis_module()
    data = {name: helper.first8(axis, path) for name, path in CASES.items()}

    helper.export_filament_metrics(data)
    helper.export_axis_metrics(axis, data)
    helper.plot_profiles(
        data, list(CASES), helper.OUTPUT / "first8_profiles_500us.png"
    )

    rows: list[dict] = []
    for case, item in data.items():
        z = item["z"]
        active = (z >= -0.060) & (z <= -0.020)
        for field in FIELDS:
            residual = np.abs(helper.radial_extremum_residual(
                item["fields"][field], FLOORS[field]
            ))
            for local_ring in range(residual.shape[0]):
                ring = local_ring + 1
                values = residual[local_ring, active]
                rows.append({
                    "case": case,
                    "field": field,
                    "ring": ring,
                    "r_mm": f"{item['r'][ring] * 1.0e3:.12g}",
                    "mean_abs_residual": f"{np.mean(values):.17g}",
                    "p95_abs_residual": f"{np.quantile(values, 0.95):.17g}",
                    "peak_abs_residual": f"{np.max(values):.17g}",
                })
    write_csv(helper.OUTPUT / "near_nozzle_minus60_to_minus20mm.csv", rows)

    plot_fields = ("pressure", "density", "temperature", "mach")
    names = list(CASES)
    figure, axes = plt.subplots(
        len(plot_fields), len(names), figsize=(15.5, 10.0),
        sharex=True, sharey=True, squeeze=False,
    )
    for column, name in enumerate(names):
        item = data[name]
        z_mm = item["z"] * 1.0e3
        active = (z_mm >= -80.0) & (z_mm <= -5.0)
        for row, field in enumerate(plot_fields):
            residual = np.abs(helper.radial_extremum_residual(
                item["fields"][field], FLOORS[field]
            ))
            image = axes[row, column].imshow(
                residual[:, active], origin="lower", aspect="auto",
                extent=(z_mm[active][0], z_mm[active][-1], 0.5, 6.5),
                vmin=0.0, vmax=0.15, cmap="magma",
            )
            axes[row, column].set_yticks(range(1, 7))
            if column == 0:
                axes[row, column].set_ylabel(f"{field}\ncentre ring")
            if row == 0:
                axes[row, column].set_title(name.replace("_", "\n"), fontsize=9)
            if row == len(plot_fields) - 1:
                axes[row, column].set_xlabel("axial z [mm]")
    colorbar = figure.colorbar(image, ax=axes, fraction=0.02, pad=0.015)
    colorbar.set_label("absolute radial local-extremum residual")
    figure.suptitle("Run165 500 us: fallback and hard-collar sweep")
    figure.savefig(helper.OUTPUT / "collar_sweep_heatmap_500us.png", dpi=180)
    plt.close(figure)

    for case in CASES:
        print(case)
        for field in FIELDS:
            selected = [r for r in rows if r["case"] == case and r["field"] == field]
            worst = max(selected, key=lambda r: float(r["p95_abs_residual"]))
            print(
                f"  {field}: ring={worst['ring']} r={float(worst['r_mm']):.3f}mm "
                f"p95={float(worst['p95_abs_residual']):.6g}"
            )


if __name__ == "__main__":
    main()
