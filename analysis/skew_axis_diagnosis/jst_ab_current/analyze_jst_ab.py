#!/usr/bin/env python3
"""Matched-time first-eight-ring audit for the central/JST Run165 A/B."""

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
FIELD_AUDIT = REPO / "analysis/skew_axis_diagnosis/field_audit/compare_first8_rings.py"

CASES = {
    "central_K4_fallback_on_80us": HERE / "central_k4_fallback_on_80us/plt00200",
    "jst_K4_fallback_on_80us": HERE / "jst_k4_fallback_on_80us/plt00200",
}


def load_helpers():
    spec = importlib.util.spec_from_file_location("first8_helpers", FIELD_AUDIT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FIELD_AUDIT}")
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
    axis = helper.load_axis_module()
    data = {name: helper.first8(axis, path) for name, path in CASES.items()}

    rows: list[dict] = []
    floors = {
        "pressure": 1.0,
        "density": 1.0e-12,
        "temperature": 1.0e-12,
        "radial_velocity": 1.0,
        "mach": 1.0e-12,
    }
    for case, item in data.items():
        active = (item["z"] >= -0.080) & (item["z"] <= -0.005)
        for field in helper.FIELDS:
            residual = np.abs(
                helper.radial_extremum_residual(item["fields"][field], floors[field])
            )
            for local_ring in range(residual.shape[0]):
                ring = local_ring + 1
                values = residual[local_ring, active]
                above = values > 0.02
                rows.append(
                    {
                        "case": case,
                        "time_s": f"{float(item['dataset'].current_time):.17g}",
                        "level": item["level"],
                        "dr_mm": f"{item['dr'] * 1.0e3:.12g}",
                        "field": field,
                        "ring": ring,
                        "r_mm": f"{item['r'][ring] * 1.0e3:.12g}",
                        "peak_abs_radial_extremum": f"{np.max(values):.17g}",
                        "p95_abs_radial_extremum": f"{np.quantile(values, 0.95):.17g}",
                        "fraction_over_2pct": f"{np.mean(above):.17g}",
                        "longest_over_2pct_mm": (
                            f"{helper.longest_true_run(above, item['dz']) * 1.0e3:.17g}"
                        ),
                    }
                )
    write_csv(HERE / "radial_filament_metrics_80us.csv", rows)

    names = list(CASES)
    figure, axes = plt.subplots(
        len(helper.FIELDS), len(names), figsize=(9.5, 12.5), sharex=True,
        squeeze=False,
    )
    colors = plt.cm.viridis(np.linspace(0.03, 0.97, 8))
    for column, name in enumerate(names):
        item = data[name]
        z_mm = item["z"] * 1.0e3
        active = (z_mm >= -80.0) & (z_mm <= 0.5)
        for row, field in enumerate(helper.FIELDS):
            values = item["fields"][field] / helper.NORMALIZER[field]
            for ring in range(8):
                axes[row, column].plot(
                    z_mm[active], values[ring, active], color=colors[ring], lw=1.0,
                    label=f"ring {ring}",
                )
            axes[row, column].grid(True, alpha=0.2)
            if column == 0:
                axes[row, column].set_ylabel(helper.LABELS[field])
            if row == 0:
                axes[row, column].set_title(name.replace("_", "\n"), fontsize=9)
            if row == len(helper.FIELDS) - 1:
                axes[row, column].set_xlabel("axial z [mm]")
    handles, labels = axes[0, -1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="center right", bbox_to_anchor=(0.995, 0.5))
    figure.suptitle("Run165 at 80 microseconds: first eight uniform-L1 rings")
    figure.tight_layout(rect=(0.0, 0.0, 0.91, 0.98))
    figure.savefig(HERE / "first8_profiles_80us.png", dpi=180)
    plt.close(figure)

    fields = ("pressure", "density", "temperature", "mach")
    figure, axes = plt.subplots(len(fields), len(names), figsize=(9.5, 9.0),
                                sharex=True, sharey=True, squeeze=False)
    for column, name in enumerate(names):
        item = data[name]
        z_mm = item["z"] * 1.0e3
        active = (z_mm >= -80.0) & (z_mm <= -5.0)
        for row, field in enumerate(fields):
            residual = np.abs(helper.radial_extremum_residual(
                item["fields"][field], floors[field]
            ))
            image = axes[row, column].imshow(
                residual[:, active], origin="lower", aspect="auto",
                extent=(z_mm[active][0], z_mm[active][-1], 0.5, 6.5),
                vmin=0.0, vmax=0.20, cmap="magma",
            )
            axes[row, column].set_yticks(range(1, 7))
            if column == 0:
                axes[row, column].set_ylabel(f"{field}\ncentre ring")
            if row == 0:
                axes[row, column].set_title(name.replace("_", "\n"), fontsize=9)
            if row == len(fields) - 1:
                axes[row, column].set_xlabel("axial z [mm]")
    colorbar = figure.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    colorbar.set_label("absolute radial local-extremum residual")
    figure.suptitle("Run165 at 80 microseconds: grid-filament diagnostic")
    figure.savefig(HERE / "radial_filament_heatmaps_80us.png", dpi=180)
    plt.close(figure)

    for case in names:
        print(case)
        for field in fields:
            candidates = [r for r in rows if r["case"] == case and r["field"] == field]
            worst = max(candidates, key=lambda r: float(r["p95_abs_radial_extremum"]))
            print(
                f"  {field}: ring={worst['ring']} "
                f"p95={float(worst['p95_abs_radial_extremum']):.6g} "
                f"longest>2%={float(worst['longest_over_2pct_mm']):.6g} mm"
            )


if __name__ == "__main__":
    main()
