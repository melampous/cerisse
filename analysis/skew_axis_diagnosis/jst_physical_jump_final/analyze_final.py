#!/usr/bin/env python3
"""First-eight-ring audit for the final physical-jump R-Z JST repair."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[3]
ANALYSIS = REPO / "analysis/skew_axis_diagnosis"
RUN165 = REPO / "exm/underexpanded_jet/2d/run165_analytic_exit_matrix"
HELPERS = ANALYSIS / "field_audit/compare_first8_rings.py"

CASES = {
    "skew_final_physical_jump_500us": HERE / "fallback_on_500us/plt01250",
    "skew_old_baseline_500us": ANALYSIS / "baseline/plt01250",
    "skew_old_axis_pair_500us": ANALYSIS / "fixed_axis_pair_500us/plt01250",
    "hllc_muscl_500us": RUN165 / "results_hllc_muscl_o2_rz_cfl02_500us/plt01178",
}


def load_helpers():
    spec = importlib.util.spec_from_file_location("first8_helpers", HELPERS)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {HELPERS}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_heatmap(helper, data: dict[str, dict]) -> None:
    names = list(CASES)
    fields = ("pressure", "density", "temperature", "mach")
    floors = {
        "pressure": 1.0,
        "density": 1.0e-12,
        "temperature": 1.0e-12,
        "mach": 1.0e-12,
    }
    figure, axes = plt.subplots(
        len(fields), len(names), figsize=(15.5, 10.5), sharex=True,
        sharey=True, squeeze=False,
    )
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
    colorbar = figure.colorbar(image, ax=axes, fraction=0.02, pad=0.015)
    colorbar.set_label("absolute radial local-extremum residual")
    figure.suptitle("Run165 at 500 microseconds: first-eight-ring filament audit")
    figure.savefig(HERE / "radial_filament_heatmaps_500us.png", dpi=180)
    plt.close(figure)


def main() -> None:
    helper = load_helpers()
    helper.OUTPUT = HERE
    axis = helper.load_axis_module()
    data = {name: helper.first8(axis, path) for name, path in CASES.items()}
    helper.export_raw(data)
    helper.export_filament_metrics(data)
    helper.export_axis_metrics(axis, data)
    helper.export_gradient_locus(data)
    helper.plot_profiles(data, list(CASES), HERE / "first8_profiles_500us.png")
    make_heatmap(helper, data)

    import csv

    with (HERE / "radial_filament_metrics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    for case in CASES:
        print(case)
        for field in helper.FIELDS:
            selected = [r for r in rows if r["case"] == case and r["field"] == field]
            worst = max(selected, key=lambda r: float(r["p95_abs_radial_extremum"]))
            print(
                f"  {field}: ring={worst['ring']} "
                f"p95={float(worst['p95_abs_radial_extremum']):.6g} "
                f"longest>2%={float(worst['longest_over_2pct_mm']):.6g} mm"
            )


if __name__ == "__main__":
    main()
