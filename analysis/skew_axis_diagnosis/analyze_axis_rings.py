#!/usr/bin/env python3
"""Quantify one-cell R-Z axis defects in uniform-finest AMReX plotfiles."""

from __future__ import annotations

import argparse
import csv
import inspect
from pathlib import Path
import textwrap

import numpy as np
import yt
from yt.frontends.boxlib import data_structures as boxlib_data


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


def load_fields(plotfile: Path):
    dataset = yt.load(str(plotfile))
    level = int(dataset.index.max_level)
    dims = np.asarray(dataset.domain_dimensions, dtype=int).copy()
    dims[:2] *= 2**level
    dims[2] = 1
    grid = dataset.covering_grid(level, dataset.domain_left_edge, dims)

    def field(name: str) -> np.ndarray:
        return grid[("boxlib", name)].to_ndarray()[:, :, 0]

    rho = field("Density")
    pressure = field("pressure")
    ur = field("x_velocity")
    uz = field("y_velocity")
    temperature = field("temperature")
    sound = np.sqrt(1.4 * pressure / rho)
    mach = np.hypot(ur, uz) / sound
    return dataset, level, {
        "pressure": pressure,
        "density": rho,
        "temperature": temperature,
        "radial_velocity": ur,
        "axial_velocity": uz,
        "mach": mach,
        "sound_speed": sound,
    }


def even_axis_prediction(values: np.ndarray) -> np.ndarray:
    # Fit q(r)=a+b r^2+c r^4 through rings r/dr=1.5,2.5,3.5 and evaluate
    # at the first cell centre r/dr=0.5.
    return 1.8 * values[1] - values[2] + 0.2 * values[3]


def odd_axis_prediction(values: np.ndarray) -> np.ndarray:
    # Fit u_r/r as an even function through rings 1..3, then evaluate u_r
    # at r/dr=0.5.
    scaled = np.stack(
        (values[1] / 1.5, values[2] / 2.5, values[3] / 3.5), axis=0
    )
    return 0.5 * (1.8 * scaled[0] - scaled[1] + 0.2 * scaled[2])


def scalar_metrics(values: np.ndarray, floor: float):
    prediction = even_axis_prediction(values)
    residual = values[0] - prediction
    local_scale = np.maximum.reduce(
        [np.abs(values[0]), np.abs(values[1]), np.abs(prediction),
         np.full_like(prediction, floor)]
    )
    relative = np.abs(residual) / local_scale
    return prediction, residual, relative


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfiles", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    patch_yt_cylindrical_readonly_edge()

    rows = []
    for plotfile in sorted(args.plotfiles):
        dataset, level, fields = load_fields(plotfile)
        nr, nz = fields["pressure"].shape
        if nr < 4:
            raise RuntimeError(f"{plotfile}: need at least four radial rings")
        z_lo = float(dataset.domain_left_edge[1])
        z_hi = float(dataset.domain_right_edge[1])
        z = z_lo + (np.arange(nz) + 0.5) * (z_hi - z_lo) / nz

        p_pred, p_res, p_rel = scalar_metrics(fields["pressure"], 1.0)
        rho_pred, rho_res, rho_rel = scalar_metrics(fields["density"], 1.0e-12)
        t_pred, t_res, t_rel = scalar_metrics(fields["temperature"], 1.0e-12)
        uz_pred, uz_res, uz_rel = scalar_metrics(fields["axial_velocity"], 1.0)
        ur_pred = odd_axis_prediction(fields["radial_velocity"])
        ur_res = fields["radial_velocity"][0] - ur_pred
        ur_rel_sound = np.abs(ur_res) / np.maximum(fields["sound_speed"][0], 1.0)

        # Exclude the two axial boundary-adjacent cells from the maximisation.
        active = slice(2, nz - 2)
        j0 = 2 + int(np.argmax(p_rel[active]))
        ju = 2 + int(np.argmax(ur_rel_sound[active]))
        row = {
            "plotfile": str(plotfile),
            "time_s": f"{float(dataset.current_time):.17g}",
            "level": level,
            "nr": nr,
            "nz": nz,
            "p_axis_defect_max": f"{p_rel[j0]:.17g}",
            "p_axis_defect_z_m": f"{z[j0]:.17g}",
            "p_ring0_Pa": f"{fields['pressure'][0, j0]:.17g}",
            "p_even_prediction_Pa": f"{p_pred[j0]:.17g}",
            "p_residual_Pa": f"{p_res[j0]:.17g}",
            "p_ring1_Pa": f"{fields['pressure'][1, j0]:.17g}",
            "p_ring2_Pa": f"{fields['pressure'][2, j0]:.17g}",
            "p_ring3_Pa": f"{fields['pressure'][3, j0]:.17g}",
            "rho_axis_defect_max": f"{np.max(rho_rel[active]):.17g}",
            "temperature_axis_defect_max": f"{np.max(t_rel[active]):.17g}",
            "uz_axis_defect_max": f"{np.max(uz_rel[active]):.17g}",
            "ur_odd_defect_over_c_max": f"{ur_rel_sound[ju]:.17g}",
            "ur_odd_defect_z_m": f"{z[ju]:.17g}",
            "ur_ring0_m_s": f"{fields['radial_velocity'][0, ju]:.17g}",
            "ur_odd_prediction_m_s": f"{ur_pred[ju]:.17g}",
            "mach_ring0_max": f"{np.max(fields['mach'][0, active]):.17g}",
        }
        rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
