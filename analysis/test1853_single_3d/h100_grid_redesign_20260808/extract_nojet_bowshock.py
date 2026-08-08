#!/usr/bin/env python3
"""Extract bow-shock locations from streamwise rays in a no-jet plotfile."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt


BODY_NOSE_X_M = 0.00244127
RAYS_M = {
    "axis": (0.0, 0.0),
    "y_plus_30mm": (0.030, 0.0),
    "y_minus_30mm": (-0.030, 0.0),
    "z_plus_30mm": (0.0, 0.030),
    "z_minus_30mm": (0.0, -0.030),
}


def named_field(dataset, name: str):
    matches = [field for field in dataset.field_list if field[1] == name]
    if len(matches) != 1:
        raise RuntimeError(f"expected one field named {name!r}, found {matches}")
    return matches[0]


def ray_shock(dataset, y: float, z: float, fields: dict[str, tuple]) -> dict:
    ray = dataset.ortho_ray(0, (y, z))
    x = np.asarray(ray["x"], dtype=np.float64)
    rho = np.asarray(ray[fields["rho"]], dtype=np.float64)
    pressure = np.asarray(ray[fields["pressure"]], dtype=np.float64)
    temperature = np.asarray(ray[fields["temperature"]], dtype=np.float64)
    velocity = np.asarray(ray[fields["velocity"]], dtype=np.float64)
    sld = np.asarray(ray[fields["sld"]], dtype=np.float64)
    ghs = np.asarray(ray[fields["ghs"]], dtype=np.float64)

    selected = (
        np.isfinite(x)
        & np.isfinite(rho)
        & (x >= -0.150)
        & (x <= -0.002)
        & (sld < 0.5)
        & (ghs < 0.5)
    )
    order = np.argsort(x[selected])
    x = x[selected][order]
    rho = rho[selected][order]
    pressure = pressure[selected][order]
    temperature = temperature[selected][order]
    velocity = velocity[selected][order]
    if x.size < 12:
        raise RuntimeError(f"insufficient fluid samples for ray y={y}, z={z}")

    gradient = np.gradient(rho, x)
    shock_index = int(np.argmax(np.abs(gradient)))
    if shock_index < 4 or shock_index + 4 >= x.size:
        raise RuntimeError(f"shock detector reached search boundary for y={y}, z={z}")

    def state(lo: int, hi: int) -> dict[str, float]:
        return {
            "density_kg_m3": float(np.mean(rho[lo:hi])),
            "pressure_Pa": float(np.mean(pressure[lo:hi])),
            "temperature_K": float(np.mean(temperature[lo:hi])),
            "x_velocity_m_s": float(np.mean(velocity[lo:hi])),
        }

    return {
        "y_m": y,
        "z_m": z,
        "shock_x_m": float(x[shock_index]),
        "shock_density_gradient_kg_m4": float(gradient[shock_index]),
        "pre_shock_left_state": state(shock_index - 4, shock_index - 1),
        "post_shock_right_state": state(shock_index + 2, shock_index + 5),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    dataset = yt.load(str(args.plotfile))
    fields = {
        "rho": named_field(dataset, "Density"),
        "pressure": named_field(dataset, "pressure"),
        "temperature": named_field(dataset, "temperature"),
        "velocity": named_field(dataset, "x_velocity"),
        "sld": named_field(dataset, "sld"),
        "ghs": named_field(dataset, "ghs"),
    }
    rays = {
        name: ray_shock(dataset, y, z, fields)
        for name, (y, z) in RAYS_M.items()
    }
    axis_x = rays["axis"]["shock_x_m"]
    result = {
        "plotfile": str(args.plotfile),
        "time_s": float(dataset.current_time),
        "method": "maximum absolute streamwise density gradient in [-0.150,-0.002] m",
        "body_nose_x_m": BODY_NOSE_X_M,
        "axis_bowshock_x_m": axis_x,
        "axis_standoff_m": BODY_NOSE_X_M - axis_x,
        "rays": rays,
        "symmetry": {
            "y_pair_shock_x_difference_m": abs(
                rays["y_plus_30mm"]["shock_x_m"]
                - rays["y_minus_30mm"]["shock_x_m"]
            ),
            "z_pair_shock_x_difference_m": abs(
                rays["z_plus_30mm"]["shock_x_m"]
                - rays["z_minus_30mm"]["shock_x_m"]
            ),
            "mean_y_vs_mean_z_shock_x_difference_m": abs(
                0.5
                * (
                    rays["y_plus_30mm"]["shock_x_m"]
                    + rays["y_minus_30mm"]["shock_x_m"]
                )
                - 0.5
                * (
                    rays["z_plus_30mm"]["shock_x_m"]
                    + rays["z_minus_30mm"]["shock_x_m"]
                )
            ),
        },
    }
    rendered = json.dumps(result, indent=2)
    print(rendered)
    args.output.write_text(rendered + "\n", encoding="ascii")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
