#!/usr/bin/env python3
"""Extract surface-recovery fields nearest the Run165 low-temperature cell."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


TARGET_M = np.array([0.00290625, 0.005181818182, -0.007363636364])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtp", nargs="+", type=Path)
    parser.add_argument("--label", action="append", default=[])
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def data_array(text: str, name: str, dtype=float) -> np.ndarray:
    pattern = (
        r'<DataArray\b[^>]*\bName="'
        + re.escape(name)
        + r'"[^>]*>(.*?)</DataArray>'
    )
    match = re.search(pattern, text, flags=re.DOTALL)
    if match is None:
        raise RuntimeError(f"missing VTP DataArray {name!r}")
    return np.fromstring(match.group(1), sep=" ", dtype=dtype)


def extract(path: Path, label: str) -> dict:
    text = path.read_text(encoding="ascii")
    points = data_array(text, "Points").reshape((-1, 3))
    connectivity = data_array(text, "connectivity", np.int64).reshape((-1, 3))
    triangles = points[connectivity]
    centroids = np.mean(triangles, axis=1)
    nearest_face = int(np.argmin(np.linalg.norm(centroids - TARGET_M, axis=1)))
    triangle = triangles[nearest_face]
    normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
    normal /= np.linalg.norm(normal)
    result = {
        "label": label,
        "file": str(path),
        "target_m": TARGET_M.tolist(),
        "face_index": nearest_face,
        "centroid_m": centroids[nearest_face].tolist(),
        "centroid_distance_m": float(
            np.linalg.norm(centroids[nearest_face] - TARGET_M)
        ),
        "normal": normal.tolist(),
    }
    for name in (
        "Pressure",
        "Temperature",
        "TauNormal",
        "Tau1",
        "Tau2",
        "dTdn",
        "PressureShockSensor",
        "PressureFallbackFraction",
        "Level",
        "IP_quality",
        "IP_fit_order",
    ):
        dtype = np.int64 if name in {"Level", "IP_quality", "IP_fit_order"} else float
        values = data_array(text, name, dtype)
        result[name] = (
            int(values[nearest_face])
            if np.issubdtype(values.dtype, np.integer)
            else float(values[nearest_face])
        )
    return result


def main() -> None:
    args = parse_args()
    if args.label and len(args.label) != len(args.vtp):
        raise SystemExit("--label must be supplied once per VTP or omitted")
    labels = args.label or [path.stem for path in args.vtp]
    output = {
        "analysis": "nearest_surface_recovery_history",
        "results": [extract(path, label) for path, label in zip(args.vtp, labels)],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
