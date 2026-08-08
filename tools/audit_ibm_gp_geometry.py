#!/usr/bin/env python3
"""Audit the geometric and polynomial contracts of IBM ghost-point CSV data."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def vector(row: dict[str, str], prefix: str) -> np.ndarray:
    return np.array([float(row[f"{prefix}_x"]), float(row[f"{prefix}_y"])])


def maximum(values: list[float]) -> float:
    return max(values) if values else float("nan")


def percentile(values: list[float], q: float) -> float:
    return float(np.percentile(np.asarray(values), q)) if values else float("nan")


def audit_geometry(
    rows: list[dict[str, str]], center: np.ndarray, radius: float, spacing: float
) -> None:
    normal_norm_error: list[float] = []
    gp_normal_error: list[float] = []
    gp_normal_angle: list[float] = []
    radial_error: list[float] = []
    radial_normal_angle: list[float] = []
    image_errors: dict[int, list[float]] = defaultdict(list)
    support_counts: dict[int, Counter[int]] = defaultdict(Counter)
    weight_sum_errors: dict[int, list[float]] = defaultdict(list)
    weight_l1: dict[int, list[float]] = defaultdict(list)

    image_count = sum(key.startswith("disIM") for key in rows[0])
    for row in rows:
        gp = np.array([float(row["x"]), float(row["y"])])
        boundary = vector(row, "ib")
        normal = vector(row, "normal")
        distance = float(row["disGP"])

        normal_norm_error.append(abs(float(np.linalg.norm(normal)) - 1.0))
        residual = gp - (boundary - distance * normal)
        gp_normal_error.append(float(np.linalg.norm(residual)) / spacing)

        gp_direction = boundary - gp
        gp_direction_norm = float(np.linalg.norm(gp_direction))
        if gp_direction_norm > 0.0:
            cosine = float(np.dot(gp_direction, normal) / gp_direction_norm)
            gp_normal_angle.append(
                math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
            )

        radial = boundary - center
        radial_norm = float(np.linalg.norm(radial))
        radial_error.append(abs(radial_norm - radius) / spacing)
        if radial_norm > 0.0:
            cosine = float(np.dot(radial / radial_norm, normal))
            radial_normal_angle.append(
                math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
            )

        for image in range(image_count):
            image_point = np.array(
                [float(row[f"ip{image}_x"]), float(row[f"ip{image}_y"])]
            )
            image_distance = float(row[f"disIM{image}"])
            expected = boundary + image_distance * normal
            image_errors[image].append(float(np.linalg.norm(image_point - expected)) / spacing)
            support_counts[image][int(row[f"ninterp{image}"])] += 1
            weight_sum_errors[image].append(abs(float(row[f"weight_sum{image}"]) - 1.0))
            weight_l1[image].append(float(row[f"weight_l1{image}"]))

    print(f"GP count: {len(rows)}")
    print(f"max |norm(n)-1|: {maximum(normal_norm_error):.3e}")
    print(
        "GP-BI normal-line residual / h: "
        f"max={maximum(gp_normal_error):.3e}, "
        f"p99={percentile(gp_normal_error, 99.0):.3e}"
    )
    print(
        "angle(BI-GP,n) [deg]: "
        f"max={maximum(gp_normal_angle):.6g}, "
        f"p99={percentile(gp_normal_angle, 99.0):.6g}"
    )
    print(
        "boundary-intercept radial error / h: "
        f"max={maximum(radial_error):.3e}, "
        f"p99={percentile(radial_error, 99.0):.3e}"
    )
    print(
        "angle(n, exact radial normal) [deg]: "
        f"max={maximum(radial_normal_angle):.6g}, "
        f"p99={percentile(radial_normal_angle, 99.0):.6g}"
    )
    for image in range(image_count):
        print(
            f"IP{image} normal-line residual / h: "
            f"max={maximum(image_errors[image]):.3e}; "
            f"ninterp={dict(sorted(support_counts[image].items()))}; "
            f"max |sum(w)-1|={maximum(weight_sum_errors[image]):.3e}; "
            f"max sum|w|={maximum(weight_l1[image]):.6g}"
        )

    selected = {
        "upstream stagnation": min(
            rows, key=lambda row: (float(row["ib_x"]), abs(float(row["ib_y"])))
        ),
        "downstream stagnation": max(
            rows, key=lambda row: (float(row["ib_x"]), -abs(float(row["ib_y"])))
        ),
        "upper symmetry": max(
            rows, key=lambda row: (float(row["ib_y"]), -abs(float(row["ib_x"])))
        ),
    }
    for label, row in selected.items():
        print(
            f"{label}: gp={row['gp_index']} cell=({row['i']},{row['j']}), "
            f"GP=({float(row['x']):.9g},{float(row['y']):.9g}), "
            f"BI=({float(row['ib_x']):.9g},{float(row['ib_y']):.9g}), "
            f"n=({float(row['normal_x']):.9g},{float(row['normal_y']):.9g}), "
            f"dGP/h={float(row['disGP']) / spacing:.6g}"
        )


def polynomial_basis(coordinate: np.ndarray, order: int) -> np.ndarray:
    x, y = coordinate
    if order == 0:
        return np.array([1.0])
    if order == 1:
        return np.array([1.0, x, y])
    return np.array([1.0, x, y, x * x, x * y, y * y])


def audit_moments(
    rows: list[dict[str, str]], prob_lo: np.ndarray, spacing: float
) -> None:
    groups: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(int(row["target"]), int(row["image"]))].append(row)

    errors: dict[int, list[float]] = defaultdict(list)
    active_weights: dict[int, list[int]] = defaultdict(list)
    for group in groups.values():
        order = int(group[0]["applied_fit_order"])
        if order < 0:
            continue
        image = np.array([float(group[0]["image_x"]), float(group[0]["image_y"])])
        moment = np.zeros({0: 1, 1: 3, 2: 6}[order])
        nonzero = 0
        for row in group:
            weight = float(row["interpolation_weight"])
            if weight != 0.0:
                nonzero += 1
            support = prob_lo + spacing * np.array(
                [int(row["support_i"]) + 0.5, int(row["support_j"]) + 0.5]
            )
            coordinate = (support - image) / spacing
            moment += weight * polynomial_basis(coordinate, order)
        target = np.zeros_like(moment)
        target[0] = 1.0
        errors[order].append(float(np.max(np.abs(moment - target))))
        active_weights[order].append(nonzero)

    print(f"WLS image-point groups: {len(groups)}")
    for order in sorted(errors):
        print(
            f"order {order} discrete-moment error: count={len(errors[order])}, "
            f"max={maximum(errors[order]):.3e}, "
            f"p99={percentile(errors[order], 99.0):.3e}, "
            f"nonzero weights={dict(sorted(Counter(active_weights[order]).items()))}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gp_csv", type=Path)
    parser.add_argument("support_csv", type=Path)
    parser.add_argument("--prob-lo", type=float, nargs=2, default=(-0.5, -0.5))
    parser.add_argument("--spacing", type=float, required=True)
    parser.add_argument("--center", type=float, nargs=2, default=(0.0, 0.0))
    parser.add_argument("--radius", type=float, default=0.1)
    args = parser.parse_args()

    gp_rows = read_csv(args.gp_csv)
    support_rows = read_csv(args.support_csv)
    if not gp_rows:
        raise SystemExit("GP CSV is empty")
    audit_geometry(
        gp_rows, np.asarray(args.center, dtype=float), args.radius, args.spacing
    )
    audit_moments(
        support_rows, np.asarray(args.prob_lo, dtype=float), args.spacing
    )


if __name__ == "__main__":
    main()
