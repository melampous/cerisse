#!/usr/bin/env python3
"""Verify the polynomial contract of the quadratic IBM image-point WLS.

The implementation mirrors ``ibm_spd_solve_e0`` and the quadratic -> linear
-> constant demotion ladder in ``src/ibm/ibm_solver_interp.h``.  Two-dimensional
3x3 masks are exhaustive; three-dimensional 3x3x3 masks combine deterministic
degenerate patterns with a seeded random campaign.
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import dataclass

import numpy as np


PIVOT_REL_TOL = 1.0e-11
MOMENT_TOL = 2.0e-10


@dataclass
class TestStats:
    masks: int = 0
    invalid: int = 0
    constant: int = 0
    linear: int = 0
    quadratic: int = 0
    max_moment_error: float = 0.0

    def add(self, order: int, error: float) -> None:
        self.masks += 1
        self.max_moment_error = max(self.max_moment_error, error)
        if order < 0:
            self.invalid += 1
        elif order == 0:
            self.constant += 1
        elif order == 1:
            self.linear += 1
        elif order == 2:
            self.quadratic += 1
        else:
            raise AssertionError(f"unexpected fit order {order}")


def basis(points: np.ndarray) -> np.ndarray:
    """Return the exact monomial ordering used by the C++ IBM WLS."""
    dim = points.shape[1]
    one = np.ones(points.shape[0])
    if dim == 2:
        x, y = points.T
        return np.column_stack((one, x, y, x * x, x * y, y * y))
    if dim == 3:
        x, y, z = points.T
        return np.column_stack(
            (one, x, y, z, x * x, y * y, z * z, x * y, x * z, y * z)
        )
    raise ValueError(f"unsupported dimension {dim}")


def solve_e0(matrix: np.ndarray) -> np.ndarray | None:
    """Mirror ``ibm_spd_solve_e0`` including its relative pivot guard."""
    size = matrix.shape[0]
    tolerance = PIVOT_REL_TOL * max(float(np.max(np.diag(matrix))), 1.0e-300)
    factor = np.zeros_like(matrix)
    for row in range(size):
        diagonal = matrix[row, row] - np.dot(
            factor[row, :row], factor[row, :row]
        )
        if not np.isfinite(diagonal) or not diagonal > tolerance:
            return None
        factor[row, row] = np.sqrt(diagonal)
        for lower_row in range(row + 1, size):
            value = matrix[lower_row, row] - np.dot(
                factor[lower_row, :row], factor[row, :row]
            )
            factor[lower_row, row] = value / factor[row, row]

    work = np.zeros(size)
    for row in range(size):
        value = (1.0 if row == 0 else 0.0) - np.dot(
            factor[row, :row], work[:row]
        )
        work[row] = value / factor[row, row]

    solution = np.zeros(size)
    for row in range(size - 1, -1, -1):
        value = work[row] - np.dot(
            factor[row + 1 :, row], solution[row + 1 :]
        )
        solution[row] = value / factor[row, row]
    return solution


def interpolation_weights(
    all_basis: np.ndarray, mask: np.ndarray, dim: int, minimum_points: int
) -> tuple[int, np.ndarray]:
    """Apply the C++ fit ladder and return ``(fit_order, weights)``."""
    point_count = int(np.count_nonzero(mask))
    weights = np.zeros(mask.size)
    if point_count < minimum_points:
        return -1, weights

    active_basis = all_basis[mask]
    normal_matrix = active_basis.T @ active_basis
    quadratic_size = 6 if dim == 2 else 10
    linear_size = dim + 1

    used_size = 0
    solution: np.ndarray | None = None
    if point_count >= quadratic_size:
        solution = solve_e0(normal_matrix)
        if solution is not None:
            used_size = quadratic_size
    if used_size == 0 and point_count >= linear_size:
        solution = solve_e0(normal_matrix[:linear_size, :linear_size])
        if solution is not None:
            used_size = linear_size

    if used_size == 0:
        weights[mask] = 1.0 / point_count
        return 0, weights

    assert solution is not None
    weights[mask] = active_basis[:, :used_size] @ solution
    return (2 if used_size == quadratic_size else 1), weights


def verify_mask(
    all_basis: np.ndarray,
    mask: np.ndarray,
    dim: int,
    minimum_points: int,
) -> tuple[int, float]:
    order, weights = interpolation_weights(
        all_basis, mask, dim, minimum_points
    )
    if order < 0:
        if np.any(weights):
            raise AssertionError("invalid stencil has nonzero weights")
        return order, 0.0
    if np.any(weights[~mask]):
        raise AssertionError("masked support has nonzero weight")

    used_size = {0: 1, 1: dim + 1, 2: (6 if dim == 2 else 10)}[order]
    target = np.zeros(used_size)
    target[0] = 1.0
    moments = all_basis[:, :used_size].T @ weights
    error = float(np.max(np.abs(moments - target)))
    if not np.isfinite(error) or error > MOMENT_TOL:
        raise AssertionError(
            f"dim={dim} order={order} points={mask.sum()} "
            f"moment_error={error:.6e} moments={moments}"
        )

    # A second, independent statement of the same contract: every polynomial
    # in the accepted basis evaluates at the image point to its constant term.
    coefficients = np.linspace(-0.7, 1.1, used_size)
    values = all_basis[:, :used_size] @ coefficients
    value_error = abs(float(weights @ values) - float(coefficients[0]))
    if value_error > 5.0 * MOMENT_TOL:
        raise AssertionError(
            f"polynomial reproduction failed: dim={dim} order={order} "
            f"error={value_error:.6e}"
        )
    return order, max(error, value_error)


def grid_basis(dim: int, image_coordinate: tuple[float, ...]) -> np.ndarray:
    # Cell centres in the candidate block are at 0.5, 1.5, and 2.5 in units
    # of grid spacing.  A centred image point lies in [1,2)^D relative to the
    # block anchor selected by ibm_interp_base_index.
    points = np.array(list(itertools.product((0.5, 1.5, 2.5), repeat=dim)))
    return basis(points - np.asarray(image_coordinate))


def exhaustive_2d() -> TestStats:
    stats = TestStats()
    image_coordinates = (
        (1.01, 1.01),
        (1.23, 1.77),
        (1.50, 1.50),
        (1.99, 1.99),
    )
    for image_coordinate in image_coordinates:
        all_basis = grid_basis(2, image_coordinate)
        for bits in range(1 << 9):
            mask = np.array([(bits >> point) & 1 for point in range(9)], dtype=bool)
            order, error = verify_mask(all_basis, mask, 2, minimum_points=2)
            stats.add(order, error)
    return stats


def structured_3d_masks() -> list[np.ndarray]:
    offsets = np.array(list(itertools.product(range(3), repeat=3)))
    masks = [np.ones(27, dtype=bool)]
    for direction in range(3):
        for index in range(3):
            masks.append(offsets[:, direction] == index)
            masks.append(offsets[:, direction] <= index)
            masks.append(offsets[:, direction] >= index)
    masks.extend(
        [
            np.count_nonzero(offsets == 1, axis=1) >= 2,
            np.sum(offsets, axis=1) <= 2,
            np.sum(offsets, axis=1) >= 4,
            np.all(offsets != 1, axis=1),
        ]
    )
    return masks


def sampled_3d(samples: int, seed: int) -> TestStats:
    stats = TestStats()
    rng = np.random.default_rng(seed)
    image_coordinates = (
        (1.01, 1.01, 1.01),
        (1.19, 1.51, 1.83),
        (1.50, 1.50, 1.50),
        (1.99, 1.99, 1.99),
    )
    probabilities = (0.12, 0.2, 0.35, 0.5, 0.7, 0.9)
    for image_coordinate in image_coordinates:
        all_basis = grid_basis(3, image_coordinate)
        for mask in structured_3d_masks():
            order, error = verify_mask(all_basis, mask, 3, minimum_points=3)
            stats.add(order, error)
        for sample in range(samples):
            probability = probabilities[sample % len(probabilities)]
            mask = rng.random(27) < probability
            order, error = verify_mask(all_basis, mask, 3, minimum_points=3)
            stats.add(order, error)
    return stats


def print_stats(name: str, stats: TestStats) -> None:
    print(
        f"{name}: masks={stats.masks} invalid={stats.invalid} "
        f"constant={stats.constant} linear={stats.linear} "
        f"quadratic={stats.quadratic} "
        f"max_moment_error={stats.max_moment_error:.3e}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-3d", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260713)
    args = parser.parse_args()
    if args.samples_3d < 1:
        parser.error("--samples-3d must be positive")

    stats_2d = exhaustive_2d()
    stats_3d = sampled_3d(args.samples_3d, args.seed)
    print_stats("2-D exhaustive", stats_2d)
    print_stats("3-D sampled", stats_3d)
    print("IBM WLS polynomial reproduction: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
