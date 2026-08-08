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
DUAL_TARGET_TOL = 2.0e-9
WALL_DERIVATIVE_TOL = 2.0e-8
CONSTRAINED_FACE_TOL = 1.0e-7


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


def dual_target_reproduction(samples: int, seed: int) -> float:
    """Verify face/reflected-face evaluation from one support polynomial.

    The face-local IBM stores two weight vectors over the same visible fluid
    support: one evaluates at the Cartesian flux-face centre, while the other
    evaluates at that point reflected through the owning wall plane.  This
    campaign exercises non-symmetric distances and randomly oriented planes.
    """
    rng = np.random.default_rng(seed)
    maximum_error = 0.0
    for dim in (2, 3):
        points = np.array(
            list(itertools.product((0.5, 1.5, 2.5), repeat=dim)),
            dtype=float,
        )
        mask = np.ones(points.shape[0], dtype=bool)
        global_basis = basis(points)
        quadratic_size = 6 if dim == 2 else 10
        for _ in range(samples):
            normal = rng.normal(size=dim)
            normal /= np.linalg.norm(normal)
            wall_point = rng.uniform(1.2, 1.8, size=dim)
            face_point = rng.uniform(0.9, 2.1, size=dim)
            signed_distance = float((face_point - wall_point) @ normal)
            reflected_point = face_point - 2.0 * signed_distance * normal

            face_order, face_weights = interpolation_weights(
                basis(points - face_point), mask, dim, minimum_points=dim
            )
            reflected_order, reflected_weights = interpolation_weights(
                basis(points - reflected_point),
                mask,
                dim,
                minimum_points=dim,
            )
            if face_order != 2 or reflected_order != 2:
                raise AssertionError(
                    "complete 3^D support unexpectedly lost quadratic rank"
                )

            coefficients = rng.normal(size=quadratic_size)
            support_values = global_basis @ coefficients
            face_exact = (basis(face_point[None, :]) @ coefficients).item()
            reflected_exact = (
                basis(reflected_point[None, :]) @ coefficients
            ).item()
            face_error = abs(float(face_weights @ support_values) - face_exact)
            reflected_error = abs(
                float(reflected_weights @ support_values) - reflected_exact
            )
            error = max(face_error, reflected_error)
            maximum_error = max(maximum_error, error)
            if not np.isfinite(error) or error > DUAL_TARGET_TOL:
                raise AssertionError(
                    f"dual-target reproduction failed: dim={dim} "
                    f"distance={signed_distance:.6e} error={error:.6e}"
                )
    return maximum_error


def polynomial_value_gradient(
    point: np.ndarray, coefficients: np.ndarray
) -> tuple[float, np.ndarray]:
    """Evaluate the C++ quadratic basis and its Cartesian gradient."""
    dim = point.size
    value = (basis(point[None, :]) @ coefficients).item()
    if dim == 2:
        x, y = point
        gradient = np.array(
            (
                coefficients[1] + 2.0 * coefficients[3] * x
                + coefficients[4] * y,
                coefficients[2] + coefficients[4] * x
                + 2.0 * coefficients[5] * y,
            )
        )
    else:
        x, y, z = point
        gradient = np.array(
            (
                coefficients[1] + 2.0 * coefficients[4] * x
                + coefficients[7] * y + coefficients[8] * z,
                coefficients[2] + 2.0 * coefficients[5] * y
                + coefficients[7] * x + coefficients[9] * z,
                coefficients[3] + 2.0 * coefficients[6] * z
                + coefficients[8] * x + coefficients[9] * y,
            )
        )
    return value, gradient


def wall_normal_derivative_reproduction(samples: int, seed: int) -> float:
    """Verify the direct wall derivative on non-grid-aligned normals.

    The direct closure reconstructs q(d) and q(2d) from one complete fluid
    WLS block, then applies (4q(d)-q(2d)-3q(0))/(2d).  A quadratic Cartesian
    polynomial restricts to a quadratic along every wall normal, so this
    composition must reproduce its wall-normal derivative to roundoff.
    """
    rng = np.random.default_rng(seed)
    maximum_error = 0.0
    for dim in (2, 3):
        support_points = np.array(
            list(itertools.product((0.5, 1.5, 2.5), repeat=dim)),
            dtype=float,
        )
        mask = np.ones(support_points.shape[0], dtype=bool)
        quadratic_size = 6 if dim == 2 else 10
        for _ in range(samples):
            normal = rng.normal(size=dim)
            normal /= np.linalg.norm(normal)
            wall_point = rng.uniform(0.8, 2.2, size=dim)
            distance = rng.uniform(0.45, 1.05)
            first_target = wall_point + distance * normal
            second_target = wall_point + 2.0 * distance * normal
            first_order, first_weights = interpolation_weights(
                basis(support_points - first_target), mask, dim,
                minimum_points=dim,
            )
            second_order, second_weights = interpolation_weights(
                basis(support_points - second_target), mask, dim,
                minimum_points=dim,
            )
            if first_order != 2 or second_order != 2:
                raise AssertionError(
                    "complete direct-wall support lost quadratic rank"
                )

            coefficients = rng.normal(size=quadratic_size)
            support_values = basis(support_points) @ coefficients
            wall_value, wall_gradient = polynomial_value_gradient(
                wall_point, coefficients
            )
            first_value = float(first_weights @ support_values)
            second_value = float(second_weights @ support_values)
            numerical = (
                4.0 * first_value - second_value - 3.0 * wall_value
            ) / (2.0 * distance)
            exact = float(wall_gradient @ normal)
            error = abs(numerical - exact)
            maximum_error = max(maximum_error, error)
            if not np.isfinite(error) or error > WALL_DERIVATIVE_TOL:
                raise AssertionError(
                    f"wall-normal derivative reproduction failed: dim={dim} "
                    f"distance={distance:.6e} error={error:.6e}"
                )
    return maximum_error


def constrained_dirichlet_basis(
    points: np.ndarray,
    wall_point: np.ndarray,
    normal: np.ndarray,
    spacing: np.ndarray,
    normal_scale: float,
) -> np.ndarray:
    delta = points - wall_point
    xi = delta / spacing
    normal_coordinate = (delta @ normal) / normal_scale
    if points.shape[1] == 2:
        x, y = xi.T
        return np.column_stack(
            (
                normal_coordinate,
                x * x,
                x * y,
                y * y,
                x**3,
                x * x * y,
                x * y * y,
                y**3,
            )
        )
    x, y, z = xi.T
    return np.column_stack(
        (
            normal_coordinate,
            x * x,
            y * y,
            z * z,
            x * y,
            x * z,
            y * z,
            x**3,
            y**3,
            z**3,
            x * x * y,
            x * x * z,
            x * y * y,
            y * y * z,
            x * z * z,
            y * z * z,
            x * y * z,
        )
    )


def constrained_monomial_exponents(dim: int) -> np.ndarray:
    """Return exponents for the degree-two/three constrained basis terms."""
    if dim == 2:
        # Entry zero is the special wall-normal linear coordinate.
        return np.asarray(
            (
                (0, 0),
                (2, 0), (1, 1), (0, 2),
                (3, 0), (2, 1), (1, 2), (0, 3),
            ),
            dtype=int,
        )
    return np.asarray(
        (
            (0, 0, 0),
            (2, 0, 0), (0, 2, 0), (0, 0, 2),
            (1, 1, 0), (1, 0, 1), (0, 1, 1),
            (3, 0, 0), (0, 3, 0), (0, 0, 3),
            (2, 1, 0), (2, 0, 1), (1, 2, 0),
            (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1),
        ),
        dtype=int,
    )


def finite_volume_target_functionals(
    point_value: np.ndarray,
    point_gradient: np.ndarray,
    face_xi: np.ndarray,
    spacing: np.ndarray,
    face_direction: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the exact fourth-order FV modified target to a cubic basis."""
    value = point_value.copy()
    gradient = point_gradient.copy()
    exponents = constrained_monomial_exponents(face_xi.size)
    h_normal_sq = spacing[face_direction] ** 2
    for basis_index in range(1, exponents.shape[0]):
        exponent = exponents[basis_index]
        if exponent[face_direction] < 2:
            continue

        derivative_counts = np.zeros(face_xi.size, dtype=int)
        derivative_counts[face_direction] = 2

        def derivative(counts: np.ndarray) -> float:
            result = 1.0
            for axis, count in enumerate(counts):
                power = int(exponent[axis])
                if power < count:
                    return 0.0
                for derivative_index in range(count):
                    result *= (power - derivative_index) / spacing[axis]
                result *= face_xi[axis] ** (power - count)
            return result

        value[basis_index] -= (
            h_normal_sq * derivative(derivative_counts) / 24.0
        )
        for gradient_direction in range(face_xi.size):
            mixed_counts = derivative_counts.copy()
            mixed_counts[gradient_direction] += 1
            gradient[basis_index, gradient_direction] -= (
                h_normal_sq * derivative(mixed_counts) / 24.0
            )
    return value, gradient


def discrete_fourth_order_face_targets(
    wall_point: np.ndarray,
    face_point: np.ndarray,
    normal: np.ndarray,
    spacing: np.ndarray,
    normal_scale: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the independent ``diff_ops.H`` four-point face operators."""
    dim = wall_point.size
    direction = int(np.argmax(np.abs(face_point - wall_point)))
    normal_offsets = np.asarray((-1.5, -0.5, 0.5, 1.5))
    interpolation = np.asarray((-1.0, 7.0, 7.0, -1.0)) / 12.0
    normal_derivative = np.asarray((1.0, -15.0, 15.0, -1.0)) / 12.0
    tangent_offsets = np.asarray((-2.0, -1.0, 1.0, 2.0))
    tangent_derivative = np.asarray((1.0, -8.0, 8.0, -1.0)) / 12.0

    row_points = np.repeat(face_point[None, :], 4, axis=0)
    row_points[:, direction] += normal_offsets * spacing[direction]
    row_basis = constrained_dirichlet_basis(
        row_points, wall_point, normal, spacing, normal_scale
    )
    value_target = interpolation @ row_basis
    gradient_target = np.empty((row_basis.shape[1], dim))
    gradient_target[:, direction] = (
        normal_derivative @ row_basis / spacing[direction]
    )

    for tangent_direction in range(dim):
        if tangent_direction == direction:
            continue
        row_derivatives = np.empty_like(row_basis)
        for row, row_point in enumerate(row_points):
            tangent_points = np.repeat(row_point[None, :], 4, axis=0)
            tangent_points[:, tangent_direction] += (
                tangent_offsets * spacing[tangent_direction]
            )
            tangent_basis = constrained_dirichlet_basis(
                tangent_points, wall_point, normal, spacing, normal_scale
            )
            row_derivatives[row] = (
                tangent_derivative @ tangent_basis
                / spacing[tangent_direction]
            )
        gradient_target[:, tangent_direction] = (
            interpolation @ row_derivatives
        )
    return value_target, gradient_target


def constrained_dirichlet_face_reproduction(
    samples: int, seed: int
) -> tuple[float, float]:
    """Verify the FV value/gradient contract at a Cartesian crossing face."""
    rng = np.random.default_rng(seed)
    maximum_value_error = 0.0
    maximum_gradient_error = 0.0
    for dim in (2, 3):
        offsets = np.array(
            list(itertools.product((0.5, 1.5, 2.5), repeat=dim)),
            dtype=float,
        )
        for _ in range(samples):
            spacing = rng.uniform(0.7, 1.3, size=dim)
            support_points = offsets * spacing
            normal = rng.normal(size=dim)
            normal /= np.linalg.norm(normal)
            wall_point = rng.uniform(0.8, 2.2, size=dim) * spacing
            normal_scale = rng.uniform(0.55, 1.1) * np.min(spacing)
            direction = int(rng.integers(0, dim))
            face_offset = rng.uniform(-0.7, 0.7) * spacing[direction]
            face_point = wall_point.copy()
            face_point[direction] += face_offset

            support_basis = constrained_dirichlet_basis(
                support_points, wall_point, normal, spacing, normal_scale
            )
            face_basis = constrained_dirichlet_basis(
                face_point[None, :], wall_point, normal, spacing,
                normal_scale,
            )[0]
            coefficients = rng.normal(size=support_basis.shape[1])
            wall_value = rng.uniform(250.0, 450.0)
            support_values = wall_value + support_basis @ coefficients

            normal_matrix = support_basis.T @ support_basis
            if np.linalg.cond(normal_matrix, p=1) > 1.0e8:
                continue
            # The first basis gradient is n/d. Polynomial gradients use the
            # same anisotropic Cartesian scaling as the C++ implementation.
            delta = face_point - wall_point
            xi = delta / spacing
            gradient_basis = np.zeros((support_basis.shape[1], dim))
            gradient_basis[0] = normal / normal_scale
            if dim == 2:
                gradient_basis[1, 0] = 2.0 * xi[0] / spacing[0]
                gradient_basis[2, 0] = xi[1] / spacing[0]
                gradient_basis[2, 1] = xi[0] / spacing[1]
                gradient_basis[3, 1] = 2.0 * xi[1] / spacing[1]
                gradient_basis[4, 0] = 3.0 * xi[0] ** 2 / spacing[0]
                gradient_basis[5, 0] = (
                    2.0 * xi[0] * xi[1] / spacing[0]
                )
                gradient_basis[5, 1] = xi[0] ** 2 / spacing[1]
                gradient_basis[6, 0] = xi[1] ** 2 / spacing[0]
                gradient_basis[6, 1] = (
                    2.0 * xi[0] * xi[1] / spacing[1]
                )
                gradient_basis[7, 1] = 3.0 * xi[1] ** 2 / spacing[1]
            else:
                gradient_basis[1, 0] = 2.0 * xi[0] / spacing[0]
                gradient_basis[2, 1] = 2.0 * xi[1] / spacing[1]
                gradient_basis[3, 2] = 2.0 * xi[2] / spacing[2]
                gradient_basis[4, 0] = xi[1] / spacing[0]
                gradient_basis[4, 1] = xi[0] / spacing[1]
                gradient_basis[5, 0] = xi[2] / spacing[0]
                gradient_basis[5, 2] = xi[0] / spacing[2]
                gradient_basis[6, 1] = xi[2] / spacing[1]
                gradient_basis[6, 2] = xi[1] / spacing[2]
                gradient_basis[7, 0] = 3.0 * xi[0] ** 2 / spacing[0]
                gradient_basis[8, 1] = 3.0 * xi[1] ** 2 / spacing[1]
                gradient_basis[9, 2] = 3.0 * xi[2] ** 2 / spacing[2]
                gradient_basis[10, 0] = (
                    2.0 * xi[0] * xi[1] / spacing[0]
                )
                gradient_basis[10, 1] = xi[0] ** 2 / spacing[1]
                gradient_basis[11, 0] = (
                    2.0 * xi[0] * xi[2] / spacing[0]
                )
                gradient_basis[11, 2] = xi[0] ** 2 / spacing[2]
                gradient_basis[12, 0] = xi[1] ** 2 / spacing[0]
                gradient_basis[12, 1] = (
                    2.0 * xi[0] * xi[1] / spacing[1]
                )
                gradient_basis[13, 1] = (
                    2.0 * xi[1] * xi[2] / spacing[1]
                )
                gradient_basis[13, 2] = xi[1] ** 2 / spacing[2]
                gradient_basis[14, 0] = xi[2] ** 2 / spacing[0]
                gradient_basis[14, 2] = (
                    2.0 * xi[0] * xi[2] / spacing[2]
                )
                gradient_basis[15, 1] = xi[2] ** 2 / spacing[1]
                gradient_basis[15, 2] = (
                    2.0 * xi[1] * xi[2] / spacing[2]
                )
                gradient_basis[16, 0] = xi[1] * xi[2] / spacing[0]
                gradient_basis[16, 1] = xi[0] * xi[2] / spacing[1]
                gradient_basis[16, 2] = xi[0] * xi[1] / spacing[2]

            fv_value_basis, fv_gradient_basis = (
                finite_volume_target_functionals(
                    face_basis, gradient_basis, xi, spacing, direction
                )
            )
            discrete_value_basis, discrete_gradient_basis = (
                discrete_fourth_order_face_targets(
                    wall_point, face_point, normal, spacing, normal_scale
                )
            )
            contract_error = max(
                float(np.max(np.abs(fv_value_basis - discrete_value_basis))),
                float(
                    np.max(
                        np.abs(fv_gradient_basis - discrete_gradient_basis)
                    )
                ),
            )
            if contract_error > 2.0e-12:
                raise AssertionError(
                    "finite-volume target contract failed: "
                    f"dim={dim} error={contract_error:.6e}"
                )

            value_weights = support_basis @ np.linalg.solve(
                normal_matrix, fv_value_basis
            )
            reconstructed_value = wall_value + value_weights @ (
                support_values - wall_value
            )
            exact_value = wall_value + discrete_value_basis @ coefficients
            value_error = abs(float(reconstructed_value - exact_value))

            reconstructed_gradient = np.empty(dim)
            for direction_index in range(dim):
                gradient_weights = support_basis @ np.linalg.solve(
                    normal_matrix, fv_gradient_basis[:, direction_index]
                )
                if normal_scale * np.sum(np.abs(gradient_weights)) > 32.0:
                    reconstructed_gradient = None
                    break
                reconstructed_gradient[direction_index] = (
                    gradient_weights @ (support_values - wall_value)
                )
            if np.sum(np.abs(value_weights)) > 32.0 or (
                reconstructed_gradient is None
            ):
                continue
            exact_gradient = discrete_gradient_basis.T @ coefficients
            gradient_error = float(
                np.max(np.abs(reconstructed_gradient - exact_gradient))
            )
            maximum_value_error = max(maximum_value_error, value_error)
            maximum_gradient_error = max(
                maximum_gradient_error, gradient_error
            )
            if value_error > CONSTRAINED_FACE_TOL or (
                gradient_error > CONSTRAINED_FACE_TOL
            ):
                raise AssertionError(
                    "constrained face reproduction failed: "
                    f"dim={dim} value={value_error:.6e} "
                    f"gradient={gradient_error:.6e}"
                )
    return maximum_value_error, maximum_gradient_error


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
    dual_error = dual_target_reproduction(
        max(1000, args.samples_3d // 2), args.seed + 1
    )
    wall_derivative_error = wall_normal_derivative_reproduction(
        max(1000, args.samples_3d // 2), args.seed + 2
    )
    constrained_value_error, constrained_gradient_error = (
        constrained_dirichlet_face_reproduction(
            max(1000, args.samples_3d // 2), args.seed + 3
        )
    )
    print_stats("2-D exhaustive", stats_2d)
    print_stats("3-D sampled", stats_3d)
    print(f"dual-target face/reflection max_error={dual_error:.3e}")
    print(
        "direct wall-normal quadratic derivative "
        f"max_error={wall_derivative_error:.3e}"
    )
    print(
        "cubic constrained Cartesian-face value/gradient max_error="
        f"{constrained_value_error:.3e}/{constrained_gradient_error:.3e}"
    )
    print("IBM WLS polynomial reproduction: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
