#!/usr/bin/env python3
"""Polynomial oracles for the finite-volume BI and entropy-jet functionals.

The tests use exact Cartesian cell averages of local quadratic polynomials for
the Dirichlet functional and point samples of a smooth discrete auxiliary field
for the one-sided jet.  Random wall orientation, anisotropic cells, wall phase,
and ghost distance exercise both two- and three-dimensional formulas.
"""

from __future__ import annotations

import argparse
import itertools

import numpy as np


def orthonormal_frame(dim: int, rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(dim, dim))
    q, _ = np.linalg.qr(matrix)
    if np.linalg.det(q) < 0.0:
        q[:, -1] *= -1.0
    # Rows are (normal, tangent_1[, tangent_2]), matching the C++ storage.
    return q.T


def coordinates(points: np.ndarray, wall: np.ndarray, frame: np.ndarray,
                scale: float) -> np.ndarray:
    return (points - wall) @ frame.T / scale


def full_basis(local: np.ndarray) -> np.ndarray:
    if local.shape[1] == 2:
        s, t = local.T
        return np.column_stack((np.ones(local.shape[0]), s, t,
                                s * s, s * t, t * t))
    s, t, r = local.T
    return np.column_stack((np.ones(local.shape[0]), s, t, r,
                            s * s, s * t, s * r,
                            t * t, t * r, r * r))


def dirichlet_cell_average_basis(
    local: np.ndarray, covariance: np.ndarray
) -> np.ndarray:
    if local.shape[1] == 2:
        s, t = local.T
        return np.column_stack((s, t,
                                s * s + covariance[0, 0],
                                s * t + covariance[0, 1],
                                t * t + covariance[1, 1]))
    s, t, r = local.T
    return np.column_stack((s, t, r,
                            s * s + covariance[0, 0],
                            s * t + covariance[0, 1],
                            s * r + covariance[0, 2],
                            t * t + covariance[1, 1],
                            t * r + covariance[1, 2],
                            r * r + covariance[2, 2]))


def target_weights(support_basis: np.ndarray,
                   target_basis: np.ndarray) -> np.ndarray:
    normal = support_basis.T @ support_basis
    if np.linalg.matrix_rank(normal, tol=1.0e-11) != normal.shape[0]:
        raise AssertionError("test support unexpectedly rank deficient")
    return support_basis @ np.linalg.solve(normal, target_basis)


def run_campaign(samples: int, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    max_fv_error = 0.0
    max_jet_error = 0.0
    max_linear_jet_error = 0.0
    max_point_semantic_defect = 0.0

    for dim in (2, 3):
        offsets = np.asarray(
            list(itertools.product((0.5, 1.5, 2.5), repeat=dim)),
            dtype=float,
        )
        for _ in range(samples):
            spacing = rng.uniform(0.35, 1.4, size=dim)
            centres = offsets * spacing
            frame = orthonormal_frame(dim, rng)
            scale = float(np.linalg.norm(spacing))

            # Put the BI immediately behind the entire Cartesian support block,
            # so every selected cell belongs to the same one-sided fluid jet.
            projection = centres @ frame[0]
            wall = rng.uniform(-0.2, 0.2, size=dim)
            wall += (float(np.min(projection - wall @ frame[0])) -
                     rng.uniform(0.15, 0.55) * scale) * frame[0]
            ghost = wall - rng.uniform(0.2, 1.1) * scale * frame[0]

            support_local = coordinates(centres, wall, frame, scale)
            ghost_local = coordinates(ghost[None, :], wall, frame, scale)
            covariance = np.zeros((dim, dim))
            for a in range(dim):
                for b in range(dim):
                    covariance[a, b] = np.sum(
                        spacing**2 * frame[a] * frame[b]
                    ) / (12.0 * scale**2)

            support_fv = dirichlet_cell_average_basis(
                support_local, covariance
            )
            target_fv = dirichlet_cell_average_basis(
                ghost_local, covariance
            )[0]
            fv_weights = target_weights(support_fv, target_fv)
            boundary_weight = 1.0 - float(np.sum(fv_weights))
            coefficients = rng.normal(size=support_fv.shape[1])
            boundary_value = float(rng.normal())
            support_values = boundary_value + support_fv @ coefficients
            exact_target = boundary_value + float(target_fv @ coefficients)
            reconstructed = (
                float(fv_weights @ support_values)
                + boundary_weight * boundary_value
            )
            fv_error = abs(reconstructed - exact_target)
            max_fv_error = max(max_fv_error, fv_error)

            support_jet = full_basis(support_local)
            target_jet = full_basis(ghost_local)[0]
            jet_weights = target_weights(support_jet, target_jet)
            jet_coefficients = rng.normal(size=support_jet.shape[1])
            jet_error = abs(
                float(jet_weights @ (support_jet @ jet_coefficients))
                - float(target_jet @ jet_coefficients)
            )
            max_jet_error = max(max_jet_error, jet_error)

            # The embedded P1 jet uses the leading {1,s,t[,r]} columns of
            # exactly the same support matrix.  It must reproduce every
            # linear field independently of the P2 solve.
            linear_columns = dim + 1
            support_linear = support_jet[:, :linear_columns]
            target_linear = target_jet[:linear_columns]
            linear_weights = target_weights(
                support_linear, target_linear
            )
            linear_coefficients = rng.normal(size=linear_columns)
            linear_jet_error = abs(
                float(
                    linear_weights
                    @ (support_linear @ linear_coefficients)
                )
                - float(target_linear @ linear_coefficients)
            )
            max_linear_jet_error = max(
                max_linear_jet_error, linear_jet_error
            )

            # Quantify why a centre-value Dirichlet functional is not a valid
            # replacement for the finite-volume functional.
            support_point = dirichlet_cell_average_basis(
                support_local, np.zeros_like(covariance)
            )
            target_point = dirichlet_cell_average_basis(
                ghost_local, np.zeros_like(covariance)
            )[0]
            point_weights = target_weights(support_point, target_point)
            point_boundary_weight = 1.0 - float(np.sum(point_weights))
            point_reconstruction = (
                float(point_weights @ support_values)
                + point_boundary_weight * boundary_value
            )
            max_point_semantic_defect = max(
                max_point_semantic_defect,
                abs(point_reconstruction - exact_target),
            )

    return {
        "cell_average_dirichlet_max_error": max_fv_error,
        "one_sided_jet_max_error": max_jet_error,
        "embedded_linear_jet_max_error": max_linear_jet_error,
        "point_functional_on_cell_averages_max_defect":
            max_point_semantic_defect,
    }


def limiter_campaign() -> dict[str, float]:
    """Verify smooth-extremum asymptotics and discontinuity rejection.

    The production hierarchy sensor regularises the local support range with
    sqrt(|P2-P1| q_scale).  It is therefore O(h) both at ordinary smooth
    points, where Delta q=O(h), and at smooth extrema, where Delta q=O(h^2).
    The independent non-smoothness sensor still uses the support range.  At
    an unresolved jump both ratios remain O(1).
    """

    tolerance = 0.10
    # Use h << tolerance so this is an asymptotic-order test rather than a
    # transition-region calibration of the limiter.
    spacing = np.asarray([1.0 / (2**level) for level in range(5, 11)])

    hierarchy_difference = spacing**2
    smooth_extremum_range = spacing**2
    hierarchy_ratio = hierarchy_difference / (
        smooth_extremum_range + np.sqrt(hierarchy_difference)
    )
    hierarchy_power = (hierarchy_ratio / tolerance) ** 4
    hierarchy_fraction = 1.0 / (1.0 + hierarchy_power)
    hierarchy_perturbation = spacing**2 * (
        hierarchy_power / (1.0 + hierarchy_power)
    )

    variation_ratio = spacing
    variation_power = (variation_ratio / tolerance) ** 4
    variation_fraction = 1.0 / (1.0 + variation_power)
    variation_perturbation = spacing * (
        variation_power / (1.0 + variation_power)
    )

    def fitted_order(error: np.ndarray) -> float:
        return float(np.polyfit(np.log(spacing), np.log(error), 1)[0])

    shock_hierarchy_ratio = 0.5 / (1.0 + np.sqrt(0.5))
    shock_hierarchy_fraction = 1.0 / (
        1.0 + (shock_hierarchy_ratio / tolerance) ** 4
    )
    shock_variation_fraction = 1.0 / (1.0 + (1.0 / tolerance) ** 4)
    return {
        "smooth_hierarchy_limiter_perturbation_order": fitted_order(
            hierarchy_perturbation
        ),
        "smooth_variation_limiter_perturbation_order": fitted_order(
            variation_perturbation
        ),
        "shock_hierarchy_high_order_fraction": shock_hierarchy_fraction,
        "shock_variation_high_order_fraction": shock_variation_fraction,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260719)
    args = parser.parse_args()
    result = run_campaign(args.samples, args.seed)
    result.update(limiter_campaign())
    for key, value in result.items():
        print(f"{key}={value:.16e}")
    if result["cell_average_dirichlet_max_error"] > 2.0e-9:
        raise SystemExit("finite-volume Dirichlet reproduction failed")
    if result["one_sided_jet_max_error"] > 2.0e-9:
        raise SystemExit("one-sided jet reproduction failed")
    if result["embedded_linear_jet_max_error"] > 2.0e-9:
        raise SystemExit("embedded linear jet reproduction failed")
    if result["point_functional_on_cell_averages_max_defect"] < 1.0e-7:
        raise SystemExit("test did not expose the point/cell-average mismatch")
    if result["smooth_hierarchy_limiter_perturbation_order"] < 5.0:
        raise SystemExit("hierarchy limiter is not high-order inactive")
    if result["smooth_variation_limiter_perturbation_order"] < 4.9:
        raise SystemExit("variation limiter is not high-order inactive")
    if result["shock_hierarchy_high_order_fraction"] > 2.0e-2:
        raise SystemExit("hierarchy limiter did not reject a jump")
    if result["shock_variation_high_order_fraction"] > 2.0e-4:
        raise SystemExit("variation limiter did not reject a jump")


if __name__ == "__main__":
    main()
