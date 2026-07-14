#!/usr/bin/env python3
"""Regression tests for IBM support visibility and filtered WLS diagnostics.

The geometry routines mirror the open-segment predicates in ``ibm_bvh.h``.
The WLS checks reuse the independent reference implementation in
``test_ibm_wls_reproduction.py``.  This is a fast mathematical gate; full
Cerisse integration tests still exercise the compiled BVH and IBM containers.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from test_ibm_wls_reproduction import (
    basis,
    grid_basis,
    interpolation_weights,
    verify_mask,
)


MACHINE_EPSILON = np.finfo(float).eps
ENDPOINT_TOLERANCE = 1.0e-10
RANK_REL_TOLERANCE = 1.0e-11


def cross_2d(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def segment_edge_first_fraction(
    start: np.ndarray,
    end: np.ndarray,
    edge_start: np.ndarray,
    edge_end: np.ndarray,
    endpoint_tolerance: float = ENDPOINT_TOLERANCE,
) -> float | None:
    """Mirror ``segment_edge_first_fraction`` from the C++ BVH."""
    direction = end - start
    edge = edge_end - edge_start
    offset = edge_start - start
    direction_norm2 = float(direction @ direction)
    edge_norm2 = float(edge @ edge)
    if direction_norm2 <= 1.0e-28 or edge_norm2 <= 1.0e-28:
        return None

    denominator = cross_2d(direction, edge)
    parallel_tolerance = (
        64.0 * MACHINE_EPSILON * math.sqrt(direction_norm2 * edge_norm2)
    )
    if abs(denominator) > parallel_tolerance:
        fraction = cross_2d(offset, edge) / denominator
        edge_fraction = cross_2d(offset, direction) / denominator
        side_tolerance = 32.0 * MACHINE_EPSILON
        if (
            endpoint_tolerance < fraction < 1.0 - endpoint_tolerance
            and -side_tolerance <= edge_fraction <= 1.0 + side_tolerance
        ):
            return fraction
        return None

    if abs(cross_2d(offset, direction)) > parallel_tolerance:
        return None
    axis = 0 if abs(direction[0]) >= abs(direction[1]) else 1
    first = (edge_start[axis] - start[axis]) / direction[axis]
    second = (edge_end[axis] - start[axis]) / direction[axis]
    overlap_begin = max(min(first, second), endpoint_tolerance)
    overlap_end = min(max(first, second), 1.0 - endpoint_tolerance)
    return overlap_begin if overlap_begin < overlap_end else None


def first_segment_hit(
    start: tuple[float, float],
    end: tuple[float, float],
    indexed_edges: list[tuple[int, tuple[float, float], tuple[float, float]]],
) -> tuple[int, float] | None:
    """Return the first hit, breaking equal-fraction ties by primitive id."""
    start_array = np.asarray(start, dtype=float)
    end_array = np.asarray(end, dtype=float)
    best_primitive = -1
    best_fraction = math.inf
    tie_tolerance = 64.0 * MACHINE_EPSILON
    for primitive, edge_start, edge_end in indexed_edges:
        fraction = segment_edge_first_fraction(
            start_array,
            end_array,
            np.asarray(edge_start, dtype=float),
            np.asarray(edge_end, dtype=float),
        )
        if fraction is None:
            continue
        if (
            fraction < best_fraction - tie_tolerance
            or (
                abs(fraction - best_fraction) <= tie_tolerance
                and primitive < best_primitive
            )
        ):
            best_primitive = primitive
            best_fraction = fraction
    return None if best_primitive < 0 else (best_primitive, best_fraction)


def polygon_edges(
    vertices: list[tuple[float, float]],
) -> list[tuple[int, tuple[float, float], tuple[float, float]]]:
    return [
        (index, vertex, vertices[(index + 1) % len(vertices)])
        for index, vertex in enumerate(vertices)
    ]


def rigid_transform(
    point: tuple[float, float], angle: float, translation: np.ndarray
) -> tuple[float, float]:
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle)],
         [math.sin(angle), math.cos(angle)]]
    )
    transformed = rotation @ np.asarray(point) + translation
    return float(transformed[0]), float(transformed[1])


def test_geometry_predicates() -> None:
    square = polygon_edges([(0.0, 0.0), (1.0, 0.0),
                            (1.0, 1.0), (0.0, 1.0)])
    hit = first_segment_hit((-1.0, 0.5), (2.0, 0.5), square)
    assert hit is not None and hit[0] == 3
    assert abs(hit[1] - 1.0 / 3.0) < 1.0e-14

    # A support exactly on the surface endpoint is not self-occluded.
    assert first_segment_hit((-1.0, 0.5), (0.0, 0.5), square) is None

    # Two edges meet at this vertex; traversal order must not change ownership.
    vertex_hit = first_segment_hit((-1.0, -1.0), (0.5, 0.5), square)
    reverse_hit = first_segment_hit((-1.0, -1.0), (0.5, 0.5), square[::-1])
    assert vertex_hit is not None and reverse_hit is not None
    assert vertex_hit[0] == reverse_hit[0] == 0
    assert abs(vertex_hit[1] - reverse_hit[1]) < 1.0e-15

    # Collinear overlap is an obstruction, not a parallel-line miss.
    overlap = first_segment_hit((-0.5, 0.0), (0.5, 0.0), square)
    assert overlap is not None and overlap[0] == 0

    # Concave L: both endpoints lie outside the solid polygon, but their direct
    # support segment crosses the horizontal arm.
    concave = polygon_edges(
        [(0.0, 0.0), (3.0, 0.0), (3.0, 1.0),
         (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)]
    )
    concave_hit = first_segment_hit((2.0, 2.0), (2.0, -1.0), concave)
    assert concave_hit is not None and concave_hit[0] == 2

    # Translation/rotation of all geometry and query points preserves hit id
    # and segment fraction.
    angle = 0.731
    translation = np.array([2.3, -4.7])
    transformed_edges = [
        (primitive,
         rigid_transform(edge_start, angle, translation),
         rigid_transform(edge_end, angle, translation))
        for primitive, edge_start, edge_end in concave
    ]
    transformed_hit = first_segment_hit(
        rigid_transform((2.0, 2.0), angle, translation),
        rigid_transform((2.0, -1.0), angle, translation),
        transformed_edges,
    )
    assert transformed_hit is not None
    assert transformed_hit[0] == concave_hit[0]
    assert abs(transformed_hit[1] - concave_hit[1]) < 2.0e-14


def numerical_rank(normal_matrix: np.ndarray) -> int:
    tolerance = RANK_REL_TOLERANCE * float(np.max(np.abs(normal_matrix)))
    singular_values = np.linalg.svd(normal_matrix, compute_uv=False)
    return int(np.count_nonzero(singular_values > tolerance))


def test_filtered_wls() -> None:
    for dim, minimum_points in ((2, 2), (3, 3)):
        all_basis = grid_basis(dim, tuple([1.37] * dim))
        full_mask = np.ones(all_basis.shape[0], dtype=bool)
        order, error = verify_mask(
            all_basis, full_mask, dim, minimum_points=minimum_points
        )
        assert order == 2 and error < 2.0e-10

        normal_matrix = all_basis.T @ all_basis
        expected_size = 6 if dim == 2 else 10
        assert numerical_rank(normal_matrix) == expected_size
        condition = np.linalg.cond(normal_matrix, p=1)
        assert np.isfinite(condition) and condition >= 1.0

    # Six 2-D points are a sufficient count for a quadratic, but two y-levels
    # make y^2 dependent on {1,y}; rank/conditioning must force linear demotion.
    offsets_2d = np.array(list(itertools.product(range(3), repeat=2)))
    two_rows = offsets_2d[:, 1] <= 1
    basis_2d = grid_basis(2, (1.27, 1.63))
    matrix_2d = basis_2d[two_rows].T @ basis_2d[two_rows]
    order_2d, _ = interpolation_weights(basis_2d, two_rows, 2, 2)
    assert int(two_rows.sum()) == 6
    assert numerical_rank(matrix_2d) < 6
    assert order_2d == 1
    verify_mask(basis_2d, two_rows, 2, minimum_points=2)

    # Eighteen 3-D points on two z-levels similarly make z^2 dependent while
    # retaining a complete affine basis.
    offsets_3d = np.array(list(itertools.product(range(3), repeat=3)))
    two_planes = offsets_3d[:, 2] <= 1
    basis_3d = grid_basis(3, (1.21, 1.47, 1.72))
    matrix_3d = basis_3d[two_planes].T @ basis_3d[two_planes]
    order_3d, _ = interpolation_weights(basis_3d, two_planes, 3, 3)
    assert int(two_planes.sum()) == 18
    assert numerical_rank(matrix_3d) < 10
    assert order_3d == 1
    verify_mask(basis_3d, two_planes, 3, minimum_points=3)

    # A sparse but geometrically complete six-point 2-D stencil must retain
    # exact quadratic reproduction.
    sparse_basis = grid_basis(2, (1.19, 1.81))
    complete_mask = None
    for selected in itertools.combinations(range(9), 6):
        candidate = np.zeros(9, dtype=bool)
        candidate[list(selected)] = True
        if np.linalg.matrix_rank(sparse_basis[candidate]) == 6:
            complete_mask = candidate
            break
    assert complete_mask is not None
    sparse_order, sparse_error = verify_mask(
        sparse_basis, complete_mask, 2, minimum_points=2
    )
    assert sparse_order == 2 and sparse_error < 2.0e-10


def main() -> int:
    test_geometry_predicates()
    test_filtered_wls()
    print("IBM visibility geometry and filtered WLS diagnostics: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
