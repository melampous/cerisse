#!/usr/bin/env python3
"""Validate Cerisse Cartesian IBM boundary-face audit records independently."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from test_ibm_visibility_audit import first_segment_hit, polygon_edges


def read_polygon(path: Path) -> list[tuple[float, float]]:
    vertices: list[tuple[float, float]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        vertices.append((float(fields[0]), float(fields[1])))
    if len(vertices) < 3:
        raise ValueError(f"{path} does not contain a 2-D closed polygon")
    return vertices


def row_vector(row: dict[str, str], prefix: str, dim: int) -> np.ndarray:
    axes = ("x", "y", "z")[:dim]
    return np.asarray([float(row[f"{prefix}_{axis}"]) for axis in axes])


def records_agree(
    left: dict[str, str], right: dict[str, str], dim: int
) -> bool:
    integer_fields = ("first_hit_geometry", "first_hit_element")
    if any(left[field] != right[field] for field in integer_fields):
        return False
    if int(left["first_hit_element"]) < 0:
        return True
    scalar_tolerance = 2.0e-12
    if not math.isclose(
        float(left["first_hit_fraction"]),
        float(right["first_hit_fraction"]),
        rel_tol=scalar_tolerance,
        abs_tol=scalar_tolerance,
    ):
        return False
    return bool(
        np.allclose(
            row_vector(left, "hit", dim),
            row_vector(right, "hit", dim),
            rtol=scalar_tolerance,
            atol=scalar_tolerance,
        )
        and np.allclose(
            row_vector(left, "normal", dim),
            row_vector(right, "normal", dim),
            rtol=scalar_tolerance,
            atol=scalar_tolerance,
        )
    )


def cpp_numerical_rank(matrix: np.ndarray) -> int:
    """Reproduce ``ibm_matrix_numerical_rank`` and its relative tolerance."""
    work = matrix.copy()
    scale = float(np.max(np.abs(work)))
    if not scale > 0.0 or not math.isfinite(scale):
        return 0
    tolerance = 1.0e-11 * scale
    rank = 0
    size = work.shape[0]
    for column in range(size):
        if rank >= size:
            break
        pivot = rank + int(np.argmax(np.abs(work[rank:, column])))
        pivot_abs = abs(float(work[pivot, column]))
        if not pivot_abs > tolerance:
            continue
        if pivot != rank:
            work[[rank, pivot], column:] = work[[pivot, rank], column:]
        inverse = 1.0 / work[rank, column]
        for row in range(rank + 1, size):
            factor = work[row, column] * inverse
            work[row, column:] -= factor * work[rank, column:]
        rank += 1
    return rank


def guarded_cholesky(matrix: np.ndarray) -> np.ndarray | None:
    """Reproduce the guarded small-SPD factorisation used by IBM WLS."""
    size = matrix.shape[0]
    tolerance = 1.0e-11 * max(float(np.max(np.diag(matrix))), 1.0e-300)
    factor = np.zeros_like(matrix)
    for row in range(size):
        diagonal = matrix[row, row] - float(
            factor[row, :row] @ factor[row, :row]
        )
        if not diagonal > tolerance or not math.isfinite(diagonal):
            return None
        factor[row, row] = math.sqrt(diagonal)
        for lower_row in range(row + 1, size):
            value = matrix[lower_row, row] - float(
                factor[lower_row, :row] @ factor[row, :row]
            )
            factor[lower_row, row] = value / factor[row, row]
    return factor


def cpp_condition_1(matrix: np.ndarray) -> float:
    """Reproduce ``ibm_spd_condition_1`` for a small normal matrix."""
    factor = guarded_cholesky(matrix)
    if factor is None:
        return math.inf
    size = matrix.shape[0]
    inverse = np.zeros_like(matrix)
    for column in range(size):
        rhs = np.zeros(size)
        rhs[column] = 1.0
        work = np.zeros(size)
        for row in range(size):
            work[row] = (
                rhs[row] - factor[row, :row] @ work[:row]
            ) / factor[row, row]
        solution = np.zeros(size)
        for row in range(size - 1, -1, -1):
            solution[row] = (
                work[row]
                - factor[row + 1 :, row] @ solution[row + 1 :]
            ) / factor[row, row]
        inverse[:, column] = solution
    condition = float(np.linalg.norm(matrix, 1) * np.linalg.norm(inverse, 1))
    return condition if math.isfinite(condition) else math.inf


def polynomial_basis(offsets: np.ndarray) -> np.ndarray:
    one = np.ones(offsets.shape[0])
    if offsets.shape[1] == 2:
        x, y = offsets.T
        return np.column_stack((one, x, y, x * x, x * y, y * y))
    x, y, z = offsets.T
    return np.column_stack(
        (one, x, y, z, x * x, y * y, z * z, x * y, x * z, y * z)
    )


def floating_values_agree(left: float, right: float, tolerance: float) -> bool:
    if math.isinf(left) or math.isinf(right):
        return left == right
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def reconstruction_records_agree(
    left: dict[str, str],
    right: dict[str, str],
    left_supports: list[dict[str, str]],
    right_supports: list[dict[str, str]],
    dim: int,
) -> bool:
    integer_fields = (
        "n_candidates",
        "n_visible",
        "requested_order",
        "final_order",
        "normal_rank",
        "accepted_rank",
        "reconstruction_valid",
    )
    for field in integer_fields:
        if field in left and field in right and left[field] != right[field]:
            return False
    tolerance = 8.0e-12
    for field in ("solid_normal_distance", "normal_cond1", "accepted_cond1"):
        if field in left and field in right and not floating_values_agree(
            float(left[field]), float(right[field]), tolerance
        ):
            return False
    if not np.allclose(
        row_vector(left, "mirror", dim),
        row_vector(right, "mirror", dim),
        rtol=tolerance,
        atol=tolerance,
    ):
        return False
    if len(left_supports) != len(right_supports):
        return False
    support_integer_fields = (
        "support_slot",
        "support_i",
        "support_j",
        "support_k",
        "candidate",
        "visible",
        "first_intersection_element",
    )
    for left_support, right_support in zip(left_supports, right_supports):
        if any(
            left_support[field] != right_support[field]
            for field in support_integer_fields
        ):
            return False
        if not floating_values_agree(
            float(left_support["weight"]),
            float(right_support["weight"]),
            tolerance,
        ):
            return False
    return True


def validate_reconstruction(
    row: dict[str, str],
    supports: list[dict[str, str]],
    dim: int,
    prob_lo: np.ndarray,
    spacing: np.ndarray,
    solid_point: np.ndarray,
    hit_point: np.ndarray,
    normal: np.ndarray,
    condition_max: float,
    moment_tolerance: float,
    label: str,
) -> list[str]:
    failures: list[str] = []
    signed_distance = float((solid_point - hit_point) @ normal)
    if not signed_distance < 0.0:
        failures.append(f"{label}: solid centre is not behind first-hit plane")
    expected_mirror = solid_point - 2.0 * signed_distance * normal
    mirror = row_vector(row, "mirror", dim)
    if not np.allclose(
        mirror, expected_mirror, rtol=3.0e-12, atol=3.0e-12
    ):
        failures.append(f"{label}: reflected solid centre is inconsistent")
    reported_distance = float(row["solid_normal_distance"])
    if not math.isclose(
        reported_distance,
        -signed_distance,
        rel_tol=3.0e-12,
        abs_tol=3.0e-12,
    ):
        failures.append(f"{label}: solid normal distance is inconsistent")

    supports = sorted(supports, key=lambda item: int(item["support_slot"]))
    expected_slots = list(range(3**dim))
    slots = [int(item["support_slot"]) for item in supports]
    if slots != expected_slots:
        failures.append(f"{label}: support slots are not exactly {expected_slots}")
        return failures

    face_fields = ("face_i", "face_j", "face_k", "direction", "fluid_side")
    support_indices: list[list[int]] = []
    candidates: list[bool] = []
    visible: list[bool] = []
    weights: list[float] = []
    for support in supports:
        if support["kind"] != "cartesian_boundary_support":
            failures.append(
                f"{label}: unexpected support kind={support['kind']}"
            )
        if any(support[field] != row[field] for field in face_fields):
            failures.append(f"{label}: support row belongs to another face")
        if "local_fab" in row and "local_fab" in support:
            if support["local_fab"] != row["local_fab"]:
                failures.append(f"{label}: support local FAB is inconsistent")
        candidate = int(support["candidate"])
        is_visible = int(support["visible"])
        if candidate not in (0, 1) or is_visible not in (0, 1):
            failures.append(f"{label}: support mask is not binary")
        if is_visible and not candidate:
            failures.append(f"{label}: visible support is not a fluid candidate")
        first_obstruction = int(support["first_intersection_element"])
        if is_visible and first_obstruction >= 0:
            failures.append(f"{label}: visible support reports an obstruction")
        if candidate and not is_visible and first_obstruction < 0:
            failures.append(f"{label}: rejected fluid support lacks obstruction")
        support_indices.append(
            [int(support[f"support_{axis}"]) for axis in ("i", "j", "k")[:dim]]
        )
        candidates.append(bool(candidate))
        visible.append(bool(is_visible))
        weights.append(float(support["weight"]))

    candidate_mask = np.asarray(candidates, dtype=bool)
    visible_mask = np.asarray(visible, dtype=bool)
    weight_array = np.asarray(weights)
    if int(row["n_candidates"]) != int(np.count_nonzero(candidate_mask)):
        failures.append(f"{label}: candidate count disagrees with support mask")
    if int(row["n_visible"]) != int(np.count_nonzero(visible_mask)):
        failures.append(f"{label}: visible count disagrees with support mask")
    if np.any(np.abs(weight_array[~visible_mask]) > moment_tolerance):
        failures.append(f"{label}: rejected support has nonzero WLS weight")

    final_order = int(row["final_order"])
    if int(row["requested_order"]) != 2:
        failures.append(f"{label}: face-local requested order is not quadratic")
    if int(row["reconstruction_valid"]) != int(final_order >= 1):
        failures.append(f"{label}: reconstruction-valid flag is inconsistent")

    if final_order < 0:
        if np.any(np.abs(weight_array) > moment_tolerance):
            failures.append(f"{label}: invalid reconstruction has nonzero weights")
        return failures

    centres = prob_lo + (np.asarray(support_indices) + 0.5) * spacing
    basis = polynomial_basis((centres - mirror) / spacing)
    active_basis = basis[visible_mask]
    requested_matrix = active_basis.T @ active_basis
    requested_rank = cpp_numerical_rank(requested_matrix)
    requested_condition = (
        cpp_condition_1(requested_matrix)
        if requested_rank == requested_matrix.shape[0]
        else math.inf
    )
    if int(row["normal_rank"]) != requested_rank:
        failures.append(
            f"{label}: requested rank {row['normal_rank']} != {requested_rank}"
        )
    if not floating_values_agree(
        float(row["normal_cond1"]), requested_condition, 2.0e-10
    ):
        failures.append(
            f"{label}: requested condition {row['normal_cond1']} "
            f"!= {requested_condition:.17g}"
        )

    linear_size = dim + 1
    quadratic_size = 6 if dim == 2 else 10
    visible_count = int(np.count_nonzero(visible_mask))
    quadratic_ok = (
        visible_count >= quadratic_size
        and math.isfinite(requested_condition)
        and requested_condition <= condition_max
    )
    linear_matrix = requested_matrix[:linear_size, :linear_size]
    linear_rank = cpp_numerical_rank(linear_matrix)
    linear_condition = (
        cpp_condition_1(linear_matrix)
        if linear_rank == linear_size
        else math.inf
    )
    linear_ok = (
        visible_count >= linear_size
        and math.isfinite(linear_condition)
        and linear_condition <= condition_max
    )
    minimum_points = 2 if dim == 2 else 3
    expected_order = (
        2
        if quadratic_ok
        else 1
        if linear_ok
        else 0
        if visible_count >= minimum_points
        else -1
    )
    if final_order != expected_order:
        failures.append(
            f"{label}: final order {final_order} != fit-ladder order "
            f"{expected_order}"
        )

    accepted_size = {0: 1, 1: linear_size, 2: quadratic_size}[final_order]
    accepted_matrix = requested_matrix[:accepted_size, :accepted_size]
    accepted_rank = cpp_numerical_rank(accepted_matrix)
    accepted_condition = cpp_condition_1(accepted_matrix)
    if final_order == 0:
        accepted_rank = 1
        accepted_condition = 1.0
    if "accepted_rank" in row and int(row["accepted_rank"]) != accepted_rank:
        failures.append(
            f"{label}: accepted rank {row['accepted_rank']} != {accepted_rank}"
        )
    if "accepted_cond1" in row and not floating_values_agree(
        float(row["accepted_cond1"]), accepted_condition, 2.0e-10
    ):
        failures.append(
            f"{label}: accepted condition {row['accepted_cond1']} "
            f"!= {accepted_condition:.17g}"
        )
    if accepted_rank != accepted_size or accepted_condition > condition_max:
        failures.append(f"{label}: accepted WLS system is rank deficient/ill-conditioned")

    moments = basis[:, :accepted_size].T @ weight_array
    target = np.zeros(accepted_size)
    target[0] = 1.0
    moment_error = float(np.max(np.abs(moments - target)))
    if not math.isfinite(moment_error) or moment_error > moment_tolerance:
        failures.append(
            f"{label}: polynomial reproduction error={moment_error:.6e}"
        )
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv_files",
        type=Path,
        nargs="+",
        help="One or more per-rank boundary-face CSV files.",
    )
    parser.add_argument("--dim", type=int, choices=(2, 3), required=True)
    parser.add_argument("--prob-lo", type=float, nargs="+", required=True)
    parser.add_argument("--prob-hi", type=float, nargs="+", required=True)
    parser.add_argument("--n-cell", type=int, nargs="+", required=True)
    parser.add_argument(
        "--geometry",
        type=Path,
        help="Optional 2-D polygon for an independent first-hit check.",
    )
    parser.add_argument(
        "--support-files",
        type=Path,
        nargs="+",
        help=(
            "Per-rank support CSV files in the same order as csv_files. "
            "By default, derive each name by replacing boundary_faces with "
            "boundary_supports."
        ),
    )
    parser.add_argument(
        "--skip-supports",
        action="store_true",
        help="Skip reconstruction/support checks and validate first-hit only.",
    )
    parser.add_argument(
        "--condition-max",
        type=float,
        default=1.0e10,
        help="Condition-number ceiling compiled into the audited case.",
    )
    parser.add_argument(
        "--moment-tol",
        type=float,
        default=2.0e-10,
        help="Maximum accepted polynomial reproduction residual.",
    )
    parser.add_argument("--require-owner-difference", action="store_true")
    args = parser.parse_args()

    dim = args.dim
    if not (len(args.prob_lo) == len(args.prob_hi) == len(args.n_cell) == dim):
        parser.error("domain bounds and n-cell must each contain --dim values")
    if args.geometry is not None and dim != 2:
        parser.error("--geometry independent first-hit validation is 2-D only")
    if args.skip_supports and args.support_files is not None:
        parser.error("--skip-supports and --support-files are mutually exclusive")
    if args.support_files is not None and len(args.support_files) != len(
        args.csv_files
    ):
        parser.error("--support-files must match the number of face CSV files")
    if not args.condition_max > 1.0 or not math.isfinite(args.condition_max):
        parser.error("--condition-max must be finite and greater than one")
    if not args.moment_tol > 0.0 or not math.isfinite(args.moment_tol):
        parser.error("--moment-tol must be finite and positive")

    prob_lo = np.asarray(args.prob_lo, dtype=float)
    prob_hi = np.asarray(args.prob_hi, dtype=float)
    n_cell = np.asarray(args.n_cell, dtype=int)
    if np.any(n_cell <= 0) or np.any(prob_hi <= prob_lo):
        parser.error("invalid domain or cell counts")
    spacing = (prob_hi - prob_lo) / n_cell
    edges = (
        polygon_edges(read_polygon(args.geometry))
        if args.geometry is not None
        else None
    )

    rows: list[dict[str, str]] = []
    supports_by_copy: dict[str, list[dict[str, str]]] = {}
    support_file_count = 0
    for file_index, csv_file in enumerate(args.csv_files):
        file_rows: list[dict[str, str]] = []
        with csv_file.open(newline="", encoding="ascii") as stream:
            for copy_index, row in enumerate(csv.DictReader(stream)):
                row["_source"] = str(csv_file)
                copy_slot = row.get("face_copy_slot", str(copy_index))
                token = f"{csv_file.resolve()}::{row['mpi_rank']}::{copy_slot}"
                row["_copy_token"] = token
                file_rows.append(row)
                rows.append(row)
        if args.skip_supports:
            continue
        if args.support_files is not None:
            support_file = args.support_files[file_index]
        else:
            support_name = csv_file.name.replace(
                "_boundary_faces_", "_boundary_supports_"
            )
            if support_name == csv_file.name:
                raise RuntimeError(
                    f"cannot derive support CSV name from {csv_file}"
                )
            support_file = csv_file.with_name(support_name)
        if not support_file.is_file():
            raise RuntimeError(f"missing boundary support CSV: {support_file}")
        support_file_count += 1
        with support_file.open(newline="", encoding="ascii") as stream:
            support_rows = list(csv.DictReader(stream))

        number_of_supports = 3**dim
        if support_rows and "face_copy_slot" in support_rows[0]:
            face_by_slot = {
                (row["mpi_rank"], row["face_copy_slot"]): row
                for row in file_rows
            }
            grouped_supports: dict[
                tuple[str, str], list[dict[str, str]]
            ] = defaultdict(list)
            for support in support_rows:
                grouped_supports[
                    (support["mpi_rank"], support["face_copy_slot"])
                ].append(support)
            unknown = set(grouped_supports) - set(face_by_slot)
            if unknown:
                raise RuntimeError(
                    f"{support_file}: support rows reference unknown face copies "
                    f"{sorted(unknown)[:4]}"
                )
            for slot_key, face_row in face_by_slot.items():
                supports = grouped_supports.get(slot_key, [])
                if len(supports) != number_of_supports:
                    raise RuntimeError(
                        f"{support_file}: face copy {slot_key} has "
                        f"{len(supports)} supports, expected {number_of_supports}"
                    )
                supports_by_copy[face_row["_copy_token"]] = supports
        else:
            expected = len(file_rows) * number_of_supports
            if len(support_rows) != expected:
                raise RuntimeError(
                    f"{support_file}: found {len(support_rows)} support rows, "
                    f"expected {expected}"
                )
            for copy_index, face_row in enumerate(file_rows):
                begin = copy_index * number_of_supports
                supports_by_copy[face_row["_copy_token"]] = support_rows[
                    begin : begin + number_of_supports
                ]
    if not rows:
        raise RuntimeError("input CSV files contain no boundary faces")

    groups: dict[tuple[int, ...], list[dict[str, str]]] = defaultdict(list)
    unsided: dict[tuple[int, ...], int] = {}
    failures: list[str] = []
    owner_differences = 0
    unavailable_legacy_owners = 0
    axes = ("i", "j", "k")[:dim]
    for row_number, row in enumerate(rows, start=2):
        if row["kind"] != "cartesian_boundary_face":
            failures.append(f"row {row_number}: unexpected kind={row['kind']}")
        face = np.asarray([int(row[f"face_{axis}"]) for axis in axes])
        direction = int(row["direction"])
        fluid_side = int(row["fluid_side"])
        if not 0 <= direction < dim:
            failures.append(f"row {row_number}: invalid direction={direction}")
            continue
        if fluid_side not in (-1, 1):
            failures.append(f"row {row_number}: invalid fluid_side={fluid_side}")
            continue

        key = (*face.tolist(), direction, fluid_side)
        groups[key].append(row)
        unsided_key = (*face.tolist(), direction)
        previous_side = unsided.setdefault(unsided_key, fluid_side)
        if previous_side != fluid_side:
            failures.append(
                f"row {row_number}: one numerical face has opposing fluid sides"
            )

        hit_element = int(row["first_hit_element"])
        hit_geometry = int(row["first_hit_geometry"])
        fraction = float(row["first_hit_fraction"])
        if hit_element < 0 or hit_geometry < 0:
            failures.append(f"row {row_number}: missing first hit")
            continue
        if not (1.0e-10 < fraction < 1.0 - 1.0e-10):
            failures.append(
                f"row {row_number}: hit fraction is outside open segment: {fraction}"
            )

        fluid_cell = face.copy()
        if fluid_side < 0:
            fluid_cell[direction] -= 1
        ghost_cell = fluid_cell.copy()
        ghost_cell[direction] -= fluid_side
        fluid_point = prob_lo + (fluid_cell + 0.5) * spacing
        ghost_point = prob_lo + (ghost_cell + 0.5) * spacing
        expected_hit = fluid_point + fraction * (ghost_point - fluid_point)
        hit_point = row_vector(row, "hit", dim)
        if not np.allclose(hit_point, expected_hit, rtol=2.0e-12, atol=2.0e-12):
            failures.append(f"row {row_number}: hit point disagrees with fraction")

        normal = row_vector(row, "normal", dim)
        normal_norm = float(np.linalg.norm(normal))
        if not math.isfinite(normal_norm) or abs(normal_norm - 1.0) > 2.0e-12:
            failures.append(
                f"row {row_number}: non-unit/non-finite normal norm={normal_norm}"
            )
        alignment = float(normal @ (fluid_point - hit_point))
        if not alignment > 0.0:
            failures.append(
                f"row {row_number}: normal does not point to fluid side"
            )
        if "local_fab" in row and row["local_fab"] != row["gp_fab"]:
            failures.append(
                f"row {row_number}: local_fab and reconstruction FAB disagree"
            )

        if not args.skip_supports:
            supports = supports_by_copy.get(row["_copy_token"])
            if supports is None:
                failures.append(f"row {row_number}: no support rows for face copy")
            else:
                failures.extend(
                    validate_reconstruction(
                        row,
                        supports,
                        dim,
                        prob_lo,
                        spacing,
                        ghost_point,
                        hit_point,
                        normal,
                        args.condition_max,
                        args.moment_tol,
                        f"row {row_number}",
                    )
                )

        legacy_element = int(row["legacy_element"])
        legacy_available = legacy_element >= 0
        if "legacy_owner_available" in row:
            if int(row["legacy_owner_available"]) != int(legacy_available):
                failures.append(
                    f"row {row_number}: legacy-availability flag is inconsistent"
                )
        expected_match = int(legacy_available and hit_element == legacy_element)
        if int(row["legacy_matches_first_hit"]) != expected_match:
            failures.append(f"row {row_number}: legacy-match flag is inconsistent")
        if legacy_available:
            owner_differences += 1 - expected_match
        else:
            unavailable_legacy_owners += 1

        if edges is not None:
            independent = first_segment_hit(
                tuple(fluid_point), tuple(ghost_point), edges
            )
            if independent is None:
                failures.append(f"row {row_number}: independent query found no hit")
            else:
                independent_element, independent_fraction = independent
                if hit_geometry != 0 or hit_element != independent_element:
                    failures.append(
                        f"row {row_number}: C++ hit=({hit_geometry},{hit_element}) "
                        f"independent=(0,{independent_element})"
                    )
                if not math.isclose(
                    fraction,
                    independent_fraction,
                    rel_tol=2.0e-11,
                    abs_tol=2.0e-11,
                ):
                    failures.append(
                        f"row {row_number}: C++ fraction={fraction} "
                        f"independent={independent_fraction}"
                    )

    # C++ duplicate counters are deliberately rank-local.  Validate those
    # declarations per rank, then independently compare every global-key copy
    # across all supplied rank files.
    local_groups: dict[
        tuple[str, tuple[int, ...]], list[dict[str, str]]
    ] = defaultdict(list)
    for key, group in groups.items():
        for row in group:
            local_groups[(row["mpi_rank"], key)].append(row)
    for (rank, key), group in local_groups.items():
        reference = group[0]
        consistent = all(records_agree(reference, row, dim) for row in group[1:])
        reconstruction_consistent = True
        if not args.skip_supports:
            reconstruction_consistent = all(
                reconstruction_records_agree(
                    reference,
                    row,
                    supports_by_copy[reference["_copy_token"]],
                    supports_by_copy[row["_copy_token"]],
                    dim,
                )
                for row in group[1:]
            )
        available_legacy = [
            row["legacy_element"]
            for row in group
            if int(row["legacy_element"]) >= 0
        ]
        legacy_consistent = len(set(available_legacy)) <= 1
        for row in group:
            if int(row["duplicate_copies"]) != len(group):
                failures.append(
                    f"rank {rank} face {key}: duplicate count is inconsistent"
                )
                break
            if int(row["duplicate_first_hit_consistent"]) != int(consistent):
                failures.append(
                    f"rank {rank} face {key}: first-hit consistency flag is wrong"
                )
                break
            if "duplicate_reconstruction_consistent" in row and int(
                row["duplicate_reconstruction_consistent"]
            ) != int(reconstruction_consistent):
                failures.append(
                    f"rank {rank} face {key}: reconstruction consistency "
                    "flag is wrong"
                )
                break
            if int(row["duplicate_legacy_owner_consistent"]) != int(
                legacy_consistent
            ):
                failures.append(
                    f"rank {rank} face {key}: legacy consistency flag is wrong"
                )
                break

    duplicate_copies = 0
    for key, group in groups.items():
        duplicate_copies += len(group) - 1
        reference = group[0]
        consistent = all(records_agree(reference, row, dim) for row in group[1:])
        if not consistent:
            ranks = sorted({row["mpi_rank"] for row in group})
            failures.append(
                f"face {key}: global first-hit copies disagree across ranks {ranks}"
            )
        if not args.skip_supports:
            reconstruction_consistent = all(
                reconstruction_records_agree(
                    reference,
                    row,
                    supports_by_copy[reference["_copy_token"]],
                    supports_by_copy[row["_copy_token"]],
                    dim,
                )
                for row in group[1:]
            )
            if not reconstruction_consistent:
                ranks = sorted({row["mpi_rank"] for row in group})
                failures.append(
                    f"face {key}: global reconstruction copies disagree "
                    f"across ranks {ranks}"
                )

    if args.require_owner_difference and owner_differences == 0:
        failures.append("fixture did not distinguish first-hit from legacy owner")
    if failures:
        preview = "\n".join(failures[:24])
        raise RuntimeError(
            f"boundary-face validation failed with {len(failures)} errors:\n"
            f"{preview}"
        )

    print(
        "IBM boundary-face CSV: PASS "
        f"files={len(args.csv_files)} rows={len(rows)} unique={len(groups)} "
        f"duplicates={duplicate_copies} "
        f"owner_differences={owner_differences} "
        f"without_local_gp={unavailable_legacy_owners} "
        f"independent_2d={int(edges is not None)} "
        f"support_files={support_file_count}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
