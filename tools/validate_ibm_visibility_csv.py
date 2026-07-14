#!/usr/bin/env python3
"""Validate a 2-D Cerisse IBM visibility-audit CSV independently."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from test_ibm_visibility_audit import first_segment_hit, polygon_edges


def read_polygon(path: Path) -> list[tuple[float, float]]:
    vertices = []
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values = stripped.split()
        vertices.append((float(values[0]), float(values[1])))
    if len(vertices) < 3:
        raise ValueError(f"{path} does not contain a closed 2-D polygon")
    return vertices


def parse_pair(values: list[float], name: str) -> np.ndarray:
    if len(values) != 2:
        raise ValueError(f"{name} requires exactly two values")
    return np.asarray(values, dtype=float)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--prob-lo", type=float, nargs=2, required=True)
    parser.add_argument("--prob-hi", type=float, nargs=2, required=True)
    parser.add_argument("--n-cell", type=int, nargs=2, required=True)
    parser.add_argument("--require-rejections", action="store_true")
    args = parser.parse_args()

    prob_lo = parse_pair(args.prob_lo, "--prob-lo")
    prob_hi = parse_pair(args.prob_hi, "--prob-hi")
    n_cell = np.asarray(args.n_cell, dtype=int)
    if np.any(n_cell <= 0) or np.any(prob_hi <= prob_lo):
        raise ValueError("invalid domain or cell counts")
    spacing = (prob_hi - prob_lo) / n_cell
    edges = polygon_edges(read_polygon(args.geometry))

    with args.csv_file.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"{args.csv_file} contains no candidate supports")

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    mismatches: list[str] = []
    rejected = 0
    for row_number, row in enumerate(rows, start=2):
        key = (row["kind"], row["target"], row["image"])
        grouped[key].append(row)
        support = int(row["support_slot"])
        candidate_bit = (int(row["candidate_mask"]) >> support) & 1
        visible_bit = (int(row["visible_mask"]) >> support) & 1
        candidate = int(row["candidate"])
        visible = int(row["visible"])
        if candidate != 1 or candidate_bit != candidate:
            mismatches.append(f"row {row_number}: candidate mask mismatch")
        if visible_bit != visible:
            mismatches.append(f"row {row_number}: visible mask mismatch")

        image = (float(row["image_x"]), float(row["image_y"]))
        support_index = np.array(
            [int(row["support_i"]), int(row["support_j"])], dtype=float
        )
        support_point = prob_lo + (support_index + 0.5) * spacing
        independent_hit = first_segment_hit(
            image, tuple(support_point), edges
        )
        independently_visible = independent_hit is None
        if bool(visible) != independently_visible:
            mismatches.append(
                f"row {row_number}: C++ visible={visible}, "
                f"independent hit={independent_hit}"
            )
        if visible:
            if int(row["first_hit_element"]) >= 0:
                mismatches.append(f"row {row_number}: clear segment has hit id")
        else:
            rejected += 1
            if int(row["first_hit_element"]) < 0:
                mismatches.append(f"row {row_number}: blocked segment lacks hit id")

    for key, group in grouped.items():
        expected_candidates = int(group[0]["n_candidates"])
        expected_visible = int(group[0]["n_visible"])
        if len(group) != expected_candidates:
            mismatches.append(
                f"{key}: rows={len(group)} != n_candidates={expected_candidates}"
            )
        if sum(int(row["visible"]) for row in group) != expected_visible:
            mismatches.append(f"{key}: n_visible is inconsistent with rows")
        if any(int(row["n_candidates"]) != expected_candidates for row in group):
            mismatches.append(f"{key}: n_candidates changes within one image")
        if any(int(row["n_visible"]) != expected_visible for row in group):
            mismatches.append(f"{key}: n_visible changes within one image")

    if args.require_rejections and rejected == 0:
        mismatches.append("fixture did not produce any rejected supports")
    if mismatches:
        preview = "\n".join(mismatches[:20])
        raise RuntimeError(
            f"visibility audit validation failed with {len(mismatches)} errors:\n"
            f"{preview}"
        )

    print(
        f"IBM visibility CSV: PASS rows={len(rows)} "
        f"images={len(grouped)} rejected={rejected}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
