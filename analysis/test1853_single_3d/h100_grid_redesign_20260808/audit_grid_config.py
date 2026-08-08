#!/usr/bin/env python3
"""Static safety and GPU-layout audit for the Test-1853 AMReX inputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


PROHIBITED_SWITCHES = (
    "cns.llf_ibm_face_local_crossing",
    "cns.llf_ibm_wall_matched_crossing",
    "cns.llf_ibm_conservative_crossing",
    "cns.llf_ibm_fluid_shell",
    "cns.llf_ibm_shifted_shell",
    "cns.afd_ibm_face_local_crossing",
    "cns.afd_ibm_crossing_wls",
    "cns.afd_ibm_conservative_crossing",
    "cns.afd_ibm_patch_conservative_crossing",
    "cns.afd_ibm_fluid_shell",
    "cns.afd_ibm_shifted_shell",
)


def parse_inputs(path: Path) -> dict[str, list[str]]:
    values: dict[str, list[str]] = {}
    for raw_line in path.read_text(encoding="ascii").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" in line:
            key, value = line.split("=", 1)
        else:
            key, value = line.split(maxsplit=1)
        values[key.strip()] = value.split()
    return values


def ints(values: dict[str, list[str]], key: str) -> list[int]:
    return [int(item) for item in values[key]]


def floats(values: dict[str, list[str]], key: str) -> list[float]:
    return [float(item) for item in values[key]]


def level_values(raw: list[int], count: int) -> list[int]:
    if len(raw) == 1:
        return raw * count
    if len(raw) != count:
        raise ValueError(f"expected 1 or {count} values, got {len(raw)}")
    return raw


def audit(path: Path) -> dict[str, object]:
    data = parse_inputs(path)
    max_level = int(data["amr.max_level"][0])
    level_count = max_level + 1
    n_cell = ints(data, "amr.n_cell")
    prob_lo = floats(data, "geometry.prob_lo")
    prob_hi = floats(data, "geometry.prob_hi")
    blocking = level_values(ints(data, "amr.blocking_factor"), level_count)
    max_grid = level_values(ints(data, "amr.max_grid_size"), level_count)
    error_buf = level_values(ints(data, "amr.n_error_buf"), max_level)
    ref_ratio = level_values(ints(data, "amr.ref_ratio"), max_level)

    dx0 = [(hi - lo) / n for lo, hi, n in zip(prob_lo, prob_hi, n_cell)]
    level_dx: list[list[float]] = [dx0]
    cumulative_ratio = 1
    for ratio in ref_ratio:
        cumulative_ratio *= ratio
        level_dx.append([dx / cumulative_ratio for dx in dx0])

    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": passed, "detail": detail})

    check(
        "cartesian_coordinate_system",
        data.get("geometry.coord_sys") == ["0"],
        f"geometry.coord_sys={data.get('geometry.coord_sys')}",
    )
    check(
        "base_domain_divisible_by_l0_blocking",
        all(n % blocking[0] == 0 for n in n_cell),
        f"n_cell={n_cell}, L0 blocking={blocking[0]}",
    )
    check(
        "max_grid_multiple_of_blocking",
        all(size % block == 0 for size, block in zip(max_grid, blocking)),
        f"max_grid_size={max_grid}, blocking_factor={blocking}",
    )
    check(
        "sfc_distribution",
        data.get("DistributionMapping.strategy") == ["SFC"],
        f"strategy={data.get('DistributionMapping.strategy')}",
    )
    check(
        "host_staged_mpi",
        data.get("amrex.use_gpu_aware_mpi") == ["0"],
        f"amrex.use_gpu_aware_mpi={data.get('amrex.use_gpu_aware_mpi')}",
    )
    check(
        "unmanaged_arena",
        data.get("amrex.the_arena_is_managed") == ["0"],
        f"amrex.the_arena_is_managed={data.get('amrex.the_arena_is_managed')}",
    )

    prohibited_active = [
        key for key in PROHIBITED_SWITCHES if data.get(key) != ["0"]
    ]
    check(
        "pure_gp_prohibited_switches_disabled",
        not prohibited_active,
        "all disabled" if not prohibited_active else str(prohibited_active),
    )

    # L2 is tagged only through x=-10 mm.  Error buffering occurs on its L1
    # parent.  Block alignment can extend the L2 valid region by at most one
    # L2 blocking unit minus one cell.  The STL minimum x is measured directly
    # from test1853_run165_single_master_SI.stl.
    l2_tag_x_hi = -0.010
    body_x_min = 0.00244127
    l2_parent_buffer = error_buf[1] * level_dx[1][0]
    l2_alignment_extension = (blocking[2] - 1) * level_dx[2][0]
    l2_worst_x_hi = l2_tag_x_hi + l2_parent_buffer + l2_alignment_extension
    l2_body_clearance = body_x_min - l2_worst_x_hi
    l2_clearance_cells = l2_body_clearance / level_dx[2][0]
    check(
        "l2_worst_case_stays_ahead_of_body",
        l2_body_clearance > 0.0,
        f"worst L2 x_hi={l2_worst_x_hi:.9f} m, "
        f"body x_min={body_x_min:.9f} m, clearance={l2_body_clearance:.9f} m",
    )
    check(
        "l2_clearance_at_least_four_cells",
        l2_clearance_cells >= 4.0,
        f"clearance={l2_clearance_cells:.3f} L2 cells",
    )

    l0_splits = [math.ceil(n / max_grid[0]) for n in n_cell]
    l0_box_lower_bound = math.prod(l0_splits)
    l0_cells = math.prod(n_cell)
    report: dict[str, object] = {
        "input": str(path),
        "pass": all(bool(item["passed"]) for item in checks),
        "n_cell": n_cell,
        "domain_lo_m": prob_lo,
        "domain_hi_m": prob_hi,
        "dx_by_level_m": level_dx,
        "blocking_factor": blocking,
        "max_grid_size": max_grid,
        "n_error_buf": error_buf,
        "regrid_int": ints(data, "amr.regrid_int"),
        "l0_cells": l0_cells,
        "l0_minimum_box_count": l0_box_lower_bound,
        "l0_average_cells_per_box_at_minimum": l0_cells / l0_box_lower_bound,
        "l2_worst_case_x_hi_m": l2_worst_x_hi,
        "l2_body_clearance_m": l2_body_clearance,
        "l2_body_clearance_cells": l2_clearance_cells,
        "checks": checks,
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    reports = [audit(path) for path in args.inputs]
    payload = {"pass": all(report["pass"] for report in reports), "reports": reports}
    rendered = json.dumps(payload, indent=2)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="ascii")
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

