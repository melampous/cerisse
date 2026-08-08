#!/usr/bin/env python3
"""Summarize Test-D timing, hierarchy, support audits, and failure markers."""

from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path


STEP_TIME = re.compile(r"\[STEP (\d+)\] Coarse TimeStep time: ([0-9.eE+-]+)")
STEP_STATE = re.compile(
    r"^STEP = (\d+) TIME = ([0-9.eE+-]+) DT = ([0-9.eE+-]+)$"
)
REGRID_BASE = re.compile(r"Now regridding at level lbase = (\d+)")
LEVEL_GRID = re.compile(r"Level (\d+)\s+(\d+) grids\s+(\d+) cells")
ADVANCED_CELLS = re.compile(r"\[Level (\d+) step \d+\] Advanced (\d+) cells")
SUPPORT = re.compile(
    r"\[IBM-AMR-Support\] (GP|surface) level=(\d+) targets=(\d+) "
    r"active-supports=(\d+) physical-BC-supports=(\d+) "
    r"coarse-fine-missing-targets=(\d+) coarse-fine-missing-supports=(\d+)"
)

FAILURE_PATTERNS = {
    "state_check_failed": re.compile(r"STATE CHECK FAILED"),
    "inadmissible_before": re.compile(r"INADMISSIBLE_BEFORE"),
    "soft_positivity": re.compile(r"\[soft_positivity\]"),
    "nan_token": re.compile(r"(?<![A-Za-z])nan(?![A-Za-z])", re.IGNORECASE),
    "inf_token": re.compile(r"(?<![A-Za-z])[-+]?inf(?![A-Za-z])", re.IGNORECASE),
}


def distribution(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"count": 0, "mean_s": None, "median_s": None, "min_s": None, "max_s": None}
    ordered = sorted(values)
    return {
        "count": len(values),
        "mean_s": statistics.fmean(values),
        "median_s": statistics.median(values),
        "min_s": ordered[0],
        "max_s": ordered[-1],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tail-steps", type=int, default=200)
    args = parser.parse_args()

    timings: list[dict[str, object]] = []
    pending_regrid: set[int] = set()
    last_state: dict[str, float | int] | None = None
    hierarchy: dict[str, dict[str, int | None]] = {}
    support: dict[str, dict[str, int]] = {}
    failure_counts = {name: 0 for name in FAILURE_PATTERNS}

    for line in args.log.read_text(encoding="utf-8", errors="replace").splitlines():
        regrid_match = REGRID_BASE.search(line)
        if regrid_match:
            pending_regrid.add(int(regrid_match.group(1)))

        state_match = STEP_STATE.match(line)
        if state_match:
            last_state = {
                "step": int(state_match.group(1)),
                "time_s": float(state_match.group(2)),
                "dt_s": float(state_match.group(3)),
            }

        timing_match = STEP_TIME.search(line)
        if timing_match:
            timings.append(
                {
                    "step": int(timing_match.group(1)),
                    "wall_s": float(timing_match.group(2)),
                    "regrid_lbase": sorted(pending_regrid),
                }
            )
            pending_regrid.clear()

        level_match = LEVEL_GRID.search(line)
        if level_match:
            level = level_match.group(1)
            hierarchy[level] = {
                "grids": int(level_match.group(2)),
                "cells": int(level_match.group(3)),
            }

        advanced_match = ADVANCED_CELLS.search(line)
        if advanced_match:
            level = advanced_match.group(1)
            hierarchy.setdefault(level, {"grids": None, "cells": 0})
            hierarchy[level]["cells"] = int(advanced_match.group(2))

        support_match = SUPPORT.search(line)
        if support_match:
            kind, level = support_match.group(1), support_match.group(2)
            support[f"{kind}_L{level}"] = {
                "targets": int(support_match.group(3)),
                "active_supports": int(support_match.group(4)),
                "physical_bc_supports": int(support_match.group(5)),
                "coarse_fine_missing_targets": int(support_match.group(6)),
                "coarse_fine_missing_supports": int(support_match.group(7)),
            }

        for name, pattern in FAILURE_PATTERNS.items():
            if pattern.search(line):
                failure_counts[name] += 1

    tail = timings[-args.tail_steps :]
    tail_regular = [float(item["wall_s"]) for item in tail if not item["regrid_lbase"]]
    tail_regrid = [float(item["wall_s"]) for item in tail if item["regrid_lbase"]]
    by_lbase = {
        str(level): distribution(
            [
                float(item["wall_s"])
                for item in tail
                if level in item["regrid_lbase"]
            ]
        )
        for level in (0, 1)
    }

    valid_cells = sum(int(item["cells"]) for item in hierarchy.values())
    weighted_cell_updates = sum(
        int(item["cells"]) * (2 ** int(level))
        for level, item in hierarchy.items()
    )
    report = {
        "log": str(args.log),
        "last_state": last_state,
        "failure_counts": failure_counts,
        "hierarchy_last_seen": hierarchy,
        "valid_cells_last_seen": valid_cells,
        "subcycled_cell_updates_per_coarse_step": weighted_cell_updates,
        "ssprk43_cell_stages_per_coarse_step": 4 * weighted_cell_updates,
        "support_last_seen": support,
        "timing_all": distribution([float(item["wall_s"]) for item in timings]),
        "timing_tail": distribution([float(item["wall_s"]) for item in tail]),
        "timing_tail_regular": distribution(tail_regular),
        "timing_tail_any_regrid": distribution(tail_regrid),
        "timing_tail_by_regrid_lbase": by_lbase,
        "tail_step_count": len(tail),
    }
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="ascii")

    support_failed = any(
        item["coarse_fine_missing_targets"] != 0
        or item["coarse_fine_missing_supports"] != 0
        for item in support.values()
    )
    return 1 if any(failure_counts.values()) or support_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
