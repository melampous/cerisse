#!/usr/bin/env python3
"""Parse the deterministic Run-165 regrid-prolongation failure audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


BAD_CELL_RE = re.compile(
    r"^\[IBM-FE-INPUT-BAD-CELL #(?P<id>\d+)\] rank=(?P<rank>\d+) "
    r"level=(?P<level>\d+) index=\((?P<index>[^)]+)\) "
    r"coord=\((?P<coord>[^)]+)\).* coarse_fine=(?P<coarse_fine>YES|no) "
    r"covered_by_finer=(?P<covered>YES|no).* provenance=(?P<provenance>\S+) "
    r"parent_index=\((?P<parent>[^)]+)\)$"
)
STATE_RE = re.compile(
    r"U=\[(?P<U>[^]]+)\] rho=(?P<rho>\S+) rhoe=(?P<rhoe>\S+) "
    r"p=(?P<p>\S+) T=(?P<T>\S+)$"
)
CHILD_RE = re.compile(
    r"^  fine_child\[offset=\((?P<offset>[^)]+)\)\] "
    r"index=\((?P<index>[^)]+)\) provenance=(?P<provenance>\S+) "
)
SUMMARY_RE = re.compile(
    r"^\[IBM-FE-INPUT-BAD-CELLS\].* bad_total=(?P<total>\d+) "
    r"bad_new_from_coarse=(?P<new>\d+) "
    r"bad_old_fine_overlap=(?P<overlap>\d+) "
    r"bad_unknown_origin=(?P<unknown>\d+)"
)
CONTEXT_RE = re.compile(
    r"^\[AMR-Regrid-Prolongation-Context\] level=(?P<level>\d+).* "
    r"fine_valid=(?P<valid>\d+) old_overlap=(?P<overlap>\d+) "
    r"new_from_coarse=(?P<new>\d+)"
)


def integer_tuple(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(","))


def real_tuple(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(","))


def parse_state(line: str) -> dict[str, object]:
    match = STATE_RE.search(line)
    if match is None:
        raise ValueError(f"cannot parse state line: {line}")
    values = real_tuple(match.group("U"))
    if len(values) != 5:
        raise ValueError(f"expected five conservative components: {line}")
    return {
        "U": values,
        "mom_x": values[0],
        "mom_y": values[1],
        "mom_z": values[2],
        "rho_E": values[3],
        "rho": float(match.group("rho")),
        "rhoe": float(match.group("rhoe")),
        "pressure": float(match.group("p")),
        "temperature": float(match.group("T")),
    }


def rhoe(state: tuple[float, ...]) -> float:
    mx, my, mz, rho_E, density = state
    return rho_E - 0.5 * (mx * mx + my * my + mz * mz) / density


def admissible_theta(
    parent: tuple[float, ...], child: tuple[float, ...], rhoe_floor: float
) -> float:
    def state_at(theta: float) -> tuple[float, ...]:
        return tuple(p + theta * (c - p) for p, c in zip(parent, child))

    high = state_at(1.0)
    if high[4] > 1.0e-19 and rhoe(high) > rhoe_floor:
        return 1.0
    if not (parent[4] > 1.0e-19 and rhoe(parent) > rhoe_floor):
        raise ValueError("coarse parent is not admissible")
    lower = 0.0
    upper = 1.0
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        candidate = state_at(midpoint)
        if candidate[4] > 1.0e-19 and rhoe(candidate) > rhoe_floor:
            lower = midpoint
        else:
            upper = midpoint
    return lower


def parse_log(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    contexts: dict[int, dict[str, int]] = {}
    summary: dict[str, int] | None = None
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    for line in lines:
        context_match = CONTEXT_RE.match(line)
        if context_match:
            level = int(context_match.group("level"))
            contexts[level] = {
                "fine_valid": int(context_match.group("valid")),
                "old_overlap": int(context_match.group("overlap")),
                "new_from_coarse": int(context_match.group("new")),
            }
            continue

        summary_match = SUMMARY_RE.match(line)
        if summary_match:
            summary = {key: int(value) for key, value in summary_match.groupdict().items()}
            continue

        bad_match = BAD_CELL_RE.match(line)
        if bad_match:
            if current is not None:
                records.append(current)
            coord = real_tuple(bad_match.group("coord"))
            index = integer_tuple(bad_match.group("index"))
            parent_index = integer_tuple(bad_match.group("parent"))
            current = {
                "id": int(bad_match.group("id")),
                "rank": int(bad_match.group("rank")),
                "level": int(bad_match.group("level")),
                "index": index,
                "coord": coord,
                "radius_yz": math.hypot(coord[1], coord[2]),
                "azimuth_deg": math.degrees(math.atan2(coord[2], coord[1])) % 360.0,
                "coarse_fine_adjacent": bad_match.group("coarse_fine") == "YES",
                "covered_by_finer": bad_match.group("covered") == "YES",
                "provenance": bad_match.group("provenance"),
                "parent_index": parent_index,
                "child_offset": tuple(i - 2 * p for i, p in zip(index, parent_index)),
                "children": [],
            }
            continue

        if current is None:
            continue
        if line.startswith("  center:"):
            current["center"] = parse_state(line)
        elif line.startswith("  coarse_parent index="):
            current["parent"] = parse_state(line)
        else:
            child_match = CHILD_RE.match(line)
            if child_match:
                child = parse_state(line)
                child["offset"] = integer_tuple(child_match.group("offset"))
                child["index"] = integer_tuple(child_match.group("index"))
                child["provenance"] = child_match.group("provenance")
                current["children"].append(child)

    if current is not None:
        records.append(current)
    if summary is None:
        raise ValueError("missing bad-cell summary")
    if len(records) != summary["total"]:
        raise ValueError(f"parsed {len(records)} records, expected {summary['total']}")

    rhoe_floor = 1.0e-8 / 0.4
    maximum_conservation_residual = 0.0
    for record in records:
        parent = record["parent"]["U"]
        children = record["children"]
        if len(children) != 8:
            raise ValueError(f"bad cell {record['id']} has {len(children)} child states")
        child_thetas = [admissible_theta(parent, child["U"], rhoe_floor) for child in children]
        record["theta_parent_max"] = min(child_thetas)
        average = tuple(sum(child["U"][n] for child in children) / 8.0 for n in range(5))
        residual = tuple(avg - coarse for avg, coarse in zip(average, parent))
        record["child_average_minus_parent"] = residual
        maximum_conservation_residual = max(
            maximum_conservation_residual, max(abs(value) for value in residual)
        )

    return {
        "source_log": str(path),
        "rhoe_floor": rhoe_floor,
        "contexts": contexts,
        "bad_cell_summary": summary,
        "maximum_child_average_conservation_residual": maximum_conservation_residual,
        "radius_yz_range": [
            min(record["radius_yz"] for record in records),
            max(record["radius_yz"] for record in records),
        ],
        "theta_parent_max_range": [
            min(record["theta_parent_max"] for record in records),
            max(record["theta_parent_max"] for record in records),
        ],
        "records": records,
    }


def write_csv(path: Path, result: dict[str, object]) -> None:
    fieldnames = [
        "id", "rank", "i", "j", "k", "x_m", "y_m", "z_m", "radius_yz_m",
        "azimuth_deg", "parent_i", "parent_j", "parent_k", "child_offset",
        "coarse_fine_adjacent", "provenance", "rho", "rhoe", "pressure_Pa",
        "temperature_K", "parent_rho", "parent_rhoe", "parent_pressure_Pa",
        "parent_temperature_K", "theta_parent_max",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for record in result["records"]:
            center = record["center"]
            parent = record["parent"]
            writer.writerow({
                "id": record["id"],
                "rank": record["rank"],
                "i": record["index"][0],
                "j": record["index"][1],
                "k": record["index"][2],
                "x_m": record["coord"][0],
                "y_m": record["coord"][1],
                "z_m": record["coord"][2],
                "radius_yz_m": record["radius_yz"],
                "azimuth_deg": record["azimuth_deg"],
                "parent_i": record["parent_index"][0],
                "parent_j": record["parent_index"][1],
                "parent_k": record["parent_index"][2],
                "child_offset": "".join(str(value) for value in record["child_offset"]),
                "coarse_fine_adjacent": record["coarse_fine_adjacent"],
                "provenance": record["provenance"],
                "rho": center["rho"],
                "rhoe": center["rhoe"],
                "pressure_Pa": center["pressure"],
                "temperature_K": center["temperature"],
                "parent_rho": parent["rho"],
                "parent_rhoe": parent["rhoe"],
                "parent_pressure_Pa": parent["pressure"],
                "parent_temperature_K": parent["temperature"],
                "theta_parent_max": record["theta_parent_max"],
            })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    result = parse_log(args.log)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_csv(args.csv, result)
    print(json.dumps({key: value for key, value in result.items() if key != "records"}, indent=2))


if __name__ == "__main__":
    main()
