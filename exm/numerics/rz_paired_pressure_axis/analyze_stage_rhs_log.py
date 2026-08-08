#!/usr/bin/env python3
"""Summarize the default-off R-Z per-stage RHS diagnostic log."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_record(line: str) -> tuple[str, dict[str, str]] | None:
    if line.startswith("[RZ-STAGE-RHS]"):
        kind = "rhs"
    elif line.startswith("[RZ-STAGE-RMOM-SPLIT]"):
        kind = "radial_momentum_split"
    else:
        return None
    fields: dict[str, str] = {}
    for token in line.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return kind, fields


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--tolerance", type=float, default=5.0e-13)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    records: list[dict[str, object]] = []
    for line in args.log.read_text(encoding="utf-8", errors="replace").splitlines():
        parsed = parse_record(line)
        if parsed is None:
            continue
        kind, fields = parsed
        records.append(
            {
                "kind": kind,
                "label": fields.get("label", "unknown"),
                "ring_i": int(fields["ring_i"]),
                "values": {
                    key: float(value)
                    for key, value in fields.items()
                    if key.endswith("_absmax")
                },
            }
        )

    if not records:
        raise SystemExit("No [RZ-STAGE-*] records found")

    rhs_keys = ("Rrho_absmax", "Rrho_ur_absmax", "Rrho_uz_absmax", "RE_absmax")
    rhs_records = [record for record in records if record["kind"] == "rhs"]
    split_records = [
        record for record in records if record["kind"] == "radial_momentum_split"
    ]
    worst_rhs = max(
        float(record["values"].get(key, 0.0))
        for record in rhs_records
        for key in rhs_keys
    )
    worst_closure = max(
        float(record["values"].get("Rpost_flux_remainder_absmax", 0.0))
        for record in split_records
    )
    labels = sorted({str(record["label"]) for record in rhs_records})
    rings = sorted({int(record["ring_i"]) for record in rhs_records})
    result = {
        "log": str(args.log),
        "tolerance": args.tolerance,
        "labels": labels,
        "rings": rings,
        "rhs_record_count": len(rhs_records),
        "split_record_count": len(split_records),
        "worst_conserved_rhs_absmax": worst_rhs,
        "worst_split_closure_absmax": worst_closure,
        "pass": worst_rhs <= args.tolerance and worst_closure <= args.tolerance,
        "records": records,
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(encoded + "\n", encoding="utf-8")
    raise SystemExit(0 if result["pass"] else 1)


if __name__ == "__main__":
    main()
