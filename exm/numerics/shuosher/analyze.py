#!/usr/bin/env python3
"""Compare Shu--Osher density profiles with a fine-grid numerical reference."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np

try:
    import yt
except ImportError as exc:  # pragma: no cover - dependency diagnostic
    raise SystemExit("analyze.py requires yt (python3 -m pip install yt)") from exc


def latest_plotfile(path: Path) -> Path:
    if path.is_dir() and (path / "Header").is_file():
        return path
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        p for p in path.glob("plt*") if pattern.fullmatch(p.name) and (p / "Header").is_file()
    ]
    if not candidates:
        candidates = [
            p
            for p in path.glob("**/plt*")
            if pattern.fullmatch(p.name) and (p / "Header").is_file()
        ]
    if not candidates:
        raise FileNotFoundError(f"no AMReX plotfile found below {path}")
    return max(candidates, key=lambda candidate: (int(pattern.fullmatch(candidate.name)[1]), str(candidate)))


def field_key(ds, name: str):
    wanted = name.lower()
    for key in ds.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} is unavailable; fields are {ds.field_list}")


def read_profile(path: Path):
    plotfile = latest_plotfile(path)
    ds = yt.load(str(plotfile))
    data = ds.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    rho = np.asarray(data[field_key(ds, "Density")].d, dtype=float)
    order = np.argsort(x)
    return plotfile, ds, x[order], rho[order]


def require_matching_uniform_domains(case_ds, reference_ds, case_plotfile: Path) -> None:
    case_level = int(case_ds.index.max_level)
    reference_level = int(reference_ds.index.max_level)
    if case_level != 0 or reference_level != 0:
        raise ValueError(
            "Shu--Osher reference comparisons require uniform level-0 data; "
            f"got levels {case_level} and {reference_level} for {case_plotfile}"
        )

    case_lo = np.asarray(case_ds.domain_left_edge.d, dtype=float)
    case_hi = np.asarray(case_ds.domain_right_edge.d, dtype=float)
    reference_lo = np.asarray(reference_ds.domain_left_edge.d, dtype=float)
    reference_hi = np.asarray(reference_ds.domain_right_edge.d, dtype=float)
    scale = max(1.0, float(np.max(np.abs(np.concatenate((reference_lo, reference_hi))))))
    if not (
        np.allclose(case_lo, reference_lo, rtol=0.0, atol=1.0e-12 * scale)
        and np.allclose(case_hi, reference_hi, rtol=0.0, atol=1.0e-12 * scale)
    ):
        raise ValueError(
            f"domain mismatch: {case_plotfile} spans [{case_lo[0]}, {case_hi[0]}], "
            f"whereas the reference spans [{reference_lo[0]}, {reference_hi[0]}]"
        )


def advance_time(path: Path) -> float:
    logs = [path / "run.log", path.parent / "run.log"]
    pattern = re.compile(r"Run Time advance\s*=\s*([0-9.eE+-]+)")
    for log in logs:
        if not log.is_file():
            continue
        match = pattern.search(log.read_text(errors="replace"))
        if match:
            return float(match.group(1))
    return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path, help="N=6400 reference directory or plotfile")
    parser.add_argument("cases", nargs="+", type=Path, help="coarser result directories")
    parser.add_argument("--csv", type=Path, default=Path("shuosher_errors.csv"))
    parser.add_argument("--figure", type=Path, help="optional density-profile figure")
    args = parser.parse_args()

    reference_plotfile, ref_ds, x_ref, rho_ref = read_profile(args.reference)
    if int(ref_ds.index.max_level) != 0:
        raise ValueError(f"reference must be a uniform level-0 result: {reference_plotfile}")
    ref_time = float(ref_ds.current_time)
    rows = []
    profiles = []
    for case in args.cases:
        plotfile, ds, x, rho = read_profile(case)
        require_matching_uniform_domains(ds, ref_ds, plotfile)
        time = float(ds.current_time)
        if abs(time - ref_time) > 1.0e-10 * max(1.0, abs(ref_time)):
            raise ValueError(f"time mismatch: {plotfile} has t={time}, reference has t={ref_time}")
        if x_ref.size <= x.size:
            raise ValueError(
                f"reference must be finer than every case: Nref={x_ref.size}, "
                f"Ncase={x.size} for {plotfile}"
            )
        tolerance = 32.0 * np.finfo(float).eps * max(1.0, np.max(np.abs(x_ref)))
        if np.min(x) < np.min(x_ref) - tolerance or np.max(x) > np.max(x_ref) + tolerance:
            raise ValueError(
                f"case cell centres in {plotfile} extend beyond the reference-centre range; "
                "interpolation would require extrapolation"
            )
        rho_reference = np.interp(x, x_ref, rho_ref)
        error = rho - rho_reference
        rows.append(
            {
                "plotfile": str(plotfile),
                "cells": int(x.size),
                "time": time,
                "advance_seconds": advance_time(case),
                "density_L1": float(np.mean(np.abs(error))),
                "density_L2": float(np.sqrt(np.mean(error * error))),
                "density_Linf": float(np.max(np.abs(error))),
                "density_TV": float(np.sum(np.abs(np.diff(rho)))),
                "density_min": float(np.min(rho)),
                "density_max": float(np.max(rho)),
            }
        )
        profiles.append((x, rho, f"N={x.size}"))

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"N={row['cells']:5d}  L1={row['density_L1']:.6e}  "
            f"L2={row['density_L2']:.6e}  Linf={row['density_Linf']:.6e}  "
            f"advance={row['advance_seconds']:.6g} s"
        )
    print(f"wrote {args.csv}")

    if args.figure:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.0, 3.8))
        ax.plot(x_ref, rho_ref, color="black", linewidth=1.0, label="N=6400 reference")
        for x, rho, label in profiles:
            ax.plot(x, rho, linewidth=0.8, label=label)
        ax.set(xlabel=r"$x$", ylabel=r"$\rho$", xlim=(-5.0, 5.0))
        ax.legend(frameon=False)
        fig.tight_layout()
        args.figure.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.figure, dpi=300)
        print(f"wrote {args.figure}")


if __name__ == "__main__":
    yt.set_log_level(40)
    main()
