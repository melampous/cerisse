#!/usr/bin/env python3
"""Measure native Cerisse BC runtime for the acoustic validation case."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "exm/numerics/bc_native"
EXE = CASE_DIR / "main2d.gnu.MPI.ex"


@dataclass(frozen=True)
class BcSpec:
    name: str
    hi_bc: str
    bc_mode: int


BCS = {
    "char": BcSpec("char", hi_bc="7 -1", bc_mode=1),
    "lodi": BcSpec("lodi", hi_bc="7 -1", bc_mode=3),
    "plane": BcSpec("plane", hi_bc="7 -1", bc_mode=2),
    "foextrap": BcSpec("foextrap", hi_bc="2 -1", bc_mode=0),
}


def run_command(cmd: list[str], cwd: Path, stdout: Path, timeout: int) -> float:
    start = time.perf_counter()
    with stdout.open("w") as out:
        subprocess.run(
            cmd,
            cwd=cwd,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
        )
    return time.perf_counter() - start


def build_case() -> None:
    run_command(["make", "-j4"], CASE_DIR, CASE_DIR / "build_native_bc_performance.out", 900)


def parse_runtime(text: str, label: str) -> float:
    match = re.search(rf"Run Time {label}\s*=\s*([0-9.eE+-]+)", text)
    if not match:
        return float("nan")
    return float(match.group(1))


def run_case(bc: BcSpec, nx: int, ny: int, steps: int,
             outdir: Path, repeat: int, timeout: int) -> dict[str, object]:
    run_dir = outdir / "runs" / f"{bc.name}_r{repeat}"
    run_dir.mkdir(parents=True, exist_ok=True)
    stdout = run_dir / "run.out"
    cmd = [
        "timeout", str(timeout),
        "mpirun", "--oversubscribe", "-np", "1",
        str(EXE), "inputs",
        f"max_step={steps}",
        "stop_time=1000.0",
        f"amr.n_cell={nx} {ny}",
        "geometry.prob_lo=0.0 0.0",
        "geometry.prob_hi=1.0 1.0",
        "geometry.is_periodic=0 1",
        "amr.max_grid_size=512",
        "amr.blocking_factor=2",
        "amr.plot_files_output=0",
        "amr.checkpoint_files_output=0",
        "cns.lo_bc=1 -1",
        f"cns.hi_bc={bc.hi_bc}",
        "cns.nstep_screen_output=1000000",
        "prob.mode=2",
        f"prob.bc_mode={bc.bc_mode}",
        "prob.mach=0.30",
        "prob.amp=1.0e-4",
        "prob.theta_deg=30.0",
        "prob.x0=0.30",
        "prob.y0=0.50",
        "prob.sigma_x=0.18",
        "prob.sigma_y=1.0",
        "prob.carrier_y=4",
    ]
    wall = run_command(cmd, CASE_DIR, stdout, timeout + 30)
    text = stdout.read_text(errors="replace")
    return {
        "bc": bc.name,
        "repeat": repeat,
        "nx": nx,
        "ny": ny,
        "steps": steps,
        "wall_s": wall,
        "run_total_s": parse_runtime(text, "total"),
        "run_advance_s": parse_runtime(text, "advance"),
    }


def write_outputs(rows: list[dict[str, object]], outdir: Path) -> tuple[Path, Path]:
    csv_path = outdir / "native_bc_performance.csv"
    fieldnames = ["bc", "repeat", "nx", "ny", "steps", "wall_s", "run_total_s", "run_advance_s"]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    md_path = outdir / "native_bc_performance.md"
    by_bc: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_bc.setdefault(str(row["bc"]), []).append(row)

    lines = [
        "# Native BC performance",
        "",
        "Case: 2D oblique acoustic wave train, no plot output, one MPI rank.",
        "",
        "| BC | median advance [s] | median total [s] | median wall [s] | rel. to char |",
        "|---|---:|---:|---:|---:|",
    ]
    char_med = median(float(r["run_advance_s"]) for r in by_bc.get("char", []) if r["run_advance_s"] == r["run_advance_s"])
    for bc_name in sorted(by_bc):
        rr = by_bc[bc_name]
        adv = median(float(r["run_advance_s"]) for r in rr if r["run_advance_s"] == r["run_advance_s"])
        total = median(float(r["run_total_s"]) for r in rr if r["run_total_s"] == r["run_total_s"])
        wall = median(float(r["wall_s"]) for r in rr)
        rel = adv / char_med if char_med > 0.0 else float("nan")
        lines.append(f"| {bc_name} | {adv:.6g} | {total:.6g} | {wall:.6g} | {rel:.3f} |")
    lines.extend([
        "",
        f"Raw data: `{csv_path.name}`",
    ])
    md_path.write_text("\n".join(lines) + "\n")
    return csv_path, md_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(ROOT / "temp/bc_performance"))
    parser.add_argument("--bcs", default="char,lodi,foextrap")
    parser.add_argument("--nx", type=int, default=192)
    parser.add_argument("--ny", type=int, default=192)
    parser.add_argument("--steps", type=int, default=160)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--timeout", type=int, default=240)
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    if outdir.exists():
        if not args.overwrite:
            raise SystemExit(f"{outdir} exists; pass --overwrite to replace it")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    if args.build or not EXE.exists():
        build_case()

    bc_names = [name.strip() for name in args.bcs.split(",") if name.strip()]
    unknown = sorted(set(bc_names) - set(BCS))
    if unknown:
        raise SystemExit(f"unknown BCs: {','.join(unknown)}")

    rows: list[dict[str, object]] = []
    for repeat in range(args.repeats):
        for bc_name in bc_names:
            print(f"running {bc_name:9s} repeat={repeat + 1}/{args.repeats}")
            rows.append(
                run_case(BCS[bc_name], args.nx, args.ny, args.steps,
                         outdir, repeat + 1, args.timeout)
            )

    csv_path, md_path = write_outputs(rows, outdir)
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
