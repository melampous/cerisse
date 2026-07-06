#!/usr/bin/env python3
"""Run or prepare cylinder lateral-boundary/domain-width sensitivity checks.

This is intentionally separate from the native acoustic validation runner.  The
cylinder case is a production IBM restart test and requires a mature checkpoint;
when the checkpoint is not present the script emits a reproducible command
matrix instead of inventing a result.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "exm/ibm/validation/canonical/cylinder_bctest"
EXE = CASE_DIR / "main2d.gnu.MPI.ex"
DEFAULT_RESTART = Path(
    "/shared/cerisse/exm/ibm/validation/canonical/"
    "cylinder_M1p7_Re1e5/plot_l3_suth/chk295000"
)


@dataclass(frozen=True)
class BcSpec:
    name: str
    lo_bc: str
    hi_bc: str


BCS = {
    "char": BcSpec("char", lo_bc="1 7", hi_bc="7 7"),
    "foextrap": BcSpec("foextrap", lo_bc="1 2", hi_bc="7 2"),
}


def run_command(cmd: list[str], cwd: Path, stdout: Path, timeout: int) -> None:
    with stdout.open("w") as out:
        subprocess.run(
            cmd,
            cwd=cwd,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
        )


def build_case() -> None:
    run_command(["make", "-j4"], CASE_DIR, CASE_DIR / "build_cylinder_side.out", 1800)


def command_for(args: argparse.Namespace, bc: BcSpec, half_width: float,
                outdir: Path) -> list[str]:
    base_half_width = 5.0
    base_ny = 800
    ny = max(32, int(round(base_ny * half_width / base_half_width)))
    run_dir = outdir / f"{bc.name}_y{half_width:g}"
    plot_dir = run_dir / "plt"
    surf_dir = run_dir / "surf"
    chk_dir = run_dir / "chk"
    max_step = args.restart_step + args.steps

    return [
        "timeout", str(args.timeout),
        "mpirun", "--oversubscribe", "-np", str(args.np),
        str(EXE), args.inputs,
        f"amr.restart={args.restart}",
        f"max_step={max_step}",
        "stop_time=120.0",
        f"geometry.prob_lo=-4.5 {-half_width}",
        f"geometry.prob_hi=22.0 {half_width}",
        f"amr.n_cell=2120 {ny}",
        f"cns.lo_bc={bc.lo_bc}",
        f"cns.hi_bc={bc.hi_bc}",
        f"amr.plot_file={plot_dir}",
        f"ib.surf_file={surf_dir}",
        f"amr.check_file={chk_dir}",
        f"amr.plot_int={args.steps}",
        "amr.plot_files_output=1",
        "amr.checkpoint_files_output=0",
        "cns.nstep_screen_output=100",
    ]


def write_manifest(commands: list[tuple[str, float, list[str]]], outdir: Path,
                   executed: bool, reason: str) -> Path:
    report = outdir / "cylinder_side_domain_sensitivity.md"
    lines = [
        "# Cylinder Side-Domain BC Sensitivity",
        "",
        "Case: M=1.7 cylinder restart, lateral y-boundary comparison.",
        "BCs: `char` uses y code 7; `foextrap` uses y code 2.",
        "Widths: half-domain `|y| <= H`; x-domain and base dx are kept fixed.",
        "",
        f"Executed: `{executed}`",
    ]
    if reason:
        lines.append(f"Reason: {reason}")
    lines.extend(["", "Commands:", ""])
    for bc_name, half_width, cmd in commands:
        lines.append(f"## {bc_name}, H={half_width:g}")
        lines.append("")
        lines.append("```bash")
        lines.append(" ".join(cmd))
        lines.append("```")
        lines.append("")
    report.write_text("\n".join(lines))
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(ROOT / "temp/cylinder_side_domain_sensitivity"))
    parser.add_argument("--restart", default=str(DEFAULT_RESTART))
    parser.add_argument("--inputs", default="inputs.char")
    parser.add_argument("--bcs", default="char,foextrap")
    parser.add_argument("--half-widths", default="4,5,6")
    parser.add_argument("--restart-step", type=int, default=295000)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    bc_names = [x.strip() for x in args.bcs.split(",") if x.strip()]
    half_widths = [float(x) for x in args.half_widths.split(",") if x.strip()]
    unknown = sorted(set(bc_names) - set(BCS))
    if unknown:
        raise SystemExit(f"unknown BCs: {','.join(unknown)}")

    commands: list[tuple[str, float, list[str]]] = []
    for bc_name in bc_names:
        for half_width in half_widths:
            commands.append((bc_name, half_width,
                             command_for(args, BCS[bc_name], half_width, outdir)))

    restart = Path(args.restart)
    if args.dry_run or not restart.is_dir():
        reason = "dry run requested" if args.dry_run else f"restart missing: {restart}"
        report = write_manifest(commands, outdir, executed=False, reason=reason)
        print(f"wrote {report}")
        return

    if args.build or not EXE.exists():
        build_case()

    start = time.time()
    for bc_name, half_width, cmd in commands:
        run_dir = outdir / f"{bc_name}_y{half_width:g}"
        run_dir.mkdir(parents=True, exist_ok=True)
        print(f"running {bc_name} H={half_width:g}")
        run_command(cmd, CASE_DIR, run_dir / "run.out", args.timeout + 60)

    report = write_manifest(commands, outdir, executed=True, reason="")
    print(f"wrote {report}")
    print(f"elapsed {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
