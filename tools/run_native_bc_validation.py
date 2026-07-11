#!/usr/bin/env python3
"""Run native Cerisse BC validation cases with N/2N refinement.

The native case lives in exm/numerics/bc_native.  Each metric compares a short
domain against a longer-domain reference solved by the same Cerisse binary and
same grid spacing.  This validates the compiled FillPatch/bcnormal/WENO path,
not just the algebra used in a lightweight model.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("temp/.mplconfig")))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import yt


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "exm/numerics/bc_native"
EXE = CASE_DIR / "main2d.gnu.MPI.ex"

GAMMA = 1.4
RHO0 = 1.0
P0 = 1.0
C0 = math.sqrt(GAMMA * P0 / RHO0)
MACH = 0.30
U0 = MACH * C0
AMP = 1.0e-4
X0 = 0.30
Y0 = 0.50
SIGMA_X = 0.060
SIGMA_Y = 0.120
THETA_DEG = 30.0
THETA = math.radians(THETA_DEG)
LONG_LX = 2.50


@dataclass(frozen=True)
class CaseSpec:
    name: str
    mode: int
    ly: float
    stop_time: float
    base_n: int
    sigma_x: float = SIGMA_X
    sigma_y: float = SIGMA_Y
    carrier_y: int = 4
    theta_deg: float = THETA_DEG
    x0: float = X0
    y0: float = Y0
    amp: float = AMP


@dataclass(frozen=True)
class BcSpec:
    name: str
    hi_bc: str
    bc_mode: int
    lodi_transverse: float | None = None
    lodi_relax: float | None = None
    extras: tuple[str, ...] = ()


CASES = {
    "normal": CaseSpec("normal", mode=0, ly=0.125, stop_time=0.72, base_n=96),
    "oblique": CaseSpec("oblique", mode=1, ly=1.0, stop_time=1.05, base_n=96),
    "oblique_train": CaseSpec(
        "oblique_train", mode=2, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.18, sigma_y=1.0, carrier_y=4,
    ),
    "acoustic15": CaseSpec(
        "acoustic15", mode=2, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.18, sigma_y=1.0, carrier_y=3, theta_deg=15.0,
    ),
    "acoustic30": CaseSpec(
        "acoustic30", mode=2, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.18, sigma_y=1.0, carrier_y=4, theta_deg=30.0,
    ),
    "acoustic45": CaseSpec(
        "acoustic45", mode=2, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.18, sigma_y=1.0, carrier_y=4, theta_deg=45.0,
    ),
    "acoustic60": CaseSpec(
        "acoustic60", mode=2, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.18, sigma_y=1.0, carrier_y=5, theta_deg=60.0,
    ),
    "entropy": CaseSpec(
        "entropy", mode=3, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.16, sigma_y=1.0, carrier_y=4, x0=0.62,
    ),
    "vortex": CaseSpec(
        "vortex", mode=4, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.16, sigma_y=1.0, carrier_y=4, x0=0.62,
    ),
    "mixed": CaseSpec(
        "mixed", mode=5, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.16, sigma_y=1.0, carrier_y=4, x0=0.62,
    ),
    "shock_front": CaseSpec(
        "shock_front", mode=6, ly=0.5, stop_time=0.30, base_n=96,
        sigma_x=0.015, sigma_y=0.5, carrier_y=1, theta_deg=0.0,
        x0=0.55, amp=0.02 * C0,
    ),
    "wake": CaseSpec(
        "wake", mode=7, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.14, sigma_y=0.08, x0=0.62, y0=0.50, amp=0.005 * C0,
    ),
    "pressure_pulse": CaseSpec(
        "pressure_pulse", mode=8, ly=1.0, stop_time=1.05, base_n=96,
        sigma_x=0.18, sigma_y=0.18, x0=0.60, y0=0.50, amp=0.001 * C0,
    ),
}

BCS = {
    "char": BcSpec("char", hi_bc="7 -1", bc_mode=1),
    "persistent_lodi": BcSpec(
        "persistent_lodi", hi_bc="2 -1", bc_mode=0,
        extras=(
            "cns.nscbc_lo=0 0",
            "cns.nscbc_hi=2 0",
            "cns.nscbc_order=2",
            "cns.nscbc_use_transverse=1",
            "cns.nscbc_transverse_relax=0.25",
        ),
    ),
    "persistent_lodi_t0": BcSpec(
        "persistent_lodi_t0", hi_bc="2 -1", bc_mode=0,
        extras=(
            "cns.nscbc_lo=0 0", "cns.nscbc_hi=2 0",
            "cns.nscbc_order=2", "cns.nscbc_use_transverse=1",
            "cns.nscbc_transverse_relax=0.0",
        ),
    ),
    "persistent_lodi_t025": BcSpec(
        "persistent_lodi_t025", hi_bc="2 -1", bc_mode=0,
        extras=(
            "cns.nscbc_lo=0 0", "cns.nscbc_hi=2 0",
            "cns.nscbc_order=2", "cns.nscbc_use_transverse=1",
            "cns.nscbc_transverse_relax=0.25",
        ),
    ),
    "persistent_lodi_t05": BcSpec(
        "persistent_lodi_t05", hi_bc="2 -1", bc_mode=0,
        extras=(
            "cns.nscbc_lo=0 0", "cns.nscbc_hi=2 0",
            "cns.nscbc_order=2", "cns.nscbc_use_transverse=1",
            "cns.nscbc_transverse_relax=0.5",
        ),
    ),
    "lodi": BcSpec("lodi", hi_bc="7 -1", bc_mode=3),
    "rhs_lodi": BcSpec("rhs_lodi", hi_bc="7 -1", bc_mode=4),
    "rhs_lodi_t0": BcSpec(
        "rhs_lodi_t0", hi_bc="7 -1", bc_mode=4, lodi_transverse=0.0,
    ),
    "rhs_lodi_raw": BcSpec(
        "rhs_lodi_raw", hi_bc="7 -1", bc_mode=4,
        extras=(
            "prob.lodi_use_shock_sensor=0",
            "prob.lodi_use_convective_sensor=0",
            "prob.lodi_use_pressure_relax_sensor=0",
        ),
    ),
    "lodi_raw": BcSpec(
        "lodi_raw", hi_bc="7 -1", bc_mode=3,
        extras=(
            "prob.lodi_use_shock_sensor=0",
            "prob.lodi_use_convective_sensor=0",
            "prob.lodi_use_pressure_relax_sensor=0",
        ),
    ),
    "lodi_convective": BcSpec(
        "lodi_convective", hi_bc="7 -1", bc_mode=3,
        extras=("prob.lodi_use_convective_extrapolation=1",),
    ),
    "lodi_sponge": BcSpec(
        "lodi_sponge", hi_bc="7 -1", bc_mode=3,
        extras=(
            "prob.sponge_strength=3.0",
            "prob.sponge_width=0.20",
            "prob.sponge_power=2.0",
            "prob.sponge_x_hi=1",
        ),
    ),
    "rhs_lodi_sponge": BcSpec(
        "rhs_lodi_sponge", hi_bc="7 -1", bc_mode=4,
        extras=(
            "prob.sponge_strength=3.0",
            "prob.sponge_width=0.20",
            "prob.sponge_power=2.0",
            "prob.sponge_x_hi=1",
        ),
    ),
    "foextrap_sponge": BcSpec(
        "foextrap_sponge", hi_bc="2 -1", bc_mode=0,
        extras=(
            "prob.sponge_strength=3.0",
            "prob.sponge_width=0.20",
            "prob.sponge_power=2.0",
            "prob.sponge_x_hi=1",
        ),
    ),
    "lodi_amp": BcSpec("lodi_amp", hi_bc="7 -1", bc_mode=3),
    "lodi_t0": BcSpec("lodi_t0", hi_bc="7 -1", bc_mode=3, lodi_transverse=0.0),
    "lodi_t025": BcSpec("lodi_t025", hi_bc="7 -1", bc_mode=3, lodi_transverse=0.25),
    "lodi_t05": BcSpec("lodi_t05", hi_bc="7 -1", bc_mode=3, lodi_transverse=0.5),
    "lodi_t1": BcSpec("lodi_t1", hi_bc="7 -1", bc_mode=3, lodi_transverse=1.0),
    "plane": BcSpec("plane", hi_bc="7 -1", bc_mode=2),
    "lodi_plane": BcSpec("lodi_plane", hi_bc="7 -1", bc_mode=2),
    "foextrap": BcSpec("foextrap", hi_bc="2 -1", bc_mode=0),
    "dirichlet": BcSpec("dirichlet", hi_bc="1 -1", bc_mode=0),
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
    run_command(["make", "-j4"], CASE_DIR, CASE_DIR / "build_native_bc_validation.out", 900)


def latest_plot(run_dir: Path) -> Path:
    plots = sorted(p for p in run_dir.glob("plt*") if (p / "Header").is_file())
    if not plots:
        raise FileNotFoundError(f"no plotfile under {run_dir}")
    return plots[-1]


def mesh_for(case: CaseSpec, nx: int, lx: float) -> tuple[int, int]:
    dx = 1.0 / nx
    ny = max(4, int(round(case.ly / dx)))
    nx_domain = int(round(lx / dx))
    return nx_domain, ny


def run_cerisse(case: CaseSpec, bc: BcSpec, nx: int, lx: float,
                outdir: Path, label: str, timeout: int) -> Path:
    nx_domain, ny = mesh_for(case, nx, lx)
    run_dir = outdir / "runs" / f"{case.name}_{bc.name}_n{nx}_{label}"
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "timeout", str(timeout),
        "mpirun", "--oversubscribe", "-np", "1",
        str(EXE), "inputs",
        f"stop_time={case.stop_time}",
        "max_step=1000000",
        f"amr.n_cell={nx_domain} {ny}",
        f"geometry.prob_hi={lx} {case.ly}",
        "geometry.prob_lo=0.0 0.0",
        "geometry.is_periodic=0 1",
        "amr.max_grid_size=512",
        "amr.blocking_factor=2",
        f"amr.plot_file={run_dir / 'plt'}",
        "amr.plot_int=1000000",
        "amr.plot_files_output=1",
        "amr.checkpoint_files_output=0",
        "cns.lo_bc=1 -1",
        f"cns.hi_bc={bc.hi_bc}",
        "cns.nstep_screen_output=1000000",
        f"prob.mode={case.mode}",
        f"prob.bc_mode={bc.bc_mode}",
        f"prob.mach={MACH}",
        f"prob.amp={case.amp}",
        f"prob.theta_deg={case.theta_deg}",
        f"prob.x0={case.x0}",
        f"prob.y0={case.y0}",
        f"prob.sigma_x={case.sigma_x}",
        f"prob.sigma_y={case.sigma_y}",
        f"prob.carrier_y={case.carrier_y}",
    ]
    if bc.lodi_transverse is not None:
        cmd.append(f"prob.lodi_transverse={bc.lodi_transverse}")
    if bc.lodi_relax is not None:
        cmd.append(f"prob.lodi_relax={bc.lodi_relax}")
    cmd.extend(bc.extras)
    run_command(cmd, CASE_DIR, run_dir / "run.out", timeout + 30)
    return latest_plot(run_dir)


def covering_fields(plot: Path) -> dict[str, np.ndarray]:
    ds = yt.load(str(plot))
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                          dims=ds.domain_dimensions)

    def arr(name: str) -> np.ndarray:
        return np.asarray(cg[("boxlib", name)]).squeeze()

    return {
        "rho": arr("Density"),
        "p": arr("pressure"),
        "u": arr("x_velocity"),
        "v": arr("y_velocity"),
        "time": float(ds.current_time),
    }


def initial_perturbations(case: CaseSpec, nx: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nx_domain, ny = mesh_for(case, nx, 1.0)
    dx = 1.0 / nx_domain
    dy = case.ly / ny
    x = (np.arange(nx_domain) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(x, y, indexing="ij")

    theta = math.radians(case.theta_deg)
    if case.mode == 0:
        envelope = np.exp(-0.5 * ((xx - case.x0) / case.sigma_x) ** 2)
        nxw = 1.0
        nyw = 0.0
    elif case.mode == 1:
        dx0 = xx - case.x0
        dy0 = yy - case.y0
        dy0 = np.where(dy0 > 0.5 * case.ly, dy0 - case.ly, dy0)
        dy0 = np.where(dy0 < -0.5 * case.ly, dy0 + case.ly, dy0)
        nxw = math.cos(theta)
        nyw = math.sin(theta)
        xi = dx0 * nxw + dy0 * nyw
        eta = -dx0 * nyw + dy0 * nxw
        envelope = np.exp(
            -0.5 * ((xi / case.sigma_x) ** 2 + (eta / case.sigma_y) ** 2)
        )
    else:
        nxw = math.cos(theta)
        nyw = math.sin(theta)
        ky = 2.0 * math.pi * case.carrier_y / case.ly
        tan_theta = max(math.tan(theta), 1.0e-12)
        kx = ky / tan_theta
        phase = kx * (xx - case.x0) + ky * (yy - case.y0)
        envelope = np.exp(-0.5 * ((xx - case.x0) / case.sigma_x) ** 2) * np.cos(phase)

    amp = case.amp
    eps = amp / C0
    drho = np.zeros_like(envelope)
    dp = np.zeros_like(envelope)
    du = np.zeros_like(envelope)
    dv = np.zeros_like(envelope)

    if case.mode <= 2:
        dvel = amp * envelope
        dp = RHO0 * C0 * dvel
        drho = dp / (C0 * C0)
        du = dvel * nxw
        dv = dvel * nyw
    elif case.mode == 3:
        drho = RHO0 * eps * envelope
    elif case.mode == 4:
        du = -amp * nyw * envelope
        dv = amp * nxw * envelope
    elif case.mode == 5:
        acoustic = 0.50 * amp * envelope
        dp = RHO0 * C0 * acoustic
        drho = dp / (C0 * C0) + 0.35 * RHO0 * eps * envelope
        du = acoustic * nxw - 0.35 * amp * nyw * envelope
        dv = acoustic * nyw + 0.35 * amp * nxw * envelope
    elif case.mode == 6:
        width = max(case.sigma_x, 2.0e-3)
        front = 0.5 * (1.0 - np.tanh((xx - case.x0) / width))
        dp = P0 * eps * front
        drho = dp / (C0 * C0)
        du = dp / (RHO0 * C0)
        dv = np.zeros_like(du)
    elif case.mode == 7:
        dx0 = (xx - case.x0) / case.sigma_x
        dy0 = yy - case.y0
        dy0 = np.where(dy0 > 0.5 * case.ly, dy0 - case.ly, dy0)
        dy0 = np.where(dy0 < -0.5 * case.ly, dy0 + case.ly, dy0)
        wake = np.exp(-0.5 * (dx0 * dx0 + (dy0 / case.sigma_y) ** 2))
        du = -amp * wake
        drho = 0.20 * RHO0 * eps * wake
    elif case.mode == 8:
        dx0 = (xx - case.x0) / case.sigma_x
        dy0 = yy - case.y0
        dy0 = np.where(dy0 > 0.5 * case.ly, dy0 - case.ly, dy0)
        dy0 = np.where(dy0 < -0.5 * case.ly, dy0 + case.ly, dy0)
        pulse = np.exp(-0.5 * (dx0 * dx0 + (dy0 / case.sigma_y) ** 2))
        dp = RHO0 * C0 * amp * pulse
        drho = dp / (C0 * C0)
        du = 0.25 * dp / (RHO0 * C0)
    return drho, du, dv, dp


def initial_energy(case: CaseSpec, nx: int, metric_xmax: float = 1.0) -> float:
    nx_domain, ny = mesh_for(case, nx, 1.0)
    dx = 1.0 / nx_domain
    dy = case.ly / ny
    drho, du, dv, dp = initial_perturbations(case, nx)
    nx_metric = max(1, min(nx_domain, int(round(metric_xmax * nx_domain))))
    drho = drho[:nx_metric, :]
    du = du[:nx_metric, :]
    dv = dv[:nx_metric, :]
    dp = dp[:nx_metric, :]
    edens = 0.5 * RHO0 * (du * du + dv * dv)
    edens += 0.5 * dp * dp / (GAMMA * P0)
    if case.mode >= 3:
        edens += 0.5 * C0 * C0 * drho * drho / RHO0
    return float(np.sum(edens) * dx * dy)


def perturbation_energy(fields: dict[str, np.ndarray], case: CaseSpec, nx: int,
                        metric_xmax: float = 1.0) -> float:
    nx_domain, ny = mesh_for(case, nx, 1.0)
    nx_metric = max(1, min(nx_domain, int(round(metric_xmax * nx_domain))))
    dx = 1.0 / nx_domain
    dy = case.ly / ny
    u = fields["u"][:nx_metric, :ny] - U0
    v = fields["v"][:nx_metric, :ny]
    p = fields["p"][:nx_metric, :ny] - P0
    rho = fields["rho"][:nx_metric, :ny] - RHO0
    edens = 0.5 * RHO0 * (u * u + v * v) + 0.5 * p * p / (GAMMA * P0)
    if case.mode >= 3:
        edens += 0.5 * C0 * C0 * rho * rho / RHO0
    return float(np.sum(edens) * dx * dy)


def compare_plots(test_plot: Path, ref_plot: Path, case: CaseSpec, nx: int,
                  metric_xmax: float = 1.0) -> dict[str, float]:
    test = covering_fields(test_plot)
    ref = covering_fields(ref_plot)
    nx_domain, ny = mesh_for(case, nx, 1.0)
    nx_metric = max(1, min(nx_domain, int(round(metric_xmax * nx_domain))))
    dx = 1.0 / nx_domain
    dy = case.ly / ny

    du = test["u"][:nx_metric, :ny] - ref["u"][:nx_metric, :ny]
    dv = test["v"][:nx_metric, :ny] - ref["v"][:nx_metric, :ny]
    dp = test["p"][:nx_metric, :ny] - ref["p"][:nx_metric, :ny]
    drho = test["rho"][:nx_metric, :ny] - ref["rho"][:nx_metric, :ny]
    err_edens = 0.5 * RHO0 * (du * du + dv * dv) + 0.5 * dp * dp / (GAMMA * P0)
    if case.mode >= 3:
        err_edens += 0.5 * C0 * C0 * drho * drho / RHO0
    err_energy = float(np.sum(err_edens) * dx * dy)
    rho_l2 = float(math.sqrt(np.sum(drho * drho) * dx * dy))

    init_e = initial_energy(case, nx, metric_xmax)
    return {
        "metric": math.sqrt(err_energy / init_e),
        "rho_l2": rho_l2,
        "test_residual": math.sqrt(perturbation_energy(test, case, nx, metric_xmax) / init_e),
        "ref_residual": math.sqrt(perturbation_energy(ref, case, nx, metric_xmax) / init_e),
        "test_time": test["time"],
        "ref_time": ref["time"],
    }


def plot_results(rows: list[dict[str, object]], outdir: Path) -> tuple[Path, Path]:
    csv_path = outdir / "native_bc_metrics.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case", "bc", "nx", "ny", "metric", "rho_l2",
                "metric_xmax", "test_residual", "ref_residual", "test_time", "ref_time",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    fig_path = outdir / "native_bc_convergence.png"
    case_names = sorted({str(r["case"]) for r in rows})
    ncols = max(1, len(case_names))
    fig, axes = plt.subplots(1, ncols, figsize=(5.3 * ncols, 4.2), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, case_name in zip(axes, case_names):
        case_rows = [r for r in rows if r["case"] == case_name]
        for bc_name in sorted({str(r["bc"]) for r in case_rows}):
            rr = sorted([r for r in case_rows if r["bc"] == bc_name],
                        key=lambda x: int(x["nx"]))
            ax.loglog(
                [int(r["nx"]) for r in rr],
                [float(r["metric"]) for r in rr],
                marker="o",
                label=bc_name,
            )
        ax.grid(True, which="both", alpha=0.3)
        ax.set_title(case_name)
        ax.set_xlabel("N in x over short domain")
        ax.set_ylabel(r"$\sqrt{E_{test-ref}/E_{initial}}$")
        if case_rows:
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=220)
    plt.close(fig)
    return csv_path, fig_path


def write_report(rows: list[dict[str, object]], csv_path: Path,
                 fig_path: Path, outdir: Path) -> tuple[Path, Path]:
    md = outdir / "native_bc_validation_report.md"
    pdf = outdir / "native_bc_validation_report.pdf"

    by_key: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        by_key.setdefault((str(row["case"]), str(row["bc"])), []).append(row)

    lines = [
        "# Native Cerisse BC validation",
        "",
        "Solver path: compiled Cerisse `FillPatch + bcnormal + WENO-Z5 + RK3/4`, no IBM, no AMR.",
        "Metric: short-domain final solution compared against a long-domain reference at identical grid spacing on `x in [0,1]`.",
        "",
        "| Case | BC | N | 2N | metric(N) | metric(2N) | rate | ref residual(2N) | verdict |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for key in sorted(by_key):
        rr = sorted(by_key[key], key=lambda x: int(x["nx"]))
        if len(rr) != 2:
            continue
        m0 = float(rr[0]["metric"])
        m1 = float(rr[1]["metric"])
        rate = math.log(m0 / m1, 2.0) if m0 > 0.0 and m1 > 0.0 else float("nan")
        ref_resid = float(rr[1]["ref_residual"])
        if key[0] == "normal" and m1 < 1.0e-5 and rate > 0.5:
            verdict = "pass"
        elif key[0].startswith("oblique") and m1 > 1.0e-2 and abs(rate) < 0.2:
            verdict = "residual floor"
        elif rate > 0.5:
            verdict = "decreasing"
        else:
            verdict = "not converged"
        lines.append(
            f"| {key[0]} | {key[1]} | {rr[0]['nx']} | {rr[1]['nx']} | "
            f"{m0:.3e} | {m1:.3e} | {rate:.2f} | {ref_resid:.3e} | {verdict} |"
        )

    present_cases = sorted({str(r["case"]) for r in rows})
    lines.extend(["", "Interpretation:", ""])
    if "normal" in present_cases:
        lines.append(
            "- The quasi-1D normal acoustic case checks the compiled solver path for a normally outgoing packet."
        )
    if "oblique" in present_cases:
        lines.append(
            "- The oblique blob case is intentionally harder: its finite transverse envelope contains slow/non-radiating content, so the long-domain reference residual is not expected to vanish."
        )
    if "oblique_train" in present_cases:
        lines.append(
            "- The oblique wave-train case is a cleaner acoustic-branch test: it is y-periodic and uses a carrier wave aligned with the intended propagation angle."
        )
    if any(c.startswith("acoustic") for c in present_cases):
        lines.append(
            "- The acoustic angle-sweep cases test clean outgoing wave trains at fixed propagation angles."
        )
    if "entropy" in present_cases:
        lines.append(
            "- The entropy case checks whether a density/contact perturbation is spuriously converted into sound at the boundary."
        )
    if "vortex" in present_cases:
        lines.append(
            "- The vortex case checks tangential/rotational content, which is not a pure acoustic branch."
        )
    if "mixed" in present_cases:
        lines.append(
            "- The mixed case combines acoustic, entropy, and vortical content and is closer to a real side-boundary signal."
        )
    if "wake" in present_cases:
        lines.append(
            "- The wake case is a slow velocity-deficit packet; foextrap or a sponge may be preferable if this dominates a lateral boundary."
        )
    if "shock_front" in present_cases:
        lines.append(
            "- The shock_front case is a weak compression front crossing the boundary; it is a robustness check, not a proof of shock non-reflection."
        )
    lines.extend([
        "- A decreasing N/2N metric is stronger evidence than a single small number; a flat metric indicates a boundary-model residual floor for this setup.",
        "- These native tests isolate boundary behavior; production cylinder side-boundary sensitivity still requires a restart-based flow case.",
        "",
        "Files:",
        "",
        f"- `{csv_path.name}`",
        f"- `{fig_path.name}`",
    ])
    md.write_text("\n".join(lines) + "\n")

    with PdfPages(pdf) as pages:
        fig = plt.figure(figsize=(8.5, 5.2))
        img = plt.imread(fig_path)
        ax = fig.add_axes([0.05, 0.08, 0.90, 0.84])
        ax.imshow(img)
        ax.axis("off")
        pages.savefig(fig)
        plt.close(fig)

        fig = plt.figure(figsize=(8.5, 11.0))
        ax = fig.add_axes([0.06, 0.06, 0.88, 0.88])
        ax.axis("off")
        ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left",
                fontsize=8, family="monospace")
        pages.savefig(fig)
        plt.close(fig)

    return md, pdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(ROOT / "temp/bc_native_validation"))
    parser.add_argument("--cases", default="normal,oblique")
    parser.add_argument("--bcs", default="char,foextrap,dirichlet")
    parser.add_argument("--res", default="")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument(
        "--metric-xmax", type=float, default=1.0,
        help="compare only x in [0,metric_xmax] of the short-domain length",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    if outdir.exists():
        if not args.overwrite:
            raise SystemExit(f"{outdir} exists; pass --overwrite to replace it")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    if args.build or not EXE.exists():
        build_case()

    case_names = [c.strip() for c in args.cases.split(",") if c.strip()]
    bc_names = [b.strip() for b in args.bcs.split(",") if b.strip()]
    unknown_cases = sorted(set(case_names) - set(CASES))
    unknown_bcs = sorted(set(bc_names) - set(BCS))
    if unknown_cases:
        raise SystemExit(f"unknown cases: {','.join(unknown_cases)}")
    if unknown_bcs:
        raise SystemExit(f"unknown BCs: {','.join(unknown_bcs)}")
    if args.res:
        res = [int(x) for x in args.res.split(",")]
    elif args.quick:
        res = [48, 96]
    else:
        res = [96, 192]
    if len(res) != 2 or res[1] != 2 * res[0]:
        raise SystemExit("--res must contain N,2N")

    rows: list[dict[str, object]] = []
    start = time.time()
    for case_name in case_names:
        case = CASES[case_name]
        for bc_name in bc_names:
            bc = BCS[bc_name]
            for nx in res:
                print(f"running {case.name:7s} {bc.name:9s} N={nx} short")
                test_plot = run_cerisse(case, bc, nx, 1.0, outdir, "short", args.timeout)
                print(f"running {case.name:7s} {bc.name:9s} N={nx} long")
                ref_plot = run_cerisse(case, bc, nx, LONG_LX, outdir, "long", args.timeout)
                metrics = compare_plots(test_plot, ref_plot, case, nx, args.metric_xmax)
                _, ny = mesh_for(case, nx, 1.0)
                row = {
                    "case": case.name,
                    "bc": bc.name,
                    "nx": nx,
                    "ny": ny,
                    "metric_xmax": args.metric_xmax,
                    **metrics,
                }
                rows.append(row)
                print(f"  metric={metrics['metric']:.3e}, "
                      f"test_residual={metrics['test_residual']:.3e}")

    csv_path, fig_path = plot_results(rows, outdir)
    md, pdf = write_report(rows, csv_path, fig_path, outdir)
    elapsed = time.time() - start
    print(f"wrote {csv_path}")
    print(f"wrote {md}")
    print(f"wrote {pdf}")
    print(f"elapsed {elapsed:.1f} s")


if __name__ == "__main__":
    main()
