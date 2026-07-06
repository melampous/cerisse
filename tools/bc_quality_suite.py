#!/usr/bin/env python3
"""Fast boundary-condition quality suite for open compressible-flow BCs.

The suite combines the existing 1D nonlinear HLLC acoustic reflection test with
three 2D linearized-Euler outflow tests:

1. normal acoustic pulse reflection,
2. oblique acoustic wave packet exiting the right boundary,
3. convected weak vortex exiting the right boundary,
4. convected entropy spot exiting the right boundary.

The 2D tests are deliberately light-weight.  They isolate boundary behaviour
without IBM geometry, AMR, or long transient statistics.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path("temp/bc_quality_suite/.mplconfig")),
)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bc_reflection_study as bc1d


GAMMA = 1.4
RHO0 = 1.0
P0 = 1.0
C0 = math.sqrt(GAMMA * P0 / RHO0)
MACH = 0.30
U0 = MACH * C0

ACTIVE_BCS = ["foextrap", "dirichlet", "char_farfield"]
LEGACY_BCS = ["foextrap", "dirichlet", "nscbc_pressure", "char_farfield"]
BC_LABELS = {
    "foextrap": "foextrap",
    "dirichlet": "hard Dirichlet",
    "nscbc_pressure": "old pressure-clamped NSCBC",
    "char_farfield": "fixed characteristic NSCBC",
}


@dataclass
class TestResult:
    case: str
    bc: str
    metric: float
    metric_name: str
    notes: str


def energy(q: np.ndarray, dx: float, dy: float) -> float:
    rho, u, v, p = q
    edens = 0.5 * RHO0 * (u * u + v * v) + 0.5 * p * p / (GAMMA * P0)
    return float(np.sum(edens) * dx * dy)


def flux_x(q: np.ndarray) -> np.ndarray:
    rho, u, v, p = q
    return np.stack(
        [
            U0 * rho + RHO0 * u,
            U0 * u + p / RHO0,
            U0 * v,
            U0 * p + GAMMA * P0 * u,
        ],
        axis=0,
    )


def flux_y(q: np.ndarray) -> np.ndarray:
    rho, u, v, p = q
    return np.stack(
        [
            RHO0 * v,
            np.zeros_like(u),
            p / RHO0,
            GAMMA * P0 * v,
        ],
        axis=0,
    )


def right_bc(qi: np.ndarray, bc: str) -> np.ndarray:
    if bc == "foextrap":
        return qi.copy()
    if bc == "dirichlet":
        return np.zeros_like(qi)

    rho_i, un_i, vt_i, p_i = qi
    entropy_i = p_i - C0 * C0 * rho_i
    jplus_i = un_i + p_i / (RHO0 * C0)

    if bc == "nscbc_pressure":
        p_g = np.zeros_like(p_i)
        un_g = jplus_i
        rho_g = (p_g - entropy_i) / (C0 * C0)
        return np.stack([rho_g, un_g, vt_i, p_g], axis=0)

    if bc == "char_farfield":
        jminus_inf = np.zeros_like(jplus_i)
        un_g = 0.5 * (jplus_i + jminus_inf)
        p_g = 0.5 * RHO0 * C0 * (jplus_i - jminus_inf)
        rho_g = (p_g - entropy_i) / (C0 * C0)
        return np.stack([rho_g, un_g, vt_i, p_g], axis=0)

    raise ValueError(f"unknown BC: {bc}")


def rhs(q: np.ndarray, dx: float, dy: float, bc: str) -> np.ndarray:
    # q shape: (4, nx, ny)
    qx = np.empty((4, q.shape[1] + 2, q.shape[2]))
    qx[:, 1:-1, :] = q
    qx[:, 0, :] = np.zeros_like(q[:, 0, :])
    qx[:, -1, :] = right_bc(q[:, -1, :], bc)

    qy = np.concatenate([q[:, :, -1:], q, q[:, :, :1]], axis=2)

    ax = abs(U0) + C0
    ay = C0
    fx_l = flux_x(qx[:, :-1, :])
    fx_r = flux_x(qx[:, 1:, :])
    fhat_x = 0.5 * (fx_l + fx_r) - 0.5 * ax * (qx[:, 1:, :] - qx[:, :-1, :])

    fy_l = flux_y(qy[:, :, :-1])
    fy_r = flux_y(qy[:, :, 1:])
    fhat_y = 0.5 * (fy_l + fy_r) - 0.5 * ay * (qy[:, :, 1:] - qy[:, :, :-1])

    return -((fhat_x[:, 1:, :] - fhat_x[:, :-1, :]) / dx +
             (fhat_y[:, :, 1:] - fhat_y[:, :, :-1]) / dy)


def rk3_step(q: np.ndarray, dt: float, dx: float, dy: float, bc: str) -> np.ndarray:
    q1 = q + dt * rhs(q, dx, dy, bc)
    q2 = 0.75 * q + 0.25 * (q1 + dt * rhs(q1, dx, dy, bc))
    return (q + 2.0 * (q2 + dt * rhs(q2, dx, dy, bc))) / 3.0


def run_2d_case(kind: str, bc: str, nx: int, ny: int, quick: bool) -> tuple[TestResult, dict[str, np.ndarray]]:
    lx, ly = 1.0, 0.5
    dx, dy = lx / nx, ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(x, y, indexing="ij")
    lx_ref = 1.65
    nx_ref = int(math.ceil(lx_ref / dx))
    x_ref = (np.arange(nx_ref) + 0.5) * dx
    xx_ref, yy_ref = np.meshgrid(x_ref, y, indexing="ij")

    amp = 1.0e-4
    x0, y0 = 0.25, 0.25
    sigma = 0.045 if quick else 0.055

    def initial_state(xm: np.ndarray, ym: np.ndarray) -> np.ndarray:
        rr2 = (xm - x0) ** 2 + (ym - y0) ** 2
        g = np.exp(-0.5 * rr2 / (sigma * sigma))
        q0 = np.zeros((4, xm.shape[0], xm.shape[1]))
        if kind == "oblique_acoustic":
            theta = math.radians(30.0)
            nxw, nyw = math.cos(theta), math.sin(theta)
            p = amp * g
            q0[0] = p / (C0 * C0)
            q0[1] = nxw * p / (RHO0 * C0)
            q0[2] = nyw * p / (RHO0 * C0)
            q0[3] = p
        elif kind == "weak_vortex":
            psi = amp * g
            q0[1] = -(ym - y0) / (sigma * sigma) * psi
            q0[2] = (xm - x0) / (sigma * sigma) * psi
        elif kind == "entropy_spot":
            q0[0] = amp * g
        else:
            raise ValueError(kind)
        return q0

    q = initial_state(xx, yy)
    q_ref = initial_state(xx_ref, yy_ref)

    if kind == "oblique_acoustic":
        theta = math.radians(30.0)
        exit_speed = U0 + C0 * math.cos(theta)
        t_end = (lx - x0) / exit_speed + 0.22
        metric_name = "sqrt(boundary-error energy / initial energy)"
        notes = "30 deg acoustic packet; compared with longer-domain reference"
    elif kind == "weak_vortex":
        t_end = (lx - x0) / U0 + 0.20
        metric_name = "sqrt(boundary-error kinetic energy / initial kinetic energy)"
        notes = "linear vortical mode; compared with longer-domain reference"
    elif kind == "entropy_spot":
        t_end = (lx - x0) / U0 + 0.20
        metric_name = "sqrt(boundary-error density L2 / initial density L2)"
        notes = "entropy/contact mode; compared with longer-domain reference"
    else:
        raise ValueError(kind)

    e0 = energy(q, dx, dy)
    rho0_l2 = float(np.sum(q[0] * q[0]) * dx * dy)
    cfl = 0.42
    dt = cfl * min(dx / (abs(U0) + C0), dy / C0)
    t = 0.0
    ts = [0.0]
    vals = [0.0]
    denom = max(e0 if kind != "entropy_spot" else rho0_l2, 1.0e-300)

    while t < t_end - 1.0e-14:
        step = min(dt, t_end - t)
        q = rk3_step(q, step, dx, dy, bc)
        q_ref = rk3_step(q_ref, step, dx, dy, "foextrap")
        t += step
        if len(ts) < 180 and (len(ts) == 1 or t - ts[-1] > t_end / 160.0):
            qerr = q - q_ref[:, :nx, :]
            ts.append(t)
            if kind == "entropy_spot":
                vals.append(float(np.sum(qerr[0] * qerr[0]) * dx * dy) / denom)
            else:
                vals.append(energy(qerr, dx, dy) / denom)

    qerr = q - q_ref[:, :nx, :]
    if kind == "entropy_spot":
        residual = float(np.sum(qerr[0] * qerr[0]) * dx * dy)
        metric = math.sqrt(max(residual, 0.0) / denom)
    else:
        residual = energy(qerr, dx, dy)
        metric = math.sqrt(max(residual, 0.0) / denom)

    data = {
        "time": np.asarray(ts),
        "signal": np.sqrt(np.maximum(np.asarray(vals), 0.0)),
        "final_pressure": qerr[3].copy(),
        "final_density": qerr[0].copy(),
    }
    result = TestResult(kind, bc, metric, metric_name, notes)
    return result, data


def run_1d_suite(quick: bool, bcs: list[str]) -> list[dict[str, object]]:
    machs = [0.0, 0.3, 0.8]
    ncell = 700 if quick else 1100
    rows = []
    for mach in machs:
        for bc in bcs:
            rows.append(bc1d.run_case(mach, bc, ncell=ncell))
    return rows


def plot_1d(rows: list[dict[str, object]], outdir: Path, bcs: list[str]) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for bc in bcs:
        xs = [float(r["mach"]) for r in rows if r["bc"] == bc]
        ys = [abs(float(r["reflection"])) for r in rows if r["bc"] == bc]
        ax.semilogy(xs, ys, marker="o", label=BC_LABELS[bc])
    ax.set_xlabel("background normal Mach number")
    ax.set_ylabel("|R| from 1D acoustic pulse")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = outdir / "case1_1d_reflection.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_2d_summary(results: list[TestResult], outdir: Path, bcs: list[str]) -> Path:
    cases = ["oblique_acoustic", "weak_vortex", "entropy_spot"]
    labels = ["2D oblique acoustic", "weak vortex", "entropy spot"]
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6), sharey=True)
    for ax, case, label in zip(axes, cases, labels):
        vals = [r.metric for r in results if r.case == case]
        ax.bar(range(len(bcs)), vals)
        ax.set_title(label)
        ax.set_xticks(range(len(bcs)))
        ax.set_xticklabels([BC_LABELS[b].replace(" ", "\n") for b in bcs], fontsize=7)
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].set_ylabel("residual amplitude ratio")
    fig.tight_layout()
    path = outdir / "case234_2d_residuals.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_traces(
    trace_data: dict[tuple[str, str], dict[str, np.ndarray]],
    outdir: Path,
    bcs: list[str],
) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4), sharey=True)
    cases = ["oblique_acoustic", "weak_vortex", "entropy_spot"]
    labels = ["oblique acoustic", "weak vortex", "entropy spot"]
    for ax, case, label in zip(axes, cases, labels):
        for bc in bcs:
            d = trace_data[(case, bc)]
            ax.plot(d["time"], d["signal"], label=BC_LABELS[bc], lw=1.4)
        ax.set_title(label)
        ax.set_xlabel("time")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("domain residual amplitude")
    axes[-1].legend(fontsize=7)
    fig.tight_layout()
    path = outdir / "case234_time_traces.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_csv(outdir: Path, rows1d: list[dict[str, object]], rows2d: list[TestResult]) -> None:
    with (outdir / "case1_1d_reflection.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "mach", "bc", "reflection_signed", "reflection_abs"])
        for r in rows1d:
            writer.writerow(["normal_acoustic_1d", r["mach"], r["bc"], r["reflection"], abs(float(r["reflection"]))])

    with (outdir / "case234_2d_metrics.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "bc", "metric", "metric_name", "notes"])
        for r in rows2d:
            writer.writerow([r.case, r.bc, r.metric, r.metric_name, r.notes])


def write_report(
    outdir: Path,
    rows1d: list[dict[str, object]],
    rows2d: list[TestResult],
    fig1: Path,
    fig2: Path,
    fig3: Path,
    nx: int,
    ny: int,
    quick: bool,
    bcs: list[str],
    include_legacy: bool,
) -> tuple[Path, Path]:
    md = outdir / "bc_quality_report.md"
    pdf = outdir / "bc_quality_report.pdf"

    def fmt(x: float) -> str:
        return f"{x:.3e}"

    lines = [
        "# Boundary-condition quality suite",
        "",
        f"Mode: {'quick' if quick else 'standard'}; 2D mesh: {nx} x {ny}; base Mach: {MACH}.",
        f"BC set: {'active + legacy diagnostic' if include_legacy else 'active production BCs only'}.",
        "The 2D metrics compare the test-domain solution against a longer-domain reference solved with the same scheme, so the reported number is a boundary-induced error rather than raw numerical diffusion.",
        "",
        "## Cases",
        "",
        "1. 1D nonlinear Euler/HLLC normal acoustic pulse reflection.",
        "2. 2D linearized-Euler oblique acoustic packet at 30 deg.",
        "3. 2D linearized-Euler weak vortical packet advected through outflow.",
        "4. 2D linearized-Euler entropy/contact spot advected through outflow.",
        "",
        "## Boundary conditions audited",
        "",
        "- foextrap: first-order extrapolated ghost state, U_g = U_i.  Use for supersonic outflow and as a robust lateral/open boundary when vortical or entropy content dominates.",
        "- hard Dirichlet: prescribed freestream ghost state, U_g = U_inf.  Use for true inflow, especially supersonic inflow.  It is not a general non-reflecting outlet.",
        "- fixed characteristic NSCBC: incoming acoustic invariant uses freestream; outgoing acoustic/entropy/tangential data use the interior.  This is a ghost-cell characteristic far-field approximation, not full Poinsot-Lele LODI.",
        "",
        "## 1D acoustic reflection",
        "",
        "| M_n | " + " | ".join(BC_LABELS[bc] for bc in bcs) + " |",
        "|---:|" + "|".join("---:" for _ in bcs) + "|",
    ]
    if include_legacy:
        lines.insert(
            lines.index("- fixed characteristic NSCBC: incoming acoustic invariant uses freestream; outgoing acoustic/entropy/tangential data use the interior.  This is a ghost-cell characteristic far-field approximation, not full Poinsot-Lele LODI."),
            "- old pressure-clamped NSCBC: legacy diagnostic only.  It clamps p_g = p_inf while keeping other outgoing information and is not a production BC.",
        )
    for mach in sorted({float(r["mach"]) for r in rows1d}):
        vals = []
        for bc in bcs:
            r = next(rr for rr in rows1d if float(rr["mach"]) == mach and rr["bc"] == bc)
            vals.append(fmt(abs(float(r["reflection"]))))
        lines.append(f"| {mach:.1f} | " + " | ".join(vals) + " |")

    lines.extend([
        "",
        "## 2D residual metrics",
        "",
        "Lower is better. Values are amplitude-like residual ratios after the outgoing packet has crossed the right boundary.",
        "",
        "| Case | " + " | ".join(BC_LABELS[bc] for bc in bcs) + " |",
        "|---|" + "|".join("---:" for _ in bcs) + "|",
    ])
    for case in ["oblique_acoustic", "weak_vortex", "entropy_spot"]:
        vals = []
        for bc in bcs:
            r = next(rr for rr in rows2d if rr.case == case and rr.bc == bc)
            vals.append(fmt(r.metric))
        lines.append(f"| {case} | " + " | ".join(vals) + " |")

    lines.extend([
        "",
        "## Interpretation",
        "",
        "- The fixed characteristic NSCBC fixes the acoustic inconsistency: incoming acoustic information is far-field, while outgoing acoustic/entropy information is taken from the interior.  It is not universally best for vortical packets with appreciable normal-velocity perturbation.",
        "- foextrap is the cleanest low-risk option for convected vortical/entropy content in these tests and remains appropriate for supersonic outflow.  Its limitation is that it has no explicit incoming acoustic control.",
        "- Hard Dirichlet can look acceptable in the 1D HLLC pulse test, but the 2D residual tests show it is less robust for general outgoing packets because it forces all perturbation components to zero in the ghost state.",
        "",
        "## Practical recommendation",
        "",
        "- Supersonic inflow: use hard freestream Dirichlet.",
        "- Supersonic outflow: foextrap and characteristic far-field are effectively equivalent; use foextrap unless a unified input policy requires code 7.",
        "- Lateral/open boundaries with M_n near zero and wake/vorticity content: prefer foextrap for the present Cerisse ghost-cell implementation.",
        "- Acoustic far-field dominated subsonic outlet: use the fixed characteristic NSCBC, but treat it as an approximation until a true LODI/NSCBC equation-level boundary is implemented.",
        "",
        "## Files",
        "",
        f"- {fig1.name}",
        f"- {fig2.name}",
        f"- {fig3.name}",
        "- case1_1d_reflection.csv",
        "- case234_2d_metrics.csv",
    ])
    if include_legacy:
        marker = lines.index("- The fixed characteristic NSCBC fixes the acoustic inconsistency: incoming acoustic information is far-field, while outgoing acoustic/entropy information is taken from the interior.  It is not universally best for vortical packets with appreciable normal-velocity perturbation.")
        lines.insert(
            marker,
            "- The old pressure-clamped algebraic NSCBC is the diagnostic bad case: it removes pressure perturbations while keeping outgoing velocity/entropy information, which creates an inconsistent ghost state for acoustic waves.",
        )
        lines.insert(
            marker + 3,
            "- The old pressure-clamped NSCBC can look good on pure entropy/vorticity packets because those packets carry little pressure perturbation.  That does not rescue it for acoustic/open-farfield use.",
        )
        lines.insert(
            lines.index("## Files") - 1,
            "- The old pressure-clamped NSCBC branch is not present in production code 7; it is only retained here as a legacy diagnostic.",
        )
    md.write_text("\n".join(lines) + "\n")

    with PdfPages(pdf) as pages:
        for title, image in [
            ("Case 1: 1D acoustic reflection", fig1),
            ("Cases 2-4: residual amplitude ratios", fig2),
            ("Cases 2-4: residual time traces", fig3),
        ]:
            fig = plt.figure(figsize=(8.5, 6.0))
            fig.suptitle(title, fontsize=14)
            img = plt.imread(image)
            ax = fig.add_axes([0.06, 0.08, 0.88, 0.82])
            ax.imshow(img)
            ax.axis("off")
            pages.savefig(fig)
            plt.close(fig)

        fig = plt.figure(figsize=(8.5, 11.0))
        ax = fig.add_axes([0.08, 0.06, 0.84, 0.88])
        ax.axis("off")
        text = "\n".join(lines[:70])
        ax.text(0.0, 1.0, text, va="top", ha="left", fontsize=8, family="monospace")
        pages.savefig(fig)
        plt.close(fig)

    return md, pdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="temp/bc_quality_suite")
    parser.add_argument("--nx", type=int, default=144)
    parser.add_argument("--ny", type=int, default=72)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--include-legacy",
        action="store_true",
        help="include old pressure-clamped NSCBC as a diagnostic comparison",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    nx = 112 if args.quick else args.nx
    ny = 56 if args.quick else args.ny

    bcs = LEGACY_BCS if args.include_legacy else ACTIVE_BCS

    print("Running 1D HLLC acoustic reflection cases...")
    rows1d = run_1d_suite(args.quick, bcs)
    print("Running 2D linearized-Euler boundary cases...")
    rows2d: list[TestResult] = []
    trace_data: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for case in ["oblique_acoustic", "weak_vortex", "entropy_spot"]:
        for bc in bcs:
            result, data = run_2d_case(case, bc, nx, ny, args.quick)
            rows2d.append(result)
            trace_data[(case, bc)] = data
            print(f"  {case:18s} {bc:16s} {result.metric:.3e}")

    write_csv(outdir, rows1d, rows2d)
    fig1 = plot_1d(rows1d, outdir, bcs)
    fig2 = plot_2d_summary(rows2d, outdir, bcs)
    fig3 = plot_traces(trace_data, outdir, bcs)
    md, pdf = write_report(
        outdir, rows1d, rows2d, fig1, fig2, fig3, nx, ny,
        args.quick, bcs, args.include_legacy,
    )
    print(f"Wrote {md}")
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()
