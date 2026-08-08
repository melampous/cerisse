#!/usr/bin/env python3
"""Fourier analysis for the periodic advected shear-wave method matrix."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from pathlib import Path

import numpy as np

try:
    import yt
except ImportError as exc:  # pragma: no cover
    raise SystemExit("analyze.py requires yt") from exc


os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")

GAMMA = 1.4
RHO0 = 1.0
P0 = 1.0
MACH = 0.5
AMPLITUDE_OVER_C0 = 1.0e-5
MODE = 1
X_LO = 0.0
X_HI = 1.0
C2 = 1.5
C4 = 0.016

SCHEMES = (
    "llf-wenoz5",
    "llf-teno5",
    "afd-hllc-wenoz5",
    "skew-jst-o4",
)
PPW = (8, 12, 16, 24, 32, 48)


def latest_plotfile(path: Path) -> Path:
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        candidate
        for candidate in path.glob("plt*")
        if pattern.fullmatch(candidate.name)
        and (candidate / "Header").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"no plotfile found in {path}")
    return max(
        candidates,
        key=lambda candidate: int(pattern.fullmatch(candidate.name)[1]),
    )


def field_key(ds, name: str):
    wanted = name.lower()
    for key in ds.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"{name!r} is absent from {ds.field_list}")


def wrapped_difference(angle: float, reference: float) -> float:
    return math.atan2(
        math.sin(angle - reference),
        math.cos(angle - reference),
    )


def skew_semidiscrete_prediction(ppw: int, time: float) -> tuple[float, float]:
    """Return the analytic C4 amplitude ratio and equivalent viscosity.

    In this exact shear mode, rho and p are uniform.  The JST sensor is zero.
    The order-four damping flux therefore has Fourier decay rate

      sigma = 16*C4*lambda/dx*sin(theta/2)^4.
    """

    length = X_HI - X_LO
    dx = length / ppw
    theta = 2.0 * math.pi * MODE / ppw
    c0 = math.sqrt(GAMMA * P0 / RHO0)
    ux = MACH * c0
    spectral_radius = abs(ux) + c0
    wave_number = 2.0 * math.pi * MODE / length
    sigma = (
        16.0
        * C4
        * spectral_radius
        / dx
        * math.sin(0.5 * theta) ** 4
    )
    amplitude_ratio = math.exp(-sigma * time)
    equivalent_viscosity = sigma / wave_number**2
    return amplitude_ratio, equivalent_viscosity


def analyse_case(
    case_dir: Path,
    scheme: str,
    ppw: int,
    cfl: float,
) -> dict[str, object]:
    plotfile = latest_plotfile(case_dir)
    ds = yt.load(str(plotfile))
    if int(ds.index.max_level) != 0:
        raise ValueError(f"{plotfile} is not a Level-0-only calculation")
    if int(ds.domain_dimensions[0]) != ppw:
        raise ValueError(
            f"{plotfile} has N={int(ds.domain_dimensions[0])}, expected {ppw}"
        )

    data = ds.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    volume = np.asarray(data[("index", "cell_volume")].d, dtype=float)
    rho = np.asarray(data[field_key(ds, "Density")].d, dtype=float)
    ymom = np.asarray(data[field_key(ds, "Ymom")].d, dtype=float)
    uy = ymom / rho
    time = float(ds.current_time)

    length = X_HI - X_LO
    c0 = math.sqrt(GAMMA * P0 / RHO0)
    ux = MACH * c0
    amplitude0 = AMPLITUDE_OVER_C0 * c0
    wave_number = 2.0 * math.pi * MODE / length
    total_volume = float(np.sum(volume))

    # The sine and cosine projections retain the correct cell-centre phase.
    sine_coefficient = float(
        2.0
        * np.sum(uy * np.sin(wave_number * (x - X_LO)) * volume)
        / total_volume
    )
    cosine_coefficient = float(
        2.0
        * np.sum(uy * np.cos(wave_number * (x - X_LO)) * volume)
        / total_volume
    )
    amplitude = math.hypot(sine_coefficient, cosine_coefficient)
    amplitude_ratio = amplitude / amplitude0

    # For uy=A*sin(k*x-phi), S-i*C=A*exp(i*phi).
    measured_phase_modulo = math.atan2(
        -cosine_coefficient,
        sine_coefficient,
    )
    exact_phase = wave_number * ux * time
    phase_error = wrapped_difference(measured_phase_modulo, exact_phase)
    phase_speed_error_over_ux = phase_error / exact_phase

    if amplitude_ratio <= 0.0:
        equivalent_viscosity = float("nan")
    else:
        equivalent_viscosity = (
            -math.log(amplitude_ratio) / (wave_number**2 * time)
        )

    uy_mean = float(np.sum(uy * volume) / total_volume)
    perturbation_ke = float(
        np.sum(0.5 * rho * (uy - uy_mean) ** 2 * volume)
    )
    exact_initial_ke = 0.25 * RHO0 * amplitude0**2 * total_volume
    ke_ratio = perturbation_ke / exact_initial_ke

    skew_amp, skew_nu = skew_semidiscrete_prediction(ppw, time)
    if scheme == "skew-jst-o4":
        skew_relative_difference = amplitude_ratio / skew_amp - 1.0
    else:
        skew_amp = float("nan")
        skew_nu = float("nan")
        skew_relative_difference = float("nan")

    return {
        "scheme": scheme,
        "cfl": cfl,
        "ppw": ppw,
        "plotfile": str(plotfile),
        "time": time,
        "amplitude": amplitude,
        "amplitude_ratio": amplitude_ratio,
        "ke_ratio": ke_ratio,
        "phase_error_rad": phase_error,
        "phase_speed_error_over_ux": phase_speed_error_over_ux,
        "nu_num": equivalent_viscosity,
        "skew_semidiscrete_amplitude_ratio": skew_amp,
        "skew_semidiscrete_nu": skew_nu,
        "skew_amplitude_relative_difference": skew_relative_difference,
        "density_max_error": float(np.max(np.abs(rho - RHO0))),
        "transverse_velocity_mean": uy_mean,
    }


def load_matrix(root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scheme in SCHEMES:
        for ppw in PPW:
            path = root / "runs" / scheme / "cfl0p05" / f"N{ppw}"
            rows.append(analyse_case(path, scheme, ppw, 0.05))
        repeat = root / "runs" / scheme / "cfl0p025" / "N24"
        rows.append(analyse_case(repeat, scheme, 24, 0.025))
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_skew_check(path: Path, rows: list[dict[str, object]]) -> None:
    skew_rows = [
        row
        for row in rows
        if row["scheme"] == "skew-jst-o4" and row["cfl"] == 0.05
    ]
    maximum = max(
        abs(float(row["skew_amplitude_relative_difference"]))
        for row in skew_rows
    )
    with path.open("w") as stream:
        stream.write(
            "Analytic curve: exp[-16 C4 lambda T sin(theta/2)^4 / dx]\n"
        )
        stream.write("C2=1.5\n")
        stream.write("C4=0.016\n")
        stream.write("lambda=abs(ux)+c0\n")
        stream.write(f"max_relative_amplitude_difference={maximum:.17e}\n")
        stream.write(
            "The comparison includes spatial C4 damping, SSPRK33 time error, "
            "roundoff, and any nonzero sensor response.\n"
        )


def plot_results(path: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    matplotlib.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 9.0,
            "axes.labelsize": 9.0,
            "legend.fontsize": 8.0,
        }
    )

    styles = {
        "llf-wenoz5": ("#0072B2", "o", "LLF WENO-Z5"),
        "llf-teno5": ("#009E73", "s", "LLF TENO5"),
        "afd-hllc-wenoz5": ("#D55E00", "^", "AFD HLLC WENO-Z5"),
        "skew-jst-o4": ("#CC79A7", "D", "Skew JST"),
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.25))
    for scheme in SCHEMES:
        colour, marker, label = styles[scheme]
        selected = [
            row
            for row in rows
            if row["scheme"] == scheme and row["cfl"] == 0.05
        ]
        x = np.array([int(row["ppw"]) for row in selected])
        amplitude_deficit = np.abs(
            1.0 - np.array([float(row["amplitude_ratio"]) for row in selected])
        )
        phase_error = np.abs(
            np.array(
                [float(row["phase_speed_error_over_ux"]) for row in selected]
            )
        )
        nu_num = np.abs(
            np.array([float(row["nu_num"]) for row in selected])
        )
        axes[0].loglog(x, amplitude_deficit, color=colour, marker=marker,
                       label=label)
        axes[1].loglog(x, phase_error, color=colour, marker=marker)
        axes[2].loglog(x, nu_num, color=colour, marker=marker)

    skew = [
        row
        for row in rows
        if row["scheme"] == "skew-jst-o4" and row["cfl"] == 0.05
    ]
    axes[0].loglog(
        [int(row["ppw"]) for row in skew],
        [
            abs(1.0 - float(row["skew_semidiscrete_amplitude_ratio"]))
            for row in skew
        ],
        color="black",
        linestyle="--",
        linewidth=1.0,
        label="Skew C4 analysis",
    )

    axes[0].set_ylabel(r"$|1-A(T)/A(0)|$")
    axes[1].set_ylabel(r"$|c_{\mathrm{num}}/u_0-1|$")
    axes[2].set_ylabel(r"$|\nu_{\mathrm{num}}|$")
    for axis in axes:
        axis.set_xlabel("Points per wavelength")
        axis.set_xlim(7.5, 50.5)
        axis.set_xticks((8, 16, 32, 48))
        axis.set_xticklabels(("8", "16", "32", "48"))
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.grid(True, which="both", color="0.88", linewidth=0.6)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()

    rows = load_matrix(args.root)
    csv_path = args.csv or args.root / "metrics.csv"
    figure_path = args.figure or args.root / "dissipation_phase.png"
    write_csv(csv_path, rows)
    write_skew_check(args.root / "analytic_skew_check.txt", rows)
    plot_results(figure_path, rows)

    for row in rows:
        print(
            f"{row['scheme']:22s} CFL={float(row['cfl']):.3f} "
            f"PPW={int(row['ppw']):2d} "
            f"A/A0={float(row['amplitude_ratio']):.10f} "
            f"KE/KE0={float(row['ke_ratio']):.10f} "
            f"dc/u={float(row['phase_speed_error_over_ux']):+.3e} "
            f"nu_num={float(row['nu_num']):+.3e}"
        )


if __name__ == "__main__":
    main()
