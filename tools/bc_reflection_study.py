#!/usr/bin/env python3
"""One-dimensional acoustic boundary-reflection study for Euler ghost-cell BCs.

The script runs small-amplitude right-running acoustic wave packets through a
finite-volume Euler solver with HLLC interface fluxes.  It compares several
right-boundary ghost-cell closures and writes figures, CSV data, and a PDF
report to temp/bc_reflection_study by default.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path("temp/bc_reflection_study/.mplconfig")),
)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


GAMMA = 1.4
RHO0 = 1.0
P0 = 1.0
C0 = math.sqrt(GAMMA * P0 / RHO0)
GM1 = GAMMA - 1.0


BC_LABELS = {
    "foextrap": "foextrap",
    "dirichlet": "hard Dirichlet",
    "nscbc_pressure": "old pressure-clamped NSCBC",
    "char_farfield": "fixed characteristic NSCBC",
}
ACTIVE_BCS = ["foextrap", "dirichlet", "char_farfield"]
LEGACY_BCS = ["foextrap", "dirichlet", "nscbc_pressure", "char_farfield"]


def prim_to_cons(rho: np.ndarray, u: np.ndarray, p: np.ndarray) -> np.ndarray:
    e = p / GM1 + 0.5 * rho * u * u
    return np.stack([rho, rho * u, e], axis=-1)


def cons_to_prim(ucons: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rho = np.maximum(ucons[..., 0], 1.0e-14)
    u = ucons[..., 1] / rho
    p = GM1 * (ucons[..., 2] - 0.5 * rho * u * u)
    p = np.maximum(p, 1.0e-14)
    return rho, u, p


def flux(ucons: np.ndarray) -> np.ndarray:
    rho, u, p = cons_to_prim(ucons)
    return np.stack(
        [rho * u, rho * u * u + p, (ucons[..., 2] + p) * u],
        axis=-1,
    )


def hllc_flux(ul: np.ndarray, ur: np.ndarray) -> np.ndarray:
    rhol, vl, pl = cons_to_prim(ul)
    rhor, vr, pr = cons_to_prim(ur)
    cl = np.sqrt(GAMMA * pl / rhol)
    cr = np.sqrt(GAMMA * pr / rhor)

    sl = np.minimum(vl - cl, vr - cr)
    sr = np.maximum(vl + cl, vr + cr)
    denom = rhol * (sl - vl) - rhor * (sr - vr)
    sm = (
        pr
        - pl
        + rhol * vl * (sl - vl)
        - rhor * vr * (sr - vr)
    ) / denom

    fl = flux(ul)
    fr = flux(ur)

    def star_state(u: np.ndarray, rho: np.ndarray, vel: np.ndarray,
                   p: np.ndarray, sk: np.ndarray) -> np.ndarray:
        factor = rho * (sk - vel) / (sk - sm)
        ek = u[..., 2] / rho
        estar = ek + (sm - vel) * (sm + p / (rho * (sk - vel)))
        return np.stack([factor, factor * sm, factor * estar], axis=-1)

    ustar_l = star_state(ul, rhol, vl, pl, sl)
    ustar_r = star_state(ur, rhor, vr, pr, sr)

    out = np.empty_like(fl)
    mask_l = 0.0 <= sl
    mask_lstar = (sl <= 0.0) & (0.0 <= sm)
    mask_rstar = (sm <= 0.0) & (0.0 <= sr)
    mask_r = sr <= 0.0

    out[mask_l] = fl[mask_l]
    out[mask_lstar] = fl[mask_lstar] + sl[mask_lstar, None] * (
        ustar_l[mask_lstar] - ul[mask_lstar]
    )
    out[mask_rstar] = fr[mask_rstar] + sr[mask_rstar, None] * (
        ustar_r[mask_rstar] - ur[mask_rstar]
    )
    out[mask_r] = fr[mask_r]
    return out


def farfield_state(mach: float) -> np.ndarray:
    u0 = mach * C0
    return prim_to_cons(np.array(RHO0), np.array(u0), np.array(P0))


def nscbc_pressure_ghost(uin: np.ndarray, mach: float) -> np.ndarray:
    """The old pressure-clamped bc_nscbc_farfield logic for a right x-boundary."""
    rhoi, ui, pi = cons_to_prim(uin)
    ci = np.sqrt(GAMMA * pi / rhoi)
    rhof, uf, pf = RHO0, mach * C0, P0
    cf = C0

    if ui <= -ci:
        return farfield_state(mach)
    if ui >= ci:
        return uin.copy()

    if ui < 0.0:
        jplus_i = ui + 2.0 * ci / GM1
        jminus_f = uf - 2.0 * cf / GM1
        cb = max(0.25 * GM1 * (jplus_i - jminus_f), 1.0e-14)
        sfar = pf / rhof ** GAMMA
        un = 0.5 * (jplus_i + jminus_f)
        rho = (cb * cb / (GAMMA * sfar)) ** (1.0 / GM1)
        p = rho * cb * cb / GAMMA
    else:
        jplus_i = ui + 2.0 * ci / GM1
        sint = pi / rhoi ** GAMMA
        p = pf
        rho = (p / sint) ** (1.0 / GAMMA)
        cb = math.sqrt(GAMMA * p / rho)
        un = jplus_i - 2.0 * cb / GM1
    return prim_to_cons(np.array(rho), np.array(un), np.array(p))


def characteristic_farfield_ghost(uin: np.ndarray, mach: float) -> np.ndarray:
    """Riemann-invariant far-field closure for a right x-boundary.

    Incoming acoustic data use the far-field J_minus. Outgoing acoustic data,
    entropy, and mean advective data use the interior state for subsonic outflow.
    This is a ghost-cell version of a transparent local characteristic far-field,
    not a fixed-pressure outlet.
    """
    rhoi, ui, pi = cons_to_prim(uin)
    ci = math.sqrt(GAMMA * pi / rhoi)
    uf = mach * C0

    if ui <= -ci:
        return farfield_state(mach)
    if ui >= ci:
        return uin.copy()

    jplus_i = ui + 2.0 * ci / GM1
    jminus_f = uf - 2.0 * C0 / GM1
    cb = max(0.25 * GM1 * (jplus_i - jminus_f), 1.0e-14)
    un = 0.5 * (jplus_i + jminus_f)
    entropy = pi / rhoi ** GAMMA if ui >= 0.0 else P0 / RHO0 ** GAMMA
    rho = (cb * cb / (GAMMA * entropy)) ** (1.0 / GM1)
    p = rho * cb * cb / GAMMA
    return prim_to_cons(np.array(rho), np.array(un), np.array(p))


def right_ghost(uin: np.ndarray, mach: float, bc: str) -> np.ndarray:
    if bc == "foextrap":
        return uin.copy()
    if bc == "dirichlet":
        return farfield_state(mach)
    if bc == "nscbc_pressure":
        return nscbc_pressure_ghost(uin, mach)
    if bc == "char_farfield":
        return characteristic_farfield_ghost(uin, mach)
    raise ValueError(f"unknown BC {bc}")


def run_case(mach: float, bc: str, ncell: int = 1400) -> dict[str, object]:
    length = 1.0
    dx = length / ncell
    x = (np.arange(ncell) + 0.5) * dx
    u0 = mach * C0
    amp = 1.0e-4
    x0 = 0.25
    sigma = 0.035
    shape = np.exp(-0.5 * ((x - x0) / sigma) ** 2)
    du = amp * shape
    dp = RHO0 * C0 * du
    drho = dp / (C0 * C0)

    state = prim_to_cons(RHO0 + drho, u0 + du, P0 + dp)
    base = farfield_state(mach)

    sensor_x = 0.72
    sensor_i = int(sensor_x / dx)
    speed_plus = u0 + C0
    speed_minus = max(C0 - u0, 1.0e-8)
    t_inc = (sensor_x - x0) / speed_plus
    t_hit = (length - x0) / speed_plus
    t_ref = t_hit + (length - sensor_x) / speed_minus
    t_end = min(max(t_ref + 6.0 * sigma / speed_minus, 1.1), 2.6)

    time = []
    wp = []
    wm = []
    rho_trace = []
    cfl = 0.45
    t = 0.0

    while t < t_end:
        rho, vel, p = cons_to_prim(state)
        max_speed = float(np.max(np.abs(vel) + np.sqrt(GAMMA * p / rho)))
        dt = cfl * dx / max_speed
        if t + dt > t_end:
            dt = t_end - t

        ext = np.empty((ncell + 2, 3))
        ext[1:-1] = state
        ext[0] = base
        ext[-1] = right_ghost(state[-1], mach, bc)
        f = hllc_flux(ext[:-1], ext[1:])
        state = state - (dt / dx) * (f[1:] - f[:-1])
        t += dt

        rs, us, ps = cons_to_prim(state[sensor_i])
        time.append(t)
        wp.append((us - u0) + (ps - P0) / (RHO0 * C0))
        wm.append((us - u0) - (ps - P0) / (RHO0 * C0))
        rho_trace.append(rs - RHO0)

    time = np.asarray(time)
    wp = np.asarray(wp)
    wm = np.asarray(wm)
    rho_trace = np.asarray(rho_trace)

    inc_hw = max(5.0 * sigma / speed_plus, 0.03)
    ref_hw = max(5.0 * sigma / speed_minus, 0.03)
    inc_mask = (time > t_inc - inc_hw) & (time < t_inc + inc_hw)
    ref_mask = (time > t_ref - ref_hw) & (time < t_ref + ref_hw)
    inc_amp = float(np.max(np.abs(wp[inc_mask])))
    if np.any(ref_mask):
        idx = np.argmax(np.abs(wm[ref_mask]))
        ref_vals = wm[ref_mask]
        ref_amp_signed = float(ref_vals[idx])
    else:
        ref_amp_signed = float("nan")
    reflection = ref_amp_signed / inc_amp

    return {
        "mach": mach,
        "bc": bc,
        "reflection": reflection,
        "incident_amplitude": inc_amp,
        "reflected_amplitude": ref_amp_signed,
        "time": time,
        "wplus": wp,
        "wminus": wm,
        "rho_trace": rho_trace,
        "t_inc": t_inc,
        "t_ref": t_ref,
    }


def write_csv(rows: list[dict[str, object]], outdir: Path) -> None:
    with (outdir / "reflection_coefficients.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["M_bg", "bc", "reflection", "incident_amp", "reflected_amp"])
        for row in rows:
            writer.writerow(
                [
                    row["mach"],
                    row["bc"],
                    f"{row['reflection']:.8e}",
                    f"{row['incident_amplitude']:.8e}",
                    f"{row['reflected_amplitude']:.8e}",
                ]
            )


def plot_coefficients(rows: list[dict[str, object]], outdir: Path,
                      bcs: list[str]) -> None:
    machs = sorted({float(r["mach"]) for r in rows})
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for bc in bcs:
        vals = [
            abs(float(next(r for r in rows if r["mach"] == m and r["bc"] == bc)["reflection"]))
            for m in machs
        ]
        ax.plot(machs, vals, marker="o", label=BC_LABELS[bc])
    ax.set_xlabel(r"background normal Mach number $M_n$")
    ax.set_ylabel(r"$|R| = |\hat w_-|_{\rm refl}/|\hat w_+|_{\rm inc}$")
    ax.set_title("1D HLLC acoustic pulse reflection at right boundary")
    ax.set_ylim(bottom=-0.02)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "reflection_coefficients.png", dpi=220)
    plt.close(fig)


def plot_sensor(rows: list[dict[str, object]], outdir: Path, mach: float,
                bcs: list[str]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.4), sharex=True)
    for bc in bcs:
        row = next(r for r in rows if r["mach"] == mach and r["bc"] == bc)
        t = row["time"]
        axes[0].plot(t, row["wplus"], label=BC_LABELS[bc], linewidth=1.1)
        axes[1].plot(t, row["wminus"], label=BC_LABELS[bc], linewidth=1.1)
    axes[0].set_ylabel(r"$w_+$ at sensor")
    axes[1].set_ylabel(r"$w_-$ at sensor")
    axes[1].set_xlabel("time")
    axes[0].set_title(f"Sensor characteristic traces, M_n={mach:g}")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / f"sensor_traces_M{str(mach).replace('.', 'p')}.png", dpi=220)
    plt.close(fig)


def make_report(rows: list[dict[str, object]], outdir: Path,
                bcs: list[str], include_legacy: bool) -> Path:
    table_lines = []
    machs = sorted({float(r["mach"]) for r in rows})
    for m in machs:
        vals = []
        for bc in bcs:
            row = next(r for r in rows if r["mach"] == m and r["bc"] == bc)
            vals.append(f"{float(row['reflection']):+.3f}")
        table_lines.append(f"{m:.2f} & " + " & ".join(vals) + r" \\")

    table_cols = " ".join(["c"] * (len(bcs) + 1))
    table_header = r"$M_n$ & " + " & ".join(BC_LABELS[bc] for bc in bcs) + r" \\"
    purpose = (
        r"""This report audits the active production ghost-cell open-boundary
closures, including the fixed characteristic-invariant far-field closure now used
by \texttt{bc\_nscbc\_farfield}.  The test is limited to subsonic normal Mach
numbers; supersonic normal outflow is characteristically determined by the
interior state."""
        if not include_legacy
        else
        r"""This report audits the pressure-clamped algebraic ghost-cell NSCBC
previously used in \texttt{bc\_nscbc\_farfield} and compares it with the fixed
characteristic-invariant closure.  The test is limited to subsonic normal Mach
numbers; supersonic normal outflow is characteristically determined by the
interior state."""
    )
    legacy_bc_section = (
        r"""
\paragraph{Old pressure-clamped algebraic NSCBC.}
The previous helper used Riemann invariants for some branches but, for subsonic
outflow, imposed
\[
  p_g=p_\infty,\qquad
  s_g=s_i,\qquad
  J_{+,g}=J_{+,i},
\]
where
\[
  J_+ = u + \frac{2c}{\gamma-1}.
\]
This is a fixed-pressure outlet translated into a ghost cell, not a transparent
characteristic far-field condition.
"""
        if include_legacy
        else ""
    )
    legacy_usage_row = (
        r"""\texttt{old form} & pressure clamp; audit only &
This was the faulty subsonic-outflow algebra \(p_g=p_\infty\),
\(J_{+,g}=J_{+,i}\), \(s_g=s_i\).  It should not be used as a non-reflecting
boundary. \\
"""
        if include_legacy
        else ""
    )
    diagnosis = (
        r"""The poor behaviour of the old pressure-clamped NSCBC is not evidence that
characteristic boundary conditions are inherently bad.  That implementation mixes
a pressure outlet constraint with a ghost-cell Riemann solve.  For a small
outgoing acoustic wave at subsonic outflow, \(p_g=p_\infty\) suppresses the
pressure perturbation in the ghost cell, while \(J_{+,g}=J_{+,i}\) retains the
outgoing velocity/acoustic invariant.  The resulting left and right states at
the HLLC interface are not connected by a single outgoing acoustic wave.  The
Riemann solver therefore decomposes the mismatch into both outgoing and incoming
waves, producing a measurable reflected \(w_-\) pulse.

The fixed characteristic NSCBC test removes that mismatch by prescribing
the incoming acoustic invariant, not the pressure itself.  This recovers near
zero reflection for the same one-dimensional acoustic benchmark."""
        if include_legacy
        else
        r"""The fixed characteristic NSCBC prescribes only the incoming acoustic
invariant from the far field.  It does not clamp the ghost-cell pressure.  This
keeps the interface states much closer to a single outgoing acoustic wave in this
one-dimensional benchmark.  The removed pressure-clamped form can still be
enabled from this script with \texttt{--include-legacy} for historical
comparison, but it is not part of production code 7."""
    )

    tex = r"""
\documentclass[10pt]{{article}}
\usepackage[margin=0.75in]{{geometry}}
\usepackage{{amsmath,amssymb,booktabs,graphicx,siunitx,array}}
\usepackage{{hyperref}}
\title{{One-Dimensional Boundary-Condition Reflection Study}}
\author{{Cerisse boundary-condition audit}}
\date{{Generated by tools/bc\_reflection\_study.py}}
\begin{{document}}
\maketitle

\section{{Purpose}}
@@PURPOSE@@

\section{{Numerical Method}}
The one-dimensional Euler equations,
\[
  \partial_t U + \partial_x F(U)=0,\qquad
  U=(\rho,\rho u,\rho E)^T,
\]
are advanced with a first-order finite-volume scheme and an HLLC numerical
flux.  The initial condition is a small-amplitude right-running acoustic
Gaussian wave packet on a uniform background state:
\[
  \delta u = \epsilon \exp\!\left[-\frac{(x-x_0)^2}{2\sigma^2}\right],\quad
  \delta p = \rho_0 c_0 \delta u,\quad
  \delta \rho = \frac{\delta p}{c_0^2}.
\]
The reflection coefficient is measured at a fixed interior sensor using the
linear acoustic characteristic variables
\[
  w_+ = \delta u + \frac{\delta p}{\rho_0 c_0},\qquad
  w_- = \delta u - \frac{\delta p}{\rho_0 c_0}.
\]
The reported quantity is
\[
  R = \frac{\max_{\rm reflected} w_-}{\max_{\rm incident} |w_+|}.
\]

\section{{Boundary Conditions Tested}}
\paragraph{{First-order extrapolation.}}
The ghost state is copied from the nearest interior cell:
\[
  U_g = U_i .
\]
In Cerisse input files this is \texttt{{cns.*\_bc = 2}}.  It does not enforce a
far-field pressure, but it is transparent to a purely outgoing acoustic wave in
this one-dimensional test.

\paragraph{{Hard Dirichlet far field.}}
The ghost state is fixed to the far-field state:
\[
  U_g = U_\infty .
\]
In a finite-volume HLLC discretisation this is not automatically a rigid wall:
the interface Riemann problem can transmit a characteristic jump if the interior
disturbance lies on an outgoing acoustic integral curve.  It is still a strong
constraint for entropy, vortical, and nonlinear disturbances.

@@LEGACY_BC_SECTION@@

\paragraph{{Fixed characteristic NSCBC.}}
For a right boundary with subsonic outflow, the incoming acoustic invariant is
taken from the far field and the outgoing invariant from the interior:
\[
  J_{+,g}=J_{+,i},\qquad
  J_{-,g}=J_{-,\infty},\qquad
  s_g=s_i,
\]
with
\[
  J_- = u - \frac{2c}{\gamma-1}.
\]
This does not instantaneously clamp pressure to \(p_\infty\); pressure relaxes
only through whatever additional sponge or LODI relaxation term is supplied.

\section{{Cerisse Usage}}
\begin{{table}}[h]
\centering
\small
\begin{{tabular}}{{>{{\raggedright\arraybackslash}}p{{0.12\linewidth}}
                  >{{\raggedright\arraybackslash}}p{{0.20\linewidth}}
                  >{{\raggedright\arraybackslash}}p{{0.56\linewidth}}}}
\toprule
Input value & Closure & Use in a case file \\
\midrule
\texttt{{1}} & hard Dirichlet / \texttt{{ext\_dir}} &
Set \texttt{{cns.lo\_bc}} or \texttt{{cns.hi\_bc}} to \texttt{{1}} and fill
the external conservative state in \texttt{{bcnormal}}.  For a far-field test
this means assigning \(U_g=U_\infty\). \\
\texttt{{2}} & first-order extrapolation &
Set the face to \texttt{{2}}.  AMReX copies the nearest valid-cell state into
the ghost cell, \(U_g=U_i\).  This is simple and robust when no far-field
pressure control is needed. \\
\texttt{{7}} & fixed char. far field &
Set the face to \texttt{{7}} and call the far-field helper from
\texttt{{bcnormal}} with an outward normal and freestream state.  Subsonic
outflow uses \(J_{+,i}\), \(J_{-,\infty}\), and \(s_i\). \\
@@LEGACY_USAGE_ROW@@
\bottomrule
\end{{tabular}}
\caption{{How the tested closures map to Cerisse input files and \texttt{{prob.h}}.}}
\end{{table}}

\section{{Results}}
\begin{{table}}[h]
\centering
\begin{{tabular}}{{@@TABLE_COLS@@}}
\toprule
@@TABLE_HEADER@@
\midrule
@@REFLECTION_TABLE_ROWS@@
\bottomrule
\end{{tabular}}
\caption{{Signed reflection coefficients from the HLLC wave-packet test.}}
\end{{table}}

\begin{{figure}}[h]
\centering
\includegraphics[width=0.82\linewidth]{{reflection_coefficients.png}}
\caption{{Absolute reflection coefficient versus background normal Mach number.}}
\end{{figure}}

\begin{{figure}}[h]
\centering
\includegraphics[width=0.82\linewidth]{{sensor_traces_M0p0.png}}
\caption{{Sensor characteristic traces for \(M_n=0\).}}
\end{{figure}}

\begin{{figure}}[h]
\centering
\includegraphics[width=0.82\linewidth]{{sensor_traces_M0p3.png}}
\caption{{Sensor characteristic traces for \(M_n=0.3\).}}
\end{{figure}}

\section{{Diagnosis}}
@@DIAGNOSIS@@

\section{{Practical Recommendation}}
For the present Cerisse ghost-cell framework:
\begin{{itemize}}
  \item Use \texttt{{2}} for simple supersonic outflow and lateral boundaries
        when no far-field pressure control is needed.
  \item Do not use the old pressure-clamped algebraic NSCBC as a lateral
        low-\(M_n\) non-reflecting boundary.
  \item Use the fixed incoming-invariant branch in
        \texttt{{bc\_nscbc\_farfield}} for characteristic far-field behaviour,
        and add a sponge layer if pressure relaxation to \(p_\infty\) is
        required over a finite distance.
\end{{itemize}}

\end{{document}}
"""
    tex = (
        tex.replace("@@REFLECTION_TABLE_ROWS@@", "\n".join(table_lines))
        .replace("@@PURPOSE@@", purpose)
        .replace("@@LEGACY_BC_SECTION@@", legacy_bc_section)
        .replace("@@LEGACY_USAGE_ROW@@", legacy_usage_row)
        .replace("@@TABLE_COLS@@", table_cols)
        .replace("@@TABLE_HEADER@@", table_header)
        .replace("@@DIAGNOSIS@@", diagnosis)
        .replace("{{", "{")
        .replace("}}", "}")
    )
    tex_path = outdir / "bc_reflection_report.tex"
    tex_path.write_text(tex)
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", tex_path.name],
        cwd=outdir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", tex_path.name],
        cwd=outdir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return outdir / "bc_reflection_report.pdf"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir",
        default="temp/bc_reflection_study",
        help="output directory",
    )
    parser.add_argument("--ncell", type=int, default=1400)
    parser.add_argument(
        "--include-legacy",
        action="store_true",
        help="include old pressure-clamped NSCBC as a diagnostic comparison",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(outdir / ".mplconfig"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    machs = [0.0, 0.1, 0.3, 0.5, 0.8]
    bcs = LEGACY_BCS if args.include_legacy else ACTIVE_BCS
    rows: list[dict[str, object]] = []
    for mach in machs:
        for bc in bcs:
            print(f"running M={mach:g}, bc={bc}")
            rows.append(run_case(mach, bc, ncell=args.ncell))

    write_csv(rows, outdir)
    plot_coefficients(rows, outdir, bcs)
    plot_sensor(rows, outdir, 0.0, bcs)
    plot_sensor(rows, outdir, 0.3, bcs)
    pdf = make_report(rows, outdir, bcs, args.include_legacy)
    print(f"wrote {pdf}")


if __name__ == "__main__":
    main()
