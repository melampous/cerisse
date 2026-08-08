#!/usr/bin/env python3
"""Locate and decompose the first inadmissible LLF-WENO RK stage."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_muscl_vs_llfweno_q2 import (
    DT_FINE, ET, GAMMA, MX, RHO, face_flux,
)


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "compare_muscl_vs_llfweno_q2"
AUDIT = BASE / "llfweno" / "stage_audit"
OUT = BASE / "analysis"
TARGET = (56, 133)
DR = DZ = 0.0010612966666666667

BEGIN_RE = re.compile(
    r"\[STAGE-STENCIL-BEGIN\].*level_step=(\d+).*label=([^ ]+) "
    r"time=([^ ]+).*center=\((\d+),(\d+)\)"
)
POINT_RE = re.compile(
    r"\[STAGE-STENCIL\] label=([^ ]+) axis=([^ ]+) offset=(-?\d+) "
    r"cell=\((\d+),(\d+)\) rho=([^ ]+) mx=([^ ]+) my=([^ ]+) "
    r"mz=([^ ]+) E=([^ ]+) rhoe_raw=([^ ]+) p_raw=([^ ]+) "
    r"ur=([^ ]+) uz=([^ ]+) rhoe_used_by_cons2prims=([^ ]+) "
    r"p_used_by_cons2prims=([^ ]+)"
)
STAGE_RE = re.compile(
    r"\[ALL-FLUID-STAGE\] level=(\d+) level_step=(\d+) label=([^ ]+) "
    r"time=([^ ]+) rhoe_min=([^ ]+) cell=\((\d+),(\d+)\) "
    r"r=([^ ]+) z=([^ ]+)"
)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def parse_trace(path: Path) -> tuple[list[dict], dict[str, float]]:
    rows, times = [], {}
    for line in path.read_text().splitlines():
        b = BEGIN_RE.search(line)
        if b:
            times[b.group(2)] = float(b.group(3))
            continue
        m = POINT_RE.search(line)
        if not m:
            continue
        (label, axis, offset, i, j, rho, mx, my, mz, energy, rhoe,
         pressure, ur, uz, used_rhoe, used_p) = m.groups()
        rows.append({
            "label": label, "time_us": 0.0, "axis": axis,
            "offset": int(offset), "i": int(i), "j": int(j),
            "rho": float(rho), "mx": float(mx), "my": float(my),
            "mz": float(mz), "energy_J_m3": float(energy),
            "rhoe_raw_J_m3": float(rhoe), "p_raw_Pa": float(pressure),
            "ur_m_s": float(ur), "uz_m_s": float(uz),
            "rhoe_used_by_cons2prims_J_m3": float(used_rhoe),
            "p_used_by_cons2prims_Pa": float(used_p),
        })
    for row in rows:
        row["time_us"] = times[row["label"]] * 1.0e6
    return rows, times


def state_from_rows(rows: list[dict], label: str) -> np.ndarray:
    U = np.full((64, 192, 5), np.nan)
    for r in rows:
        if r["label"] != label:
            continue
        U[r["i"], r["j"]] = [r["mx"], r["my"], r["mz"],
                               r["energy_J_m3"], r["rho"]]
    return U


def rhoe(U: np.ndarray) -> float:
    return float(U[ET] - 0.5*np.dot(U[:3],U[:3])/U[RHO])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    first_rows, _ = parse_trace(AUDIT / "first_event_trace.log")
    write_csv(OUT / "first_violation_stencil_all_stages.csv", first_rows)

    # Full stage-minimum history exposes violations hidden by the final RK check.
    stage_rows = []
    for line in (AUDIT / "run.log").read_text().splitlines():
        m = STAGE_RE.search(line)
        if not m:
            continue
        level, step, label, time, minimum, i, j, r, z = m.groups()
        if int(level) != 1:
            continue
        stage_rows.append({
            "level_step": int(step), "label": label,
            "time_us": float(time)*1.0e6, "min_rhoe_J_m3": float(minimum),
            "i": int(i), "j": int(j), "r_mm": float(r)*1.0e3,
            "z_mm": float(z)*1.0e3,
        })
    write_csv(OUT / "llfweno_rk_stage_minimum_history.csv", stage_rows)

    # Follow the first bad cell through the four preceding fine steps.
    point_history = []
    trace_paths = [AUDIT / f"trace_{s}_cell56_133.log" for s in (1080,1081,1082)]
    trace_paths.append(AUDIT / "first_event_trace.log")
    for step, path in zip((1080,1081,1082,1083), trace_paths):
        rows, _ = parse_trace(path)
        for row in rows:
            if row["axis"] == "axial" and row["offset"] == 0:
                point_history.append({"level_step": step, **row})
    write_csv(OUT / "first_bad_cell_raw_history.csv", point_history)

    # Reconstruct all four faces from the exact t=216.6 us stage-1 input.
    U = state_from_rows(first_rows, "heun_stage1_input")
    faces = {
        "radial_low": face_flux(U, TARGET[0], 0, DR, center=TARGET),
        "radial_high": face_flux(U, TARGET[0]+1, 0, DR, center=TARGET),
        "axial_low": face_flux(U, TARGET[1], 1, DR, center=TARGET),
        "axial_high": face_flux(U, TARGET[1]+1, 1, DR, center=TARGET),
    }
    i, j = TARGET
    rlo, rhi = i*DR, (i+1)*DR
    volr = (rlo+rhi)*DR
    contributions = {}
    for name, sign, radius in (("radial_low",1.0,rlo),
                               ("radial_high",-1.0,rhi)):
        F, P = faces[name]["flux"], faces[name]["pressure_face"]
        c = sign*2.0*radius*F/volr
        c[MX] = sign*(2.0*radius*(F[MX]-P)/volr + P/DR)
        contributions[name] = c
    contributions["axial_low"] = faces["axial_low"]["flux"]/DZ
    contributions["axial_high"] = -faces["axial_high"]["flux"]/DZ
    rhs = sum(contributions.values())
    U0 = U[TARGET]
    Ufe_field = state_from_rows(first_rows, "heun_fe_stage1_candidate")
    Ufe = Ufe_field[TARGET]
    inferred = (Ufe-U0)/DT_FINE
    component_names = ["mx","my","mz","energy","rho"]
    rhs_rows = []
    for n,name in enumerate(component_names):
        rhs_rows.append({
            "component": name,
            **{f"{face}_rhs": contributions[face][n] for face in contributions},
            "offline_total_rhs": rhs[n], "inferred_FE_rhs": inferred[n],
            "relative_error": abs(rhs[n]-inferred[n])/max(abs(inferred[n]),1.0),
        })
    write_csv(OUT / "first_violation_face_rhs_decomposition.csv", rhs_rows)

    attribution = []
    for name,c in {**contributions, "all_four":rhs,
                   "all_except_axial_high":rhs-contributions["axial_high"]}.items():
        candidate = U0+DT_FINE*c
        attribution.append({
            "case": name, "candidate_rhoe_J_m3": rhoe(candidate),
            "candidate_pressure_Pa": (GAMMA-1.0)*rhoe(candidate),
        })
    write_csv(OUT / "first_violation_internal_energy_attribution.csv", attribution)

    weights = []
    for face_name,face in faces.items():
        for row in face["rows"]:
            weights.append({"face":face_name, **row})
    write_csv(OUT / "first_violation_weno_q2_weights.csv", weights)

    # Compact plot of the hidden FE violation and accepted Heun state.
    fig,ax = plt.subplots(figsize=(9,4.8), constrained_layout=True)
    for label,name,marker in (("heun_stage1_input","stage input","o"),
                              ("heun_fe_stage1_candidate","FE candidate","x"),
                              ("heun_final","Heun final","s")):
        rr=sorted((r for r in point_history if r["label"]==label),
                  key=lambda r:r["time_us"])
        ax.plot([r["time_us"] for r in rr],
                [r["rhoe_raw_J_m3"] for r in rr], marker=marker,label=name)
    ax.axhline(0,color="k",lw=1)
    ax.set_xlabel("time [μs]")
    ax.set_ylabel(r"cell (56,133) raw $\rho e$ [J m$^{-3}$]")
    ax.legend()
    fig.savefig(OUT/"first_bad_cell_rk_stage_history.png",dpi=190)
    plt.close(fig)

    print("first violation:", next(r for r in stage_rows if r["min_rhoe_J_m3"]<0))
    print("RHS relative errors:", [r["relative_error"] for r in rhs_rows])
    print("attribution:", attribution)


if __name__ == "__main__":
    main()
