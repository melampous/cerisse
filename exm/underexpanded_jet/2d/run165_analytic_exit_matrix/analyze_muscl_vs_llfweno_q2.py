#!/usr/bin/env python3
"""Matched-time HLLC-MUSCL / q=2 LLF-WENO Run165 RZ diagnosis."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import yt

from analyze_hllc_muscl_afd_weno_pointwise import (
    load_fields,
    patch_yt_cylindrical_readonly_edge,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_muscl_vs_llfweno_q2"
OUT = RESULTS / "analysis"
GAMMA = 1.4
P_INF = 534.5806702256042
RHOE_INF = P_INF / (GAMMA - 1.0)
TARGET = (57, 131)
DT_FINE = 2.0e-7
RZ_WENO_EPS_REL = 1.0e-6

# Conservative component ordering in src/set/Index.h.
MX, MY, MZ, ET, RHO = range(5)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def path_for(scheme: str, step: int) -> Path:
    if step <= 550:
        return RESULTS / scheme / "evolution" / f"plt{step:05d}"
    return RESULTS / scheme / "replay" / f"plt{step:05d}"


def cons_array(fields: dict[str, np.ndarray]) -> np.ndarray:
    return np.stack(
        [fields["mx"], fields["my"], fields["mz"],
         fields["energy"], fields["rho"]], axis=-1
    )


def primitive(U: np.ndarray) -> dict[str, np.ndarray]:
    rho = U[..., RHO]
    vel = U[..., :3] / rho[..., None]
    kinetic = 0.5 * np.sum(U[..., :3] ** 2, axis=-1) / rho
    rhoe = U[..., ET] - kinetic
    pressure = (GAMMA - 1.0) * rhoe
    return {
        "rho": rho,
        "vel": vel,
        "rhoe": rhoe,
        "eint": rhoe / rho,
        "p": pressure,
        "c": np.sqrt(GAMMA * pressure / rho),
    }


def roe_state(Ul: np.ndarray, Ur: np.ndarray, direction: int) -> dict[str, float | int]:
    ql, qr = primitive(Ul), primitive(Ur)
    sl, sr = np.sqrt(ql["rho"]), np.sqrt(qr["rho"])
    wl = sl / (sl + sr)
    cn, ct, ctt = direction, (direction + 1) % 3, (direction + 2) % 3
    uvw = ql["vel"] * wl + qr["vel"] * (1.0 - wl)
    u, v, w = uvw[cn], uvw[ct], uvw[ctt]
    q2 = float(np.dot(uvw, uvw))
    El = ql["eint"] + 0.5 * float(np.dot(ql["vel"], ql["vel"]))
    Er = qr["eint"] + 0.5 * float(np.dot(qr["vel"], qr["vel"]))
    H = (El + ql["p"] / ql["rho"]) * wl + (Er + qr["p"] / qr["rho"]) * (1.0 - wl)
    h = H - 0.5 * q2
    c = np.sqrt((GAMMA - 1.0) * h)
    return {"cn": cn, "ct": ct, "ctt": ctt, "u": u, "v": v,
            "w": w, "q2": q2, "H": H, "h": h, "c": c}


def cons2char(r: dict, values: np.ndarray) -> np.ndarray:
    f = np.asarray(values, dtype=float)
    out = np.empty(5)
    un, ut, utt = r["cn"], r["ct"], r["ctt"]
    u, v, w, h, c, q2 = r["u"], r["v"], r["w"], r["h"], r["c"], r["q2"]
    out[0] = f[ut] - v * f[RHO]
    out[1] = f[utt] - w * f[RHO]
    out[2] = (
        0.5 * ((-h - c*u) / c * f[un] - v*f[ut] - w*f[utt] + f[ET])
        + (2*h*u + c*q2) / (4*c) * f[RHO]
    ) / h
    out[3] = (
        0.5 * ((h - c*u) / c * f[un] - v*f[ut] - w*f[utt] + f[ET])
        + (-2*h*u + c*q2) / (4*c) * f[RHO]
    ) / h
    out[4] = (u*f[un] + v*f[ut] + w*f[utt] - f[ET] + (h-0.5*q2)*f[RHO]) / h
    return out


def char2cons(r: dict, values: np.ndarray) -> np.ndarray:
    f = np.asarray(values, dtype=float)
    out = np.empty(5)
    un, ut, utt = r["cn"], r["ct"], r["ctt"]
    u, v, w, H, c, q2 = r["u"], r["v"], r["w"], r["H"], r["c"], r["q2"]
    out[un] = f[2]*(u-c) + f[3]*(u+c) + u*f[4]
    out[ut] = f[0] + v*(f[2]+f[3]+f[4])
    out[utt] = f[1] + w*(f[2]+f[3]+f[4])
    out[ET] = v*f[0] + w*f[1] + (H-u*c)*f[2] + (H+u*c)*f[3] + 0.5*q2*f[4]
    out[RHO] = f[2] + f[3] + f[4]
    return out


def physical_flux(U: np.ndarray, direction: int) -> np.ndarray:
    q = primitive(U)
    un = q["vel"][direction]
    flux = np.empty(5)
    flux[:3] = U[:3] * un
    flux[direction] += q["p"]
    flux[ET] = (U[ET] + q["p"]) * un
    flux[RHO] = U[RHO] * un
    return flux


def weno_details(stencil: np.ndarray, epsilon: float) -> dict[str, np.ndarray | float]:
    s = np.asarray(stencil, dtype=float)
    beta = np.empty(3)
    beta[2] = 13/12*(s[4]-2*s[3]+s[2])**2 + 0.25*(s[4]-4*s[3]+3*s[2])**2
    beta[1] = 13/12*(s[3]-2*s[2]+s[1])**2 + 0.25*(s[3]-s[1])**2
    beta[0] = 13/12*(s[2]-2*s[1]+s[0])**2 + 0.25*(3*s[2]-4*s[1]+s[0])**2
    cand6 = np.array([
        2*s[2] + 5*s[1] - s[0],
        -s[3] + 5*s[2] + 2*s[1],
        11*s[2] - 7*s[3] + 2*s[4],
    ])
    tau = abs(beta[2] - beta[0])
    alpha = np.array([3.0, 6.0, 1.0]) * (1.0 + (tau / (epsilon + beta)) ** 2)
    weight = alpha / alpha.sum()
    recon = float(np.dot(weight, cand6) / 6.0)
    return {"samples": s, "beta": beta, "cand": cand6/6.0,
            "alpha": alpha, "weight": weight, "recon": recon}


def face_flux(U: np.ndarray, face: int, direction: int,
              dr: float, r0: float = 0.0,
              center: tuple[int, int] = TARGET) -> dict:
    """Port src/rhs/Weno.h for one all-fluid q=2 WENO-Z5 face."""
    if direction == 0:
        samples = np.stack([U[face + off, center[1]] for off in range(-3, 3)])
    else:
        samples = np.stack([U[center[0], face + off] for off in range(-3, 3)])
    q = primitive(samples)
    alpha_llf = float(np.max(np.abs(q["vel"][:, direction]) + q["c"]))
    roe = roe_state(samples[2], samples[3], direction)
    radial = direction == 0
    r_face = r0 + face * dr if radial else 1.0
    scales = np.ones(6)
    if radial:
        indices = np.arange(face-3, face+3)
        r_cell = r0 + (indices + 0.5) * dr
        scales = r_cell / r_face if r_face > 0.0 else r_cell

    pos, neg = np.empty((6, 5)), np.empty((6, 5))
    pressure_char = np.zeros((6, 5))
    for m in range(6):
        F = physical_flux(samples[m], direction)
        pos[m] = cons2char(roe, 0.5 * scales[m] * (samples[m] + F / alpha_llf))
        neg[m] = cons2char(roe, 0.5 * scales[m] * (samples[m] - F / alpha_llf))
        if radial:
            p_only = np.zeros(5)
            p_only[MX] = 0.5 * q["p"][m] / alpha_llf
            pressure_char[m] = cons2char(roe, p_only)

    epsilon = RZ_WENO_EPS_REL * max(np.max(np.abs(pos)), np.max(np.abs(neg)))**2
    full_char = np.empty(5)
    adv_char = np.empty(5)
    pressure_flux_char = np.empty(5)
    rows: list[dict] = []
    for cidx in range(5):
        p_stencil = pos[4::-1, cidx]
        n_stencil = neg[1:6, cidx]
        pd = weno_details(p_stencil, epsilon)
        nd = weno_details(n_stencil, epsilon)
        full_char[cidx] = alpha_llf * (pd["recon"] - nd["recon"])
        for branch, details in (("positive", pd), ("negative", nd)):
            for cand in range(3):
                rows.append({
                    "direction": "radial" if radial else "axial",
                    "face_index": face,
                    "char_component": cidx,
                    "branch": branch,
                    "candidate": cand,
                    "sample0": details["samples"][0],
                    "sample1": details["samples"][1],
                    "sample2": details["samples"][2],
                    "sample3": details["samples"][3],
                    "sample4": details["samples"][4],
                    "beta": details["beta"][cand],
                    "candidate_value": details["cand"][cand],
                    "alpha_weight": details["alpha"][cand],
                    "normalized_weight": details["weight"][cand],
                    "branch_reconstruction": details["recon"],
                    "llf_alpha": alpha_llf,
                    "face_epsilon": epsilon,
                })

        if radial:
            pos_adv = pos[:, cidx] - scales * pressure_char[:, cidx]
            neg_adv = neg[:, cidx] + scales * pressure_char[:, cidx]
            adv_p = np.dot(pd["weight"], np.array([
                (2*pos_adv[2]+5*pos_adv[3]-pos_adv[4])/6,
                (-pos_adv[1]+5*pos_adv[2]+2*pos_adv[3])/6,
                (11*pos_adv[2]-7*pos_adv[1]+2*pos_adv[0])/6,
            ]))
            adv_n = np.dot(nd["weight"], np.array([
                (2*neg_adv[3]+5*neg_adv[2]-neg_adv[1])/6,
                (-neg_adv[4]+5*neg_adv[3]+2*neg_adv[2])/6,
                (11*neg_adv[3]-7*neg_adv[4]+2*neg_adv[5])/6,
            ]))
            # The arrays above are exactly gather-positive [4,3,2,1,0]
            # and gather-negative [1,2,3,4,5] with frozen weights.
            adv_char[cidx] = alpha_llf * (adv_p - adv_n)
            ppos = pressure_char[4::-1, cidx]
            pneg = -pressure_char[1:6, cidx]
            ppos_cand = np.array([
                (2*ppos[2]+5*ppos[1]-ppos[0])/6,
                (-ppos[3]+5*ppos[2]+2*ppos[1])/6,
                (11*ppos[2]-7*ppos[3]+2*ppos[4])/6,
            ])
            pneg_cand = np.array([
                (2*pneg[2]+5*pneg[1]-pneg[0])/6,
                (-pneg[3]+5*pneg[2]+2*pneg[1])/6,
                (11*pneg[2]-7*pneg[3]+2*pneg[4])/6,
            ])
            pressure_flux_char[cidx] = alpha_llf * (
                np.dot(pd["weight"], ppos_cand) - np.dot(nd["weight"], pneg_cand)
            )

    flux = char2cons(roe, full_char)
    pressure_face = 0.0
    if radial:
        adv_flux = char2cons(roe, adv_char)
        pressure_flux = char2cons(roe, pressure_flux_char)
        pressure_face = float(pressure_flux[MX])
        flux[MX] = adv_flux[MX] + pressure_face if r_face > 0.0 else 0.0
    return {"flux": flux, "pressure_face": pressure_face, "rows": rows,
            "alpha": alpha_llf, "epsilon": epsilon, "samples": samples}


def failed_fe_state_from_log(path: Path) -> np.ndarray:
    text = path.read_text()
    match = re.search(
        r"offending-cell state: rho = ([^ ]+)\s+mom = \(([^,]+), ([^,]+), ([^)]+)\)\s+E = ([^\n]+)",
        text,
    )
    if not match:
        raise RuntimeError(f"cannot parse failed state from {path}")
    rho, mx, my, mz, energy = map(float, match.groups())
    return np.array([mx, my, mz, energy, rho])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    patch_yt_cylindrical_readonly_edge()
    yt.funcs.mylog.setLevel(50)

    broad_steps = list(range(0, 551, 25))
    cached: dict[tuple[str, int], tuple] = {}
    metric_rows: list[dict] = []
    target_rows: list[dict] = []
    for step in broad_steps:
        states = {}
        for scheme in ("llfweno", "muscl"):
            data = load_fields(path_for(scheme, step))
            cached[(scheme, step)] = data
            fields, time_us, _, _, rc, zc = data
            states[scheme] = fields
            target_rows.append({
                "step": step, "time_us": time_us, "scheme": scheme,
                "rho": fields["rho"][TARGET], "pressure_Pa": fields["pressure"][TARGET],
                "rhoe_J_m3": fields["rhoe"][TARGET], "ur_m_s": fields["ur"][TARGET],
                "uz_m_s": fields["uz"][TARGET],
            })
        dp = states["llfweno"]["pressure"] - states["muscl"]["pressure"]
        de = states["llfweno"]["rhoe"] - states["muscl"]["rhoe"]
        arg = np.unravel_index(np.argmax(np.abs(dp)), dp.shape)
        _, time_us, _, _, rc, zc = cached[("llfweno", step)]
        metric_rows.append({
            "step": step, "time_us": time_us,
            "pressure_l1_over_pinf": np.mean(np.abs(dp))/P_INF,
            "pressure_linf_over_pinf": np.max(np.abs(dp))/P_INF,
            "rhoe_l1_over_freestream": np.mean(np.abs(de))/RHOE_INF,
            "rhoe_linf_over_freestream": np.max(np.abs(de))/RHOE_INF,
            "dp_argmax_i": arg[0], "dp_argmax_j": arg[1],
            "dp_argmax_r_mm": rc[arg[0]], "dp_argmax_z_mm": zc[arg[1]],
            "llf_min_pressure_Pa": np.min(states["llfweno"]["pressure"]),
            "muscl_min_pressure_Pa": np.min(states["muscl"]["pressure"]),
            "llf_min_rhoe_J_m3": np.min(states["llfweno"]["rhoe"]),
            "muscl_min_rhoe_J_m3": np.min(states["muscl"]["rhoe"]),
        })

    fine_steps = list(range(551, 558))
    for step in fine_steps:
        for scheme in ("llfweno", "muscl"):
            data = load_fields(path_for(scheme, step))
            cached[(scheme, step)] = data
            f, time_us, _, _, _, _ = data
            target_rows.append({
                "step": step, "time_us": time_us, "scheme": scheme,
                "rho": f["rho"][TARGET], "pressure_Pa": f["pressure"][TARGET],
                "rhoe_J_m3": f["rhoe"][TARGET], "ur_m_s": f["ur"][TARGET],
                "uz_m_s": f["uz"][TARGET],
            })

    write_csv(OUT / "matched_time_metrics.csv", metric_rows)
    write_csv(OUT / "target_cell_history.csv", target_rows)

    stencil_steps = [400, 425, 450, 475, 500, 525, 550] + fine_steps
    stencil_rows: list[dict] = []
    stencil_points = [("radial", off, TARGET[0]+off, TARGET[1]) for off in range(-3,4)]
    stencil_points += [("axial", off, TARGET[0], TARGET[1]+off) for off in range(-3,4)]
    for step in stencil_steps:
        for scheme in ("llfweno", "muscl"):
            if (scheme, step) not in cached:
                cached[(scheme, step)] = load_fields(path_for(scheme, step))
            f, time_us, _, _, rc, zc = cached[(scheme, step)]
            for axis, off, i, j in stencil_points:
                stencil_rows.append({
                    "step": step, "time_us": time_us, "scheme": scheme,
                    "axis": axis, "offset": off, "i": i, "j": j,
                    "r_mm": rc[i], "z_mm": zc[j], "rho": f["rho"][i,j],
                    "pressure_Pa": f["pressure"][i,j], "ur_m_s": f["ur"][i,j],
                    "uz_m_s": f["uz"][i,j], "energy_J_m3": f["energy"][i,j],
                    "rhoe_J_m3": f["rhoe"][i,j],
                })
    write_csv(OUT / "failure_stencil_evolution.csv", stencil_rows)

    # Exact q=2 LLF-WENO face reconstruction at the last admissible state.
    f557, _, r_edges, z_edges, rc, zc = cached[("llfweno", 557)]
    U = cons_array(f557)
    dr = (r_edges[1] - r_edges[0]) * 1.0e-3
    dz = (z_edges[1] - z_edges[0]) * 1.0e-3
    faces = {
        "radial_low": face_flux(U, TARGET[0], 0, dr),
        "radial_high": face_flux(U, TARGET[0]+1, 0, dr),
        "axial_low": face_flux(U, TARGET[1], 1, dr),
        "axial_high": face_flux(U, TARGET[1]+1, 1, dr),
    }
    weight_rows = []
    for label, face in faces.items():
        for row in face["rows"]:
            row = dict(row)
            row["face_label"] = label
            weight_rows.append(row)
    write_csv(OUT / "failure_face_weno_q2.csv", weight_rows)

    i, j = TARGET
    rlo, rhi = i*dr, (i+1)*dr
    volr = (rhi+rlo)*dr
    contributions = {}
    for name, sign, radius in (("radial_low", 1.0, rlo),
                               ("radial_high", -1.0, rhi)):
        F = faces[name]["flux"]
        P = faces[name]["pressure_face"]
        c = sign * 2.0 * radius * F / volr
        c[MX] = sign * (2.0*radius*(F[MX]-P)/volr + P/dr)
        contributions[name] = c
    contributions["axial_low"] = faces["axial_low"]["flux"] / dz
    contributions["axial_high"] = -faces["axial_high"]["flux"] / dz
    rhs = sum(contributions.values())
    U0 = U[TARGET]
    Ufe_reported = failed_fe_state_from_log(
        RESULTS / "llfweno" / "stage_probe" / "rk1.log"
    )
    rhs_inferred = (Ufe_reported - U0) / DT_FINE
    flux_rows = []
    component_names = ["mx", "my", "mz", "energy", "rho"]
    for n, comp in enumerate(component_names):
        flux_rows.append({
            "component": comp,
            "radial_low_rhs": contributions["radial_low"][n],
            "radial_high_rhs": contributions["radial_high"][n],
            "axial_low_rhs": contributions["axial_low"][n],
            "axial_high_rhs": contributions["axial_high"][n],
            "offline_total_rhs": rhs[n], "inferred_FE_rhs": rhs_inferred[n],
            "relative_error": abs(rhs[n]-rhs_inferred[n]) / max(abs(rhs_inferred[n]), 1.0),
            "U_at_222p8us": U0[n], "FE_candidate_at_223us": Ufe_reported[n],
        })
    write_csv(OUT / "failure_face_rhs_decomposition.csv", flux_rows)

    def point_rhoe(state: np.ndarray) -> float:
        return float(state[ET] - 0.5*np.dot(state[:3], state[:3])/state[RHO])

    attribution_rows = []
    cases = {name: contribution for name, contribution in contributions.items()}
    cases["radial_pair"] = contributions["radial_low"] + contributions["radial_high"]
    cases["axial_pair"] = contributions["axial_low"] + contributions["axial_high"]
    cases["all_four_faces"] = rhs
    cases["all_except_axial_high"] = rhs - contributions["axial_high"]
    for name, contribution in cases.items():
        candidate = U0 + DT_FINE*contribution
        attribution_rows.append({
            "case": name, "candidate_rhoe_J_m3": point_rhoe(candidate),
            "candidate_pressure_Pa": (GAMMA-1.0)*point_rhoe(candidate),
            "delta_rhoe_J_m3": point_rhoe(candidate)-point_rhoe(U0),
        })

    # A local low-order control: replace only the dangerous upper axial WENO
    # face by the adjacent two-state LLF flux.  This is an offline diagnostic,
    # not a modification of the simulation.
    Ul, Ur = U[i,j], U[i,j+1]
    qlr = primitive(np.stack([Ul, Ur]))
    alpha_two_state = float(np.max(np.abs(qlr["vel"][:,1]) + qlr["c"]))
    llf_two_state_flux = (
        0.5*(physical_flux(Ul,1) + physical_flux(Ur,1))
        - 0.5*alpha_two_state*(Ur-Ul)
    )
    upper_llf_rhs = -llf_two_state_flux/dz
    replacement = U0 + DT_FINE*(rhs-contributions["axial_high"]+upper_llf_rhs)
    attribution_rows.append({
        "case": "replace_axial_high_by_two_state_LLF",
        "candidate_rhoe_J_m3": point_rhoe(replacement),
        "candidate_pressure_Pa": (GAMMA-1.0)*point_rhoe(replacement),
        "delta_rhoe_J_m3": point_rhoe(replacement)-point_rhoe(U0),
    })
    write_csv(OUT / "failure_internal_energy_attribution.csv", attribution_rows)

    face_summary = []
    for label, face in faces.items():
        face_summary.append({
            "face": label, "llf_alpha_m_s": face["alpha"],
            "weno_epsilon": face["epsilon"], "pressure_face_Pa": face["pressure_face"],
            **{f"flux_{name}": face["flux"][n] for n, name in enumerate(component_names)},
        })
    write_csv(OUT / "failure_face_fluxes.csv", face_summary)

    # Matched field evolution.
    map_steps = [250, 450, 550, 557]
    fig, axes = plt.subplots(2, len(map_steps), figsize=(16, 7), constrained_layout=True)
    norm = LogNorm(vmin=50.0, vmax=2.0e5)
    for col, step in enumerate(map_steps):
        for row, scheme in enumerate(("muscl", "llfweno")):
            if (scheme, step) not in cached:
                cached[(scheme, step)] = load_fields(path_for(scheme, step))
            f, t, re, ze, _, _ = cached[(scheme, step)]
            im = axes[row,col].pcolormesh(ze, re, f["pressure"], shading="auto",
                                         norm=norm, cmap="turbo")
            axes[row,col].plot(zc[TARGET[1]], rc[TARGET[0]], "wo", ms=3, mec="k")
            axes[row,col].set_title(f"{t:.1f} μs")
            axes[row,col].set_xlabel("z [mm]")
            if col == 0:
                axes[row,col].set_ylabel(("HLLC-MUSCL\n" if row == 0 else "LLF-WENO q=2\n") + "r [mm]")
    fig.colorbar(im, ax=axes, label="pressure [Pa]")
    fig.savefig(OUT / "pressure_evolution_muscl_vs_llfweno_q2.png", dpi=190)
    plt.close(fig)

    # Target-cell collapse and comparison.
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    for scheme, label in (("muscl", "HLLC-MUSCL"), ("llfweno", "LLF-WENO-Z5 q=2")):
        rows = sorted((r for r in target_rows if r["scheme"] == scheme), key=lambda r:r["time_us"])
        t = np.array([r["time_us"] for r in rows])
        axes[0].semilogy(t, [r["pressure_Pa"] for r in rows], "o-", ms=3, label=label)
        axes[1].semilogy(t, [r["rhoe_J_m3"] for r in rows], "o-", ms=3, label=label)
    axes[0].axvline(223.0, color="k", ls="--", lw=1)
    axes[1].axvline(223.0, color="k", ls="--", lw=1)
    axes[0].set_ylabel("target pressure [Pa]")
    axes[1].set_ylabel(r"target $\rho e$ [J m$^{-3}$]")
    axes[1].set_xlabel("time [μs]")
    axes[0].legend()
    fig.savefig(OUT / "target_cell_collapse_history.png", dpi=190)
    plt.close(fig)

    # Fine-time local pressure maps around the failing cell.
    local_steps = [550, 553, 555, 557]
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), constrained_layout=True)
    for ax, step in zip(axes, local_steps):
        f, t, re, ze, _, _ = cached[("llfweno", step)]
        slr, slz = slice(48,64), slice(120,145)
        im = ax.pcolormesh(ze[120:146], re[48:65], f["pressure"][slr,slz],
                           shading="auto", norm=LogNorm(vmin=10, vmax=1e5), cmap="turbo")
        ax.plot(zc[TARGET[1]], rc[TARGET[0]], "wo", ms=4, mec="k")
        ax.set_title(f"LLF {t:.1f} μs")
        ax.set_xlabel("z [mm]")
    axes[0].set_ylabel("r [mm]")
    fig.colorbar(im, ax=axes, label="pressure [Pa]")
    fig.savefig(OUT / "llfweno_failure_local_pressure_evolution.png", dpi=190)
    plt.close(fig)

    print(f"analysis written to {OUT}")
    print("offline/inferred RHS relative errors:", [r["relative_error"] for r in flux_rows])


if __name__ == "__main__":
    main()
