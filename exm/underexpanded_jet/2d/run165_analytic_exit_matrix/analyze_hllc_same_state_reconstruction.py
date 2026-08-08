#!/usr/bin/env python3
"""Reconstruct one HLLC face with MUSCL and AFD-WENO on the same state.

The source state is the stage-1 Forward-Euler stencil immediately before the
AFD-HLLC-WENO failure at L1 face (4,171), dir=1.  Values are parsed from the
existing all-fluid stage trace, so this diagnostic does not evolve a new flow
solution and does not involve IBM.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_hllc_muscl_vs_afd_weno"
TRACE = RESULTS / "afd_fallback_off" / "cell_trace_radius3_fail.log"
OUTDIR = RESULTS / "same_state_reconstruction"
GAMMA = 1.4
FACE_J = 171
RADIAL_I = 4
DT_FINE = 2.0e-7
DX = 1.0612966666666667e-3

CELL_RE = re.compile(
    r"label=rk2_stage1_fe.*cell=\((?P<i>-?\d+),(?P<j>-?\d+)\).*"
    r"values=\(rho=(?P<rho>[-+0-9.eE]+),mx=(?P<mx>[-+0-9.eE]+),"
    r"my=(?P<my>[-+0-9.eE]+),mz=(?P<mz>[-+0-9.eE]+),"
    r"E=(?P<E>[-+0-9.eE]+)\)"
)


def load_stencil() -> dict[int, dict[str, float]]:
    cells: dict[int, dict[str, float]] = {}
    for line in TRACE.read_text().splitlines():
        match = CELL_RE.search(line)
        if not match or int(match["i"]) != RADIAL_I:
            continue
        j = int(match["j"])
        if not 168 <= j <= 174 or j in cells:
            continue
        values = {name: float(match[name]) for name in ("rho", "mx", "my", "mz", "E")}
        rho = values["rho"]
        kinetic = 0.5 * (values["mx"] ** 2 + values["my"] ** 2 + values["mz"] ** 2) / rho
        pressure = (GAMMA - 1.0) * (values["E"] - kinetic)
        cells[j] = {
            **values,
            "ur": values["mx"] / rho,
            "un": values["my"] / rho,
            "ut2": values["mz"] / rho,
            "p": pressure,
            "c": math.sqrt(GAMMA * pressure / rho),
        }
    expected = set(range(168, 175))
    if set(cells) != expected:
        raise RuntimeError(f"incomplete stencil: got {sorted(cells)}, expected {sorted(expected)}")
    return cells


def mc_limiter(left: float, right: float) -> float:
    centered = 0.5 * (left + right)
    sign = math.copysign(1.0, centered)
    limited = 2.0 * min(abs(left), abs(right)) if left * right >= 0.0 else 0.0
    return sign * min(limited, abs(centered))


def muscl_slope(cells: dict[int, dict[str, float]], j: int) -> list[float]:
    qm, q, qp = cells[j - 1], cells[j], cells[j + 1]
    c = q["c"]
    rho = q["rho"]
    dpl = q["p"] - qm["p"]
    dpr = qp["p"] - q["p"]
    dul = q["un"] - qm["un"]
    dur = qp["un"] - q["un"]
    d0 = mc_limiter(0.5 * dpl / c - 0.5 * rho * dul,
                    0.5 * dpr / c - 0.5 * rho * dur)
    d1 = mc_limiter(q["rho"] - qm["rho"] - dpl / c**2,
                    qp["rho"] - q["rho"] - dpr / c**2)
    d2 = mc_limiter(0.5 * dpl / c + 0.5 * rho * dul,
                    0.5 * dpr / c + 0.5 * rho * dur)
    d3 = mc_limiter(q["ur"] - qm["ur"], qp["ur"] - q["ur"])
    d4 = mc_limiter(q["ut2"] - qm["ut2"], qp["ut2"] - q["ut2"])
    return [d0, d1, d2, d3, d4]


def muscl_state(q: dict[str, float], slope: list[float], side: str) -> dict[str, float]:
    factor = 0.5 if side == "left" else -0.5
    d0, d1, d2, d3, d4 = slope
    rho = q["rho"] + factor * ((d0 + d2) / q["c"] + d1)
    un = q["un"] + factor * ((d2 - d0) / q["rho"])
    pressure = q["p"] + factor * (d0 + d2) * q["c"]
    return {
        "rho": max(rho, 1.0e-40),
        "un": un,
        "ut1": q["ur"] + factor * d3,
        "ut2": q["ut2"] + factor * d4,
        "p": max(pressure, 1.0e-40),
    }


def characteristic(q: dict[str, float]) -> list[float]:
    scale = math.sqrt(GAMMA * q["rho"] * q["p"])
    return [
        q["rho"] * (1.0 - 1.0 / GAMMA),
        0.5 * (q["p"] + scale * q["un"]),
        0.5 * (q["p"] - scale * q["un"]),
        q["ur"],
        q["ut2"],
        GAMMA,
    ]


def weno_right(samples: list[float]) -> dict[str, object]:
    s0, s1, s2, s3, s4 = samples
    beta = [
        13.0 / 12.0 * (s0 - 2.0 * s1 + s2) ** 2
        + 0.25 * (s0 - 4.0 * s1 + 3.0 * s2) ** 2,
        13.0 / 12.0 * (s1 - 2.0 * s2 + s3) ** 2
        + 0.25 * (s1 - s3) ** 2,
        13.0 / 12.0 * (s2 - 2.0 * s3 + s4) ** 2
        + 0.25 * (3.0 * s2 - 4.0 * s3 + s4) ** 2,
    ]
    candidate = [
        (3.0 * s0 - 10.0 * s1 + 15.0 * s2) / 8.0,
        (-s1 + 6.0 * s2 + 3.0 * s3) / 8.0,
        (3.0 * s2 + 6.0 * s3 - s4) / 8.0,
    ]
    tau = abs(beta[0] - beta[2])
    optimal = [1.0, 10.0, 5.0]
    alpha = [optimal[k] * (1.0 + (tau / (1.0e-40 + beta[k])) ** 2) for k in range(3)]
    total = sum(alpha)
    weight = [value / total for value in alpha]
    value = sum(weight[k] * candidate[k] for k in range(3))
    return {"samples": samples, "beta": beta, "candidate": candidate, "weight": weight, "value": value}


def weno_left(samples: list[float]) -> dict[str, object]:
    return weno_right(list(reversed(samples)))


def weno_states(cells: dict[int, dict[str, float]]) -> tuple[dict[str, float], dict[str, float], list[dict], list[dict]]:
    wc = {j: characteristic(cells[j]) for j in cells}
    left_details: list[dict] = []
    right_details: list[dict] = []
    wl: list[float] = []
    wr: list[float] = []
    for component in range(6):
        left = weno_right([wc[j][component] for j in range(168, 173)])
        right = weno_left([wc[j][component] for j in range(169, 174)])
        wl.append(float(left["value"]))
        wr.append(float(right["value"]))
        left_details.append(left)
        right_details.append(right)

    def inverse(w: list[float]) -> dict[str, float]:
        gamma = w[5]
        pressure = w[1] + w[2]
        rho = w[0] / (1.0 - 1.0 / gamma)
        un = (w[1] - w[2]) / math.sqrt(gamma * rho * pressure)
        return {"rho": rho, "un": un, "ut1": w[3], "ut2": w[4], "p": pressure}

    return inverse(wl), inverse(wr), left_details, right_details


def frozen_characteristic_states(
    cells: dict[int, dict[str, float]],
) -> tuple[dict[str, float], dict[str, float]]:
    """Diagnostic frozen linear primitive-characteristic reconstruction."""
    adjacent_left = cells[FACE_J - 1]
    adjacent_right = cells[FACE_J]
    rho_ref = 0.5 * (adjacent_left["rho"] + adjacent_right["rho"])
    c_ref = 0.5 * (adjacent_left["c"] + adjacent_right["c"])
    transformed = {
        j: [
            0.5 * q["p"] / c_ref - 0.5 * rho_ref * q["un"],
            q["rho"] - q["p"] / c_ref**2,
            0.5 * q["p"] / c_ref + 0.5 * rho_ref * q["un"],
            q["ur"],
            q["ut2"],
        ]
        for j, q in cells.items()
    }

    def reconstruct(side: str) -> dict[str, float]:
        indices = range(168, 173) if side == "left" else range(169, 174)
        values = []
        for component in range(5):
            samples = [transformed[j][component] for j in indices]
            detail = weno_right(samples) if side == "left" else weno_left(samples)
            values.append(float(detail["value"]))
        pressure = c_ref * (values[0] + values[2])
        return {
            "rho": values[1] + (values[0] + values[2]) / c_ref,
            "un": (values[2] - values[0]) / rho_ref,
            "ut1": values[3],
            "ut2": values[4],
            "p": pressure,
        }

    return reconstruct("left"), reconstruct("right")


def hllc(
    left: dict[str, float],
    right: dict[str, float],
    sound_left: float | None = None,
    sound_right: float | None = None,
) -> dict[str, float]:
    def augmented(state: dict[str, float]) -> dict[str, float]:
        rho, p = state["rho"], state["p"]
        specific = p / ((GAMMA - 1.0) * rho)
        specific += 0.5 * (state["un"] ** 2 + state["ut1"] ** 2 + state["ut2"] ** 2)
        return {**state, "c": math.sqrt(GAMMA * p / rho), "e": specific}

    l, r = augmented(left), augmented(right)
    if sound_left is not None:
        l["c"] = sound_left
    if sound_right is not None:
        r["c"] = sound_right
    sl = min(l["un"] - l["c"], r["un"] - r["c"])
    sr = max(l["un"] + l["c"], r["un"] + r["c"])
    ratio = math.sqrt(r["rho"] / l["rho"])
    uroe = (l["un"] + r["un"] * ratio) / (1.0 + ratio)
    croe = (l["c"] + r["c"] * ratio) / (1.0 + ratio)
    sl = min(sl, uroe - croe)
    sr = max(sr, uroe + croe)
    denom = l["rho"] * (sl - l["un"]) - r["rho"] * (sr - r["un"])
    sstar = (
        r["p"] - l["p"]
        + l["rho"] * l["un"] * (sl - l["un"])
        - r["rho"] * r["un"] * (sr - r["un"])
    ) / denom

    if sl > 0.0:
        q = l
        mass = q["rho"] * q["un"]
        mn = q["rho"] * q["un"] ** 2 + q["p"]
        energy = q["un"] * (q["rho"] * q["e"] + q["p"])
    elif sr < 0.0:
        q = r
        mass = q["rho"] * q["un"]
        mn = q["rho"] * q["un"] ** 2 + q["p"]
        energy = q["un"] * (q["rho"] * q["e"] + q["p"])
    else:
        q, speed = (l, sl) if sstar >= 0.0 else (r, sr)
        frac = (speed - q["un"]) / (speed - sstar) - 1.0
        rho_star_ratio = frac + 1.0
        mass = q["rho"] * q["un"] + speed * q["rho"] * frac
        mn = (
            q["rho"] * q["un"] ** 2
            + q["p"]
            + speed * q["rho"] * (rho_star_ratio * sstar - q["un"])
        )
        energy = q["un"] * (q["rho"] * q["e"] + q["p"]) + speed * q["rho"] * (
            frac * q["e"]
            + rho_star_ratio
            * (sstar - q["un"])
            * (sstar + q["p"] / q["rho"] / (speed - q["un"]))
        )
    return {
        "mass": mass,
        "normal_momentum": mn,
        "energy": energy,
        "sl": sl,
        "sr": sr,
        "sstar": sstar,
        "amax": max(abs(sl), abs(sr)),
        "face_cfl": DT_FINE * max(abs(sl), abs(sr)) / DX,
    }


def main() -> None:
    cells = load_stencil()
    muscl_left = muscl_state(cells[170], muscl_slope(cells, 170), "left")
    muscl_right = muscl_state(cells[171], muscl_slope(cells, 171), "right")
    weno_left, weno_right_state, left_details, right_details = weno_states(cells)
    frozen_left, frozen_right = frozen_characteristic_states(cells)
    # Riemann.h passes the adjacent cell-centred sound speeds to its MUSCL
    # HLLC call.  In this supersonic branch the flux itself is nevertheless
    # exactly the right reconstructed physical flux.
    muscl_flux = hllc(muscl_left, muscl_right, cells[170]["c"], cells[171]["c"])
    weno_flux = hllc(weno_left, weno_right_state)
    frozen_flux = hllc(frozen_left, frozen_right)

    OUTDIR.mkdir(exist_ok=True)
    with (OUTDIR / "six_point_raw_and_acoustic_state.csv").open("w", newline="") as stream:
        fieldnames = [
            "j", "rho", "p", "un", "ur", "sound", "w0_entropy_density",
            "w_plus", "w_minus", "u_radial", "u_theta", "gamma",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for j in range(168, 174):
            q = cells[j]
            w = characteristic(q)
            writer.writerow({
                "j": j,
                "rho": f"{q['rho']:.17g}",
                "p": f"{q['p']:.17g}",
                "un": f"{q['un']:.17g}",
                "ur": f"{q['ur']:.17g}",
                "sound": f"{q['c']:.17g}",
                "w0_entropy_density": f"{w[0]:.17g}",
                "w_plus": f"{w[1]:.17g}",
                "w_minus": f"{w[2]:.17g}",
                "u_radial": f"{w[3]:.17g}",
                "u_theta": f"{w[4]:.17g}",
                "gamma": f"{w[5]:.17g}",
            })
    with (OUTDIR / "same_state_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "scheme", "side", "rho", "un", "ut1", "p", "wave_speed", "face_cfl", "mass_flux", "momentum_flux", "energy_flux"
        ])
        writer.writeheader()
        for scheme, left, right, flux in (
            ("HLLC-MUSCL", muscl_left, muscl_right, muscl_flux),
            ("AFD-HLLC-WENO-Z5", weno_left, weno_right_state, weno_flux),
            ("frozen-linear-WENO diagnostic", frozen_left, frozen_right, frozen_flux),
        ):
            for side, state in (("left", left), ("right", right)):
                writer.writerow({
                    "scheme": scheme,
                    "side": side,
                    **{name: f"{state[name]:.17g}" for name in ("rho", "un", "ut1", "p")},
                    "wave_speed": f"{flux['amax']:.17g}",
                    "face_cfl": f"{flux['face_cfl']:.17g}",
                    "mass_flux": f"{flux['mass']:.17g}",
                    "momentum_flux": f"{flux['normal_momentum']:.17g}",
                    "energy_flux": f"{flux['energy']:.17g}",
                })

    names = ["entropy_density", "w_plus", "w_minus", "u_radial", "u_theta", "gamma"]
    for side, details in (("left", left_details), ("right", right_details)):
        with (OUTDIR / f"weno_{side}_candidates.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=[
                "component", "sample0", "sample1", "sample2", "sample3", "sample4",
                "candidate0", "candidate1", "candidate2", "beta0", "beta1", "beta2",
                "weight0", "weight1", "weight2", "reconstructed",
            ])
            writer.writeheader()
            for name, detail in zip(names, details):
                writer.writerow({
                    "component": name,
                    **{f"sample{k}": f"{detail['samples'][k]:.17g}" for k in range(5)},
                    **{f"candidate{k}": f"{detail['candidate'][k]:.17g}" for k in range(3)},
                    **{f"beta{k}": f"{detail['beta'][k]:.17g}" for k in range(3)},
                    **{f"weight{k}": f"{detail['weight'][k]:.17g}" for k in range(3)},
                    "reconstructed": f"{detail['value']:.17g}",
                })

    paired_candidates = []
    for candidate in range(3):
        w0 = float(right_details[0]["candidate"][candidate])
        wplus = float(right_details[1]["candidate"][candidate])
        wminus = float(right_details[2]["candidate"][candidate])
        density = w0 / (1.0 - 1.0 / GAMMA)
        pressure = wplus + wminus
        velocity = (wplus - wminus) / math.sqrt(GAMMA * density * pressure)
        paired_candidates.append((density, pressure, velocity))
    with (OUTDIR / "weno_right_paired_candidate_states.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["candidate", "rho", "p", "un"])
        writer.writeheader()
        for candidate, (density, pressure, velocity) in enumerate(paired_candidates):
            writer.writerow({
                "candidate": candidate,
                "rho": f"{density:.17g}",
                "p": f"{pressure:.17g}",
                "un": f"{velocity:.17g}",
            })

    js = list(range(168, 175))
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.6), constrained_layout=True)
    face_x = FACE_J - 0.5
    axes[0, 0].plot(js, [cells[j]["p"] for j in js], "o-k", label="stage-1 cell centres")
    for label, state, marker, colour in (
        ("MUSCL right", muscl_right, "s", "tab:blue"),
        ("current WENO right", weno_right_state, "X", "tab:red"),
        ("frozen-linear WENO right", frozen_right, "D", "tab:green"),
    ):
        axes[0, 0].scatter(face_x, state["p"], marker=marker, s=75, color=colour, label=label, zorder=4)
    axes[0, 0].set_ylabel("pressure [Pa]")
    axes[0, 0].set_xlabel("axial cell index j (face at 170.5)")
    axes[0, 0].set_yscale("symlog", linthresh=100.0)
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.25)

    axes[0, 1].plot(js, [cells[j]["un"] for j in js], "o-k", label="stage-1 cell centres")
    for label, state, marker, colour in (
        ("MUSCL right", muscl_right, "s", "tab:blue"),
        ("current WENO right", weno_right_state, "X", "tab:red"),
        ("frozen-linear WENO right", frozen_right, "D", "tab:green"),
    ):
        axes[0, 1].scatter(face_x, state["un"], marker=marker, s=75, color=colour, label=label, zorder=4)
    axes[0, 1].set_ylabel("axial velocity [m/s]")
    axes[0, 1].set_xlabel("axial cell index j (face at 170.5)")
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(alpha=0.25)

    candidates = range(3)
    axes[1, 0].bar([k - 0.18 for k in candidates], [v[1] for v in paired_candidates], width=0.36,
                   color="tab:orange", label="paired-candidate pressure")
    weight_axis = axes[1, 0].twinx()
    weight_axis.plot(candidates, right_details[1]["weight"], "o-", color="tab:blue", label=r"$w_+$ weight")
    weight_axis.plot(candidates, right_details[2]["weight"], "s-", color="tab:red", label=r"$w_-$ weight")
    axes[1, 0].set_xlabel("right-state candidate stencil")
    axes[1, 0].set_ylabel("candidate pressure [Pa]")
    weight_axis.set_ylabel("nonlinear weight")
    axes[1, 0].grid(alpha=0.25)
    handles1, labels1 = axes[1, 0].get_legend_handles_labels()
    handles2, labels2 = weight_axis.get_legend_handles_labels()
    axes[1, 0].legend(handles1 + handles2, labels1 + labels2, fontsize=8, loc="upper left")

    names = ["MUSCL", "current\nWENO", "frozen-linear\nWENO"]
    fluxes = [muscl_flux, weno_flux, frozen_flux]
    axes[1, 1].bar(names, [abs(flux["energy"]) for flux in fluxes],
                   color=["tab:blue", "tab:red", "tab:green"])
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel(r"$|F_{\rho E}|$ [W/m$^2$]")
    axes[1, 1].grid(axis="y", alpha=0.25)
    for position, flux in enumerate(fluxes):
        axes[1, 1].text(position, abs(flux["energy"]) * 1.15,
                        f"CFL$_f$={flux['face_cfl']:.3f}", ha="center", fontsize=8)
    fig.suptitle("Same-state HLLC input diagnosis at face (4,171), RK2 stage 2")
    fig.savefig(OUTDIR / "hllc_same_state_reconstruction.png", dpi=180)
    plt.close(fig)

    print("same AFD stage-1 FE stencil, face=(4,171), dir=1")
    for name, left, right, flux in (
        ("HLLC-MUSCL", muscl_left, muscl_right, muscl_flux),
        ("AFD-HLLC-WENO-Z5", weno_left, weno_right_state, weno_flux),
        ("frozen-linear-WENO diagnostic", frozen_left, frozen_right, frozen_flux),
    ):
        print(name)
        print("  left ", left)
        print("  right", right)
        print("  flux ", flux)
    print("WENO/MUSCL |energy flux| ratio =", abs(weno_flux["energy"] / muscl_flux["energy"]))
    print("WENO right acoustic candidates/weights:")
    for component in (1, 2):
        detail = right_details[component]
        print(component, "candidate=", detail["candidate"], "beta=", detail["beta"], "weight=", detail["weight"])


if __name__ == "__main__":
    main()
