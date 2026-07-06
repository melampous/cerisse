#!/usr/bin/env python3
"""Driver: 1-D BC reflection benchmark -> all figures + JSON for the report.
Uses the vectorised core bc_1d_fast."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import bc_1d_fast as B

FIG = "figures"
out = {}
bcs = [("FOE", "foextrap"), ("DIR", "Dirichlet"), ("NSC", "NSCBC-algebraic"),
       ("LODI", "NSCBC-LODI ($\\sigma$=0.25)")]
colors = {"FOE": "tab:green", "DIR": "tab:orange", "NSC": "tab:red", "LODI": "tab:blue"}

# 1. x-t diagrams (M=0) ------------------------------------------------------
print("[1] x-t diagrams")
xt_store = {}
fig, axes = plt.subplots(1, 4, figsize=(20, 5.2), sharey=True)
for ax, (bc, name) in zip(axes, bcs):
    res = B.run_pulse(bc, M=0.0, record_xt=True, t_end=13.0)
    xt_store[bc] = res
    vlim = res["amp"]
    im = ax.pcolormesh(res["x"], res["xt_t"], res["xt"], cmap="RdBu_r",
                       vmin=-vlim, vmax=vlim, shading="auto", rasterized=True)
    ax.set_xlabel("x"); ax.set_title(name, fontsize=12); ax.set_xlim(0, 10)
    ax.text(0.05, 0.95, f"R = {res['R_energy']:.3f}", transform=ax.transAxes,
            va="top", fontsize=12, fontweight="bold",
            bbox=dict(boxstyle="round", fc="white", alpha=0.85))
axes[0].set_ylabel("t  (acoustic units, $c_0$=1)")
fig.suptitle("Acoustic-pulse x–t diagram, quiescent medium (M=0): incident pulse → right boundary → reflection",
             fontsize=13)
fig.colorbar(im, ax=axes, shrink=0.85, pad=0.01).set_label("$\\delta p$")
plt.savefig(f"{FIG}/xt_diagrams_M0.png", dpi=140, bbox_inches="tight"); plt.close()

# 2. sensor traces (M=0) -----------------------------------------------------
print("[2] sensor traces")
fig, ax = plt.subplots(figsize=(11, 5))
R_M0 = {}
for bc, name in bcs:
    res = xt_store[bc]; R_M0[bc] = res["R_energy"]
    ax.plot(res["t"], res["p_sensor"] / res["amp"], color=colors[bc], lw=1.6,
            label=f"{name}  (R={res['R_energy']:.3f})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("t  (acoustic units)"); ax.set_ylabel("$\\delta p/\\delta p_0$ at sensor x=6")
ax.set_title("Sensor pressure history (M=0): first lobe = incident, later lobe = reflected")
ax.legend(); ax.grid(alpha=0.3)
plt.savefig(f"{FIG}/sensor_traces_M0.png", dpi=150, bbox_inches="tight"); plt.close()

# 3. R vs Mach ---------------------------------------------------------------
print("[3] R vs Mach")
machs = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
Rtab = {bc: [] for bc, _ in bcs}
for M in machs:
    for bc, _ in bcs:
        Rtab[bc].append(B.run_pulse(bc, M=M, t_end=13.0)["R_energy"])
    print(f"   M={M:.2f}  " + "  ".join(f"{bc}={Rtab[bc][-1]:.3f}" for bc, _ in bcs))
fig, ax = plt.subplots(figsize=(10, 6))
for bc, name in bcs:
    ax.plot(machs, Rtab[bc], "o-", color=colors[bc], lw=1.8, ms=5, label=name)
ax.set_xlabel("background Mach number  $M_n$ (normal to boundary)")
ax.set_ylabel("reflection coefficient  $R$ (energy)")
ax.set_title("Acoustic reflection vs normal Mach number at a subsonic outflow")
ax.legend(); ax.grid(alpha=0.3); ax.axhline(0, color="k", lw=0.5)
ax.annotate("lateral far-field of all 3\nvalidation cases ($M_n\\approx0$)",
            xy=(0.01, Rtab["NSC"][0]), xytext=(0.25, 0.5),
            arrowprops=dict(arrowstyle="->", color="red"), color="red", fontsize=10)
plt.savefig(f"{FIG}/R_vs_mach.png", dpi=150, bbox_inches="tight"); plt.close()

# 4. ROOT CAUSE: R vs sigma --------------------------------------------------
print("[4] R vs sigma")
sigmas = [0.0, 0.05, 0.1, 0.15, 0.25, 0.5, 1.0, 2.0, 5.0, 15.0, 50.0, 200.0]
R_sigma = [B.run_pulse("LODI", M=0.1, sigma=s, t_end=13.0)["R_energy"] for s in sigmas]
R_nsc_M01 = B.run_pulse("NSC", M=0.1, t_end=13.0)["R_energy"]
for s, r in zip(sigmas, R_sigma):
    print(f"   sigma={s:7.2f}  R={r:.3f}")
print(f"   NSC algebraic R={R_nsc_M01:.3f}")
fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogx(np.array(sigmas) + 1e-3, R_sigma, "o-", color="tab:blue", lw=1.8, ms=6,
            label="NSCBC-LODI (proper, relaxation $\\sigma$)")
ax.axhline(R_nsc_M01, color="tab:red", ls="--", lw=2,
           label=f"Cerisse algebraic NSCBC code (R={R_nsc_M01:.3f})")
ax.axvline(0.25, color="g", ls=":", lw=1.5, label="recommended $\\sigma$=0.25")
ax.set_xlabel("relaxation coefficient  $\\sigma$"); ax.set_ylabel("reflection coefficient $R$ (M=0.1)")
ax.set_title("ROOT CAUSE: algebraic code = the $\\sigma\\to\\infty$ (hard $p=p_\\infty$) limit of NSCBC")
ax.legend(); ax.grid(alpha=0.3, which="both")
plt.savefig(f"{FIG}/root_cause_sigma.png", dpi=150, bbox_inches="tight"); plt.close()

# 5. pressure drift ----------------------------------------------------------
print("[5] pressure drift")
def run_drift(bc, M=0.1, dp0=5e-3, L=10.0, N=400, t_end=40.0, sigma=0.25):
    r0 = 1.0; p0 = 1.0 / B.GAM; c0 = B.csound(r0, p0); u0 = M * c0; dx = L / N
    rho = np.full(N, r0); u = np.full(N, u0); p = np.full(N, p0 + dp0)
    U = B.prim2cons_arr(rho, u, p)
    Uinf = B.prim2cons_arr(np.array([r0]), np.array([u0]), np.array([p0]))[:, 0]
    dt = 0.4 * dx / (c0 + abs(u0)); ts, pm = [], []; t = 0.0
    while t < t_end:
        if t + dt > t_end: dt = t_end - t
        U = B.step(U, dt, dx, bc, Uinf, sigma, L); t += dt
        _, _, pf = B.cons2prim(U); ts.append(t); pm.append(np.mean(pf) - p0)
    return np.array(ts), np.array(pm)
fig, ax = plt.subplots(figsize=(11, 5)); drift = {}
for bc, name in bcs:
    t, pm = run_drift(bc); drift[bc] = float(pm[-1] / 5e-3)
    ax.plot(t, pm / 5e-3, color=colors[bc], lw=1.7, label=f"{name} (residual={drift[bc]:.2f})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("t"); ax.set_ylabel("domain-mean $\\delta p/\\delta p_0$")
ax.set_title("Mean-pressure drift (M=0.1): does the BC relax the box back to $p_\\infty$?")
ax.legend(); ax.grid(alpha=0.3)
plt.savefig(f"{FIG}/pressure_drift.png", dpi=150, bbox_inches="tight"); plt.close()

# 6. entropy advection -------------------------------------------------------
print("[6] entropy advection")
ent_store = {}; spurious = {}
fig, axes = plt.subplots(1, 4, figsize=(20, 5.2), sharey=True)
for ax, (bc, name) in zip(axes, bcs):
    res = B.run_pulse(bc, M=0.5, kind="entropy", amp=2e-2, t_end=9.0, record_xt=True)
    ent_store[bc] = res
    sp = float(np.max(np.abs(res["p_sensor"])) / (res["amp"] * res["p0"]))
    spurious[bc] = sp
    vlim = res["amp"] * res["p0"] * 0.4
    im = ax.pcolormesh(res["x"], res["xt_t"], res["xt"], cmap="RdBu_r",
                       vmin=-vlim, vmax=vlim, shading="auto", rasterized=True)
    ax.set_xlabel("x"); ax.set_title(f"{name}\nspurious $\\delta p$={sp:.3f}", fontsize=11)
    ax.set_xlim(0, 10)
    print(f"   {bc:5s} spurious={sp:.4f}")
axes[0].set_ylabel("t")
fig.suptitle("Entropy-wave advection (M=0.5): $p$ field should stay 0 — any $\\delta p$ is spurious outflow noise",
             fontsize=13)
fig.colorbar(im, ax=axes, shrink=0.85, pad=0.01).set_label("$\\delta p$ (spurious)")
plt.savefig(f"{FIG}/entropy_xt.png", dpi=140, bbox_inches="tight"); plt.close()

# summary --------------------------------------------------------------------
out["R_M0"] = R_M0
out["R_vs_mach"] = {"mach": machs, **{bc: Rtab[bc] for bc, _ in bcs}}
out["sigma"] = sigmas; out["R_sigma"] = R_sigma; out["R_nsc_M01"] = R_nsc_M01
out["drift_residual"] = drift; out["spurious_entropy"] = spurious
json.dump(out, open("bc_study_summary.json", "w"), indent=2)
print("\nDONE wrote bc_study_summary + 6 figures")
print("R(M=0):  " + "  ".join(f"{bc}={R_M0[bc]:.3f}" for bc, _ in bcs))
print("drift:   " + "  ".join(f"{bc}={drift[bc]:+.2f}" for bc, _ in bcs))
print("entropy: " + "  ".join(f"{bc}={spurious[bc]:.3f}" for bc, _ in bcs))
