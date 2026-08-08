#!/usr/bin/env python3
"""Discrete inflow-flux audit, m142 3D (faithful to prob.h lines 146-175):
  r0 = r_jet - 3*delta;  f = 0.5*(1 - tanh((r-r0)/delta));  u = U_e*f
  T = T0 - u^2/(2*Cp);   rho = p_jet/(R*T);   hole mask r < r_jet.
Cell-centre quadrature at each AMR level vs near-exact radial quadrature.
Metrics: mdot, Jx_mom = int rho u^2 dA (pressure term excluded: identical
A*p_jet for a given mask), Cd_num = discrete/analytic at SAME delta, and
inter-delta discriminability at each level. fig54 + table."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
R, gam = 287.0, 1.4
Cp = gam*R/(gam-1.0)
T0 = 297.15
T_e = 247.625
p_jet = 99780.0*3.273446/1.89293
r_jet = 0.0127
D = 2*r_jet
U_e = np.sqrt(gam*R*T_e)

def props(f):
    u = U_e*f
    T = T0 - 0.5*u*u/Cp
    rho = p_jet/(R*T)
    return u, rho

def fprof(rr, delta):
    if delta == 0.0:
        return np.ones_like(rr)
    r0 = r_jet - 3.0*delta
    return 0.5*(1.0 - np.tanh((rr - r0)/delta))

def discrete(h, delta):
    N = int(np.ceil(1.2*r_jet/h)) + 2
    i = np.arange(-N, N+1)
    yc = (i + 0.5)*h
    Y, Z = np.meshgrid(yc, yc, indexing="ij")
    rr = np.sqrt(Y**2 + Z**2)
    m = rr < r_jet
    f = fprof(rr[m], delta)
    u, rho = props(f)
    dA = h*h
    return np.sum(rho*u)*dA, np.sum(rho*u*u)*dA, np.count_nonzero(m)*dA

def analytic(delta, n=200000):
    r = (np.arange(n) + 0.5)*(r_jet/n)
    f = fprof(r, delta)
    u, rho = props(f)
    dr = r_jet/n
    md = np.sum(rho*u*2*np.pi*r)*dr
    jx = np.sum(rho*u*u*2*np.pi*r)*dr
    # momentum thickness at the lip (planar def., f(1-f) weight)
    th = np.sum(f*(1.0-f))*dr
    return md, jx, th

DELTAS = [0.0, 25e-6, 50e-6, 100e-6, 150e-6, 254e-6]
LEVELS = [("L0", D/8), ("L1", D/16), ("L2", D/32), ("L3", D/64),
          ("L4", D/128), ("L5", D/256)]
ana = {d: analytic(d) for d in DELTAS}

print("Analytic reference (per delta):")
print("delta[um]   mdot[kg/s]   Jx_mom[N]   theta[um]   r_half=r0[mm]")
for d in DELTAS:
    md, jx, th = ana[d]
    print(f"{d*1e6:8.0f} {md:12.5f} {jx:11.3f} {th*1e6:10.2f} {(r_jet-3*d)*1e3:12.4f}")

print("\nDiscrete/analytic ratios (Cd_num for mdot; same-delta):")
hdr = "level  h[um] " + "".join(f"  d={d*1e6:>3.0f}" for d in DELTAS)
print(hdr)
res = {}
for name, h in LEVELS:
    row = []
    for d in DELTAS:
        md, jx, A = discrete(h, d)
        res[(name, d)] = (md, jx, A)
        row.append(md/ana[d][0])
    print(f"{name:5s} {h*1e6:6.0f} " + "".join(f" {v:6.4f}" for v in row))

print("\nInter-delta discriminability: (mdot(d)-mdot(0))/mdot(0) [%]")
print(hdr)
for src in ["ana"] + [n for n, _ in LEVELS]:
    row = []
    for d in DELTAS:
        if src == "ana":
            v = (ana[d][0]-ana[0.0][0])/ana[0.0][0]*100
        else:
            v = (res[(src, d)][0]-res[(src, 0.0)][0])/res[(src, 0.0)][0]*100
        row.append(v)
    print(f"{src:12s} " + "".join(f" {v:6.2f}" for v in row))

print("\nSame for momentum flux Jx_mom: (Jx(d)-Jx(0))/Jx(0) [%]")
for src in ["ana"] + [n for n, _ in LEVELS]:
    row = []
    for d in DELTAS:
        if src == "ana":
            v = (ana[d][1]-ana[0.0][1])/ana[0.0][1]*100
        else:
            v = (res[(src, d)][1]-res[(src, 0.0)][1])/res[(src, 0.0)][1]*100
        row.append(v)
    print(f"{src:12s} " + "".join(f" {v:6.2f}" for v in row))

# ---- fig54 ----
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
dd = np.array(DELTAS)*1e6
cols = plt.cm.Blues(np.linspace(0.35, 1.0, len(LEVELS)))
for (name, h), c in zip(LEVELS, cols):
    axs[0].plot(dd, [res[(name, d)][0]/ana[d][0] for d in DELTAS],
                "o-", color=c, lw=1.4, ms=4, label=f"{name} ({h*1e6:.0f} µm)")
axs[0].axhline(1.0, color="k", lw=0.8, ls=":")
axs[0].set_xlabel("$\\delta$ [µm]"); axs[0].set_ylabel("$\\dot m_{disc}/\\dot m_{exact}$")
axs[0].set_title("(a) quadrature fidelity at same $\\delta$")
axs[0].legend(fontsize=8); axs[0].grid(alpha=0.25, lw=0.4)
axs[1].plot(dd, [(ana[d][1]/ana[0.0][1]-1)*100 for d in DELTAS],
            "k^-", lw=1.8, ms=6, label="analytic")
for (name, h), c in zip([LEVELS[1], LEVELS[5]], [cols[1], cols[5]]):
    axs[1].plot(dd, [(res[(name, d)][1]/res[(name, 0.0)][1]-1)*100 for d in DELTAS],
                "o--", color=c, lw=1.4, ms=5, label=f"{name} discrete")
axs[1].set_xlabel("$\\delta$ [µm]")
axs[1].set_ylabel("$[J_x(\\delta)-J_x(0)]/J_x(0)$ [%]")
axs[1].set_title("(b) momentum-flux deficit vs $\\delta$: is the effect representable?")
axs[1].legend(fontsize=8); axs[1].grid(alpha=0.25, lw=0.4)
fig.suptitle("m142 inflow BC audit — cell-centre sampling of the tanh profile (prob.h model)")
fig.tight_layout()
fig.savefig(HERE / "fig54_inflow_flux_audit.png", dpi=160)
print("\nwrote fig54_inflow_flux_audit.png")
