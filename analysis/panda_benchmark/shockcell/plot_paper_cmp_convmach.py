#!/usr/bin/env python3
"""Convective Mach number comparison with Panda & Seasholtz (1999) fig10d.

Experimental M_c: phase speed of the SCREECH-frequency (5.4 kHz, helical)
component, from the phase gradient of laser-scattering signals cross-
correlated with a fixed lip microphone; M_c = U_c/c_amb. Strongly modulated
by the standing wave (KH wave + upstream acoustic).

Simulation (3D L5, lip-line probes y = 0.47 D_e, z ~ 0, x/D_e = 0.22-6.0,
dt ~ 5.5 us over the statistics window 6.5-10.16 ms):
(i) broadband eddy convection velocity from the time-domain cross-
    correlation peak of adjacent-probe pressure fluctuations (parabolic lag
    refinement; error bars = spread over 3 time blocks);
(ii) band phase speed at 5.4 +- 1 kHz from the cross-spectral phase,
    U_c(f) = 2*pi*f*dx/dphi (no screech lock exists in the computation).
The two definitions are NOT equivalent to the experimental screech-locked
phase speed; the comparison is qualitative by construction."""
import numpy as np
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
SCRATCH = Path("/tmp/claude-1000/-home-qiaoj-testcerisse-cerisse/"
               "efba2e60-081b-46b5-af00-eb8de7ed5170/scratchpad")
D = 0.0254
C_AMB = np.sqrt(1.4 * 287.0 * 300.0)      # 347.2 m/s (T_amb = 300 K)
DX1 = 1.5875e-3                           # level-1 spacing [m]

# ---- experimental fig10d ----
exp = []
for l in open(HERE.parent.parent.parent / "U_jet/analysis/panda_benchmark/"
              "epaps_raw/fig10d.dat.txt"):
    toks = l.split()
    if len(toks) == 2:
        try:
            exp.append((float(toks[0]), float(toks[1])))
        except ValueError:
            pass
exp = np.array(exp)

# ---- simulation probes ----
P = SCRATCH / "probes_statistics.log"
lines = open(P).read().strip().split("\n")
cols = re.findall(r"pressure\(\((\d+),(\d+),(\d+)\)", lines[0])
rows = [l for l in lines if not l.startswith("time")]
data = np.array([[float(v) for v in l.split(",")] for l in rows])
t = data[:, 0]
# lip-line axial probes: (j,k) = (135,127), unique i ascending
lip = {}
for c, (i, j, k) in enumerate(cols):
    if (j, k) == ("135", "127"):
        lip.setdefault(int(i), c + 1)
II = sorted(lip)
XD = [(i + 0.5) * DX1 / D for i in II]
# uniform resampling
dt = np.median(np.diff(t))
tu = np.arange(t[0], t[-1], dt)
sig = {i: np.interp(tu, t, data[:, lip[i]]) for i in II}
for i in II:
    sig[i] -= sig[i].mean()

def xc_uc(s1, s2, dx, tmin_v=40.0):
    """Broadband: lag of cross-correlation max -> Uc."""
    n = len(s1)
    c = np.correlate(s2, s1, "full")          # positive lag = s2 delayed
    lags = np.arange(-n + 1, n) * dt
    m = (lags > 0) & (lags < dx / tmin_v)
    if not m.any():
        return np.nan
    li = np.where(m)[0]
    k = li[np.argmax(c[li])]
    if 0 < k < len(c) - 1:
        d2 = c[k-1] - 2*c[k] + c[k+1]
        off = 0.5 * (c[k-1] - c[k+1]) / d2 if d2 != 0 else 0.0
    else:
        off = 0.0
    tau = lags[k] + off * dt
    return dx / tau if tau > 0 else np.nan

def band_phase_uc(s1, s2, dx, f0=5400.0, hb=1000.0):
    """Cross-spectral phase in [f0-hb, f0+hb] -> phase speed."""
    n = len(s1)
    w = np.hanning(n)
    F1 = np.fft.rfft(s1 * w); F2 = np.fft.rfft(s2 * w)
    fr = np.fft.rfftfreq(n, dt)
    S12 = np.conj(F1) * F2
    m = (fr > f0 - hb) & (fr < f0 + hb)
    if not m.any():
        return np.nan
    # amplitude-weighted mean of per-frequency phase speeds
    ph = np.angle(S12[m])
    good = ph > 0.05
    if not good.any():
        return np.nan
    uc = 2 * np.pi * fr[m][good] * dx / ph[good]
    wgt = np.abs(S12[m][good])
    uc = uc[(uc > 30) & (uc < 500)]
    wgt = wgt[:len(uc)] if len(wgt) != len(uc) else wgt
    if len(uc) == 0:
        return np.nan
    return float(np.average(uc, weights=np.abs(wgt[:len(uc)])))

pairs = [(a, b) for a, b in zip(II[:-1], II[1:])]
xm, uc_bb, uc_er, uc_ph = [], [], [], []
nb = 3
for a, b in pairs:
    dx = (b - a) * DX1
    xm.append(0.5 * (lipx := ((a + 0.5) * DX1 / D)) + 0.5 * ((b + 0.5) * DX1 / D))
    uc_bb.append(xc_uc(sig[a], sig[b], dx))
    n = len(tu); blk = []
    for q in range(nb):
        sl = slice(q * n // nb, (q + 1) * n // nb)
        s1 = sig[a][sl] - sig[a][sl].mean()
        s2 = sig[b][sl] - sig[b][sl].mean()
        v = xc_uc(s1, s2, dx)
        if np.isfinite(v):
            blk.append(v)
    uc_er.append(np.std(blk, ddof=1) / np.sqrt(len(blk)) if len(blk) > 1
                 else np.nan)
    uc_ph.append(band_phase_uc(sig[a], sig[b], dx))
xm = np.array(xm); uc_bb = np.array(uc_bb)
uc_er = np.array(uc_er); uc_ph = np.array(uc_ph)
for x_, u_, e_, p_ in zip(xm, uc_bb, uc_er, uc_ph):
    print("pair mid x/D=%.2f  Uc_bb=%.0f+-%.0f m/s (Mc=%.2f)  Uc_5.4k=%s"
          % (x_, u_, e_, u_/C_AMB,
             "%.0f (Mc=%.2f)" % (p_, p_/C_AMB) if np.isfinite(p_) else "n/a"))

fig, ax = plt.subplots(figsize=(7.0, 2.9))
ax.plot(exp[:, 0], exp[:, 1], "o", color="k", ms=3.2, mfc="none", mew=0.8,
        label="Exp. screech-locked phase speed (fig10d)")
ax.axhline(np.mean(exp[:, 1]), color="k", lw=0.6, ls=":",
           label="Exp. mean")
ax.errorbar(xm, uc_bb / C_AMB, yerr=uc_er / C_AMB, color="#b3251e",
            marker="o", ms=4, ls="-", lw=0.9, capsize=2, elinewidth=0.7,
            label="3D L5 broadband cross-correlation")
ax.set_xlim(0, 7)
ax.set_ylim(0, 1.8)
ax.set_xlabel("$x/D_e$", labelpad=2)
ax.set_ylabel("$M_c = U_c/c_\\infty$", labelpad=2)
ax.tick_params(length=2.5)
ax.legend(fontsize=6.6, frameon=False, loc="upper right", handlelength=1.5)
fig.tight_layout()
fig.savefig(HERE / "paperfig_cmp_convmach.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_convmach.pdf")
print("wrote paperfig_cmp_convmach")
