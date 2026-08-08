#!/usr/bin/env python3
"""Numerical-dissipation evidence set for the three schemes (M_j = 1.42,
3D L5, identical grids): baseline WENO-Z5 NS, TENO5 NS, Euler (WENO-Z5).

Figures
  paperfig_schemes_inst     : instantaneous |omega| De/Uj and numerical
                              schlieren, z=0 plane, 0<x<4, |r|<1.2;
                              teno/euler at t=9.300 ms, baseline at
                              t=9.424 ms (nearest stored frame); common
                              colour scales, common schlieren reference.
  paperfig_schemes_shear    : (a) vorticity thickness delta_omega(x) and
                              half-velocity radius; (b) M=1 supersonic
                              core width. From <u_x>, <u>, <T> (z=0).
  paperfig_schemes_tke      : resolved TKE/Uj^2 fields, common scale.
  paperfig_schemes_spectra  : premultiplied lip-line pressure spectra
                              f*E_pp at x/De = 0.22, 0.47, 0.97, 1.97
                              (y = 0.47 De), Welch; all windows cut to
                              the common [6.5, 9.3] ms; azimuthal probe
                              duplicates averaged where available.
Table: schemes_centreline_metrics.md (x_fc from raw <rho>; x_n and S_n from
<p>; per-cell peak/trough amplitudes from <rho>)."""
import re
import numpy as np
from scipy.signal import find_peaks, welch
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
SCR = Path("/tmp/claude-1000/-home-qiaoj-testcerisse-cerisse/"
           "efba2e60-081b-46b5-af00-eb8de7ed5170/scratchpad")
D = 0.0254
UJ = 414.2
RHO_J = 1.6413
P_AMB = 99780.0
DX1 = 1.5875e-3
COL = {"baseline": "#b3251e", "teno": "#167c72", "euler": "#6a51a3"}
LAB = {"baseline": "WENO-Z5 NS", "teno": "TENO5 NS",
       "euler": "Euler (WENO-Z5)"}
SCHEMES = ["baseline", "teno", "euler"]

# ================= Fig 1: instantaneous pair =========================
XI, RI = 4.0, 1.2
inst = {}
for s in SCHEMES:
    d = np.load(HERE / ("SCHEME_inst_%s.npz" % s))
    om, g, dx = d["om"], d["grho"], float(d["dx"])
    nx, ny = om.shape
    x = (np.arange(nx) + 0.5) * dx / D
    y = (np.arange(ny) - ny // 2 + 0.5) * dx / D
    mi = x <= XI
    mj = np.abs(y) <= RI
    inst[s] = (x[mi], y[mj], om[np.ix_(mi, mj)] * D / UJ,
               g[np.ix_(mi, mj)], float(d["time"]) * 1e3)
gref = np.nanpercentile(np.concatenate(
    [inst[s][3].ravel() for s in SCHEMES]), 99.9)

fig, axs = plt.subplots(3, 2, figsize=(7.2, 6.4), sharex=True,
                        sharey=True)
for r, s in enumerate(SCHEMES):
    x, y, om, g, tms = inst[s]
    ext = [x[0], x[-1], y[0], y[-1]]
    im0 = axs[r, 0].imshow(om.T, origin="lower", extent=ext, vmin=0,
                           vmax=12, cmap="viridis", aspect="equal",
                           interpolation="nearest", rasterized=True)
    im1 = axs[r, 1].imshow(np.exp(-8.0 * g.T / gref), origin="lower",
                           extent=ext, vmin=0, vmax=1, cmap="gray",
                           aspect="equal", interpolation="nearest",
                           rasterized=True)
    for c, tcol in ((0, "w"), (1, "k")):
        axs[r, c].text(0.02, 0.94,
                       "(%s) %s, $t=%.3f$ ms" % (chr(97 + 2 * r + c),
                                                 LAB[s], tms),
                       transform=axs[r, c].transAxes, ha="left",
                       va="top", fontsize=7.5, color=tcol,
                       bbox=dict(facecolor="black" if tcol == "w"
                                 else "white", edgecolor="none",
                                 alpha=0.45, pad=1.2))
    axs[r, 0].set_ylabel("$r/D_e$", labelpad=1)
for ax in axs[-1]:
    ax.set_xlabel("$x/D_e$", labelpad=2)
for ax in axs.flat:
    ax.set_xlim(0, XI); ax.set_ylim(-RI, RI)
    ax.tick_params(length=2.5, labelsize=8)
fig.subplots_adjust(left=0.07, right=0.985, top=0.94, bottom=0.16,
                    wspace=0.06, hspace=0.10)
for x0, im, lab in ((0.10, im0, "$|\\omega|\\,D_e/U_j$"),
                    (0.56, im1,
                     "$\\exp(-8\\,|\\nabla\\rho|/|\\nabla\\rho|_{99.9})$")):
    cax = fig.add_axes([x0, 0.065, 0.35, 0.018])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label(lab, fontsize=8.5, labelpad=2)
    cb.ax.tick_params(length=2.5, labelsize=8)
fig.savefig(HERE / "paperfig_schemes_inst.png", dpi=300)
fig.savefig(HERE / "paperfig_schemes_inst.pdf", dpi=300)
print("wrote paperfig_schemes_inst")

# ================= Fig 2: shear-layer scales + core width ============
def ysmooth(F, n=3):
    k = np.ones(n) / n
    return np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, F)

shear = {}
for s in SCHEMES:
    d = np.load(HERE / ("SCHEME_mean_%s.npz" % s))
    x, y = d["x"], d["y"]
    ux = ysmooth(d["uxm"].astype(np.float64))
    M = d["M"].astype(np.float64)
    dy = (y[1] - y[0]) * D
    dudy = np.gradient(ux, dy, axis=1)
    umax = ux.max(1); umin = ux.min(1)
    gmax = np.abs(dudy).max(1)
    dw = (umax - umin) / np.maximum(gmax, 1e-9) / D
    # half-velocity radius, averaged over the two sides
    rh = np.full(len(x), np.nan)
    jc = np.argmin(np.abs(y))
    for i in range(len(x)):
        uh = 0.5 * (umax[i] + umin[i])
        pos = np.where(ux[i, jc:] < uh)[0]
        neg = np.where(ux[i, :jc][::-1] < uh)[0]
        if len(pos) and len(neg):
            rh[i] = 0.5 * (y[jc + pos[0]] + abs(y[jc - 1 - neg[0]]))
    # M = 1 supersonic core width
    wc = np.array([(np.ptp(y[(M[i] >= 1.0) & (np.abs(y) <= 1.0)])
                    if ((M[i] >= 1.0) & (np.abs(y) <= 1.0)).any() else 0.0)
                   for i in range(len(x))])
    shear[s] = (x, dw, rh, wc)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.9))
for ax in (a1, a2):
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87")
    ax.tick_params(length=2.5)
for s in SCHEMES:
    x, dw, rh, wc = shear[s]
    m = (x > 0.15) & (x <= 4.0)
    a1.plot(x[m], dw[m], color=COL[s], lw=1.2, label=LAB[s])
    a1.plot(x[m], rh[m], color=COL[s], lw=0.9, ls="--")
    a2.plot(x[m], wc[m], color=COL[s], lw=1.2, label=LAB[s])
a1.set_xlabel("$x/D_e$", labelpad=2)
a1.set_ylabel("$\\delta_\\omega/D_e$ (solid), $r_{1/2}/D_e$ (dashed)",
              labelpad=2)
a1.set_xlim(0, 4); a1.set_ylim(0, 0.8)
a1.legend(fontsize=6.8, frameon=False, loc="upper left",
          handlelength=1.5)
a1.text(0.02, 0.03, "(a)", transform=a1.transAxes, fontsize=9)
a2.set_xlabel("$x/D_e$", labelpad=2)
a2.set_ylabel("$M=1$ core width $/D_e$", labelpad=2)
a2.set_xlim(0, 7); a2.set_ylim(0, 1.2)
a2.legend(fontsize=6.8, frameon=False, loc="lower left",
          handlelength=1.5)
a2.text(0.02, 0.93, "(b)", transform=a2.transAxes, fontsize=9)
fig.tight_layout(w_pad=1.2)
fig.savefig(HERE / "paperfig_schemes_shear.png", dpi=300)
fig.savefig(HERE / "paperfig_schemes_shear.pdf")
print("wrote paperfig_schemes_shear")

# ================= Fig 3: resolved TKE composite =====================
fig, axes = plt.subplots(3, 1, figsize=(8.2, 6.2), sharex=True,
                         gridspec_kw={"hspace": 0.12})
tkes = []
for s in SCHEMES:
    d = np.load(HERE / ("SCHEME_mean_%s.npz" % s))
    tkes.append((d["x"], d["y"], d["tke"] / UJ**2))
VMAX = float(np.nanpercentile(
    np.concatenate([t[2].ravel() for t in tkes]), 99.5))
im = None
for ax, s, (sx, sy, tk) in zip(axes, SCHEMES, tkes):
    im = ax.pcolormesh(sx, sy, tk.T, cmap="magma", vmin=0, vmax=VMAX,
                       shading="auto", rasterized=True)
    ax.set_xlim(0, 7.2); ax.set_ylim(-1.0, 1.0)
    ax.set_ylabel(r"$r/D_e$")
    ax.text(0.012, 0.06, LAB[s], transform=ax.transAxes, fontsize=9,
            va="bottom", color="w",
            bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.55,
                      boxstyle="round,pad=0.25"))
axes[-1].set_xlabel(r"$x/D_e$")
cax = fig.add_axes([0.25, 0.058, 0.5, 0.015])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label(r"resolved TKE $/U_j^2$", fontsize=9)
cb.ax.tick_params(labelsize=8)
fig.subplots_adjust(left=0.075, right=0.985, top=0.985, bottom=0.15)
fig.savefig(HERE / "paperfig_schemes_tke.png", dpi=220)
fig.savefig(HERE / "paperfig_schemes_tke.pdf", dpi=220)
print("wrote paperfig_schemes_tke  vmax=%.4f" % VMAX)

# ================= Table: centreline metrics =========================
def boxfilt(a, x, wD):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w//2] = a[:w//2]; out[-(w//2):] = a[-(w//2):]
    return out

def refine(x, y, i):
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    d2 = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if d2 == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/d2*(x[1]-x[0])

lines = ["| scheme | x_fc | x1 | x2 | x3 | x4 | S1 | S2 | S3 | A1 | A2 | A3 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for s in SCHEMES:
    d = np.load(HERE / ("SCHEME_mean_%s.npz" % s))
    x, y = d["x"], d["y"]
    jc = np.argsort(np.abs(y))[:2]
    p = boxfilt(d["pm"][:, jc].mean(1) / P_AMB, x, 0.05)
    dr = np.load(HERE / ("SIM_%s_rhoMEAN_z0.npz"
                         % ("L5" if s == "baseline" else s)))
    jr = np.argsort(np.abs(dr["y"]))[:2]
    rho_raw = dr["rho"][:, jr].mean(1) / RHO_J
    rho = boxfilt(rho_raw, dr["x"], 0.05)
    xr = dr["x"]
    # x_fc from the steepest rise of the raw mean density between the first
    # expansion minimum and the subsequent principal compression maximum.
    drhodx = np.gradient(rho_raw, xr)
    expansion = np.where((xr >= 0.35) & (xr <= 1.0))[0]
    imin = expansion[np.nanargmin(rho_raw[expansion])]
    compression = np.where((np.arange(xr.size) > imin) & (xr <= 1.65))[0]
    imax = compression[np.nanargmax(rho_raw[compression])]
    m0 = np.arange(imin, imax + 1)
    xfc = refine(xr, drhodx, m0[np.nanargmax(drhodx[m0])])
    dist = max(1, int(round(0.5 / (x[1] - x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, p, -np.inf), height=1.0,
                       prominence=0.03, distance=dist)
    xn = [refine(x, p, i) for i in pk[:4]]
    Sn = np.diff(xn)
    # per-cell peak/trough amplitude from <rho>
    distr = max(1, int(round(0.5 / (xr[1] - xr[0]))))
    rpk, _ = find_peaks(np.where(xr > 0.5, rho, -np.inf),
                        prominence=0.02, distance=distr)
    rtr, _ = find_peaks(np.where(xr > 0.5, -rho, -np.inf),
                        prominence=0.02, distance=distr)
    rpk = rpk[xr[rpk] > 0.6]        # drop mask-boundary artefact
    xpk = [refine(xr, rho, i) for i in rpk[:4]]
    A = []
    for n in range(min(3, len(xpk))):
        prev = [xr[i] for i in rtr if xr[i] < xpk[n]]
        A.append(np.interp(xpk[n], xr, rho)
                 - (np.interp(prev[-1], xr, rho) if prev else np.nan))
    lines.append("| %s | %.3f | %s | %s | %s |"
                 % (LAB[s], xfc,
                    " | ".join("%.3f" % v for v in (xn + [np.nan]*4)[:4]),
                    " | ".join("%.3f" % v for v in (list(Sn) + [np.nan]*3)[:3]),
                    " | ".join("%.3f" % v for v in (A + [np.nan]*3)[:3])))
open(HERE / "schemes_centreline_metrics.md", "w").write(
    "\n".join(lines) + "\n")
print("\n".join(lines))

# ================= Fig 4: premultiplied lip-line spectra =============
T0, T1 = 6.5e-3, 9.3e-3
STA = [3, 7, 15, 31]          # level-1 i-index -> x/De 0.22 0.47 0.97 1.97
def load_probes(path):
    txt = open(path).read().strip().split("\n")
    cols = re.findall(r"pressure\(\((\d+),(\d+),(\d+)\)", txt[0])
    rows = [l for l in txt[1:] if not l.startswith("time")]
    dat = np.array([[float(v) for v in l.split(",")] for l in rows])
    t = dat[:, 0]
    m = (t >= T0) & (t <= T1)
    series = {}
    for c, (i, j, k) in enumerate(cols):
        series.setdefault(int(i) if (j, k) == ("135", "127") else None, [])
    out = {}
    for c, (i, j, k) in enumerate(cols):
        ii = int(i)
        if ii in STA:
            out.setdefault(ii, []).append(dat[m, c + 1])
    return t[m], out

PROBES = {"baseline": SCR / "probes_statistics.log",
          "teno": SCR / "probesB_teno.log",
          "euler": SCR / "probesB_euler.log"}
# resolved cutoff: convective wavelength 8*dx at L5 (dx=99.2um),
# Uc ~ 0.6 Uj -> f_cut = Uc/(8 dx). Above this, do not interpret.
DXL5 = 99.2e-6
FCUT = 0.6 * UJ / (8 * DXL5) / 1e3     # kHz
FHI = 30.0                              # broadband-tail band start [kHz]

fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.7), sharey=True)
for ax in axs:
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", lw=0.5, color="0.87", which="both")
    ax.tick_params(length=2.5, labelsize=8)
tail = {s: [] for s in SCHEMES}
for s in SCHEMES:
    t, series = load_probes(PROBES[s])
    dt = np.median(np.diff(t))
    fs = 1.0 / dt
    for ax, ii in zip(axs, STA):
        if ii not in series:
            tail[s].append(np.nan); continue
        Ps = []
        for sig in series[ii]:
            sig = sig - sig.mean()
            f, P = welch(sig, fs=fs, window="hann", nperseg=256,
                         noverlap=128)
            Ps.append(P)
        P = np.mean(Ps, 0)
        m = f > 0
        ax.loglog(f[m] / 1e3, P[m] / P_AMB**2, color=COL[s], lw=1.0,
                  label=LAB[s])
        fk = f / 1e3
        band = (fk >= FHI) & (fk <= min(FCUT, fk.max()))
        tail[s].append((np.trapezoid if hasattr(np,"trapezoid") else np.trapz)(P[band], f[band]) / P_AMB**2)
for ax, ii in zip(axs, STA):
    ax.axvspan(min(FCUT, 90), 95, color="0.85", alpha=0.5, lw=0)
    ax.axvline(FCUT, color="0.4", lw=0.7, ls=":")
    ax.set_xlabel("$f$ [kHz]", labelpad=2)
    ax.set_title("$x/D_e=%.2f$" % ((ii + 0.5) * DX1 / D), fontsize=8.5)
    ax.set_xlim(1, 95); ax.set_ylim(1e-10, 3e-6)
axs[0].set_ylabel("$E_{pp}/p_\\infty^2$ [Hz$^{-1}$]", labelpad=2)
axs[0].legend(fontsize=6.2, frameon=False, loc="lower left",
              handlelength=1.3)
fig.tight_layout(w_pad=0.5)
fig.savefig(HERE / "paperfig_schemes_spectra.png", dpi=300)
fig.savefig(HERE / "paperfig_schemes_spectra.pdf")
print("wrote paperfig_schemes_spectra")
print("f_cut (8dx, L5) = %.0f kHz; broadband-tail band [%.0f, %.0f] kHz"
      % (FCUT, FHI, min(FCUT, 90)))
print("high-freq tail energy (x/De=0.22/0.47/0.97/1.97):")
for s in SCHEMES:
    print("  %-16s %s" % (LAB[s],
          "  ".join("%.2e" % v for v in tail[s])))
