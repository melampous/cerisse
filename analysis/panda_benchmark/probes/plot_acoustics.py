#!/usr/bin/env python3
"""Pressure spectra and azimuthal content from the 210-probe records.

Three formulations have a usable record: skew-symmetric+JST, AFD-HLLC+TENO5
and AFD-HLLC+WENO-Z5, each 7.2 ms long. MUSCL+HLLC has 3.6 ms and LLF+TENO5
only 1.25 ms; both are printed in the table but kept out of the figures,
because a Welch estimate over 1.25 ms has a 1.4 kHz bin and no useful number of
averages.

Three things have to be handled before any spectrum is meaningful.

Repeated headers and restarts. The data log is reopened at every restart, so
the header line recurs 6 to 22 times and a few samples straddle a restart with
a non-increasing time. Every line beginning with "time" is dropped, the record
is sorted, and non-increasing steps are removed.

Non-uniform sampling. The probe interval is ten level-1 steps, but the step is
adaptive, so the sampling interval varies by 5 to 50 per cent about its median
depending on the case. The record is resampled onto a uniform axis at the
median interval before any transform. Linear resampling attenuates the top of
the band; at the frequencies of interest here, a few kHz against a Nyquist of
about 90 kHz, that attenuation is negligible.

Spatial averaging. Each probe reports the mean over a 2x2x2 block of level-1
cells, 0.125 D_e on a side. Waves shorter than about twice that are attenuated
by the probe itself, which is a soft ceiling near 40 kHz - well above the
frequencies examined, but it is the reason no claim is made above 20 kHz.

Azimuthal decomposition. The rings are snapped to the Cartesian grid: on the
r = 0.5 rings the spacing runs from 19.25 to 27.0 degrees and the radius from
0.442 to 0.534 D_e. A discrete Fourier transform over probe index would assume
uniform spacing and leak between modes. The modes are therefore obtained by
least squares on the actual angles,

    p(phi_k, t) = sum_{m=-4}^{4} a_m(t) exp(i m phi_k),

solved once per time sample. The range is limited to |m| <= 4 rather than the
nominal 7: with sixteen unevenly spaced points the normal matrix for the full
set is poorly conditioned, and the condition number is printed so the choice
can be checked rather than taken on trust.

The Strouhal axis uses U_j = 414.2 m/s and D_e = 0.0254 m, so St = 1
corresponds to 16.31 kHz. Panda & Seasholtz report screech for the M_j = 1.42
jet at 5.4 kHz, St = 0.33, in a helical m = +-1 mode; that line is drawn for
reference, not as a claim that these records reproduce it.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch

HERE = Path(__file__).resolve().parent
D, U_J = 0.0254, 414.2
FST = U_J/D                       # St = 1 at 16.31 kHz
PREF = 20e-6
SCREECH = 5400.0

CASES = [("s3d_skewjst", "Skew-symmetric + JST", "#009E73", "-", 2.0),
         ("s3d_afdhllc", "AFD-HLLC + TENO5", "#D55E00", "--", 1.5),
         ("s3d_afdhllc_weno", "AFD-HLLC + WENO-Z5", "#D55E00", "-", 1.5)]
SHORT = [("s3d_musclhllc", "MUSCL + HLLC"), ("s3d_teno5", "LLF + TENO5")]

Z = np.load(HERE/"probe_coords.npz")
NAMES = [str(s) for s in Z["names"]]
PHI = np.radians(Z["phi"])
RAD = Z["r"]
XPOS = Z["xyz"][:, 0]/D
IDX = {}
for i, n in enumerate(NAMES):
    IDX.setdefault(n.split("_")[0] if "_" in n else
                   (n.rstrip("abcdefgh") if n[-1].isalpha() and
                    any(c.isdigit() for c in n) else n.rstrip("0123456789")),
                   []).append(i)


def series(tag):
    z = np.load(HERE/("TS_%s.npz" % tag))
    t, p = z["t"], z["p"].astype(np.float64)
    dt = np.median(np.diff(t))
    tu = np.arange(t[0], t[-1], dt)
    pu = np.empty((tu.size, p.shape[1]))
    for k in range(p.shape[1]):
        pu[:, k] = np.interp(tu, t, p[:, k])
    return tu, pu, dt


def psd(p, dt, nseg):
    f, P = welch(p - p.mean(axis=0), fs=1.0/dt, nperseg=min(nseg, p.shape[0]),
                 noverlap=None, detrend="linear", axis=0)
    return f, P


def modes(p, phi, mmax=4):
    """Least-squares azimuthal modes at the true angles."""
    M = np.arange(-mmax, mmax+1)
    A = np.exp(1j*np.outer(phi, M))
    cond = np.linalg.cond(A.conj().T @ A)
    coef = np.linalg.lstsq(A, p.T, rcond=None)[0]      # (2mmax+1, nt)
    return M, coef.T, cond


print("record inventory")
print("%-22s %8s %10s %9s %9s %8s"
      % ("scheme", "samples", "span [ms]", "dt [us]", "df [Hz]", "usable"))
for tag, name in [(c[0], c[1]) for c in CASES] + SHORT:
    z = np.load(HERE/("TS_%s.npz" % tag))
    t = z["t"]
    dt = np.median(np.diff(t))
    n = 512 if t.size > 1000 else 256 if t.size > 500 else 128
    print("%-22s %8d %10.3f %9.3f %9.0f %8s"
          % (name, t.size, (t[-1]-t[0])*1e3, dt*1e6, 1.0/(n*dt),
             "yes" if t.size > 1000 else "no"))

Q = {}
print("\noverall sound pressure level, 0.5 - 20 kHz  [dB re 20 uPa]")
print("%-22s %9s %9s %9s %9s %9s"
      % ("scheme", "s1 r=0.5", "s2 r=0.5", "o2 r=1.0", "n15 r=1.5", "n30 r=3.0"))
for tag, name, col, ls, lw in CASES:
    t, p, dt = series(tag)
    f, P = psd(p, dt, 512)
    Q[tag] = dict(t=t, p=p, dt=dt, f=f, P=P, name=name, col=col, ls=ls, lw=lw)
    band = (f >= 500) & (f <= 20000)
    row = []
    for g in ("s1", "s2", "o2", "n15", "n30"):
        i = IDX[g]
        ms = np.trapezoid(P[band][:, i], f[band], axis=0).mean()
        row.append(10*np.log10(ms/PREF**2))
    print("%-22s %9.1f %9.1f %9.1f %9.1f %9.1f" % (name, *row))

print("\npeak of the ring-averaged spectrum, 1 - 20 kHz")
print("%-22s %-6s %10s %8s %12s" % ("scheme", "ring", "f [kHz]", "St", "dB"))
for tag, *_ in CASES:
    q = Q[tag]
    for g in ("s1", "o2", "n30"):
        i = IDX[g]
        S = q["P"][:, i].mean(axis=1)
        b = (q["f"] >= 1000) & (q["f"] <= 20000)
        k = np.flatnonzero(b)[np.argmax(S[b])]
        print("%-22s %-6s %10.2f %8.3f %12.1f"
              % (q["name"], g, q["f"][k]/1e3, q["f"][k]/FST,
                 10*np.log10(S[k]/PREF**2)))

print("\nazimuthal decomposition on ring s1 (x = 1 D_e, r = 0.5)")
print("%-22s %8s %s" % ("scheme", "cond", "share of 0.5-20 kHz energy by |m|"))
MODE = {}
for tag, *_ in CASES:
    q = Q[tag]
    i = IDX["s1"]
    M, a, cond = modes(q["p"][:, i], PHI[i])
    fm, Pm = psd(a, q["dt"], 512)
    band = (fm >= 500) & (fm <= 20000)
    e = np.trapezoid(np.abs(Pm[band]), fm[band], axis=0)
    tot = e.sum()
    sh = {}
    for mm in range(0, 5):
        sh[mm] = sum(e[k] for k, v in enumerate(M) if abs(v) == mm)/tot
    MODE[tag] = dict(M=M, f=fm, P=np.abs(Pm), cond=cond, share=sh)
    print("%-22s %8.1f %s" % (q["name"], cond,
                              "  ".join("m=%d %4.1f%%" % (k, 100*v)
                                        for k, v in sh.items())))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": False, "ytick.right": True,
                     "legend.frameon": False})

# ---------------------------------------------------- spectra at three radii
fig, AX = plt.subplots(1, 3, figsize=(12.6, 4.3), sharey=True)
fig.subplots_adjust(left=0.062, right=0.986, top=0.845, bottom=0.135,
                    wspace=0.065)
for ax, g, ttl in zip(AX, ("s1", "o2", "n30"),
                      (r"lip line, $x=1$, $r\simeq0.5$",
                       r"$x=2$, $r\simeq1.0$", r"$x=2$, $r\simeq3.0$")):
    for tag, *_ in CASES:
        q = Q[tag]
        S = q["P"][:, IDX[g]].mean(axis=1)
        m = q["f"] > 300
        ax.semilogx(q["f"][m]/1e3, 10*np.log10(S[m]/PREF**2), color=q["col"],
                    ls=q["ls"], lw=q["lw"], zorder=5, label=q["name"])
    ax.axvline(SCREECH/1e3, color="0.35", lw=1.0, ls=":", zorder=3)
    ax.text(SCREECH/1e3*1.06, 60, "Panda screech\n5.4 kHz, $St$ = 0.33",
            fontsize=7.8, color="0.30", rotation=90, va="bottom")
    ax.set_xlim(0.3, 40)
    ax.set_xlabel("$f$ [kHz]")
    ax.set_title(ttl, fontsize=9.8, pad=4)
    ax.grid(True, which="both", ls="--", lw=0.4, color="0.93")
    ax.set_axisbelow(True)
    sec = ax.secondary_xaxis("top", functions=(lambda v: v*1e3/FST,
                                               lambda s: s*FST/1e3))
    sec.set_xlabel(r"$St = fD_e/U_j$", fontsize=9.2, labelpad=3)
    sec.tick_params(labelsize=8.2)
AX[0].set_ylabel(r"PSD  [dB re 20 $\mu$Pa$^2$/Hz]")
AX[0].set_ylim(40, 135)
AX[0].legend(loc="lower left", fontsize=8.4, handlelength=2.2)
fig.savefig(HERE/"acoustic_psd.png", dpi=250)
fig.savefig(HERE/"acoustic_psd.pdf")
plt.close(fig)
print("\nwrote acoustic_psd.png/pdf")

# ------------------------------------------------------- azimuthal modes
fig, AX = plt.subplots(1, 3, figsize=(12.6, 4.3), sharey=True)
fig.subplots_adjust(left=0.062, right=0.986, top=0.845, bottom=0.135,
                    wspace=0.065)
CM = {0: "#000000", 1: "#D55E00", 2: "#0072B2", 3: "#009E73", 4: "#CC79A7"}
for ax, (tag, *_) in zip(AX, CASES):
    q, mo = Q[tag], MODE[tag]
    m = mo["f"] > 300
    for mm in range(0, 5):
        k = [j for j, v in enumerate(mo["M"]) if abs(v) == mm]
        S = mo["P"][:, k].sum(axis=1)
        ax.semilogx(mo["f"][m]/1e3, 10*np.log10(S[m]/PREF**2), color=CM[mm],
                    lw=1.5 if mm < 2 else 1.1, zorder=6-mm,
                    label="$|m|$ = %d   %.1f %%" % (mm, 100*mo["share"][mm]))
    ax.axvline(SCREECH/1e3, color="0.35", lw=1.0, ls=":", zorder=3)
    ax.set_xlim(0.3, 40)
    ax.set_xlabel("$f$ [kHz]")
    ax.set_title(q["name"], fontsize=9.8, pad=4)
    ax.grid(True, which="both", ls="--", lw=0.4, color="0.93")
    ax.set_axisbelow(True)
    ax.legend(loc="lower left", fontsize=8.0, handlelength=1.6)
    sec = ax.secondary_xaxis("top", functions=(lambda v: v*1e3/FST,
                                               lambda s: s*FST/1e3))
    sec.set_xlabel(r"$St = fD_e/U_j$", fontsize=9.2, labelpad=3)
    sec.tick_params(labelsize=8.2)
AX[0].set_ylabel(r"PSD  [dB re 20 $\mu$Pa$^2$/Hz]")
AX[0].set_ylim(40, 135)
fig.suptitle(r"azimuthal modes on the lip-line ring $x=1\,D_e$, "
             r"least squares at the true probe angles", fontsize=10.2, y=0.975)
fig.savefig(HERE/"acoustic_modes.png", dpi=250)
fig.savefig(HERE/"acoustic_modes.pdf")
print("wrote acoustic_modes.png/pdf")
