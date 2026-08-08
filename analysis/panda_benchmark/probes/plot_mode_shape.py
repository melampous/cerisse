#!/usr/bin/env python3
"""The shape of the mode that carries the 5.6-5.8 kHz peak.

The band-integrated decomposition in plot_acoustics.py is dominated by m = 0
simply because the broadband floor is. Restricting to a band around the peak
reverses it completely: the peak is carried by m = +-1. This figure asks what
kind of m = +-1 it is, which the sign of m answers and the modulus does not.

    a_{+1} and a_{-1} of equal strength   ->  a standing, flapping mode; the
                                             pattern oscillates in one plane
    one of them dominant                  ->  a spinning helix, rotating in a
                                             fixed sense

The mode coefficients come from least squares on the true probe angles, which
matters here more than for the modulus: the rings are grid-snapped, spacing
19.25 to 27.0 degrees, and an assumed-uniform transform mixes +m into -m and
would manufacture a false balance between them.

The time-mean and the rms fields cannot show any of this. A statistically
stationary round jet is axisymmetric in every moment however strong its
instantaneous helical content, so the mean is blind to the mode by
construction; the information lives in the cross-spectral phase around the
ring, which is what is plotted here.

The record is 7.2 ms, so at 5.7 kHz there are about 41 cycles. That is enough
to establish which mode carries the band and the sense of rotation; it is not
enough to quote a stable amplitude for either.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch, csd

HERE = Path(__file__).resolve().parent
D, U_J = 0.0254, 414.2
FST, PREF = U_J/D, 20e-6
NSEG = 1024
FPK = {"s3d_skewjst": 5618.0, "s3d_afdhllc": 5779.0,
       "s3d_afdhllc_weno": 5706.0}
CASES = [("s3d_skewjst", "Skew-symmetric + JST", "#009E73"),
         ("s3d_afdhllc", "AFD-HLLC + TENO5", "#D55E00"),
         ("s3d_afdhllc_weno", "AFD-HLLC + WENO-Z5", "#E8A33D")]

Z = np.load(HERE/"probe_coords.npz")
NAMES = [str(s) for s in Z["names"]]
PHI, RAD = np.radians(Z["phi"]), Z["r"]
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
    return tu, np.stack([np.interp(tu, t, p[:, k])
                         for k in range(p.shape[1])], 1), dt


def coeffs(p, phi, mmax=4):
    M = np.arange(-mmax, mmax+1)
    A = np.exp(1j*np.outer(phi, M))
    return M, np.linalg.lstsq(A, p.T, rcond=None)[0].T


print("signed azimuthal content in the peak band (peak +- 700 Hz), ring o2")
print("%-22s %8s %8s %8s %8s %10s"
      % ("scheme", "m=-1", "m=0", "m=+1", "+1/-1", "verdict"))
Q = {}
for tag, name, col in CASES:
    t, p, dt = series(tag)
    i = IDX["o2"]
    M, a = coeffs(p[:, i], PHI[i])
    f, P = welch(a, fs=1.0/dt, nperseg=NSEG, detrend="linear", axis=0)
    P = np.abs(P)
    b = np.abs(f - FPK[tag]) <= 700.0
    e = P[b].sum(axis=0)
    ep = {int(v): e[k] for k, v in enumerate(M)}
    ratio = ep[1]/ep[-1]
    Q[tag] = dict(t=t, p=p, dt=dt, f=f, P=P, M=M, a=a, name=name, col=col,
                  e=ep, ratio=ratio)
    tot = ep[1]+ep[-1]+ep[0]
    print("%-22s %7.1f%% %7.1f%% %7.1f%% %8.2f %10s"
          % (name, 100*ep[-1]/tot, 100*ep[0]/tot, 100*ep[1]/tot, ratio,
             "spinning" if ratio > 2 or ratio < 0.5 else "standing/flapping"))

# Coherence needs averaging: with nperseg = 1024 there is one segment in a
# 1300-sample record and the estimate is identically 1 by construction. The
# phase is taken at the fine resolution, where it is well defined, and the
# coherence separately at nperseg = 256, which gives about nine averages and a
# 707 Hz bin - coarse, but an honest number instead of a tautological one.
NCOH = 256
print("\ncross-spectral phase around the ring at the peak, ring o2")
print("(a pure m mode advances by 360*m degrees over one turn;"
      " coherence at nperseg=%d, ~%.0f averages)" % (NCOH, 2*1300/NCOH-1))
print("%-22s %14s %10s %12s"
      % ("scheme", "slope [deg/turn]", "implied m", "coherence"))
for tag, *_ in CASES:
    q = Q[tag]
    i = IDX["o2"]
    ref = i[int(np.argmin(np.abs(PHI[i])))]
    ph, co = [], []
    for j in i:
        f, Pxy = csd(q["p"][:, ref], q["p"][:, j], fs=1.0/q["dt"],
                     nperseg=NSEG, detrend="linear")
        k = int(np.argmin(np.abs(f - FPK[tag])))
        ph.append(np.angle(Pxy[k]))
        fc, Pc = csd(q["p"][:, ref], q["p"][:, j], fs=1.0/q["dt"],
                     nperseg=NCOH, detrend="linear")
        _, Pxx = welch(q["p"][:, ref], fs=1.0/q["dt"], nperseg=NCOH,
                       detrend="linear")
        _, Pyy = welch(q["p"][:, j], fs=1.0/q["dt"], nperseg=NCOH,
                       detrend="linear")
        kc = int(np.argmin(np.abs(fc - FPK[tag])))
        co.append(np.abs(Pc[kc])**2/(Pxx[kc]*Pyy[kc]))
    ph = np.unwrap(np.array(ph)[np.argsort(PHI[i])])
    ang = np.sort(PHI[i])
    sl = np.polyfit(ang, ph, 1)[0]
    q["ph"], q["co"], q["ang"] = ph, np.array(co)[np.argsort(PHI[i])], ang
    cc = np.array(co)
    print("%-22s %14.0f %10.2f   %.2f (min %.2f)"
          % (q["name"], np.degrees(sl*2*np.pi), sl, cc.mean(), cc.min()))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": False, "ytick.right": True,
                     "legend.frameon": False})
fig = plt.figure(figsize=(12.8, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.255,
                      left=0.055, right=0.988, top=0.885, bottom=0.135)
a1 = fig.add_subplot(gs[0, 0])
a2 = fig.add_subplot(gs[0, 1])
a3 = fig.add_subplot(gs[0, 2], projection="polar")

for tag, *_ in CASES:
    q = Q[tag]
    m = (q["f"] >= 2000) & (q["f"] <= 12000)
    for sgn, lsty in ((+1, "-"), (-1, "--")):
        k = int(np.flatnonzero(q["M"] == sgn)[0])
        a1.plot(q["f"][m]/1e3, 10*np.log10(q["P"][m, k]/PREF**2),
                color=q["col"], ls=lsty, lw=1.6 if sgn > 0 else 1.2, zorder=5,
                label=("%s,  $m=%+d$" % (q["name"], sgn)))
    a1.axvline(FPK[tag]/1e3, color=q["col"], lw=0.8, ls=":", zorder=3)
a1.set_xlim(2, 12)
a1.set_xlabel("$f$ [kHz]")
a1.set_ylabel(r"PSD  [dB re 20 $\mu$Pa$^2$/Hz]")
a1.set_title(r"(a)  $m=+1$ against $m=-1$, ring $o2$", fontsize=9.8, pad=4)
a1.legend(loc="lower left", fontsize=7.4, handlelength=2.0, ncol=1)

for tag, *_ in CASES:
    q = Q[tag]
    a2.plot(np.degrees(q["ang"]), np.degrees(q["ph"] - q["ph"][0]), "o-",
            color=q["col"], ms=5, mec="k", mew=0.5, lw=1.3, zorder=5,
            label="%s  (coh %.2f)" % (q["name"], q["co"].mean()))
a2.plot([-180, 180], [-360, 360], "-", color="0.55", lw=1.2, zorder=3,
        label="$m=+1$")
a2.plot([-180, 180], [360, -360], "--", color="0.55", lw=1.2, zorder=3,
        label="$m=-1$")
a2.set_xlim(-180, 180)
a2.set_ylim(-420, 420)
a2.set_xticks([-180, -90, 0, 90, 180])
a2.set_yticks([-360, -180, 0, 180, 360])
a2.set_xlabel(r"$\varphi$ [deg]")
a2.set_ylabel("cross-spectral phase [deg]")
a2.set_title("(b)  phase around the ring at the peak", fontsize=9.8, pad=4)
a2.legend(loc="upper left", fontsize=7.4, handlelength=1.8)

th = np.linspace(0, 2*np.pi, 361)
for tag, *_ in CASES:
    q = Q[tag]
    kp = int(np.flatnonzero(q["M"] == +1)[0])
    km = int(np.flatnonzero(q["M"] == -1)[0])
    b = np.abs(q["f"] - FPK[tag]) <= 700.0
    Ap = np.sqrt(q["P"][b, kp].sum())
    Am = np.sqrt(q["P"][b, km].sum())
    amp = np.abs(Ap*np.exp(1j*th) + Am*np.exp(-1j*th))
    a3.plot(th, amp/amp.max(), color=q["col"], lw=1.8, zorder=5,
            label=q["name"])
a3.set_rlim(0, 1.15)
a3.set_rticks([0.5, 1.0])
a3.set_title(r"(c)  $|p'|(\varphi)$ in the $y$-$z$ plane" "\n"
             r"from $m=\pm1$, normalised", fontsize=9.6, pad=14)
a3.tick_params(labelsize=7.8)
a3.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30), fontsize=7.4)

for ax in (a1, a2):
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"mode_shape.png", dpi=250)
fig.savefig(HERE/"mode_shape.pdf")
print("\nwrote mode_shape.png/pdf")
