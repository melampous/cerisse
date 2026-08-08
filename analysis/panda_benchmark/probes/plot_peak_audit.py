#!/usr/bin/env python3
"""What the 5-6 kHz peak is, and what it is not.

A narrow rise near St = 0.35 in a screeching-jet spectrum has two competing
explanations and they are separable with the probe array as laid out.

  screech          a feedback tone. Frequency set by the shock-cell spacing
                   through Powell's relation, the SAME at every observation
                   angle, phase-locked, narrow, with harmonics.
  broadband shock  the peak frequency shifts with observation angle, because
  associated noise  it comes from a wavepacket-shock interference whose
                   condition contains the observer direction.

Four tests are applied.

1. Powell scaling. f_s = U_c / [L_s (1 + U_c/a_inf)] with U_c = 0.7 U_j. Each
   formulation has its own measured shock-cell spacing L_s, taken from the
   centreline mean density in plot_scheme_metrics.py, so the prediction is
   different for each and the test has content: if the measured peak tracks
   each formulation's own L_s, the peak is locked to the shock cells.

2. Angle independence. The n15 and n30 rings sit at the same x = 2 D_e but at
   r = 1.5 and 3.0, so they view the source at 36.9 and 56.3 degrees from the
   axis. A screech tone keeps its frequency; broadband shock noise does not.

3. Tone-to-background ratio and bandwidth. Measured against a median-filtered
   background over a decade around the peak. A tone is narrow and stands well
   clear; a wavepacket hump is broad and low.

4. Which azimuthal mode carries the peak, evaluated in a band around the peak
   rather than over the whole spectrum.

Resolution is the binding constraint and is stated with every number. The
record is 7.2 ms. Welch with 1024-sample segments gives a 177 Hz bin and about
three averages, so the tone can be located to roughly a bin but its width and
its peak level carry a large uncertainty. Nothing below is a claim that these
records reproduce screech; the tests say which explanation the data is
consistent with, at this record length.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch, csd, medfilt

HERE = Path(__file__).resolve().parent
D, U_J, A_INF = 0.0254, 414.2, 340.0
FST = U_J/D
PREF = 20e-6
SCREECH_EXP = 5400.0

# shock-cell spacing of each formulation, from the centreline mean density
LS = {"s3d_skewjst": 1.147, "s3d_afdhllc": 1.127, "s3d_afdhllc_weno": 1.115}
CASES = [("s3d_skewjst", "Skew-symmetric + JST", "#009E73", "-", 2.0),
         ("s3d_afdhllc", "AFD-HLLC + TENO5", "#D55E00", "--", 1.5),
         ("s3d_afdhllc_weno", "AFD-HLLC + WENO-Z5", "#D55E00", "-", 1.5)]

Z = np.load(HERE/"probe_coords.npz")
NAMES = [str(s) for s in Z["names"]]
PHI = np.radians(Z["phi"])
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
    pu = np.stack([np.interp(tu, t, p[:, k]) for k in range(p.shape[1])], 1)
    return tu, pu, dt


def spec(p, dt, nseg):
    f, P = welch(p, fs=1.0/dt, nperseg=min(nseg, p.shape[0]),
                 detrend="linear", axis=0)
    return f, P


def peak(f, S, lo=3000.0, hi=9000.0):
    """Peak against a median-filtered background, with a parabolic refinement."""
    bg = medfilt(10*np.log10(np.maximum(S, 1e-30)), 31)
    b = (f >= lo) & (f <= hi)
    k = np.flatnonzero(b)[np.argmax(S[b])]
    db = 10*np.log10(S) - bg
    fpk = f[k]
    if 0 < k < len(f)-1:
        y0, y1, y2 = np.log(S[k-1]), np.log(S[k]), np.log(S[k+1])
        d2 = y0 - 2*y1 + y2
        if d2 < 0:
            fpk = f[k] + 0.5*(y0-y2)/d2*(f[1]-f[0])
    # half-power width about the peak
    half = S[k]/2.0
    a = k
    while a > 0 and S[a] > half:
        a -= 1
    c = k
    while c < len(S)-1 and S[c] > half:
        c += 1
    return fpk, db[k], f[c]-f[a], 10*np.log10(S[k]/PREF**2)


print("resolution of this estimate")
NSEG = 1024
for tag, name, *_ in CASES:
    t, p, dt = series(tag)
    print("   %-22s %d samples, %.3f ms, bin %.0f Hz, ~%.1f averages"
          % (name, t.size, (t[-1]-t[0])*1e3, 1.0/(NSEG*dt),
             2.0*t.size/NSEG - 1))

Q = {}
print("\n1. Powell scaling   f_s = U_c / [L_s (1 + U_c/a_inf)],  U_c = 0.7 U_j")
print("%-22s %8s %11s %11s %10s %9s"
      % ("scheme", "L_s [D_e]", "f_Powell", "f_measured", "diff", "f*L_s/U_j"))
for tag, name, col, ls, lw in CASES:
    t, p, dt = series(tag)
    f, P = spec(p, dt, NSEG)
    S = P[:, IDX["o2"]].mean(axis=1)
    fpk, tbr, bw, lev = peak(f, S)
    ls_m = LS[tag]*D
    uc = 0.7*U_J
    fpow = uc/(ls_m*(1.0 + uc/A_INF))
    Q[tag] = dict(t=t, p=p, dt=dt, f=f, P=P, name=name, col=col, lsty=ls,
                  lw=lw, fpk=fpk, tbr=tbr, bw=bw, lev=lev, fpow=fpow)
    print("%-22s %8.3f %10.0f Hz %10.0f Hz %+9.1f%% %9.4f"
          % (name, LS[tag], fpow, fpk, 100*(fpk/fpow - 1), fpk*ls_m/U_J))

print("\n   experiment: Panda & Seasholtz 5400 Hz;  Powell with L_s = 1.160 D_e"
      " gives %.0f Hz" % (0.7*U_J/(1.160*D*(1+0.7*U_J/A_INF))))

print("\n2. Angle independence   n15 at 36.9 deg, n30 at 56.3 deg from the axis")
print("%-22s %12s %12s %10s" % ("scheme", "f(n15)", "f(n30)", "shift"))
for tag, *_ in CASES:
    q = Q[tag]
    fs = []
    for g in ("n15", "n30"):
        S = q["P"][:, IDX[g]].mean(axis=1)
        fs.append(peak(q["f"], S)[0])
    print("%-22s %9.0f Hz %9.0f Hz %+9.1f%%"
          % (q["name"], fs[0], fs[1], 100*(fs[1]/fs[0]-1)))
    q["fn15"], q["fn30"] = fs

print("\n3. Tone-to-background and width, ring o2 (bin is %.0f Hz)"
      % (1.0/(NSEG*Q["s3d_skewjst"]["dt"])))
print("%-22s %10s %12s %12s %10s"
      % ("scheme", "peak dB", "above bg", "-3dB width", "Q = f/dw"))
for tag, *_ in CASES:
    q = Q[tag]
    print("%-22s %10.1f %10.1f dB %9.0f Hz %10.1f"
          % (q["name"], q["lev"], q["tbr"], q["bw"], q["fpk"]/max(q["bw"], 1)))

print("\n4. azimuthal content in a band around the peak (peak +- 700 Hz),"
      " ring o2")


def modes(p, phi, mmax=4):
    M = np.arange(-mmax, mmax+1)
    A = np.exp(1j*np.outer(phi, M))
    return M, np.linalg.lstsq(A, p.T, rcond=None)[0].T


print("%-22s %s" % ("scheme", "share of band energy by |m|"))
for tag, *_ in CASES:
    q = Q[tag]
    i = IDX["o2"]
    M, a = modes(q["p"][:, i], PHI[i])
    fm, Pm = spec(a, q["dt"], NSEG)
    b = np.abs(fm - q["fpk"]) <= 700.0
    e = np.abs(Pm[b]).sum(axis=0)
    sh = {mm: sum(e[k] for k, v in enumerate(M) if abs(v) == mm)/e.sum()
          for mm in range(5)}
    q["share"] = sh
    print("%-22s %s" % (q["name"],
                        "  ".join("m=%d %4.1f%%" % (k, 100*v)
                                  for k, v in sh.items())))

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": False, "ytick.right": True,
                     "legend.frameon": False})
fig, AX = plt.subplots(1, 3, figsize=(12.8, 4.4))
fig.subplots_adjust(left=0.058, right=0.986, top=0.845, bottom=0.135,
                    wspace=0.235)

a1, a2, a3 = AX
for tag, *_ in CASES:
    q = Q[tag]
    S = q["P"][:, IDX["o2"]].mean(axis=1)
    m = (q["f"] >= 1500) & (q["f"] <= 14000)
    a1.plot(q["f"][m]/1e3, 10*np.log10(S[m]/PREF**2), color=q["col"],
            ls=q["lsty"], lw=q["lw"], zorder=5, label=q["name"])
    a1.plot([q["fpk"]/1e3], [q["lev"]], "v", ms=7, mfc=q["col"], mec="k",
            mew=0.7, zorder=7)
a1.axvline(SCREECH_EXP/1e3, color="0.35", lw=1.0, ls=":", zorder=3)
a1.text(SCREECH_EXP/1e3*0.985, 101, "Panda 5.4 kHz", fontsize=8.0,
        color="0.30", rotation=90, ha="right", va="bottom")
a1.set_xlim(1.5, 14)
a1.set_ylim(100, 122)
a1.set_xlabel("$f$ [kHz]")
a1.set_ylabel(r"PSD  [dB re 20 $\mu$Pa$^2$/Hz]")
a1.set_title(r"(a)  ring $o2$, $x=2$, $r\simeq1.0$", fontsize=9.8, pad=4)
a1.legend(loc="lower left", fontsize=8.2, handlelength=2.2)

for tag, *_ in CASES:
    q = Q[tag]
    a2.plot([q["fpow"]/1e3], [q["fpk"]/1e3], "o", ms=11, mfc=q["col"],
            mec="k", mew=0.8, zorder=6)
    a2.annotate(q["name"].replace("Skew-symmetric", "skew").replace(" + ", "+"),
                (q["fpow"]/1e3, q["fpk"]/1e3), textcoords="offset points",
                xytext=(0, 13), fontsize=8.0, ha="center")
lo, hi = 5.0, 6.6
a2.plot([lo, hi], [lo, hi], "-", color="0.55", lw=1.2, zorder=3,
        label="measured = Powell")
fp_exp = 0.7*U_J/(1.160*D*(1+0.7*U_J/A_INF))
a2.plot([fp_exp/1e3], [SCREECH_EXP/1e3], "*", ms=17, mfc="w", mec="k",
        mew=1.0, zorder=7)
a2.annotate("experiment", (fp_exp/1e3, SCREECH_EXP/1e3),
            textcoords="offset points", xytext=(0, -18), fontsize=8.4,
            ha="center")
a2.set_xlim(lo, hi)
a2.set_ylim(lo, hi)
a2.set_xlabel(r"Powell from each $L_s$  [kHz]")
a2.set_ylabel(r"measured peak  [kHz]")
a2.set_title("(b)  does the peak follow the shock cells?", fontsize=9.8, pad=4)
a2.legend(loc="upper left", fontsize=8.2, handlelength=1.8)

w = 0.26
for k, (tag, *_) in enumerate(CASES):
    q = Q[tag]
    a3.bar(np.arange(5) + (k-1)*w, [100*q["share"][m] for m in range(5)],
           width=w, color=q["col"], edgecolor="k", linewidth=0.5,
           alpha=1.0 if q["lsty"] == "-" else 0.55, zorder=5,
           label=q["name"])
a3.set_xticks(range(5))
a3.set_xticklabels(["0", "1", "2", "3", "4"])
a3.set_xlabel(r"azimuthal mode $|m|$")
a3.set_ylabel("share of band energy [%]")
a3.set_title(r"(c)  who carries the peak  (peak $\pm$ 700 Hz)",
             fontsize=9.8, pad=4)
a3.legend(loc="upper right", fontsize=8.0)
for ax in AX:
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)

fig.savefig(HERE/"peak_audit.png", dpi=250)
fig.savefig(HERE/"peak_audit.pdf")
print("\nwrote peak_audit.png/pdf")
