#!/usr/bin/env python3
"""Reduce each centreline density profile to a few physical scalars by
fitting a damped, linearly-chirped cosine (sine x decay) on top of a
slowly decaying mean.

  rho(x)/rho_j = trend(x) + A * exp(-(x-x1)/Ld) * cos(phi(x))
  trend(x)     = c0 + c1 * exp(-(x-x1)/Lm)
  phi(x)       = (2*pi/(Ls0*beta)) * ln(1 + beta*(x-x1)) + phi0
                 (instantaneous cell spacing Ls(x) = Ls0*(1+beta*(x-x1)))

Fitted over x in [x1, 8] D_e (the cell train; the first diamond is not
sinusoidal). Each curve collapses to the scalar fingerprint
(A, Ld, Ls0, beta, c0, Lm): first-cell shock strength, envelope decay
length, initial cell spacing, spacing growth rate, far-field mean and
mean-decay length. Curve-to-curve comparison then reduces to comparing
these six numbers."""
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RHO_J = 1.6413
CASES = [("baseline", "Baseline", "#aab0b7"),
         ("tophat", "Top hat", "#0b6ef5"),
         ("walltanh", "Wall tanh", "#2ca02c"),
         ("shift100", "Thin shifted tanh", "#ff6a00")]


def boxfilt(a, x, wD=0.05):
    w = max(3, int(round(wD / (x[1] - x[0]))))
    if w % 2 == 0:
        w += 1
    out = np.convolve(a, np.ones(w) / w, mode="same")
    out[:w // 2] = a[:w // 2]; out[-(w // 2):] = a[-(w // 2):]
    return out


def first_peak(x, y):
    exp = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = exp[np.nanargmin(y[exp])]
    comp = np.where((np.arange(x.size) > imin) & (x <= 1.7))[0]
    return x[comp[np.nanargmax(y[comp])]]


def model(x, x1, A, Ld, Ls0, beta, c0, c1, Lm, phi0):
    xr = x - x1
    trend = c0 + c1 * np.exp(-xr / Lm)
    phase = (2 * np.pi / (Ls0 * beta)) * np.log1p(beta * xr) + phi0
    return trend + A * np.exp(-xr / Ld) * np.cos(phase)


def fit_case(v):
    d = np.load(HERE / f"AXIS8_{v}.npz")
    x = d["x"]
    rho = boxfilt(d["DensityMEAN"] / RHO_J, x)
    x1 = first_peak(x, rho)
    m = (x >= x1) & (x <= 8.0)
    xf, yf = x[m], rho[m]
    # initial guesses from peaks
    pk, _ = find_peaks(yf, prominence=0.02,
                       distance=max(1, int(0.5 / (x[1] - x[0]))))
    xn = xf[pk]
    An = yf[pk]
    Ls0 = float(np.diff(xn)[0]) if len(xn) > 1 else 1.3
    beta = 0.02
    if len(xn) > 2:
        beta = float(np.clip((np.diff(xn)[-1] - np.diff(xn)[0])
                             / (Ls0 * (xn[-1] - xn[0])), 1e-3, 0.5))
    A0 = float(An[0] - yf.mean())
    Ld0 = 3.0
    p0 = [x1, A0, Ld0, Ls0, beta, 0.75, 0.3, 3.0, 0.0]
    lo = [x1 - 1e-6, 0.05, 0.5, 0.8, 1e-3, 0.5, -0.5, 0.5, -np.pi]
    hi = [x1 + 1e-6, 0.8, 30.0, 2.0, 0.5, 1.1, 1.0, 20.0, np.pi]
    try:
        popt, _ = curve_fit(model, xf, yf, p0=p0, bounds=(lo, hi),
                            maxfev=40000)
    except Exception as e:
        print("  fit failed for %s: %s" % (v, e)); popt = p0
    resid = yf - model(xf, *popt)
    rms = float(np.sqrt(np.mean(resid**2)))
    return x, rho, xf, yf, popt, rms


PARN = ["x1", "A", "Ld", "Ls0", "beta", "c0", "c1", "Lm", "phi0"]
rows = []
fig, axs = plt.subplots(2, 2, figsize=(8.4, 5.4), sharex=True, sharey=True)
plt.rcParams.update({"font.size": 9})
for ax, (v, lab, c) in zip(axs.flat, CASES):
    x, rho, xf, yf, p, rms = fit_case(v)
    rows.append((lab, p, rms))
    m = (x > 0.05) & (x <= 8.0)
    ax.plot(x[m], rho[m], color="0.4", lw=1.4, label="data")
    ax.plot(xf, model(xf, *p), color=c, lw=1.5, ls="--", label="fit")
    ax.set_axisbelow(True); ax.grid(True, ls="--", lw=0.5, color="0.88")
    ax.set_title(lab, fontsize=9)
    ax.text(0.03, 0.05,
            "$L_{s0}$=%.2f $\\beta$=%.3f\n$L_d$=%.2f $A$=%.2f\nrms=%.3f"
            % (p[3], p[4], p[2], p[1], rms),
            transform=ax.transAxes, fontsize=7.5, va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))
    ax.set_ylim(0.3, 1.55)
axs[0, 0].legend(fontsize=7.5, frameon=False, loc="upper right")
for ax in axs[1]:
    ax.set_xlabel("$x/D_e$")
for ax in axs[:, 0]:
    ax.set_ylabel(r"$\langle\rho\rangle/\rho_j$")
fig.tight_layout()
fig.savefig(HERE / "fit_shockcell_model.png", dpi=250)
fig.savefig(HERE / "fit_shockcell_model.pdf")
print("wrote fit_shockcell_model.png/pdf\n")
print("%-18s %6s %6s %6s %6s %6s %6s %6s"
      % ("case", "A", "Ld", "Ls0", "beta", "c0", "Lm", "rms"))
for lab, p, rms in rows:
    print("%-18s %6.3f %6.2f %6.3f %6.3f %6.3f %6.2f %6.4f"
          % (lab, p[1], p[2], p[3], p[4], p[5], p[7], rms))
