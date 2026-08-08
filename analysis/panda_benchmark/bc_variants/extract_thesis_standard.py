#!/usr/bin/env python3
"""Re-extract every case using the thesis definition verbatim.

Section 05 (sec:uej-verification) of the chapter states:

  - the centreline profile is built from the finest available AMR level at
    each axial location;
  - x_fc is the maximum positive density gradient between the first expansion
    minimum and the first downstream principal density maximum, taken on the
    RAW mean density;
  - x_n are the successive principal maxima of the numerical mean density
    after a moving box filter of nominal width 0.05 D_e, each with prominence
    at least 0.02 in <rho>/rho_j and separated by at least 0.5 D_e;
  - every selected gradient or density maximum is refined by a local
    quadratic fit through the sample and its two neighbours;
  - S_n = x_{n+1} - x_n.

The same definition is applied to calculation and experiment. The sampling
uncertainty reported alongside is not part of the thesis definition; it is
obtained by applying that same definition to disjoint sub-windows of each
case's own statistics record.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parent))
from design_2d_sweep import scalars, D

HERE = Path(__file__).resolve().parent
RHO_J = 1.6413
T0 = 4.0e-4                    # statistics phase start for the new 2D runs

SPEC = [("A0_tophat",    0.0, 3.0, "shifted"),
        ("C1_wall200", 200.0, 3.0, "wall"),
        ("A1_a100",    100.0, 3.0, "shifted"),
        ("B1_wall680", 679.5, 3.0, "wall"),
        ("A2_a170",    170.0, 3.0, "shifted"),
        ("S4_a170_s4", 170.0, 4.0, "shifted"),
        ("A3_a254",    254.0, 3.0, "shifted"),
        ("S5_a170_s5", 170.0, 5.0, "shifted"),
        ("S6_a170_s6", 170.0, 6.0, "shifted"),
        ("S7_a170_s7", 170.0, 7.0, "shifted")]


def boxfilt(a, x, wD=0.05):
    """Moving box filter of nominal width wD in D_e, as in section 05."""
    n = max(3, int(round(wD / (x[1] - x[0]))))
    if n % 2 == 0:
        n += 1
    o = np.convolve(a, np.ones(n) / n, mode="same")
    o[:n//2] = a[:n//2]
    o[-(n//2):] = a[-(n//2):]
    return o


def refine(x, y, i):
    """Local quadratic fit through the sample and its two neighbours."""
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + 0.5*(y[i-1]-y[i+1])/den*(x[1]-x[0])


def thesis_metrics(x, rho, nmax=5):
    """x_fc from the raw profile; x_n from the box-filtered profile."""
    e = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = e[np.nanargmin(rho[e])]
    c = np.where((np.arange(x.size) > imin) & (x <= 1.7))[0]
    imax = c[np.nanargmax(rho[c])]
    g = np.gradient(rho, x)                      # raw, unsmoothed
    seg = np.arange(imin, imax + 1)
    xfc = refine(x, g, seg[np.nanargmax(g[seg])])

    sm = boxfilt(rho, x)
    dist = max(1, int(round(0.5 / (x[1] - x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, sm, -np.inf),
                       prominence=0.02, distance=dist)
    xn = np.array([refine(x, sm, i) for i in pk[:nmax]])
    return xfc, xn


def load(nm):
    f = HERE / "sweep_out" / ("XFC_%s.npz" % nm)
    if not f.exists():
        return None
    d = np.load(f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda s: int("".join(c for c in s if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    t = np.array([float(d[k+"_t"][0]) for k in ks])
    M = np.array([d[k+"_rho"] for k in ks]) / RHO_J
    return (np.arange(M.shape[1])+0.5)*dx/D, t, M


print("thesis definition: raw-gradient maximum for x_fc, 0.05 D_e box filter "
      "for x_n,\nprominence 0.02, separation 0.5 D_e, quadratic refinement\n")
hdr = ("%-13s %6s %7s | %7s %7s %7s %7s %7s | %7s %7s %7s"
       % ("case", "a[um]", "d*[um]", "x_fc", "x_1", "x_2", "x_3", "x_4",
          "S_1", "S_2", "S_3"))
print(hdr); print("-" * len(hdr))
rows = []
for nm, a, s, form in SPEC:
    got = load(nm)
    if got is None:
        print("%-13s %6.0f %7s | MISSING" % (nm, a, "-")); continue
    x, t, M = got
    xfc, xn = thesis_metrics(x, M[-1])
    S = np.diff(xn)          # S_n = x_{n+1} - x_n, over the x_n series only
    sc = scalars(a*1e-6, form, s)
    # sampling uncertainty: same definition on disjoint sub-windows
    W = t - T0
    sub = [(W[i+1]*M[i+1]-W[i]*M[i])/(t[i+1]-t[i])
           for i in range(len(t)-1) if t[i+1]-t[i] > 1.5e-4]
    vv = np.array([thesis_metrics(x, q)[0] for q in sub])
    sig = np.nanstd(vv, ddof=1)/np.sqrt(len(vv))
    rows.append((nm, a, s, form, sc, xfc, sig, xn, S))
    print("%-13s %6.0f %7.1f | %7.4f %7.4f %7.4f %7.4f %7.4f | %7.4f %7.4f %7.4f"
          % (nm, a, sc["dstar_c"]*1e6, xfc,
             *(list(xn[:4]) + [np.nan]*(4-len(xn[:4]))),
             *(list(S[:3]) + [np.nan]*(3-len(S[:3])))))

print("\nsampling uncertainty on x_fc (thesis definition, disjoint sub-windows)")
for nm, a, s, form, sc, xfc, sig, xn, S in rows:
    print("   %-13s x_fc = %.4f +- %.4f" % (nm, xfc, sig))

lad = [r for r in rows if r[1] == 170.0]
lad.sort(key=lambda r: r[2])
if len(lad) >= 3:
    print("\nshifted tanh series, a = 170 um fixed "
          "(delta_omega = 340 um, theta = 75 um frozen)")
    print("   %4s %9s %9s %9s" % ("s", "d*[um]", "x_fc", "+-"))
    for r in lad:
        print("   %4.0f %9.1f %9.4f %9.4f"
              % (r[2], r[4]["dstar_c"]*1e6, r[5], r[6]))
    X = np.array([r[4]["dstar_c"]*1e6 for r in lad])
    Y = np.array([r[5] for r in lad]); Sg = np.array([r[6] for r in lad])
    m, c = np.polyfit(X, Y, 1, w=1/Sg)
    r2 = 1 - ((Y-(m*X+c))**2).sum()/((Y-Y.mean())**2).sum()
    print("   weighted slope = %+.4f per 1000 um   R2 = %.3f" % (m*1000, r2))
    print("   3D four-case correlation = -0.0766 per 1000 um")
np.savez(HERE/"sweep_out"/"THESIS_STD.npz",
         **{r[0]: np.array([r[5], r[6], r[4]["dstar_c"]]) for r in rows})
