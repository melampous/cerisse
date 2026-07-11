#!/usr/bin/env python3
"""
cvp_rigorous.py -- Rigorous, publication-grade far-field-BC corruption metric
for the radiated sound of a co-rotating Lamb-Oseen vortex pair (CVP).

Replaces a crude single-snapshot p' difference with TWO drift-free, phase-aware,
time-windowed metrics evaluated in the acoustic far-field annulus and referenced
to the reflection-free WIDE truth:

  (A) Time-windowed per-point RMS field  -- broadband amplitude of the unsteady
      (acoustic) pressure, drift removed by subtracting each point's OWN time mean
      over the analysis window.  Phase-insensitive by construction (an amplitude),
      so instantaneous phase misalignment between runs cannot create spurious error.

  (B) Quadrupole Fourier amplitude |P_hat| at the acoustic fundamental
      f_ac = Omega/pi = 0.12732 (= 2*Omega/(2*pi); the m=2 quadrupole radiates at
      TWICE the orbital frequency).  Single-bin DFT, Hann-windowed in time so a
      non-integer-period window does not leak the DC/drift mode into the tone bin.
      Cleanest pure-acoustic measure (rejects near-field hydrodynamics at other f).

The single corruption number per BC is the area-weighted annulus-L2 ratio of the
(narrow - wide) field to the wide field, in percent.  A noise floor from a
split-half window on the wide run, a per-BC mean-pressure drift, a 4-panel |P_hat|
figure, and a directivity curve |P_hat|(theta) on r=9 are also produced.

Usage:
    python3 cvp_rigorous.py  <case_dir>

    <case_dir> must contain:  plot_wide_char/   plot_narrow_foe/
                              plot_narrow_char/ plot_narrow_lodi/
    each holding yt-loadable plotfiles  plt?????  with fields
    Density, pressure, x_velocity, y_velocity.

Outputs (written into <case_dir>):
    cvp_pHat_quadrupole.png   -- 4-panel |P_hat| amplitude field
    cvp_directivity.png       -- |P_hat|(theta) on r=9, all 4 BCs overlaid
    cvp_rigorous_summary.json -- machine-readable table
    + a clean stdout table.

Verified non-dim constants (from prob.h):
    gamma=1.4, rho0=1, c0=1, p0=1/gamma=0.7142857,
    Omega=Gamma/(pi*b^2)=0.40, T_rot=2*pi/Omega=15.70796,
    f_ac=Omega/pi=0.1273240, T_ac=lambda=1/f_ac=7.853982 (since c0=1).
"""

import os
import sys
import json
import glob
import math

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import yt

yt.set_log_level(50)  # quiet

# numpy renamed trapz -> trapezoid in 2.0+; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# ----------------------------------------------------------------------------
# VERIFIED CONSTANTS (non-dimensional)
# ----------------------------------------------------------------------------
GAMMA = 1.4
RHO0  = 1.0
C0    = 1.0
P0    = 1.0 / GAMMA                 # 0.7142857  quiescent background pressure
OMEGA = 0.40                        # orbital angular velocity
T_ROT = 2.0 * math.pi / OMEGA       # 15.70796
F_AC  = OMEGA / math.pi             # 0.1273240   acoustic fundamental (= 2*Omega/2pi)
T_AC  = 1.0 / F_AC                  # 7.853982   = lambda (c0=1)
W_AC  = 2.0 * math.pi * F_AC        # 0.80       = 2*Omega

# Sampling window (the common +/-12 comparison box, inscribed circle r=12)
L_WIN = 12.0
N_GRID = 480                        # dx = 2*L/N = 24/480 = 0.05

# Far-field annulus: R_IN >= 1.15*lambda (=9.03) keeps the WHOLE annulus in the
# radiation zone (quadrupole near-field ~1/r^2,1/r^3 terms are negligible),
# R_OUT=10.5 leaves >=1.5 inside the r=12 narrow boundary.
R_IN  = 9.0
R_OUT = 10.5

# Analysis window upper bound (stop_time of the AWS runs)
T_HI_TARGET = 55.0
# We choose t_lo so [t_lo,55] spans an INTEGER number of T_ac with >=2 periods.
# Preference order: 5,4,3 whole acoustic periods.  Whichever fits the common
# frame coverage of all 4 runs is used; if it cannot be made integer-period a
# Hann time window is applied so the DFT/mean are still leakage-clean.
MIN_PERIODS = 3                     # >=3 T_ac requested (>=2 hard minimum)

# Directivity sampling circle
R_DIR = 9.0
N_THETA = 180

# Runs:  label -> plot subdirectory name
RUNS = [
    ("wide", "plot_wide_char"),
    ("foe",  "plot_narrow_foe"),
    ("char", "plot_narrow_char"),
    ("lodi", "plot_narrow_lodi"),
]
PRETTY = {
    "wide": "wide_char (TRUTH)",
    "foe":  "narrow_foe",
    "char": "narrow_char",
    "lodi": "narrow_lodi (relax=0.25)",
}


# ----------------------------------------------------------------------------
# Frame inventory
# ----------------------------------------------------------------------------
def list_frames(plot_dir):
    """Return sorted list of (path, time) for yt-loadable plotfiles in plot_dir.

    Robust to a few different plotfile naming conventions and to unreadable
    frames (skipped with a warning).
    """
    if not os.path.isdir(plot_dir):
        return []
    candidates = []
    for pat in ("plt?????", "plt*", "plot?????", "*plt*"):
        candidates += glob.glob(os.path.join(plot_dir, pat))
    # de-dup, keep directories (AMReX plotfiles are directories) or Header files
    seen = set()
    plts = []
    for c in sorted(candidates):
        if c in seen:
            continue
        seen.add(c)
        if os.path.isdir(c) or os.path.isfile(os.path.join(c, "Header")):
            plts.append(c)
    out = []
    for p in plts:
        try:
            ds = yt.load(p)
            t = float(ds.current_time)
            out.append((p, t))
        except Exception as e:  # noqa: BLE001
            sys.stderr.write(f"  [warn] cannot load {p}: {e}\n")
            continue
    out.sort(key=lambda z: z[1])
    return out


# ----------------------------------------------------------------------------
# Field sampling onto the common +/-12 grid
# ----------------------------------------------------------------------------
def sample_pressure(ds_path):
    """Sample 'pressure' on the common +/-12 uniform grid via arbitrary_grid.

    Returns a (N_GRID, N_GRID) float64 array, or None on failure.
    """
    try:
        ds = yt.load(ds_path)
        ds.force_periodicity()  # required: fill the edge band of the sampling box
        left = ds.arr([-L_WIN, -L_WIN, 0.0], "code_length")
        right = ds.arr([L_WIN, L_WIN, 1.0], "code_length")
        g = ds.arbitrary_grid(left, right, dims=[N_GRID, N_GRID, 1])
        # field-tuple access, robust to yt naming
        try:
            p = np.asarray(g["boxlib", "pressure"])[:, :, 0]
        except Exception:  # noqa: BLE001
            p = np.asarray(g["pressure"])[:, :, 0]
        return np.ascontiguousarray(p, dtype=np.float64)
    except Exception as e:  # noqa: BLE001
        sys.stderr.write(f"  [warn] sample failed {ds_path}: {e}\n")
        return None


def build_geometry():
    """Cell-centered x,y axes and derived R, THETA on the common grid."""
    # cell centers of the +/-12 box at N_GRID samples
    dx = (2.0 * L_WIN) / N_GRID
    xc = -L_WIN + (np.arange(N_GRID) + 0.5) * dx
    yc = xc.copy()
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    R = np.hypot(X, Y)
    TH = np.arctan2(Y, X)
    return xc, yc, X, Y, R, TH, dx


# ----------------------------------------------------------------------------
# Window selection: integer T_ac if possible, else Hann
# ----------------------------------------------------------------------------
def choose_window(times_by_run):
    """Choose a common [t_lo, t_hi] integer-T_ac window covering all runs.

    times_by_run: dict label -> np.array of frame times (sorted).
    Returns (t_lo, t_hi, n_periods, integer_window: bool).
    """
    # common upper time = min over runs of the max available time, capped at 55
    t_hi = min(T_HI_TARGET, min(t[-1] for t in times_by_run.values()))
    # common lower bound on availability = max over runs of the min time
    t_min_common = max(t[0] for t in times_by_run.values())
    # earliest physically meaningful start: 1 rotation (transient cleared)
    t_lo_floor = max(T_ROT, t_min_common)

    # try the largest integer number of whole T_ac that fits [t_lo_floor, t_hi]
    span = t_hi - t_lo_floor
    n_max = int(math.floor(span / T_AC + 1e-9))

    def feasible(n):
        t_lo = t_hi - n * T_AC
        if t_lo < t_lo_floor - 1e-9:
            return None
        for t in times_by_run.values():
            k = int(np.count_nonzero((t >= t_lo - 1e-6) & (t <= t_hi + 1e-6)))
            if k < 6:
                return None
        return t_lo

    feas = [n for n in range(n_max, MIN_PERIODS - 1, -1) if feasible(n) is not None]
    if feas:
        # PREFER the largest EVEN n>=4 (lets us form a disjoint integer-period
        # split-half noise floor); else the largest feasible n.
        even = [n for n in feas if n >= 4 and n % 2 == 0]
        n = even[0] if even else feas[0]
        return feasible(n), t_hi, n, True

    # fallback: non-integer window from t_lo_floor (Hann will be applied)
    n_frac = span / T_AC
    return t_lo_floor, t_hi, n_frac, False


def frames_in_window(times, t_lo, t_hi):
    """Boolean index of frames inside [t_lo, t_hi] (small eps slack)."""
    return (times >= t_lo - 1e-6) & (times <= t_hi + 1e-6)


# ----------------------------------------------------------------------------
# Metric core
# ----------------------------------------------------------------------------
def trapz_time_mean(P, t):
    """Per-point trapezoidal time mean over [t0,t1]. P shape (K, N, N)."""
    t0, t1 = t[0], t[-1]
    denom = max(t1 - t0, 1e-30)
    return _trapz(P, t, axis=0) / denom


def rms_field(P, t):
    """Drift-free per-point temporal RMS (trapezoidal). P shape (K,N,N).

    Removes per-point DC + LINEAR drift (least-squares fit of [1, t-tbar]),
    NOT just the time-mean.  A constant-mean removal leaves a ramping mean
    (the un-anchored BC drift, e.g. LODI p_relax small) in the residual and
    inflates the RMS spuriously; detrending kills it so the RMS is the true
    unsteady (acoustic + unsteady near-field) amplitude only.
    """
    K = t.shape[0]
    tb = t.mean()
    G = np.empty((K, 2), dtype=np.float64)
    G[:, 0] = 1.0
    G[:, 1] = t - tb
    Pp = P.reshape(K, -1)
    coef, *_ = np.linalg.lstsq(G, Pp, rcond=None)     # (2, N*N): DC + slope
    res = (Pp - G @ coef).reshape(P.shape)            # detrended residual
    t0, t1 = t[0], t[-1]
    denom = max(t1 - t0, 1e-30)
    msq = _trapz(res * res, t, axis=0) / denom
    return np.sqrt(np.maximum(msq, 0.0))


def hann_weights(t, t_lo, t_hi):
    """Hann window evaluated at the actual (non-uniform) frame times."""
    xi = (t - t_lo) / max(t_hi - t_lo, 1e-30)
    xi = np.clip(xi, 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(2.0 * math.pi * xi))


def quad_fourier_amp(P, t, t_lo, t_hi, integer_window):
    """Complex quadrupole phasor field P_hat at f_ac via a per-point
    LEAST-SQUARES single-tone fit with an explicit DC + linear-drift basis.

    This is the verified-robust estimator: on non-uniform (CFL-driven) snapshot
    times with DC offset + mean drift, a naive single-bin DFT recovers |P_hat|
    to ~1% but gets the COMPLEX value wrong by >100% (the DC/drift mode is not
    orthogonal to the f_ac kernel off the uniform integer-period grid).  Fitting
    p'(t) = c0 + c1*cos(w_ac t) + c2*sin(w_ac t) + c3*(t - tbar) makes the tone
    EXACTLY orthogonal to mean + linear drift for ANY sampling, so neither the
    DC offset nor a ramping mean (foextrap's un-anchored drift) can pollute the
    tone.  Returns complex (N,N) P_hat = c1 + 1j*c2  with the physical
    convention p'(t) ~ Re{P_hat * exp(+i w_ac t)} = c1 cos - c2 sin; |P_hat| is
    phase-invariant (used for the amplitude metric) and arg(P_hat) carries the
    BC-induced phase shift (used for the complex metric).

    `integer_window` only governs the cross-check DFT; the LS fit is exact
    regardless, so the argument is accepted for interface symmetry.
    """
    K = t.shape[0]
    tbar = t.mean()
    # design matrix (K,4): [1, cos, sin, drift]
    G = np.empty((K, 4), dtype=np.float64)
    G[:, 0] = 1.0
    G[:, 1] = np.cos(W_AC * t)
    G[:, 2] = np.sin(W_AC * t)
    G[:, 3] = t - tbar
    # subtract base pressure to keep magnitudes O(1e-3) (the DC column absorbs
    # any residual; this is cosmetic for conditioning, not for drift removal)
    Pp = (P - P0).reshape(K, -1)                    # (K, N*N)
    coef, *_ = np.linalg.lstsq(G, Pp, rcond=None)   # (4, N*N)
    a1 = coef[1].reshape(P.shape[1], P.shape[2])
    b1 = coef[2].reshape(P.shape[1], P.shape[2])
    return a1 + 1j * b1


def quad_fourier_amp_dft(P, t, t_lo, t_hi, integer_window):
    """Hann-windowed single-bin DFT estimate of |P_hat| -- MAGNITUDE-ONLY
    cross-check on the LS fit (do not use its phase on non-uniform t)."""
    pbar = trapz_time_mean(P, t)
    pp = P - pbar[None, :, :]
    w = np.ones_like(t) if integer_window else hann_weights(t, t_lo, t_hi)
    cg = np.mean(w)
    if cg <= 0:
        cg = 1.0
    kern = np.exp(-1j * W_AC * t)
    wk = w * kern
    num = _trapz(pp * wk[:, None, None], t, axis=0)
    span = max(t[-1] - t[0], 1e-30)
    return (2.0 / span) * num / cg


def annulus_l2_ratio(field_narrow, field_wide, mask):
    """Area-weighted (uniform grid) annulus-L2 of (narrow-wide)/wide, percent.

    Works for real or complex fields; for complex it is the full complex L2.
    """
    diff = field_narrow[mask] - field_wide[mask]
    num = np.sqrt(np.sum(np.abs(diff) ** 2))
    den = np.sqrt(np.sum(np.abs(field_wide[mask]) ** 2))
    if den <= 0:
        return float("nan")
    return 100.0 * num / den


# ----------------------------------------------------------------------------
# Directivity
# ----------------------------------------------------------------------------
def directivity_curve(P_hat, xc, yc):
    """|P_hat|(theta) sampled on the r=R_DIR circle by bilinear interpolation
    of the real and imaginary parts separately. Returns (theta, |P_hat|)."""
    th = np.linspace(0.0, 2.0 * math.pi, N_THETA, endpoint=False)
    xq = R_DIR * np.cos(th)
    yq = R_DIR * np.sin(th)

    def bilinear(F, xq, yq):
        # map query coords to fractional grid indices on cell-centered axes
        dx = xc[1] - xc[0]
        fi = (xq - xc[0]) / dx
        fj = (yq - yc[0]) / dx
        i0 = np.clip(np.floor(fi).astype(int), 0, N_GRID - 2)
        j0 = np.clip(np.floor(fj).astype(int), 0, N_GRID - 2)
        ti = fi - i0
        tj = fj - j0
        v00 = F[i0, j0]
        v10 = F[i0 + 1, j0]
        v01 = F[i0, j0 + 1]
        v11 = F[i0 + 1, j0 + 1]
        return (v00 * (1 - ti) * (1 - tj) + v10 * ti * (1 - tj)
                + v01 * (1 - ti) * tj + v11 * ti * tj)

    ar = bilinear(P_hat.real, xq, yq)
    ai = bilinear(P_hat.imag, xq, yq)
    return th, np.hypot(ar, ai)


# ----------------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        sys.stderr.write("usage: python3 cvp_rigorous.py <case_dir>\n")
        sys.exit(2)
    case_dir = os.path.abspath(sys.argv[1])
    if not os.path.isdir(case_dir):
        sys.stderr.write(f"error: not a directory: {case_dir}\n")
        sys.exit(2)

    print("=" * 78)
    print("CVP rigorous far-field-BC acoustic corruption metric")
    print("=" * 78)
    print(f"case_dir = {case_dir}")
    print(f"constants: gamma={GAMMA}, p0={P0:.7f}, Omega={OMEGA}, "
          f"T_rot={T_ROT:.5f}, f_ac={F_AC:.6f}, T_ac=lambda={T_AC:.5f}")
    print(f"grid: +/-{L_WIN} window, N={N_GRID} (dx={2*L_WIN/N_GRID:.4f}); "
          f"annulus r in [{R_IN},{R_OUT}]")
    print()

    # --- 1. inventory frames per run -----------------------------------------
    frames = {}
    times_by_run = {}
    for label, sub in RUNS:
        pdir = os.path.join(case_dir, sub)
        fl = list_frames(pdir)
        if not fl:
            sys.stderr.write(f"  [warn] no frames for {label} in {pdir}\n")
            continue
        frames[label] = fl
        times_by_run[label] = np.array([t for _, t in fl])
        print(f"  {label:5s}: {len(fl):3d} frames, "
              f"t in [{fl[0][1]:.3f}, {fl[-1][1]:.3f}]  ({sub})")
    if "wide" not in frames:
        sys.stderr.write("error: wide (truth) run missing; cannot reference.\n")
        sys.exit(1)
    print()

    # --- 2. choose common window ---------------------------------------------
    t_lo, t_hi, n_per, integer_window = choose_window(times_by_run)
    if integer_window:
        print(f"window = [{t_lo:.4f}, {t_hi:.4f}]  = {int(round(n_per))} x T_ac "
              f"(integer-period, no Hann)")
    else:
        print(f"window = [{t_lo:.4f}, {t_hi:.4f}]  = {n_per:.3f} x T_ac "
              f"(NON-integer -> Hann time window applied)")
    print()

    # --- 3. sample fields, build P stacks ------------------------------------
    xc, yc, X, Y, R, TH, dx = build_geometry()
    mask = (R >= R_IN) & (R <= R_OUT)
    n_cells = int(np.count_nonzero(mask))
    print(f"annulus cells in mask = {n_cells}")
    print()

    rms_fields = {}
    phat_fields = {}
    drift = {}
    Pstacks = {}      # keep wide stack for split-half noise floor
    times_used = {}
    xcheck_methodfloor = float("nan")   # LS-vs-DFT |P_hat| disagreement (method noise)

    for label, _ in RUNS:
        if label not in frames:
            continue
        sel = frames_in_window(times_by_run[label], t_lo, t_hi)
        chosen = [frames[label][i] for i in np.where(sel)[0]]
        if len(chosen) < 6:
            sys.stderr.write(f"  [warn] {label}: only {len(chosen)} frames in "
                             f"window; skipping.\n")
            continue
        Pk = []
        tk = []
        for path, t in chosen:
            p = sample_pressure(path)
            if p is None:
                continue
            if not np.all(np.isfinite(p)):
                sys.stderr.write(f"  [warn] non-finite pressure in {path}; "
                                 f"skipping frame.\n")
                continue
            Pk.append(p)
            tk.append(t)
        if len(Pk) < 6:
            sys.stderr.write(f"  [warn] {label}: <6 valid frames; skipping.\n")
            continue
        P = np.stack(Pk, axis=0)            # (K,N,N) absolute pressure
        tarr = np.array(tk)
        times_used[label] = tarr

        # sanity: quiescent far corners ~ p0
        corner = P[:, :8, :8].mean()
        print(f"  {label:5s}: K={len(tk):3d} frames used; "
              f"far-corner <p>={corner:.5f} (p0={P0:.5f})")

        # metric A: drift-free RMS field
        rms_fields[label] = rms_field(P, tarr)
        # metric B: quadrupole Fourier amplitude (complex phasor, LS estimator)
        phat_fields[label] = quad_fourier_amp(P, tarr, t_lo, t_hi, integer_window)
        # mean-p drift over the window in the annulus (annulus-mean of p-p0)
        pbar_t = np.array([(P[k][mask]).mean() for k in range(P.shape[0])]) - P0
        drift[label] = float(pbar_t[-1] - pbar_t[0])
        if label == "wide":
            Pstacks["wide"] = (P, tarr)
            # DFT cross-check (magnitude only): LS and Hann-DFT |P_hat| should
            # agree to a few % in the far-field annulus; large divergence flags
            # leakage / a windowing problem.
            ph_dft = quad_fourier_amp_dft(P, tarr, t_lo, t_hi, integer_window)
            xcheck = annulus_l2_ratio(np.abs(ph_dft),
                                      np.abs(phat_fields[label]), mask)
            xcheck_methodfloor = xcheck
            print(f"         [xcheck] wide LS-vs-Hann-DFT |P_hat| "
                  f"annulus disagreement = {xcheck:.2f}%  (method-noise floor)")
    print()

    # --- 4. noise floor (DISJOINT integer-period split-half on wide) ----------
    # A split-half floor is only meaningful if EACH half is itself a whole
    # number of acoustic periods (>=2), i.e. the analysis window spans >=4 whole
    # T_ac.  On a shorter (e.g. 3-period) window the half-records alias the
    # single tone and the "floor" balloons to a meaningless >100% -- so we GATE
    # it and otherwise fall back to the LS-vs-DFT method-noise proxy (xcheck).
    floor_rms = float("nan")
    floor_fou = float("nan")
    floor_kind = "none"
    n_whole = int(round(n_per)) if integer_window else int(math.floor(n_per))
    if "wide" in Pstacks and integer_window and n_whole >= 4 and n_whole % 2 == 0:
        Pw, tw = Pstacks["wide"]
        half = n_whole // 2
        t_split = t_lo + half * T_AC                 # boundary between two halves
        i1 = tw <= t_split + 1e-6                     # first  `half` periods
        i2 = tw >  t_split - 1e-6                     # second `half` periods (disjoint)
        if np.count_nonzero(i1) >= 6 and np.count_nonzero(i2) >= 6:
            rms_h1 = rms_field(Pw[i1], tw[i1])
            rms_h2 = rms_field(Pw[i2], tw[i2])
            floor_rms = annulus_l2_ratio(rms_h1, rms_h2, mask) / 100.0
            ph_h1 = quad_fourier_amp(Pw[i1], tw[i1], tw[i1][0], tw[i1][-1], True)
            ph_h2 = quad_fourier_amp(Pw[i2], tw[i2], tw[i2][0], tw[i2][-1], True)
            floor_fou = annulus_l2_ratio(np.abs(ph_h1), np.abs(ph_h2), mask) / 100.0
            floor_kind = f"split-half ({half}+{half} T_ac, disjoint)"
    if not np.isfinite(floor_fou):
        # window too short for a clean split-half -> use method-noise proxy
        if np.isfinite(xcheck_methodfloor):
            floor_fou = xcheck_methodfloor / 100.0
            floor_kind = "LS-vs-DFT method-noise proxy (window <4 T_ac)"

    # --- 5. corruption numbers vs wide ---------------------------------------
    rms_wide = rms_fields.get("wide")
    phat_wide = phat_fields.get("wide")

    results = {}
    for label, _ in RUNS:
        if label not in rms_fields:
            continue
        if label == "wide":
            results[label] = dict(rms_corruption_pct=0.0,
                                  fourier_amp_corruption_pct=0.0,
                                  fourier_cplx_corruption_pct=0.0,
                                  mean_drift=drift.get(label, float("nan")))
            continue
        rms_err = annulus_l2_ratio(rms_fields[label], rms_wide, mask)
        # amplitude-only Fourier (robust, phase-free)
        fou_amp = annulus_l2_ratio(np.abs(phat_fields[label]),
                                   np.abs(phat_wide), mask)
        # complex Fourier (amplitude + phase, stricter)
        fou_cplx = annulus_l2_ratio(phat_fields[label], phat_wide, mask)
        results[label] = dict(rms_corruption_pct=rms_err,
                              fourier_amp_corruption_pct=fou_amp,
                              fourier_cplx_corruption_pct=fou_cplx,
                              mean_drift=drift.get(label, float("nan")))

    # floor-subtracted (in quadrature) clean values
    for label in results:
        if label == "wide":
            continue
        r = results[label]
        if np.isfinite(floor_rms):
            r["rms_clean_pct"] = 100.0 * math.sqrt(
                max((r["rms_corruption_pct"] / 100.0) ** 2 - floor_rms ** 2, 0.0))
        if np.isfinite(floor_fou):
            r["fourier_amp_clean_pct"] = 100.0 * math.sqrt(
                max((r["fourier_amp_corruption_pct"] / 100.0) ** 2
                    - floor_fou ** 2, 0.0))

    # --- 6. print the table ---------------------------------------------------
    print("=" * 78)
    print("CORRUPTION TABLE  (annulus-L2 vs WIDE truth, percent)")
    print("=" * 78)
    hdr = (f"{'BC':<24}{'fourier_amp_%':>14}{'fourier_cplx_%':>16}"
           f"{'rms_%':>10}{'mean_drift':>14}")
    print(hdr)
    print("-" * len(hdr))
    for label, _ in RUNS:
        if label not in results:
            continue
        r = results[label]
        print(f"{PRETTY[label]:<24}"
              f"{r['fourier_amp_corruption_pct']:>14.3f}"
              f"{r['fourier_cplx_corruption_pct']:>16.3f}"
              f"{r['rms_corruption_pct']:>10.3f}"
              f"{r['mean_drift']:>14.3e}")
    print("-" * len(hdr))
    if np.isfinite(floor_fou):
        rms_str = f"{100*floor_rms:>10.3f}" if np.isfinite(floor_rms) else f"{'n/a':>10}"
        print(f"{'NOISE FLOOR':<24}"
              f"{100*floor_fou:>14.3f}{'':>16}{rms_str}")
        print(f"  floor kind: {floor_kind}")
    print()
    print("Notes:")
    print("  * fourier_amp = |P_hat| amplitude error at f_ac (phase-free, robust).")
    print("  * fourier_cplx = full complex P_hat error (amplitude+phase, stricter).")
    print("  * rms = broadband per-point temporal-RMS error (drift removed).")
    print("  * mean_drift = annulus-mean (p-p0) change over the window;")
    print("    the real BC discriminator at near-normal incidence (anchoring).")
    print("  * E_clean = sqrt(E^2 - floor^2) reported in JSON.")
    print()

    # --- 7. figure 1: 4-panel |P_hat| ----------------------------------------
    order = [l for l, _ in RUNS if l in phat_fields]
    vmax = 0.0
    for l in order:
        amp = np.abs(phat_fields[l])
        vmax = max(vmax, float(amp[mask].max()) if np.any(mask) else amp.max())
    if vmax <= 0:
        vmax = 1e-6

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.ravel()
    extent = [-L_WIN, L_WIN, -L_WIN, L_WIN]
    for ax, l in zip(axes, ["wide", "foe", "char", "lodi"]):
        if l not in phat_fields:
            ax.set_visible(False)
            continue
        amp = np.abs(phat_fields[l]).T  # transpose: imshow wants [row=y, col=x]
        im = ax.imshow(amp, origin="lower", extent=extent, vmin=0, vmax=vmax,
                       cmap="inferno", aspect="equal")
        # annulus + boundary overlays
        for rr, ls in ((R_IN, "--"), (R_OUT, "--")):
            ax.add_patch(plt.Circle((0, 0), rr, fill=False, color="cyan",
                                    ls=ls, lw=1.0))
        ax.add_patch(plt.Circle((0, 0), T_AC, fill=False, color="white",
                                ls=":", lw=1.0))  # lambda circle
        ax.add_patch(plt.Rectangle((-L_WIN, -L_WIN), 2 * L_WIN, 2 * L_WIN,
                                   fill=False, color="lime", lw=0.8))
        r = results.get(l, {})
        sub = (f"{r.get('fourier_amp_corruption_pct', 0.0):.2f}%"
               if l != "wide" else "TRUTH")
        ax.set_title(f"{PRETTY[l]}\n|P_hat| corruption = {sub}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(r"Quadrupole acoustic amplitude $|\hat P|$ at $f_{ac}=%.4f$  "
                 r"(annulus dashed, $\lambda$ dotted)" % F_AC)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out1 = os.path.join(case_dir, "cvp_pHat_quadrupole.png")
    fig.savefig(out1, dpi=140)
    plt.close(fig)
    print(f"wrote {out1}")

    # --- 8. figure 2: directivity |P_hat|(theta) on r=9 ----------------------
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    colors = {"wide": "k", "foe": "tab:blue", "char": "tab:green",
              "lodi": "tab:red"}
    dir_curves = {}
    for l in ["wide", "foe", "char", "lodi"]:
        if l not in phat_fields:
            continue
        th, amp = directivity_curve(phat_fields[l], xc, yc)
        dir_curves[l] = (th.tolist(), amp.tolist())
        lw = 2.5 if l == "wide" else 1.6
        ax2.plot(np.degrees(th), amp, color=colors[l], lw=lw,
                 label=PRETTY[l])
    ax2.set_xlabel(r"$\theta$ (deg)")
    ax2.set_ylabel(r"$|\hat P|$ on $r=%.0f$" % R_DIR)
    ax2.set_title(r"Directivity of the radiated quadrupole "
                  r"(4-lobe $\cos 2\theta$ pattern)")
    ax2.set_xlim(0, 360)
    ax2.set_xticks(range(0, 361, 45))
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best", fontsize=9)
    fig2.tight_layout()
    out2 = os.path.join(case_dir, "cvp_directivity.png")
    fig2.savefig(out2, dpi=140)
    plt.close(fig2)
    print(f"wrote {out2}")

    # --- 9. JSON summary ------------------------------------------------------
    summary = dict(
        case_dir=case_dir,
        constants=dict(gamma=GAMMA, p0=P0, Omega=OMEGA, T_rot=T_ROT,
                       f_ac=F_AC, T_ac=T_AC),
        window=dict(t_lo=t_lo, t_hi=t_hi, n_periods=n_per,
                    integer_window=bool(integer_window)),
        annulus=dict(r_in=R_IN, r_out=R_OUT, n_cells=n_cells),
        grid=dict(L=L_WIN, N=N_GRID, dx=2 * L_WIN / N_GRID),
        noise_floor=dict(rms_pct=100 * floor_rms if np.isfinite(floor_rms)
                         else None,
                         fourier_amp_pct=100 * floor_fou
                         if np.isfinite(floor_fou) else None,
                         kind=floor_kind),
        results={l: results[l] for l in results},
        directivity={l: dict(theta=dir_curves[l][0], amp=dir_curves[l][1])
                     for l in dir_curves},
        frames_used={l: int(len(times_used[l])) for l in times_used},
    )
    out3 = os.path.join(case_dir, "cvp_rigorous_summary.json")
    with open(out3, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {out3}")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
