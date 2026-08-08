#!/usr/bin/env python3
"""Mesh-sensitivity metrics with block-based statistical uncertainty.

Data: blockax{2d,3d}_L*_rho.npz — cumulative centreline MEAN/SQR profiles
of pressure, axial velocity and density at every statistics snapshot.
Block statistics are obtained by tau-weighted differencing
(tau_k = t_k - t_s); blocks shorter than 0.2 ms are merged.

Figures:
  paperfig_shockpos_vs_res : x_fc, x1 and S1-S3 vs D_e/dx_min
  paperfig_rmsheat_2x2     : rho_rms/rho_j and u_rms/U_j heatmaps, x vs level
  paperfig_rmsint_vs_res   : A_p^[0.05,2], A_p^[2,6] and xi_c vs D_e/dx_min

Metric definitions: x_fc is the maximum positive gradient of the unsmoothed
time-averaged centreline density between the first expansion minimum and the
first downstream principal density maximum. The locations x_n are successive
principal maxima of the mean centreline density after smoothing over 0.05 D
(prominence 0.02 rho_j, minimum separation 0.5 D, parabolic refinement).
S_n = x_{n+1}-x_n is evaluated between successive density maxima;
A_p^[a,b] is the root mean square of p_rms/p_inf over [a,b]D; xi_c is the
q_rms^2-weighted axial centre over [0.05,8]D."""
import numpy as np
trapz_f = getattr(np, 'trapezoid', None) or np.trapz
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

plt.rcParams.update({
    "font.size": 9.4, "axes.linewidth": 0.65,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254; P_AMB = 99780.0; RHO_J = 1.6413; UJ = 414.2
TBLK_MIN = 0.2e-3

CASES = [("2d", t) for t in "0123456"] + [("3d", t) for t in "12345"]

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
    d = (y[i-1] - 2*y[i] + y[i+1])
    if d == 0:
        return x[i]
    return x[i] + 0.5*(y[i-1] - y[i+1])/d * (x[1]-x[0])

def density_xfc(mode, tag):
    """Density-based closure location from the full mean."""
    if mode == "2d":
        fname = {"0": "AXT_L0lo.npz", "1": "AXT_L1lo.npz",
                 "2": "AXT_L2lo.npz"}.get(tag, "AXT_L%s.npz" % tag)
    else:
        fname = "AX3T_L%s.npz" % tag
    d = np.load(HERE / fname)
    x = (d["z"] if "z" in d else d["x"]) / D
    rho = (d["stat_DensityMEAN"] if "stat_DensityMEAN" in d
           else d["ax_DensityMEAN"]) / RHO_J
    drho = np.gradient(rho, x)
    expansion = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = expansion[np.nanargmin(rho[expansion])]
    compression = np.where((np.arange(x.size) > imin) & (x <= 1.65))[0]
    imax = compression[np.nanargmax(rho[compression])]
    ii = np.arange(imin, imax + 1)
    xfc = refine(x, drho, ii[np.nanargmax(drho[ii])])
    return xfc

def metrics(x, pM, pS, uM, uS, rhoM, rhoS, npk=4):
    """All shock/rms metrics from one set of axis profiles. The shock-cell
    positions x_fc and x_n are density-based; the r.m.s. measures use the
    pressure and axial-velocity moments."""
    rho = rhoM / RHO_J
    # density-based x_fc: maximum positive gradient of the unsmoothed mean
    # density across the first recompression interval
    drho = np.gradient(rho, x)
    exp = np.where((x >= 0.35) & (x <= 1.0))[0]
    imin = exp[np.nanargmin(rho[exp])]
    comp = np.where((np.arange(x.size) > imin) & (x <= 1.65))[0]
    imax = comp[np.nanargmax(rho[comp])]
    ii = np.arange(imin, imax + 1)
    xfc = refine(x, drho, ii[np.nanargmax(drho[ii])])
    # density-based x_n: successive principal maxima of the smoothed mean
    # density downstream of x_fc
    rs = boxfilt(rho, x, 0.05)
    dist = max(1, int(round(0.5 / (x[1] - x[0]))))
    pk, _ = find_peaks(np.where(x > xfc, rs, -np.inf),
                       prominence=0.02, distance=dist)
    xn = [refine(x, rs, i) for i in pk[:npk]]
    prms = np.sqrt(np.maximum(pS - pM**2, 0)) / P_AMB
    urms = np.sqrt(np.maximum(uS - uM**2, 0)) / UJ
    rhorms = np.sqrt(np.maximum(rhoS - rhoM**2, 0)) / RHO_J
    def segment(q, a, b):
        """Interpolate to common integration bounds for every hierarchy."""
        mm = (x > a) & (x < b)
        xx = np.concatenate(([a], x[mm], [b]))
        qq = np.interp(xx, x, q)
        return xx, qq
    def Ap(a, b):
        xx, qq = segment(prms, a, b)
        return np.sqrt(trapz_f(qq**2, xx) / (b-a))
    x8p, p8 = segment(prms, 0.05, 8.0)
    x8u, u8 = segment(urms, 0.05, 8.0)
    xcp = trapz_f(x8p*p8**2, x8p) / trapz_f(p8**2, x8p)
    xcu = trapz_f(x8u*u8**2, x8u) / trapz_f(u8**2, x8u)
    return dict(xfc=xfc, x1=xn[0] if xn else np.nan, xn=xn,
                Ap02=Ap(0.05, 2.0), Ap26=Ap(2.0, 6.0),
                xcp=xcp, xcu=xcu, prms=prms, urms=urms, rhorms=rhorms)

RES = {}
for mode, tag in CASES:
    d = np.load(HERE / ("blockax%s_L%s_rho.npz" % (mode, tag)))
    prof, times, ts, dx = d["prof"], d["times"], float(d["t_s"]), float(d["dx"])
    x = (np.arange(prof.shape[2]) + 0.5) * dx / D
    tau = times - ts
    # full window (x_fc and x_n now density-based inside metrics())
    full = metrics(x, *prof[-1])
    # block edges: greedy >= TBLK_MIN, always keep last snapshot
    edges = [-1]                       # -1 = the tau=0 zero anchor
    last_t = 0.0
    for k in range(len(times)):
        if tau[k] - last_t >= TBLK_MIN or k == len(times) - 1:
            edges.append(k); last_t = tau[k]
    blocks = []
    for a, b in zip(edges[:-1], edges[1:]):
        Ma = prof[a] if a >= 0 else 0.0
        ta = tau[a] if a >= 0 else 0.0
        Mb, tb = prof[b], tau[b]
        blocks.append(metrics(x, *((tb*Mb - ta*Ma) / (tb - ta))))
    def se(key):
        v = np.array([bl[key] for bl in blocks], dtype=float)
        v = v[np.isfinite(v)]
        return np.std(v, ddof=1)/np.sqrt(len(v)) if len(v) > 1 else np.nan
    RES[(mode, tag)] = dict(full=full, x=x, nominal_resolution=int(round(D/dx)),
                            nblk=len(blocks),
                            se={k: se(k) for k in ("xfc", "x1", "Ap02",
                                                   "Ap26", "xcp", "xcu")})
    print("%s L%s D/dx=%4d nblk=%d xfc=%.3f+-%.3f x1=%.3f+-%.3f "
          "Ap02=%.4f+-%.4f Ap26=%.4f+-%.4f xcp=%.2f+-%.2f xcu=%.2f+-%.2f"
          % (mode, tag, RES[(mode, tag)]["nominal_resolution"], len(blocks),
             full["xfc"], se("xfc"), full["x1"], se("x1"),
             full["Ap02"], se("Ap02"), full["Ap26"], se("Ap26"),
             full["xcp"], se("xcp"), full["xcu"], se("xcu")))

def fam(mode):
    return dict(marker="o", mfc=None, ls="-") if mode == "2d" \
        else dict(marker="s", mfc="none", ls="--")

def fam_shockpos(mode):
    return dict(marker="o", mfc="none", ls="--") if mode == "2d" \
        else dict(marker="s", mfc=None, ls="-")

def mode_label(mode):
    """Publication label for the two coordinate formulations."""
    return "RZ" if mode == "2d" else "Cartesian"

# ---- Figure 1: shock positions vs resolution -------------------------------
# Experimental reference: the mean of the 37 phase-conditioned M_j = 1.42
# centreline density profiles in Panda & Seasholtz (1999), EPAPS m142den.
# The reference locations and spacings are extracted from the unsmoothed
# cycle-mean density field with three-point parabolic refinement. The
# M_j = 1.42 survey has an axial sampling interval of 0.08 D, giving a nominal
# position-sampling range of +-0.04 D. The conservative spacing range
# is +-0.08 D. These are sampling-based ranges, not experimental
# confidence intervals.
EXP_XFC, EXP_X1 = 0.888, 1.199
EXP_S = (1.359, 1.186, 1.130)
EB_POS = 0.04
EB_S = 0.08
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
for mode in ("2d", "3d"):
    tags = "0123456" if mode == "2d" else "12345"
    Rv = [RES[(mode, t)]["nominal_resolution"] for t in tags]
    st = fam_shockpos(mode)
    for key, col, lab in (("xfc", "#1f4e9c", "$x_{fc}$"), ("x1", "#b3251e", "$x_1$")):
        y = [RES[(mode, t)]["full"][key] for t in tags]
        axs[0].plot(Rv, y, color=col, ms=4, lw=0.9,
                    label="%s %s" % (lab, mode_label(mode)), **st)
    for n, (ordinal, col) in enumerate((("First", "#1f4e9c"),
                                         ("Second", "#2e8b57"),
                                         ("Third", "#c8781e"))):
        y = [(RES[(mode, t)]["full"]["xn"][n+1] - RES[(mode, t)]["full"]["xn"][n])
             if len(RES[(mode, t)]["full"]["xn"]) > n+1 else np.nan
             for t in tags]
        # The 2D spacing estimates are retained only for illustration. Use a
        # common light grey and distinguish S1--S3 by their open markers.
        plot_col = "#b8b8b8" if mode == "2d" else col
        plot_style = dict(st)
        if mode == "2d":
            plot_style["marker"] = ("o", "s", "^")[n]
        axs[1].plot(Rv, y, color=plot_col, ms=4, lw=0.9,
                    label="$S_%d$ %s" % (n+1, mode_label(mode)),
                    **plot_style)
for y, col in zip((EXP_XFC, EXP_X1), ("#1f4e9c", "#b3251e")):
    axs[0].axhspan(y-EB_POS, y+EB_POS, color=col, alpha=0.10, lw=0,
                   zorder=0.5)
    axs[0].axhline(y, color=col, alpha=0.65, lw=0.6, ls="-.", zorder=0.7)
for y, col in zip(EXP_S, ("#1f4e9c", "#2e8b57", "#c8781e")):
    axs[1].axhspan(y-EB_S, y+EB_S, color=col, alpha=0.10, lw=0,
                   zorder=0.5)
    axs[1].axhline(y, color=col, alpha=0.65, lw=0.6, ls="-.", zorder=0.7)
for ax in axs:
    ax.set_xscale("log", base=2)
    ax.set_xticks([16, 32, 64, 128, 256, 512, 1024])
    ax.set_xticklabels(["16", "32", "64", "128", "256", "512", "1024"])
    ax.set_xlabel("$D_e/\\Delta x_{\\min}$", labelpad=2)
    ax.tick_params(length=2.5)
    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="both", color="0.88", alpha=0.70,
            linewidth=0.5, zorder=0)
axs[0].set_ylabel("Axial position, $x/D_e$", labelpad=2)
axs[1].set_ylabel("$S_n/D_e$", labelpad=2)
handles, labels = axs[0].get_legend_handles_labels()
legend_order = (0, 2, 1, 3)
axs[0].legend([handles[i] for i in legend_order],
              [labels[i] for i in legend_order],
              fontsize=7.2, frameon=False, ncol=2, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), handlelength=1.6,
              columnspacing=0.8, borderaxespad=0)
handles, labels = axs[1].get_legend_handles_labels()
legend_order = (0, 3, 1, 4, 2, 5)
axs[1].legend([handles[i] for i in legend_order],
              [labels[i] for i in legend_order],
              fontsize=6.9, frameon=False, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), handlelength=1.4,
              columnspacing=0.7, borderaxespad=0)
for k, ax in enumerate(axs):
    ax.text(0.02, 0.96, "(%s)" % chr(97+k), transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
fig.tight_layout(w_pad=1.4)
fig.savefig(HERE / "paperfig_shockpos_vs_res.png", dpi=300)
fig.savefig(HERE / "paperfig_shockpos_vs_res.pdf")
print("wrote paperfig_shockpos_vs_res")

# ---- Figure 2: rms heatmaps 2x2 --------------------------------------------
XH = 8.0
fig, axs = plt.subplots(2, 2, figsize=(7.2, 4.4), sharex=True)

# Build all arrays first so that the two coordinate formulations use the same
# colour limit for a given variable.
HEAT = {}
for mode in ("2d", "3d"):
    tags = "0123456" if mode == "2d" else "12345"
    for qk in ("rho", "u"):
        NXC = 800
        xg = np.linspace(0, XH, NXC+1)
        xc = 0.5*(xg[:-1]+xg[1:])
        Z = np.full((len(tags), NXC), np.nan)
        for i, t in enumerate(tags):
            r = RES[(mode, t)]
            q = r["full"]["rhorms" if qk == "rho" else "urms"]
            qs = boxfilt(q, r["x"], 0.05)
            Z[i] = np.interp(xc, r["x"], qs)
        HEAT[(mode, qk)] = (tags, xg, Z)
VMAX = {
    qk: float(np.ceil(np.nanpercentile(
        np.concatenate([HEAT[(mode, qk)][2].ravel()
                        for mode in ("2d", "3d")]), 99.5) * 50) / 50)
    for qk in ("rho", "u")
}

for row, mode in enumerate(("2d", "3d")):
    for colk, (qk, cmap) in enumerate((("rho", "magma"), ("u", "viridis"))):
        ax = axs[row, colk]
        tags, xg, Z = HEAT[(mode, qk)]
        im = ax.pcolormesh(xg, np.arange(len(tags)+1), Z, cmap=cmap,
                           vmin=0, vmax=VMAX[qk], rasterized=True)
        for i, t in enumerate(tags):
            result = RES[(mode, t)]["full"]
            positions = [result["xfc"]] + result["xn"][:3]
            colours = ("#1f4e9c", "#b3251e", "#2e8b57", "#c8781e")
            for xpos, colour in zip(positions, colours):
                ax.scatter(xpos, i + 0.5, marker="o", s=7, color=colour,
                           edgecolor="white", linewidth=0.25, zorder=5)
        ax.set_yticks(np.arange(len(tags)) + 0.5)
        ax.set_yticklabels([str(RES[(mode, t)]["nominal_resolution"])
                            for t in tags],
                           fontsize=7.5)
        ax.tick_params(length=2.5)
        cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.05)
        cb.set_label("$\\rho_{\\mathrm{rms}}/\\rho_j$" if qk == "rho"
                     else "$u_{x,\\mathrm{rms}}/U_j$", fontsize=8, labelpad=2)
        cb.ax.tick_params(length=2, labelsize=7.5)
        tagl = "(%s) %s" % (chr(97 + 2*row + colk), mode_label(mode))
        ax.text(0.985, 0.94, tagl, transform=ax.transAxes, ha="right",
                va="top", fontsize=8.5, color="w",
                bbox=dict(facecolor="black", edgecolor="none", alpha=0.35,
                          pad=1.0))
        if colk == 0:
            ax.set_ylabel("$D_e/\\Delta x_{\\min}$", labelpad=2)
for ax in axs[1]:
    ax.set_xlabel("$x/D_e$", labelpad=2)
fig.legend(
    [Line2D([0], [0], marker="o", ms=3.5, color="none",
            markerfacecolor=colour, markeredgecolor="white",
            markeredgewidth=0.25)
     for colour in ("#1f4e9c", "#b3251e", "#2e8b57", "#c8781e")],
    ["$x_{fc}$", "$x_1$", "$x_2$", "$x_3$"],
    loc="upper center", ncol=4,
    frameon=False, fontsize=8, handlelength=1.0, columnspacing=0.9,
    bbox_to_anchor=(0.5, 1.005),
)
fig.tight_layout(rect=(0, 0, 1, 0.96), w_pad=1.0, h_pad=0.6)
fig.savefig(HERE / "paperfig_rmsheat_2x2.png", dpi=300)
fig.savefig(HERE / "paperfig_rmsheat_2x2.pdf")
print("wrote paperfig_rmsheat_2x2")

# ---- Figure 3: integral fluctuation metrics vs resolution ------------------
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7))
for mode in ("2d", "3d"):
    tags = "0123456" if mode == "2d" else "12345"
    Rv = [RES[(mode, t)]["nominal_resolution"] for t in tags]
    st = fam(mode)
    for key, col, lab in (("Ap02", "#7b1fa2", "$A_p^{[0.05,2]}$"),
                          ("Ap26", "#e08214", "$A_p^{[2,6]}$")):
        y = [RES[(mode, t)]["full"][key] for t in tags]
        axs[0].plot(Rv, y, color=col, ms=4, lw=0.9,
                    label="%s %s" % (lab, mode_label(mode)), **st)
    for key, col, lab in (("xcp", "#7b1fa2", "$\\xi_c^{(p)}$"),
                          ("xcu", "#2e8b57", "$\\xi_c^{(u_x)}$")):
        y = [RES[(mode, t)]["full"][key] for t in tags]
        axs[1].plot(Rv, y, color=col, ms=4, lw=0.9,
                    label="%s %s" % (lab, mode_label(mode)), **st)
for ax in axs:
    ax.set_xscale("log", base=2)
    ax.set_xticks([16, 32, 64, 128, 256, 512, 1024])
    ax.set_xticklabels(["16", "32", "64", "128", "256", "512", "1024"])
    ax.set_xlabel("$D_e/\\Delta x_{\\min}$", labelpad=2)
    ax.tick_params(length=2.5)
axs[0].set_ylabel("$A_p$", labelpad=2)
axs[1].set_ylabel("$\\xi_c$", labelpad=2)
axs[0].legend(fontsize=6.5, frameon=False, loc="lower right", ncol=2,
              handlelength=1.6, columnspacing=0.8)
axs[1].legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2,
              handlelength=1.6, columnspacing=0.8)
for k, ax in enumerate(axs):
    ax.text(0.02, 0.96, "(%s)" % chr(97+k), transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
fig.tight_layout(w_pad=1.4)
fig.savefig(HERE / "paperfig_rmsint_vs_res.png", dpi=300)
fig.savefig(HERE / "paperfig_rmsint_vs_res.pdf")
print("wrote paperfig_rmsint_vs_res")
