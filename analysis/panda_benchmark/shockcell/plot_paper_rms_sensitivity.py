"""Paper figures: grid sensitivity of second-order statistics (centreline rho_rms).

Replaces the mean+rms 4-variable overlays for the grid-sensitivity section.
Design: (a) smoothed rms profiles — fine-branch min-max envelope + finest as
reference + coarse levels as individual labeled lines; (b) deviation from the
finest solution (fine branch only); (c) near-field peak amplitude and location
versus cell size (the scalar convergence view).
2D family: L0(D/16), L1(D/32), L2..L6 (D/64..D/1024, envelope).
3D family: L1(D/16)..L3(D/64) lines, L4-L5 (D/128, D/256) envelope.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/shockcell"
D, RHOJ = 0.0254, 1.6413
SMOOTH = 0.25          # boxcar width in D_e units
XMAX = 10.0

plt.rcParams.update({"font.size": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True})


def smooth(y, dx):
    w = max(3, int(SMOOTH / dx) | 1)
    return np.convolve(y, np.ones(w) / w, mode="same")


def load_rms(fname, zkey, mkey, skey):
    d = np.load(f"{BASE}/{fname}")
    z = d[zkey] / D
    rms = np.sqrt(np.clip(d[skey] - d[mkey]**2, 0.0, None)) / RHOJ
    m = z <= XMAX + 1.0
    return z[m], rms[m]


FAMS = {
    "2d": {
        "coarse": [("AXT_L0lo.npz", 16, r"$\ell_{\max}=0$ ($D_e/16$)"),
                   ("AXT_L1lo.npz", 32, r"$\ell_{\max}=1$ ($D_e/32$)")],
        "fine":   [("AXT_L2lo.npz", 64), ("AXT_L3.npz", 128), ("AXT_L4.npz", 256),
                   ("AXT_L5.npz", 512), ("AXT_L6.npz", 1024)],
        "keys": ("z", "stat_DensityMEAN", "stat_DensitySQR"),
        "scale": lambda d: (d["stat_DensityMEAN"], d["stat_DensitySQR"]),
        "env_label": r"envelope $\ell_{\max}=2$–$6$ ($D_e/64$–$D_e/1024$)",
        "ref_label": r"$\ell_{\max}=6$ ($D_e/1024$, reference)",
        "title": "axisymmetric (RZ) family",
    },
    "3d": {
        "coarse": [("AX3T_L1.npz", 16, r"$\ell_{\max}=1$ ($D_e/16$)"),
                   ("AX3T_L2.npz", 32, r"$\ell_{\max}=2$ ($D_e/32$)"),
                   ("AX3T_L3.npz", 64, r"$\ell_{\max}=3$ ($D_e/64$)")],
        "fine":   [("AX3T_L4.npz", 128), ("AX3T_L5.npz", 256)],
        "keys": ("x", "ax_DensityMEAN", "ax_DensitySQR"),
        "env_label": r"envelope $\ell_{\max}=4$–$5$ ($D_e/128$–$D_e/256$)",
        "ref_label": r"$\ell_{\max}=5$ ($D_e/256$, reference)",
        "title": "three-dimensional family",
    },
}

CO_COL = ["#e6a35c", "#d95f02", "#a63603"]          # warm = under-resolved
CO_LS = [":", "-.", "--"]
ENV_FC, ENV_EC = "#c6dbef", "#2171b5"                # Blues family
REF_C = "#08306b"

for fam, cfg in FAMS.items():
    zk, mk, sk = cfg["keys"]

    # ---- load everything on its own grid
    coarse = [(load_rms(f, zk, mk, sk), n, lab) for f, n, lab in cfg["coarse"]]
    fine = [(load_rms(f, zk, mk, sk), n) for f, n in cfg["fine"]]

    # ---- common grid for envelope/deviation (coarsest fine-branch spacing)
    zc = np.arange(0.0, XMAX, 1.0 / 64)
    fine_s = []
    for (z, r), n in fine:
        rs = smooth(r, z[1] - z[0])
        fine_s.append((np.interp(zc, z, rs), n))
    F = np.array([a for a, _ in fine_s])
    env_lo, env_hi = F.min(0), F.max(0)
    ref = F[-1]

    # ---- scalar metrics: near-field peak (x<=5) amplitude + location
    mets = []
    for (z, r), n, _ in coarse:
        rs = smooth(r, z[1] - z[0]); m = z <= 5
        mets.append((n, rs[m].max(), z[m][np.argmax(rs[m])]))
    for (z, r), n in fine:
        rs = smooth(r, z[1] - z[0]); m = z <= 5
        mets.append((n, rs[m].max(), z[m][np.argmax(rs[m])]))
    mets.sort(key=lambda a: a[0])

    # ---- figure
    fig = plt.figure(figsize=(9.6, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0],
                          hspace=0.30, wspace=0.27,
                          left=0.085, right=0.955, top=0.94, bottom=0.085)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    # (a) envelope + reference + coarse lines
    axA.fill_between(zc, env_lo, env_hi, color=ENV_FC, edgecolor=ENV_EC,
                     lw=0.6, alpha=0.9, label=cfg["env_label"])
    axA.plot(zc, ref, color=REF_C, lw=1.8, label=cfg["ref_label"])
    for ((z, r), n, lab), c, ls in zip(coarse, CO_COL, CO_LS):
        axA.plot(z, smooth(r, z[1] - z[0]), color=c, ls=ls, lw=1.6, label=lab)
    axA.set_xlim(0, XMAX); axA.set_ylim(0, None)
    axA.set_xlabel(r"$x/D_e$"); axA.set_ylabel(r"$\rho_{\mathrm{rms}}/\rho_j$")
    axA.legend(fontsize=8.5, loc="upper left", ncol=2, framealpha=0.95)
    axA.text(0.015, 0.04, "(a)", transform=axA.transAxes, fontsize=12)
    axA.set_title(cfg["title"], fontsize=11)

    # (b) fine-branch deviation from the reference
    for (a, n), shade in zip(fine_s[:-1], np.linspace(0.45, 0.8, len(fine_s) - 1)):
        axB.plot(zc, a - ref, color=plt.cm.Blues(shade), lw=1.2,
                 label=rf"$D_e/{n}$")
    axB.axhline(0, color="0.3", lw=0.8)
    axB.set_xlim(0, XMAX)
    axB.set_xlabel(r"$x/D_e$")
    axB.set_ylabel(r"$\rho_{\mathrm{rms}}/\rho_j - (\rho_{\mathrm{rms}}/\rho_j)_{\mathrm{ref}}$")
    axB.legend(fontsize=8, loc="upper right", ncol=2, framealpha=0.95)
    axB.text(0.03, 0.05, "(b)  fine branch $-$ reference", transform=axB.transAxes, fontsize=10)

    # (c) peak amplitude + location vs cell size
    ns = np.array([m[0] for m in mets], dtype=float)
    amp = [m[1] for m in mets]; loc = [m[2] for m in mets]
    axC.plot(ns, amp, "o-", color=REF_C, lw=1.3, ms=6, label="peak amplitude")
    axC.set_xscale("log", base=2)
    axC.set_xticks(ns); axC.set_xticklabels([rf"$D_e/{int(n)}$" for n in ns],
                                            rotation=45, fontsize=8.5)
    axC.set_ylabel(r"$\max_{x\leq 5D_e}\ \rho_{\mathrm{rms}}/\rho_j$", color=REF_C)
    axC.tick_params(axis="y", colors=REF_C)
    axC2 = axC.twinx()
    axC2.plot(ns, loc, "s--", color="#d95f02", lw=1.2, ms=6, mfc="none",
              label="peak location")
    axC2.set_ylabel(r"$x_{\mathrm{peak}}/D_e$", color="#d95f02", labelpad=4)
    axC2.tick_params(axis="y", colors="#d95f02")
    axC.invert_xaxis()
    axC.text(0.06, 0.88, "(c)", transform=axC.transAxes, fontsize=12)

    for ext in ("png", "pdf"):
        fig.savefig(f"{BASE}/paperfig_rms_sensitivity_{fam}.{ext}", dpi=250)
    plt.close(fig)
    print(fam, "saved;", " ".join(f"{int(n)}:{a:.3f}@{l:.2f}" for n, a, l in mets))
