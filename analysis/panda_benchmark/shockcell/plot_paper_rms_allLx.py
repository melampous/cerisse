"""Grid-sensitivity rms figure, single-axes variant (user spec 2026-07-20):
all levels overlaid (L0..L6 for 2D, L1..L5 for 3D), 0.25 D_e smoothed,
with an in-axes table over the x/D_e > 5 region listing, per level, the two
largest peaks of rho_rms (location and value). Level labels are L0..L6 only.
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/shockcell"
D, RHOJ = 0.0254, 1.6413
SMOOTH = 0.25
XMAX = 10.0

plt.rcParams.update({"font.size": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True})


def load_rms(fname, zkey, mkey, skey):
    d = np.load(f"{BASE}/{fname}")
    z = d[zkey] / D
    rms = np.sqrt(np.clip(d[skey] - d[mkey]**2, 0.0, None)) / RHOJ
    m = z <= XMAX + 0.5
    return z[m], rms[m]


def smooth(y, dx):
    w = max(3, int(SMOOTH / dx) | 1)
    return np.convolve(y, np.ones(w) / w, mode="same")


def top2_peaks(z, s):
    pk, props = find_peaks(s, prominence=0.008)
    if len(pk) == 0:
        return [(np.nan, np.nan), (np.nan, np.nan)]
    order = np.argsort(s[pk])[::-1][:2]
    out = [(z[pk[i]], s[pk[i]]) for i in order]
    while len(out) < 2:
        out.append((np.nan, np.nan))
    return out


FAMS = {
    "2d": {
        "files": [("AXT_L0lo.npz", "L0"), ("AXT_L1lo.npz", "L1"), ("AXT_L2lo.npz", "L2"),
                  ("AXT_L3.npz", "L3"), ("AXT_L4.npz", "L4"), ("AXT_L5.npz", "L5"),
                  ("AXT_L6.npz", "L6")],
        "keys": ("z", "stat_DensityMEAN", "stat_DensitySQR"),
        "coarse_n": 2,                       # L0, L1 drawn warm
        "title": "axisymmetric (RZ) family",
        "ylim": 0.355,
        "tab_y": 0.985,
    },
    "3d": {
        "files": [("AX3T_L1.npz", "L1"), ("AX3T_L2.npz", "L2"), ("AX3T_L3.npz", "L3"),
                  ("AX3T_L4.npz", "L4"), ("AX3T_L5.npz", "L5")],
        "keys": ("x", "ax_DensityMEAN", "ax_DensitySQR"),
        "coarse_n": 1,                       # L1 warm
        "title": "three-dimensional family",
        "ylim": 0.20,
        "tab_y": 0.985,
    },
}

WARM = ["#d95f02", "#a63603"]
WARM_LS = [":", "-."]

for fam, cfg in FAMS.items():
    zk, mk, sk = cfg["keys"]
    fig, ax = plt.subplots(figsize=(11.2, 6.0))
    rows = []
    n_co = cfg["coarse_n"]
    n_fi = len(cfg["files"]) - n_co
    blues = [plt.cm.Blues(v) for v in np.linspace(0.45, 1.0, n_fi)]
    for i, (f, lab) in enumerate(cfg["files"]):
        z, r = load_rms(f, zk, mk, sk)
        s = smooth(r, z[1] - z[0])
        mplt = z <= XMAX - SMOOTH
        z, s = z[mplt], s[mplt]
        if i < n_co:
            ax.plot(z, s, color=WARM[i], ls=WARM_LS[i], lw=1.7, label=lab)
        else:
            j = i - n_co
            ax.plot(z, s, color=blues[j], lw=1.5 if j < n_fi - 1 else 2.2, label=lab)
        (xa, va), (xb, vb) = top2_peaks(z, s)
        rows.append([lab, f"{xa:.2f}", f"{va:.3f}", f"{xb:.2f}", f"{vb:.3f}"])

    ax.set_xlim(0, XMAX); ax.set_ylim(0, cfg["ylim"])
    ax.set_xlabel(r"$x/D_e$"); ax.set_ylabel(r"$\rho_{\mathrm{rms}}/\rho_j$")
    ax.legend(fontsize=9, loc="upper left", ncol=1, framealpha=0.95)
    ax.set_title(cfg["title"], fontsize=11)

    # in-axes table over the x/De > 5 region
    col_labels = ["level", r"$x^{(1)}/D_e$", r"$\rho^{(1)}_{\mathrm{rms}}/\rho_j$",
                  r"$x^{(2)}/D_e$", r"$\rho^{(2)}_{\mathrm{rms}}/\rho_j$"]
    tab = ax.table(cellText=rows, colLabels=col_labels,
                   cellLoc="center", colLoc="center",
                   bbox=[0.505, 0.545, 0.487, 0.435])
    tab.auto_set_font_size(False)
    tab.set_fontsize(8.6)
    for (ri, ci), cell in tab.get_celld().items():
        cell.set_edgecolor("0.55"); cell.set_linewidth(0.5)
        cell.set_facecolor("white")
        if ri == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#eef3fa")
    tab.set_zorder(6)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{BASE}/paperfig_rms_allLx_{fam}.{ext}", dpi=250)
    plt.close(fig)
    print(fam, "rows:", rows)
