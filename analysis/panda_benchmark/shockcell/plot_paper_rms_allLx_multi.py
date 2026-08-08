"""Grid-sensitivity rms figures, flat single-axes variant, three quantities:
rho_rms/rho_j, p_rms/p_inf, M_rms = u'_rms/sqrt(gamma R <T>) (fluctuating
Mach number: axial-velocity rms over local mean sound speed). All levels
overlaid (L0..L6 for 2D, L1..L5 for 3D), 0.25 D_e smoothing, in-axes table
over x/D_e > 5 with each level's two largest peaks. Labels L0..L6 only.
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/shockcell"
D, RHOJ, PAMB = 0.0254, 1.6413, 99780.0
GAM, RGAS = 1.4, 8.31446 / 0.02896
SMOOTH = 0.25
XMAX = 10.0

plt.rcParams.update({"font.size": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True})

def rms_of(d, mkey, skey):
    return np.sqrt(np.clip(d[skey] - d[mkey]**2, 0.0, None))

def profile(fname, fam, qty):
    d = np.load(f"{BASE}/{fname}")
    if fam == "2d":
        z = d["z"] / D; pre = "stat_"; uax = "y_velocity"   # axial = y in 2D RZ
    else:
        z = d["x"] / D; pre = "ax_"; uax = "x_velocity"     # axial = x in 3D
    if qty == "rho":
        y = rms_of(d, pre+"DensityMEAN", pre+"DensitySQR") / RHOJ
    elif qty == "p":
        y = rms_of(d, pre+"pressureMEAN", pre+"pressureSQR") / PAMB
    else:
        urms = rms_of(d, pre+uax+"MEAN", pre+uax+"SQR")
        c = np.sqrt(GAM * RGAS * d[pre+"temperatureMEAN"])
        y = urms / c
    m = z <= XMAX + 0.5
    return z[m], y[m]

def smooth(y, dx):
    w = max(3, int(SMOOTH/dx) | 1)
    return np.convolve(y, np.ones(w)/w, mode="same")

def top2(z, s):
    pk, _ = find_peaks(s, prominence=0.008)
    if len(pk) == 0: return [(np.nan,np.nan)]*2
    order = np.argsort(s[pk])[::-1][:2]
    out = [(z[pk[i]], s[pk[i]]) for i in order]
    while len(out) < 2: out.append((np.nan, np.nan))
    return out

FAMS = {
    "2d": {"files": [("AXT_L0lo.npz","L0"),("AXT_L1lo.npz","L1"),("AXT_L2lo.npz","L2"),
                     ("AXT_L3.npz","L3"),("AXT_L4.npz","L4"),("AXT_L5.npz","L5"),
                     ("AXT_L6.npz","L6")],
           "coarse_n": 2, "title": "axisymmetric (RZ) family"},
    "3d": {"files": [("AX3T_L1.npz","L1"),("AX3T_L2.npz","L2"),("AX3T_L3.npz","L3"),
                     ("AX3T_L4.npz","L4"),("AX3T_L5.npz","L5")],
           "coarse_n": 1, "title": "three-dimensional family"},
}
QTYS = {"rho": r"$\rho_{\mathrm{rms}}/\rho_j$",
        "p":   r"$p_{\mathrm{rms}}/p_\infty$",
        "M":   r"$M_{\mathrm{rms}}$"}
QSYM = {"rho": r"\rho_{\mathrm{rms}}/\rho_j",
        "p":   r"p_{\mathrm{rms}}/p_\infty",
        "M":   r"M_{\mathrm{rms}}"}
WARM = ["#d95f02", "#a63603"]; WARM_LS = [":", "-."]

for fam, cfg in FAMS.items():
    for qty, ylab in QTYS.items():
        fig, ax = plt.subplots(figsize=(11.2, 4.1))
        rows, ymax = [], 0.0
        n_co = cfg["coarse_n"]; n_fi = len(cfg["files"]) - n_co
        blues = [plt.cm.Blues(v) for v in np.linspace(0.45, 1.0, n_fi)]
        for i, (f, lab) in enumerate(cfg["files"]):
            z, r = profile(f, fam, qty)
            s = smooth(r, z[1]-z[0])
            mplt = z <= XMAX - SMOOTH
            z, s = z[mplt], s[mplt]
            if i < n_co:
                ax.plot(z, s, color=WARM[i], ls=WARM_LS[i], lw=1.7, label=lab)
            else:
                j = i - n_co
                ax.plot(z, s, color=blues[j], lw=1.4 if j < n_fi-1 else 2.1, label=lab)
            (xa,va),(xb,vb) = top2(z, s)
            rows.append([lab, f"{xa:.2f}", f"{va:.3f}", f"{xb:.2f}", f"{vb:.3f}"])
            ymax = max(ymax, s.max())
        ax.set_xlim(0, XMAX); ax.set_ylim(0, 1.42*ymax)
        ax.set_xlabel(r"$x/D_e$"); ax.set_ylabel(ylab)
        ax.legend(fontsize=8.5, loc="upper left", ncol=1, framealpha=0.95,
                  handlelength=1.6, borderaxespad=0.4, labelspacing=0.3)
        ax.set_title(cfg["title"], fontsize=11)
        cols = ["level", r"$x^{(1)}/D_e$", rf"$({QSYM[qty]})^{{(1)}}$",
                r"$x^{(2)}/D_e$", rf"$({QSYM[qty]})^{{(2)}}$"]
        tab = ax.table(cellText=rows, colLabels=cols, cellLoc="center",
                       colLoc="center", bbox=[0.505, 0.50, 0.487, 0.48])
        tab.auto_set_font_size(False); tab.set_fontsize(8.2)
        for (ri,ci), cell in tab.get_celld().items():
            cell.set_edgecolor("0.55"); cell.set_linewidth(0.5)
            cell.set_facecolor("white")
            if ri == 0:
                cell.set_text_props(weight="bold"); cell.set_facecolor("#eef3fa")
        tab.set_zorder(6)
        fig.tight_layout()
        for ext in ("png","pdf"):
            fig.savefig(f"{BASE}/paperfig_rms_allLx_{fam}_{qty}.{ext}", dpi=250)
        plt.close(fig)
        print(fam, qty, "L-peaks:", [(r[0], r[1], r[2]) for r in rows])
