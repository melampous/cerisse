#!/usr/bin/env python3
"""m142 shock-cell train analysis: CFD (time-mean, chk10083 stats) vs
Panda & Seasholtz 1999 fig4 (M=1.43 centerline density).

Outputs figures + a numbers report:
  1) centerline profiles (rho, p, Mach) + detected peaks + AMR level
     boundaries + experimental peaks
  2) E_n (cumulative position error), L_n (cell length), eps_n
  3) jet-scale diagnostics: shear-layer radius r_s(x), vorticity thickness
     delta_omega(x), M=1 core radius, centerline Mach, mean total-pressure loss
  4) Prandtl-Pack check with nominal D and sqrt(Cd)*D
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter

BASE = "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/shockcell"
D = 0.0254
GAM, RGAS = 1.4, 287.0
MJ = 1.42

d = np.load(f"{BASE}/m142_stats_axis_plane.npz")
ax_x = d["axis_x"] / D
A = d["axis"]                       # comps: [0]=u [1]=v [2]=w [3]=u2 [9->idx4]=p ...
CA = list(d["comps_axis"])          # [0,1,2,3,9,10,11,12,14]
def ac(c): return A[CA.index(c)]
u, p, T, rho = ac(0), ac(9), ac(10), ac(11)
mach = u / np.sqrt(GAM * RGAS * T)
lev_bounds = sorted(set(round(b / D, 3) for l, a, b in d["cover"] if b / D < 23.9))

exp = np.loadtxt(f"{BASE}/exp_m142_centerline.csv", delimiter=",", skiprows=1)
EXP_PEAKS = np.array([1.20, 2.60, 3.80, 4.90, 5.90, 6.60, 7.30])

# ---- peak detection: 3 quantities, mild smoothing, parabolic refinement ----
def detect(sig, x, kind="max", xmin=0.5, xmax=9.5, win=9):
    s = savgol_filter(sig, win, 3)
    if kind == "min":
        s = -s
    idx, _ = find_peaks(s, prominence=0.02 * (s.max() - s.min()))
    out = []
    for i in idx:
        if not (xmin <= x[i] <= xmax):
            continue
        i0, i1 = max(i - 3, 1), min(i + 4, len(x) - 1)
        c = np.polyfit(x[i0:i1], s[i0:i1], 2)
        xp = -c[1] / (2 * c[0]) if c[0] != 0 else x[i]
        out.append(xp if abs(xp - x[i]) < 0.2 else x[i])
    return np.array(out)

pk_rho = detect(rho, ax_x)
pk_p = detect(p, ax_x)
pk_M = detect(mach, ax_x, kind="min")

def consensus(plists, tol=0.25):
    ref = plists[0]
    cons = []
    for x0 in ref:
        grp = [x0]
        for pl in plists[1:]:
            if len(pl):
                j = np.argmin(np.abs(pl - x0))
                if abs(pl[j] - x0) < tol:
                    grp.append(pl[j])
        if len(grp) >= 2:                      # seen by >=2 methods
            cons.append((np.median(grp), np.ptp(grp), len(grp)))
    return cons

cons = consensus([pk_rho, pk_p, pk_M])
print("=== CFD consensus peaks (median [spread] nmethods) ===")
for i, (xm, sp, nm) in enumerate(cons):
    print(f"  raw {i+1}: x/D = {xm:.3f}  [{sp:.3f}]  ({nm}/3)")
all_cfd = np.array([c[0] for c in cons])
all_spread = np.array([c[1] for c in cons])

# pair each experimental peak with the NEAREST CFD peak (within 0.5D);
# CFD peaks left unpaired are substructure / beyond-exp-range, listed separately
cfd = np.full(len(EXP_PEAKS), np.nan)
spread = np.zeros(len(EXP_PEAKS))
used = set()
for i, xe in enumerate(EXP_PEAKS):
    j = int(np.argmin(np.abs(all_cfd - xe)))
    if abs(all_cfd[j] - xe) < 0.5 and j not in used:
        cfd[i] = all_cfd[j]; spread[i] = all_spread[j]; used.add(j)
extras = [x for k, x in enumerate(all_cfd) if k not in used]
print("unpaired CFD peaks (substructure / beyond exp range):",
      [round(x, 3) for x in extras])

n = np.arange(1, len(EXP_PEAKS) + 1)
E = cfd - EXP_PEAKS
Lc = np.diff(cfd)
Le = np.diff(EXP_PEAKS)
eps = (Lc - Le) / Le
EXP_ERR = np.array([0.05, 0.05, 0.05, 0.05, 0.08, 0.10, 0.10])  # quantization + falling prominence
print("=== paired E_n / L_n / eps_n ===")
for i in range(len(EXP_PEAKS)):
    ln = f"  n={i+1}: x_cfd={cfd[i]:.3f}  x_exp={EXP_PEAKS[i]:.2f}  E_n={E[i]:+.3f}D (exp +-{EXP_ERR[i]:.2f})"
    if i < len(Lc):
        ln += f"   L_n: cfd={Lc[i]:.3f} exp={Le[i]:.2f}  eps={eps[i]*100:+.1f}%"
    print(ln)
scale = np.polyfit(EXP_PEAKS[:4], cfd[:4], 1)[0]
print(f"linear-fit scale factor over cells 1-4: x_cfd ~ {scale:.4f} * x_exp"
      f"  -> implied Cd_eff = {scale**2:.3f}")

# ---- Prandtl-Pack ----
DJ = ((1 + 0.2 * MJ**2) / 1.2) ** 1.5 / np.sqrt(MJ)   # Dj/D fully expanded
LS_pp = 1.306 * np.sqrt(MJ**2 - 1) * DJ
for cd, tag in ((1.0, "nominal D"), (0.905, "Cd=0.905 (D/8 audit)"), (0.8789, "Cd=0.879 (fine audit)")):
    print(f"  Prandtl-Pack Ls1 [{tag}]: {LS_pp*np.sqrt(cd):.3f} D  "
          f"(CFD L1={Lc[0]:.3f}, exp L1={Le[0]:.2f})")

# ---- jet-scale diagnostics from mid-plane ----
P = d["plane"]                       # comps [0,1,2,9,10,11] at L4 res
CPL = list(d["comps_plane"])
NXP, NZP = P.shape[1], P.shape[2]
xp_ = (np.arange(NXP) + 0.5) * 24.0 / NXP
zp_ = -8.0 + (np.arange(NZP) + 0.5) * 16.0 / NZP
up = P[CPL.index(0)]
pp_ = P[CPL.index(9)]
Tp = P[CPL.index(10)]
Mp = np.abs(up) / np.sqrt(GAM * RGAS * Tp)
kc = NZP // 2
zpos = zp_[kc:]
rs, dom, rM1 = np.zeros(NXP), np.zeros(NXP), np.zeros(NXP)
for i in range(NXP):
    prof = up[i, kc:]
    g = np.abs(np.gradient(prof, zpos * D))
    j = np.argmax(g)
    rs[i] = zpos[j]
    du = prof.max() - 0.0
    dom[i] = du / max(g[j], 1e-9) / D
    sup = np.where(Mp[i, kc:] > 1.0)[0]
    rM1[i] = zpos[sup].max() if sup.size else 0.0
p0 = p * (1 + 0.2 * np.clip(mach, 0, None) ** 2) ** 3.5
p0_in = p0[np.argmin(np.abs(ax_x - 0.2))]

# ================= figures =================
fig, axs = plt.subplots(3, 1, figsize=(15, 12), sharex=True, constrained_layout=True)
for a in axs:
    for b in lev_bounds:
        a.axvline(b, color="0.6", ls=":", lw=1)
    for xe in EXP_PEAKS:
        a.axvline(xe, color="crimson", ls="--", lw=0.7, alpha=0.5)
axs[0].plot(ax_x, rho / 1.1698, "k-", lw=1.2, label="CFD mean rho/rho_amb")
axs[0].plot(exp[:, 0], exp[:, 1], "o-", color="crimson", ms=3, lw=0.8, label="Panda fig4 (M=1.43)")
axs[0].plot(cfd, np.interp(cfd, ax_x, rho / 1.1698), "kv", ms=7, label="CFD peaks (consensus)")
axs[0].set_ylabel("rho / rho_amb"); axs[0].legend(); axs[0].set_xlim(0, 10)
axs[1].plot(ax_x, p / 1e3, "b-", lw=1.2)
axs[1].set_ylabel("mean p [kPa]")
axs[2].plot(ax_x, mach, "g-", lw=1.2)
axs[2].set_ylabel("mean centerline Mach"); axs[2].set_xlabel("x/D")
axs[0].set_title("m142 centerline time-mean vs experiment — dotted: AMR level boundaries "
                 f"(x/D = {', '.join(str(b) for b in lev_bounds)}); dashed red: exp peaks")
fig.savefig(f"{BASE}/fig1_centerline_peaks_amr.png", dpi=140)

fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
axs[0].axhline(0, color="0.5")
axs[0].errorbar(n, E, yerr=np.sqrt(np.maximum(spread / 2, 0.02)**2 + EXP_ERR**2),
                fmt="ko-", capsize=3)
axs[0].set_xlabel("peak index n"); axs[0].set_ylabel("E_n = x_cfd - x_exp  [D]")
axs[0].set_title("cumulative position error")
axs[1].plot(n[:-1], Le, "o-", color="crimson", label="exp")
axs[1].plot(n[:-1], Lc, "ks-", label="CFD")
axs[1].set_xlabel("cell index n"); axs[1].set_ylabel("L_n [D]"); axs[1].legend()
axs[1].set_title("cell length")
axs[2].axhline(0, color="0.5")
axs[2].plot(n[:-1], eps * 100, "ko-")
axs[2].set_xlabel("cell index n"); axs[2].set_ylabel("eps_n [%]")
axs[2].set_title("cell-length error")
for a in axs[1:]:
    pass
fig.savefig(f"{BASE}/fig2_En_Ln_eps.png", dpi=140)

fig, axs = plt.subplots(3, 1, figsize=(15, 11), sharex=True, constrained_layout=True)
for a in axs:
    for b in lev_bounds:
        a.axvline(b, color="0.6", ls=":", lw=1)
    for xc in cfd:
        a.axvline(xc, color="k", ls="--", lw=0.6, alpha=0.5)
axs[0].plot(xp_, rs, "b-", label="shear-layer radius r_s (max |dU/dz|)")
axs[0].plot(xp_, rM1, "g-", label="M=1 core radius (mean field)")
axs[0].set_ylabel("r / D"); axs[0].legend(); axs[0].set_ylim(0, 2); axs[0].set_xlim(0, 10)
axs[1].plot(xp_, dom, "m-")
axs[1].set_ylabel("vorticity thickness / D")
axs[2].plot(ax_x, p0 / p0_in, "r-")
axs[2].set_ylabel("mean p0 / p0(0.2D)"); axs[2].set_xlabel("x/D"); axs[2].set_ylim(0, 1.1)
axs[0].set_title("jet-scale diagnostics — dotted: AMR boundaries; dashed: CFD shock peaks")
fig.savefig(f"{BASE}/fig3_jet_scales.png", dpi=140)
print("figures written to", BASE)
