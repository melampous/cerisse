#!/usr/bin/env python3
"""Overlay centerline Mach(z) of several plotfiles (AMR vs uniform) to show
physics invariance under the AMR-tagging fix. Each arg is label=plotfile.

Usage: python3 compare_centerlines.py UNI-steady=plot_uni/plt10000 \
                                       UNI-t2.1e-4=plot_uni/plt02000 \
                                       AMR-t2.4e-4=plot_amrlong/plt01200
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level(50)
GAMMA = 1.4
DE = 0.01
AS_MD_MM = 0.67 * np.sqrt(3.0) * DE * 1000.0

def centerline(plotpath):
    ds = yt.load(plotpath)
    ds.force_periodicity()
    lo = ds.domain_left_edge.d            # code units == metres (NOT .to('m'))
    hi = ds.domain_right_edge.d
    lev = int(getattr(ds.index, "max_level", 0))
    dims = np.array(ds.domain_dimensions) * (2 ** lev)
    nr, nz = int(dims[0]), int(dims[1])
    cg = ds.covering_grid(level=lev, left_edge=ds.domain_left_edge, dims=dims)
    rho = np.array(cg["Density"][:, :, 0])
    p   = np.array(cg["pressure"][:, :, 0])
    ux  = np.array(cg["x_velocity"][:, :, 0])
    uy  = np.array(cg["y_velocity"][:, :, 0])
    mach = np.sqrt(ux**2 + uy**2) / np.sqrt(GAMMA * np.maximum(p, 1e-6) / np.maximum(rho, 1e-9))
    z = np.linspace(lo[1], hi[1], nz, endpoint=False) + 0.5 * (hi[1] - lo[1]) / nz
    mc = mach[0, :]                       # r=0 centerline column
    t = float(ds.current_time.d)
    return z * 1000.0, mc, t

def first_md(zmm, mc):
    mins = []
    for k in range(2, len(zmm) - 2):
        if zmm[k] < 3.0 or zmm[k] > 40.0:
            continue
        if mc[k] < mc[k-1] and mc[k] <= mc[k+1] and mc[k] < 1.6:
            mins.append(zmm[k])
    return mins[0] if mins else float("nan")

fig, ax = plt.subplots(figsize=(7.5, 8.5))
colors = ["k", "tab:blue", "tab:red", "tab:green", "tab:purple"]
for i, spec in enumerate(sys.argv[1:]):
    label, path = spec.split("=", 1)
    zmm, mc, t = centerline(path)
    md = first_md(zmm, mc)
    sel = (zmm >= 0) & (zmm <= 60)
    ls = "-" if "AMR" in label else "--"
    ax.plot(mc[sel], zmm[sel], ls, color=colors[i % len(colors)], lw=1.7,
            label=f"{label}  (t={t:.2e}s, 1st MD={md:.1f}mm)")
    print(f"{label:16s} t={t:.3e}  1st_MD={md:.2f}mm  path={path}")

ax.axvline(1.0, color="gray", ls=":", lw=1)
ax.axhline(AS_MD_MM, color="cyan", ls=":", lw=1.5, label=f"Ashkenas-Sherman {AS_MD_MM:.1f}mm")
ax.set_xlabel("centerline Mach")
ax.set_ylabel("z (mm)")
ax.set_xlim(0, 2.6)
ax.set_ylim(0, 60)
ax.grid(alpha=0.3)
ax.legend(loc="upper right", fontsize=8)
ax.set_title("Centerline Mach(z): AMR (solid) vs uniform (dashed)\n"
             "physics invariance under the AMR-tagging fix")
fig.savefig("ibm_jet_npr3_AMRvsUNI.png", dpi=130, bbox_inches="tight")
print("saved ibm_jet_npr3_AMRvsUNI.png")
