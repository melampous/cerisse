#!/usr/bin/env python3
"""Reliable under-expanded-jet visualization for the IBM-solid NPR=3 case.

Three panels:
  (1) numerical schlieren  exp(-k|grad rho|/max)  -- crisp shock structure
  (2) Mach with shock-cell contour lines, axis-mirrored full jet
  (3) centerline Mach(z) with first-Mach-disk + Ashkenas-Sherman markers

Usage: python3 render_jet.py [plotfile] [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level(50)

PLT = sys.argv[1] if len(sys.argv) > 1 else "plot_uni/plt10000"
OUT = sys.argv[2] if len(sys.argv) > 2 else "ibm_jet_npr3_schlieren.png"

GAMMA = 1.4
DE = 0.01           # jet exit diameter = 2*r_jet = 2*5mm  (m)
NPR = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
MVMAX = max(2.2, 1.1 * np.sqrt(NPR) + 1.4)   # Mach colorbar / xlim scales with NPR
MD_ZMAX = 65.0                                # MD search upper bound (mm), covers high NPR
# Ashkenas-Sherman: x_MD/De = 0.67*sqrt(NPR)  (NPR = p0_jet/p_amb here ~ p_jet/p_amb*... use static-ratio form)
# We use the same target the case was validated against: 11.6 mm
AS_MD_MM = 0.67 * np.sqrt(NPR) * DE * 1000.0   # ~11.6 mm

ds = yt.load(PLT)
ds.force_periodicity()   # avoid covering_grid edge round-off error
# NOTE: code_length IS already in metres; .to("m") wrongly rescales by 0.01, so use raw .d
dom_lo = ds.domain_left_edge.d
dom_hi = ds.domain_right_edge.d
# render at the finest available level (uniform -> 0; AMR -> max_level)
LEVEL = int(getattr(ds.index, "max_level", 0))
dims = np.array(ds.domain_dimensions) * (2 ** LEVEL)
nr, nz = int(dims[0]), int(dims[1])

cg = ds.covering_grid(level=LEVEL, left_edge=ds.domain_left_edge, dims=dims)
rho = np.array(cg["Density"][:, :, 0])
p   = np.array(cg["pressure"][:, :, 0])
ux  = np.array(cg["x_velocity"][:, :, 0])
uy  = np.array(cg["y_velocity"][:, :, 0])
sld = np.array(cg["sld"][:, :, 0])    # 1 in solid, 0 in fluid (or fractional)

# coordinate axes (cell centers)
r = np.linspace(dom_lo[0], dom_hi[0], nr, endpoint=False) + 0.5 * (dom_hi[0] - dom_lo[0]) / nr
z = np.linspace(dom_lo[1], dom_hi[1], nz, endpoint=False) + 0.5 * (dom_hi[1] - dom_lo[1]) / nz
dr = (dom_hi[0] - dom_lo[0]) / nr
dz = (dom_hi[1] - dom_lo[1]) / nz

# arrays are [nr, nz]; transpose to [nz, nr] for imshow (z vertical, r horizontal)
def T(a):
    return a.T

rhoT, pT, uxT, uyT, sldT = T(rho), T(p), T(ux), T(uy), T(sld)
mach = np.sqrt(uxT**2 + uyT**2) / np.sqrt(GAMMA * np.maximum(pT, 1e-6) / np.maximum(rhoT, 1e-9))

# numerical schlieren: gradient magnitude of density
grad_z, grad_r = np.gradient(rhoT, dz, dr)
gmag = np.sqrt(grad_z**2 + grad_r**2)
solid_mask = sldT > 0.5
gmag_fluid = gmag.copy()
gmag_fluid[solid_mask] = 0.0
gnorm = gmag_fluid / (np.percentile(gmag_fluid[~solid_mask], 99.5) + 1e-30)
schlieren = np.exp(-8.0 * np.clip(gnorm, 0, 1))   # 1 = quiescent (white), dark = shock

mach_masked = np.ma.masked_where(solid_mask, mach)
schlieren_masked = np.ma.masked_where(solid_mask, schlieren)

# axis-mirror (RZ -> full jet)  : stack r<0 (flip) | r>=0
def mirror(field):
    return np.concatenate([field[:, ::-1], field], axis=1)
r_full = np.concatenate([-r[::-1], r])
extent_full = [r_full[0] * 1000, r_full[-1] * 1000, z[0] * 1000, z[-1] * 1000]  # mm

# near-field window for plotting: include the solid plate (z<0) up to z=80 mm
ZMAX_MM = 80.0
ZMIN_MM = -16.0     # show the solid plate (z in [-15, 0] mm) + sliver of cavity below
zsel = z * 1000 <= ZMAX_MM
iz = np.argmax(z * 1000 > ZMAX_MM) if np.any(z * 1000 > ZMAX_MM) else nz
PLATE_Z0, PLATE_Z1 = -15.0, 0.0    # plate spans z in [-15,0] mm, r in [0,50] mm
PLATE_RMAX = 50.0
RJET_MM = 0.5 * DE * 1000.0        # 5 mm (r_jet = De/2)

# centerline Mach (innermost fluid column, r ~ dr/2)
mc = mach[:, 0].copy()
# first Mach-disk = first local minimum of centerline Mach above z>2mm in the supersonic train
zc_mm = z * 1000
mask_search = zc_mm > 2.0
# find local minima
mins = []
for k in range(2, nz - 2):
    if zc_mm[k] < 3.0 or zc_mm[k] > MD_ZMAX:
        continue
    if mc[k] < mc[k - 1] and mc[k] <= mc[k + 1] and mc[k] < 1.6:
        mins.append((zc_mm[k], mc[k]))
md_mm = mins[0][0] if mins else float("nan")

# ---------------- plot ----------------
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

EXTENT = [r_full[0]*1000, r_full[-1]*1000, z[0]*1000, ZMAX_MM]
SOLID_COLOR = "#c2b280"   # sand: the IBM solid plate

def draw_solid(ax, label=False):
    """Overlay the IBM solid plate + outline + jet-exit patch on a contour panel."""
    solid_ov = np.where(mirror(solid_mask.astype(float))[:iz] > 0.5, 1.0, np.nan)
    ax.imshow(solid_ov, origin="lower", extent=EXTENT,
              cmap=ListedColormap([SOLID_COLOR]), vmin=0, vmax=1, aspect="auto", zorder=3)
    # crisp plate outline (mirrored rectangle), broken at the jet aperture
    ax.add_patch(Rectangle((-PLATE_RMAX, PLATE_Z0), 2*PLATE_RMAX, PLATE_Z1-PLATE_Z0,
                           fill=False, ec="black", lw=1.0, zorder=4))
    # jet-exit aperture on the top face: |r| < r_jet at z=0
    ax.plot([-RJET_MM, RJET_MM], [0, 0], color="magenta", lw=3.0, zorder=5)
    for rr in (-3.2, 0.0, 3.2):
        ax.annotate("", xy=(rr, 7), xytext=(rr, 0.4),
                    arrowprops=dict(arrowstyle="-|>", color="magenta", lw=1.4), zorder=6)
    if label:
        ax.text(0, -7.5, "solid plate (IBM, slip wall)", ha="center", va="center",
                fontsize=8, color="black", zorder=7)
        ax.text(RJET_MM+1.5, 9, "jet exit\nM=1 (from solid)", ha="left", va="center",
                fontsize=8, color="magenta", zorder=7)

fig = plt.figure(figsize=(15, 9.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.05, 0.9], wspace=0.32)

# Panel 1: schlieren (mirrored)
ax1 = fig.add_subplot(gs[0, 0])
im1 = ax1.imshow(mirror(schlieren_masked)[:iz], origin="lower", extent=EXTENT,
                 cmap="gray", vmin=0, vmax=1, aspect="auto")
draw_solid(ax1, label=True)
ax1.axhline(md_mm, color="red", lw=1.2, ls="--", label=f"1st Mach disk {md_mm:.1f} mm")
ax1.axhline(AS_MD_MM, color="cyan", lw=1.2, ls=":", label=f"A-S {AS_MD_MM:.1f} mm")
ax1.set_title("Numerical schlieren  |∇ρ|\n(barrel shock + Mach disks + diamonds)")
ax1.set_xlabel("r (mm)")
ax1.set_ylabel("z (mm)")
ax1.set_ylim(ZMIN_MM, ZMAX_MM)
ax1.legend(loc="upper right", fontsize=8)

# Panel 2: Mach (mirrored) with contour lines
ax2 = fig.add_subplot(gs[0, 1])
machM = mirror(mach_masked)[:iz]
im2 = ax2.imshow(machM, origin="lower", extent=EXTENT,
                 cmap="turbo", vmin=0, vmax=MVMAX, aspect="auto")
RR, ZZ = np.meshgrid(r_full * 1000, z[:iz] * 1000)
try:
    ax2.contour(RR, ZZ, np.ma.filled(machM, 0), levels=[1.0], colors="white", linewidths=0.8)
except Exception:
    pass
draw_solid(ax2)
ax2.axhline(md_mm, color="red", lw=1.0, ls="--")
ax2.axhline(AS_MD_MM, color="cyan", lw=1.0, ls=":")
ax2.set_title("Mach (axis-mirrored full jet)\nwhite line = sonic M=1")
ax2.set_xlabel("r (mm)")
ax2.set_ylim(ZMIN_MM, ZMAX_MM)
cb = fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.02)
cb.set_label("Mach")

# Panel 3: centerline Mach(z)
ax3 = fig.add_subplot(gs[0, 2])
sel3 = (zc_mm <= ZMAX_MM) & (zc_mm >= ZMIN_MM)
ax3.axhspan(PLATE_Z0, PLATE_Z1, color=SOLID_COLOR, alpha=0.6, zorder=0)
ax3.text(1.3, PLATE_Z0/2, "solid plate", ha="center", va="center", fontsize=8, zorder=1)
ax3.axhline(0.0, color="magenta", lw=2.0, zorder=2)
ax3.plot(mc[sel3], zc_mm[sel3], "b-", lw=1.6, zorder=3)
ax3.axvline(1.0, color="gray", ls=":", lw=1)
ax3.axhline(md_mm, color="red", lw=1.4, ls="--", label=f"1st MD {md_mm:.1f} mm")
ax3.axhline(AS_MD_MM, color="cyan", lw=1.4, ls=":", label=f"A-S {AS_MD_MM:.1f} mm")
ax3.set_title("centerline Mach(z)\nshock-cell train")
ax3.set_xlabel("centerline Mach")
ax3.set_ylabel("z (mm)")
ax3.set_xlim(0, MVMAX + 0.2)
ax3.set_ylim(ZMIN_MM, ZMAX_MM)
ax3.grid(alpha=0.3)
ax3.legend(loc="upper right", fontsize=8)

err = (md_mm - AS_MD_MM) / AS_MD_MM * 100.0
fig.suptitle(
    f"IBM-solid under-expanded jet  NPR={NPR:.0f}  (M_exit=1, sonic)   "
    f"1st Mach disk {md_mm:.1f} mm vs Ashkenas-Sherman {AS_MD_MM:.1f} mm  ({err:+.0f}%)   "
    f"[{PLT.split('/')[-1]}]",
    fontsize=12, y=0.99)
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"saved {OUT}   x_MD={md_mm:.2f}mm  A-S={AS_MD_MM:.2f}mm  err={err:+.1f}%  (mins={mins[:4]})")
