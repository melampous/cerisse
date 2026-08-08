"""F5: Mach-disk region structure, Edgington-Mitchell et al. (2014) style —
mean axial velocity, axial-velocity rms, centreline profiles, mean vorticity.
Source: 2D RZ NPR=5.746 plateau statistics (D_e/256), MEAN+SQR fields.

F6: numerical schlieren triptych, Zapryagaev (2015) Fig. 3 style —
(a) RZ mean field |grad rho_mean| (long-exposure analogue),
(b) RZ instantaneous, (c) m180 3D instantaneous meridional slice.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from lit_common import (D, GAM, RGAS, PINF, BASE, load2d, axis_state,
                        mach_field, load_m180_slice)

plt.rcParams.update({"font.size": 11, "xtick.direction": "in",
                     "ytick.direction": "in"})

c = load2d("LIT2D_sweep_5.746.npz")
q = np.load(f"{BASE}/data/LIT2D_sweep_5.746_sqr.npz")
s = load_m180_slice()
dx = c["r"][1] - c["r"][0]

uz = c["y_velocityMEAN"]; ur = c["x_velocityMEAN"]
uz_rms = np.sqrt(np.clip(q["y_velocitySQR"] - uz**2, 0, None))
Mf = mach_field(c)
# exit velocity scale (sonic, from the axis exit state)
z, p_ax, T_ax, rho_ax, uz_ax, M_ax = axis_state(c)
i0 = np.argmin(np.abs(z - 0.05))
UE = float(uz_ax[i0])
print(f"U_e = {UE:.1f} m/s")

# mean azimuthal vorticity  omega_theta = d(ur)/dz - d(uz)/dr   (r,z grids)
dur_dz = np.gradient(ur, dx * D, axis=1)
duz_dr = np.gradient(uz, dx * D, axis=0)
omt = (dur_dz - duz_dr) * D / UE

XM, RM = 4.0, 1.6
jm = int(XM / dx); im = int(RM / dx)
ext = [0, XM, 0, RM]


def mirror(A):
    return A[:im, :jm]


fig = plt.figure(figsize=(12.8, 7.6))
gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.18,
                      left=0.06, right=0.98, top=0.95, bottom=0.07)

axA = fig.add_subplot(gs[0, 0])
h = axA.imshow(mirror(uz) / UE, origin="lower", extent=ext, aspect="equal",
               cmap="RdYlBu_r", vmin=-0.2, vmax=2.8)
axA.contour(c["z"][:jm], c["r"][:im], mirror(Mf), levels=[1.0],
            colors="k", linewidths=0.9)
plt.colorbar(h, ax=axA, pad=0.01, label=r"$\bar u_x/U_e$")
axA.set_title("(a) mean axial velocity + sonic line", fontsize=10)
axA.set_ylabel(r"$r/D_e$")

axB = fig.add_subplot(gs[0, 1])
h = axB.imshow(mirror(uz_rms) / UE, origin="lower", extent=ext, aspect="equal",
               cmap="magma", vmin=0, vmax=0.5)
axB.contour(c["z"][:jm], c["r"][:im], mirror(Mf), levels=[1.0],
            colors="w", linewidths=0.7)
plt.colorbar(h, ax=axB, pad=0.01, label=r"$u_{x,\mathrm{rms}}/U_e$")
axB.set_title("(b) axial-velocity rms — outer + internal shear layers",
              fontsize=10)

axC = fig.add_subplot(gs[1, 0])
axC.plot(z, uz_ax / UE, color="#08306b", lw=1.8, label=r"$\bar u_x/U_e$")
rms_ax = np.sqrt(np.clip(q["y_velocitySQR"][0] - uz_ax**2, 0, None)) / UE
axC.plot(z, rms_ax, color="#d95f02", lw=1.5, label=r"$u_{x,\mathrm{rms}}/U_e$")
axC.axhline(0, color="0.5", lw=0.7)
sel_w = (z > 1.0) & (z < 6.0)
umin = uz_ax[sel_w].min(); xmin_ = z[sel_w][uz_ax[sel_w].argmin()]
note = (rf"$\min \bar u_x = {umin:.0f}$ m/s at $x/D_e={xmin_:.2f}$" "\n"
        + ("mean recirculation (RZ); Edgington-Mitchell exp. (NPR 4.2): none"
           if umin < 0 else "no mean reverse flow"))
axC.text(0.98, 0.55, note, transform=axC.transAxes, ha="right",
         fontsize=8.5, color="0.25")
axC.set_xlim(0, 6); axC.set_xlabel(r"$x/D_e$")
axC.set_ylabel(r"centreline velocity $/U_e$")
axC.legend(fontsize=9, loc="upper right")
axC.set_title("(c) centreline mean and rms (Edgington-Mitchell Fig. 3 analogue)",
              fontsize=10)

axD = fig.add_subplot(gs[1, 1])
h = axD.imshow(mirror(omt), origin="lower", extent=ext, aspect="equal",
               cmap="RdBu_r", vmin=-8, vmax=8)
axD.contour(c["z"][:jm], c["r"][:im], mirror(Mf), levels=[1.0],
            colors="k", linewidths=0.7)
plt.colorbar(h, ax=axD, pad=0.01,
             label=r"$\bar\omega_\theta D_e/U_e$")
axD.set_title("(d) mean azimuthal vorticity — slip line / annular layer",
              fontsize=10)
axD.set_xlabel(r"$x/D_e$"); axD.set_ylabel(r"$r/D_e$")
axA.set_xlabel(r"$x/D_e$"); axB.set_xlabel(r"$x/D_e$")
axB.set_ylabel(r"$r/D_e$")

fig.suptitle("Mach-disk region, NPR = 5.746 — RZ plateau statistics "
             r"($D_e/256$, 1.25 ms window)", fontsize=11)
for ext_ in ("png", "pdf"):
    fig.savefig(f"{BASE}/litfig_m180_diskregion.{ext_}", dpi=250)
plt.close(fig)

# ================= F6: schlieren triptych =================
def schl(rho, dxm):
    gx = np.gradient(rho, dxm, axis=1)
    gy = np.gradient(rho, dxm, axis=0)
    g = np.sqrt(gx**2 + gy**2)
    return g


XS, RS = 6.0, 1.8
jm2 = int(XS / dx); im2 = int(RS / dx)
g_mean = schl(c["DensityMEAN"], dx * D)[:im2, :jm2]
g_inst = schl(c["Density"], dx * D)[:im2, :jm2]
dx3 = s["x"][1] - s["x"][0]
i3m = int(XS / dx3)
g3 = schl(np.nan_to_num(s["Density"], nan=np.nanmedian(s["Density"])),
          dx3 * D)
jsel = (s["y"] > -RS) & (s["y"] < RS)
g3c = g3[jsel][:, :i3m]

fig, axes = plt.subplots(3, 1, figsize=(11.6, 8.4), sharex=True)
titles = [r"(a) $|\nabla\bar\rho|$, RZ mean ($D_e/256$) — long-exposure analogue",
          r"(b) $|\nabla\rho|$, RZ instantaneous — short-exposure analogue",
          r"(c) $|\nabla\rho|$, 3D instantaneous meridional slice ($D_e/128$)"]
full_r = np.concatenate([-c["r"][:im2][::-1], c["r"][:im2]])
for ax, G, ttl in zip(axes[:2], [g_mean, g_inst], titles[:2]):
    Gm = np.vstack([G[::-1], G])
    ax.imshow(Gm, origin="lower", extent=[0, XS, -RS, RS], aspect="equal",
              cmap="gray_r", norm=LogNorm(vmin=2, vmax=2000))
    ax.set_title(ttl, fontsize=10)
    ax.set_ylabel(r"$r/D_e$")
axes[2].imshow(g3c, origin="lower", extent=[0, XS, -RS, RS], aspect="equal",
               cmap="gray_r", norm=LogNorm(vmin=2, vmax=2000))
axes[2].set_title(titles[2], fontsize=10)
axes[2].set_ylabel(r"$y/D_e$")
axes[2].set_xlabel(r"$x/D_e$")
fig.suptitle("NPR = 5.746 numerical schlieren — Zapryagaev et al. (2015) "
             "Fig. 3 layout (exp. NPR = 5.0)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
for ext_ in ("png", "pdf"):
    fig.savefig(f"{BASE}/litfig_m180_schlieren.{ext_}", dpi=250)
plt.close(fig)
print("saved litfig_m180_diskregion + litfig_m180_schlieren")
