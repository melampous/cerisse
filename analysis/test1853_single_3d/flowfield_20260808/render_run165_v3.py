#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

PLT, OUT = sys.argv[1], sys.argv[2]
ds = yt.load(PLT)
ad = ds.all_data()
print("global sld max:", float(ad[("boxlib","sld")].max()),
      " ghs max:", float(ad[("boxlib","ghs")].max()))
for xp in [0.15, 0.25, 0.35]:
    try:
        pt = ds.point([xp, 0.0, 0.0])
        print("probe x=%+.3f m: T=%7.2f K  u=%+8.1f  p=%9.1f  sld=%.1f"
              % (xp, float(pt[("boxlib","temperature")][0]),
                 float(pt[("boxlib","x_velocity")][0]),
                 float(pt[("boxlib","pressure")][0]),
                 float(pt[("boxlib","sld")][0])))
    except Exception as e:
        print("probe x=%+.3f m: no data (%s)" % (xp, type(e).__name__))

xlo, ylo, zlo = [float(v) for v in ds.domain_left_edge.d]
xhi, yhi, zhi = [float(v) for v in ds.domain_right_edge.d]
NX, NZ = 1152, 768
slc = ds.slice("y", 0.0)
frb = slc.to_frb(width=(zhi-zlo, "code_length"), resolution=(NZ, NX),
                 center=[(xlo+xhi)/2, 0.0, 0.0],
                 height=(xhi-xlo, "code_length"))
def g(name):
    a = np.asarray(frb[("boxlib", name)], dtype=np.float64)
    return a
rho = g("Density")
print("frb shape:", rho.shape)   # expect (NX, NZ) rows=x or (NZ,NX); disambiguate
T = g("temperature"); p = g("pressure")
u = g("x_velocity"); v = g("y_velocity"); w = g("z_velocity")
sld = g("sld"); ghs = g("ghs")
if rho.shape == (NZ, NX):        # rows = image-y = x per yt for normal 'y'? enforce rows=x
    rho, T, p, u, v, w, sld, ghs = [a.T for a in (rho, T, p, u, v, w, sld, ghs)]
# now arrays are (NX, NZ): index [x, z]
mach = np.sqrt(u*u+v*v+w*w)/np.sqrt(1.4*287.1*np.maximum(T, 1.0))
gx, gz = np.gradient(rho, (xhi-xlo)/NX, (zhi-zlo)/NZ)
schl = np.log10(1.0+np.hypot(gx, gz))
ext = [xlo*1e3, xhi*1e3, zlo*1e3, zhi*1e3]
plt.rcParams.update({"font.size": 9})
fig, axs = plt.subplots(2, 2, figsize=(13.6, 9.2), sharex=True, sharey=True)
fig.subplots_adjust(left=0.05, right=0.97, top=0.93, bottom=0.06,
                    hspace=0.10, wspace=0.14)
panels = [(schl, "gray_r", None, "numerical schlieren"),
          (mach, "turbo", (0, 5.0), "Mach number"),
          (T, "inferno", (10, 360), "temperature  [K]"),
          (np.log10(np.maximum(p, 1.0)), "viridis", (2.0, 5.6), "$\\log_{10}\\,p$  [Pa]")]
smask = (sld >= 0.5)
for ax, (F, cmap, clim, title) in zip(axs.flat, panels):
    im = ax.imshow(F.T, origin="lower", extent=ext, cmap=cmap,
                   aspect="equal", interpolation="nearest")
    if clim: im.set_clim(*clim)
    ax.contourf(np.linspace(ext[0], ext[1], NX), np.linspace(ext[2], ext[3], NZ),
                smask.T.astype(float), levels=[0.5, 1.5], colors=["#909090"])
    ax.contour(np.linspace(ext[0], ext[1], NX), np.linspace(ext[2], ext[3], NZ),
               smask.T.astype(float), levels=[0.5], colors="k", linewidths=0.6)
    fig.colorbar(im, ax=ax, fraction=0.036, pad=0.015)
    ax.set_title(title, fontsize=10)
for ax in axs[1]: ax.set_xlabel("$x$  [mm]")
for ax in axs[:, 0]: ax.set_ylabel("$z$  [mm]")
fig.suptitle("Run 165 (Test 1853 SRP)  %s   t = %.4f ms   y = 0 midplane  "
             "(grey: solid body, sld)" %
             (PLT.split("/")[-1], float(ds.current_time)*1e3), fontsize=11)
fig.savefig(OUT + "_fields.png", dpi=165)
print("wrote", OUT + "_fields.png")
