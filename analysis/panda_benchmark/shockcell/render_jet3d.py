#!/usr/bin/env python3
"""Perspective 3D view of the underexpanded jet from an extracted block.

Combines three robust elements, each chosen to survive the resolution of
the extracted block:
  (i)   a DENSITY isosurface  rho/rho_j = const  -> the barrel shock and
        shock-cell surfaces (smooth, resolution-robust);
  (ii)  a VORTEX-CORE isosurface  max(-lambda_2,0) = const, coloured by
        |omega| De/Uj  -> the shear-layer structures (needs fine data);
  (iii) a numerical-schlieren slice on z = 0 -> the shock system in a cut.
Plus the nozzle-lip circle for scale.

Usage: python3 render_jet3d.py BLOCK3D_L4_lite.npz [rho_iso] [core_pct]
"""
import sys
from pathlib import Path
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True
HERE = Path(__file__).parent
FN = sys.argv[1] if len(sys.argv) > 1 else "BLOCK3D_L4_lite.npz"
RHO_ISO = float(sys.argv[2]) if len(sys.argv) > 2 else 0.55
CORE_PCT = float(sys.argv[3]) if len(sys.argv) > 3 else 97.5
TAG = Path(FN).stem.replace("BLOCK3D_", "")
D, UJ, RHO_J = 0.0254, 414.2, 1.6413

d = np.load(HERE / FN)
x, y, z = np.asarray(d["x"], float), np.asarray(d["y"], float), np.asarray(d["z"], float)
dx = float(d["dx"])
if "core" in d:                      # lite file: already non-dimensional
    core = np.nan_to_num(np.asarray(d["core"], np.float32))
    om = np.nan_to_num(np.asarray(d["om"], np.float32))
else:
    core = np.nan_to_num(np.maximum(-d["lam2"], 0.0)) * (D/UJ)**2
    om = np.nan_to_num(d["omag"]) * (D/UJ)
rho = np.nan_to_num(np.asarray(d["rho"], np.float32), nan=RHO_J) / RHO_J
print("block %s  t=%.3f ms  dx=%.1f um" % (core.shape, float(d["time"])*1e3, dx*1e6))

grid = pv.ImageData()
grid.dimensions = np.array(core.shape) + 1
grid.origin = (x[0], y[0], z[0])
grid.spacing = (dx/D, dx/D, dx/D)
grid.cell_data["core"] = core.ravel(order="F")
grid.cell_data["om"] = om.ravel(order="F")
grid.cell_data["rho"] = rho.ravel(order="F")
pts = grid.cell_data_to_point_data()

# (i) density isosurface = the shock-cell / barrel-shock envelope
iso_r = pts.contour([RHO_ISO], scalars="rho").clip("x", value=2.3, invert=True)
print("density isosurface rho/rho_j=%.2f : %d cells" % (RHO_ISO, iso_r.n_cells))

# (ii) vortex cores
pos = core[core > 0]
lev = float(np.percentile(pos, CORE_PCT)) if pos.size else 1.0
iso_c = pts.contour([lev], scalars="core")
if iso_c.n_cells:
    iso_c = iso_c.smooth(n_iter=30, relaxation_factor=0.1)
print("vortex isosurface max(-lam2)=%.0f (p%.0f) : %d cells" % (lev, CORE_PCT, iso_c.n_cells))

# (iii) schlieren on z = 0
kz = int(np.argmin(np.abs(z)))
gy, gx = np.gradient(rho[:, :, kz], dx/D, dx/D)
gm = np.sqrt(gx**2 + gy**2)
shl = np.exp(-6.0*gm/np.percentile(gm, 99.5))
ny = shl.shape[1]; jhalf = int(np.searchsorted(y, 0.0))
sl = pv.ImageData()
sl.dimensions = (shl.shape[0]+1, jhalf+1, 1)
sl.origin = (x[0], y[0], 0.0)
sl.spacing = (dx/D, dx/D, 1.0)
sl.cell_data["schlieren"] = shl[:, :jhalf].ravel(order="F")

th = np.linspace(0, 2*np.pi, 400)
lip = pv.lines_from_points(np.c_[np.zeros_like(th), 0.5*np.cos(th), 0.5*np.sin(th)])

p = pv.Plotter(off_screen=True, window_size=(1800, 1100))
p.set_background("white")
p.add_mesh(sl, scalars="schlieren", cmap="gray", clim=(0, 1),
           show_scalar_bar=False, opacity=0.85, lighting=False)
p.add_mesh(iso_r, color="#8fa8c8", opacity=0.22, smooth_shading=True,
           specular=0.2, show_scalar_bar=False)
if iso_c.n_cells:
    p.add_mesh(iso_c, scalars="om", cmap="turbo", clim=(2, 22), opacity=1.0,
               smooth_shading=True, specular=0.4, specular_power=20,
               scalar_bar_args=dict(title="|w| De/Uj  ", vertical=True,
                                    position_x=0.905, position_y=0.30,
                                    height=0.42, width=0.028,
                                    title_font_size=16, label_font_size=14,
                                    color="black"))
p.add_mesh(lip, color="black", line_width=5)
ax_line = pv.lines_from_points(np.c_[[0.0, float(x[-1])], [0.0,0.0], [0.0,0.0]])
p.add_mesh(ax_line, color="#777777", line_width=2)
p.camera_position = [(-1.8, -5.4, 3.4), (2.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
p.camera.zoom(1.55)
p.add_text("t = %.2f ms" % (float(d["time"])*1e3), position="upper_left",
           font_size=11, color="black")
out = HERE / ("paperfig_jet3d_%s.png" % TAG)
p.screenshot(str(out))
print("wrote", out.name)
