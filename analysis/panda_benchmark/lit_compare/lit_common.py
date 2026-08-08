"""Shared helpers for the literature-comparison figures (Muraoka/Li/
Zapryagaev/Edgington-Mitchell/Cheng/Duronio benchmarks vs the NPR campaign).

2D RZ arrays: shape (nr=512, nz=2560), D_e/256 spacing, r,z cell-centred.
3D m180 slice: shape (ny=768, nx=1024), D_e/128, y in [-3D,3D), x in [0,8D).
"""
import os
import numpy as np

D = 0.0254
GAM = 1.4
RGAS = 287.0
PINF = 101325.0
BASE = os.environ.get(
    "LITBASE",
    "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/lit_compare")
AUX = os.environ.get("LITAUX", BASE + "/..")

CASES_2D = [  # (label, npr, filename)
    ("sweep", 1.89293, "LIT2D_sweep_1.893.npz"),
    ("sweep", 2.394, "LIT2D_sweep_2.394.npz"),
    ("sweep", 2.800, "LIT2D_sweep_2.800.npz"),
    ("sweep", 3.273446, "LIT2D_sweep_3.273.npz"),
    ("sweep", 3.671, "LIT2D_sweep_3.671.npz"),
    ("sweep", 4.200, "LIT2D_sweep_4.200.npz"),
    ("sweep", 4.900, "LIT2D_sweep_4.900.npz"),
    ("sweep", 5.746, "LIT2D_sweep_5.746.npz"),
    ("nprfix", 7.5, "LIT2D_nprfix_7.5.npz"),
    ("nprfix", 10.0, "LIT2D_nprfix_10.npz"),
    ("nprfix", 15.0, "LIT2D_nprfix_15.npz"),
    ("nprfix", 20.0, "LIT2D_nprfix_20.npz"),
]


def load2d(fn):
    d = np.load(f"{BASE}/data/{fn}")
    dx = float(d["dx"])
    out = {k: d[k] for k in d.files}
    nr, nz = d["DensityMEAN"].shape
    out["r"] = (np.arange(nr) + 0.5) * dx / D
    out["z"] = (np.arange(nz) + 0.5) * dx / D
    return out


def axis_state(c):
    """Axis (first radial cell) mean profiles -> z/D, p, T, rho, uz, M."""
    p = c["pressureMEAN"][0]
    T = c["temperatureMEAN"][0]
    rho = c["DensityMEAN"][0]
    uz = c["y_velocityMEAN"][0]
    M = uz / np.sqrt(GAM * RGAS * T)
    return c["z"], p, T, rho, uz, M


def mach_field(c):
    """|u|-based mean Mach field (nr, nz)."""
    q = np.sqrt(c["x_velocityMEAN"] ** 2 + c["y_velocityMEAN"] ** 2)
    return q / np.sqrt(GAM * RGAS * c["temperatureMEAN"])


def pocket(z, M, mpeak_min=1.25, msub=0.95, zmax=6.0):
    """Subsonic pocket behind the first shock on the axis.
    Mach-disk classification requires a genuine deceleration (Mmin < msub),
    not a near-sonic wiggle. Returns (x_in, x_out, Mmin, censored) or None;
    censored=True when the pocket still extends at z=zmax (x_out unresolved)."""
    sel = z < zmax
    zs, Ms = z[sel], M[sel]
    ipk = int(np.argmax(Ms))
    if Ms[ipk] < mpeak_min:
        return None
    sub = np.where(Ms[ipk:] < 1.0)[0]
    if len(sub) == 0:
        return None
    i0 = ipk + sub[0]
    up = np.where(Ms[i0:] > 1.0)[0]
    censored = len(up) == 0
    i1 = len(Ms) - 1 if censored else i0 + up[0]
    Mmin = float(Ms[i0:i1 + 1].min())
    if Mmin > msub:
        return None
    return zs[i0], zs[i1], Mmin, censored


def disk_radius(c, x_in):
    """Radial half-extent of the subsonic region on the disk plane,
    taken just behind the axis M=1 downcross (x_in + 0.05 D)."""
    M = mach_field(c)
    j = int(np.argmin(np.abs(c["z"] - (x_in + 0.05))))
    col = M[:, j]
    sup = np.where(col > 1.0)[0]
    if len(sup) == 0 or sup[0] == 0:
        return 0.0
    return c["r"][sup[0] - 1]


def pitot_ratio(p, M):
    """Rayleigh Pitot-tube pressure (gamma=1.4): probe stagnation pressure."""
    M = np.clip(M, 0.0, None)
    sub = p * (1 + 0.2 * M ** 2) ** 3.5
    M2 = np.maximum(M ** 2, 1.0)
    sup = p * ((5.76 * M2) / (5.6 * M2 - 0.8)) ** 3.5 * (2.8 * M2 - 0.4) / 2.4
    return np.where(M <= 1.0, sub, sup)


def probe_avg(y, dx_over_D, probe_D=0.0236):
    """Boxcar over the Zapryagaev probe outer diameter (0.6 mm = 0.0236 D)."""
    w = max(1, int(round(probe_D / dx_over_D)) | 1)
    if w < 3:
        return y
    return np.convolve(y, np.ones(w) / w, mode="same")


def load_m180_slice():
    d = np.load(f"{BASE}/data/LIT3D_m180_slice.npz")
    dx = float(d["dx"])
    out = {k: d[k] for k in d.files}
    ny, nx = d["Density"].shape
    out["x"] = (np.arange(nx) + 0.5) * dx / D
    out["y"] = (float(d["ylo"]) + (np.arange(ny) + 0.5) * dx) / D
    return out


def m180_axis_inst(s):
    """Instantaneous axis line from the 3D slice (average of the two
    cell rows straddling y=0). Axial velocity is x_velocity in 3D."""
    j0 = np.argmin(np.abs(s["y"] + 0.0))  # nearest row
    rows = [j0 - 1, j0] if s["y"][j0] > 0 else [j0, j0 + 1]
    p = s["pressure"][rows].mean(0)
    T = s["temperature"][rows].mean(0)
    u = s["x_velocity"][rows].mean(0)
    M = u / np.sqrt(GAM * RGAS * T)
    return s["x"], p, T, u, M
