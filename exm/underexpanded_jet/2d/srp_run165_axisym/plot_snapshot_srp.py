#!/usr/bin/env python3
"""Quick schlieren + Mach number snapshot for 2D axisymmetric SRP Run-165.

Usage:  python plot_snapshot_srp.py [plt_dir]
"""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level(40)

GAM   = 1.4
R_GAS = 287.06

def main():
    if len(sys.argv) > 1:
        plt_dir = sys.argv[1]
    else:
        # find newest plt?????
        cands = sorted([d for d in os.listdir('.') if d.startswith('plt') and d[3:].isdigit()])
        if not cands:
            print("No plt directory found"); sys.exit(1)
        plt_dir = cands[-1]

    print(f"Loading {plt_dir} ...")
    ds = yt.load(plt_dir)

    # Use covering_grid at finest level so AMR is fully resolved
    lmax = ds.index.max_level
    dims = ds.domain_dimensions * ds.refine_by**lmax
    cg   = ds.covering_grid(level=lmax, left_edge=ds.domain_left_edge, dims=dims)

    rho = np.asarray(cg[("boxlib","Density")])
    mx  = np.asarray(cg[("boxlib","Xmom")])
    my  = np.asarray(cg[("boxlib","Ymom")])
    E   = np.asarray(cg[("boxlib","Energy")])
    # covering_grid may be 3D with trailing 1 ; take [:,:,0] or reshape to 2D
    while rho.ndim > 2:
        rho = rho[..., 0]; mx = mx[..., 0]; my = my[..., 0]; E = E[..., 0]
    rho = np.maximum(rho, 1e-30)
    ur  = mx / rho
    uz  = my / rho
    ke  = 0.5 * (ur**2 + uz**2)
    p   = (GAM - 1.0) * (E - rho * ke)
    T   = p / (rho * R_GAS)
    c   = np.sqrt(GAM * R_GAS * T)
    M   = np.sqrt(ur**2 + uz**2) / c

    # In AMReX 2D RZ: axis 0 is r (i), axis 1 is z (j). yt returns data as [i, j, 1].
    # For a plot with r on x-axis and z on y-axis, transpose so rows=z, cols=r.
    rho_T = rho.T
    M_T   = M.T

    # Numerical schlieren (|∇ρ|/ρ), normalised (exp scaling)
    gz, gr = np.gradient(rho_T)
    dx_r = (ds.domain_right_edge[0].v - ds.domain_left_edge[0].v) / rho_T.shape[1]
    dx_z = (ds.domain_right_edge[1].v - ds.domain_left_edge[1].v) / rho_T.shape[0]
    gr /= dx_r; gz /= dx_z
    gmag = np.sqrt(gr**2 + gz**2) / rho_T
    # Exp-scaled schlieren
    s = np.exp(-5.0 * gmag / max(gmag.max(), 1e-30))
    # flip to light-on-dark convention used in the thesis figs
    schlieren = 1.0 - s

    # Mirror about the axis (r -> -r) for visualisation
    rho_mirror = np.concatenate([rho_T[:, ::-1], rho_T], axis=1)
    M_mirror   = np.concatenate([M_T[:, ::-1],   M_T  ], axis=1)
    sch_mirror = np.concatenate([schlieren[:, ::-1], schlieren], axis=1)

    rlo = -float(ds.domain_right_edge[0].v) * 1000  # mm
    rhi = -rlo
    zlo = float(ds.domain_left_edge[1].v) * 1000
    zhi = float(ds.domain_right_edge[1].v) * 1000

    fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharey=True)
    extent = [rlo, rhi, zlo, zhi]

    im0 = axes[0].imshow(sch_mirror, origin="lower", extent=extent,
                          cmap="gray", aspect="equal")
    axes[0].set_title(f"Schlieren  —  {plt_dir}  t={float(ds.current_time):.3e} s")
    axes[0].set_xlabel("r [mm]"); axes[0].set_ylabel("z [mm]")
    plt.colorbar(im0, ax=axes[0], shrink=0.7, label="1 - exp(-5|∇ρ|/ρ / max)")

    im1 = axes[1].imshow(M_mirror, origin="lower", extent=extent,
                          cmap="jet", aspect="equal", vmin=0, vmax=5)
    axes[1].set_title(f"Mach number  —  {plt_dir}")
    axes[1].set_xlabel("r [mm]")
    plt.colorbar(im1, ax=axes[1], shrink=0.7, label="M")

    # jet exit marker
    for ax in axes:
        ax.axhline(0, color="w", lw=0.3, alpha=0.3)
        ax.axvline( 6.3, color="w", lw=0.3, ls="--", alpha=0.5)
        ax.axvline(-6.3, color="w", lw=0.3, ls="--", alpha=0.5)

    out = f"snapshot_{plt_dir}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Saved {out}  (shape={rho.shape}, t={float(ds.current_time):.3e} s)")

if __name__ == "__main__":
    main()
