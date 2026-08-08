#!/usr/bin/env python3
"""Track the main SRP bow-shock ridge and compare overview windows.

This QA helper consumes the instantaneous XZ cache written by
postprocess_srp_tri_nozzle_reference.py.  It traces upper and lower
high-gradient branches outward from a physically constrained seed near the
jet/body interaction.  Three transverse half-heights are then compared.  For
each candidate, the upstream display edge is placed a prescribed distance
ahead of the tracked branch where it crosses the upper/lower frame.

The ridge and crop guides are QA overlays only; they are never placed on the
final publication figure.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import numpy as np
from scipy.ndimage import gaussian_filter, median_filter


BODY_DIAMETER = 0.127
BODY_CENTRE_X = 0.300
BODY_CENTRE_Z = 0.275
DEFAULT_HALF_HEIGHTS = (1.65, 1.85, 2.05)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--half-heights", nargs=3, type=float,
                        default=DEFAULT_HALF_HEIGHTS,
                        metavar=("H1", "H2", "H3"))
    parser.add_argument("--upstream-margin", type=float, default=0.25,
                        help="Target upstream margin in body diameters")
    parser.add_argument("--right-x-m", type=float, default=0.60)
    parser.add_argument("--schlieren-gref", type=float,
                        default=72.18947715759369)
    parser.add_argument("--schlieren-k", type=float, default=4.0)
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def load_cache(path):
    with np.load(str(path), allow_pickle=False) as data:
        return {
            "u": data["u"].astype(np.float64),
            "v": data["v"].astype(np.float64),
            "grad": data["instantaneous_grad_rho_magnitude"].astype(np.float64),
            "solid": data["solid"].astype(bool),
            "time_s": float(data["time_s"]),
            "plotfile": str(data["plotfile"]),
        }


def trace_branch(u, v, grad, sign):
    """Dynamic-programming ridge from |Z|=0.30 to 2.15 D_b.

    The emission score is the smoothed log density-gradient.  A tight seed at
    X=-1.45 selects the body/jet bow-shock branch rather than the unrelated
    upstream-boundary/startup structure.  Subsequent states are data driven,
    with only a smoothness and weak outward-monotonicity penalty.
    """
    oriented_z = sign * v
    z_indices = np.flatnonzero((oriented_z >= 0.30) & (oriented_z <= 2.15))
    z_indices = z_indices[np.argsort(oriented_z[z_indices])]
    x_indices = np.flatnonzero((u >= -2.32) & (u <= -0.90))
    if z_indices.size < 5 or x_indices.size < 5:
        raise RuntimeError("Cache does not cover the ridge-tracking window")

    smooth = gaussian_filter(np.nan_to_num(grad, nan=0.0), sigma=(1.7, 1.25))
    emission = np.log1p(np.maximum(smooth[np.ix_(x_indices, z_indices)], 0.0)).T
    scale = np.percentile(emission[np.isfinite(emission)], 95.0)
    emission = emission / max(float(scale), np.finfo(float).tiny)
    x_values = u[x_indices]
    z_values = oriented_z[z_indices]
    nx = len(x_values)
    nz = len(z_values)
    score = np.full((nz, nx), -np.inf, dtype=np.float64)
    back = np.full((nz, nx), -1, dtype=np.int32)

    seed_penalty = 8.0 * ((x_values + 1.45) / 0.20) ** 2
    score[0] = emission[0] - seed_penalty
    dx_grid = float(np.median(np.diff(x_values)))
    max_jump = max(2, int(round(0.055 / dx_grid)))
    for row in range(1, nz):
        for current in range(nx):
            lo = max(0, current - max_jump)
            hi = min(nx, current + max_jump + 1)
            previous = np.arange(lo, hi)
            delta = x_values[current] - x_values[previous]
            transition = 2.8 * (delta / 0.055) ** 2
            # As |Z| grows, a detached shock normally moves upstream.  Allow
            # downstream motion, but penalize it weakly rather than imposing it.
            transition += 1.4 * np.maximum(delta, 0.0) / 0.055
            candidates = score[row - 1, previous] - transition
            best_local = int(np.argmax(candidates))
            best_previous = int(previous[best_local])
            score[row, current] = emission[row, current] + candidates[best_local]
            back[row, current] = best_previous

    path = np.empty(nz, dtype=np.int32)
    path[-1] = int(np.argmax(score[-1]))
    for row in range(nz - 1, 0, -1):
        path[row - 1] = back[row, path[row]]
    ridge_x = x_values[path]
    ridge_x = median_filter(ridge_x, size=7, mode="nearest")
    ridge_strength = np.asarray([
        smooth[x_indices[path[row]], z_indices[row]] for row in range(nz)
    ])
    return {
        "z_abs": z_values,
        "z_signed": sign * z_values,
        "x": ridge_x,
        "strength": ridge_strength,
    }


def interpolate_branch(branch, half_height):
    return float(np.interp(half_height, branch["z_abs"], branch["x"]))


def candidate_windows(u, upper, lower, half_heights, target_margin, right_x_m):
    result = []
    for half_height in half_heights:
        upper_x = interpolate_branch(upper, half_height)
        lower_x = interpolate_branch(lower, half_height)
        crossing_x = min(upper_x, lower_x)
        left_x = max(float(u[0]), crossing_x - target_margin)
        actual_margin = crossing_x - left_x
        result.append({
            "half_height_Db": float(half_height),
            "upper_crossing_X": upper_x,
            "lower_crossing_X": lower_x,
            "limiting_crossing_X": crossing_x,
            "display_left_X": left_x,
            "display_right_X": (right_x_m - BODY_CENTRE_X) / BODY_DIAMETER,
            "actual_upstream_margin_Db": actual_margin,
            "display_x_m": [BODY_CENTRE_X + BODY_DIAMETER * left_x,
                            float(right_x_m)],
            "display_z_m": [BODY_CENTRE_Z - BODY_DIAMETER * half_height,
                            BODY_CENTRE_Z + BODY_DIAMETER * half_height],
        })
    return result


def draw_solid(ax, u, v, solid):
    overlay = np.ma.masked_where(~solid.T, np.ones(solid.T.shape))
    ax.pcolormesh(u, v, overlay, shading="nearest",
                  cmap=ListedColormap(["0.84"]), vmin=0.0, vmax=1.0,
                  rasterized=True, zorder=3)
    ax.contour(u, v, solid.T.astype(float), levels=[0.5], colors="0.15",
               linewidths=0.55, zorder=4)


def make_montage(data, upper, lower, candidates, args, output_dir):
    schlieren = np.exp(-args.schlieren_k * data["grad"] / args.schlieren_gref)
    schlieren[~np.isfinite(data["grad"])] = np.nan
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad("white")
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.35))
    mesh = None
    for index, (ax, candidate) in enumerate(zip(axes, candidates)):
        mesh = ax.pcolormesh(data["u"], data["v"],
                             np.ma.masked_invalid(schlieren.T),
                             shading="nearest", cmap=cmap,
                             norm=Normalize(0.0, 1.0), rasterized=True)
        draw_solid(ax, data["u"], data["v"], data["solid"])
        for branch, color in ((upper, "#d62728"), (lower, "#1f77b4")):
            selected = branch["z_abs"] <= candidate["half_height_Db"]
            ax.plot(branch["x"][selected], branch["z_signed"][selected],
                    color=color, linewidth=1.0, linestyle="--", zorder=6)
        ax.axvline(candidate["display_left_X"], color="#2ca02c",
                   linewidth=0.9, linestyle=":", zorder=7)
        ax.set_xlim(candidate["display_left_X"], candidate["display_right_X"])
        height = candidate["half_height_Db"]
        ax.set_ylim(-height, height)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$X=(x-x_0)/D_b$")
        if index == 0:
            ax.set_ylabel(r"$Z=(z-z_0)/D_b$")
        else:
            ax.tick_params(labelleft=False)
        ax.set_title(
            r"$|Z|\leq{:.2f}$; upstream margin ${:.2f}D_b$".format(
                height, candidate["actual_upstream_margin_Db"]), fontsize=9.0)
        ax.text(0.02, 0.98, "({})".format(chr(97 + index)),
                transform=ax.transAxes, ha="left", va="top",
                fontsize=9.0, fontstyle="italic",
                bbox={"facecolor": "white", "edgecolor": "none",
                      "alpha": 0.75, "pad": 0.5})
    colorbar = fig.colorbar(mesh, ax=axes, fraction=0.024, pad=0.02,
                            ticks=(0.0, 0.5, 1.0))
    colorbar.set_label(r"Numerical schlieren, $S$")
    fig.suptitle(
        "Bow-shock window QA (dashed: tracked main gradient ridge; "
        "green: proposed upstream frame)", fontsize=10.0, y=0.98)
    fig.subplots_adjust(left=0.06, right=0.91, bottom=0.15, top=0.85,
                        wspace=0.13)
    output = output_dir / "nozzle_overview_window_candidates.png"
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    return output


def main():
    args = parse_args()
    if not (args.upstream_margin > 0.0 and args.schlieren_gref > 0.0
            and args.schlieren_k > 0.0):
        raise ValueError("Margin and schlieren parameters must be positive")
    half_heights = sorted(float(value) for value in args.half_heights)
    data = load_cache(Path(args.cache))
    upper = trace_branch(data["u"], data["v"], data["grad"], +1.0)
    lower = trace_branch(data["u"], data["v"], data["grad"], -1.0)
    candidates = candidate_windows(data["u"], upper, lower, half_heights,
                                   args.upstream_margin, args.right_x_m)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    montage = make_montage(data, upper, lower, candidates, args, output_dir)
    payload = {
        "purpose": "quantitative QA only; ridge overlays are excluded from final figure",
        "source_cache": str(Path(args.cache).resolve()),
        "plotfile": data["plotfile"],
        "time_ms": 1.0e3 * data["time_s"],
        "ridge_method": {
            "quantity": "Gaussian-smoothed log instantaneous |grad(rho)|",
            "seed": "X=-1.45 at |Z|=0.30, selecting the body/jet shock branch",
            "tracking": "dynamic programming with local smoothness and weak outward monotonicity",
            "warning": "upstream-boundary-connected startup wave deliberately excluded by seed",
        },
        "target_upstream_margin_Db": float(args.upstream_margin),
        "candidates": candidates,
        "upper_ridge": {
            "z_abs_Db": upper["z_abs"].tolist(),
            "x_over_Db": upper["x"].tolist(),
            "grad_rho_kg_m-4": upper["strength"].tolist(),
        },
        "lower_ridge": {
            "z_abs_Db": lower["z_abs"].tolist(),
            "x_over_Db": lower["x"].tolist(),
            "grad_rho_kg_m-4": lower["strength"].tolist(),
        },
        "montage": montage.name,
    }
    json_path = output_dir / "nozzle_overview_window_candidates.json"
    with json_path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    print(montage)
    print(json_path)


if __name__ == "__main__":
    main()
