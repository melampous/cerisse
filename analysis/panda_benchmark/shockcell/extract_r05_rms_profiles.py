#!/usr/bin/env python3
"""Extract off-axis RMS profiles at r/D_e = 0.5 from AMReX plotfiles.

The two-dimensional profiles are linearly interpolated in radius to the
orifice lip.  For the three-dimensional cases, the local temporal variance is
interpolated to a 256-point circle at r/D_e = 0.5 and then averaged in the
azimuthal direction.  AMR data are composited with the finest available level
at each axial and angular location.

This script is intended to run on aws-hpc, where the original plotfiles are
available.  It writes compact NPZ profiles to /tmp/codex_r05_profiles.
"""
import glob
import os
import re
from pathlib import Path

import numpy as np


D = 0.0254
R_SAMPLE = 0.5 * D
XMAX_D = 8.2
NTHETA = 256
OUT = Path("/tmp/codex_r05_profiles")

B2 = re.compile(r"\(\((-?\d+),(-?\d+)\) \((-?\d+),(-?\d+)\)")
B3 = re.compile(
    r"\(\((-?\d+),(-?\d+),(-?\d+)\) "
    r"\((-?\d+),(-?\d+),(-?\d+)\)"
)

VARS2D = ["pressureMEAN", "pressureSQR", "y_velocityMEAN", "y_velocitySQR"]
VARS3D = ["pressureMEAN", "pressureSQR", "x_velocityMEAN", "x_velocitySQR"]


def read_header(plotfile, pattern):
    header = (Path(plotfile) / "Header").read_text().split("\n")
    i = 1
    nvar = int(header[i])
    i += 1
    names = header[i:i + nvar]
    i += nvar
    i += 1
    time = float(header[i])
    i += 1
    finest = int(header[i])
    i += 1
    prob_lo = [float(v) for v in header[i].split()]
    i += 1
    prob_hi = [float(v) for v in header[i].split()]
    i += 1
    i += 1
    match = pattern.search(header[i])
    if match is None:
        raise RuntimeError(f"Cannot parse level-zero box in {plotfile}/Header")
    return names, nvar, time, finest, prob_lo, prob_hi, [int(v) for v in match.groups()]


def level_boxes(plotfile, level, pattern):
    cell_header = (Path(plotfile) / f"Level_{level}" / "Cell_H").read_text().split("\n")
    boxes = [
        [int(v) for v in match.groups()]
        for match in (pattern.search(line) for line in cell_header)
        if match is not None
    ]
    fabs = [
        (line.split()[1], int(line.split()[2]))
        for line in cell_header
        if line.strip().startswith("FabOnDisk:")
    ]
    if len(boxes) > len(fabs):
        boxes = boxes[len(boxes) - len(fabs):]
    return boxes, fabs


def linear_cells(coordinate, lower, spacing, ncells):
    """Return bracketing cell-centre indices and linear weights."""
    s = (coordinate - lower) / spacing - 0.5
    i0 = int(np.floor(s))
    fraction = s - i0
    pairs = [(i0, 1.0 - fraction), (i0 + 1, fraction)]
    return [(i, w) for i, w in pairs if 0 <= i < ncells and w > 1.0e-12]


def final_plotfile(case_directory):
    candidates = sorted(
        path for path in glob.glob(str(Path(case_directory) / "plt?????*"))
        if os.path.isdir(path) and not path.endswith(".temp")
    )
    if not candidates:
        raise RuntimeError(f"No complete plotfile found in {case_directory}")
    return candidates[-1]


def extract_2d(plotfile):
    names, _, time, finest, prob_lo, prob_hi, level0 = read_header(plotfile, B2)
    indices = [names.index(name) for name in VARS2D]
    nr0 = level0[2] + 1
    nx0 = level0[3] + 1
    dx_finest = (prob_hi[1] - prob_lo[1]) / (nx0 * 2**finest)
    nx = int(round(XMAX_D * D / dx_finest))
    pressure_variance = np.full(nx, np.nan)
    velocity_variance = np.full(nx, np.nan)

    for level in range(finest + 1):
        ratio = 2**(finest - level)
        nx_level = -(-nx // ratio)
        dr = (prob_hi[0] - prob_lo[0]) / (nr0 * 2**level)
        radial_cells = linear_cells(R_SAMPLE, prob_lo[0], dr, nr0 * 2**level)
        lines = {cell: np.full((4, nx_level), np.nan) for cell, _ in radial_cells}

        boxes, fabs = level_boxes(plotfile, level, B2)
        handles = {}
        for (ilo, jlo, ihi, jhi), (filename, offset) in zip(boxes, fabs):
            wanted = [cell for cell, _ in radial_cells if ilo <= cell <= ihi]
            if not wanted or jlo >= nx_level:
                continue
            ni, nj = ihi - ilo + 1, jhi - jlo + 1
            handle = handles.get(filename)
            if handle is None:
                handle = open(Path(plotfile) / f"Level_{level}" / filename, "rb")
                handles[filename] = handle
            handle.seek(offset)
            handle.readline()
            base = handle.tell()
            x1 = min(jhi + 1, nx_level)
            for variable_index, component in enumerate(indices):
                handle.seek(base + component * ni * nj * 8)
                array = np.fromfile(handle, dtype=np.float64, count=ni * nj).reshape(nj, ni)
                for cell in wanted:
                    lines[cell][variable_index, jlo:x1] = array[:x1 - jlo, cell - ilo]
        for handle in handles.values():
            handle.close()

        pvar_level = np.zeros(nx_level)
        uvar_level = np.zeros(nx_level)
        valid_p = np.ones(nx_level, dtype=bool)
        valid_u = np.ones(nx_level, dtype=bool)
        for cell, weight in radial_cells:
            line = lines[cell]
            pvar_cell = line[1] - line[0] ** 2
            uvar_cell = line[3] - line[2] ** 2
            valid_p &= np.isfinite(pvar_cell)
            valid_u &= np.isfinite(uvar_cell)
            pvar_level += weight * np.maximum(pvar_cell, 0.0)
            uvar_level += weight * np.maximum(uvar_cell, 0.0)

        p_fine = np.repeat(pvar_level, ratio)[:nx]
        u_fine = np.repeat(uvar_level, ratio)[:nx]
        valid_p_fine = np.repeat(valid_p, ratio)[:nx]
        valid_u_fine = np.repeat(valid_u, ratio)[:nx]
        pressure_variance[valid_p_fine] = p_fine[valid_p_fine]
        velocity_variance[valid_u_fine] = u_fine[valid_u_fine]

    return pressure_variance, velocity_variance, time, dx_finest


def extract_3d(plotfile):
    names, _, time, finest, prob_lo, prob_hi, level0 = read_header(plotfile, B3)
    indices = [names.index(name) for name in VARS3D]
    nx0, ny0, nz0 = level0[3] + 1, level0[4] + 1, level0[5] + 1
    dx_finest = (prob_hi[0] - prob_lo[0]) / (nx0 * 2**finest)
    nx = int(round(XMAX_D * D / dx_finest))
    pvar_theta = np.full((NTHETA, nx), np.nan)
    uvar_theta = np.full((NTHETA, nx), np.nan)
    theta = 2.0 * np.pi * np.arange(NTHETA) / NTHETA

    for level in range(finest + 1):
        ratio = 2**(finest - level)
        nx_level = -(-nx // ratio)
        ny = ny0 * 2**level
        nz = nz0 * 2**level
        dy = (prob_hi[1] - prob_lo[1]) / ny
        dz = (prob_hi[2] - prob_lo[2]) / nz

        samples = []
        needed = set()
        for angle in theta:
            y = R_SAMPLE * np.cos(angle)
            z = R_SAMPLE * np.sin(angle)
            yc = linear_cells(y, prob_lo[1], dy, ny)
            zc = linear_cells(z, prob_lo[2], dz, nz)
            weights = [((j, k), wy * wz) for j, wy in yc for k, wz in zc]
            samples.append(weights)
            needed.update(cell for cell, _ in weights)

        lines = {cell: np.full((4, nx_level), np.nan) for cell in needed}
        boxes, fabs = level_boxes(plotfile, level, B3)
        handles = {}
        for (ilo, jlo, klo, ihi, jhi, khi), (filename, offset) in zip(boxes, fabs):
            wanted = [cell for cell in needed if jlo <= cell[0] <= jhi and klo <= cell[1] <= khi]
            if not wanted or ilo >= nx_level:
                continue
            ni, nj, nk = ihi - ilo + 1, jhi - jlo + 1, khi - klo + 1
            handle = handles.get(filename)
            if handle is None:
                handle = open(Path(plotfile) / f"Level_{level}" / filename, "rb")
                handles[filename] = handle
            handle.seek(offset)
            handle.readline()
            base = handle.tell()
            x1 = min(ihi + 1, nx_level)
            for cell in wanted:
                j, k = cell
                for variable_index, component in enumerate(indices):
                    position = base + ((component * nk + (k - klo)) * nj + (j - jlo)) * ni * 8
                    handle.seek(position)
                    array = np.fromfile(handle, dtype=np.float64, count=ni)
                    lines[cell][variable_index, ilo:x1] = array[:x1 - ilo]
        for handle in handles.values():
            handle.close()

        for itheta, weights in enumerate(samples):
            pvar_level = np.zeros(nx_level)
            uvar_level = np.zeros(nx_level)
            valid_p = np.ones(nx_level, dtype=bool)
            valid_u = np.ones(nx_level, dtype=bool)
            for cell, weight in weights:
                line = lines[cell]
                pvar_cell = line[1] - line[0] ** 2
                uvar_cell = line[3] - line[2] ** 2
                valid_p &= np.isfinite(pvar_cell)
                valid_u &= np.isfinite(uvar_cell)
                pvar_level += weight * np.maximum(pvar_cell, 0.0)
                uvar_level += weight * np.maximum(uvar_cell, 0.0)
            p_fine = np.repeat(pvar_level, ratio)[:nx]
            u_fine = np.repeat(uvar_level, ratio)[:nx]
            valid_p_fine = np.repeat(valid_p, ratio)[:nx]
            valid_u_fine = np.repeat(valid_u, ratio)[:nx]
            pvar_theta[itheta, valid_p_fine] = p_fine[valid_p_fine]
            uvar_theta[itheta, valid_u_fine] = u_fine[valid_u_fine]

    pressure_variance = np.nanmean(pvar_theta, axis=0)
    velocity_variance = np.nanmean(uvar_theta, axis=0)
    return pressure_variance, velocity_variance, time, dx_finest


CASES = [
    ("2d", "0", "/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L012/stat_L0"),
    ("2d", "1", "/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L012/stat_L1"),
    ("2d", "2", "/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L012/stat_L2"),
    ("2d", "3", "/shared/cerisse/exm/underexpanded_jet/2d/panda_m142_L3start/statT"),
    ("2d", "4", "/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L4/statT"),
    ("2d", "5", "/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L5/statT"),
    ("2d", "6", "/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L6stat/statT"),
    ("3d", "1", "/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_L1_sens/stat"),
    ("3d", "2", "/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_L2_sens/stat"),
    ("3d", "3", "/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_L3_sens/stat"),
    ("3d", "4", "/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_L4_sens/stat"),
    ("3d", "5", "/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_3d_gpu/l5_prod_cb_20260715"),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for mode, tag, directory in CASES:
        plotfile = final_plotfile(directory)
        extractor = extract_2d if mode == "2d" else extract_3d
        pvar, uvar, time, dx = extractor(plotfile)
        x = (np.arange(len(pvar)) + 0.5) * dx / D
        output = OUT / f"r05_{mode}_L{tag}.npz"
        np.savez_compressed(
            output,
            x=x,
            pressure_variance=pvar,
            axial_velocity_variance=uvar,
            time=time,
            dx=dx,
            radius_over_D=0.5,
            ntheta=NTHETA if mode == "3d" else 1,
            source_plotfile=plotfile,
        )
        coverage = np.mean(np.isfinite(pvar))
        print(
            f"{mode} L{tag}: {Path(plotfile).name}, t={time*1e3:.5f} ms, "
            f"D/dx={D/dx:.0f}, pressure coverage={coverage:.3f}, wrote {output}",
            flush=True,
        )


if __name__ == "__main__":
    main()
