#!/usr/bin/env python3
"""Extract the scaled azimuthal-vorticity diagnostic from native AMR data.

The dimensionless quantity is evaluated on each native AMR level before
coarse data are repeated onto the display raster:

    q* = (rho_j D_e^2 / U_j) omega_theta / (rho r).

This order matters because both ``r`` and the axis mask must use the native
cell spacing.  The 2D source is the axisymmetric L4 checkpoint at 8.00 ms.
The 3D source is the Cartesian L5 plotfile at 10.16 ms.  In 3D, fields from
the two cell-centred planes adjacent to z=0 are averaged before omega_z and
the local cylindrical projection are formed.
"""

from pathlib import Path
import re

import numpy as np


D_E = 0.0254
U_J = 414.2
RHO_J = 1.6413
GAMMA = 1.4
XMAX_D = 4.2
RMAX_D = 1.25

CHECKPOINT_2D = Path(
    "/shared/cerisse_reflux/exm/underexpanded_jet/2d/"
    "m142_2d_L4/statT/chk16237"
)
PLOTFILE_3D = Path(
    "/shared/cerisse_reflux/exm/underexpanded_jet/3d/"
    "m142_3d_gpu/l5_prod_cb_20260715/plt09170"
)
OUTPUT = Path("/tmp/omega_rhor_native.npz")

BOX_2D = re.compile(
    r"\(\((-?\d+),(-?\d+)\) \((-?\d+),(-?\d+)\)"
)
BOX_3D = re.compile(
    r"\(\((-?\d+),(-?\d+),(-?\d+)\) "
    r"\((-?\d+),(-?\d+),(-?\d+)\)"
)


def dilate(mask):
    expanded = mask.copy()
    expanded[1:, :] |= mask[:-1, :]
    expanded[:-1, :] |= mask[1:, :]
    expanded[:, 1:] |= mask[:, :-1]
    expanded[:, :-1] |= mask[:, 1:]
    return expanded


def header_boxes(header_path, pattern):
    lines = header_path.read_text().splitlines()
    boxes = [
        tuple(int(value) for value in match.groups())
        for match in (pattern.search(line) for line in lines)
        if match is not None
    ]
    fabs = [
        (line.split()[1], int(line.split()[2]))
        for line in lines
        if line.strip().startswith("FabOnDisk:")
    ]
    if len(boxes) > len(fabs):
        boxes = boxes[-len(fabs):]
    if len(boxes) != len(fabs):
        raise RuntimeError(
            f"{header_path}: found {len(boxes)} boxes and {len(fabs)} FABs"
        )
    return boxes, fabs


def checkpoint_metadata():
    lines = (CHECKPOINT_2D / "Header").read_text().splitlines()
    time = float(lines[2])
    finest = int(lines[3])
    match = BOX_2D.search(lines[6])
    if match is None:
        raise RuntimeError("Cannot parse axisymmetric checkpoint domain")
    n0 = [int(match.group(3)) + 1, int(match.group(4)) + 1]
    prob_match = re.search(
        r"RealBox\s+([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s+"
        r"([+\-0-9.eE]+)\s+([+\-0-9.eE]+)",
        lines[6],
    )
    if prob_match is None:
        raise RuntimeError("Cannot parse axisymmetric physical domain")
    values = [float(value) for value in prob_match.groups()]
    prob_lo = np.array([values[0], values[2]])
    prob_hi = np.array([values[1], values[3]])
    return time, finest, n0, prob_lo, prob_hi


def read_checkpoint_state(level_dir, box, fab):
    """Read the five conservative components over one valid 2D box."""
    filename, offset = fab
    with (level_dir / filename).open("rb") as handle:
        handle.seek(offset)
        line = handle.readline().decode("ascii")
        ghost_match = BOX_2D.search(line)
        if ghost_match is None:
            raise RuntimeError(f"Cannot parse FAB header in {filename}")
        gi0, gj0, gi1, gj1 = (
            int(value) for value in ghost_match.groups()
        )
        ni = gi1 - gi0 + 1
        nj = gj1 - gj0 + 1
        count = ni * nj
        base = handle.tell()
        components = []
        for component in range(5):
            handle.seek(base + component * count * 8)
            values = np.fromfile(handle, dtype=np.float64, count=count)
            components.append(values.reshape(nj, ni))
    i0, j0, i1, j1 = box
    radial = slice(i0 - gi0, i1 - gi0 + 1)
    axial = slice(j0 - gj0, j1 - gj0 + 1)
    return np.stack(
        [component[axial, radial].T for component in components],
        axis=0,
    )


def even_radial_gradient(field, spacing):
    """Second-order radial derivative with even symmetry at the axis."""
    result = np.empty_like(field)
    result[0] = (field[1] - field[0]) / (2.0 * spacing)
    result[1:-1] = (field[2:] - field[:-2]) / (2.0 * spacing)
    result[-1] = (field[-1] - field[-2]) / spacing
    return result


def place_native(field, level, finest, target, level_map):
    refinement = 2 ** (finest - level)
    repeated = field
    if refinement > 1:
        repeated = np.repeat(
            np.repeat(repeated, refinement, axis=0),
            refinement,
            axis=1,
        )
    repeated = repeated[:target.shape[0], :target.shape[1]]
    good = np.isfinite(repeated)
    target_view = target[:repeated.shape[0], :repeated.shape[1]]
    level_view = level_map[:repeated.shape[0], :repeated.shape[1]]
    target_view[good] = repeated[good]
    level_view[good] = level


def extract_axisymmetric():
    time, finest, n0, prob_lo, prob_hi = checkpoint_metadata()
    finest_spacing = (prob_hi[0] - prob_lo[0]) / (n0[0] * 2**finest)
    nr = int(np.ceil(RMAX_D * D_E / finest_spacing))
    nx = int(np.ceil(XMAX_D * D_E / finest_spacing))
    q_composite = np.full((nr, nx), np.nan, dtype=np.float64)
    omega_composite = np.full_like(q_composite, np.nan)
    rho_composite = np.full_like(q_composite, np.nan)
    level_map = np.full((nr, nx), -1, dtype=np.int8)

    for level in range(finest + 1):
        refinement = 2 ** (finest - level)
        spacing = finest_spacing * refinement
        nr_level = int(np.ceil(RMAX_D * D_E / spacing))
        nx_level = int(np.ceil(XMAX_D * D_E / spacing))
        state = np.full((5, nr_level, nx_level), np.nan)
        level_dir = CHECKPOINT_2D / f"Level_{level}"
        boxes, fabs = header_boxes(
            level_dir / "SD_0_New_MF_H",
            BOX_2D,
        )
        for box, fab in zip(boxes, fabs):
            i0, j0, i1, j1 = box
            if i0 >= nr_level or j0 >= nx_level:
                continue
            values = read_checkpoint_state(level_dir, box, fab)
            i_end = min(i1 + 1, nr_level)
            j_end = min(j1 + 1, nx_level)
            state[:, i0:i_end, j0:j_end] = values[
                :,
                :i_end - i0,
                :j_end - j0,
            ]

        mx, mr, _, energy, rho = state
        ur = mx / rho
        ux = mr / rho
        pressure = (GAMMA - 1.0) * (
            energy - 0.5 * (mx**2 + mr**2) / rho
        )
        omega_theta = (
            np.gradient(ur, spacing, axis=1)
            - even_radial_gradient(ux, spacing)
        )
        radius = (np.arange(nr_level) + 0.5) * spacing
        q_scaled = (
            RHO_J * D_E**2 * omega_theta
            / (rho * radius[:, None] * U_J)
        )
        bad = dilate(
            ~np.isfinite(ur + ux + rho + pressure)
        )
        omega_theta[bad] = np.nan
        rho = rho.copy()
        rho[bad] = np.nan
        q_scaled[bad] = np.nan
        place_native(
            omega_theta,
            level,
            finest,
            omega_composite,
            level_map,
        )
        place_native(rho, level, finest, rho_composite, level_map)
        place_native(q_scaled, level, finest, q_composite, level_map)
        print(f"2D level {level} done", flush=True)

    r = (np.arange(nr) + 0.5) * finest_spacing / D_E
    x = (np.arange(nx) + 0.5) * finest_spacing / D_E
    return {
        "r2": r,
        "x2": x,
        "q2": q_composite,
        "omega2": omega_composite,
        "rho2": rho_composite,
        "level2": level_map,
        "time2": time,
        "dx2": finest_spacing,
    }


def plotfile_metadata():
    lines = (PLOTFILE_3D / "Header").read_text().splitlines()
    index = 1
    nvar = int(lines[index])
    index += 1
    names = lines[index:index + nvar]
    index += nvar + 1
    time = float(lines[index])
    index += 1
    finest = int(lines[index])
    index += 1
    prob_lo = np.array([float(value) for value in lines[index].split()])
    index += 1
    prob_hi = np.array([float(value) for value in lines[index].split()])
    index += 2
    match = BOX_3D.search(lines[index])
    if match is None:
        raise RuntimeError("Cannot parse Cartesian plotfile domain")
    n0 = [
        int(match.group(4)) + 1,
        int(match.group(5)) + 1,
        int(match.group(6)) + 1,
    ]
    variables = ["x_velocity", "y_velocity", "Density"]
    components = [names.index(variable) for variable in variables]
    return time, finest, n0, prob_lo, prob_hi, components


def extract_cartesian():
    time, finest, n0, prob_lo, prob_hi, components = plotfile_metadata()
    finest_spacing = (prob_hi[0] - prob_lo[0]) / (n0[0] * 2**finest)
    nx = int(np.ceil(XMAX_D * D_E / finest_spacing))
    half_width = int(np.ceil(RMAX_D * D_E / finest_spacing))
    ny = 2 * half_width
    q_composite = np.full((nx, ny), np.nan, dtype=np.float64)
    omega_composite = np.full_like(q_composite, np.nan)
    rho_composite = np.full_like(q_composite, np.nan)
    level_map = np.full((nx, ny), -1, dtype=np.int8)

    for level in range(finest + 1):
        refinement = 2 ** (finest - level)
        spacing = finest_spacing * refinement
        nx_level = int(np.ceil(XMAX_D * D_E / spacing))
        half_level = int(np.ceil(RMAX_D * D_E / spacing))
        y_centre = (n0[1] * 2**level) // 2
        z_centre = (n0[2] * 2**level) // 2
        y0, y1 = y_centre - half_level, y_centre + half_level
        z_indices = [z_centre - 1, z_centre]
        arrays = {
            name: np.full((2, nx_level, 2 * half_level), np.nan)
            for name in ("ux", "uy", "rho")
        }
        level_dir = PLOTFILE_3D / f"Level_{level}"
        boxes, fabs = header_boxes(level_dir / "Cell_H", BOX_3D)
        handles = {}
        for box, (filename, offset) in zip(boxes, fabs):
            i0, j0, k0, i1, j1, k1 = box
            if (
                i0 >= nx_level
                or j1 < y0
                or j0 >= y1
                or k1 < z_indices[0]
                or k0 > z_indices[1]
            ):
                continue
            ni, nj, nk = i1 - i0 + 1, j1 - j0 + 1, k1 - k0 + 1
            handle = handles.get(filename)
            if handle is None:
                handle = (level_dir / filename).open("rb")
                handles[filename] = handle
            handle.seek(offset)
            handle.readline()
            base = handle.tell()
            x_end = min(i1 + 1, nx_level)
            j_low, j_high = max(j0, y0), min(j1 + 1, y1)
            for plane_index, global_k in enumerate(z_indices):
                if global_k < k0 or global_k > k1:
                    continue
                for name, component in zip(
                    ("ux", "uy", "rho"),
                    components,
                ):
                    handle.seek(
                        base
                        + (
                            component * nk + global_k - k0
                        ) * nj * ni * 8
                    )
                    values = np.fromfile(
                        handle,
                        dtype=np.float64,
                        count=nj * ni,
                    ).reshape(nj, ni)
                    arrays[name][
                        plane_index,
                        i0:x_end,
                        j_low - y0:j_high - y0,
                    ] = values[
                        j_low - j0:j_high - j0,
                        :x_end - i0,
                    ].T
        for handle in handles.values():
            handle.close()

        omega_planes = []
        for plane in range(2):
            d_uy_dx = np.gradient(
                arrays["uy"][plane],
                spacing,
                axis=0,
            )
            d_ux_dy = np.gradient(
                arrays["ux"][plane],
                spacing,
                axis=1,
            )
            omega_planes.append(d_uy_dx - d_ux_dy)
        omega_theta = 0.5 * (omega_planes[0] + omega_planes[1])
        rho = 0.5 * (arrays["rho"][0] + arrays["rho"][1])
        y = (np.arange(-half_level, half_level) + 0.5) * spacing
        radius = np.abs(y)
        omega_theta *= np.sign(y)[None, :]
        q_scaled = (
            RHO_J * D_E**2 * omega_theta
            / (rho * radius[None, :] * U_J)
        )
        bad = dilate(
            ~np.isfinite(
                arrays["ux"][0]
                + arrays["ux"][1]
                + arrays["uy"][0]
                + arrays["uy"][1]
                + rho
            )
        )
        omega_theta[bad] = np.nan
        rho[bad] = np.nan
        q_scaled[bad] = np.nan
        place_native(
            omega_theta,
            level,
            finest,
            omega_composite,
            level_map,
        )
        place_native(rho, level, finest, rho_composite, level_map)
        place_native(q_scaled, level, finest, q_composite, level_map)
        print(f"3D level {level} done", flush=True)

    x = (np.arange(nx) + 0.5) * finest_spacing / D_E
    y = (
        np.arange(-half_width, half_width) + 0.5
    ) * finest_spacing / D_E
    return {
        "x3": x,
        "y3": y,
        "q3": q_composite,
        "omega3": omega_composite,
        "rho3": rho_composite,
        "level3": level_map,
        "time3": time,
        "dx3": finest_spacing,
    }


def main():
    data = {}
    data.update(extract_axisymmetric())
    data.update(extract_cartesian())
    data.update(
        {
            "D_e": D_E,
            "U_j": U_J,
            "rho_j": RHO_J,
            "source2": str(CHECKPOINT_2D),
            "source3": str(PLOTFILE_3D),
        }
    )
    np.savez_compressed(
        OUTPUT,
        **{
            key: value.astype(np.float32)
            if isinstance(value, np.ndarray)
            and np.issubdtype(value.dtype, np.floating)
            else value
            for key, value in data.items()
        },
    )
    print(
        f"WROTE {OUTPUT}: q2={data['q2'].shape}, "
        f"q3={data['q3'].shape}",
        flush=True,
    )


if __name__ == "__main__":
    main()
