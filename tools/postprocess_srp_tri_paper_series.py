#!/usr/bin/env python3
"""Paper-grade three-view and time-series post-processing for SRP tri-nozzle.

The script reads only selected AMReX slices.  It deliberately keeps all
display limits and the numerical-schlieren reference fixed across time.
Three modes make it practical to run on CX3 without repeatedly scanning the
large plotfiles:

  calibrate  read/cache the latest state and determine the global G_ref;
  extract    read/cache/render a strided subset of the plotfiles;
  assemble   build time montages, the data archive, and the manifest;
  all        perform all three stages sequentially.

The plotted planes are interpolated to their requested physical locations.
For schlieren, the full three-dimensional density-gradient magnitude is first
computed at the bracketing cell-centre planes and then interpolated.
"""

import argparse
import copy
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm, Normalize
from matplotlib.ticker import LogFormatterMathtext, MaxNLocator
import numpy as np


GAMMA = 1.4
P_INF = 574.56
BODY_DIAMETER = 0.127
BODY_CENTRE = (0.300, 0.275, 0.275)
PLANE_LOCATIONS = {"xy": 0.275, "xz": 0.259, "yz": 0.200}
AXIAL_CROP = (0.0, 0.56)
TRANSVERSE_CROP = (0.08, 0.47)
CACHE_SCHEMA_VERSION = 3
SAMPLER_DESCRIPTION = "yt smoothed_covering_grid for primitives; covering_grid for IBM IDs"
VIEW_ORDER = ("xy", "xz", "yz")
VARIABLES = (
    "Density", "pressure", "x_velocity", "y_velocity", "z_velocity",
    "sld", "ghs",
)
FIELD_ORDER = ("mach", "pressure_ratio", "schlieren")
EARLY_STEPS = (0, 500, 1000, 1500, 2000, 2500)
REPRESENTATIVE_STEPS = (0, 1500, 3000, 4500, 6000, 7500, 9000)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("calibrate", "extract", "assemble", "all"),
                        default="all")
    parser.add_argument("--plot-root", required=True)
    parser.add_argument("--reader", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--worker-index", type=int, default=0)
    parser.add_argument("--worker-count", type=int, default=1)
    parser.add_argument("--paper-dpi", type=int, default=600)
    parser.add_argument("--frame-dpi", type=int, default=400)
    parser.add_argument("--schlieren-k", type=float, default=4.0)
    parser.add_argument("--schlieren-percentile", type=float, default=99.5)
    parser.add_argument("--fixed-gref", type=float, default=None,
                        help="Use this fixed schlieren G_ref instead of recalibrating")
    parser.add_argument("--axial-crop", type=float, nargs=2, default=AXIAL_CROP,
                        metavar=("X_MIN", "X_MAX"))
    parser.add_argument("--transverse-crop", type=float, nargs=2,
                        default=TRANSVERSE_CROP, metavar=("R_MIN", "R_MAX"))
    parser.add_argument("--overwrite-cache", action="store_true")
    return parser.parse_args()


def configure_crops(args):
    global AXIAL_CROP, TRANSVERSE_CROP
    axial = tuple(float(value) for value in args.axial_crop)
    transverse = tuple(float(value) for value in args.transverse_crop)
    if not (axial[0] < axial[1] and transverse[0] < transverse[1]):
        raise ValueError("Crop lower bounds must be smaller than upper bounds")
    AXIAL_CROP = axial
    TRANSVERSE_CROP = transverse


def configure_matplotlib():
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9.0,
        "axes.labelsize": 9.0,
        "axes.titlesize": 8.0,
        "axes.linewidth": 0.6,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def load_reader(path):
    spec = importlib.util.spec_from_file_location("amrex_slice_reader", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load AMReX reader: {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_yt_module():
    try:
        import yt
    except ImportError as error:
        raise RuntimeError(
            "yt is required for ghost-aware AMR sampling; load the CX3 yt/4.3 module"
        ) from error
    yt.funcs.mylog.setLevel(30)
    return yt


def discover_plotfiles(root):
    found = []
    for path in Path(root).glob("plt[0-9]*"):
        match = re.fullmatch(r"plt(\d+)", path.name)
        if path.is_dir() and match and (path / "Header").is_file():
            found.append((int(match.group(1)), str(path)))
    found.sort()
    if not found:
        raise RuntimeError("No complete pltNNNNN directories found under {}".format(root))
    return found


def bracket_stencil(coords, value):
    """Return bracketing indices, interpolation weight, and gradient stencil."""
    hi = int(np.searchsorted(coords, value, side="left"))
    hi = max(1, min(hi, len(coords) - 1))
    lo = hi - 1
    weight = float((value - coords[lo]) / (coords[hi] - coords[lo]))
    weight = max(0.0, min(1.0, weight))
    indices = sorted(set(max(0, min(len(coords) - 1, i))
                         for i in (lo - 1, lo, hi, hi + 1)))
    return lo, hi, weight, indices


def coordinates_from_header(header):
    level = int(header["finest_level"])
    box = header["all_boxes"][level]
    nx = int(box[3]) - int(box[0]) + 1
    ny = int(box[4]) - int(box[1]) + 1
    nz = int(box[5]) - int(box[2]) + 1
    lo = np.asarray(header["lo"], dtype=float)
    hi = np.asarray(header["hi"], dtype=float)
    spacing = (hi - lo) / np.asarray([nx, ny, nz], dtype=float)
    coords = tuple(lo[d] + (np.arange((nx, ny, nz)[d]) + 0.5) * spacing[d]
                   for d in range(3))
    return coords, tuple(float(v) for v in spacing)


def interp_planes(planes, lo, hi, weight):
    return (1.0 - weight) * planes[lo] + weight * planes[hi]


def plane_valid(rho, sld, ghs):
    return (np.isfinite(rho) & np.isfinite(sld) & np.isfinite(ghs) &
            (rho > 0.0) & (sld < 0.5) & (ghs < 0.5))


def masked_derivative(field, valid, spacing, axis):
    """Centred derivative with one-sided fallback, never crossing invalid cells."""
    forward = np.full(field.shape, np.nan, dtype=np.float64)
    backward = np.full(field.shape, np.nan, dtype=np.float64)
    forward_valid = np.zeros(field.shape, dtype=bool)
    backward_valid = np.zeros(field.shape, dtype=bool)
    if axis == 0:
        forward[:-1, :] = field[1:, :]
        backward[1:, :] = field[:-1, :]
        forward_valid[:-1, :] = valid[1:, :]
        backward_valid[1:, :] = valid[:-1, :]
    elif axis == 1:
        forward[:, :-1] = field[:, 1:]
        backward[:, 1:] = field[:, :-1]
        forward_valid[:, :-1] = valid[:, 1:]
        backward_valid[:, 1:] = valid[:, :-1]
    else:
        raise ValueError("2-D derivative axis must be 0 or 1")

    result = np.full(field.shape, np.nan, dtype=np.float64)
    central = valid & forward_valid & backward_valid
    result[central] = (forward[central] - backward[central]) / (2.0 * spacing)
    use_forward = valid & forward_valid & ~backward_valid
    result[use_forward] = (forward[use_forward] - field[use_forward]) / spacing
    use_backward = valid & backward_valid & ~forward_valid
    result[use_backward] = (field[use_backward] - backward[use_backward]) / spacing
    return result


def normal_derivative_at(planes, valid_planes, index, spacing):
    centre = planes[index]
    valid = valid_planes[index]
    prev_index = index - 1 if index - 1 in planes else None
    next_index = index + 1 if index + 1 in planes else None
    result = np.full(centre.shape, np.nan, dtype=np.float64)
    if prev_index is not None and next_index is not None:
        central = valid & valid_planes[prev_index] & valid_planes[next_index]
        result[central] = ((planes[next_index][central] - planes[prev_index][central]) /
                           (2.0 * spacing))
    else:
        central = np.zeros(centre.shape, dtype=bool)
    if next_index is not None:
        forward = valid & valid_planes[next_index] & ~central
        result[forward] = (planes[next_index][forward] - centre[forward]) / spacing
    if prev_index is not None:
        backward = valid & valid_planes[prev_index] & ~central & ~np.isfinite(result)
        result[backward] = (centre[backward] - planes[prev_index][backward]) / spacing
    return result


def gradient_components_at(planes, sld_planes, ghs_planes, index,
                           normal_spacing, inplane_spacings):
    valid_planes = {key: plane_valid(planes[key], sld_planes[key], ghs_planes[key])
                    for key in planes}
    valid = valid_planes[index]
    g0 = masked_derivative(planes[index], valid, inplane_spacings[0], axis=0)
    g1 = masked_derivative(planes[index], valid, inplane_spacings[1], axis=1)
    gn = normal_derivative_at(planes, valid_planes, index, normal_spacing)
    for component in (g0, g1, gn):
        component[~valid] = np.nan
    return g0, g1, gn


def crop_indices(coords, lower, upper):
    indices = np.flatnonzero((coords >= lower) & (coords <= upper))
    if indices.size < 2:
        raise RuntimeError("Requested crop [{}, {}] contains too few cells".format(lower, upper))
    return indices


def make_view(slices, direction, coords, spacings, location):
    """Interpolate primitives and full |grad rho| onto one requested view."""
    axis = {"yz": 0, "xz": 1, "xy": 2}[direction]
    normal_coords = coords[axis]
    lo, hi, weight, _ = bracket_stencil(normal_coords, location)
    interpolated = {}
    for name in VARIABLES:
        interpolated[name] = interp_planes(slices[name][axis], lo, hi, weight)

    rho_planes = slices["Density"][axis]
    sld_planes = slices["sld"][axis]
    ghs_planes = slices["ghs"][axis]
    inplane_axes = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[axis]
    components_lo = gradient_components_at(
        rho_planes, sld_planes, ghs_planes, lo, spacings[axis],
        (spacings[inplane_axes[0]], spacings[inplane_axes[1]]),
    )
    components_hi = gradient_components_at(
        rho_planes, sld_planes, ghs_planes, hi, spacings[axis],
        (spacings[inplane_axes[0]], spacings[inplane_axes[1]]),
    )
    # Interpolate signed gradient components before taking the magnitude.
    # This is essential on symmetry planes, where the normal component has
    # opposite signs on the two bracketing cell-centre planes.
    components = tuple((1.0 - weight) * lo_component + weight * hi_component
                       for lo_component, hi_component in zip(components_lo, components_hi))
    grad = np.sqrt(sum(component * component for component in components))

    rho = interpolated["Density"]
    pressure = interpolated["pressure"]
    sld = interpolated["sld"]
    ghs = interpolated["ghs"]
    valid = (plane_valid(rho, sld, ghs) & np.isfinite(pressure) & (pressure > 0.0) &
             np.isfinite(interpolated["x_velocity"]) &
             np.isfinite(interpolated["y_velocity"]) &
             np.isfinite(interpolated["z_velocity"]))
    speed = np.sqrt(interpolated["x_velocity"] ** 2 +
                    interpolated["y_velocity"] ** 2 +
                    interpolated["z_velocity"] ** 2)
    sound = np.sqrt(GAMMA * np.maximum(pressure, 0.0) /
                    np.maximum(rho, np.finfo(float).tiny))
    mach = speed / np.maximum(sound, np.finfo(float).tiny)
    pressure_ratio = pressure / P_INF
    solid = (sld >= 0.5) & np.isfinite(sld)
    for array in (mach, pressure_ratio, grad):
        array[~valid] = np.nan

    if direction == "xy":
        u_raw, v_raw = coords[0], coords[1]
        u_index = crop_indices(u_raw, AXIAL_CROP[0], AXIAL_CROP[1])
        v_index = crop_indices(v_raw, TRANSVERSE_CROP[0], TRANSVERSE_CROP[1])
        u = (u_raw[u_index] - BODY_CENTRE[0]) / BODY_DIAMETER
        v = (v_raw[v_index] - BODY_CENTRE[1]) / BODY_DIAMETER
    elif direction == "xz":
        u_raw, v_raw = coords[0], coords[2]
        u_index = crop_indices(u_raw, AXIAL_CROP[0], AXIAL_CROP[1])
        v_index = crop_indices(v_raw, TRANSVERSE_CROP[0], TRANSVERSE_CROP[1])
        u = (u_raw[u_index] - BODY_CENTRE[0]) / BODY_DIAMETER
        v = (v_raw[v_index] - BODY_CENTRE[2]) / BODY_DIAMETER
    else:
        u_raw, v_raw = coords[1], coords[2]
        u_index = crop_indices(u_raw, TRANSVERSE_CROP[0], TRANSVERSE_CROP[1])
        v_index = crop_indices(v_raw, TRANSVERSE_CROP[0], TRANSVERSE_CROP[1])
        u = (u_raw[u_index] - BODY_CENTRE[1]) / BODY_DIAMETER
        v = (v_raw[v_index] - BODY_CENTRE[2]) / BODY_DIAMETER

    selection = np.ix_(u_index, v_index)
    return {
        "u": u.astype(np.float64),
        "v": v.astype(np.float64),
        "mach": mach[selection].astype(np.float32),
        "pressure_ratio": pressure_ratio[selection].astype(np.float32),
        "grad_rho": grad[selection].astype(np.float32),
        "solid": solid[selection].astype(np.uint8),
        "source_level": np.zeros(mach[selection].shape, dtype=np.uint8),
        "location": float(location),
        "bracket_lo": int(lo),
        "bracket_hi": int(hi),
        "weight": float(weight),
    }


def process_plotfile(reader, plotfile):
    header = reader._parse_header_3d(plotfile)
    coords, spacings = coordinates_from_header(header)
    x_info = bracket_stencil(coords[0], PLANE_LOCATIONS["yz"])
    y_info = bracket_stencil(coords[1], PLANE_LOCATIONS["xz"])
    z_info = bracket_stencil(coords[2], PLANE_LOCATIONS["xy"])
    slices, xc, yc, zc, dx, dy, dz, meta = reader.read_amrex_3d_multi_slices_multi_ix(
        plotfile, VARIABLES,
        ix_list=x_info[3], iy_list=y_info[3], iz_list=z_info[3],
        out_level=-1, prolong_method="linear",
    )
    read_coords = (xc, yc, zc)
    read_spacings = (dx, dy, dz)
    views = {name: make_view(slices, name, read_coords, read_spacings,
                             PLANE_LOCATIONS[name]) for name in VIEW_ORDER}
    return {
        "time_s": float(meta["sim_time"]),
        "finest_level": int(meta["finest_level"]),
        "output_level": int(meta["output_level"]),
        "shape": (int(meta["Nx"]), int(meta["Ny"]), int(meta["Nz"])),
        "spacing": tuple(float(x) for x in read_spacings),
        "views": views,
    }


def _yt_array(values):
    """Return a unit-stripped NumPy view for yt arrays and plain arrays."""
    if hasattr(values, "to_value"):
        return np.asarray(values.to_value(), dtype=np.float64)
    if hasattr(values, "d"):
        return np.asarray(values.d, dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


def _inclusive_cell_range(coords, bounds):
    selected = np.flatnonzero((coords >= bounds[0]) & (coords <= bounds[1]))
    if selected.size < 2:
        raise RuntimeError("Crop {} contains too few cells".format(bounds))
    return int(selected[0]), int(selected[-1]) + 1


def _take_plane(array, axis, index):
    return np.take(array, int(index), axis=axis)


def _source_level_at(ds, u_coords, v_coords, u_axis, v_axis,
                     normal_axis, normal_position):
    """Map the native AMR source level at one cell-centre plane."""
    levels = np.zeros((len(u_coords), len(v_coords)), dtype=np.uint8)
    grids = sorted(ds.index.grids, key=lambda grid: int(grid.Level))
    tolerance = 1.0e-12 * float(np.max(_yt_array(ds.domain_width)))
    for grid in grids:
        left = _yt_array(grid.LeftEdge)
        right = _yt_array(grid.RightEdge)
        if not (left[normal_axis] - tolerance <= normal_position <
                right[normal_axis] - tolerance):
            continue
        i0 = int(np.searchsorted(u_coords, left[u_axis] - tolerance, side="left"))
        i1 = int(np.searchsorted(u_coords, right[u_axis] - tolerance, side="left"))
        j0 = int(np.searchsorted(v_coords, left[v_axis] - tolerance, side="left"))
        j1 = int(np.searchsorted(v_coords, right[v_axis] - tolerance, side="left"))
        i0, i1 = max(0, i0), min(len(u_coords), i1)
        j0, j1 = max(0, j0), min(len(v_coords), j1)
        if i0 < i1 and j0 < j1:
            levels[i0:i1, j0:j1] = np.maximum(
                levels[i0:i1, j0:j1], np.uint8(grid.Level))
    return levels


def _yt_view(ds, direction, fine_coords, fine_spacings):
    """Build one ghost-aware composite view using thin yt covering grids."""
    normal_axis = {"yz": 0, "xz": 1, "xy": 2}[direction]
    inplane_axes = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[normal_axis]
    location = PLANE_LOCATIONS[direction]
    lo, hi, weight, stencil = bracket_stencil(fine_coords[normal_axis], location)

    starts = np.zeros(3, dtype=np.int64)
    stops = np.zeros(3, dtype=np.int64)
    for axis in range(3):
        if axis == normal_axis:
            starts[axis] = int(stencil[0])
            stops[axis] = int(stencil[-1]) + 1
        else:
            bounds = AXIAL_CROP if axis == 0 else TRANSVERSE_CROP
            starts[axis], stops[axis] = _inclusive_cell_range(fine_coords[axis], bounds)
    # yt's smoothed covering grid expands by two cells at every recursive AMR
    # level.  Keep the level-0-equivalent halo inside non-periodic domains;
    # this removes only 6.67 mm (0.0525 D_b) at the x inflow for this L2 case
    # and avoids inventing periodic/ghost data there.
    fine_dims = np.asarray([len(values) for values in fine_coords], dtype=np.int64)
    interpolation_halo = 2 * int(ds.refine_by) ** int(ds.max_level)
    starts = np.maximum(starts, interpolation_halo)
    stops = np.minimum(stops, fine_dims - interpolation_halo)
    dims = stops - starts
    domain_left = _yt_array(ds.domain_left_edge)
    left_edge = domain_left + starts * np.asarray(fine_spacings)

    continuous_fields = [("boxlib", name) for name in VARIABLES[:5]]
    geometry_fields = [("boxlib", "sld"), ("boxlib", "ghs")]
    smooth = ds.smoothed_covering_grid(
        int(ds.max_level), left_edge=left_edge, dims=dims,
        fields=continuous_fields, num_ghost_zones=0, use_pbar=False,
    )
    categorical = ds.covering_grid(
        int(ds.max_level), left_edge=left_edge, dims=dims,
        fields=geometry_fields, num_ghost_zones=0, use_pbar=False,
    )
    arrays = {name: _yt_array(smooth[("boxlib", name)]) for name in VARIABLES[:5]}
    arrays["sld"] = _yt_array(categorical[("boxlib", "sld")])
    arrays["ghs"] = _yt_array(categorical[("boxlib", "ghs")])

    local_lo = int(lo - starts[normal_axis])
    local_hi = int(hi - starts[normal_axis])
    normal_indices = range(int(dims[normal_axis]))
    rho_planes = {index: _take_plane(arrays["Density"], normal_axis, index)
                  for index in normal_indices}
    sld_planes = {index: _take_plane(arrays["sld"], normal_axis, index)
                  for index in normal_indices}
    ghs_planes = {index: _take_plane(arrays["ghs"], normal_axis, index)
                  for index in normal_indices}
    components_lo = gradient_components_at(
        rho_planes, sld_planes, ghs_planes, local_lo,
        fine_spacings[normal_axis],
        (fine_spacings[inplane_axes[0]], fine_spacings[inplane_axes[1]]),
    )
    components_hi = gradient_components_at(
        rho_planes, sld_planes, ghs_planes, local_hi,
        fine_spacings[normal_axis],
        (fine_spacings[inplane_axes[0]], fine_spacings[inplane_axes[1]]),
    )
    components = tuple((1.0 - weight) * component_lo + weight * component_hi
                       for component_lo, component_hi in zip(components_lo, components_hi))
    grad = np.sqrt(sum(component * component for component in components))

    plane_values = {}
    for name, array in arrays.items():
        lower = _take_plane(array, normal_axis, local_lo)
        upper = _take_plane(array, normal_axis, local_hi)
        plane_values[name] = (lower, upper, (1.0 - weight) * lower + weight * upper)

    rho_lo, rho_hi, rho = plane_values["Density"]
    p_lo, p_hi, pressure = plane_values["pressure"]
    sld_lo, sld_hi, _ = plane_values["sld"]
    ghs_lo, ghs_hi, _ = plane_values["ghs"]
    fluid_lo = plane_valid(rho_lo, sld_lo, ghs_lo) & np.isfinite(p_lo) & (p_lo > 0.0)
    fluid_hi = plane_valid(rho_hi, sld_hi, ghs_hi) & np.isfinite(p_hi) & (p_hi > 0.0)
    valid = fluid_lo & fluid_hi
    for name in ("x_velocity", "y_velocity", "z_velocity"):
        lower, upper, _ = plane_values[name]
        valid &= np.isfinite(lower) & np.isfinite(upper)

    ux = plane_values["x_velocity"][2]
    uy = plane_values["y_velocity"][2]
    uz = plane_values["z_velocity"][2]
    speed = np.sqrt(ux * ux + uy * uy + uz * uz)
    sound = np.sqrt(GAMMA * np.maximum(pressure, 0.0) /
                    np.maximum(rho, np.finfo(float).tiny))
    mach = speed / np.maximum(sound, np.finfo(float).tiny)
    pressure_ratio = pressure / P_INF
    solid = ((sld_lo >= 0.5) | (sld_hi >= 0.5))
    for array in (mach, pressure_ratio, grad):
        array[~valid] = np.nan

    u_axis, v_axis = inplane_axes
    u_raw = fine_coords[u_axis][starts[u_axis]:stops[u_axis]]
    v_raw = fine_coords[v_axis][starts[v_axis]:stops[v_axis]]
    source_lo = _source_level_at(ds, u_raw, v_raw, u_axis, v_axis,
                                 normal_axis, fine_coords[normal_axis][lo])
    source_hi = _source_level_at(ds, u_raw, v_raw, u_axis, v_axis,
                                 normal_axis, fine_coords[normal_axis][hi])
    source_level = np.minimum(source_lo, source_hi)

    if direction == "xy":
        u = (u_raw - BODY_CENTRE[0]) / BODY_DIAMETER
        v = (v_raw - BODY_CENTRE[1]) / BODY_DIAMETER
    elif direction == "xz":
        u = (u_raw - BODY_CENTRE[0]) / BODY_DIAMETER
        v = (v_raw - BODY_CENTRE[2]) / BODY_DIAMETER
    else:
        u = (u_raw - BODY_CENTRE[1]) / BODY_DIAMETER
        v = (v_raw - BODY_CENTRE[2]) / BODY_DIAMETER

    return {
        "u": u.astype(np.float64),
        "v": v.astype(np.float64),
        "mach": mach.astype(np.float32),
        "pressure_ratio": pressure_ratio.astype(np.float32),
        "grad_rho": grad.astype(np.float32),
        "solid": solid.astype(np.uint8),
        "source_level": source_level.astype(np.uint8),
        "location": float(location),
        "bracket_lo": int(lo),
        "bracket_hi": int(hi),
        "weight": float(weight),
    }


def process_plotfile_yt(yt_module, plotfile):
    """Read three thin, globally interpolated AMR slabs with yt."""
    dataset = yt_module.load(plotfile)
    base_dims = np.asarray(dataset.domain_dimensions, dtype=np.int64)
    fine_dims = base_dims * int(dataset.refine_by) ** int(dataset.max_level)
    left = _yt_array(dataset.domain_left_edge)
    right = _yt_array(dataset.domain_right_edge)
    spacings = tuple(float(value) for value in (right - left) / fine_dims)
    fine_coords = tuple(left[axis] + (np.arange(fine_dims[axis]) + 0.5) * spacings[axis]
                        for axis in range(3))
    views = {name: _yt_view(dataset, name, fine_coords, spacings) for name in VIEW_ORDER}
    time_value = float(_yt_array(dataset.current_time))
    result = {
        "time_s": time_value,
        "finest_level": int(dataset.max_level),
        "output_level": int(dataset.max_level),
        "shape": tuple(int(value) for value in fine_dims),
        "spacing": spacings,
        "sampler": SAMPLER_DESCRIPTION,
        "views": views,
    }
    del dataset
    return result


def cache_path(output_dir, step):
    return Path(output_dir) / "cache" / ("plt{:05d}.npz".format(step))


def save_cache(path, step, plotfile, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_schema_version": np.asarray(CACHE_SCHEMA_VERSION, dtype=np.int32),
        "step": np.asarray(step, dtype=np.int64),
        "plotfile": np.asarray(plotfile),
        "sampler": np.asarray(result.get("sampler", "legacy selective-slice reader")),
        "time_s": np.asarray(result["time_s"], dtype=np.float64),
        "finest_level": np.asarray(result["finest_level"], dtype=np.int32),
        "output_level": np.asarray(result["output_level"], dtype=np.int32),
        "shape": np.asarray(result["shape"], dtype=np.int32),
        "spacing": np.asarray(result["spacing"], dtype=np.float64),
    }
    for view_name, view in result["views"].items():
        for key in ("u", "v", "mach", "pressure_ratio", "grad_rho", "solid",
                    "source_level"):
            payload[view_name + "_" + key] = view[key]
        payload[view_name + "_location"] = np.asarray(view["location"], dtype=np.float64)
        payload[view_name + "_bracket"] = np.asarray(
            [view["bracket_lo"], view["bracket_hi"], view["weight"]], dtype=np.float64)
    temporary = str(path) + ".tmp.npz"
    np.savez_compressed(temporary, **payload)
    os.replace(temporary, str(path))


def load_cache(path):
    with np.load(str(path), allow_pickle=False) as data:
        if "cache_schema_version" not in data or int(data["cache_schema_version"]) != CACHE_SCHEMA_VERSION:
            raise RuntimeError(
                "Stale cache schema in {} (expected version {})".format(
                    path, CACHE_SCHEMA_VERSION))
        views = {}
        for view_name in VIEW_ORDER:
            views[view_name] = {
                key: data[view_name + "_" + key]
                for key in ("u", "v", "mach", "pressure_ratio", "grad_rho", "solid",
                            "source_level")
            }
            views[view_name]["location"] = float(data[view_name + "_location"])
        return {
            "step": int(data["step"]),
            "plotfile": str(data["plotfile"]),
            "sampler": str(data["sampler"]),
            "time_s": float(data["time_s"]),
            "finest_level": int(data["finest_level"]),
            "output_level": int(data["output_level"]),
            "shape": tuple(int(v) for v in data["shape"]),
            "spacing": tuple(float(v) for v in data["spacing"]),
            "views": views,
        }


def field_style(field):
    if field == "mach":
        return {
            "label": r"Mach number, $M$", "cmap": "cividis",
            "norm": Normalize(0.0, 13.0), "ticks": np.arange(0.0, 13.0, 2.0),
            "extend": "max", "formatter": None,
        }
    if field == "pressure_ratio":
        return {
            "label": r"Pressure ratio, $p/p_\infty$", "cmap": "magma",
            "norm": LogNorm(1.0e-2, 1.0e4),
            "ticks": [1.0e-2, 1.0e-1, 1.0, 10.0, 1.0e2, 1.0e3, 1.0e4],
            "extend": "both", "formatter": LogFormatterMathtext(),
        }
    return {
        "label": r"Numerical schlieren, $S$", "cmap": "gray",
        "norm": Normalize(0.0, 1.0), "ticks": [0.0, 0.5, 1.0],
        "extend": "neither", "formatter": None,
    }


def view_data(view, field, gref, schlieren_k):
    if field == "schlieren":
        data = np.exp(-schlieren_k * view["grad_rho"].astype(np.float64) / gref)
        data[~np.isfinite(view["grad_rho"])] = np.nan
        return data
    return view[field]


def view_title(view_name):
    if view_name == "xy":
        return r"$z-z_0=0$ (nozzle 0 plane)"
    if view_name == "xz":
        return r"$y-y_0=-0.126D_b$ (nozzles 1, 2)"
    return r"$x-x_0=-0.787D_b$ (plume section)"


def axis_labels(view_name):
    if view_name == "xy":
        return r"$X$", r"$Y$"
    if view_name == "xz":
        return r"$X$", r"$Z$"
    return r"$Y$", r"$Z$"


def draw_panel(ax, view, view_name, field, gref, schlieren_k,
               panel_label=None, title=True, labels=True):
    style = field_style(field)
    cmap = copy.copy(plt.get_cmap(style["cmap"]))
    cmap.set_bad("white")
    u, v = view["u"], view["v"]
    values = view_data(view, field, gref, schlieren_k)
    mesh = ax.pcolormesh(u, v, np.ma.masked_invalid(values.T), shading="nearest",
                         cmap=cmap, norm=style["norm"], rasterized=True)

    solid = view["solid"].astype(bool)
    if np.any(solid):
        overlay = np.ma.masked_where(~solid.T, np.ones(solid.T.shape, dtype=float))
        ax.pcolormesh(u, v, overlay, shading="nearest",
                      cmap=ListedColormap(["0.88"]), vmin=0.0, vmax=1.0,
                      rasterized=True)
        if np.any(~solid):
            ax.contour(u, v, solid.T.astype(float), levels=[0.5], colors="black",
                       linewidths=0.65)

    xlabel, ylabel = axis_labels(view_name)
    if labels:
        ax.set_xlabel(xlabel, labelpad=1.5)
        ax.set_ylabel(ylabel, labelpad=1.5)
    if title:
        ax.set_title(view_title(view_name), pad=3.0, fontsize=7.5)
    if panel_label:
        ax.text(0.025, 0.965, panel_label, transform=ax.transAxes,
                ha="left", va="top", fontsize=9.0, fontstyle="italic",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72,
                      "pad": 0.8})
    ax.set_xlim(float(u[0]), float(u[-1]))
    ax.set_ylim(float(v[0]), float(v[-1]))
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    return mesh


def add_colorbar(fig, mesh, cax, field):
    style = field_style(field)
    colorbar = fig.colorbar(mesh, cax=cax, ticks=style["ticks"],
                            extend=style["extend"], format=style["formatter"])
    colorbar.set_label(style["label"], labelpad=3.0)
    colorbar.ax.tick_params(labelsize=7.0)
    return colorbar


def save_static(fig, stem, dpi):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(stem) + ".pdf", dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(str(stem) + ".png", dpi=dpi, bbox_inches="tight", pad_inches=0.02)


def render_threeview(cache, field, gref, schlieren_k, stem, dpi,
                     static=True, frame=False):
    panel_ratios = [
        (cache["views"][name]["u"][-1] - cache["views"][name]["u"][0]) /
        (cache["views"][name]["v"][-1] - cache["views"][name]["v"][0])
        for name in VIEW_ORDER
    ]
    fig = plt.figure(figsize=(7.08, 2.55))
    grid = fig.add_gridspec(1, 4, width_ratios=tuple(panel_ratios + [0.055]),
                            left=0.075, right=0.955, bottom=0.19, top=0.82, wspace=0.24)
    axes = [fig.add_subplot(grid[0, col]) for col in range(3)]
    cax = fig.add_subplot(grid[0, 3])
    mesh = None
    for col, view_name in enumerate(VIEW_ORDER):
        mesh = draw_panel(axes[col], cache["views"][view_name], view_name, field,
                          gref, schlieren_k, panel_label="({})".format(chr(97 + col)))
    add_colorbar(fig, mesh, cax, field)
    fig.suptitle(r"$t={:.6f}\ \mathrm{{ms}}$   (plt{:05d})".format(
        1.0e3 * cache["time_s"], cache["step"]), y=0.985, fontsize=8.5)
    if static:
        save_static(fig, stem, dpi)
    if frame:
        path = Path(stem)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(path), dpi=dpi)
    plt.close(fig)


def render_allfields_contact(cache, gref, schlieren_k, stem, dpi):
    panel_ratios = [
        (cache["views"][name]["u"][-1] - cache["views"][name]["u"][0]) /
        (cache["views"][name]["v"][-1] - cache["views"][name]["v"][0])
        for name in VIEW_ORDER
    ]
    fig = plt.figure(figsize=(7.08, 5.75))
    grid = fig.add_gridspec(3, 4, width_ratios=tuple(panel_ratios + [0.055]),
                            left=0.075, right=0.955, bottom=0.075, top=0.92,
                            hspace=0.29, wspace=0.24)
    letter = 0
    for row, field in enumerate(FIELD_ORDER):
        mesh = None
        for col, view_name in enumerate(VIEW_ORDER):
            ax = fig.add_subplot(grid[row, col])
            mesh = draw_panel(ax, cache["views"][view_name], view_name, field,
                              gref, schlieren_k,
                              panel_label="({})".format(chr(97 + letter)),
                              title=(row == 0), labels=True)
            letter += 1
        add_colorbar(fig, mesh, fig.add_subplot(grid[row, 3]), field)
    fig.suptitle(r"$t={:.6f}\ \mathrm{{ms}}$   (plt{:05d})".format(
        1.0e3 * cache["time_s"], cache["step"]), y=0.988, fontsize=9.0)
    save_static(fig, stem, dpi)
    plt.close(fig)


def render_time_contact(caches, view_name, gref, schlieren_k, stem, dpi, title):
    """Three fields by six early times on one fixed physical plane."""
    ncols = len(caches)
    fig = plt.figure(figsize=(7.08, 4.35))
    grid = fig.add_gridspec(3, ncols + 1,
                            width_ratios=tuple([1.0] * ncols + [0.055]),
                            left=0.055, right=0.965, bottom=0.10, top=0.91,
                            hspace=0.24, wspace=0.12)
    for row, field in enumerate(FIELD_ORDER):
        mesh = None
        for col, cache in enumerate(caches):
            ax = fig.add_subplot(grid[row, col])
            mesh = draw_panel(ax, cache["views"][view_name], view_name, field, gref,
                              schlieren_k, panel_label=None, title=False,
                              labels=False)
            if row == 0:
                ax.set_title(r"${:.3f}\ \mathrm{{ms}}$".format(1.0e3 * cache["time_s"]),
                             fontsize=7.3, pad=2.0)
            if row == 2:
                ax.set_xlabel(axis_labels(view_name)[0], labelpad=1.0)
            if col == 0:
                ax.set_ylabel((r"$M$" if field == "mach" else
                               (r"$p/p_\infty$" if field == "pressure_ratio" else r"$S$")),
                              labelpad=1.0)
            else:
                ax.tick_params(axis="y", labelleft=False)
            ax.tick_params(labelsize=6.2)
        add_colorbar(fig, mesh, fig.add_subplot(grid[row, ncols]), field)
    fig.suptitle(title, y=0.985, fontsize=9.0)
    save_static(fig, stem, dpi)
    plt.close(fig)


def render_field_montage(caches, field, gref, schlieren_k, stem, dpi):
    """Seven representative XY states in a compact 2 x 4 publication montage."""
    fig = plt.figure(figsize=(7.08, 4.05))
    grid = fig.add_gridspec(2, 5, width_ratios=(1.0, 1.0, 1.0, 1.0, 0.06),
                            left=0.06, right=0.96, bottom=0.105, top=0.91,
                            hspace=0.25, wspace=0.17)
    mesh = None
    for index, cache in enumerate(caches):
        row, col = divmod(index, 4)
        ax = fig.add_subplot(grid[row, col])
        mesh = draw_panel(ax, cache["views"]["xy"], "xy", field, gref,
                          schlieren_k, panel_label="({})".format(chr(97 + index)),
                          title=False, labels=False)
        ax.set_title(r"$t={:.3f}\ \mathrm{{ms}}$".format(1.0e3 * cache["time_s"]),
                     fontsize=7.5, pad=2.0)
        ax.set_xlabel(r"$X$", labelpad=1.0)
        if col == 0:
            ax.set_ylabel(r"$Y$", labelpad=1.0)
        else:
            ax.tick_params(axis="y", labelleft=False)
    if len(caches) < 8:
        fig.add_subplot(grid[1, 3]).axis("off")
    add_colorbar(fig, mesh, fig.add_subplot(grid[:, 4]), field)
    fig.suptitle("Flow-field development on the nozzle-0 plane", y=0.985,
                 fontsize=9.0)
    save_static(fig, stem, dpi)
    plt.close(fig)


def render_amr_level_contact(caches, stem, dpi):
    """Show the native AMR source level underlying the early XY sequence."""
    ncols = 3
    nrows = int(np.ceil(len(caches) / float(ncols)))
    fig = plt.figure(figsize=(7.08, 3.55))
    grid = fig.add_gridspec(nrows, ncols + 1,
                            width_ratios=(1.0, 1.0, 1.0, 0.045),
                            left=0.07, right=0.955, bottom=0.105, top=0.90,
                            wspace=0.20, hspace=0.30)
    axes = np.empty((nrows, ncols), dtype=object)
    for row in range(nrows):
        for col in range(ncols):
            axes[row, col] = fig.add_subplot(grid[row, col])
    cmap = ListedColormap(["#d9d9d9", "#80b1d3", "#fb8072"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    mesh = None
    for index, cache in enumerate(caches):
        row, col = divmod(index, ncols)
        ax = axes[row, col]
        view = cache["views"]["xy"]
        mesh = ax.pcolormesh(view["u"], view["v"], view["source_level"].T,
                             shading="nearest", cmap=cmap, norm=norm, rasterized=True)
        solid = view["solid"].astype(bool)
        if np.any(solid):
            ax.contour(view["u"], view["v"], solid.T.astype(float), levels=[0.5],
                       colors="black", linewidths=0.65)
        ax.set_title(r"$t={:.3f}\ \mathrm{{ms}}$".format(1.0e3 * cache["time_s"]),
                     fontsize=8.0, pad=2.0)
        ax.set_xlabel(r"$X$", labelpad=1.0)
        if col == 0:
            ax.set_ylabel(r"$Y$", labelpad=1.0)
        else:
            ax.tick_params(axis="y", labelleft=False)
        ax.set_xlim(float(view["u"][0]), float(view["u"][-1]))
        ax.set_ylim(float(view["v"][0]), float(view["v"][-1]))
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(4))
    for index in range(len(caches), nrows * ncols):
        axes.flat[index].axis("off")
    colorbar = fig.colorbar(mesh, cax=fig.add_subplot(grid[:, ncols]), ticks=[0, 1, 2])
    colorbar.set_label("Native AMR source level")
    fig.suptitle("Resolution map for the early startup sequence", y=0.985,
                 fontsize=9.0)
    save_static(fig, stem, dpi)
    plt.close(fig)


def calibration_path(output_dir):
    return Path(output_dir) / "schlieren_calibration.json"


def write_calibration(output_dir, cache, gref, args):
    if args.fixed_gref is None:
        reference_source = "plt{:05d}, pooled cropped fluid cells in XY/XZ/YZ".format(
            cache["step"])
        reference_mode = "computed from this output crop"
    else:
        reference_source = "fixed external G_ref supplied on command line"
        reference_mode = "reused for direct comparison with the near-field paper series"
    payload = {
        "definition": "S = exp(-k * |grad(rho)| / G_ref)",
        "gradient": "full 3-D magnitude at bracketing cell centres, then plane interpolation",
        "reference_source": reference_source,
        "reference_mode": reference_mode,
        "percentile": float(args.schlieren_percentile),
        "G_ref_kg_m-4": float(gref),
        "k": float(args.schlieren_k),
        "fixed_for_all_times": True,
        "axial_crop_m": list(AXIAL_CROP),
        "transverse_crop_m": list(TRANSVERSE_CROP),
    }
    path = calibration_path(output_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def read_calibration(output_dir):
    with calibration_path(output_dir).open() as handle:
        payload = json.load(handle)
    return float(payload["G_ref_kg_m-4"]), float(payload["k"]), payload


def calibrate(args, yt_module, plotfiles):
    step, plotfile = plotfiles[-1]
    path = cache_path(args.output_dir, step)
    result = None
    if path.exists() and not args.overwrite_cache:
        try:
            result = load_cache(path)
            print("Using cached calibration state {}".format(path), flush=True)
        except RuntimeError as error:
            print("Ignoring stale calibration cache: {}".format(error), flush=True)
    if result is None:
        print("Reading calibration state {}".format(plotfile), flush=True)
        result = process_plotfile_yt(yt_module, plotfile)
        save_cache(path, step, plotfile, result)
        result = load_cache(path)
    samples = []
    for view_name in VIEW_ORDER:
        values = result["views"][view_name]["grad_rho"]
        finite = values[np.isfinite(values)]
        if finite.size:
            samples.append(finite.astype(np.float64))
    if not samples:
        raise RuntimeError("No finite fluid density gradients in calibration state")
    if args.fixed_gref is None:
        gref = float(np.percentile(np.concatenate(samples), args.schlieren_percentile))
    else:
        gref = float(args.fixed_gref)
    if not np.isfinite(gref) or gref <= 0.0:
        raise RuntimeError("Invalid schlieren G_ref: {}".format(gref))
    write_calibration(args.output_dir, result, gref, args)
    latest_dir = Path(args.output_dir) / "latest_threeview"
    for field in FIELD_ORDER:
        render_threeview(result, field, gref, args.schlieren_k,
                         latest_dir / (field + "_threeview"), args.paper_dpi)
    render_allfields_contact(result, gref, args.schlieren_k,
                             latest_dir / "allfields_threeview_supplement",
                             args.paper_dpi)
    print("Calibrated G_ref={:.9g} kg m^-4".format(gref), flush=True)
    return gref


def extract(args, yt_module, plotfiles):
    if args.worker_count < 1 or not (0 <= args.worker_index < args.worker_count):
        raise ValueError("Require 0 <= worker-index < worker-count")
    gref, calibrated_k, _ = read_calibration(args.output_dir)
    if abs(calibrated_k - args.schlieren_k) > 1.0e-12:
        raise RuntimeError("Requested schlieren k differs from calibration file")
    assigned = [item for index, item in enumerate(plotfiles)
                if index % args.worker_count == args.worker_index]
    print("Worker {}/{} assigned {} states".format(
        args.worker_index, args.worker_count, len(assigned)), flush=True)
    for item_index, (step, plotfile) in enumerate(assigned, 1):
        path = cache_path(args.output_dir, step)
        result = None
        if path.exists() and not args.overwrite_cache:
            try:
                result = load_cache(path)
                print("[{}/{}] cache plt{:05d}".format(
                    item_index, len(assigned), step), flush=True)
            except RuntimeError as error:
                print("[{}/{}] stale plt{:05d}: {}".format(
                    item_index, len(assigned), step, error), flush=True)
        if result is None:
            print("[{}/{}] read plt{:05d}".format(item_index, len(assigned), step), flush=True)
            result = process_plotfile_yt(yt_module, plotfile)
            save_cache(path, step, plotfile, result)
            result = load_cache(path)

        for field in FIELD_ORDER:
            frame = Path(args.output_dir) / "frames" / field / ("plt{:05d}.png".format(step))
            render_threeview(result, field, gref, args.schlieren_k, frame,
                             args.frame_dpi, static=False, frame=True)
        if step in EARLY_STEPS:
            stem = Path(args.output_dir) / "early_snapshots" / (
                "plt{:05d}_allfields_threeview".format(step))
            render_allfields_contact(result, gref, args.schlieren_k,
                                     stem, args.paper_dpi)
        print("[{}/{}] done plt{:05d}, t={:.6f} ms".format(
            item_index, len(assigned), step, 1.0e3 * result["time_s"]), flush=True)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def assemble(args, plotfiles):
    gref, calibrated_k, calibration = read_calibration(args.output_dir)
    if abs(calibrated_k - args.schlieren_k) > 1.0e-12:
        raise RuntimeError("Requested schlieren k differs from calibration file")
    caches = []
    for step, _ in plotfiles:
        path = cache_path(args.output_dir, step)
        if not path.exists():
            raise RuntimeError("Missing cache: {}".format(path))
        caches.append(load_cache(path))
    steps = [cache["step"] for cache in caches]
    times = [cache["time_s"] for cache in caches]
    if steps != sorted(steps) or np.any(np.diff(times) <= 0.0):
        raise RuntimeError("Cache steps/times are not strictly increasing")

    by_step = {cache["step"]: cache for cache in caches}
    early = [by_step[step] for step in EARLY_STEPS if step in by_step]
    representative = [by_step[step] for step in REPRESENTATIVE_STEPS if step in by_step]
    development_dir = Path(args.output_dir) / "development"
    render_time_contact(
        early, "xy", gref, args.schlieren_k,
        development_dir / "early_allfields_xy", args.paper_dpi,
        "Early flow-field development on the nozzle-0 plane",
    )
    render_time_contact(
        early, "xz", gref, args.schlieren_k,
        development_dir / "early_allfields_xz", args.paper_dpi,
        r"Early flow-field development on the two-nozzle offset plane "
        r"($y-y_0=-0.126D_b$)",
    )
    render_amr_level_contact(early, development_dir / "early_amr_source_level_xy",
                             args.paper_dpi)
    for field in FIELD_ORDER:
        render_field_montage(representative, field, gref, args.schlieren_k,
                             development_dir / (field + "_representative_xy"),
                             args.paper_dpi)

    first_xy = caches[0]["views"]["xy"]
    archive = {
        "steps": np.asarray(steps, dtype=np.int32),
        "times_ms": (1.0e3 * np.asarray(times, dtype=np.float64)),
        "x_over_Db": first_xy["u"],
        "y_over_Db": first_xy["v"],
        "mach_xy": np.stack([cache["views"]["xy"]["mach"] for cache in caches]),
        "pressure_ratio_xy": np.stack(
            [cache["views"]["xy"]["pressure_ratio"] for cache in caches]),
        "pressure_Pa_xy": P_INF * np.stack(
            [cache["views"]["xy"]["pressure_ratio"] for cache in caches]),
        "schlieren_xy": np.stack([
            view_data(cache["views"]["xy"], "schlieren", gref, args.schlieren_k)
            for cache in caches
        ]).astype(np.float32),
        "solid_xy": np.stack([cache["views"]["xy"]["solid"] for cache in caches]),
        "source_level_xy": np.stack(
            [cache["views"]["xy"]["source_level"] for cache in caches]),
        "G_ref_kg_m-4": np.asarray(gref, dtype=np.float64),
        "schlieren_k": np.asarray(args.schlieren_k, dtype=np.float64),
        "p_inf_Pa": np.asarray(P_INF, dtype=np.float64),
        "body_diameter_m": np.asarray(BODY_DIAMETER, dtype=np.float64),
    }
    np.savez_compressed(str(development_dir / "development_xy_data.npz"), **archive)

    output_files = []
    output_root = Path(args.output_dir)
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and "cache" not in path.parts and path.name != "manifest.json":
            output_files.append(str(path.relative_to(output_root)))
    script_path = os.path.abspath(__file__)
    manifest = {
        "case": str(Path(args.plot_root).parent),
        "plot_root": os.path.abspath(args.plot_root),
        "plotfiles": [{"step": cache["step"], "time_ms": 1.0e3 * cache["time_s"],
                       "source": cache["plotfile"]} for cache in caches],
        "fields": {
            "mach": "sqrt(u^2+v^2+w^2) / sqrt(gamma*p/rho), gamma=1.4",
            "pressure_ratio": "p / p_inf, p_inf=574.56 Pa",
            "schlieren": calibration,
        },
        "planes_m": PLANE_LOCATIONS,
        "crop_m": {
            "axial_x": list(AXIAL_CROP),
            "transverse_yz": list(TRANSVERSE_CROP),
        },
        "coordinates": "normalised by D_b=0.127 m about (x0,y0,z0)=(0.300,0.275,0.275) m",
        "mask": "sld < 0.5 and ghs < 0.5, with positive finite rho and p",
        "startup_note": (
            "plt00000 is an imposed startup initial condition with a discontinuity at x=0.25 m; "
            "the early sequence is a startup transient, not a sequence of statistically steady states"
        ),
        "amr_source_level_fractions": {
            "plt{:05d}".format(cache["step"]): {
                view_name: {
                    "L{}".format(level): float(np.mean(
                        cache["views"][view_name]["source_level"] == level))
                    for level in (0, 1, 2)
                } for view_name in VIEW_ORDER
            } for cache in caches
        },
        "fixed_display_limits": {
            "mach": [0.0, 13.0], "pressure_ratio": [1.0e-2, 1.0e4],
            "schlieren": [0.0, 1.0],
        },
        "sampler": caches[0]["sampler"],
        "legacy_reader_not_used": os.path.abspath(args.reader),
        "script": script_path,
        "script_sha256": file_sha256(script_path),
        "outputs": output_files,
    }
    with (output_root / "manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    print("Assembled {} states; manifest has {} output files".format(
        len(caches), len(output_files)), flush=True)


def main():
    args = parse_args()
    configure_crops(args)
    configure_matplotlib()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    plotfiles = discover_plotfiles(args.plot_root)
    print("Discovered {} plotfiles: plt{:05d} ... plt{:05d}".format(
        len(plotfiles), plotfiles[0][0], plotfiles[-1][0]), flush=True)
    print("Crops: x={} m, y/z={} m".format(AXIAL_CROP, TRANSVERSE_CROP), flush=True)
    if args.mode == "all" and args.worker_count != 1:
        raise RuntimeError("mode=all is single-worker only; use calibrate, extract array, assemble")
    yt_module = None
    if args.mode in ("calibrate", "extract", "all"):
        yt_module = load_yt_module()
    if args.mode in ("calibrate", "all"):
        if args.worker_index != 0 or args.worker_count != 1:
            raise RuntimeError("calibrate must run as one non-array worker")
        calibrate(args, yt_module, plotfiles)
    if args.mode in ("extract", "all"):
        extract(args, yt_module, plotfiles)
    if args.mode in ("assemble", "all"):
        assemble(args, plotfiles)


if __name__ == "__main__":
    main()
