#!/usr/bin/env python3
"""Offline, finite-window means for the SRP tri-nozzle paper slices.

This is deliberately a separate pipeline from the instantaneous paper series.
It extracts the seven late snapshots (plt06000 ... plt09000) into its own
cache and then forms a non-uniform-time, trapezoidally weighted mean.  The
primary products are

* Mach number built from Favre velocity, mean density, and mean pressure;
* Reynolds/time-mean pressure;
* numerical schlieren of mean density, using the magnitude of the time-mean
  *signed* density-gradient components.

For comparison, the pipeline also writes the time mean of instantaneous Mach
and the time mean of instantaneous numerical schlieren.  These pairs are not
mathematically interchangeable and are kept under unambiguous names.

The slice extractor follows the ghost-aware yt thin-covering-grid algorithm
used by postprocess_srp_tri_paper_series.py.  It interpolates signed quantities
onto the requested physical plane and never differentiates across IBM cells.
"""

import argparse
import copy
import hashlib
import json
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm, Normalize
from matplotlib.ticker import LogFormatterMathtext, MaxNLocator
import numpy as np


GAMMA = 1.4
P_INF = 574.56
BODY_DIAMETER = 0.127
BODY_CENTRE = (0.300, 0.275, 0.275)
PLANE_LOCATIONS = {"xy": 0.275, "xz": 0.259, "yz": 0.200}
VIEW_ORDER = ("xy", "xz", "yz")
LATE_STEPS = (6000, 6500, 7000, 7500, 8000, 8500, 9000)
AXIAL_CROP = (0.01, 0.60)
TRANSVERSE_CROP = (-0.02, 0.57)
CONTINUOUS_FIELDS = (
    "Density", "pressure", "x_velocity", "y_velocity", "z_velocity",
)
GEOMETRY_FIELDS = ("sld", "ghs")
CACHE_SCHEMA_VERSION = 1
SAMPLER_DESCRIPTION = (
    "yt level-max smoothed_covering_grid thin slabs for primitives; "
    "level-max covering_grid for IBM IDs; signed gradients use centred "
    "differences with one-sided fluid-only fallback"
)
PRIMARY_FIELDS = (
    "favre_mach", "mean_pressure_ratio", "schlieren_of_mean_density",
)
DIAGNOSTIC_FIELDS = (
    "mean_instantaneous_mach", "mean_instantaneous_schlieren",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("extract", "assemble", "all"),
                        default="all")
    parser.add_argument("--plot-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--worker-index", type=int, default=0)
    parser.add_argument("--worker-count", type=int, default=1)
    parser.add_argument("--steps", type=int, nargs="+", default=list(LATE_STEPS))
    parser.add_argument("--axial-crop", type=float, nargs=2, default=AXIAL_CROP,
                        metavar=("X_MIN", "X_MAX"))
    parser.add_argument("--transverse-crop", type=float, nargs=2,
                        default=TRANSVERSE_CROP, metavar=("R_MIN", "R_MAX"))
    parser.add_argument("--schlieren-gref", type=float,
                        default=72.18947715759369)
    parser.add_argument("--schlieren-k", type=float, default=4.0)
    parser.add_argument("--paper-dpi", type=int, default=600)
    parser.add_argument("--overwrite-cache", action="store_true")
    return parser.parse_args()


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
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def load_yt_module():
    try:
        import yt
    except ImportError as error:
        raise RuntimeError(
            "yt is required; load the CX3 yt/4.3 module before extraction"
        ) from error
    yt.funcs.mylog.setLevel(30)
    return yt


def discover_selected_plotfiles(root, requested_steps):
    found = {}
    for path in Path(root).glob("plt[0-9]*"):
        match = re.fullmatch(r"plt(\d+)", path.name)
        if path.is_dir() and match and (path / "Header").is_file():
            found[int(match.group(1))] = str(path)
    steps = tuple(sorted(set(int(step) for step in requested_steps)))
    missing = [step for step in steps if step not in found]
    if missing:
        raise RuntimeError("Missing complete requested plotfiles: {}".format(missing))
    if len(steps) < 2:
        raise RuntimeError("At least two snapshots are required for a time mean")
    return [(step, found[step]) for step in steps]


def validate_args(args):
    axial = tuple(float(value) for value in args.axial_crop)
    transverse = tuple(float(value) for value in args.transverse_crop)
    if not (axial[0] < axial[1] and transverse[0] < transverse[1]):
        raise ValueError("Crop lower bounds must be smaller than upper bounds")
    if not (np.isfinite(args.schlieren_gref) and args.schlieren_gref > 0.0):
        raise ValueError("schlieren G_ref must be finite and positive")
    if not (np.isfinite(args.schlieren_k) and args.schlieren_k > 0.0):
        raise ValueError("schlieren k must be finite and positive")
    if args.worker_count < 1 or not (0 <= args.worker_index < args.worker_count):
        raise ValueError("Require 0 <= worker-index < worker-count")
    if args.mode == "all" and (args.worker_index != 0 or args.worker_count != 1):
        raise ValueError("mode=all is single-worker only")
    return axial, transverse


def _yt_array(values):
    if hasattr(values, "to_value"):
        return np.asarray(values.to_value(), dtype=np.float64)
    if hasattr(values, "d"):
        return np.asarray(values.d, dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


def bracket_stencil(coords, value):
    hi = int(np.searchsorted(coords, value, side="left"))
    hi = max(1, min(hi, len(coords) - 1))
    lo = hi - 1
    weight = float((value - coords[lo]) / (coords[hi] - coords[lo]))
    weight = max(0.0, min(1.0, weight))
    indices = sorted(set(max(0, min(len(coords) - 1, index))
                         for index in (lo - 1, lo, hi, hi + 1)))
    return lo, hi, weight, indices


def _inclusive_cell_range(coords, bounds):
    selected = np.flatnonzero((coords >= bounds[0]) & (coords <= bounds[1]))
    if selected.size < 2:
        raise RuntimeError("Crop {} contains too few cells".format(bounds))
    return int(selected[0]), int(selected[-1]) + 1


def _take_plane(array, axis, index):
    return np.take(array, int(index), axis=axis)


def plane_valid(rho, sld, ghs):
    return (np.isfinite(rho) & np.isfinite(sld) & np.isfinite(ghs) &
            (rho > 0.0) & (sld < 0.5) & (ghs < 0.5))


def masked_derivative(field, valid, spacing, axis):
    """Centred derivative with one-sided fallback, never crossing IBM cells."""
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
    previous = index - 1 if index - 1 in planes else None
    following = index + 1 if index + 1 in planes else None
    result = np.full(centre.shape, np.nan, dtype=np.float64)
    central = np.zeros(centre.shape, dtype=bool)
    if previous is not None and following is not None:
        central = valid & valid_planes[previous] & valid_planes[following]
        result[central] = ((planes[following][central] - planes[previous][central]) /
                           (2.0 * spacing))
    if following is not None:
        forward = valid & valid_planes[following] & ~central
        result[forward] = (planes[following][forward] - centre[forward]) / spacing
    if previous is not None:
        backward = valid & valid_planes[previous] & ~central & ~np.isfinite(result)
        result[backward] = (centre[backward] - planes[previous][backward]) / spacing
    return result


def gradient_components_at(rho_planes, sld_planes, ghs_planes, index,
                           normal_spacing, inplane_spacings):
    valid_planes = {
        key: plane_valid(rho_planes[key], sld_planes[key], ghs_planes[key])
        for key in rho_planes
    }
    valid = valid_planes[index]
    gu = masked_derivative(rho_planes[index], valid, inplane_spacings[0], axis=0)
    gv = masked_derivative(rho_planes[index], valid, inplane_spacings[1], axis=1)
    gn = normal_derivative_at(rho_planes, valid_planes, index, normal_spacing)
    for component in (gu, gv, gn):
        component[~valid] = np.nan
    return gu, gv, gn


def source_level_at(ds, u_coords, v_coords, u_axis, v_axis,
                    normal_axis, normal_position):
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


def extract_view(ds, direction, fine_coords, spacings, axial_crop, transverse_crop):
    """Extract one physical plane from a ghost-aware, finest-level thin slab."""
    normal_axis = {"yz": 0, "xz": 1, "xy": 2}[direction]
    inplane_axes = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[normal_axis]
    location = PLANE_LOCATIONS[direction]
    lo, hi, plane_weight, stencil = bracket_stencil(fine_coords[normal_axis], location)

    starts = np.zeros(3, dtype=np.int64)
    stops = np.zeros(3, dtype=np.int64)
    for axis in range(3):
        if axis == normal_axis:
            starts[axis] = int(stencil[0])
            stops[axis] = int(stencil[-1]) + 1
        else:
            bounds = axial_crop if axis == 0 else transverse_crop
            starts[axis], stops[axis] = _inclusive_cell_range(fine_coords[axis], bounds)

    # This is the same non-periodic-domain safeguard used by the paper series.
    fine_dims = np.asarray([len(values) for values in fine_coords], dtype=np.int64)
    interpolation_halo = 2 * int(ds.refine_by) ** int(ds.max_level)
    starts = np.maximum(starts, interpolation_halo)
    stops = np.minimum(stops, fine_dims - interpolation_halo)
    dims = stops - starts
    if np.any(dims < 2):
        raise RuntimeError("Ghost-aware crop left too few cells: {}".format(dims))
    left_edge = _yt_array(ds.domain_left_edge) + starts * np.asarray(spacings)

    smooth = ds.smoothed_covering_grid(
        int(ds.max_level), left_edge=left_edge, dims=dims,
        fields=[("boxlib", name) for name in CONTINUOUS_FIELDS],
        num_ghost_zones=0, use_pbar=False,
    )
    categorical = ds.covering_grid(
        int(ds.max_level), left_edge=left_edge, dims=dims,
        fields=[("boxlib", name) for name in GEOMETRY_FIELDS],
        num_ghost_zones=0, use_pbar=False,
    )
    arrays = {name: _yt_array(smooth[("boxlib", name)])
              for name in CONTINUOUS_FIELDS}
    arrays.update({name: _yt_array(categorical[("boxlib", name)])
                   for name in GEOMETRY_FIELDS})

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
        rho_planes, sld_planes, ghs_planes, local_lo, spacings[normal_axis],
        (spacings[inplane_axes[0]], spacings[inplane_axes[1]]),
    )
    components_hi = gradient_components_at(
        rho_planes, sld_planes, ghs_planes, local_hi, spacings[normal_axis],
        (spacings[inplane_axes[0]], spacings[inplane_axes[1]]),
    )
    components_uvn = tuple(
        (1.0 - plane_weight) * lower + plane_weight * upper
        for lower, upper in zip(components_lo, components_hi)
    )
    components_xyz = [None, None, None]
    components_xyz[inplane_axes[0]] = components_uvn[0]
    components_xyz[inplane_axes[1]] = components_uvn[1]
    components_xyz[normal_axis] = components_uvn[2]

    plane_values = {}
    for name, array in arrays.items():
        lower = _take_plane(array, normal_axis, local_lo)
        upper = _take_plane(array, normal_axis, local_hi)
        plane_values[name] = (lower, upper,
                              (1.0 - plane_weight) * lower + plane_weight * upper)

    rho_lo, rho_hi, rho = plane_values["Density"]
    p_lo, p_hi, pressure = plane_values["pressure"]
    sld_lo, sld_hi, _ = plane_values["sld"]
    ghs_lo, ghs_hi, _ = plane_values["ghs"]
    valid = (plane_valid(rho_lo, sld_lo, ghs_lo) &
             plane_valid(rho_hi, sld_hi, ghs_hi) &
             np.isfinite(p_lo) & np.isfinite(p_hi) & (p_lo > 0.0) & (p_hi > 0.0))
    velocities = []
    for name in ("x_velocity", "y_velocity", "z_velocity"):
        lower, upper, value = plane_values[name]
        valid &= np.isfinite(lower) & np.isfinite(upper)
        velocities.append(value)
    ux, uy, uz = velocities
    speed = np.sqrt(ux * ux + uy * uy + uz * uz)
    sound = np.sqrt(GAMMA * pressure / rho)
    instant_mach = speed / sound
    grad_magnitude = np.sqrt(sum(component * component
                                 for component in components_xyz))
    solid = ((sld_lo >= 0.5) | (sld_hi >= 0.5))

    flow_arrays = [rho, pressure, ux, uy, uz, instant_mach]
    for array in flow_arrays:
        array[~valid] = np.nan
    for component in components_xyz:
        component[~valid] = np.nan
    grad_magnitude[~valid] = np.nan

    u_axis, v_axis = inplane_axes
    u_raw = fine_coords[u_axis][starts[u_axis]:stops[u_axis]]
    v_raw = fine_coords[v_axis][starts[v_axis]:stops[v_axis]]
    source_lo = source_level_at(ds, u_raw, v_raw, u_axis, v_axis,
                                normal_axis, fine_coords[normal_axis][lo])
    source_hi = source_level_at(ds, u_raw, v_raw, u_axis, v_axis,
                                normal_axis, fine_coords[normal_axis][hi])
    source_level = np.minimum(source_lo, source_hi)
    offsets = np.asarray(BODY_CENTRE)
    u = (u_raw - offsets[u_axis]) / BODY_DIAMETER
    v = (v_raw - offsets[v_axis]) / BODY_DIAMETER

    return {
        "u": u.astype(np.float64),
        "v": v.astype(np.float64),
        "rho": rho.astype(np.float32),
        "pressure": pressure.astype(np.float32),
        "rho_u": (rho * ux).astype(np.float32),
        "rho_v": (rho * uy).astype(np.float32),
        "rho_w": (rho * uz).astype(np.float32),
        "instantaneous_mach": instant_mach.astype(np.float32),
        "grad_rho_x": components_xyz[0].astype(np.float32),
        "grad_rho_y": components_xyz[1].astype(np.float32),
        "grad_rho_z": components_xyz[2].astype(np.float32),
        "instantaneous_grad_rho_magnitude": grad_magnitude.astype(np.float32),
        "solid": solid.astype(np.uint8),
        "source_level": source_level.astype(np.uint8),
        "location": float(location),
        "bracket_lo": int(lo),
        "bracket_hi": int(hi),
        "plane_weight": float(plane_weight),
        "effective_u_m": [float(u_raw[0]), float(u_raw[-1])],
        "effective_v_m": [float(v_raw[0]), float(v_raw[-1])],
    }


def process_plotfile(yt_module, plotfile, axial_crop, transverse_crop):
    dataset = yt_module.load(plotfile)
    base_dims = np.asarray(dataset.domain_dimensions, dtype=np.int64)
    fine_dims = base_dims * int(dataset.refine_by) ** int(dataset.max_level)
    left = _yt_array(dataset.domain_left_edge)
    right = _yt_array(dataset.domain_right_edge)
    spacings = tuple(float(value) for value in (right - left) / fine_dims)
    fine_coords = tuple(
        left[axis] + (np.arange(fine_dims[axis]) + 0.5) * spacings[axis]
        for axis in range(3)
    )
    views = {
        name: extract_view(dataset, name, fine_coords, spacings,
                           axial_crop, transverse_crop)
        for name in VIEW_ORDER
    }
    result = {
        "time_s": float(_yt_array(dataset.current_time)),
        "finest_level": int(dataset.max_level),
        "shape": tuple(int(value) for value in fine_dims),
        "spacing": spacings,
        "sampler": SAMPLER_DESCRIPTION,
        "views": views,
    }
    del dataset
    return result


CACHE_ARRAY_KEYS = (
    "u", "v", "rho", "pressure", "rho_u", "rho_v", "rho_w",
    "instantaneous_mach", "grad_rho_x", "grad_rho_y", "grad_rho_z",
    "instantaneous_grad_rho_magnitude", "solid", "source_level",
)


def cache_path(output_dir, step):
    return Path(output_dir) / "cache" / ("plt{:05d}.npz".format(step))


def save_cache(path, step, plotfile, result, axial_crop, transverse_crop):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_schema_version": np.asarray(CACHE_SCHEMA_VERSION, dtype=np.int32),
        "step": np.asarray(step, dtype=np.int64),
        "plotfile": np.asarray(plotfile),
        "time_s": np.asarray(result["time_s"], dtype=np.float64),
        "finest_level": np.asarray(result["finest_level"], dtype=np.int32),
        "shape": np.asarray(result["shape"], dtype=np.int32),
        "spacing": np.asarray(result["spacing"], dtype=np.float64),
        "sampler": np.asarray(result["sampler"]),
        "axial_crop_m": np.asarray(axial_crop, dtype=np.float64),
        "transverse_crop_m": np.asarray(transverse_crop, dtype=np.float64),
    }
    for view_name, view in result["views"].items():
        for key in CACHE_ARRAY_KEYS:
            payload[view_name + "_" + key] = view[key]
        payload[view_name + "_location"] = np.asarray(view["location"], dtype=np.float64)
        payload[view_name + "_bracket"] = np.asarray(
            [view["bracket_lo"], view["bracket_hi"], view["plane_weight"]],
            dtype=np.float64,
        )
        payload[view_name + "_effective_u_m"] = np.asarray(
            view["effective_u_m"], dtype=np.float64)
        payload[view_name + "_effective_v_m"] = np.asarray(
            view["effective_v_m"], dtype=np.float64)
    temporary = str(path) + ".tmp.npz"
    np.savez_compressed(temporary, **payload)
    os.replace(temporary, str(path))


def load_cache(path, axial_crop, transverse_crop):
    with np.load(str(path), allow_pickle=False) as data:
        if ("cache_schema_version" not in data or
                int(data["cache_schema_version"]) != CACHE_SCHEMA_VERSION):
            raise RuntimeError("Stale mean-cache schema: {}".format(path))
        if not np.allclose(data["axial_crop_m"], axial_crop, rtol=0.0, atol=1.0e-14):
            raise RuntimeError("Mean-cache axial crop mismatch: {}".format(path))
        if not np.allclose(data["transverse_crop_m"], transverse_crop,
                           rtol=0.0, atol=1.0e-14):
            raise RuntimeError("Mean-cache transverse crop mismatch: {}".format(path))
        views = {}
        for view_name in VIEW_ORDER:
            view = {key: data[view_name + "_" + key] for key in CACHE_ARRAY_KEYS}
            view["location"] = float(data[view_name + "_location"])
            view["effective_u_m"] = data[view_name + "_effective_u_m"]
            view["effective_v_m"] = data[view_name + "_effective_v_m"]
            views[view_name] = view
        return {
            "step": int(data["step"]),
            "plotfile": str(data["plotfile"]),
            "time_s": float(data["time_s"]),
            "finest_level": int(data["finest_level"]),
            "shape": tuple(int(value) for value in data["shape"]),
            "spacing": tuple(float(value) for value in data["spacing"]),
            "sampler": str(data["sampler"]),
            "views": views,
        }


def extract(args, yt_module, plotfiles, axial_crop, transverse_crop):
    assigned = [item for index, item in enumerate(plotfiles)
                if index % args.worker_count == args.worker_index]
    print("Worker {}/{} assigned {} late states".format(
        args.worker_index, args.worker_count, len(assigned)), flush=True)
    for item_index, (step, plotfile) in enumerate(assigned, 1):
        path = cache_path(args.output_dir, step)
        result = None
        if path.exists() and not args.overwrite_cache:
            try:
                result = load_cache(path, axial_crop, transverse_crop)
                print("[{}/{}] cache plt{:05d}".format(
                    item_index, len(assigned), step), flush=True)
            except RuntimeError as error:
                print("[{}/{}] ignore {}".format(item_index, len(assigned), error),
                      flush=True)
        if result is None:
            print("[{}/{}] read plt{:05d}".format(
                item_index, len(assigned), step), flush=True)
            result = process_plotfile(yt_module, plotfile, axial_crop, transverse_crop)
            save_cache(path, step, plotfile, result, axial_crop, transverse_crop)
        print("[{}/{}] done plt{:05d}, t={:.6f} ms".format(
            item_index, len(assigned), step, 1.0e3 * result["time_s"]), flush=True)


def trapezoidal_time_weights(times):
    times = np.asarray(times, dtype=np.float64)
    if times.ndim != 1 or times.size < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("Snapshot times must be one-dimensional and increasing")
    weights = np.empty(times.size, dtype=np.float64)
    weights[0] = 0.5 * (times[1] - times[0])
    weights[-1] = 0.5 * (times[-1] - times[-2])
    weights[1:-1] = 0.5 * (times[2:] - times[:-2])
    weights /= times[-1] - times[0]
    if not np.isclose(np.sum(weights), 1.0, rtol=0.0, atol=5.0e-15):
        raise RuntimeError("Trapezoidal weights do not sum to unity")
    return weights


def complete_record_weighted_mean(arrays, weights):
    arrays = [np.asarray(array, dtype=np.float64) for array in arrays]
    common = np.logical_and.reduce([np.isfinite(array) for array in arrays])
    result = np.zeros(arrays[0].shape, dtype=np.float64)
    for weight, array in zip(weights, arrays):
        result[common] += weight * array[common]
    result[~common] = np.nan
    return result, common


def assemble_view(caches, view_name, weights, gref, schlieren_k):
    views = [cache["views"][view_name] for cache in caches]
    reference = views[0]
    for view in views[1:]:
        if (view["u"].shape != reference["u"].shape or
                view["v"].shape != reference["v"].shape or
                not np.allclose(view["u"], reference["u"], rtol=0.0, atol=1.0e-12) or
                not np.allclose(view["v"], reference["v"], rtol=0.0, atol=1.0e-12)):
            raise RuntimeError("Grid mismatch across {} late caches".format(view_name))

    means = {}
    common_masks = {}
    for key in ("rho", "pressure", "rho_u", "rho_v", "rho_w",
                "instantaneous_mach", "grad_rho_x", "grad_rho_y", "grad_rho_z",
                "instantaneous_grad_rho_magnitude"):
        means[key], common_masks[key] = complete_record_weighted_mean(
            [view[key] for view in views], weights)

    flow_valid = np.logical_and.reduce([
        common_masks[key] for key in ("rho", "pressure", "rho_u", "rho_v", "rho_w",
                                      "instantaneous_mach")
    ])
    gradient_valid = np.logical_and.reduce([
        common_masks[key] for key in ("grad_rho_x", "grad_rho_y", "grad_rho_z",
                                      "instantaneous_grad_rho_magnitude")
    ])
    rho_mean = means["rho"]
    pressure_mean = means["pressure"]
    tiny = np.finfo(np.float64).tiny
    favre_u = means["rho_u"] / np.maximum(rho_mean, tiny)
    favre_v = means["rho_v"] / np.maximum(rho_mean, tiny)
    favre_w = means["rho_w"] / np.maximum(rho_mean, tiny)
    favre_speed = np.sqrt(favre_u ** 2 + favre_v ** 2 + favre_w ** 2)
    mean_sound_speed = np.sqrt(GAMMA * pressure_mean / np.maximum(rho_mean, tiny))
    favre_mach = favre_speed / np.maximum(mean_sound_speed, tiny)
    favre_mach[~flow_valid] = np.nan
    mean_instantaneous_mach = means["instantaneous_mach"]
    mean_instantaneous_mach[~flow_valid] = np.nan

    grad_mean_magnitude = np.sqrt(
        means["grad_rho_x"] ** 2 + means["grad_rho_y"] ** 2 +
        means["grad_rho_z"] ** 2
    )
    schlieren_of_mean_density = np.exp(-schlieren_k * grad_mean_magnitude / gref)
    schlieren_of_mean_density[~gradient_valid] = np.nan
    instantaneous_schlieren = []
    for view in views:
        values = np.exp(-schlieren_k *
                        view["instantaneous_grad_rho_magnitude"].astype(np.float64) /
                        gref)
        values[~np.isfinite(view["instantaneous_grad_rho_magnitude"])] = np.nan
        instantaneous_schlieren.append(values)
    mean_instantaneous_schlieren, schlieren_record = complete_record_weighted_mean(
        instantaneous_schlieren, weights)
    mean_instantaneous_schlieren[~(gradient_valid & schlieren_record)] = np.nan

    solid = np.logical_or.reduce([view["solid"].astype(bool) for view in views])
    source_level_min = np.minimum.reduce([view["source_level"] for view in views])
    source_level_max = np.maximum.reduce([view["source_level"] for view in views])
    output = {
        "u": reference["u"].astype(np.float64),
        "v": reference["v"].astype(np.float64),
        "rho_mean": rho_mean.astype(np.float32),
        "pressure_mean_Pa": pressure_mean.astype(np.float32),
        "mean_pressure_ratio": (pressure_mean / P_INF).astype(np.float32),
        "rho_u_mean": means["rho_u"].astype(np.float32),
        "rho_v_mean": means["rho_v"].astype(np.float32),
        "rho_w_mean": means["rho_w"].astype(np.float32),
        "favre_u": favre_u.astype(np.float32),
        "favre_v": favre_v.astype(np.float32),
        "favre_w": favre_w.astype(np.float32),
        "favre_mach": favre_mach.astype(np.float32),
        "mean_instantaneous_mach": mean_instantaneous_mach.astype(np.float32),
        "grad_rho_mean_x": means["grad_rho_x"].astype(np.float32),
        "grad_rho_mean_y": means["grad_rho_y"].astype(np.float32),
        "grad_rho_mean_z": means["grad_rho_z"].astype(np.float32),
        "grad_rho_mean_magnitude": grad_mean_magnitude.astype(np.float32),
        "mean_instantaneous_grad_rho_magnitude":
            means["instantaneous_grad_rho_magnitude"].astype(np.float32),
        "schlieren_of_mean_density": schlieren_of_mean_density.astype(np.float32),
        "mean_instantaneous_schlieren": mean_instantaneous_schlieren.astype(np.float32),
        "flow_complete_record": flow_valid.astype(np.uint8),
        "gradient_complete_record": gradient_valid.astype(np.uint8),
        "solid": solid.astype(np.uint8),
        "source_level_min": source_level_min.astype(np.uint8),
        "source_level_max": source_level_max.astype(np.uint8),
        "location": float(reference["location"]),
        "effective_u_m": reference["effective_u_m"].astype(np.float64),
        "effective_v_m": reference["effective_v_m"].astype(np.float64),
    }
    for key, array in output.items():
        if isinstance(array, np.ndarray) and array.ndim == 2 and key not in (
                "solid", "source_level_min", "source_level_max",
                "flow_complete_record", "gradient_complete_record"):
            array[solid] = np.nan
    return output


def field_style(field):
    if field in ("favre_mach", "mean_instantaneous_mach"):
        label = (r"Favre-mean-state Mach, $M_{\widetilde{\mathbf{u}}}$"
                 if field == "favre_mach" else
                 r"Mean instantaneous Mach, $\overline{M}$")
        return {
            "label": label, "cmap": "cividis", "norm": Normalize(0.0, 13.0),
            "ticks": np.arange(0.0, 13.0, 2.0), "extend": "max",
            "formatter": None,
        }
    if field == "mean_pressure_ratio":
        return {
            "label": r"Mean pressure ratio, $\overline{p}/p_\infty$",
            "cmap": "magma", "norm": LogNorm(1.0e-2, 1.0e4),
            "ticks": [1.0e-2, 1.0e-1, 1.0, 10.0, 1.0e2, 1.0e3, 1.0e4],
            "extend": "both", "formatter": LogFormatterMathtext(),
        }
    if field == "schlieren_of_mean_density":
        label = r"Schlieren of mean density, $S[\overline{\rho}]$"
    elif field == "mean_instantaneous_schlieren":
        label = r"Mean instantaneous schlieren, $\overline{S[\rho]}$"
    else:
        raise KeyError(field)
    return {
        "label": label, "cmap": "gray", "norm": Normalize(0.0, 1.0),
        "ticks": [0.0, 0.5, 1.0], "extend": "neither", "formatter": None,
    }


def axis_labels(view_name):
    if view_name == "xy":
        return r"$X$", r"$Y$"
    if view_name == "xz":
        return r"$X$", r"$Z$"
    return r"$Y$", r"$Z$"


def view_title(view_name):
    if view_name == "xy":
        return r"$z-z_0=0$ (nozzle 0 plane)"
    if view_name == "xz":
        return r"$y-y_0=-0.126D_b$ (nozzles 1, 2)"
    return r"$x-x_0=-0.787D_b$ (plume section)"


def draw_panel(ax, view, view_name, field, panel_label=None, title=True):
    style = field_style(field)
    cmap = copy.copy(plt.get_cmap(style["cmap"]))
    cmap.set_bad("white")
    values = view[field]
    mesh = ax.pcolormesh(view["u"], view["v"],
                         np.ma.masked_invalid(values.T), shading="nearest",
                         cmap=cmap, norm=style["norm"], rasterized=True)
    solid = view["solid"].astype(bool)
    if np.any(solid):
        overlay = np.ma.masked_where(~solid.T, np.ones(solid.T.shape, dtype=float))
        ax.pcolormesh(view["u"], view["v"], overlay, shading="nearest",
                      cmap=ListedColormap(["0.88"]), vmin=0.0, vmax=1.0,
                      rasterized=True)
        if np.any(~solid):
            ax.contour(view["u"], view["v"], solid.T.astype(float), levels=[0.5],
                       colors="black", linewidths=0.65)
    xlabel, ylabel = axis_labels(view_name)
    ax.set_xlabel(xlabel, labelpad=1.5)
    ax.set_ylabel(ylabel, labelpad=1.5)
    if title:
        ax.set_title(view_title(view_name), pad=3.0, fontsize=7.5)
    if panel_label:
        ax.text(0.025, 0.965, panel_label, transform=ax.transAxes,
                ha="left", va="top", fontsize=9.0, fontstyle="italic",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72,
                      "pad": 0.8})
    ax.set_xlim(float(view["u"][0]), float(view["u"][-1]))
    ax.set_ylim(float(view["v"][0]), float(view["v"][-1]))
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


def save_static(fig, stem, dpi):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(stem) + ".pdf", dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(str(stem) + ".png", dpi=dpi, bbox_inches="tight", pad_inches=0.02)


def render_threeview(mean_views, field, stem, dpi, interval_ms):
    ratios = [
        (mean_views[name]["u"][-1] - mean_views[name]["u"][0]) /
        (mean_views[name]["v"][-1] - mean_views[name]["v"][0])
        for name in VIEW_ORDER
    ]
    fig = plt.figure(figsize=(7.08, 2.55))
    grid = fig.add_gridspec(1, 4, width_ratios=tuple(ratios + [0.055]),
                            left=0.075, right=0.955, bottom=0.19, top=0.82,
                            wspace=0.24)
    mesh = None
    for col, view_name in enumerate(VIEW_ORDER):
        mesh = draw_panel(fig.add_subplot(grid[0, col]), mean_views[view_name],
                          view_name, field,
                          panel_label="({})".format(chr(97 + col)))
    add_colorbar(fig, mesh, fig.add_subplot(grid[0, 3]), field)
    fig.suptitle(
        r"Offline trapezoidal time mean, $t=[{:.6f},{:.6f}]\ \mathrm{{ms}}$".format(
            interval_ms[0], interval_ms[1]), y=0.985, fontsize=8.5)
    save_static(fig, stem, dpi)
    plt.close(fig)


def render_contact(mean_views, stem, dpi, interval_ms):
    fields = PRIMARY_FIELDS + DIAGNOSTIC_FIELDS
    ratios = [
        (mean_views[name]["u"][-1] - mean_views[name]["u"][0]) /
        (mean_views[name]["v"][-1] - mean_views[name]["v"][0])
        for name in VIEW_ORDER
    ]
    fig = plt.figure(figsize=(7.08, 9.15))
    grid = fig.add_gridspec(len(fields), 4, width_ratios=tuple(ratios + [0.055]),
                            left=0.075, right=0.955, bottom=0.055, top=0.945,
                            hspace=0.30, wspace=0.24)
    letter = 0
    for row, field in enumerate(fields):
        mesh = None
        for col, view_name in enumerate(VIEW_ORDER):
            mesh = draw_panel(fig.add_subplot(grid[row, col]), mean_views[view_name],
                              view_name, field,
                              panel_label="({})".format(chr(97 + letter)),
                              title=(row == 0))
            letter += 1
        add_colorbar(fig, mesh, fig.add_subplot(grid[row, 3]), field)
    fig.suptitle(
        r"Late-window offline means, $t=[{:.6f},{:.6f}]\ \mathrm{{ms}}$".format(
            interval_ms[0], interval_ms[1]), y=0.993, fontsize=9.0)
    save_static(fig, stem, dpi)
    plt.close(fig)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def save_mean_archive(path, mean_views, steps, times, weights, args):
    payload = {
        "steps": np.asarray(steps, dtype=np.int32),
        "times_s": np.asarray(times, dtype=np.float64),
        "times_ms": 1.0e3 * np.asarray(times, dtype=np.float64),
        "trapezoidal_time_weights": np.asarray(weights, dtype=np.float64),
        "interval_s": np.asarray([times[0], times[-1]], dtype=np.float64),
        "gamma": np.asarray(GAMMA, dtype=np.float64),
        "p_inf_Pa": np.asarray(P_INF, dtype=np.float64),
        "body_diameter_m": np.asarray(BODY_DIAMETER, dtype=np.float64),
        "schlieren_G_ref_kg_m-4": np.asarray(args.schlieren_gref, dtype=np.float64),
        "schlieren_k": np.asarray(args.schlieren_k, dtype=np.float64),
        "requested_axial_crop_m": np.asarray(args.axial_crop, dtype=np.float64),
        "requested_transverse_crop_m": np.asarray(args.transverse_crop,
                                                   dtype=np.float64),
    }
    for view_name, view in mean_views.items():
        for key, value in view.items():
            if isinstance(value, np.ndarray):
                payload[view_name + "_" + key] = value
            elif key == "location":
                payload[view_name + "_location_m"] = np.asarray(value, dtype=np.float64)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(str(path), **payload)


def definitions_payload(caches, mean_views, weights, args, axial_crop,
                        transverse_crop, plotfiles):
    times = np.asarray([cache["time_s"] for cache in caches], dtype=np.float64)
    script_path = os.path.abspath(__file__)
    return {
        "product": "offline finite-window time mean on three physical slices",
        "native_solver_statistics": False,
        "stationarity_claimed": False,
        "interpretation_note": (
            "This seven-snapshot late-window mean is a diagnostic finite-window mean. "
            "It is not a native solver statistical accumulator, and statistical "
            "stationarity has not been established."
        ),
        "steps": [int(cache["step"]) for cache in caches],
        "plotfiles": [source for _, source in plotfiles],
        "times_s": times.tolist(),
        "times_ms": (1.0e3 * times).tolist(),
        "averaging_interval_s": [float(times[0]), float(times[-1])],
        "averaging_interval_ms": [float(1.0e3 * times[0]),
                                   float(1.0e3 * times[-1])],
        "quadrature": {
            "method": "composite trapezoidal rule for non-uniform snapshot times",
            "normalised_weights": [float(value) for value in weights],
            "equation": (
                "<f>_T = sum_i(w_i f_i); w_0=dt_0/(2T), "
                "w_N=dt_{N-1}/(2T), w_i=(t_{i+1}-t_{i-1})/(2T)"
            ),
            "complete_record_policy": (
                "a pixel is published only when the required quantity is finite at "
                "all seven snapshots; no temporal renormalisation around missing values"
            ),
        },
        "definitions": {
            "mean_density": "rho_bar = <rho>_T",
            "mean_pressure": "p_bar = <p>_T (Reynolds/time mean)",
            "favre_velocity": "u_tilde_j = <rho u_j>_T / <rho>_T",
            "favre_mach": (
                "M_favre = |u_tilde| / sqrt(gamma p_bar/rho_bar), gamma=1.4; "
                "this is Mach formed from the Favre velocity and mean thermodynamic state"
            ),
            "mean_instantaneous_mach": (
                "M_inst_bar = < |u|/sqrt(gamma p/rho) >_T; generally not equal to M_favre"
            ),
            "signed_mean_density_gradient": (
                "grad(rho)_bar = (<d rho/dx>_T,<d rho/dy>_T,<d rho/dz>_T); "
                "signed components are averaged before taking magnitude"
            ),
            "schlieren_of_mean_density": (
                "S[rho_bar] = exp(-k |grad(rho)_bar|/G_ref)"
            ),
            "mean_instantaneous_schlieren": (
                "S_inst_bar = <exp(-k |grad rho|/G_ref)>_T; generally not equal "
                "to S[rho_bar]"
            ),
        },
        "schlieren": {
            "k": float(args.schlieren_k),
            "G_ref_kg_m-4": float(args.schlieren_gref),
            "reference_note": (
                "fixed external reference reused from the instantaneous paper series "
                "to preserve visual comparability"
            ),
        },
        "planes_m": PLANE_LOCATIONS,
        "requested_crop_m": {
            "axial_x": list(axial_crop),
            "transverse_yz": list(transverse_crop),
        },
        "effective_crop_m": {
            name: {
                "u": [float(value) for value in mean_views[name]["effective_u_m"]],
                "v": [float(value) for value in mean_views[name]["effective_v_m"]],
            } for name in VIEW_ORDER
        },
        "coordinates": (
            "plot coordinates are normalised by D_b=0.127 m about "
            "(x0,y0,z0)=(0.300,0.275,0.275) m"
        ),
        "mask": "sld < 0.5 and ghs < 0.5 on both bracketing planes, with rho,p>0",
        "sampler": caches[0]["sampler"],
        "source_level_note": (
            "source_level_min/max record the native AMR levels contributing across "
            "the seven snapshots; all sampled arrays are represented on the Lmax grid"
        ),
        "script": script_path,
        "script_sha256": file_sha256(script_path),
        "archive": "late_mean_fields.npz",
    }


def assemble(args, plotfiles, axial_crop, transverse_crop):
    caches = []
    for step, _ in plotfiles:
        path = cache_path(args.output_dir, step)
        if not path.exists():
            raise RuntimeError("Missing extraction cache: {}".format(path))
        caches.append(load_cache(path, axial_crop, transverse_crop))
    steps = [cache["step"] for cache in caches]
    times = np.asarray([cache["time_s"] for cache in caches], dtype=np.float64)
    weights = trapezoidal_time_weights(times)
    mean_views = {
        view_name: assemble_view(caches, view_name, weights,
                                 args.schlieren_gref, args.schlieren_k)
        for view_name in VIEW_ORDER
    }
    interval_ms = (float(1.0e3 * times[0]), float(1.0e3 * times[-1]))
    figures = Path(args.output_dir) / "figures"
    for field in PRIMARY_FIELDS + DIAGNOSTIC_FIELDS:
        render_threeview(mean_views, field, figures / (field + "_threeview"),
                         args.paper_dpi, interval_ms)
    render_contact(mean_views, figures / "late_mean_allfields_threeview",
                   args.paper_dpi, interval_ms)

    save_mean_archive(Path(args.output_dir) / "late_mean_fields.npz", mean_views,
                      steps, times, weights, args)
    definitions = definitions_payload(caches, mean_views, weights, args,
                                      axial_crop, transverse_crop, plotfiles)
    with (Path(args.output_dir) / "mean_definitions.json").open("w") as handle:
        json.dump(definitions, handle, indent=2, sort_keys=True)
    print("Assembled {} late states over [{:.6f}, {:.6f}] ms".format(
        len(caches), interval_ms[0], interval_ms[1]), flush=True)
    print("Trapezoidal weights: {}".format(
        ", ".join("{:.9f}".format(value) for value in weights)), flush=True)


def main():
    args = parse_args()
    configure_matplotlib()
    axial_crop, transverse_crop = validate_args(args)
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    plotfiles = discover_selected_plotfiles(args.plot_root, args.steps)
    print("Selected plotfiles: {}".format(
        ", ".join("plt{:05d}".format(step) for step, _ in plotfiles)), flush=True)
    print("Requested crops: x={} m, y/z={} m".format(
        axial_crop, transverse_crop), flush=True)
    yt_module = None
    if args.mode in ("extract", "all"):
        yt_module = load_yt_module()
        extract(args, yt_module, plotfiles, axial_crop, transverse_crop)
    if args.mode in ("assemble", "all"):
        assemble(args, plotfiles, axial_crop, transverse_crop)


if __name__ == "__main__":
    main()
