#!/usr/bin/env python3
"""Clean VTK rendering of the plt09000 tri-nozzle jet core.

The plotfile has no passive tracer.  We therefore retain only the weak
high-total-pressure region reachable from geometric nozzle seeds, render a
much stricter high-total-pressure core, and label the result as a proxy.
The preprocessed L2 mask is cached so camera/style iterations do not reread
the 25 GiB AMReX plotfile.
"""

import argparse
import json
import os

import numpy as np
from PIL import Image
from scipy import ndimage
import yt
import vtk
from vtk.util.numpy_support import numpy_to_vtk


GAMMA = 1.4
P0_INF = 188222.90041860205
U_INF = 743.3947269116185
NOZZLE_X = 0.299
NOZZLE_R = 0.0028
NOZZLES = (
    ("jet 0", 0.307, 0.275, (1.00, 0.18, 0.12)),
    ("jet 1", 0.259, 0.302, (1.00, 0.52, 0.05)),
    ("jet 2", 0.259, 0.248, (1.00, 0.84, 0.08)),
)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plotfile", required=True)
    parser.add_argument("--stl", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--envelope-p0", type=float, default=3.5e5)
    parser.add_argument("--core-p0", type=float, default=1.5e6)
    parser.add_argument("--rebuild-cache", action="store_true")
    return parser.parse_args()


def load_field(grid, name):
    print("Loading {}".format(name), flush=True)
    result = np.asarray(grid[("boxlib", name)].d, dtype=np.float32).copy()
    grid.clear_data()
    return result


def sample(array, point, lo, dx):
    coordinates = (np.asarray(point, dtype=float) - lo) / dx - 0.5
    return float(ndimage.map_coordinates(
        array, coordinates[:, None], order=1, mode="constant", cval=np.nan,
        prefilter=False,
    )[0])


def sample_velocity(velocity, point, lo, dx):
    return np.asarray([sample(field, point, lo, dx) for field in velocity])


def trace(seed, velocity, solid, lo, hi, dx, step, maximum_steps):
    point = np.asarray(seed, dtype=float)
    path = [point.copy()]
    margin = 1.5 * dx
    for _ in range(maximum_steps):
        if np.any(point <= lo + margin) or np.any(point >= hi - margin):
            break
        if sample(solid, point, lo, dx) > 0.48:
            break
        vector = sample_velocity(velocity, point, lo, dx)
        speed = float(np.linalg.norm(vector))
        if not np.all(np.isfinite(vector)) or speed < 5.0:
            break
        midpoint = point + 0.5 * step * vector / speed
        middle_velocity = sample_velocity(velocity, midpoint, lo, dx)
        middle_speed = float(np.linalg.norm(middle_velocity))
        if not np.all(np.isfinite(middle_velocity)) or middle_speed < 5.0:
            break
        point = point + step * middle_velocity / middle_speed
        path.append(point.copy())
    return np.asarray(path, dtype=np.float32)


def make_streamlines(velocity, solid, lo, hi, dx):
    jet_paths = []
    ring = [(0.0, 0.0)] + [
        (0.00145 * np.cos(a), 0.00145 * np.sin(a))
        for a in np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)
    ]
    for nozzle_id, (_, y0, z0, _) in enumerate(NOZZLES):
        for dy, dz in ring:
            path = trace(
                (NOZZLE_X - 0.0032, y0 + dy, z0 + dz), velocity, solid,
                lo, hi, dx, step=0.00075, maximum_steps=360,
            )
            if len(path) >= 12:
                jet_paths.append((nozzle_id, path))

    free_paths = []
    for radius in (0.105, 0.128):
        for angle in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
            y0 = 0.275 + radius * np.cos(angle)
            z0 = 0.275 + radius * np.sin(angle)
            if not (lo[1] < y0 < hi[1] and lo[2] < z0 < hi[2]):
                continue
            path = trace(
                (lo[0] + 0.0025, y0, z0), velocity, solid,
                lo, hi, dx, step=0.0013, maximum_steps=260,
            )
            if len(path) >= 20:
                free_paths.append(path)
    return jet_paths, free_paths


def block_mean(array):
    nx, ny, nz = (size // 2 for size in array.shape)
    view = array[:2 * nx, :2 * ny, :2 * nz].reshape(nx, 2, ny, 2, nz, 2)
    return view.mean(axis=(1, 3, 5)).astype(np.float32)


def pack_paths(paths, with_ids=False):
    if not paths:
        return np.empty((0, 3), np.float32), np.empty(0, np.int32), np.empty(0, np.int16)
    if with_ids:
        identifiers = np.asarray([item[0] for item in paths], dtype=np.int16)
        arrays = [item[1] for item in paths]
    else:
        identifiers = np.zeros(len(paths), dtype=np.int16)
        arrays = paths
    lengths = np.asarray([len(path) for path in arrays], dtype=np.int32)
    return np.concatenate(arrays, axis=0), lengths, identifiers


def unpack_paths(points, lengths, identifiers=None):
    paths = []
    start = 0
    for index, length in enumerate(lengths):
        path = points[start:start + int(length)]
        start += int(length)
        if identifiers is None:
            paths.append(path)
        else:
            paths.append((int(identifiers[index]), path))
    return paths


def prepare_cache(args):
    yt.set_log_level(50)
    dataset = yt.load(args.plotfile)
    level = int(dataset.max_level)
    base = np.asarray(dataset.domain_dimensions, dtype=int)
    width = np.asarray(dataset.domain_right_edge.d - dataset.domain_left_edge.d)
    dx = float(width[0] / (base[0] * (2 ** level)))
    lo = np.asarray([0.060, 0.140, 0.140], dtype=float)
    target_hi = np.asarray([0.305, 0.410, 0.410], dtype=float)
    dimensions = np.rint((target_hi - lo) / dx).astype(int)
    hi = lo + dimensions * dx
    print("L{} ROI {} dx={:.9g}".format(level, dimensions.tolist(), dx), flush=True)

    grid = dataset.covering_grid(
        level=level, left_edge=dataset.arr(lo, "code_length"), dims=dimensions,
    )
    rho = load_field(grid, "Density")
    pressure = load_field(grid, "pressure")
    ux = load_field(grid, "x_velocity")
    uy = load_field(grid, "y_velocity")
    uz = load_field(grid, "z_velocity")
    solid = load_field(grid, "sld")
    ghost = load_field(grid, "ghs")
    del grid

    fluid = (
        (solid < 0.5) & (ghost < 0.5) & np.isfinite(rho) & np.isfinite(pressure)
        & (rho > 0.0) & (pressure > 0.0)
    )
    speed2 = ux * ux + uy * uy + uz * uz
    mach2 = speed2 * rho / np.maximum(GAMMA * pressure, 1.0e-30)
    mach = np.sqrt(np.maximum(mach2, 0.0)).astype(np.float32)
    total_pressure = (
        pressure * np.power(1.0 + 0.2 * mach2, 3.5)
    ).astype(np.float32)

    weak = fluid & (total_pressure > args.envelope_p0)
    seeds = np.zeros(weak.shape, dtype=bool)
    seed_counts = {}
    x_centres = lo[0] + (np.arange(dimensions[0]) + 0.5) * dx
    y_centres = lo[1] + (np.arange(dimensions[1]) + 0.5) * dx
    z_centres = lo[2] + (np.arange(dimensions[2]) + 0.5) * dx
    ix = np.where((x_centres >= NOZZLE_X - 0.006) & (x_centres <= NOZZLE_X))[0]
    for name, y0, z0, _ in NOZZLES:
        jy = np.where(np.abs(y_centres - y0) <= 1.6 * NOZZLE_R)[0]
        kz = np.where(np.abs(z_centres - z0) <= 1.6 * NOZZLE_R)[0]
        local_y = y_centres[jy][:, None]
        local_z = z_centres[kz][None, :]
        disk = (local_y - y0) ** 2 + (local_z - z0) ** 2 <= (1.5 * NOZZLE_R) ** 2
        local = np.zeros((len(ix), len(jy), len(kz)), dtype=bool)
        local[:] = disk[None, :, :]
        target = np.ix_(ix, jy, kz)
        local &= weak[target] & (total_pressure[target] > 9.0e5)
        seeds[target] |= local
        seed_counts[name] = int(np.count_nonzero(local))
    if any(value == 0 for value in seed_counts.values()):
        raise RuntimeError("Missing nozzle seed cells: {}".format(seed_counts))

    selected = ndimage.binary_propagation(
        seeds, structure=ndimage.generate_binary_structure(3, 1), mask=weak,
    )
    core = selected & (total_pressure > args.core_p0)
    if not np.any(selected) or not np.any(core):
        raise RuntimeError("Empty connected envelope/core")
    boundary_contacts = {
        "x_lo": int(np.count_nonzero(selected[0])),
        "x_hi": int(np.count_nonzero(selected[-1])),
        "y_lo": int(np.count_nonzero(selected[:, 0, :])),
        "y_hi": int(np.count_nonzero(selected[:, -1, :])),
        "z_lo": int(np.count_nonzero(selected[:, :, 0])),
        "z_hi": int(np.count_nonzero(selected[:, :, -1])),
    }
    if any(boundary_contacts.values()):
        raise RuntimeError("Connected mask touches ROI: {}".format(boundary_contacts))

    jet_paths, free_paths = make_streamlines((ux, uy, uz), solid, lo, hi, dx)
    jet_points, jet_lengths, jet_ids = pack_paths(jet_paths, with_ids=True)
    free_points, free_lengths, _ = pack_paths(free_paths)

    mach_for_average = mach.copy()
    mach_for_average[~fluid] = 0.0
    valid_average = block_mean(fluid.astype(np.float32))
    mach_surface = block_mean(mach_for_average) / np.maximum(valid_average, 1.0e-6)
    envelope_surface = block_mean(selected.astype(np.float32))
    core_surface = block_mean(core.astype(np.float32))
    surface_dx = 2.0 * dx
    surface_origin = lo + 0.5 * surface_dx

    metadata = {
        "plotfile": os.path.abspath(args.plotfile),
        "simulation_time_ms": float(dataset.current_time) * 1.0e3,
        "classification_is_exact_tracer": False,
        "classification": "6-neighbour nozzle-seeded high-local-total-pressure proxy",
        "envelope_p0_Pa": float(args.envelope_p0),
        "core_p0_Pa": float(args.core_p0),
        "freestream_p0_Pa": P0_INF,
        "level": level,
        "fine_dimensions": dimensions.tolist(),
        "fine_spacing_m": dx,
        "roi_lo_m": lo.tolist(),
        "roi_hi_m": hi.tolist(),
        "surface_origin_m": surface_origin.tolist(),
        "surface_spacing_m": surface_dx,
        "surface_dimensions": list(core_surface.shape),
        "seed_counts": seed_counts,
        "boundary_contacts": boundary_contacts,
        "envelope_cells": int(np.count_nonzero(selected)),
        "core_cells": int(np.count_nonzero(core)),
        "jet_streamlines": len(jet_paths),
        "freestream_streamlines": len(free_paths),
        "surface_warning": "2x block reconstruction is qualitative; do not use for exact area.",
    }
    np.savez_compressed(
        args.cache,
        envelope=envelope_surface,
        core=core_surface,
        mach=mach_surface.astype(np.float32),
        origin=surface_origin.astype(np.float64),
        spacing=np.asarray(surface_dx, dtype=np.float64),
        jet_points=jet_points,
        jet_lengths=jet_lengths,
        jet_ids=jet_ids,
        free_points=free_points,
        free_lengths=free_lengths,
        metadata=np.asarray(json.dumps(metadata)),
    )
    print("Cached {}".format(args.cache), flush=True)
    return metadata


def vtk_image(indicator, mach, origin, spacing, indicator_name):
    image = vtk.vtkImageData()
    image.SetDimensions(*indicator.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(spacing, spacing, spacing)
    vtk_indicator = numpy_to_vtk(
        indicator.ravel(order="F"), deep=True, array_type=vtk.VTK_FLOAT,
    )
    vtk_indicator.SetName(indicator_name)
    image.GetPointData().SetScalars(vtk_indicator)
    if mach is not None:
        vtk_mach = numpy_to_vtk(mach.ravel(order="F"), deep=True, array_type=vtk.VTK_FLOAT)
        vtk_mach.SetName("Mach")
        image.GetPointData().AddArray(vtk_mach)
    return image


def smooth_surface(image, scalar_name):
    image.GetPointData().SetActiveScalars(scalar_name)
    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputData(image)
    contour.SetValue(0, 0.5)
    contour.ComputeScalarsOn()
    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputConnection(contour.GetOutputPort())
    smoother.SetNumberOfIterations(14)
    smoother.BoundarySmoothingOff()
    smoother.FeatureEdgeSmoothingOff()
    smoother.SetPassBand(0.10)
    smoother.NonManifoldSmoothingOn()
    smoother.NormalizeCoordinatesOn()
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(smoother.GetOutputPort())
    normals.SetFeatureAngle(60.0)
    normals.ConsistencyOn()
    normals.SplittingOff()
    normals.Update()
    return normals


def turbo_lookup():
    from matplotlib import cm
    table = vtk.vtkLookupTable()
    table.SetNumberOfTableValues(256)
    table.SetRange(0.0, 11.0)
    table.Build()
    cmap = cm.get_cmap("turbo")
    for index in range(256):
        r, g, b, a = cmap(index / 255.0)
        table.SetTableValue(index, r, g, b, a)
    return table


def polyline_data(paths):
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    next_id = 0
    for path in paths:
        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(len(path))
        for local_id, point in enumerate(path):
            points.InsertNextPoint(*[float(value) for value in point])
            line.GetPointIds().SetId(local_id, next_id)
            next_id += 1
        lines.InsertNextCell(line)
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetLines(lines)
    return poly


def tube_actor(paths, radius, color, opacity=1.0):
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(polyline_data(paths))
    tube.SetRadius(radius)
    tube.SetNumberOfSides(10)
    tube.CappingOn()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().SetSpecular(0.25)
    return actor


def text_actor(text, x, y, size, color=(0.92, 0.95, 1.0)):
    actor = vtk.vtkTextActor()
    actor.SetInput(text)
    actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    actor.GetPositionCoordinate().SetValue(x, y)
    prop = actor.GetTextProperty()
    prop.SetFontFamilyToArial()
    prop.SetFontSize(size)
    prop.SetColor(*color)
    prop.SetShadow(True)
    return actor


def write_png(window, path):
    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.ReadFrontBufferOff()
    capture.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(path)
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()


def render(args):
    cache = np.load(args.cache, allow_pickle=False)
    metadata = json.loads(str(cache["metadata"]))
    envelope = cache["envelope"]
    core = cache["core"]
    mach = cache["mach"]
    origin = cache["origin"]
    spacing = float(cache["spacing"])
    jet_paths = unpack_paths(cache["jet_points"], cache["jet_lengths"], cache["jet_ids"])
    free_paths = unpack_paths(cache["free_points"], cache["free_lengths"])

    envelope_image = vtk_image(envelope, None, origin, spacing, "Envelope")
    core_image = vtk_image(core, mach, origin, spacing, "Core")
    envelope_surface = smooth_surface(envelope_image, "Envelope")
    core_surface = smooth_surface(core_image, "Core")

    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(core_surface.GetOutputPort())
    probe.SetSourceData(core_image)
    probe.Update()
    lookup = turbo_lookup()
    core_mapper = vtk.vtkPolyDataMapper()
    core_mapper.SetInputConnection(probe.GetOutputPort())
    core_mapper.SetLookupTable(lookup)
    core_mapper.SetScalarRange(0.0, 11.0)
    core_mapper.SetScalarModeToUsePointFieldData()
    core_mapper.SelectColorArray("Mach")
    core_mapper.ScalarVisibilityOn()
    core_actor = vtk.vtkActor()
    core_actor.SetMapper(core_mapper)
    core_actor.GetProperty().SetInterpolationToPhong()
    core_actor.GetProperty().SetSpecular(0.32)
    core_actor.GetProperty().SetSpecularPower(24.0)

    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(args.stl)
    body_normals = vtk.vtkPolyDataNormals()
    body_normals.SetInputConnection(stl_reader.GetOutputPort())
    body_normals.SetFeatureAngle(50.0)
    body_mapper = vtk.vtkPolyDataMapper()
    body_mapper.SetInputConnection(body_normals.GetOutputPort())
    body_actor = vtk.vtkActor()
    body_actor.SetMapper(body_mapper)
    body_actor.GetProperty().SetColor(0.63, 0.67, 0.72)
    body_actor.GetProperty().SetInterpolationToPhong()
    body_actor.GetProperty().SetSpecular(0.50)
    body_actor.GetProperty().SetSpecularPower(35.0)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.018, 0.026, 0.050)
    renderer.SetBackground2(0.075, 0.105, 0.155)
    renderer.GradientBackgroundOn()
    renderer.SetUseDepthPeeling(True)
    renderer.SetMaximumNumberOfPeels(100)
    renderer.SetOcclusionRatio(0.05)
    renderer.AddActor(body_actor)
    renderer.AddActor(core_actor)

    camera = vtk.vtkCamera()
    renderer.SetActiveCamera(camera)
    silhouette = vtk.vtkPolyDataSilhouette()
    silhouette.SetInputConnection(envelope_surface.GetOutputPort())
    silhouette.SetCamera(camera)
    silhouette.SetEnableFeatureAngle(1)
    silhouette.SetFeatureAngle(42.0)
    silhouette_mapper = vtk.vtkPolyDataMapper()
    silhouette_mapper.SetInputConnection(silhouette.GetOutputPort())
    silhouette_actor = vtk.vtkActor()
    silhouette_actor.SetMapper(silhouette_mapper)
    silhouette_actor.GetProperty().SetColor(1.0, 0.72, 0.16)
    silhouette_actor.GetProperty().SetLineWidth(2.2)
    renderer.AddActor(silhouette_actor)

    for nozzle_id, (_, _, _, color) in enumerate(NOZZLES):
        paths = [path for path_id, path in jet_paths if path_id == nozzle_id]
        if paths:
            renderer.AddActor(tube_actor(paths, 0.00048, color, 0.98))
    if free_paths:
        renderer.AddActor(tube_actor(free_paths, 0.00024, (0.10, 0.55, 1.0), 0.62))

    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lookup)
    scalar_bar.SetTitle("jet-core Mach")
    scalar_bar.SetNumberOfLabels(6)
    scalar_bar.SetPosition(0.865, 0.12)
    scalar_bar.SetWidth(0.105)
    scalar_bar.SetHeight(0.36)
    scalar_bar.GetTitleTextProperty().SetColor(0.93, 0.95, 1.0)
    scalar_bar.GetLabelTextProperty().SetColor(0.93, 0.95, 1.0)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(text_actor(
        "plt09000  t={:.3f} ms | jet core p0 > {:.2f} MPa, colored by Mach".format(
            metadata["simulation_time_ms"], metadata["core_p0_Pa"] / 1.0e6),
        0.025, 0.945, 20,
    ))
    renderer.AddActor2D(text_actor(
        "gold outline: nozzle-connected p0 > {:.2f} MPa | warm tubes: jets | blue tubes: freestream".format(
            metadata["envelope_p0_Pa"] / 1.0e6),
        0.025, 0.905, 15, (0.82, 0.86, 0.92),
    ))
    renderer.AddActor2D(text_actor(
        "proxy only: exact source separation requires a transported Cjet tracer",
        0.025, 0.030, 14, (0.70, 0.74, 0.81),
    ))

    key = vtk.vtkLight()
    key.SetPosition(-0.25, 0.05, 0.65)
    key.SetFocalPoint(0.22, 0.275, 0.275)
    key.SetIntensity(1.05)
    renderer.AddLight(key)
    fill = vtk.vtkLight()
    fill.SetPosition(0.45, 0.55, 0.18)
    fill.SetFocalPoint(0.22, 0.275, 0.275)
    fill.SetIntensity(0.45)
    renderer.AddLight(fill)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    window.SetSize(1800, 1100)
    window.AddRenderer(renderer)

    prefix = "srp_tri_L2_plt09000_jet3d_vtk"
    upstream_path = os.path.join(args.output_dir, prefix + "_upstream.png")
    side_path = os.path.join(args.output_dir, prefix + "_side.png")
    combined_path = os.path.join(args.output_dir, prefix + "_combined.png")
    core_mesh_path = os.path.join(args.output_dir, prefix + "_core.vtp")
    envelope_mesh_path = os.path.join(args.output_dir, prefix + "_envelope.vtp")

    camera.SetPosition(-0.42, 0.02, 0.53)
    camera.SetFocalPoint(0.18, 0.275, 0.275)
    camera.SetViewUp(0.0, 0.0, 1.0)
    camera.SetViewAngle(28.0)
    renderer.ResetCameraClippingRange()
    window.Render()
    write_png(window, upstream_path)

    camera.SetPosition(-0.12, 0.64, 0.49)
    camera.SetFocalPoint(0.19, 0.275, 0.275)
    camera.SetViewUp(0.0, 0.0, 1.0)
    camera.SetViewAngle(27.0)
    renderer.ResetCameraClippingRange()
    silhouette.Modified()
    window.Render()
    write_png(window, side_path)

    with Image.open(upstream_path) as first, Image.open(side_path) as second:
        combined = Image.new("RGB", (first.width + second.width, max(first.height, second.height)), "white")
        combined.paste(first.convert("RGB"), (0, 0))
        combined.paste(second.convert("RGB"), (first.width, 0))
        combined.save(combined_path, quality=95)

    core_writer = vtk.vtkXMLPolyDataWriter()
    core_writer.SetFileName(core_mesh_path)
    core_writer.SetInputConnection(probe.GetOutputPort())
    core_writer.SetDataModeToBinary()
    core_writer.Write()
    envelope_writer = vtk.vtkXMLPolyDataWriter()
    envelope_writer.SetFileName(envelope_mesh_path)
    envelope_writer.SetInputConnection(envelope_surface.GetOutputPort())
    envelope_writer.SetDataModeToBinary()
    envelope_writer.Write()

    probe.Update()
    envelope_surface.Update()
    metadata["render_outputs"] = {
        "upstream_png": upstream_path,
        "side_png": side_path,
        "combined_png": combined_path,
        "core_vtp": core_mesh_path,
        "envelope_vtp": envelope_mesh_path,
    }
    metadata["core_surface"] = {
        "points": int(probe.GetOutput().GetNumberOfPoints()),
        "polygons": int(probe.GetOutput().GetNumberOfPolys()),
    }
    metadata["envelope_surface"] = {
        "points": int(envelope_surface.GetOutput().GetNumberOfPoints()),
        "polygons": int(envelope_surface.GetOutput().GetNumberOfPolys()),
    }
    summary_path = os.path.join(args.output_dir, prefix + "_summary.json")
    with open(summary_path, "w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
    print(json.dumps(metadata, indent=2, sort_keys=True), flush=True)


def main():
    args = arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    if args.rebuild_cache or not os.path.exists(args.cache):
        prepare_cache(args)
    else:
        print("Using cache {}".format(args.cache), flush=True)
    render(args)


if __name__ == "__main__":
    main()
