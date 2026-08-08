#!/usr/bin/env python3
"""Render a camera-view 3-D density-gradient volume with the SRP STL body.

This is a schlieren-like visualization of the physical scalar
``|grad(rho)|``.  It is not a source tracer and it is not a strict optical
schlieren calculation.  A camera-side cutaway is applied only to the outer
volume so the three central retro-jet plumes remain visible.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy import ndimage
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import yt


G_REF = 72.18947715759369  # kg m^-4, frozen paper-series plt09000 P99.5
BODY_CENTRE_YZ = (0.275, 0.275)
DEFAULT_ROI_LO = (0.055, 0.135, 0.135)
DEFAULT_ROI_HI = (0.315, 0.415, 0.415)
DEFAULT_CAMERA_POSITION = (-0.42, 0.02, 0.53)
DEFAULT_CAMERA_FOCAL = (0.22, 0.275, 0.275)
DEFAULT_CAMERA_UP = (0.0, 0.0, 1.0)
OPACITY_Q = (0.0, 0.15, 0.30, 0.60, 1.0, 2.0, 4.0, 8.0)
OPACITY_ALPHA = (0.0, 0.0, 0.0020, 0.0080, 0.0200, 0.0350, 0.0120, 0.0030)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plotfile", required=True)
    parser.add_argument("--stl", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--output-name", default="srp_tri_L2_plt09000_schlieren3d_camera")
    parser.add_argument("--roi-lo", type=float, nargs=3, default=DEFAULT_ROI_LO)
    parser.add_argument("--roi-hi", type=float, nargs=3, default=DEFAULT_ROI_HI)
    parser.add_argument("--cutaway-radius", type=float, default=0.075)
    parser.add_argument("--cutaway-fade", type=float, default=0.015)
    parser.add_argument("--boundary-taper-cells", type=int, default=12)
    parser.add_argument("--mask-iterations", type=int, default=2)
    parser.add_argument("--opacity-scale", type=float, default=1.0)
    parser.add_argument("--opacity-unit-distance", type=float, default=0.0040)
    parser.add_argument("--image-size", type=int, nargs=2, default=(3200, 2000))
    parser.add_argument("--camera-position", type=float, nargs=3,
                        default=DEFAULT_CAMERA_POSITION)
    parser.add_argument("--camera-focal", type=float, nargs=3,
                        default=DEFAULT_CAMERA_FOCAL)
    parser.add_argument("--camera-up", type=float, nargs=3, default=DEFAULT_CAMERA_UP)
    parser.add_argument("--camera-view-angle", type=float, default=27.0)
    parser.add_argument("--rebuild-cache", action="store_true")
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def yt_array(values):
    if hasattr(values, "to_value"):
        return np.asarray(values.to_value(), dtype=np.float64)
    if hasattr(values, "d"):
        return np.asarray(values.d, dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


def aligned_region(dataset, target_lo, target_hi, level):
    domain_lo = yt_array(dataset.domain_left_edge)
    domain_hi = yt_array(dataset.domain_right_edge)
    base_dims = np.asarray(dataset.domain_dimensions, dtype=np.int64)
    refine = int(dataset.refine_by) ** int(level)
    dims_full = base_dims * refine
    spacing = (domain_hi - domain_lo) / dims_full
    target_lo = np.asarray(target_lo, dtype=float)
    target_hi = np.asarray(target_hi, dtype=float)
    starts = np.ceil((target_lo - domain_lo) / spacing - 1.0e-10).astype(np.int64)
    stops = np.floor((target_hi - domain_lo) / spacing + 1.0e-10).astype(np.int64)
    if np.any(starts < 0) or np.any(stops > dims_full) or np.any(stops <= starts):
        raise ValueError("Requested ROI is outside the plotfile domain")
    lo = domain_lo + starts * spacing
    dims = stops - starts
    hi = lo + dims * spacing
    return lo, hi, dims, spacing


def load_field(container, name):
    print("Loading {}".format(name), flush=True)
    result = np.asarray(container[("boxlib", name)].d, dtype=np.float32).copy()
    container.clear_data()
    return result


def gradient_magnitude(rho, spacing):
    magnitude2 = np.zeros(rho.shape, dtype=np.float32)
    for axis in range(3):
        derivative = np.gradient(rho, float(spacing[axis]), axis=axis, edge_order=2)
        np.multiply(derivative, derivative, out=derivative)
        magnitude2 += derivative
        del derivative
    np.sqrt(magnitude2, out=magnitude2)
    return magnitude2


def smoothstep(values):
    clipped = np.clip(values, 0.0, 1.0)
    return clipped * clipped * (3.0 - 2.0 * clipped)


def boundary_taper(shape, cells):
    """Separable raised-cosine-like taper that removes finite-ROI sheets."""
    weights = []
    cells = max(2, int(cells))
    for length in shape:
        index = np.arange(length, dtype=np.float32)
        distance = np.minimum(index, (length - 1) - index)
        weights.append(smoothstep(distance / float(cells)))
    return weights


def prepare_cache(args):
    yt.funcs.mylog.setLevel(30)
    dataset = yt.load(args.plotfile)
    level = int(dataset.max_level)
    lo, hi, dims, spacing = aligned_region(dataset, args.roi_lo, args.roi_hi, level)
    if not np.allclose(spacing, spacing[0], rtol=1.0e-10, atol=1.0e-15):
        raise RuntimeError("This renderer currently requires isotropic cells")
    print("L{} ROI dims={} lo={} hi={} dx={:.9g}".format(
        level, dims.tolist(), lo.tolist(), hi.tolist(), spacing[0]), flush=True)

    smooth = dataset.smoothed_covering_grid(
        level=level, left_edge=dataset.arr(lo, "code_length"), dims=dims,
        fields=[("boxlib", "Density")], use_pbar=False,
    )
    rho = load_field(smooth, "Density")
    del smooth

    categorical = dataset.covering_grid(
        level=level, left_edge=dataset.arr(lo, "code_length"), dims=dims,
        fields=[("boxlib", "sld"), ("boxlib", "ghs")], use_pbar=False,
    )
    sld = load_field(categorical, "sld")
    fluid = np.isfinite(rho) & (rho > 0.0) & np.isfinite(sld) & (sld < 0.5)
    del sld
    ghs = load_field(categorical, "ghs")
    fluid &= np.isfinite(ghs) & (ghs < 0.5)
    del ghs, categorical

    structure = ndimage.generate_binary_structure(3, 1)
    interior = ndimage.binary_erosion(
        fluid, structure=structure, iterations=max(1, args.mask_iterations),
        border_value=0,
    )
    del fluid
    grad = gradient_magnitude(rho, spacing)
    del rho
    grad[~interior] = 0.0

    # Smoothly cut away only the camera-facing OUTER field.  The full central
    # cylinder containing all three nozzles and near-field plumes is retained.
    y = lo[1] + (np.arange(dims[1]) + 0.5) * spacing[1]
    z = lo[2] + (np.arange(dims[2]) + 0.5) * spacing[2]
    dy = y[:, None] - BODY_CENTRE_YZ[0]
    dz = z[None, :] - BODY_CENTRE_YZ[1]
    radius = np.sqrt(dy * dy + dz * dz)
    camera_yz = np.asarray(args.camera_position[1:3], dtype=float) - np.asarray(
        BODY_CENTRE_YZ, dtype=float)
    camera_yz /= np.linalg.norm(camera_yz)
    camera_side_distance = dy * camera_yz[0] + dz * camera_yz[1]
    radial_gate = smoothstep((radius - args.cutaway_radius) / args.cutaway_fade)
    side_gate = smoothstep(camera_side_distance / args.cutaway_fade)
    cut_strength = radial_gate * side_gate
    grad *= (1.0 - cut_strength)[None, :, :]

    # Smoothly remove the finite ROI boundaries instead of leaving visible
    # rectangular sheets in the ray-cast image.
    taper_x, taper_y, taper_z = boundary_taper(grad.shape, args.boundary_taper_cells)
    grad *= taper_x[:, None, None]
    grad *= taper_y[None, :, None]
    grad *= taper_z[None, None, :]

    samples = grad[(grad > 0.0) & np.isfinite(grad)][::16]
    percentiles = {}
    if samples.size:
        for percentile in (50.0, 75.0, 90.0, 95.0, 99.0, 99.5, 99.9):
            percentiles["p{:g}".format(percentile)] = float(
                np.percentile(samples, percentile))

    q = np.clip(grad / G_REF, 0.0, 8.0).astype(np.float16)
    nonzero_cells = int(np.count_nonzero(q))
    del grad, interior
    origin = lo + 0.5 * spacing
    metadata = {
        "plotfile": os.path.abspath(args.plotfile),
        "simulation_time_ms": float(dataset.current_time) * 1.0e3,
        "level": level,
        "dimensions": [int(value) for value in dims],
        "origin_cell_centre_m": [float(value) for value in origin],
        "spacing_m": [float(value) for value in spacing],
        "roi_lo_m": [float(value) for value in lo],
        "roi_hi_m": [float(value) for value in hi],
        "quantity": "q = |grad(rho)| / G_ref",
        "G_ref_kg_m-4": G_REF,
        "G_ref_provenance": "paper-series plt09000 pooled XY/XZ/YZ P99.5",
        "gradient_percentiles_kg_m-4_after_cutaway": percentiles,
        "fluid_mask": "sld < 0.5 and ghs < 0.5, eroded by {} L2 cells".format(
            max(1, args.mask_iterations)),
        "cutaway": {
            "description": "camera-side outer volume smoothly faded; central cylinder retained",
            "central_radius_m": float(args.cutaway_radius),
            "fade_width_m": float(args.cutaway_fade),
            "camera_side_normal_yz": [float(value) for value in camera_yz],
        },
        "boundary_taper_cells": int(args.boundary_taper_cells),
        "nonzero_volume_cells": nonzero_cells,
        "classification_is_jet_tracer": False,
        "scientific_label": "schlieren-like density-gradient volume rendering",
        "caveat": "not strict optical schlieren and not source-resolved without Cjet",
        "sampler": "yt smoothed_covering_grid Density; categorical covering_grid IBM IDs",
    }
    Path(args.cache).parent.mkdir(parents=True, exist_ok=True)
    temporary = args.cache + ".tmp.npz"
    np.savez_compressed(
        temporary, q=q, origin=origin.astype(np.float64),
        spacing=spacing.astype(np.float64), metadata=np.asarray(json.dumps(metadata)),
    )
    os.replace(temporary, args.cache)
    print("Cached {} ({:,} nonzero cells)".format(args.cache, nonzero_cells), flush=True)
    return metadata


def vtk_image(q, origin, spacing):
    image = vtk.vtkImageData()
    image.SetDimensions(*q.shape)
    image.SetOrigin(*[float(value) for value in origin])
    image.SetSpacing(*[float(value) for value in spacing])
    scalars = numpy_to_vtk(q.ravel(order="F"), deep=True, array_type=vtk.VTK_FLOAT)
    scalars.SetName("density_gradient_ratio")
    image.GetPointData().SetScalars(scalars)
    return image


def text_actor(text, x, y, size, color=(0.10, 0.11, 0.13), bold=False):
    actor = vtk.vtkTextActor()
    actor.SetInput(text)
    actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    actor.GetPositionCoordinate().SetValue(x, y)
    prop = actor.GetTextProperty()
    prop.SetFontFamilyToArial()
    prop.SetFontSize(size)
    prop.SetColor(*color)
    prop.SetBold(bool(bold))
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
    with np.load(args.cache, allow_pickle=False) as cache:
        q = np.asarray(cache["q"], dtype=np.float32)
        origin = np.asarray(cache["origin"], dtype=np.float64)
        spacing = np.asarray(cache["spacing"], dtype=np.float64)
        metadata = json.loads(str(cache["metadata"]))

    image = vtk_image(q, origin, spacing)
    mapper = vtk.vtkSmartVolumeMapper()
    mapper.SetInputData(image)
    mapper.SetBlendModeToComposite()
    if hasattr(mapper, "SetSampleDistance"):
        mapper.SetSampleDistance(0.60 * float(spacing[0]))
    if hasattr(mapper, "AutoAdjustSampleDistancesOff"):
        mapper.AutoAdjustSampleDistancesOff()

    color = vtk.vtkColorTransferFunction()
    for scalar in OPACITY_Q:
        # Near-black ink on a white background; a subtle blue component keeps
        # depth accumulation from looking like a flat binary silhouette.
        fraction = min(1.0, scalar / 8.0)
        color.AddRGBPoint(scalar, 0.055 + 0.045 * fraction,
                         0.065 + 0.050 * fraction, 0.085 + 0.070 * fraction)
    opacity = vtk.vtkPiecewiseFunction()
    for scalar, alpha in zip(OPACITY_Q, OPACITY_ALPHA):
        opacity.AddPoint(scalar, min(0.95, args.opacity_scale * alpha))
    volume_property = vtk.vtkVolumeProperty()
    volume_property.SetColor(color)
    volume_property.SetScalarOpacity(opacity)
    volume_property.SetInterpolationTypeToLinear()
    volume_property.ShadeOff()
    volume_property.SetScalarOpacityUnitDistance(args.opacity_unit_distance)
    volume = vtk.vtkVolume()
    volume.SetMapper(mapper)
    volume.SetProperty(volume_property)

    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(args.stl)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(stl_reader.GetOutputPort())
    normals.SetFeatureAngle(50.0)
    normals.ConsistencyOn()
    normals.SplittingOff()
    body_mapper = vtk.vtkPolyDataMapper()
    body_mapper.SetInputConnection(normals.GetOutputPort())
    body_actor = vtk.vtkActor()
    body_actor.SetMapper(body_mapper)
    body_actor.GetProperty().SetColor(0.74, 0.76, 0.79)
    body_actor.GetProperty().SetInterpolationToPhong()
    body_actor.GetProperty().SetAmbient(0.48)
    body_actor.GetProperty().SetDiffuse(0.48)
    body_actor.GetProperty().SetSpecular(0.18)
    body_actor.GetProperty().SetSpecularPower(28.0)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.AddVolume(volume)
    renderer.AddActor(body_actor)
    if hasattr(renderer, "SetUseFXAA"):
        renderer.SetUseFXAA(True)

    camera = vtk.vtkCamera()
    camera.SetPosition(*args.camera_position)
    camera.SetFocalPoint(*args.camera_focal)
    camera.SetViewUp(*args.camera_up)
    camera.SetViewAngle(args.camera_view_angle)
    renderer.SetActiveCamera(camera)

    silhouette = vtk.vtkPolyDataSilhouette()
    silhouette.SetInputConnection(normals.GetOutputPort())
    silhouette.SetCamera(camera)
    silhouette.SetEnableFeatureAngle(0)
    silhouette_mapper = vtk.vtkPolyDataMapper()
    silhouette_mapper.SetInputConnection(silhouette.GetOutputPort())
    silhouette_actor = vtk.vtkActor()
    silhouette_actor.SetMapper(silhouette_mapper)
    silhouette_actor.GetProperty().SetColor(0.08, 0.09, 0.11)
    silhouette_actor.GetProperty().SetLineWidth(1.6)
    renderer.AddActor(silhouette_actor)

    key = vtk.vtkLight()
    key.SetPosition(-0.25, -0.05, 0.70)
    key.SetFocalPoint(*args.camera_focal)
    key.SetIntensity(0.95)
    renderer.AddLight(key)
    fill = vtk.vtkLight()
    fill.SetPosition(0.62, 0.60, 0.18)
    fill.SetFocalPoint(*args.camera_focal)
    fill.SetIntensity(0.42)
    renderer.AddLight(fill)

    width, height = args.image_size
    title_size = max(18, int(round(27.0 * width / 3200.0)))
    note_size = max(13, int(round(17.0 * width / 3200.0)))
    renderer.AddActor2D(text_actor(
        "3-D schlieren-like density-gradient volume  |grad(rho)|",
        0.025, 0.950, title_size, bold=True,
    ))
    renderer.AddActor2D(text_actor(
        "plt09000  t={:.6f} ms   |   complete STL geometry".format(
            metadata["simulation_time_ms"]),
        0.025, 0.914, note_size,
    ))
    renderer.AddActor2D(text_actor(
        "outer camera-side volume clipped for plume visibility; central jet volume retained",
        0.025, 0.035, note_size, color=(0.25, 0.27, 0.30),
    ))
    renderer.AddActor2D(text_actor(
        "density-gradient visualization, not a transported jet tracer",
        0.025, 0.012, max(11, note_size - 2), color=(0.38, 0.40, 0.43),
    ))

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    window.SetSize(int(width), int(height))
    window.AddRenderer(renderer)
    renderer.ResetCameraClippingRange()
    window.Render()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / (args.output_name + ".png")
    write_png(window, str(png_path))

    metadata["render"] = {
        "output_png": str(png_path),
        "image_size": [int(width), int(height)],
        "camera_position_m": [float(value) for value in args.camera_position],
        "camera_focal_m": [float(value) for value in args.camera_focal],
        "camera_up": [float(value) for value in args.camera_up],
        "camera_view_angle_deg": float(args.camera_view_angle),
        "blend_mode": "front-to-back composite DVR",
        "sample_distance_m": 0.60 * float(spacing[0]),
        "opacity_unit_distance_m": float(args.opacity_unit_distance),
        "opacity_scale": float(args.opacity_scale),
        "opacity_q_nodes": list(OPACITY_Q),
        "opacity_alpha_nodes": [float(args.opacity_scale * value)
                                for value in OPACITY_ALPHA],
        "body": "complete opaque STL with Phong shading and silhouette",
    }
    metadata["script"] = os.path.abspath(__file__)
    metadata["script_sha256"] = sha256(__file__)
    metadata["stl"] = os.path.abspath(args.stl)
    metadata["stl_sha256"] = sha256(args.stl)
    summary_path = output_dir / (args.output_name + "_summary.json")
    with summary_path.open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
    print("Wrote {}".format(png_path), flush=True)
    print("Wrote {}".format(summary_path), flush=True)


def main():
    args = arguments()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    if args.rebuild_cache or not os.path.exists(args.cache):
        prepare_cache(args)
    else:
        print("Using cache {}".format(args.cache), flush=True)
    render(args)


if __name__ == "__main__":
    main()
