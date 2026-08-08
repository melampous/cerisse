#!/usr/bin/env python3
"""Scale, roll and translate a Test 1853 STL for the raw-coordinate IBM reader.

The CAD/STL masters are in millimetres at NASA TSP x=0.  Cerisse's IBM reader
does not parse STL units and initializes static transforms to identity, so an
SI computation normally needs a pre-transformed metre STL.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import trimesh


HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"
MODELS = {
    "single": "test1853_run165_single",
    "tri_n3at120": "test1853_run262_263_tri_n3at120",
    "tri_n3at240": "test1853_run262_263_tri_n3at240",
    "quad_n3at120": "test1853_quad_n3at120",
    "quad_n3at240": "test1853_quad_n3at240",
    "ring4": "test1853_four_peripheral_ring90",
    "single_surface_normal": "test1853_single_peripheral_surface_normal",
}


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(HERE))
    except ValueError:
        return str(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def transform_matrix(scale: float, roll_deg: float, translation: list[float]) -> np.ndarray:
    phi = math.radians(roll_deg)
    c, s = math.cos(phi), math.sin(phi)
    matrix = np.array(
        [
            [scale, 0.0, 0.0, translation[0]],
            [0.0, scale * c, -scale * s, translation[1]],
            [0.0, scale * s, scale * c, translation[2]],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    return matrix


def transform_point(point_mm: list[float], matrix: np.ndarray) -> list[float]:
    homogeneous = np.array([*point_mm, 1.0], dtype=float)
    return (matrix @ homogeneous)[:3].tolist()


def write_transformed(source: Path, destination: Path, matrix: np.ndarray) -> dict[str, Any]:
    mesh = trimesh.load_mesh(source, process=False)
    if not isinstance(mesh, trimesh.Trimesh):
        raise RuntimeError(f"{source}: expected one triangular mesh")
    mesh.apply_transform(matrix)
    destination.parent.mkdir(parents=True, exist_ok=True)
    data = mesh.export(file_type="stl")
    if not isinstance(data, bytes):
        data = data.encode("utf-8")
    destination.write_bytes(data)
    checked = trimesh.load_mesh(destination, process=True)
    if not isinstance(checked, trimesh.Trimesh):
        raise RuntimeError(f"{destination}: transformed output is not one mesh")
    return {
        "source": display_path(source),
        "destination": str(destination),
        "source_sha256": sha256(source),
        "destination_sha256": sha256(destination),
        "triangles": int(len(checked.faces)),
        "bounds_m": checked.bounds.tolist(),
        "watertight": bool(checked.is_watertight),
        "winding_consistent": bool(checked.is_winding_consistent),
        "is_volume": bool(checked.is_volume),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configuration", choices=sorted(MODELS), required=True)
    parser.add_argument(
        "--input-stl",
        type=Path,
        help=(
            "closed millimetre/body-frame STL to transform; defaults to the "
            "standard outputs/stl model selected by --configuration"
        ),
    )
    parser.add_argument(
        "--patch-dir",
        type=Path,
        help=(
            "optional directory of matching millimetre/body-frame patch STLs; "
            "custom --input-stl files omit patches unless this is supplied"
        ),
    )
    parser.add_argument(
        "--audit-file",
        type=Path,
        help="build_audit.json matching a custom input STL; defaults to outputs/reports",
    )
    parser.add_argument("--run", choices=[165, 247, 262, 263, 307, 315], type=int)
    parser.add_argument(
        "--roll-deg",
        type=float,
        help="roll about +x; default is the selected run phi, otherwise zero",
    )
    parser.add_argument(
        "--tsp-m",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=[0.0, 0.0, 0.0],
        help="world-coordinate location of NASA theoretical sharp point in metres",
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS / "ibm_prepared")
    parser.add_argument("--name", help="output stem; generated automatically by default")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    parameters = json.loads((HERE / "geometry_parameters.yaml").read_text(encoding="utf-8"))
    audit_path = (
        args.audit_file.resolve()
        if args.audit_file is not None
        else OUTPUTS / "reports" / "build_audit.json"
    )
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    stem = MODELS[args.configuration]
    if args.run == 165 and args.configuration != "single":
        raise SystemExit("Run 165 requires --configuration single")
    if args.run in (247, 262, 263) and not args.configuration.startswith("tri"):
        raise SystemExit("Runs 247/262/263 require a tri configuration")
    if args.run in (307, 315) and not args.configuration.startswith("quad"):
        raise SystemExit("Runs 307/315 require a quad configuration")
    if args.configuration.startswith("quad") and args.run not in (None, 307, 315):
        raise SystemExit("Quad configurations require Run 307/315 or no --run")
    if args.configuration in ("ring4", "single_surface_normal") and args.run is not None:
        raise SystemExit(
            f"The parametric {args.configuration} configuration has no official run ID"
        )
    roll = args.roll_deg
    if roll is None:
        roll = float(parameters["runs"][str(args.run)]["roll_phi_deg"]) if args.run else 0.0
    matrix = transform_matrix(1.0e-3, roll, list(args.tsp_m))
    suffix = f"run{args.run}" if args.run else "master"
    source_stl = (
        args.input_stl.resolve()
        if args.input_stl is not None
        else (OUTPUTS / "stl" / f"{stem}.stl").resolve()
    )
    if not source_stl.is_file():
        raise SystemExit(f"Input STL does not exist: {source_stl}")
    name = args.name or f"{source_stl.stem}_{suffix}_SI"
    output_dir = args.output_dir.resolve() / name
    full = write_transformed(
        source_stl, output_dir / f"{name}.stl", matrix
    )
    if not (full["watertight"] and full["winding_consistent"] and full["is_volume"]):
        raise RuntimeError(f"Transformed closed STL failed topology checks: {full}")
    patches: dict[str, Any] = {}
    if args.patch_dir is not None:
        patch_dir = args.patch_dir.resolve()
        if not patch_dir.is_dir():
            raise SystemExit(f"Patch directory does not exist: {patch_dir}")
    elif args.input_stl is None:
        patch_dir = OUTPUTS / "patches" / stem
    else:
        patch_dir = None
    if patch_dir is not None:
        for source in sorted(patch_dir.glob("*.stl")):
            patches[source.name] = write_transformed(
                source, output_dir / "patches" / source.name, matrix
            )
    nozzles = []
    for nozzle in audit[stem]["nozzles"]:
        center_mm = list(
            nozzle.get(
                "throat_center_mm",
                [
                    float(nozzle["throat_x_mm"]),
                    float(nozzle["axis_y_mm"]),
                    float(nozzle["axis_z_mm"]),
                ],
            )
        )
        body_normal = np.asarray(
            nozzle.get("throat_plane_normal_solid_to_fluid", [-1.0, 0.0, 0.0]),
            dtype=float,
        )
        phi = math.radians(roll)
        c, s = math.cos(phi), math.sin(phi)
        rotation = np.array(
            [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float
        )
        world_normal = rotation @ body_normal
        world_normal /= np.linalg.norm(world_normal)
        nozzles.append(
            {
                "number": nozzle["number"],
                "throat_center_world_m": transform_point(center_mm, matrix),
                "throat_radius_m": float(nozzle["throat_radius_mm"]) * 1.0e-3,
                "throat_plane_normal_solid_to_fluid_world": world_normal.tolist(),
                "selection_rule": "abs(dot(x-center,normal))<=epsilon AND norm((x-center)-dot(x-center,normal)*normal)<=Rt+epsilon; do not use a 3-D sphere alone",
            }
        )
    metadata = {
        "status": "PASS",
        "units": "m",
        "configuration": args.configuration,
        "source_model": stem,
        "parameter_file_sha256": sha256(HERE / "geometry_parameters.yaml"),
        "preparer_sha256": sha256(Path(__file__).resolve()),
        "source_input_stl": str(source_stl),
        "source_input_units_and_frame": "mm in the NASA TSP body frame",
        "run": args.run,
        "run_conditions": parameters["runs"][str(args.run)] if args.run else None,
        "roll_deg": roll,
        "tsp_world_m": list(args.tsp_m),
        "body_to_world_matrix_mm_input_to_m_world": matrix.tolist(),
        "closed_ibm_stl": full,
        "patches": patches,
        "active_nozzle_throats": nozzles,
        "important": "Cerisse static IBM reads raw STL coordinates with identity transform; use this prepared file, not the millimetre master, in an SI-metre domain. For a custom high-resolution STL, pass it explicitly with --input-stl and verify closed_ibm_stl.source_sha256.",
    }
    metadata_path = output_dir / f"{name}_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"PASS: {full['destination']}")
    print(f"metadata: {metadata_path}")


if __name__ == "__main__":
    main()
