#!/usr/bin/env python3
"""Independent artifact checks for generated NASA Test 1853 geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import gmsh
import numpy as np
import trimesh

from build_test1853 import HERE, OUTPUT_ROOT, configurations, load_parameters


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_one_mesh(path: Path, process: bool = True) -> trimesh.Trimesh:
    mesh = trimesh.load_mesh(path, process=process)
    if not isinstance(mesh, trimesh.Trimesh):
        raise RuntimeError(f"{path}: expected one mesh, got {type(mesh)}")
    return mesh


def validate_stl(path: Path, expected_bounds: list[list[float]]) -> dict[str, Any]:
    raw = load_one_mesh(path, process=False)
    mesh = load_one_mesh(path, process=True)
    areas = raw.area_faces
    unique_mask = raw.unique_faces()
    checks = {
        "sha256": digest(path),
        "bytes": path.stat().st_size,
        "raw_triangles": int(len(raw.faces)),
        "processed_vertices": int(len(mesh.vertices)),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "is_volume": bool(mesh.is_volume),
        "body_count": int(mesh.body_count),
        "volume_mm3": float(mesh.volume),
        "bounds_mm": mesh.bounds.tolist(),
        "max_bounds_error_mm": float(
            np.max(np.abs(mesh.bounds - np.asarray(expected_bounds, dtype=float)))
        ),
        "degenerate_triangles_area_le_1e-12": int(np.count_nonzero(areas <= 1.0e-12)),
        "duplicate_triangles": int(len(unique_mask) - np.count_nonzero(unique_mask)),
    }
    required = (
        checks["watertight"]
        and checks["winding_consistent"]
        and checks["is_volume"]
        and checks["body_count"] == 1
        and checks["volume_mm3"] > 0.0
        and checks["max_bounds_error_mm"] <= 2.0e-3
        and checks["degenerate_triangles_area_le_1e-12"] == 0
        and checks["duplicate_triangles"] == 0
    )
    if not required:
        raise RuntimeError(f"STL validation failed for {path}: {checks}")
    return checks


def gts_self_intersection(path: Path) -> dict[str, Any]:
    if not shutil.which("stl2gts") or not shutil.which("gtscheck"):
        raise RuntimeError("stl2gts and gtscheck are required unless --skip-gts is used")
    mesh = load_one_mesh(path, process=False)
    ascii_stl = mesh.export(file_type="stl_ascii")
    converted = subprocess.run(
        ["stl2gts"], input=ascii_stl, text=True, capture_output=True, check=False
    )
    if converted.returncode != 0:
        raise RuntimeError(f"stl2gts failed for {path}: {converted.stderr}")
    checked = subprocess.run(
        ["gtscheck", "-v"],
        input=converted.stdout,
        text=True,
        capture_output=True,
        check=False,
    )
    result = {
        "return_code": checked.returncode,
        "passed": checked.returncode == 0,
        "diagnostic": (checked.stdout + checked.stderr).strip(),
        "meaning": "0=valid orientable manifold without self-intersection; 2=topology failure; 3=self-intersection",
    }
    if not result["passed"]:
        raise RuntimeError(f"GTS validation failed for {path}: {result}")
    return result


def validate_cad(path: Path, expected_bounds: list[list[float]]) -> dict[str, Any]:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
        gmsh.model.add("validation")
        gmsh.model.occ.importShapes(str(path), highestDimOnly=False)
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(3)
        if len(volumes) != 1:
            raise RuntimeError(f"{path}: expected one solid, got {volumes}")
        volume = volumes[0][1]
        bounds = np.asarray(gmsh.model.getBoundingBox(3, volume)).reshape((2, 3))
        difference = float(
            np.max(np.abs(bounds - np.asarray(expected_bounds, dtype=float)))
        )
        result = {
            "sha256": digest(path),
            "bytes": path.stat().st_size,
            "solid_count": 1,
            "surface_count": len(
                gmsh.model.getBoundary([(3, volume)], combined=False, recursive=False)
            ),
            "volume_mm3": float(gmsh.model.occ.getMass(3, volume)),
            "bounds_mm": bounds.tolist(),
            "max_bounds_error_mm": difference,
        }
        if result["volume_mm3"] <= 0.0 or difference > 2.0e-4:
            raise RuntimeError(f"CAD validation failed for {path}: {result}")
        return result
    finally:
        gmsh.finalize()


def validate_patches(
    p: dict[str, Any],
    root: Path,
    name: str,
    active: list[int],
    expected_areas_mm2: dict[int, float],
) -> dict[str, Any]:
    directory = root / "patches" / name
    wall = load_one_mesh(directory / "body_wall_excluding_throat_caps.stl", process=True)
    result: dict[str, Any] = {
        "body_wall": {
            "watertight_expected_false": not wall.is_watertight,
            "winding_consistent": bool(wall.is_winding_consistent),
            "triangles": int(len(wall.faces)),
        },
        "throat_caps": {},
    }
    if wall.is_watertight or not wall.is_winding_consistent:
        raise RuntimeError(f"{name}: open body-wall patch has wrong topology")
    found = sorted(
        int(path.stem.split("_")[1])
        for path in directory.glob("nozzle_*_throat_cap.stl")
    )
    if found != sorted(active):
        raise RuntimeError(f"{name}: expected throat caps {active}, found {found}")
    tolerance = float(p["mesh"]["mesh_area_relative_tolerance"])
    for number in found:
        path = directory / f"nozzle_{number}_throat_cap.stl"
        cap = load_one_mesh(path, process=False)
        expected = expected_areas_mm2[number]
        relative = abs(float(cap.area) - expected) / expected
        item = {
            "triangles": int(len(cap.faces)),
            "mesh_area_mm2": float(cap.area),
            "expected_area_mm2": expected,
            "relative_error": relative,
        }
        if relative > tolerance:
            raise RuntimeError(f"{name} Nozzle {number}: throat patch area failed: {item}")
        result["throat_caps"][str(number)] = item
    return result


def write_reports(root: Path, report: dict[str, Any]) -> None:
    directory = root / "reports"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "independent_validation.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Independent geometry validation",
        "",
        "| Model | STEP solid | BREP solid | STL watertight | STL volume | GTS self-intersection | Triangles |",
        "|---|---:|---:|---|---|---|---:|",
    ]
    for name, item in report["models"].items():
        gts = item.get("gts", {"passed": "SKIPPED"})
        lines.append(
            f"| `{name}` | {item['step']['solid_count']} | {item['brep']['solid_count']} | "
            f"{item['stl']['watertight']} | {item['stl']['is_volume']} | "
            f"{gts['passed']} | {item['stl']['raw_triangles']} |"
        )
    lines += [
        "",
        "All listed models passed CAD re-import, bounds, positive-volume, STL manifold/winding, duplicate/degenerate-face, patch-area, and (unless explicitly skipped) GTS self-intersection gates.",
        "",
        "SHA-256 hashes and numeric diagnostics are in `independent_validation.json`.",
    ]
    (directory / "independent_validation.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_ROOT)
    parser.add_argument(
        "--config",
        action="append",
        help="model stem to validate; repeat as needed (default: models in build audit)",
    )
    parser.add_argument("--skip-gts", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    p = load_parameters()
    cfgs = configurations(p)
    build_path = args.output / "reports" / "build_audit.json"
    build = json.loads(build_path.read_text(encoding="utf-8"))
    selected = args.config or list(build)
    unknown = sorted(set(selected) - set(cfgs))
    missing = sorted(set(selected) - set(build))
    if unknown or missing:
        raise RuntimeError(f"Unknown configurations={unknown}; missing audits={missing}")
    report: dict[str, Any] = {"status": "PASS", "models": {}}
    for name in selected:
        cfg = cfgs[name]
        expected_bounds = build[name]["occ"]["expected_bounds_mm"]
        stl_path = args.output / "stl" / f"{name}.stl"
        expected_areas: dict[int, float] = {}
        for nozzle in build[name]["nozzles"]:
            number = int(nozzle["number"])
            geometry_key = str(nozzle.get("geometry_key", number))
            expected_areas[number] = float(
                nozzle.get(
                    "expected_throat_area_mm2",
                    p["nozzles"][geometry_key]["throat_area_mm2"],
                )
            )
        item = {
            "step": validate_cad(args.output / "cad" / f"{name}.step", expected_bounds),
            "brep": validate_cad(args.output / "cad" / f"{name}.brep", expected_bounds),
            "stl": validate_stl(stl_path, expected_bounds),
            "patches": validate_patches(
                p,
                args.output,
                name,
                sorted(cfg["active"]),
                expected_areas,
            ),
        }
        if not args.skip_gts:
            item["gts"] = gts_self_intersection(stl_path)
        report["models"][name] = item
        print(f"PASS {name}")
    write_reports(args.output, report)
    print(f"PASS: independent validation report written under {args.output / 'reports'}")


if __name__ == "__main__":
    main()
