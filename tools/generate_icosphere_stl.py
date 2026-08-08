#!/usr/bin/env python3
"""Generate one fixed, outward-oriented binary icosphere STL."""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path


def normalize(point: tuple[float, float, float]) -> tuple[float, float, float]:
    magnitude = math.sqrt(sum(value * value for value in point))
    return tuple(value / magnitude for value in point)


def base_icosahedron() -> tuple[
    list[tuple[float, float, float]], list[tuple[int, int, int]]
]:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = [
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ]
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]
    return [normalize(vertex) for vertex in vertices], faces


def subdivide(
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, int, int]],
) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    vertices = list(vertices)
    cache: dict[tuple[int, int], int] = {}

    def midpoint(first: int, second: int) -> int:
        key = (min(first, second), max(first, second))
        if key not in cache:
            a, b = vertices[first], vertices[second]
            vertices.append(normalize(tuple((a[d] + b[d]) * 0.5 for d in range(3))))
            cache[key] = len(vertices) - 1
        return cache[key]

    refined = []
    for first, second, third in faces:
        ab = midpoint(first, second)
        bc = midpoint(second, third)
        ca = midpoint(third, first)
        refined.extend(((first, ab, ca), (second, bc, ab),
                        (third, ca, bc), (ab, bc, ca)))
    return vertices, refined


def cross(
    first: tuple[float, float, float], second: tuple[float, float, float]
) -> tuple[float, float, float]:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def write_binary_stl(
    output: Path,
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, int, int]],
    radius: float,
    center: tuple[float, float, float],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    header = b"Cerisse fixed validation icosphere".ljust(80, b"\0")
    with output.open("wb") as stream:
        stream.write(header)
        stream.write(struct.pack("<I", len(faces)))
        for face in faces:
            unit_points = [vertices[index] for index in face]
            edge1 = tuple(unit_points[1][d] - unit_points[0][d] for d in range(3))
            edge2 = tuple(unit_points[2][d] - unit_points[0][d] for d in range(3))
            normal = normalize(cross(edge1, edge2))
            centroid = tuple(sum(point[d] for point in unit_points) / 3.0
                             for d in range(3))
            if sum(normal[d] * centroid[d] for d in range(3)) < 0.0:
                unit_points[1], unit_points[2] = unit_points[2], unit_points[1]
                normal = tuple(-value for value in normal)
            points = [tuple(center[d] + radius * point[d] for d in range(3))
                      for point in unit_points]
            stream.write(struct.pack("<12fH", *normal, *points[0], *points[1],
                                     *points[2], 0))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--level", type=int, default=6,
                        help="Subdivision level; facets = 20*4^level (default: 6).")
    parser.add_argument("--radius", type=float, default=0.1)
    parser.add_argument("--center", type=float, nargs=3, default=(0.0, 0.0, 0.0))
    args = parser.parse_args()
    if not 0 <= args.level <= 8:
        parser.error("--level must be between 0 and 8")
    if not math.isfinite(args.radius) or args.radius <= 0.0:
        parser.error("--radius must be finite and positive")
    if not all(math.isfinite(value) for value in args.center):
        parser.error("--center values must be finite")

    vertices, faces = base_icosahedron()
    for _ in range(args.level):
        vertices, faces = subdivide(vertices, faces)
    write_binary_stl(args.output, vertices, faces, args.radius, tuple(args.center))
    print(f"{args.output}: {len(vertices)} vertices, {len(faces)} facets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
