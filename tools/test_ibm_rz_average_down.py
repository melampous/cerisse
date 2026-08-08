#!/usr/bin/env python3
"""Fast regression for IBM mixed-cell average-down metric semantics.

This is intentionally independent of a Cerisse executable: it checks the
closed-form RZ conservation identity and fails if the production source drops
the geometry-volume path.  A normal 2-D IBM build remains the compile gate.
"""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CNS = ROOT / "src" / "CNS.cpp"


def corrected_average(values, fluid, volumes, rz, coarse_before):
    active = [index for index, marker in enumerate(fluid) if marker]
    solid = [index for index, marker in enumerate(fluid) if not marker]
    if not active or not solid:
        return coarse_before, False
    weights = volumes if rz else [1.0] * len(values)
    denominator = sum(weights[index] for index in active)
    result = sum(values[index] * weights[index] for index in active) / denominator
    return result, True


def annular_volumes(radial_index: int) -> list[float]:
    # AMReX RZ uses coordinate 0 as radius.  A common axial-width factor
    # cancels.  Ordering is (i0,j0), (i1,j0), (i0,j1), (i1,j1), matching a
    # 2x2 refinement patch.
    inner = math.pi * ((radial_index + 1) ** 2 - radial_index**2)
    outer = math.pi * ((radial_index + 2) ** 2 - (radial_index + 1) ** 2)
    return [inner, outer, inner, outer]


def main() -> None:
    source = CNS.read_text(encoding="utf-8")
    required = (
        "const bool metric_weighted = fine_lev.geom.IsRZ();",
        "fine_lev.geom.GetVolume(*fine_volume",
        "metric_weighted\n                ? fvol(fi,fj,0,0) : Real(1.0)",
        "fluid_volume += weight;",
        "const Real inv = Real(1.0) / fluid_volume;",
        "if (n_fluid > 0 && n_solid > 0)",
    )
    missing = [token for token in required if token not in source]
    assert not missing, f"production RZ average-down path missing: {missing}"

    values = [1.0, 9.0, 3.0, 5.0]
    fluid = [True, False, True, True]
    prior = -17.0

    cartesian, valid = corrected_average(
        values, fluid, [99.0] * 4, False, prior
    )
    assert valid
    assert cartesian == (1.0 + 3.0 + 5.0) / 3.0

    for radial_index in (0, 4):
        volumes = annular_volumes(radial_index)
        rz_average, valid = corrected_average(values, fluid, volumes, True, prior)
        assert valid
        active = (0, 2, 3)
        fine_integral = sum(values[index] * volumes[index] for index in active)
        active_volume = sum(volumes[index] for index in active)
        assert math.isclose(
            rz_average * active_volume, fine_integral, rel_tol=2.0e-16
        )
        if radial_index == 0:
            # The axis-adjacent radial child volumes are exactly 1:3.  This
            # deliberately differs from the legacy arithmetic result (3).
            assert math.isclose(rz_average, 19.0 / 5.0, rel_tol=5.0e-16)
            assert not math.isclose(rz_average, cartesian)

    unchanged, valid = corrected_average(
        values, [True] * 4, annular_volumes(0), True, prior
    )
    assert not valid and unchanged == prior
    unchanged, valid = corrected_average(
        values, [False] * 4, annular_volumes(0), True, prior
    )
    assert not valid and unchanged == prior

    print("PASS: IBM mixed-cell average-down preserves Cartesian arithmetic")
    print("PASS: RZ mixed-cell average-down conserves active annular volume")
    print("PASS: pure-fluid/pure-solid overwrite semantics are unchanged")


if __name__ == "__main__":
    main()
