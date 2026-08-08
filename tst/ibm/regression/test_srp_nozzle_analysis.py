#!/usr/bin/env python3
"""Small regression tests for the planar SRP nozzle post-processing."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
ANALYZER = ROOT / "IBM/cases/srp_planar_70deg_nozzle_2d/analyze_nozzle.py"
SPEC = importlib.util.spec_from_file_location("srp_nozzle_analysis", ANALYZER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def uniform_section(n_aperture: int) -> dict[str, object]:
    """Return a constant-state throat represented by full Cartesian cells."""
    dy = 1.0 / n_aperture
    n_y = 4 * n_aperture
    y = (-2.0 + 0.5 * dy) + np.arange(n_y) * dy
    active_y = np.zeros(n_y, dtype=bool)
    first = (n_y - n_aperture) // 2
    active_y[first : first + n_aperture] = True
    rho = 2.0
    velocity = 3.0
    pressure = 5.0
    shape = (1, n_y)
    active = active_y[None, :]
    return {
        "x": np.asarray([0.0]),
        "y": y,
        "rho": np.full(shape, rho),
        "u": np.full(shape, velocity),
        "p": np.full(shape, pressure),
        "mach": np.full(shape, velocity / math.sqrt(MODULE.GAMMA * pressure / rho)),
        "active": active,
        "dy": dy,
        "analysis_level": 0,
    }


class SectionMetricsTest(unittest.TestCase):
    def test_constant_full_cell_mass_flow_is_exact(self) -> None:
        expected = 2.0 * 3.0 * 1.0
        for n_aperture in (16, 32, 64):
            with self.subTest(n_aperture=n_aperture):
                data = uniform_section(n_aperture)
                metrics = MODULE.section_metrics(
                    data,
                    0.0,
                    section_name="synthetic throat",
                    min_active_cells=16,
                )
                self.assertEqual(metrics["n_active"], n_aperture)
                self.assertAlmostEqual(
                    metrics["mass_flow_per_depth"],
                    expected,
                    places=14,
                )

    def test_underresolved_section_fails_closed(self) -> None:
        data = uniform_section(2)
        with self.assertRaisesRegex(RuntimeError, "only 2 active cells"):
            MODULE.section_metrics(
                data,
                0.0,
                section_name="synthetic throat",
                min_active_cells=16,
            )


if __name__ == "__main__":
    unittest.main()
