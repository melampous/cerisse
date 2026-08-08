#!/usr/bin/env python3
"""Unit tests for NPR15 primary/secondary shock identity tracking."""

from __future__ import annotations

import unittest

import numpy as np

from npr15_shock_tracking import (
    track_primary_mach_disk,
    track_secondary_subsonic_compression,
)


class ShockIdentityTest(unittest.TestCase):
    def test_stronger_downstream_peak_does_not_replace_mach_disk(self) -> None:
        z = np.linspace(2.30, 2.55, 251)
        primary_index = int(np.argmin(np.abs(z - 2.40)))
        secondary_index = int(np.argmin(np.abs(z - 2.48)))
        gradient = np.zeros((4, z.size))
        for radial_index in range(4):
            gradient[radial_index] = (
                np.exp(-((z - z[primary_index]) / 0.003) ** 2)
                + 1.2 * np.exp(-((z - z[secondary_index]) / 0.003) ** 2)
            )

        density = np.ones_like(gradient)
        density[:, primary_index + 1:] += 1.0
        density[:, secondary_index + 1:] += 0.5
        mach = np.full_like(gradient, 4.0)
        mach[:, primary_index + 1:] = 0.6
        mach[:, secondary_index + 1:] = 0.2

        primary = track_primary_mach_disk(
            gradient, density, mach, z, radial_count=4
        )
        secondary = track_secondary_subsonic_compression(
            gradient, mach, z, primary, radial_count=4
        )

        np.testing.assert_allclose(primary["position"], z[primary_index])
        np.testing.assert_allclose(secondary["position"], z[secondary_index])
        self.assertTrue(
            np.all(
                secondary["gradient_peak"] > primary["gradient_peak"]
            )
        )


if __name__ == "__main__":
    unittest.main()
