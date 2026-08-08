#!/usr/bin/env python3
"""Roundoff-level oracle for the pure full-cell R-Z limiter algebra.

The paired radial-momentum operator is

    D(F, P) = D_metric(F - P) + D_cart(P),

where ``F`` is the complete radial momentum flux and ``P`` is its pressure
companion.  This test intentionally has no geometry fractions, cut cells, EB
data, or named-case constants.  It exercises only full background annuli.
"""

from __future__ import annotations

from dataclasses import dataclass
import sys
import unittest


LOW = 0
HIGH = 1


@dataclass(frozen=True)
class RadialCell:
    """One full R-Z background cell with irrelevant common pi factors removed."""

    index: int
    dr: float

    def __post_init__(self) -> None:
        if self.index < 0 or not self.dr > 0.0:
            raise ValueError("the radial mesh must begin at r=0 with dr>0")

    @property
    def radii(self) -> tuple[float, float]:
        return self.index * self.dr, (self.index + 1) * self.dr

    @property
    def volume(self) -> float:
        r_lo, r_hi = self.radii
        return r_hi * r_hi - r_lo * r_lo

    def area(self, side: int) -> float:
        return 2.0 * self.radii[side]


def blend(low: float, high: float, theta: float) -> float:
    """Use the limiter's low + theta * (high - low) convention."""

    return low + theta * (high - low)


def paired_momentum_operator(
    cell: RadialCell,
    flux: tuple[float, float],
    pressure: tuple[float, float],
) -> float:
    """Return D_metric(F-P)+D_cart(P) for one full annular cell."""

    r_lo, _ = cell.radii
    metric_advective_lo = (
        cell.area(LOW) * (flux[LOW] - pressure[LOW])
        if r_lo > 0.0
        else 0.0
    )
    metric_advective_hi = cell.area(HIGH) * (
        flux[HIGH] - pressure[HIGH]
    )
    return (
        (metric_advective_lo - metric_advective_hi) / cell.volume
        + (pressure[LOW] - pressure[HIGH]) / cell.dr
    )


def scalar_metric_operator(
    cell: RadialCell, flux: tuple[float, float]
) -> float:
    """Return the metric divergence for a non-radial-momentum component.

    At r=0 the WENO channel stores an undivided metric h-flux.  It therefore
    has the explicit 2/volume coefficient used by the production assembler;
    multiplying it by the zero physical face area would incorrectly erase it.
    """

    r_lo, _ = cell.radii
    metric_lo = (
        2.0 * flux[LOW]
        if r_lo == 0.0
        else cell.area(LOW) * flux[LOW]
    )
    metric_hi = cell.area(HIGH) * flux[HIGH]
    return (metric_lo - metric_hi) / cell.volume


def paired_face_correction(
    cell: RadialCell,
    side: int,
    dt: float,
    delta_flux: float,
    delta_pressure: float,
) -> float:
    """Antidiffusive state correction contributed by one radial face."""

    orientation = 1.0 if side == LOW else -1.0
    area_over_volume = cell.area(side) / cell.volume
    return orientation * dt * (
        area_over_volume * delta_flux
        + (1.0 / cell.dr - area_over_volume) * delta_pressure
    )


def scalar_face_correction(
    cell: RadialCell, side: int, dt: float, delta_flux: float
) -> float:
    """One-face metric correction, including the zero-area axis channel."""

    orientation = 1.0 if side == LOW else -1.0
    r_lo, _ = cell.radii
    if side == LOW and r_lo == 0.0:
        coefficient = 2.0 / cell.volume
    else:
        coefficient = cell.area(side) / cell.volume
    return orientation * dt * coefficient * delta_flux


class RzDualChannelAlgebraTest(unittest.TestCase):
    def assert_roundoff(self, actual: float, expected: float) -> None:
        scale = max(1.0, abs(actual), abs(expected))
        tolerance = 64.0 * sys.float_info.epsilon * scale
        self.assertLessEqual(
            abs(actual - expected),
            tolerance,
            msg=(
                f"actual={actual:.17e}, expected={expected:.17e}, "
                f"error={abs(actual - expected):.3e}, tol={tolerance:.3e}"
            ),
        )

    def test_theta_endpoints_select_exact_low_and_high_pair(self) -> None:
        low_flux = (3.25, -7.5)
        high_flux = (-11.0, 19.75)
        low_pressure = (2.5, 8.0)
        high_pressure = (13.0, -4.25)

        for theta, expected_flux, expected_pressure in (
            (0.0, low_flux, low_pressure),
            (1.0, high_flux, high_pressure),
        ):
            selected_flux = tuple(
                blend(low_flux[q], high_flux[q], theta) for q in (LOW, HIGH)
            )
            selected_pressure = tuple(
                blend(low_pressure[q], high_pressure[q], theta)
                for q in (LOW, HIGH)
            )
            for actual, expected in zip(selected_flux, expected_flux):
                self.assert_roundoff(actual, expected)
            for actual, expected in zip(selected_pressure, expected_pressure):
                self.assert_roundoff(actual, expected)

            for index in (0, 4):
                cell = RadialCell(index=index, dr=0.037)
                actual = paired_momentum_operator(
                    cell, selected_flux, selected_pressure
                )
                expected = paired_momentum_operator(
                    cell, expected_flux, expected_pressure
                )
                self.assert_roundoff(actual, expected)

    def test_same_theta_pair_operator_is_affine_at_axis_and_off_axis(self) -> None:
        low_flux = (6.125, -3.75)
        high_flux = (-14.5, 27.0)
        low_pressure = (9.0, 4.25)
        high_pressure = (-2.0, 18.75)

        for index in (0, 3):
            cell = RadialCell(index=index, dr=0.021)
            low_operator = paired_momentum_operator(
                cell, low_flux, low_pressure
            )
            high_operator = paired_momentum_operator(
                cell, high_flux, high_pressure
            )
            for theta in (0.0, 0.125, 0.37, 0.875, 1.0):
                selected_flux = tuple(
                    blend(low_flux[q], high_flux[q], theta)
                    for q in (LOW, HIGH)
                )
                selected_pressure = tuple(
                    blend(low_pressure[q], high_pressure[q], theta)
                    for q in (LOW, HIGH)
                )
                direct = paired_momentum_operator(
                    cell, selected_flux, selected_pressure
                )
                affine = blend(low_operator, high_operator, theta)
                self.assert_roundoff(direct, affine)

    def test_facewise_pair_correction_matches_direct_local_operator(self) -> None:
        low_flux = (17.0, -8.0)
        high_flux = (-5.5, 31.0)
        low_pressure = (4.0, 15.0)
        high_pressure = (23.0, -6.5)
        dt = 7.0e-4

        # Different face coefficients exercise a genuinely local limiter.
        for index, theta in ((0, (0.23, 0.81)), (5, (0.64, 0.17))):
            cell = RadialCell(index=index, dr=0.013)
            selected_flux = tuple(
                blend(low_flux[q], high_flux[q], theta[q])
                for q in (LOW, HIGH)
            )
            selected_pressure = tuple(
                blend(low_pressure[q], high_pressure[q], theta[q])
                for q in (LOW, HIGH)
            )
            direct_delta = dt * (
                paired_momentum_operator(
                    cell, selected_flux, selected_pressure
                )
                - paired_momentum_operator(cell, low_flux, low_pressure)
            )
            correction_sum = sum(
                theta[q]
                * paired_face_correction(
                    cell,
                    q,
                    dt,
                    high_flux[q] - low_flux[q],
                    high_pressure[q] - low_pressure[q],
                )
                for q in (LOW, HIGH)
            )
            self.assert_roundoff(direct_delta, correction_sum)

    def test_axis_and_nonaxis_face_coefficients(self) -> None:
        dt = 3.0e-3
        delta_flux = 29.0
        delta_pressure = -7.25
        axis_cell = RadialCell(index=0, dr=0.04)
        outer_cell = RadialCell(index=4, dr=0.04)

        # Radial momentum: the complete-flux channel has zero axis area, but
        # the ordinary pressure-gradient companion remains active.
        axis_pair = paired_face_correction(
            axis_cell, LOW, dt, delta_flux, delta_pressure
        )
        self.assert_roundoff(axis_pair, dt * delta_pressure / axis_cell.dr)

        # Other conserved components use the explicit undivided metric h-flux
        # at the axis and ordinary full-face area away from it.
        axis_scalar = scalar_face_correction(
            axis_cell, LOW, dt, delta_flux
        )
        self.assert_roundoff(
            axis_scalar, 2.0 * dt * delta_flux / (axis_cell.dr**2)
        )
        nonaxis_scalar = scalar_face_correction(
            outer_cell, LOW, dt, delta_flux
        )
        self.assert_roundoff(
            nonaxis_scalar,
            dt * outer_cell.area(LOW) / outer_cell.volume * delta_flux,
        )

        # Direct one-face replacements independently verify both formulae.
        for cell, side in (
            (axis_cell, LOW),
            (axis_cell, HIGH),
            (outer_cell, LOW),
            (outer_cell, HIGH),
        ):
            low_flux = (2.0, -3.0)
            low_pressure = (5.0, 11.0)
            high_flux = list(low_flux)
            high_pressure = list(low_pressure)
            high_flux[side] += delta_flux
            high_pressure[side] += delta_pressure
            direct_pair = dt * (
                paired_momentum_operator(
                    cell, tuple(high_flux), tuple(high_pressure)
                )
                - paired_momentum_operator(cell, low_flux, low_pressure)
            )
            self.assert_roundoff(
                direct_pair,
                paired_face_correction(
                    cell, side, dt, delta_flux, delta_pressure
                ),
            )

            low_scalar = (7.0, -13.0)
            high_scalar = list(low_scalar)
            high_scalar[side] += delta_flux
            direct_scalar = dt * (
                scalar_metric_operator(cell, tuple(high_scalar))
                - scalar_metric_operator(cell, low_scalar)
            )
            self.assert_roundoff(
                direct_scalar,
                scalar_face_correction(cell, side, dt, delta_flux),
            )

    def test_one_shared_theta_drives_both_face_channels(self) -> None:
        cell = RadialCell(index=2, dr=0.025)
        low_flux = (8.0, -12.0)
        high_flux = (21.0, 4.0)
        low_pressure = (3.0, 17.0)
        high_pressure = (14.0, -9.0)
        theta = 0.43

        correct = paired_momentum_operator(
            cell,
            tuple(blend(low_flux[q], high_flux[q], theta) for q in (LOW, HIGH)),
            tuple(
                blend(low_pressure[q], high_pressure[q], theta)
                for q in (LOW, HIGH)
            ),
        )
        affine = blend(
            paired_momentum_operator(cell, low_flux, low_pressure),
            paired_momentum_operator(cell, high_flux, high_pressure),
            theta,
        )
        self.assert_roundoff(correct, affine)

        # A separately chosen pressure theta is not the same selected paired
        # operator.  Keep this as a guard against updating only F in R-Z.
        wrong_pressure_theta = 0.71
        mismatched = paired_momentum_operator(
            cell,
            tuple(blend(low_flux[q], high_flux[q], theta) for q in (LOW, HIGH)),
            tuple(
                blend(
                    low_pressure[q], high_pressure[q], wrong_pressure_theta
                )
                for q in (LOW, HIGH)
            ),
        )
        scale = max(1.0, abs(mismatched), abs(affine))
        self.assertGreater(
            abs(mismatched - affine), 1.0e6 * sys.float_info.epsilon * scale
        )

    def test_one_nonaxis_face_is_shared_by_both_adjacent_cells(self) -> None:
        left = RadialCell(index=0, dr=0.031)
        right = RadialCell(index=1, dr=0.031)
        dt = 1.7e-3
        theta = 0.38
        low_shared_flux, high_shared_flux = -4.0, 22.0
        low_shared_pressure, high_shared_pressure = 7.5, -3.0
        delta_flux = high_shared_flux - low_shared_flux
        delta_pressure = high_shared_pressure - low_shared_pressure

        # The same nodal face value is the high face of the left cell and the
        # low face of the right cell.  Verify each local operator independently.
        for cell, side in ((left, HIGH), (right, LOW)):
            base_flux = [1.25, -6.75]
            base_pressure = [2.0, 9.0]
            base_flux[side] = low_shared_flux
            base_pressure[side] = low_shared_pressure
            selected_flux = list(base_flux)
            selected_pressure = list(base_pressure)
            selected_flux[side] = blend(
                low_shared_flux, high_shared_flux, theta
            )
            selected_pressure[side] = blend(
                low_shared_pressure, high_shared_pressure, theta
            )
            direct_delta = dt * (
                paired_momentum_operator(
                    cell, tuple(selected_flux), tuple(selected_pressure)
                )
                - paired_momentum_operator(
                    cell, tuple(base_flux), tuple(base_pressure)
                )
            )
            expected = theta * paired_face_correction(
                cell, side, dt, delta_flux, delta_pressure
            )
            self.assert_roundoff(direct_delta, expected)


if __name__ == "__main__":
    unittest.main()
