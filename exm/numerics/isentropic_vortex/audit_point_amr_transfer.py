#!/usr/bin/env python3
"""Audit AMR transfer operators for cell-centred point samples.

This script is independent of the flow solver.  It compares the transfer
semantics used by AMReX's built-in quartic interpolaters with a matched pair
of point-value operators for refinement ratio two.

The point prolongation is the degree-four Lagrange interpolant through five
coarse cell centres.  The point restriction is the degree-five Lagrange
interpolant through the six fine cell centres nearest a coarse cell centre.
All coefficients are represented exactly before conversion to float.
"""

from __future__ import annotations

from fractions import Fraction
import math

import numpy as np


POINT_PROLONG_LEFT = (
    Fraction(-45, 2048),
    Fraction(105, 512),
    Fraction(945, 1024),
    Fraction(-63, 512),
    Fraction(35, 2048),
)
POINT_PROLONG_RIGHT = tuple(reversed(POINT_PROLONG_LEFT))

# Values currently written in AMReX_Interp_C.H.  Their eight-decimal
# truncation does not preserve constants exactly.
AMREX_POINT_LEFT_DECIMAL = (
    0.01708984,
    -0.12304688,
    0.92285156,
    0.20507812,
    -0.02197266,
)

# CellConservativeQuartic left-child coefficients in one dimension.  These
# act on coarse cell averages, not coarse point samples.
CONSERVATIVE_QUARTIC_LEFT = (
    Fraction(-3, 128),
    Fraction(11, 64),
    Fraction(1, 1),
    Fraction(-11, 64),
    Fraction(3, 128),
)
CONSERVATIVE_QUARTIC_RIGHT = tuple(
    2 * Fraction(int(offset == 0), 1) - coefficient
    for offset, coefficient in zip(range(-2, 3), CONSERVATIVE_QUARTIC_LEFT)
)

POINT_RESTRICT_CUBIC = (
    Fraction(-1, 16),
    Fraction(9, 16),
    Fraction(9, 16),
    Fraction(-1, 16),
)
POINT_RESTRICT_QUINTIC = (
    Fraction(3, 256),
    Fraction(-25, 256),
    Fraction(75, 128),
    Fraction(75, 128),
    Fraction(-25, 256),
    Fraction(3, 256),
)


def weighted_value(weights: tuple[Fraction, ...], values: np.ndarray) -> float:
    return float(sum(float(weight) * value for weight, value in zip(weights, values)))


def polynomial_reproduction() -> None:
    print("Polynomial reproduction at refinement ratio two")
    print("degree  point-P error  FV-P-on-points error  avg-R error  cubic-R error  quintic-R error")

    coarse_x = np.arange(-2, 3, dtype=float)
    left_x = -0.25
    fine_average_x = np.array([-0.25, 0.25])
    fine_cubic_x = np.array([-0.75, -0.25, 0.25, 0.75])
    fine_quintic_x = np.array([-1.25, -0.75, -0.25, 0.25, 0.75, 1.25])

    for degree in range(7):
        exact_left = left_x**degree
        point_left = weighted_value(POINT_PROLONG_LEFT, coarse_x**degree)
        finite_volume_left = weighted_value(
            CONSERVATIVE_QUARTIC_LEFT, coarse_x**degree
        )
        exact_centre = 0.0**degree
        average_restriction = np.mean(fine_average_x**degree)
        cubic_restriction = weighted_value(
            POINT_RESTRICT_CUBIC, fine_cubic_x**degree
        )
        quintic_restriction = weighted_value(
            POINT_RESTRICT_QUINTIC, fine_quintic_x**degree
        )
        print(
            f"{degree:6d}  {abs(point_left-exact_left):13.6e}"
            f"  {abs(finite_volume_left-exact_left):20.6e}"
            f"  {abs(average_restriction-exact_centre):11.6e}"
            f"  {abs(cubic_restriction-exact_centre):13.6e}"
            f"  {abs(quintic_restriction-exact_centre):15.6e}"
        )


def constant_reproduction() -> None:
    exact_sum = float(sum(POINT_PROLONG_LEFT))
    decimal_sum = sum(AMREX_POINT_LEFT_DECIMAL)
    print("\nConstant reproduction")
    print(f"exact point-quartic one-dimensional sum = {exact_sum:.17g}")
    print(f"AMReX decimal one-dimensional sum       = {decimal_sum:.17g}")
    for dimension in (1, 2, 3):
        error = abs(decimal_sum**dimension - 1.0)
        print(f"AMReX tensor constant error in {dimension}D       = {error:.9e}")


def smooth_convergence() -> None:
    print("\nSmooth sine transfer errors")
    print("N   avg-R rate/error        cubic-R rate/error      quintic-R rate/error")
    previous = None
    for ncell in (16, 32, 64, 128, 256):
        coarse_h = 1.0 / ncell
        centre = 0.37
        exact = math.sin(2.0 * math.pi * centre)

        def sample(offset_quarters: int) -> float:
            return math.sin(
                2.0 * math.pi * (centre + offset_quarters * coarse_h / 4.0)
            )

        average = 0.5 * (sample(-1) + sample(1))
        cubic = weighted_value(
            POINT_RESTRICT_CUBIC,
            np.array([sample(offset) for offset in (-3, -1, 1, 3)]),
        )
        quintic = weighted_value(
            POINT_RESTRICT_QUINTIC,
            np.array([sample(offset) for offset in (-5, -3, -1, 1, 3, 5)]),
        )
        errors = tuple(abs(value - exact) for value in (average, cubic, quintic))
        if previous is None:
            rates = (math.nan,) * 3
        else:
            rates = tuple(math.log(old / new, 2.0) for old, new in zip(previous, errors))
        print(
            f"{ncell:3d} "
            + "  ".join(
                f"{rate:5.2f}/{error:.6e}" for rate, error in zip(rates, errors)
            )
        )
        previous = errors


def main() -> None:
    polynomial_reproduction()
    constant_reproduction()
    smooth_convergence()


if __name__ == "__main__":
    main()
