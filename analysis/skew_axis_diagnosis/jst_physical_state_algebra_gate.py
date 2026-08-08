#!/usr/bin/env python3
"""Algebraic gates for the R-Z Skew/JST physical-state correction.

This script deliberately has no Cerisse or plotting dependency.  It mirrors
the radial face-storage convention used by ``Skew.h`` and
``assemble_rz_flux_divergence`` for a scalar conservative component.

The two gates are:

1. A constant physical state with a face-varying C2 coefficient.  Dissipating
   rU/r_f manufactures a non-zero R-Z RHS; dissipating U and assigning zero
   physical JST metric flux at the axis gives zero to roundoff.
2. Smooth axis-regular even and odd modes.  With the corrected physical-state
   C2/C4 stencils and a zero-area axis JST flux, the first-ring contribution
   decays with the expected local parity rates.  Its cylindrical-volume-
   weighted L2 contribution is at least third order, so the axis skip cannot
   lower the designed global order of Skew4+JST.
"""

from __future__ import annotations

import math


# Skew4 coefficients from src/rhs/Skew.h.  Cell offsets relative to radial
# face f are [-2, -1, 0, +1].
C2 = (-1.0 / 6.0, -5.0 / 6.0, 7.0 / 6.0, -1.0 / 6.0)
C4 = (-1.0, 3.0, -3.0, 1.0)
OFFSETS = (-2, -1, 0, 1)


def radius(cell: int, spacing: float) -> float:
    """Signed cell-centre radius for an axis face at index zero."""

    return (cell + 0.5) * spacing


def integrated_jst_flux(
    face: int,
    spacing: float,
    epsilon2: float,
    epsilon4: float,
    state,
    *,
    metric_weighted_state: bool,
    skip_axis: bool,
) -> float:
    """Return the metric-integrated scalar JST flux G_f = r_f D_f.

    At a non-axis face, Cerisse stores the per-unit-area dissipative flux D_f
    and the R-Z divergence supplies r_f.  At the axis, the array slot is
    already a metric auxiliary flux, so this function directly returns G_0.
    """

    if face == 0 and skip_axis:
        return 0.0

    r_face = face * spacing
    c2_sum = 0.0
    c4_sum = 0.0
    for coefficient2, coefficient4, offset in zip(C2, C4, OFFSETS):
        cell = face + offset
        value = state(radius(cell, spacing))
        if metric_weighted_state:
            # At r_f > 0 this is the old stored-state factor r_j/r_f.  The
            # axis slot stores metric flux directly, hence it uses r_j.
            scale = radius(cell, spacing) / r_face if face else radius(cell, spacing)
        else:
            scale = 1.0
        c2_sum += coefficient2 * scale * value
        c4_sum += coefficient4 * scale * value

    stored_or_metric_flux = -epsilon2 * c2_sum + epsilon4 * c4_sum
    return stored_or_metric_flux if face == 0 else r_face * stored_or_metric_flux


def radial_rhs(
    ring: int,
    spacing: float,
    metric_flux_lo: float,
    metric_flux_hi: float,
) -> float:
    """Scalar R-Z flux divergence with Cerisse's annular volume factor."""

    radial_volume_factor = (2 * ring + 1) * spacing * spacing
    return 2.0 * (metric_flux_lo - metric_flux_hi) / radial_volume_factor


def observed_rates(values: list[float]) -> list[float]:
    return [
        math.log(abs(values[index] / values[index + 1]), 2.0)
        for index in range(len(values) - 1)
    ]


def gate_constant_state() -> None:
    spacing = 0.01
    constant = 2.5
    epsilons = (1.0, 3.0, 2.0, 4.0)
    state = lambda _r: constant

    old_fluxes = [
        integrated_jst_flux(
            face,
            spacing,
            epsilons[face],
            0.0,
            state,
            metric_weighted_state=True,
            skip_axis=False,
        )
        for face in range(len(epsilons))
    ]
    new_fluxes = [
        integrated_jst_flux(
            face,
            spacing,
            epsilons[face],
            0.0,
            state,
            metric_weighted_state=False,
            skip_axis=True,
        )
        for face in range(len(epsilons))
    ]

    print("constant-state variable-epsilon C2 gate")
    print("ring  old_RHS          analytic_old_RHS  corrected_RHS")
    for ring in range(len(epsilons) - 1):
        old_rhs = radial_rhs(ring, spacing, old_fluxes[ring], old_fluxes[ring + 1])
        new_rhs = radial_rhs(ring, spacing, new_fluxes[ring], new_fluxes[ring + 1])
        analytic = (
            2.0
            * constant
            * (epsilons[ring + 1] - epsilons[ring])
            / ((2 * ring + 1) * spacing)
        )
        print(f"{ring:4d}  {old_rhs: .12e}  {analytic: .12e}  {new_rhs: .12e}")
        assert math.isclose(old_rhs, analytic, rel_tol=2.0e-14, abs_tol=2.0e-14)
        assert abs(new_rhs) < 2.0e-13


def gate_smooth_axis_parity() -> None:
    # Even scalar modes represent rho, rho*u_z, rho*E, ... .  The odd mode
    # represents radial momentum.  Defining them on signed r automatically
    # supplies the correct axis ghost parity.
    even_state = lambda r: 1.0 + 0.3 * r * r + 0.2 * r**4
    odd_state = lambda r: 0.7 * r - 0.2 * r**3 + 0.1 * r**5
    spacings = [1.0 / n for n in (16, 32, 64, 128, 256)]

    first_ring_rhs: dict[str, list[float]] = {"even": [], "odd": []}
    weighted_l2: dict[str, list[float]] = {"even": [], "odd": []}
    for name, state in (("even", even_state), ("odd", odd_state)):
        for spacing in spacings:
            # In a smooth field the JST pressure/density sensor is O(h^2).
            # epsilon2=h^2 captures that asymptotic scaling; epsilon4=1
            # exercises the active third-difference background damping.
            axis_flux = integrated_jst_flux(
                0,
                spacing,
                spacing * spacing,
                1.0,
                state,
                metric_weighted_state=False,
                skip_axis=True,
            )
            first_face_flux = integrated_jst_flux(
                1,
                spacing,
                spacing * spacing,
                1.0,
                state,
                metric_weighted_state=False,
                skip_axis=True,
            )
            rhs = radial_rhs(0, spacing, axis_flux, first_face_flux)
            first_ring_rhs[name].append(rhs)
            # The first annulus has radial measure proportional to h^2, so
            # its contribution to a cylindrical-volume L2 norm carries h.
            weighted_l2[name].append(abs(rhs) * spacing)

    print("\nsmooth parity gate for corrected C2+C4 axis treatment")
    print("mode  first-ring pointwise rates       cylindrical-L2 contribution rates")
    for name in ("even", "odd"):
        point_rates = observed_rates(first_ring_rhs[name])
        l2_rates = observed_rates(weighted_l2[name])
        print(
            f"{name:4s}  "
            + " ".join(f"{rate:7.4f}" for rate in point_rates)
            + "       "
            + " ".join(f"{rate:7.4f}" for rate in l2_rates)
        )

    # Skew4+JST is globally third order because of its active C4 term.  The
    # even variables are third order already pointwise in the first ring.  An
    # odd radial-momentum mode is second order pointwise there, but the first
    # annulus has O(h^2) cylindrical measure and hence contributes at third
    # order to the global L2 norm.  Most importantly, neither mode has an O(1)
    # or inverse-h axis defect.
    assert min(observed_rates(first_ring_rhs["even"])[-2:]) > 2.95
    assert min(observed_rates(first_ring_rhs["odd"])[-2:]) > 1.95
    assert min(observed_rates(weighted_l2["even"])[-2:]) > 3.95
    assert min(observed_rates(weighted_l2["odd"])[-2:]) > 2.95


def main() -> None:
    # Coefficient identities behind the constant-state proof.
    assert abs(sum(C2)) < 1.0e-15
    first_moment = sum(
        coefficient * (offset + 0.5)
        for coefficient, offset in zip(C2, OFFSETS)
    )
    assert math.isclose(first_moment, 1.0, rel_tol=0.0, abs_tol=1.0e-15)

    gate_constant_state()
    gate_smooth_axis_parity()
    print("\nPASS: physical-state R-Z JST algebra gates")


if __name__ == "__main__":
    main()
