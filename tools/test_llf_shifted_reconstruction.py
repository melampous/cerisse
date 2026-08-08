#!/usr/bin/env python3
"""Polynomial and order checks for shifted LLF-WENO/TENO flux traces."""

from __future__ import annotations

import math


LINEAR = {
    1: (3.0 / 5.0, 47.0 / 110.0, -3.0 / 110.0),
    3: (3.0 / 10.0, 3.0 / 5.0, 1.0 / 10.0),
    5: (1.0 / 10.0, 3.0 / 5.0, 3.0 / 10.0),
    7: (-3.0 / 110.0, 47.0 / 110.0, 3.0 / 5.0),
}


def pair_within_jump_limit(a: float, b: float, threshold: float) -> bool:
    if not (math.isfinite(a + b) and a > 0.0 and b > 0.0):
        return False
    return abs(b - a) / (abs(a) + abs(b)) <= threshold


def smooth_transition(metric: float, upper: float) -> float:
    if upper == 0.0:
        return 1.0 if metric == 0.0 else 0.0
    lower = 0.5 * upper
    if metric <= lower:
        return 1.0
    if metric >= upper:
        return 0.0
    coordinate = (upper - metric) / (upper - lower)
    return coordinate * coordinate * (3.0 - 2.0 * coordinate)


def polynomial(coefficients: list[float], x: float) -> float:
    return sum(value * x**power for power, value in enumerate(coefficients))


def cell_average(coefficients: list[float], center: float) -> float:
    result = 0.0
    for power, value in enumerate(coefficients):
        result += value * (
            (center + 0.5) ** (power + 1)
            - (center - 0.5) ** (power + 1)
        ) / (power + 1)
    return result


def shifted_data(values: list[float], target_twice: int, left: bool):
    target = 0.5 * target_twice
    lower = target - 1.0 if left else target
    upper = lower + 1.0
    candidates: list[float] = []
    beta: list[float] = []
    for stencil in range(3):
        z = target - stencil
        second = (
            values[stencil]
            - 2.0 * values[stencil + 1]
            + values[stencil + 2]
        )
        a2 = 0.5 * second
        a1 = values[stencil + 1] - values[stencil] - a2
        candidates.append(
            0.5 * (z - 1.0) * (z - 2.0) * values[stencil]
            - z * (z - 2.0) * values[stencil + 1]
            + 0.5 * z * (z - 1.0) * values[stencil + 2]
            - a2 / 12.0
        )
        lo = lower - stencil
        hi = upper - stencil
        first = (
            a1 * a1 * (hi - lo)
            + 2.0 * a1 * a2 * (hi * hi - lo * lo)
            + (4.0 / 3.0) * a2 * a2 * (hi**3 - lo**3)
        )
        beta.append(max(0.0, first + 4.0 * a2 * a2 * (hi - lo)))
    return candidates, beta, LINEAR[target_twice]


def wenoz(values: list[float], target_twice: int, left: bool) -> float:
    candidates, beta, linear = shifted_data(values, target_twice, left)
    tau = abs(beta[0] - beta[2])
    positive = []
    negative = []
    positive_linear = []
    negative_linear = []
    for weight, smoothness in zip(linear, beta):
        split_positive = 0.5 * (weight + 3.0 * abs(weight))
        split_negative = split_positive - weight
        sensor = 1.0 + tau / (1.0e-40 + smoothness)
        positive.append(split_positive * sensor)
        negative.append(split_negative * sensor)
        positive_linear.append(split_positive)
        negative_linear.append(split_negative)
    plus = sum(w * q for w, q in zip(positive, candidates)) / sum(positive)
    minus = sum(w * q for w, q in zip(negative, candidates)) / sum(negative)
    return sum(positive_linear) * plus - sum(negative_linear) * minus


def teno(values: list[float], target_twice: int, left: bool) -> float:
    candidates, beta, linear = shifted_data(values, target_twice, left)
    tau = abs(beta[0] - beta[2])
    gamma = [
        (1.0 + tau / (math.ulp(1.0) + smoothness)) ** 6
        for smoothness in beta
    ]
    selected = [item / sum(gamma) >= 1.0e-3 for item in gamma]
    if all(selected):
        return sum(w * q for w, q in zip(linear, candidates))
    weights = [abs(w) if keep else 0.0 for w, keep in zip(linear, selected)]
    if sum(weights) == 0.0:
        return wenoz(values, target_twice, left)
    return sum(w * q for w, q in zip(weights, candidates)) / sum(weights)


def check_polynomial_reproduction() -> None:
    coefficients = [0.7, -1.1, 0.4, 0.9, -0.3]
    values = [cell_average(coefficients, float(i)) for i in range(5)]
    for target_twice, linear in LINEAR.items():
        target = 0.5 * target_twice
        exact = polynomial(coefficients, target)
        for left in (False, True):
            candidates, _, _ = shifted_data(values, target_twice, left)
            reconstructed = sum(w * q for w, q in zip(linear, candidates))
            if not math.isclose(
                reconstructed, exact, rel_tol=3.0e-13, abs_tol=3.0e-13
            ):
                raise AssertionError(
                    f"quartic cell-average reproduction failed at "
                    f"target={target}, left={left}: {reconstructed} != {exact}"
                )


def check_standard_equivalence() -> None:
    values = [0.3, -0.8, 1.7, 0.2, -0.4]
    _, beta, _ = shifted_data(values, 5, True)
    d0 = values[0] - 2.0 * values[1] + values[2]
    d1 = values[1] - 2.0 * values[2] + values[3]
    d2 = values[2] - 2.0 * values[3] + values[4]
    expected = [
        13.0 / 12.0 * d0 * d0
        + 0.25 * (values[0] - 4.0 * values[1] + 3.0 * values[2]) ** 2,
        13.0 / 12.0 * d1 * d1 + 0.25 * (values[1] - values[3]) ** 2,
        13.0 / 12.0 * d2 * d2
        + 0.25 * (3.0 * values[2] - 4.0 * values[3] + values[4]) ** 2,
    ]
    for actual, reference in zip(beta, expected):
        if not math.isclose(actual, reference, rel_tol=2.0e-14,
                            abs_tol=2.0e-14):
            raise AssertionError(f"standard beta mismatch: {actual} != {reference}")


def check_masked_candidate_moments() -> None:
    """Reduced WENO5 candidate pairs must reproduce cubic FV data."""
    coefficients = [0.7, -1.1, 0.4, 0.9]
    values = [
        cell_average(coefficients, center)
        for center in (-2.0, -1.0, 0.0, 1.0, 2.0)
    ]
    q0 = (2.0 * values[2] + 5.0 * values[1] - values[0]) / 6.0
    q1 = (-values[3] + 5.0 * values[2] + 2.0 * values[1]) / 6.0
    q2 = (11.0 * values[2] - 7.0 * values[3] + 2.0 * values[4]) / 6.0
    exact = polynomial(coefficients, -0.5)
    reduced = ((q0 + q1) / 2.0, (3.0 * q1 + q2) / 4.0)
    for reconstructed in reduced:
        if not math.isclose(
            reconstructed, exact, rel_tol=3.0e-14, abs_tol=3.0e-14
        ):
            raise AssertionError(
                "masked WENO5 cubic moment failed: "
                f"{reconstructed} != {exact}"
            )


def check_smooth_order() -> None:
    coefficients = [1.0, 0.2, -0.3, 0.4, -0.25, 0.7]
    spacings = [0.2, 0.1, 0.05, 0.025]
    for name, reconstruction in (("WENO-Z5", wenoz), ("TENO5", teno)):
        for target_twice in LINEAR:
            target = 0.5 * target_twice
            for left in (False, True):
                errors = []
                for spacing in spacings:
                    shifted_coefficients = [
                        coefficient * spacing**power
                        for power, coefficient in enumerate(coefficients)
                    ]
                    values = [
                        cell_average(shifted_coefficients, i - target)
                        for i in range(5)
                    ]
                    errors.append(abs(
                        reconstruction(values, target_twice, left)
                        - coefficients[0]
                    ))
                aggregate = math.log(errors[0] / errors[-1], 2.0) / 3.0
                if aggregate < 4.5:
                    raise AssertionError(
                        f"{name} order failed at target={target}, left={left}: "
                        f"errors={errors}, aggregate={aggregate}"
                    )


def check_reduced_shell_selector() -> None:
    threshold = 0.15
    smooth = [1.0, 1.02, 1.04, 1.06]
    if not all(
        pair_within_jump_limit(smooth[i], smooth[i + 1], threshold)
        for i in range(len(smooth) - 1)
    ):
        raise AssertionError("smooth one-sided shell support was rejected")

    # A linearly smeared shock has zero interior second difference. The
    # adjacent-jump condition is therefore needed in addition to curvature.
    shock = [1.0, 2.0, 3.0, 4.0]
    if all(
        pair_within_jump_limit(shock[i], shock[i + 1], threshold)
        for i in range(len(shock) - 1)
    ):
        raise AssertionError("linearly smeared shock escaped shell fallback")

    weights = [smooth_transition(value, 0.15)
               for value in (0.0, 0.075, 0.10, 0.149, 0.15, 0.30)]
    if any(not (0.0 <= weight <= 1.0) for weight in weights):
        raise AssertionError(f"non-convex shell weights: {weights}")
    if any(weights[i] < weights[i + 1] for i in range(len(weights) - 1)):
        raise AssertionError(f"non-monotone shell weights: {weights}")
    if weights[0] != 1.0 or weights[-1] != 0.0:
        raise AssertionError(f"incorrect shell endpoints: {weights}")


def check_wall_normal_ghost_constraint() -> None:
    """Check the equal-distance quadratic used on the crossing line."""
    # u_n(0)=0. The formula must reproduce every quadratic exactly.
    for spacing in (0.37, 0.11, 0.025):
        a1, a2 = 1.7, -0.9
        first = a1 * spacing + a2 * spacing**2
        second = 2.0 * a1 * spacing + 4.0 * a2 * spacing**2
        ghost = second - 3.0 * first
        exact = -a1 * spacing + a2 * spacing**2
        if not math.isclose(ghost, exact, rel_tol=2.0e-14,
                            abs_tol=2.0e-14):
            raise AssertionError(
                f"quadratic wall-normal extension failed: {ghost} != {exact}"
            )

    # A smooth cubic remainder must converge as O(h^3), which is the value
    # accuracy required for an O(h^2) crossing-face divergence.
    a1, a2, a3 = 1.7, -0.9, 0.6
    errors = []
    for spacing in (0.2, 0.1, 0.05, 0.025):
        def normal_velocity(s: float) -> float:
            return a1 * s + a2 * s**2 + a3 * s**3

        ghost = normal_velocity(2.0 * spacing) - 3.0 * normal_velocity(spacing)
        errors.append(abs(ghost - normal_velocity(-spacing)))
    orders = [math.log(errors[i] / errors[i + 1], 2.0)
              for i in range(len(errors) - 1)]
    if min(orders) < 2.99:
        raise AssertionError(
            f"wall-normal constrained extension lost cubic accuracy: {orders}"
        )

    # Curved-wall regression: q(s) must be sampled on the Cartesian line
    # containing the actual first hit and solid cell.  Sampling on a separate
    # wall-normal projection line introduces an O(h) tangential offset and
    # fails this test after division by the mesh spacing in the RHS.
    radius = 0.2
    angle = 0.8
    hit = (radius * math.cos(angle), radius * math.sin(angle))
    normal = (math.cos(angle), math.sin(angle))

    def swirl_velocity(point: tuple[float, float]) -> tuple[float, float]:
        radial = math.hypot(point[0], point[1])
        return (-50.0 * point[1] / radial, 50.0 * point[0] / radial)

    def crossing_normal_velocity(s: float) -> float:
        velocity = swirl_velocity((hit[0] + s, hit[1]))
        return velocity[0] * normal[0] + velocity[1] * normal[1]

    curved_errors = []
    for spacing in (0.02, 0.01, 0.005, 0.0025):
        ghost = (
            crossing_normal_velocity(2.0 * spacing)
            - 3.0 * crossing_normal_velocity(spacing)
        )
        curved_errors.append(
            abs(ghost - crossing_normal_velocity(-spacing))
        )
    curved_orders = [
        math.log(curved_errors[i] / curved_errors[i + 1], 2.0)
        for i in range(len(curved_errors) - 1)
    ]
    if min(curved_orders) < 2.8:
        raise AssertionError(
            "curved-wall Cartesian crossing extension lost cubic accuracy: "
            f"{curved_orders}"
        )


def main() -> None:
    check_polynomial_reproduction()
    check_standard_equivalence()
    check_masked_candidate_moments()
    check_smooth_order()
    check_reduced_shell_selector()
    check_wall_normal_ghost_constraint()
    print("LLF shifted same-side WENO/TENO reconstruction: PASS")


if __name__ == "__main__":
    main()
