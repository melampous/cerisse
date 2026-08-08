#!/usr/bin/env python3
"""Polynomial/order checks for AFD same-side shifted WENO/TENO traces."""

from __future__ import annotations

import math


LINEAR = {
    1: (35.0 / 48.0, 7.0 / 24.0, -1.0 / 48.0),
    3: (5.0 / 16.0, 5.0 / 8.0, 1.0 / 16.0),
    5: (1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0),
    7: (-1.0 / 48.0, 7.0 / 24.0, 35.0 / 48.0),
}


def shifted_data(values: list[float], target_twice: int, left: bool):
    target = 0.5 * target_twice
    lower = target - 1.0 if left else target
    upper = lower + 1.0
    candidates: list[float] = []
    beta: list[float] = []
    for stencil in range(3):
        z = target - stencil
        candidates.append(
            0.5 * (z - 1.0) * (z - 2.0) * values[stencil]
            - z * (z - 2.0) * values[stencil + 1]
            + 0.5 * z * (z - 1.0) * values[stencil + 2]
        )
        second = (
            values[stencil]
            - 2.0 * values[stencil + 1]
            + values[stencil + 2]
        )
        a2 = 0.5 * second
        a1 = values[stencil + 1] - values[stencil] - a2
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
        sensor = 1.0 + (tau / (1.0e-40 + smoothness)) ** 2
        positive.append(split_positive * sensor)
        negative.append(split_negative * sensor)
        positive_linear.append(split_positive)
        negative_linear.append(split_negative)
    plus = sum(w * q for w, q in zip(positive, candidates)) / sum(positive)
    minus = sum(w * q for w, q in zip(negative, candidates)) / sum(negative)
    return sum(positive_linear) * plus - sum(negative_linear) * minus


def teno(values: list[float], target_twice: int, left: bool) -> float:
    candidates, beta, linear = shifted_data(values, target_twice, left)
    eps = math.ulp(1.0)
    tau = abs(
        abs(beta[0] - beta[2])
        - (beta[0] + 4.0 * beta[1] + beta[2]) / 6.0
    )
    gamma = [(1.0 + tau / (eps + item)) ** 6 for item in beta]
    selected = [item / sum(gamma) >= 1.0e-4 for item in gamma]
    if all(selected):
        return sum(w * q for w, q in zip(linear, candidates))
    weights = [abs(w) if keep else 0.0 for w, keep in zip(linear, selected)]
    if sum(weights) == 0.0:
        return wenoz(values, target_twice, left)
    return sum(w * q for w, q in zip(weights, candidates)) / sum(weights)


def polynomial(coefficients: list[float], x: float) -> float:
    return sum(value * x**power for power, value in enumerate(coefficients))


def check_polynomial_reproduction() -> None:
    coefficients = [0.7, -1.1, 0.4, 0.9, -0.3]
    for target_twice, linear in LINEAR.items():
        target = 0.5 * target_twice
        values = [polynomial(coefficients, float(i)) for i in range(5)]
        for left in (False, True):
            candidates, _, _ = shifted_data(values, target_twice, left)
            reconstructed = sum(w * q for w, q in zip(linear, candidates))
            exact = polynomial(coefficients, target)
            if not math.isclose(reconstructed, exact, rel_tol=2.0e-13,
                                abs_tol=2.0e-13):
                raise AssertionError(
                    f"quartic reproduction failed at target={target}, "
                    f"left={left}: {reconstructed} != {exact}"
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


def check_smooth_order() -> None:
    # The x^5 term gives a nonzero fifth-order interpolation error without
    # reaching roundoff on the tested spacings.
    coefficients = [1.0, 0.2, -0.3, 0.4, -0.25, 0.7]
    spacings = [0.2, 0.1, 0.05, 0.025]
    for name, reconstruction in (("WENO-Z5", wenoz), ("TENO5", teno)):
        for target_twice in LINEAR:
            target = 0.5 * target_twice
            for left in (False, True):
                errors = []
                for spacing in spacings:
                    values = [
                        polynomial(coefficients, (i - target) * spacing)
                        for i in range(5)
                    ]
                    errors.append(abs(reconstruction(
                        values, target_twice, left) - coefficients[0]))
                orders = [
                    math.log(errors[i] / errors[i + 1], 2.0)
                    for i in range(len(errors) - 1)
                ]
                aggregate_order = math.log(
                    errors[0] / errors[-1], 2.0
                ) / (len(errors) - 1)
                if aggregate_order < 4.5:
                    raise AssertionError(
                        f"{name} order failed at target={target}, left={left}: "
                        f"errors={errors}, orders={orders}, "
                        f"aggregate={aggregate_order}"
                    )


def main() -> None:
    check_polynomial_reproduction()
    check_standard_equivalence()
    check_smooth_order()
    print("AFD shifted same-side WENO/TENO reconstruction: PASS")


if __name__ == "__main__":
    main()
