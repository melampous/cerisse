#!/usr/bin/env python3
"""Extract the Chapter 3 shock-cell features for the exit-profile study.

The locations are density-based and use the common Chapter 3 extraction:

* ``x_fc`` is the largest positive gradient of the raw centreline mean density
  between the first expansion minimum and the first downstream maximum;
* ``x_n`` are maxima of the centreline mean density after a 0.05 D_e moving
  box filter. Accepted maxima have
  prominence at least 0.02 and separation at least 0.5 D_e;
* every selected location is refined with a three-point quadratic fit.

The experimental locations are extracted from the raw sampled mean-density
profile, as specified in Section 3.4, with the same prominence and separation
criteria.  This script is the single source for the feature markers and tables
used by the exit-profile comparison.
"""

from pathlib import Path

import numpy as np
from scipy.signal import find_peaks


HERE = Path(__file__).resolve().parent
RHO_J = 1.6413

CASES = (
    ("baseline", "Baseline"),
    ("tophat", "Top hat"),
    ("walltanh", "Wall-attached tanh"),
    ("shift100", "Thin shifted tanh"),
)


def box_filter(values, x, width=0.05):
    """Moving box filter with the nominal width used in Section 3.4."""
    count = max(3, int(round(width / (x[1] - x[0]))))
    if count % 2 == 0:
        count += 1
    result = np.convolve(values, np.ones(count) / count, mode="same")
    result[: count // 2] = values[: count // 2]
    result[-(count // 2) :] = values[-(count // 2) :]
    return result


def refine(x, values, index):
    """Three-point quadratic refinement of a selected sample location."""
    if index <= 0 or index >= values.size - 1:
        return float(x[index])
    denominator = values[index - 1] - 2.0 * values[index] + values[index + 1]
    if denominator == 0.0:
        return float(x[index])
    offset = 0.5 * (values[index - 1] - values[index + 1]) / denominator
    return float(x[index] + offset * (x[1] - x[0]))


def closure_location(x, raw_density):
    expansion = np.where((x >= 0.35) & (x <= 1.0))[0]
    first_minimum = expansion[np.nanargmin(raw_density[expansion])]
    recompression = np.where(
        (np.arange(x.size) > first_minimum) & (x <= 1.7)
    )[0]
    first_maximum = recompression[np.nanargmax(raw_density[recompression])]
    gradient = np.gradient(raw_density, x)
    interval = np.arange(first_minimum, first_maximum + 1)
    selected = interval[np.nanargmax(gradient[interval])]
    return refine(x, gradient, selected)


def density_maxima(x, values, x_fc, maximum_count=5):
    minimum_distance = max(1, int(round(0.5 / (x[1] - x[0]))))
    peaks, _ = find_peaks(
        np.where(x > x_fc, values, -np.inf),
        prominence=0.02,
        distance=minimum_distance,
    )
    return np.array([refine(x, values, i) for i in peaks[:maximum_count]])


def numerical_features(case):
    data = np.load(HERE / f"AXIS8_{case}.npz")
    x = data["x"]
    raw_density = data["DensityMEAN"] / RHO_J
    x_fc = closure_location(x, raw_density)
    filtered_density = box_filter(raw_density, x, width=0.05)
    return x_fc, density_maxima(x, filtered_density, x_fc)


def experimental_features():
    data = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
    x = data["x"]
    raw_density = data["rho"]
    x_fc = closure_location(x, raw_density)
    return x_fc, density_maxima(x, raw_density, x_fc)


features = {}
experiment_xfc, experiment_xn = experimental_features()
features["experiment_xfc"] = experiment_xfc
features["experiment_xmax"] = experiment_xn

for key, _ in CASES:
    x_fc, x_n = numerical_features(key)
    features[f"{key}_xfc"] = x_fc
    features[f"{key}_xmax"] = x_n

np.savez(HERE / "shockcell_marks.npz", **features)

experiment_spacing = np.diff(experiment_xn)
print("Density-based Chapter 3 feature extraction")
print(
    f"{'Profile/reference':22s} {'x_fc':>8s} {'x_1':>8s} "
    f"{'S_1':>8s} {'S_2':>8s} {'S_3':>8s} {'S_4':>8s} {'E_S':>8s}"
)
print("-" * 86)


def report(label, x_fc, x_n):
    spacing = np.diff(x_n)
    error = (
        np.mean(np.abs(spacing[:3] - experiment_spacing[:3]))
        if label != "Experiment"
        else np.nan
    )
    print(
        f"{label:22s} {x_fc:8.5f} {x_n[0]:8.5f} "
        + " ".join(f"{value:8.5f}" for value in spacing[:4])
        + (f" {error:8.5f}" if np.isfinite(error) else f" {'--':>8s}")
    )


report("Experiment", experiment_xfc, experiment_xn)
for key, label in CASES:
    report(label, features[f"{key}_xfc"], features[f"{key}_xmax"])
