"""Physical shock trackers for the NPR15 underexpanded-jet diagnostics."""

from __future__ import annotations

import numpy as np


def subcell_peak(
    signal: np.ndarray, coordinate: np.ndarray, index: int
) -> float:
    if index <= 0 or index + 1 >= signal.size:
        return float(coordinate[index])
    left, center, right = signal[index - 1:index + 2]
    denominator = left - 2.0 * center + right
    offset = 0.0
    if denominator != 0.0:
        offset = 0.5 * (left - right) / denominator
    offset = float(np.clip(offset, -0.5, 0.5))
    return float(
        coordinate[index]
        + offset * (coordinate[index + 1] - coordinate[index])
    )


def local_compression_peaks(
    gradient: np.ndarray,
    coordinate: np.ndarray,
    lower: float,
    upper: float,
) -> np.ndarray:
    window = (coordinate >= lower) & (coordinate <= upper)
    interior = np.flatnonzero(
        window[1:-1]
        & (gradient[1:-1] > 0.0)
        & (gradient[1:-1] >= gradient[:-2])
        & (gradient[1:-1] > gradient[2:])
    )
    return interior + 1


def side_average(
    field: np.ndarray, radial_index: int, axial_index: int,
    side: int, width: int = 2,
) -> float:
    if side < 0:
        section = slice(max(0, axial_index - width), axial_index)
    else:
        section = slice(
            axial_index + 1,
            min(field.shape[1], axial_index + width + 1),
        )
    return float(np.mean(field[radial_index, section]))


def track_primary_mach_disk(
    density_gradient: np.ndarray,
    density: np.ndarray,
    mach: np.ndarray,
    z_over_d: np.ndarray,
    radial_count: int,
    z_min: float = 2.25,
    z_max: float = 2.55,
    upstream_mach_min: float = 1.5,
    downstream_mach_max: float = 1.5,
) -> dict[str, np.ndarray]:
    """Track the first supersonic-to-subsonic compression in each radial row."""

    position = np.full(radial_count, np.nan)
    cell_index = np.full(radial_count, -1, dtype=np.int64)
    gradient_peak = np.full(radial_count, np.nan)
    mach_before = np.full(radial_count, np.nan)
    mach_after = np.full(radial_count, np.nan)
    density_before = np.full(radial_count, np.nan)
    density_after = np.full(radial_count, np.nan)

    for radial_index in range(radial_count):
        candidates = local_compression_peaks(
            density_gradient[radial_index], z_over_d, z_min, z_max
        )
        for axial_index in candidates:
            m_before = side_average(
                mach, radial_index, axial_index, -1
            )
            m_after = side_average(
                mach, radial_index, axial_index, 1
            )
            rho_before = side_average(
                density, radial_index, axial_index, -1
            )
            rho_after = side_average(
                density, radial_index, axial_index, 1
            )
            if (
                m_before >= upstream_mach_min
                and m_after <= downstream_mach_max
                and rho_after > rho_before
            ):
                position[radial_index] = subcell_peak(
                    density_gradient[radial_index],
                    z_over_d,
                    int(axial_index),
                )
                cell_index[radial_index] = int(axial_index)
                gradient_peak[radial_index] = float(
                    density_gradient[radial_index, axial_index]
                )
                mach_before[radial_index] = m_before
                mach_after[radial_index] = m_after
                density_before[radial_index] = rho_before
                density_after[radial_index] = rho_after
                break

    return {
        "position": position,
        "cell_index": cell_index,
        "gradient_peak": gradient_peak,
        "mach_before": mach_before,
        "mach_after": mach_after,
        "density_before": density_before,
        "density_after": density_after,
    }


def track_secondary_subsonic_compression(
    density_gradient: np.ndarray,
    mach: np.ndarray,
    z_over_d: np.ndarray,
    primary: dict[str, np.ndarray],
    radial_count: int,
    z_max: float = 2.65,
    minimum_separation: float = 0.025,
    mach_ceiling: float = 1.0,
) -> dict[str, np.ndarray]:
    """Track the strongest post-disk compression whose two sides are subsonic."""

    position = np.full(radial_count, np.nan)
    cell_index = np.full(radial_count, -1, dtype=np.int64)
    gradient_peak = np.full(radial_count, np.nan)
    mach_before = np.full(radial_count, np.nan)
    mach_after = np.full(radial_count, np.nan)

    for radial_index in range(radial_count):
        primary_position = primary["position"][radial_index]
        if not np.isfinite(primary_position):
            continue
        candidates = local_compression_peaks(
            density_gradient[radial_index],
            z_over_d,
            primary_position + minimum_separation,
            z_max,
        )
        accepted = []
        for axial_index in candidates:
            m_before = side_average(
                mach, radial_index, axial_index, -1
            )
            m_after = side_average(
                mach, radial_index, axial_index, 1
            )
            if m_before <= mach_ceiling and m_after <= mach_ceiling:
                accepted.append((int(axial_index), m_before, m_after))
        if not accepted:
            continue
        axial_index, m_before, m_after = max(
            accepted,
            key=lambda item: density_gradient[radial_index, item[0]],
        )
        position[radial_index] = subcell_peak(
            density_gradient[radial_index], z_over_d, axial_index
        )
        cell_index[radial_index] = axial_index
        gradient_peak[radial_index] = float(
            density_gradient[radial_index, axial_index]
        )
        mach_before[radial_index] = m_before
        mach_after[radial_index] = m_after

    return {
        "position": position,
        "cell_index": cell_index,
        "gradient_peak": gradient_peak,
        "mach_before": mach_before,
        "mach_after": mach_after,
    }
