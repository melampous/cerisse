#!/usr/bin/env python3
"""Audit the Mj=1.19 probe record for a defensible screech signature.

The script deliberately separates detection of a narrow-band near-field
component from identification of a closed screech feedback loop.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import coherence, csd, periodogram, stft, welch


ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "exm/underexpanded_jet/2d/panda_m119_L4_nscbc_4090"
SOURCE = CASE / "X_steady.npz"
OUT = ROOT / "analysis/panda_benchmark/m119_screech_audit"

NAMES = (
    [f"lp{i:02d}" for i in range(1, 11)]
    + ["ol1", "ol2", "ol3"]
    + [f"cl{i:02d}" for i in range(1, 11)]
    + ["f1a", "f1b", "f1c", "f1d", "f1e"]
    + ["f3a", "f3b", "f3c", "f3d"]
)
LP_Z_OVER_D = np.array([0.10, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])
D = 0.0254
F_EXPERIMENT = 8400.0
F_CANDIDATE = 9833.423137768435


def fixed_band_prominence(f: np.ndarray, psd: np.ndarray, fc: float) -> float:
    """Peak-to-median ratio using one fixed physical comparison band."""
    candidate = (f >= fc - 350.0) & (f <= fc + 350.0)
    floor = (f >= 7000.0) & (f <= 12500.0) & (
        (f < fc - 800.0) | (f > fc + 800.0)
    )
    return float(10.0 * np.log10(np.max(psd[candidate]) / np.median(psd[floor])))


def legacy_selected_peak(f: np.ndarray, psd: np.ndarray) -> tuple[float, float]:
    """Reproduce the moving 41-bin floor used in the original assessment."""
    search = (f >= 800.0) & (f <= 30000.0)
    fb = f[search]
    pb = psd[search]
    logp = np.log10(pb + 1.0e-300)
    local_floor = np.array(
        [
            np.median(logp[max(0, k - 41) : min(logp.size, k + 42)])
            for k in range(logp.size)
        ]
    )
    excess = 10.0 * (logp - local_floor)
    peak = int(np.argmax(excess))
    return float(fb[peak]), float(excess[peak])


def tone_amplitude(t: np.ndarray, x: np.ndarray, fc: float) -> float:
    """Hann-windowed peak sinusoidal amplitude at a fixed frequency."""
    w = np.hanning(x.size)
    x = x - np.mean(x)
    return float(2.0 * np.abs(np.sum(w * x * np.exp(-2j * np.pi * fc * t))) / np.sum(w))


def load_record() -> tuple[np.ndarray, np.ndarray, float]:
    data = np.load(SOURCE)
    t = np.asarray(data["tu"], dtype=float)
    x = np.asarray(data["X"], dtype=float)
    fs = float(data["fs"])
    if x.shape[1] != len(NAMES) or np.any(np.diff(t) <= 0.0):
        raise ValueError("Unexpected probe layout or non-monotone time vector")
    return t, x - np.mean(x, axis=0, keepdims=True), fs


def probe_metrics(t: np.ndarray, x: np.ndarray, fs: float) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    half = t[0] + 0.5 * (t[-1] - t[0])
    for i, name in enumerate(NAMES):
        f, p = periodogram(
            x[:, i], fs=fs, window="hann", detrend="linear", scaling="density"
        )
        search = (f >= 800.0) & (f <= 30000.0)
        global_peak = int(np.flatnonzero(search)[np.argmax(p[search])])
        legacy_peak_hz, legacy_peak_db = legacy_selected_peak(f, p)
        fw, pw = welch(
            x[:, i],
            fs=fs,
            window="hann",
            nperseg=4096,
            noverlap=2048,
            detrend="linear",
            scaling="density",
            average="median",
        )
        wband = (fw >= 7000.0) & (fw <= 12500.0)
        welch_peak = int(np.flatnonzero(wband)[np.argmax(pw[wband])])
        amp = tone_amplitude(t, x[:, i], F_CANDIDATE)
        early = t <= half
        late = t > half
        rows.append(
            {
                "probe": name,
                "group": name.rstrip("0123456789"),
                "global_peak_hz": float(f[global_peak]),
                "legacy_selected_peak_hz": legacy_peak_hz,
                "legacy_selected_peak_db": legacy_peak_db,
                "candidate_periodogram_prominence_db": fixed_band_prominence(
                    f, p, F_CANDIDATE
                ),
                "welch_band_peak_hz": float(fw[welch_peak]),
                "candidate_welch_prominence_db": fixed_band_prominence(
                    fw, pw, F_CANDIDATE
                ),
                "candidate_amplitude_pa": amp,
                "candidate_rms_fraction": amp / np.sqrt(2.0) / float(np.std(x[:, i])),
                "early_amplitude_pa": tone_amplitude(t[early], x[early, i], F_CANDIDATE),
                "late_amplitude_pa": tone_amplitude(t[late], x[late, i], F_CANDIDATE),
            }
        )
    return rows


def phase_audit(x: np.ndarray, fs: float) -> dict[str, object]:
    phases = []
    coherences = []
    frequency = None
    for i in range(10):
        f, spectrum = csd(
            x[:, 0],
            x[:, i],
            fs=fs,
            window="hann",
            nperseg=2048,
            noverlap=1024,
            detrend="linear",
        )
        _, coh = coherence(
            x[:, 0],
            x[:, i],
            fs=fs,
            window="hann",
            nperseg=2048,
            noverlap=1024,
            detrend="linear",
        )
        k = int(np.argmin(np.abs(f - F_CANDIDATE)))
        frequency = float(f[k])
        phases.append(float(np.angle(spectrum[k])))
        coherences.append(float(coh[k]))

    phase = np.unwrap(np.asarray(phases))
    coh = np.asarray(coherences)
    fits = []
    for count in (5, 6, 7):
        z = LP_Z_OVER_D[:count] * D
        weights = coh[:count] ** 2
        design = np.column_stack((z, np.ones(count)))
        coef = np.linalg.lstsq(
            design * np.sqrt(weights[:, None]),
            phase[:count] * np.sqrt(weights),
            rcond=None,
        )[0]
        prediction = design @ coef
        residual = float(
            np.sqrt(np.average((phase[:count] - prediction) ** 2, weights=weights))
        )
        fits.append(
            {
                "last_probe": NAMES[count - 1],
                "phase_speed_m_per_s": float(2.0 * np.pi * frequency / coef[0]),
                "phase_rmse_deg": float(np.degrees(residual)),
            }
        )
    return {
        "frequency_hz": frequency,
        "phase_rad": phase.tolist(),
        "coherence": coh.tolist(),
        "fits": fits,
    }


def coherent_mode_ranking(x: np.ndarray, fs: float) -> list[dict[str, float]]:
    transforms = []
    frequency = None
    for i in range(x.shape[1]):
        frequency, _, z = stft(
            x[:, i],
            fs=fs,
            window="hann",
            nperseg=2048,
            noverlap=1024,
            detrend="linear",
            boundary=None,
            padded=False,
        )
        transforms.append(z)
    zall = np.stack(transforms, axis=0)
    rows = []
    for k, fk in enumerate(frequency):
        if not 1000.0 <= fk <= 30000.0:
            continue
        q = zall[:, k, :]
        spectrum = q @ q.conj().T / q.shape[1]
        scale = np.sqrt(np.maximum(np.real(np.diag(spectrum)), 1.0e-300))
        normalized = spectrum / (scale[:, None] * scale[None, :])
        eigenvalues = np.linalg.eigvalsh(normalized)
        rows.append(
            {
                "frequency_hz": float(fk),
                "leading_mode_fraction": float(eigenvalues[-1] / np.sum(eigenvalues)),
            }
        )
    return sorted(rows, key=lambda row: row["leading_mode_fraction"], reverse=True)


def sliding_candidate(t: np.ndarray, x: np.ndarray, fs: float, probe: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    window = 0.002
    stride = 0.00025
    starts = np.arange(t[0], t[-1] - window + 1.0e-12, stride)
    centers = []
    ridge = []
    prominence = []
    for start in starts:
        chosen = (t >= start) & (t < start + window)
        f, p = periodogram(
            x[chosen, probe],
            fs=fs,
            window="hann",
            detrend="linear",
            scaling="density",
        )
        band = (f >= 7000.0) & (f <= 12500.0)
        peak = int(np.flatnonzero(band)[np.argmax(p[band])])
        centers.append((start + 0.5 * window) * 1000.0)
        ridge.append(f[peak] / 1000.0)
        prominence.append(fixed_band_prominence(f, p, F_CANDIDATE))
    return np.asarray(centers), np.asarray(ridge), np.asarray(prominence)


def make_figure(
    t: np.ndarray,
    x: np.ndarray,
    fs: float,
    metrics: list[dict[str, float | str]],
    phase: dict[str, object],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.4))

    ax = axes[0, 0]
    for name, color in (("lp03", "#0072B2"), ("cl02", "#D55E00"), ("f1a", "#009E73"), ("f3a", "#666666")):
        i = NAMES.index(name)
        f, p = welch(
            x[:, i], fs=fs, window="hann", nperseg=4096, noverlap=2048,
            detrend="linear", scaling="density", average="median"
        )
        keep = (f >= 1000.0) & (f <= 30000.0)
        ax.plot(f[keep] / 1000.0, 10.0 * np.log10(p[keep] / np.var(x[:, i])), label=name, color=color)
    ax.axvline(F_EXPERIMENT / 1000.0, color="black", linestyle=":", label="Panda 8.4 kHz")
    ax.axvline(F_CANDIDATE / 1000.0, color="#CC79A7", linestyle="--", label="candidate 9.83 kHz")
    ax.set(xlabel="frequency [kHz]", ylabel="variance-normalized PSD [dB/Hz]", title="Median Welch spectra (same estimator)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25)

    ax = axes[0, 1]
    for name, color in (("cl02", "#D55E00"), ("lp03", "#0072B2"), ("f1a", "#009E73"), ("f3a", "#666666")):
        center, _, prominence = sliding_candidate(t, x, fs, NAMES.index(name))
        ax.plot(center, prominence, marker="o", ms=3, label=name, color=color)
    ax.axhline(6.0, color="black", linestyle=":", linewidth=1)
    ax.set(xlabel="window centre [ms]", ylabel="9.83 kHz prominence [dB]", title="Two-millisecond moving-window persistence")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    fraction = np.array([float(row["candidate_rms_fraction"]) for row in metrics])
    colors = ["#0072B2"] * 10 + ["#56B4E9"] * 3 + ["#D55E00"] * 10 + ["#009E73"] * 5 + ["#666666"] * 4
    ax.bar(np.arange(len(NAMES)), fraction, color=colors, width=0.82)
    ax.set_xticks(np.arange(len(NAMES)))
    ax.set_xticklabels(NAMES, rotation=90, fontsize=7)
    ax.set(xlabel="probe", ylabel="fitted-tone RMS / total RMS", title="Spatial support of the 9.83 kHz component")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1, 1]
    phase_deg = np.degrees(np.asarray(phase["phase_rad"], dtype=float))
    coh = np.asarray(phase["coherence"], dtype=float)
    ax.plot(LP_Z_OVER_D, phase_deg, "o-", color="#CC79A7", label="unwrapped phase")
    fit_z = LP_Z_OVER_D[:5] * D
    weights = coh[:5] ** 2
    design = np.column_stack((fit_z, np.ones(5)))
    coef = np.linalg.lstsq(
        design * np.sqrt(weights[:, None]), phase_deg[:5] * np.sqrt(weights), rcond=None
    )[0]
    ax.plot(LP_Z_OVER_D[:5], design @ coef, linestyle="--", color="black", label="fit through lp05 only")
    ax.set(xlabel="lip-line station x/D", ylabel="phase relative to lp01 [deg]", title="The apparent upstream phase branch is local")
    ax.grid(alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(LP_Z_OVER_D, coh, "s:", color="#0072B2", label="MSC")
    ax2.set_ylabel("magnitude-squared coherence", color="#0072B2")
    ax2.set_ylim(0.0, 1.05)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], fontsize=8, loc="best")

    fig.suptitle(
        "Mj=1.19 2D probe audit: a local 9.83 kHz component is not a closed screech identification",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "m119_screech_audit.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t, x, fs = load_record()
    metrics = probe_metrics(t, x, fs)
    phase = phase_audit(x, fs)
    ranking = coherent_mode_ranking(x, fs)

    with (OUT / "probe_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)

    candidate_rank = min(
        range(len(ranking)),
        key=lambda i: abs(ranking[i]["frequency_hz"] - F_CANDIDATE),
    )
    summary = {
        "source": str(SOURCE.relative_to(ROOT)),
        "record_start_ms": float(t[0] * 1000.0),
        "record_end_ms": float(t[-1] * 1000.0),
        "record_length_ms": float((t[-1] - t[0]) * 1000.0),
        "samples": int(t.size),
        "sample_rate_hz": fs,
        "full_record_bin_width_hz": fs / t.size,
        "experimental_frequency_hz": F_EXPERIMENT,
        "candidate_frequency_hz": F_CANDIDATE,
        "candidate_strouhal_D_over_Uj": F_CANDIDATE * D / 363.0,
        "experimental_strouhal_D_over_Uj": F_EXPERIMENT * D / 363.0,
        "original_six_near_field_probes": [
            row["probe"]
            for row in metrics
            if abs(float(row["legacy_selected_peak_hz"]) - F_CANDIDATE) <= 250.0
            and float(row["legacy_selected_peak_db"]) > 10.0
        ],
        "phase_audit": phase,
        "candidate_coherent_mode_rank": candidate_rank + 1,
        "candidate_coherent_mode": ranking[candidate_rank],
        "top_coherent_modes": ranking[:10],
    }
    (OUT / "summary_metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figure(t, x, fs, metrics, phase)

    print(f"wrote {OUT / 'probe_metrics.csv'}")
    print(f"wrote {OUT / 'summary_metrics.json'}")
    print(f"wrote {OUT / 'm119_screech_audit.png'}")


if __name__ == "__main__":
    main()
