"""Shared loader for Cerisse time_probe.log files.

IMPORTANT — Cerisse probe sentinel semantics
--------------------------------------------
`src/set/Diagnosis.cpp::recordTimeProbe` writes `sum_over_level_cells / numPts`
once per step at level `cns.time_probe_lev`. When the probe's `in_box` is NOT
covered by any grid at that level (e.g. AMR has not refined there at this
step), the sum stays 0 while `numPts` is the geometric cell count, so the
output is **exactly 0.0** — this is a sentinel for "not sampled", NOT a real
pressure reading of 0 Pa. Treating it as a measurement biases mean/RMS and
causes low-f spectral leakage.

This loader:
  1. Skips repeated header rows written on restart (`comments="t"`).
  2. Returns a mask array that flags sentinel zeros as invalid.
  3. Emits a warning if any probe column has coverage < 50 %.

If you ever patch the Cerisse writer to emit NaN instead of 0, the mask here
still works (NaN is caught by `~np.isfinite`).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Cerisse writes 6-digit precision; the smallest non-zero pressure we ever see
# in the jet runs is ~50 kPa. Use 500 Pa = 0.5 kPa as the sentinel threshold.
SENTINEL_ABS_TOL = 500.0  # Pa


@dataclass
class ProbeLog:
    path: Path
    t: np.ndarray                 # shape (N,)
    values: np.ndarray            # shape (N, num_probes)
    names: list[str]              # len = num_probes
    mask: np.ndarray              # bool, shape (N, num_probes); True = valid
    coverage: np.ndarray          # per-probe fraction of valid rows

    def masked(self, probe_index: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (t, v) with sentinel zeros removed for one probe."""
        m = self.mask[:, probe_index]
        return self.t[m], self.values[m, probe_index]


def _parse_header(path: Path) -> list[str]:
    """Extract probe names from header line.

    Header example (commas inside the box indices confuse a naive split):
        time, pressure((0,64)(5,76)), pressure((0,128)(5,140)), ...
    So we split on ", " *only at depth 0* (not inside parens).
    """
    with path.open() as f:
        hdr = f.readline().strip()
    parts, depth, cur = [], 0, []
    for ch in hdr:
        if ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if cur:
        parts.append("".join(cur).strip())
    return parts[1:]  # drop leading "time"


def load_probe_log(
    path: str | Path,
    tol: float = SENTINEL_ABS_TOL,
    warn_below: float = 0.5,
) -> ProbeLog:
    """Load a Cerisse `time_probe.log` and mask sentinel zeros.

    Parameters
    ----------
    path: path to the log file.
    tol: |value| < tol treated as sentinel (not a measurement).
    warn_below: emit a warning for any probe whose valid-coverage falls below
        this fraction (default 50 %).
    """
    path = Path(path)
    names = _parse_header(path)
    # `comments="t"` skips lines starting with the literal "t" — the restart-
    # appended header rows ("time,...") — while numeric lines survive.
    data = np.loadtxt(path, delimiter=",", comments="t")
    t = data[:, 0]
    v = data[:, 1:]
    assert v.shape[1] == len(names), (
        f"Header columns ({len(names)}) != data columns ({v.shape[1]}) in {path}"
    )
    mask = np.isfinite(v) & (np.abs(v) > tol)
    coverage = mask.mean(axis=0)
    for name, cov in zip(names, coverage):
        if cov < warn_below:
            warnings.warn(
                f"{path.name}: probe '{name}' coverage = {cov:.0%} "
                f"(< {warn_below:.0%}); statistics may be unreliable.",
                stacklevel=2,
            )
    return ProbeLog(path=path, t=t, values=v, names=names,
                    mask=mask, coverage=coverage)


if __name__ == "__main__":  # quick self-check
    import sys, glob
    for f in sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else "npr*.log")):
        log = load_probe_log(f)
        print(f"{Path(f).name:12s} N={len(log.t):>6d}  coverage=",
              ", ".join(f"{c*100:3.0f}%" for c in log.coverage))
