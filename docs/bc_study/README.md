# Far-field BC reflection study

Quantitative 1-D benchmark of the Cerisse far-field boundary conditions, with
the **report in `bc_report.pdf`**.

## Files
- `bc_report.pdf`        — the report (principle + usage + benchmark + root cause)
- `bc_report.tex`        — its source
- `bc_1d_fast.py`        — vectorised 1-D Euler + HLLC core (generates all data)
- `bc_1d_reflection.py`  — readable (non-vectorised) reference implementation
- `run_study.py`         — driver: regenerates every figure + `bc_study_summary.json`
- `figures/`             — all figures
- `bc_study_summary.json`— numeric results

## Reproduce
    python3 run_study.py     # ~3 min, writes figures/ + bc_study_summary.json
    pdflatex bc_report.tex   # twice, for cross-refs

## Headline result
| BC | acoustic R (M=0) | mean-drift residual |
|----|------------------|---------------------|
| foextrap (code 2)        | 0.00 | 1.00 (no anchor) |
| Dirichlet (code 1 macro) | 0.00 | 0.50 |
| NSCBC-algebraic (code 7) | **0.78** | 0.00 (perfect anchor) |
| NSCBC-LODI σ=0.25 (n/a)  | 0.01 | 0.65 |

Root cause: the algebraic "NSCBC" (`bc_types.h:266`, `p_ghost = p_inf`) is the
σ→∞ (perfectly reflecting) limit of a proper Poinsot–Lele NSCBC. The ghost-cell
architecture cannot express the temporal relaxation that makes true NSCBC
non-reflecting. Recommendation: use `foextrap` (code 2) on lateral far-field
boundaries (the sphere cases already do).
