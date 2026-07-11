# Underexpanded Jet — central hub

All artefacts related to the underexpanded jet study (data, post-processing,
figures, animations, documentation) live under this directory.

## Layout

```
exm/underexpanded_jet/
├── README.md                  ← this file
├── 2d/                        ← 2D axisymmetric case shells
│   └── npr<NPR>_lev<L>_<mode>[_extra]/
├── 3d/                        ← 3D Cartesian case shells
│   └── npr<NPR>_lev<L>_<mode>[_extra]/
├── _docs/                     ← human-readable documentation
│   ├── underexpanded_jet_inventory.pdf   ← full case inventory (cx3+AWS+local)
│   └── work_summary.md                   ← what has been done & what was found
├── _analysis/                 ← all post-processing scripts + their outputs
│   ├── case2_powell_tam/      ← Powell–Tam mechanism investigation (Case 2 = npr3 lev4 NS)
│   ├── visc_vs_inviscid/      ← viscous vs inviscid 2D-axi comparison
│   ├── grid_convergence_2d/   ← 2D-axi grid-convergence study (coarse/medium/fine)
│   ├── psd_2d_vs_3d/          ← 2D vs 3D PSD comparison (job 64)
│   └── npz_data/              ← cached probe / mean-std .npz arrays
├── _gifs/                     ← animations (rz / cart, sym / asym, multi-panel comparisons)
├── _figures/                  ← stand-alone figures not tied to a single analysis script
├── _postproc/                 ← (reserved) reusable post-processing helpers
├── _report/                   ← (reserved) report-ready figures / LaTeX fragments
└── _org/                      ← organisation tooling (manifest, mover, README/PDF generators)
    ├── manifest.py            ← single source of truth: the 31 cases on cx3+AWS+local
    ├── gen_move_scripts.py    ← regenerates the host-local mv scripts
    ├── gen_readmes.py         ← regenerates per-case READMEs everywhere
    └── gen_inventory_pdf.py   ← regenerates _docs/underexpanded_jet_inventory.pdf
```

## Naming convention

```
<dim>/npr<NPR>_lev<L>_<mode>[_extra]/
```

| field | values | meaning |
|-------|--------|---------|
| `<dim>` | `2d` / `3d` | physical dimensionality |
| `<NPR>` | `2`, `2p5`, `3`, `4`, `5`, `7p5`, `10`, `15`, `20` | nozzle pressure ratio (decimal point → `p`) |
| `<L>` | `2`, `3`, `4` | AMR max level |
| `<mode>` | `euler` / `ns` | inviscid Euler vs viscous Navier–Stokes |
| `_extra` | `axisym`, `gridconv_<coarse|medium|fine>`, `short`, … | study qualifier |

`2d/*` cases are 2D axisymmetric (the `_axisym` suffix is kept explicit so it
is impossible to confuse with a hypothetical 2D Cartesian variant).

## Where to look first

1. **What data exists, where, and how big?**
   → `_docs/underexpanded_jet_inventory.pdf` (5 pages, landscape A4).
   The same information is encoded programmatically in `_org/manifest.py`.

2. **What has been investigated and what was found?**
   → `_docs/work_summary.md`.

3. **Per-case metadata** (NPR, level, mode, host, plt/chk index range, size,
   grid spacing, integration window): every case directory on every host has
   its own `README.md` auto-generated from the manifest.

## Hosts

| host | root |
|------|------|
| local (laptop) | `exm/underexpanded_jet/` (this tree — case shells only, no plt/chk) |
| cx3 EPHEMERAL  | `/rds/general/ephemeral/user/jq319/ephemeral/cerisse/exm/underexpanded_jet/` |
| AWS HPC        | `ssh aws-hpc:/shared/cerisse/exm/underexpanded_jet/` |

The local tree mirrors the directory structure of the heavy hosts so that
relative paths in scripts work uniformly; the actual plt/chk data lives on
cx3 and AWS only.

## Regenerating documentation

```bash
# refresh the case inventory PDF after editing _org/manifest.py
python exm/underexpanded_jet/_org/gen_inventory_pdf.py

# refresh per-case READMEs everywhere
python exm/underexpanded_jet/_org/gen_readmes.py

# regenerate host-local mover scripts (only needed if cases are added/renamed)
python exm/underexpanded_jet/_org/gen_move_scripts.py
```
