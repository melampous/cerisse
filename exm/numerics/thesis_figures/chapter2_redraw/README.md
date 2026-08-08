# Chapter 2 figure audit and redraw specification

This directory is a staging area.  Nothing here is included by the thesis
unless it is copied explicitly into `02_numerical_methods/figures`.

## Audit scope

The audit covers all eight figures referenced by the current Chapter 2.
All PDF fonts are embedded.  Seven figures are fully vector.  The field
matrix is a hybrid PDF with vector labels, six 491 by 491 pixel field layers,
and two raster colour bar gradients.  Those layers are about 300 dpi in the
source PDF and remain adequate at the current printed size.

| Figure | Asset | Source | Native size | Main issue |
|---|---|---|---|---|
| 2.1 | `rz_metric_grid.pdf` | thesis `rz_metric_grid.tex` | 154.3 by 82.6 mm | Latin Modern vector output is sound.  The axial arrow crosses the radial index labels.  Flux arrowheads can be read as physical flow directions. |
| 2.2 | `inviscid_reconstruction_stencil.pdf` | thesis `inviscid_reconstruction_stencil.tex` | 114.1 by 56.0 mm | Vector output is sound.  The split branches should be tied explicitly to \(g^+\) and \(g^-\). |
| 2.3 | `amrex_hierarchy_schematic.pdf` | thesis `amrex_hierarchy_schematic.tex` | 151.6 by 114.1 mm | Vector output and final font size are suitable.  The single coarse to fine arrow should not imply that ghost filling occurs only once per step. |
| 2.4 | `mms_convergence.pdf` | `../plot_mms_convergence.py` | 179.1 by 148.6 mm | It is reduced to about 84 per cent in the thesis.  Its effective font is about 7 pt.  The legend says “\(\rho u_y\) or \(\rho u_r\)”.  Panel (c) puts global \(L_2\) evolution errors and a first ring semidiscrete \(L_\infty\) residual on one axis. |
| 2.5 | `shu_osher_three_scheme.pdf` | `../plot_three_scheme_results.py` | 172.9 by 158.6 mm | Type 3 DejaVu fonts and an effective final font of about 7 pt.  Scheme names omit en dashes.  The plot shows only \(-1.5\leq x\leq2.5\), while the caption does not state that it is a cropped interval. |
| 2.6 | `afd_hllc_wenoz5_low_temporal_error_convergence.pdf` | `../../isentropic_vortex/plot_low_temporal_error_convergence.py` | 238.5 by 75.4 mm | The three across layout is reduced to about 63 per cent.  Its effective font is about 5.3 pt.  The annotation box in panel (c) hides part of the curve.  It also uses Type 3 DejaVu fonts. |
| 2.7 | `vortex_accuracy_three_scheme.pdf` | `../plot_three_scheme_results.py` | 174.8 by 80.8 mm | Type 3 DejaVu fonts and an effective final font of about 7 pt.  Norm notation and scheme names differ from the chapter notation. |
| 2.8 | `vortex_fields_three_scheme.pdf` | `../plot_three_scheme_results.py` | 187.3 by 112.3 mm | Effective labels are about 6 pt.  Scheme names omit en dashes.  Lower row panel labels use a different corner from the upper row.  The caption does not state that the upper row is cropped to the vortex core or that vorticity is evaluated by centred differences. |

The current TikZ figures use embedded Latin Modern Type 1 fonts.  The current
Matplotlib figures use three different font configurations.  The MMS plot
uses embedded STIX TrueType fonts.  The other four plots use Type 3 DejaVu
fonts, with DejaVu Sans math glyphs mixed with serif text in the three method
figures.  This is the main source of visible typographic inconsistency.

## Fixed style for every Chapter 2 figure

The thesis text block is 150 mm wide.  Full width figures must therefore be
generated at 5.9055 inches.  A script must not create a 7 to 9 inch canvas
and rely on LaTeX to reduce it.

- Text is Latin Modern Roman.  Mathematics uses Computer Modern or Latin
  Modern math.  All PDF fonts are embedded.  Type 3 output is rejected.
- Body labels and titles are 8.2 to 8.5 pt at the final printed size.  Tick
  labels and legends are at least 7.2 pt.  Panel labels are 8.2 pt bold.
- Panel labels are placed just outside the upper left of line-plot axes, or
  are included at the start of a centred panel title.  Field panels place the
  label inside the upper left with a small white backing.  Diagram panel
  titles remain at the upper left of their frame.
- Axes use 0.7 pt spines.  Main curves use 1.1 pt lines.  Reference slopes use
  0.8 to 0.9 pt grey dotted lines.  Markers are 4.4 pt.
- Only light grey dotted major grids are used.  Minor grids are used only on
  log plots when they improve reading.
- The colour blind safe method mapping is fixed:
  LLF–WENO-Z5 blue `#0072B2`, LLF–TENO5 green `#009E73`, and
  AFD–HLLC–WENO-Z5 vermillion `#D55E00`.
- The conserved variable mapping is fixed:
  \(\rho\) blue circle, axial or \(x\) momentum vermillion square,
  transverse or radial momentum green triangle, and \(\rho E\) magenta
  diamond.  Cartesian and RZ legends state \(\rho u_y\) and \(\rho u_r\)
  separately.
- Method names use en dashes everywhere.  Coordinate compounds use the same
  form as the prose, including “axisymmetric RZ”, “axis-aligned”, and
  “coarse–fine”.
- Error plots use the chapter notation
  \(\lVert e_q\rVert_{1,h}\), \(\lVert e_q\rVert_{2,h}\), and
  \(\lVert e_q\rVert_{\infty,h}\).  The first ring operator error remains
  \(\lVert e_{\mathcal L,\rho u_r}\rVert_{\infty,i=0}\).
- Line plots remain vector PDF.  A field layer or colour bar gradient may be
  rasterised at no less than 300 dpi at final size.  Text, axes, interface
  lines, ticks, and colour bar labels remain vector.

## Redraw decisions

1. Figure 2.1 keeps its two panel TikZ layout.  The unused axial arrow is
   removed from the one dimensional radial indexing sketch.  Weighted flux
   labels use line stubs without arrowheads.
2. Figure 2.2 keeps the mirrored six point stencil.  The two headings identify
   \(g^+\) and \(g^-\).  No mathematical stencil content changes.
3. Figure 2.3 keeps the hierarchy and subcycling panels.  The stage time
   wording is made explicit.  No AMR operation is added.
4. Figure 2.4 becomes a two by two plot.  Cartesian Euler, Cartesian
   Navier–Stokes, and global RZ errors occupy separate panels.  The first ring
   semidiscrete residual has its own fourth panel.  This removes the mixed norm
   and mixed scale axis.  Cartesian and RZ momentum notation use separate
   legends.
5. Figure 2.5 keeps the three vertically aligned resolution panels.  It is
   generated at 150 mm width with the fixed method mapping.
6. Figure 2.6 becomes two panels on the first row and one full width panel on
   the second row.  The large result box is removed.  A short arrow marks only
   the nonmonotone \(N=320\) density point.  Numerical details remain in the
   table and prose.
7. Figure 2.7 keeps its two panels.  The density ordinate becomes
   \(\lVert e_\rho\rVert_{2,h}\).  Fill state continues to distinguish uniform
   and AMR results.  The chapter must define the plotted perturbation kinetic
   energy \(K^\prime\) and its exact value before this figure.  They are used
   in the current ordinate but are not defined in the current prose.
8. Figure 2.8 keeps a two by three matrix with shared row colour scales.  All
   panel letters use the upper left.  Repeated in-panel level labels are
   removed because the dashed interface and caption define the two levels.
   Vector labels remain at least
   7.2 pt.  Raster field layers are written at 400 dpi.

The Chapter 2 captions must also be updated when the staged assets are
adopted.  Figure 2.4 must describe four panels.  Figure 2.5 must state that
only \(-1.5\leq x\leq2.5\) is shown.  Figure 2.8 must state that the upper
row shows \(0.30\leq x/L,y/L\leq0.70\) and that its vorticity is evaluated
from centred velocity differences.

## Rebuild

Run these commands from this directory.

```bash
pdflatex -interaction=nonstopmode -halt-on-error rz_metric_grid.tex
pdflatex -interaction=nonstopmode -halt-on-error inviscid_reconstruction_stencil.tex
pdflatex -interaction=nonstopmode -halt-on-error amrex_hierarchy_schematic.tex
python3 plot_mms_convergence.py
python3 plot_three_scheme_results.py --output-dir output
python3 plot_low_temporal_error_convergence.py \
  ../../isentropic_vortex/results/y9000x_afd_hllc_wenoz5_spatial_20260727 \
  ../../isentropic_vortex/results/y9000x_afd_hllc_wenoz5_spatial_angle45_20260727 \
  --output-dir output
```

The three TikZ PDFs are written beside their sources.  The five numerical
PDFs and their PNG previews are written to `output`.
