# Test-1853 H100 grid redesign (2026-08-08)

## Scope

These are case-level AMR configurations.  They do not modify GP geometry,
BI-CWLS, thermodynamic extension, Cartesian flux reconstruction, the full-cell
positivity limiter, or surface recovery.  Both inputs retain the pure shared-GP
full-Cartesian path and explicitly disable all existing crossing/cut-control
surrogates.  No EB, cut cell, open-face fraction, aggregate, or redistribution
path is used.

The launch script for any run using these inputs must still write and validate:

```text
method=pure_shared_gp_full_cartesian
cut_control=disabled
```

## Why the current mesh is slow

The repair run uses `256x176x176`, blocking factor 8, and maximum L0/L1 box size
64.  L0 regrids every four coarse steps and the L1 clock triggers an `lbase=1`
regrid every two coarse steps.  At chk03800 it has 36 L0 boxes, roughly 270-290
L1 boxes, and 56 L2 boxes.  Its approximately 15.8 million valid cells therefore
carry substantial ghost-fill, GP support, MPI metadata, kernel-launch, and
regrid overhead.

The historical H100 underexpanded-jet case used a fixed hierarchy, SFC,
`grid_eff=0.75`, and much larger boxes.  Its roughly 50.9 million valid cells
advanced faster despite doing far more cell work.  That comparison is not a
solver-only benchmark because the old case had no GP-IBM, used RK3, and had a
different AMR hierarchy, but it strongly identifies box fragmentation and
frequent dynamic regridding as avoidable costs in the present case.

## Profiles

### `inputs.chk4000_h100_large_boxes`

Use only with the existing `256x176x176` checkpoint family.  It changes L1/L2
clustering to `blocking_factor=16 16 8`, `max_grid_size=128 128 64`, SFC,
`grid_eff=0.75`, buffers `8 4`, and regrid intervals `16 16` (16 coarse steps
for `lbase=0`, about eight coarse steps for `lbase=1`).

The L0 BoxArray is checkpoint state and will remain fragmented after restart.
This profile can test coarser L1/L2 clustering and higher CFL without changing
the physical mesh, but it cannot establish the full large-box performance.

### `inputs.fresh_256x192x192_h100`

This is the preferred cold-start H100 grid qualification profile.  The domain
is `[-0.192,0.384] x [-0.216,0.216]^2` m, giving isotropic 2.25 mm L0 cells.
The 9,437,184-cell base grid is divisible by 32 in every direction and should
form eight nominal L0 boxes of about `128x96x96`, one per H100 rank.

L0 and L1 use blocking factor 32.  L2 deliberately remains at 8: its tagged
jet-core cap is only 12.4 mm ahead of the first solid point.  With the four-cell
L1 error buffer and worst block-alignment extension, the static bound leaves
about 4.00 mm, or 7.12 L2 cells, before the STL.  The runtime strict IBM support
audit remains authoritative.

The wider transverse domain is an architecture-oriented test domain, not a
completed domain-independence result.  Final SRP production boundaries still
require a separate domain-size study.

## Qualification sequence

1. Run the static audit in this directory.
2. Cold-start the fresh profile for grid construction only, with fallback-mask
   and prolongation diagnostics enabled.  Record per-level boxes/cells, HBM,
   GP target counts, surface-target counts, and coarse-fine support misses.
3. If the hierarchy is safe, run a short CFL matrix at 0.05, 0.10, 0.20, and
   0.30.  A low-CFL pass alone is not a flux/RK acceptance result.
4. Export the accepted dynamic hierarchy and switch to `amr.regrid_file` with a
   very large regrid interval, following the underexpanded-jet production case.
5. Reprofile with diagnostics disabled and a non-TPROF executable before using
   wall time as a production estimate.

No H100 job is submitted by the files in this directory.

## Current positivity limitation

These grid-qualification inputs use strict state checking, no soft positivity,
no clipping, no temperature-floor injection, and `stage_positivity=0`.  The
current 3-D conservative shared-face positivity limiter aborts when enabled,
while the historical `stage_positivity` option is not an independent SSPRK
Forward-Euler-bracket check.  This is an explicit integration-gate limitation;
these inputs are not a claim that the final production positivity requirement
has already been met.
