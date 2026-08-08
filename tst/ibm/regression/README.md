# IBM regression helpers

This directory contains test drivers and log checks. Nothing here is included
by a production executable.

`results/code_cleanup_20260723/` contains the code-layout cleanup regressions:

- matched lazy-cache and forced-audit Mach-4 circle outputs;
- frozen-protocol PM expansion and compression-corner outputs;
- axis-aligned and `(0.37h,0.23h)` shifted smooth curved-wall MMS outputs.

The valid comparisons and provenance are summarized in
`results/code_cleanup_20260723/README.md`. Generated plotfiles and logs are
test artifacts and are never included by `src/`.

MMS runners accept either relative or absolute artifact roots. They resolve
the executable and output paths before passing them to AMReX, so an absolute
path is not accidentally prefixed with the case working directory.
