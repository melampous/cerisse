# Euler-admissible conservative prolongation regression

This CPU unit regression embeds the strongest inadmissible parent block from
the deterministic Test-1853 `chk03059` regrid audit. It verifies that:

1. AMReX `lincc_interp` reproduces the logged negative internal-energy-density
   child;
2. the Cerisse parent-block mapper makes every child Euler-admissible and, in
   a separate branch, satisfies a configured specific-internal-energy floor;
3. one common child-deviation coefficient preserves the parent conservative
   average in Cartesian coordinates and the annular-volume-weighted parent
   average in R-Z; and
4. the mapper is bitwise transparent for uniform and nonuniform smooth
   admissible states, with and without the specific-energy constraint.

It exercises only full background cells. It contains no IBM geometry, cut
cell, area fraction, aggregate, or flow update.

```bash
make -j8
./main3d.gnu.ex
make -j8 DIM=2
./main2d.gnu.ex
```
