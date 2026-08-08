# R-Z Skew/JST physical-state algebra gate

Run:

```bash
python3 analysis/skew_axis_diagnosis/jst_physical_state_algebra_gate.py
```

The script mirrors the scalar part of the radial face convention in
`src/rhs/Skew.h` and `src/tim/compute_rhs.cpp`; it does not run the CFD solver.

For the Skew4 C2 coefficients,

\[
\sum_l c_l=0,\qquad
\sum_l c_l\xi_l=1.
\]

Consequently, for a constant physical state \(U\), the old metric-state
dissipation gives

\[
D_f^{old}=-\epsilon_f\sum_lc_l\frac{r_l}{r_f}U
         =-\epsilon_f\frac{h}{r_f}U,
\]

and its metric-integrated flux is \(G_f^{old}=-\epsilon_f hU\).  If the JST
sensor makes \(\epsilon_f\) vary across faces, ring \(i\) receives the spurious
source

\[
R_i^{old}=\frac{2U(\epsilon_{i+1}-\epsilon_i)}{(2i+1)h}.
\]

The corrected physical-state stencil has \(D_f\propto\sum_lc_lU=0\), and the
physical JST contribution at the zero-area axis face is also zero.  Its RHS is
therefore zero to roundoff for arbitrary facewise \(\epsilon_f\).

The second gate uses smooth even and odd axis-regular polynomials.  With the
smooth-sensor scaling \(\epsilon_2=O(h^2)\) and active third-difference C4
damping, the first-ring JST contribution is third order pointwise for even
variables and second order pointwise for odd radial momentum.  Since the first
annulus has cylindrical measure \(O(h^2)\), those contributions are fourth and
third order, respectively, in a cylindrical-volume-weighted L2 norm.  Thus the
axis skip introduces no singular defect and does not lower the expected global
third order of Skew4+JST.

This is an algebra/unit gate, not a replacement for the R-Z Euler MMS.  The
minimum solver gate remains N=32/64/128 Skew4+JST RHS or short-evolution MMS,
checking global cylindrical L2 order near three and separately recording the
first-ring radial-momentum error.
