# Pure-GP R-Z dual-channel algebra regression

This is a fast, test-only oracle for the shared-background-face positivity
limiter. It does not compile into Cerisse and does not use cut cells, embedded
boundaries, area fractions, aggregate cells, or case-specific geometry.

The test uses complete R-Z background annuli and checks:

- the paired radial-momentum operator
  `D_metric(F - P) + D_cart(P)`;
- the nonzero pressure correction at the zero-area axis face;
- the explicit undivided metric-flux correction for other components at the
  axis;
- ordinary full-face metric corrections away from the axis;
- one shared face value as seen by both adjacent cells;
- identical limiter coefficients for the complete flux `F` and pressure
  companion `P`;
- the `theta=0` all-low and `theta=1` all-high endpoints;
- direct selected-operator and face-correction identities to floating-point
  roundoff, including different limiter coefficients on a cell's two faces.

Run from the repository root:

```bash
python3 tst/ibm/positivity/rz_axis_shared_face/test_rz_dual_channel_algebra.py -v
```

This regression verifies the local algebra only. It is not a substitute for
the required pure-GP field, curved-wall MMS/BI-pressure, inverse-MOC, and AMR
reflux qualification gates.
