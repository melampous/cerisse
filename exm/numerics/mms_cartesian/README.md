# Cartesian Euler and Navier-Stokes MMS

This directory contains a two-dimensional, time-dependent manufactured
solution for the point-sample finite-difference formulation used by CERISSE.
The domain is periodic, so the measured error is not contaminated by a
physical-boundary closure.

The inviscid method is selected with `CERISSE_MMS_EULER_SCHEME`. The
supported values are `llf-wenoz5`, `llf-teno5`, `llf-teno6`,
`afd-hllc-wenoz5`, and `afd-hllc-teno5`. Explicit `-old` selectors retain the
archived implementation for controlled regression only. The default
`llf-wenoz5` selector uses the compact shared driver in `Weno.h`. The
Navier-Stokes build adds centred viscous stress and Fourier heat flux.
`CERISSE_MMS_VISCOUS_ORDER` selects order 2, 4, or 6. The default is 2.
`CERISSE_MMS_TRANSPORT_MODEL` selects `constant` or `sutherland`. The default
is `constant`. The manufactured fields activate normal stress, shear stress,
dilatation, viscous work, and heat conduction. The Sutherland option also
tests products involving variable viscosity and thermal conductivity.

For Cartesian, single-species, non-LES calculations, the fourth-order method
first constructs a fourth-order pointwise physical viscous flux at each face.
It forms the complete stress, viscous-work, and heat-conduction products before
converting that flux to the conservative finite-difference auxiliary flux.
The RZ, GP-IBM, LES, and multispecies paths retain their existing operators.

The formal study uses `N=16,32,64,128` and a fixed final time. The number of
steps increases by four after every grid refinement. Thus `dt` is
proportional to `dx^2`, and the global SSPRK(3,3) error is sixth order in
`dx`.

Run both physics cases with:

```bash
CERISSE_MMS_RESULTS_ROOT=results/formal_cartesian_mms ./run_study.sh
```

The formal Cartesian result used for the thesis is generated with:

```bash
CERISSE_MMS_RESULTS_ROOT=results/new_cartesian_afd_mms \
CERISSE_MMS_EULER_SCHEME=afd-hllc-wenoz5 \
./run_study.sh
```

A variable-transport fourth-order study is generated with:

```bash
CERISSE_MMS_RESULTS_ROOT=results/cartesian_ns_visc4_sutherland \
CERISSE_MMS_EULER_SCHEME=afd-hllc-wenoz5 \
CERISSE_MMS_PHYSICS=navier-stokes \
CERISSE_MMS_VISCOUS_ORDER=4 \
CERISSE_MMS_TRANSPORT_MODEL=sutherland \
./run_study.sh
```

For the AFD method, the runner sets and records the flux correction, shock
fallback, smoothness threshold, and TENO cut-off explicitly.

The analyzer reports pointwise conservative-variable errors in the discrete
`L1`, `L2`, and `Linf` norms. Do not compare these results with legacy MMS
data generated from separately averaged primitive variables.

## Relation to Tin and Enson

Tin's principal high-order inviscid method uses point-state TENO5
interpolation. WENO-Z5 is also available for comparison. The corresponding
CERISSE1 source provides both choices. The reconstructed states are passed to
HLLC. An alternative finite-difference correction then removes the leading
derivative errors in smooth regions. The
`afd-hllc-wenoz5` and `afd-hllc-teno5` selectors represent this method in the
present code.

The LLF selectors are a separate finite-difference formulation. They apply
characteristic local Lax--Friedrichs flux splitting before WENO or TENO
reconstruction. The LLF result must therefore be verified separately. It
must not be described as the Tin AFD--HLLC method.

Tin's thesis writes the WENO-Z exponent as \(q=1\), whereas the archived
CERISSE1 source uses the squared ratio, equivalent to \(q=2\). The present LLF
implementation defaults to \(q=1\). A face-scaled WENO epsilon can be supplied
to the runner with `CERISSE_MMS_LLF_WENO_EPSILON_RELATIVE`. It is an explicit
test option and is recorded in the result manifest.
