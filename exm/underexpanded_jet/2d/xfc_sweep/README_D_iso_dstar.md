# D series — inlet profiles at fixed displacement thickness

Four shifted-tanh inlet profiles that share a compressible displacement
thickness and a discharge coefficient while their peak shear and momentum
thickness span a factor of 2.2. The series exists to test whether the first
shock-cell location follows the mass deficit alone.

The shifted tanh has two parameters and they enter the integral scales
differently:

| scale | depends on |
|---|---|
| `delta_omega = 2a` | `a` alone |
| `theta ~ a/2` | `a` alone |
| `delta*` | both `a` and `s` |

so `delta* = const` is a curve in the `(a, s)` plane along which `a` — and with
it the shear and the momentum deficit — is free to move. The four members walk
that curve. This is the complement of the S series, which fixes `a` and lets
`delta*` run.

## Members

Anchored on the existing case `A2_a170`, which is member 3 and needs no new
directory. Numbering keeps the ordering in `a`, so the gap at D3 is deliberate.

| case | `a` [um] | `s` | `delta_omega` [um] | `theta_c` [um] | `delta*_c` [um] | `C_d` |
|---|---|---|---|---|---|---|
| `D1_a100_s5p196` | 100 | 5.195590 | 200 | 44.24 | 532.9673 | 0.9179 |
| `D2_a130_s3p966` | 130 | 3.965528 | 260 | 57.50 | 532.9672 | 0.9179 |
| `A2_a170` (= D3) | 170 | 3.000000 | 340 | 75.04 | 532.9672 | 0.9180 |
| `D4_a220_s2p284` | 220 | 2.284208 | 440 | 96.40 | 532.9672 | 0.9181 |

`delta*_c` is held to 0.0001 um across the series and `C_d` to 2e-4, so any
movement of `x_fc` across these four is not a mass-flux effect.

`delta_omega` and `theta_c` move together, both being proportional to `a`. The
series therefore does not separate those two from each other; that separation
needs the cross-family pairs of the B series.

## Why the series stops where it does

Two bounds, both hit outside the four members:

- **lip leak.** `f(R-) = (1 - tanh s)/2` grows as `s` falls. D4 sits at
  1.0e-2; below `s = 2.28` the profile no longer closes on the lip and the
  boundary condition stops representing a solid edge.
- **resolution.** The lip is resolved at 24.8 um, so `delta_omega` = 200 um in
  D1 is eight cells. Going thinner would confound a physical difference with an
  under-resolved one.

The forced-L5 band reaches 2.5 mm inside the lip. The thickest member, D4,
reaches `f = 0.99` at `R - 1.008 mm`, so all four clear the band by at least
1.49 mm and the tagging block is unchanged.

## What was edited

Each `prob.h` differs from `A2_a170/prob.h` in exactly three numeric literals,
and nothing else:

| line | template | this series |
|---|---|---|
| 45 | `Real delta_jet = 1.700000e-04;` | `a` |
| 85 | `const Real r0 = p.r_jet - 3.0*p.delta_jet;` | `s` |
| 180 | `const Real r0 = pparm.r_jet - Real(3.0) * pparm.delta_jet;` | `s` |

`GNUmakefile` and `inputs_L5lip` are byte-identical to the template.

The `_LO` variants carry one further line, the same one the existing `_LO`
cases carry: the forced-L5 tag drops the axis strip and keeps the lip band
only.

Comments were deliberately left untouched, following the precedent set by the
S and B series. Two of them are consequently stale in every case of the sweep,
including the pre-existing ones: line 45 still reads `delta/De = 0.01` (true
only for `A2_a170`; `B1_wall680` has carried the same wrong tag since it was
created) and lines 83/168/170 still say the profile is centred at `R-3*delta`.
Correct values are in the table above.

## Running

Same two-phase protocol as the rest of the sweep, driven by
`run_sweep2d_usw2.sh`:

```
make DIM=2 USE_CUDA=TRUE CUDA_ARCH=<80 for A100, 90 for H100> -j
./main2d.ex inputs_L5lip amr.max_level=5 cns.record_stats=0 stop_time=4.0e-4 ...
./main2d.ex inputs_L5lip amr.max_level=5 cns.record_stats=1 amr.restart=<chk> \
            stop_time=2.4e-3 amr.regrid_on_restart=1 amr.regrid_int=100000 ...
```

The statistics window is [0.4, 2.4] ms. `post_restart` zeroes the accumulators,
so phase 2 is a clean average with no startup transient.

The result files to compare are `XFC_<case>_LO__L5.npz`, which is the naming
the existing analysis reads.
