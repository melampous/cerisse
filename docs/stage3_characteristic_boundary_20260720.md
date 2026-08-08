# Stage 3 characteristic-boundary verification

Date: 2026-07-20

## Decision

\[
\boxed{\text{Stage 3A: PASS}}
\]

\[
\boxed{
\text{Stage 3B: operator-level PASS; boundary-system-level IN PROGRESS}
}
\]

Stage 3 as a whole remains **IN PROGRESS**.  The normal characteristic inlet
and outlet models are suitable for the straight-buffer internal-flow tests
defined below.  The Giles/GC-NSCBC implementation is a mathematically and
discretely verified single-outlet prototype; it is not yet a general external
far-field boundary.

No IBM geometry, BI-constrained reconstruction, ghost-cell extension, WENO
weight, or positivity-limiter threshold was changed in Stage 3B.

## Frozen numerical scope

The verified Stage 3B configuration is:

- two-dimensional Cartesian coordinates;
- one uniform level-0 mesh;
- single-species calorically perfect-gas Euler equations;
- LLF-WENO-Z5 and SSPRK(4,3);
- exactly one type-2 or type-3 subsonic outflow boundary;
- periodicity in the boundary-tangential direction;
- three physical ghost-cell layers;
- positive outward mean velocity;
- opt-in and default-off.

The implementation terminates for R--Z coordinates, refined AMR levels,
reverse flow, more than one active NSCBC boundary, and non-periodic tangential
boundaries.  These are explicit scope guards, not claimed capabilities.

## Stage 3A: normal characteristic models

For an outward unit normal \(\boldsymbol n\), the Euler eigenvalues are

\[
u_n-c,\qquad u_n,\qquad u_n,\qquad u_n+c.
\]

The accepted normal models are:

- subsonic inlet: prescribe \(p_t,T_t\) and flow direction while retaining the
  outgoing acoustic invariant from the interior;
- subsonic outlet: prescribe only the incoming acoustic information through
  \(p_b\), while entropy, tangential velocity, and the outgoing acoustic mode
  come from the interior;
- supersonic inlet: prescribe the full state;
- supersonic outlet: extrapolate all characteristics and ignore \(p_b\);
- linearized characteristic far field: incoming modes use the target state and
  outgoing modes use the interior state.

The ideal-gas total/static relations are

\[
\frac{T_t}{T}=1+\frac{\gamma-1}{2}M^2,
\qquad
\frac{p_t}{p}=\left(\frac{T_t}{T}\right)^{\gamma/(\gamma-1)}.
\]

The implementation does not prescribe static pressure, velocity, mass flow,
and total state simultaneously at a subsonic inlet.

### Stage 3A acceptance evidence

| Test | Result | Evidence |
|---|---|---|
| Subsonic oblique uniform flow, \(M=0.3\), \(23^\circ\) | PASS | maximum conservative RHS \(1.172\times10^{-13}\) |
| Supersonic oblique uniform flow, \(M=2\), \(23^\circ\) | PASS | maximum conservative RHS \(3.411\times10^{-13}\) |
| Subsonic \(p_t,T_t,p_b\) duct | PASS | mass-flow relative error \(-1.848\times10^{-6}\); entropy span \(1.992\times10^{-14}\) |
| Supersonic outlet back-pressure independence | PASS | \(p_b/p_0=0.2\) and 5.0 give bitwise-identical fields |
| Normal outgoing acoustic packet | PASS | incoming/outgoing characteristic ratio \(O(10^{-5})\) |
| Reverse flow | PASS, fail-closed | explicit non-success status; no state clipping |

Stage 3A permits smooth internal-nozzle verification when the inlet and outlet
are placed in straight, approximately one-dimensional buffer sections and
buffer-length sensitivity is reported.

## Stage 3B mathematical closure

For frozen local mean flow, define the dimensional characteristic variables

\[
c_2=\rho c\,\delta u_t,
\qquad
c_4=\delta p-\rho c\,\delta u_n.
\]

The incoming acoustic equation from the full two-dimensional linearized Euler
system is

\[
\partial_t c_4
+(u_n-c)\partial_n c_4
+c\,\partial_s c_2
+u_t\partial_s c_4=0.
\]

The implemented boundary closure is Giles' second-order two-dimensional
unsteady outflow condition, Eq. (174):

\[
\partial_t c_4
+u_n\partial_s c_2
+u_t\partial_s c_4=0.
\]

Equating the two relations at the boundary gives

\[
\partial_n c_4=\partial_s c_2
=\rho c\,\partial_s u_t.
\]

Thus the interior one-sided derivative supplies the outgoing acoustic,
entropy, and tangential characteristics; only the incoming acoustic normal
derivative is replaced.  Giles Eq. (175) is intentionally not used because
its higher acoustic order introduces first-order vorticity reflection, which
is undesirable for nozzle wakes and SRP shear layers.

This method is called a **Giles-type multidimensional Euler characteristic
outflow**.  It is not described as an exact non-reflecting boundary and is not
misattributed to the original one-dimensional Poinsot--Lele LODI closure.

## Ghost-cell finite-volume mapping

Each SSPRK stage executes the following sequence:

1. complete the same-level and periodic halo exchange;
2. evaluate a second-order one-sided normal derivative at the
   boundary-adjacent interior cell;
3. evaluate the centered tangential velocity divergence;
4. preserve outgoing characteristic derivatives and impose Giles Eq. (174)
   on the incoming acoustic derivative;
5. generate exactly three normal ghost layers with the GC-NSCBC recurrences
   corresponding to Motheau--Almgren--Bell Eqs. (31)--(33);
6. reconstruct an EOS-consistent conservative ghost state from
   \((\rho,p,u_n,u_t)\);
7. reject non-finite, non-positive, or reverse-flow states without clipping.

The published fourth-layer recurrence was not adopted: Cerisse requires three
layers for this LLF-WENO-Z5 path, and the reported fourth expression does not
pass the same linear-polynomial reproduction check.  The implemented three
layers reproduce constant, linear, and quadratic normal profiles to roundoff.

High- and low-order flux paths read the same stage-local physical boundary
state.  Communication overlap is rejected while NSCBC is active, avoiding a
path that would bypass the stage-local fill.

## Stage 3B verification

### Algebra and polynomial consistency

The CPU startup oracle verifies:

- total-condition inlet recovery;
- pressure-outlet preservation of \(J^+\), entropy, and tangential velocity;
- supersonic outlet independence from \(p_b\);
- explicit reverse-flow status;
- far-field uniform-state recovery;
- three-layer linear and quadratic GC-NSCBC reproduction;
- exact Giles derivative closure for a divergence-free convected vorticity
  mode.

All checks pass.  CUDA coverage is provided by the compiled stage-local kernel
and the field-level CPU/CUDA comparison.

### Boundary-local divergence-free vorticity mode

This is the decisive local transverse-consistency test.  The exact
streamfunction construction has zero pressure perturbation and satisfies the
continuum Giles closure.  The pressure-conversion error is:

| \(N\) | Pressure conversion | Incoming acoustic ratio | Adjacent pressure order |
|---:|---:|---:|---:|
| 48 | \(5.62486\times10^{-3}\) | \(3.00616\times10^{-2}\) | -- |
| 96 | \(1.18334\times10^{-3}\) | \(5.51448\times10^{-3}\) | 2.25 |
| 192 | \(2.80294\times10^{-4}\) | \(1.22742\times10^{-3}\) | 2.08 |
| 384 | \(6.73718\times10^{-5}\) | \(2.92145\times10^{-4}\) | 2.06 |

The local spurious acoustic conversion therefore converges at approximately
second order.  A tenfold perturbation-amplitude reduction gives the same
normalized result, confirming linear rather than nonlinear contamination.

The earlier long-exit pressure-conversion plateau of about \(1.2\times10^{-3}\)
is not used as the local consistency verdict: that final-time metric also
contains bulk propagation, packet departure, and characteristic-projection
errors accumulated before and after the boundary interaction.

### Oblique acoustic packets

At \(M_n=0.3\), the measured final-time incoming characteristic ratios are:

| Incidence | \(N=96\) | \(N=192\) | Giles single-frequency pressure-amplitude coefficient |
|---:|---:|---:|---:|
| \(15^\circ\) | 0.001086 | 0.005268 | 0.009447 |
| \(30^\circ\) | 0.017325 | 0.030852 | 0.040590 |
| \(45^\circ\) | 0.080464 | 0.084439 | 0.103107 |
| \(60^\circ\) | 0.171590 | 0.176664 | 0.218225 |

The angular trend and nonzero finite-angle reflection agree with the expected
behavior of a finite-order local outflow condition.  The measured quantity is
a final-time domain \(L_2\) wave-packet ratio, whereas Giles' expression is a
single-frequency boundary pressure-amplitude coefficient; their numerical
difference is not an implementation error estimate.  In particular, changing
\(M_n\) changes the reflected packet group velocity and the portion remaining
inside the sampled domain.

At \(30^\circ\), the \(N=96/192\) incoming ratios are 0.05310/0.05638 for
\(M_n=0.1\), 0.01732/0.03085 for \(M_n=0.3\), and
0.001080/0.001830 for \(M_n=0.6\).  The corresponding single-frequency Giles
coefficients are 0.06007, 0.04059, and 0.01915.  This confirms Mach-number
sensitivity but does not replace a boundary-time Fourier-amplitude test.

For \(+30^\circ\) and \(-30^\circ\), the corrected fixture gives matching
metrics to about \(3.2\times10^{-13}\) absolute at \(N=96\).  This verifies
the sign symmetry of the tangential derivative for the tested boundary.

### Normal frequency and convective modes

The normal acoustic incoming ratios at \(N=96\) are approximately
\(1.59\times10^{-8}\), \(4.29\times10^{-8}\), and \(8.19\times10^{-7}\)
for the low-, medium-, and high-frequency packets.  The transverse correction
therefore does not regress the accepted near-normal result.

For convective disturbances at \(N=96\):

- entropy-wave pressure conversion: \(1.87\times10^{-10}\);
- one-dimensional tangential-shear pressure conversion:
  \(8.32\times10^{-9}\).

These results show no material spurious acoustic conversion in the tested
one-dimensional convective modes.  The divergence-free oblique vorticity mode
is covered separately by the second-order local convergence table above.

### Parallel and failure-path checks

| Check | Result |
|---|---|
| CPU versus CUDA, \(N=96\) local vorticity test | maximum field difference \(2.44\times10^{-15}\) |
| One rank/one large FAB versus two ranks/32-cell FABs | bitwise-identical \(\rho,p,u,v\) |
| Reverse flow at the active outlet | exits with code 6 and diagnostic code 1; no clipping |
| Non-periodic tangential boundary | exits with code 6 before use; physical-corner ownership explicitly uncertified |
| R--Z or AMR refined level | explicitly rejected by the implementation |

## Interpretation

The original normal LODI result exhibited the expected quasi-one-dimensional
failure: normal waves were absorbed well, while oblique acoustic and vortical
modes produced finite reflection or mode conversion.  The Giles/GC-NSCBC
prototype removes the leading local inconsistency for a divergence-free
vorticity mode and gives approximately second-order local convergence.

It does not make a finite local boundary exactly transparent.  Reflection at
large incidence angle remains expected, and the \(60^\circ\) result must not be
reported as zero-reflection.  If the resulting reflection exceeds the SRP
error budget, the physically defensible remedies are a larger outer domain or
an independently verified sponge layer, not empirical retuning of the Giles
coefficient.

## Remaining Stage 3B gates

Stage 3B cannot be closed for external flow until all of the following pass:

1. all four Cartesian boundary orientations using positive and negative
   tangential wavenumbers;
2. deterministic ownership and one-sided tangential derivatives at physical
   corners;
3. simultaneous multi-boundary far-field operation;
4. separation of mean back-pressure relaxation from nonzero-frequency acoustic
   absorption for the stage-local type-3 outlet;
5. boundary-time Fourier extraction of pressure-amplitude reflection over the
   required Mach, frequency, incidence-angle, and cells-per-wavelength range;
6. repeated CPU/CUDA/FAB/MPI and fail-closed tests after that extension.

Consequently:

- smooth internal nozzle with straight buffer sections: **GO under Stage 3A**;
- shocked internal nozzle: **CONDITIONAL GO**, with shock-to-outlet and outlet-
  position sensitivity;
- circle/sphere quantitative external flow: **NO-GO pending Stage 3B closure**;
- underexpanded jet and SRP external far field: **NO-GO pending Stage 3B
  closure**.

## Evidence locations

- `temp/stage3_characteristic_bc/gc_giles2_vortex_boundary_48_96/`
- `temp/stage3_characteristic_bc/gc_giles2_vortex_boundary_96_192/`
- `temp/stage3_characteristic_bc/gc_giles2_vortex_boundary_192_384_cuda/`
- `temp/stage3_characteristic_bc/gc_giles2_angle_96_192/`
- `temp/stage3_characteristic_bc/gc_giles2_angle_pm30_final/`
- `temp/stage3_characteristic_bc/gc_giles2_mach01_96_192/`
- `temp/stage3_characteristic_bc/gc_giles2_mach06_96_192/`
- `temp/stage3_characteristic_bc/gc_giles2_frequency_48_96/`
- `temp/stage3_characteristic_bc/gc_giles2_entropy_shear_final/`
- `temp/stage3_characteristic_bc/gc_giles2_vortex_boundary_cpu_final/`
- `temp/stage3_characteristic_bc/gc_giles2_vortex_boundary_cuda_final/`
- `temp/stage3_characteristic_bc/gc_giles2_vortex_boundary_mpi2_final/`

## Provenance

| Item | Commit / SHA-256 |
|---|---|
| Cerisse HEAD | `dc35646d8eb90aaa362a3292a9d702e30c9d5af5` plus the documented dirty-worktree changes |
| AMReX | `bd922c6216e0a734f3b1cf0ca73e7d669b90f3ef` |
| `src/set/nscbc.h` | `7a4396f7073146674cfae404dfe1b3c86f2b7d86ba1959b415112b6cd4ef1494` |
| `src/set/nscbc_ghost.cpp` | `0a411555504c0231ed74fdec78a404cde1caa681d3573c8f06ecc2ed7312d39f` |
| `src/CNS.cpp` | `bdffbdf8a8747e9272fb04e74692b8cad70f788a3b42f8d08ae8792625a84078` |
| `src/CNS.h` | `f97fe9226c89fac5d7f272eb6fca181fed6e6d696ee5dff3fba98be32516764c` |
| `src/tim/advance.cpp` | `948205fff670ae141163f304db61e7bd62e1130e82a968ca7f8c36ab85e5917c` |
| `exm/numerics/bc_native/prob.h` | `9301d202f4f5909a3e92ba2dd5d4b7e42c44f83b1112d4eec6c43b8d4c2ef812` |
| `tools/run_native_bc_validation.py` | `385562f83641f223d61078f8059003cffe0ca310c743c96161954b550d3d1582` |
| CPU executable | `06e9433b448bee540340236d814693cb4ce62bdd3e5cbfad3f9527cef525b8bc` |
| CUDA executable | `d850609a278bf30c761c6b3bc6d195cd65a03d3cdc5b0dbd495306a5de186a8f` |

Toolchain: GCC 13.3.0, Open MPI 4.1.6, CUDA 12.0.140, NVIDIA RTX
4090 Laptop GPU, compute capability 8.9.

## References

- M. B. Giles, *Non-Reflecting Boundary Conditions for Euler Equation
  Calculations*, CFDL-TR-88-1, especially Eqs. (160), (174), and (175):
  <https://people.maths.ox.ac.uk/gilesm/files/bcs.pdf>
- E. Motheau, A. S. Almgren, and J. B. Bell, *Navier--Stokes Characteristic
  Boundary Conditions Using Ghost Cells*, AIAA Journal 55 (2017):
  <https://doi.org/10.2514/1.J055885>
- G. Lodato, P. Domingo, and L. Vervisch, multidimensional NSCBC treatment:
  <https://web.stanford.edu/group/ctr/ResBriefs/2010/12_lodato.pdf>
- T. J. Poinsot and S. K. Lele, *Boundary Conditions for Direct Simulations of
  Compressible Viscous Flows*, JCP 101 (1992):
  <https://doi.org/10.1016/0021-9991(92)90046-2>
