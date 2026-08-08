# LLF and AMR accuracy recovery

This note records the interpretation of the Cartesian MMS and AMR isentropic
vortex studies. It separates the reference methods from new CERISSE
diagnostics.

## Reference basis

| Reference | Principal inviscid method | Relevant verification evidence |
|---|---|---|
| Tin Hang Un | Point-state TENO5, HLLC, and the alternative finite-difference correction. WENO-Z5 is also available. | Smooth AFD accuracy and AMR manufactured-solution tests. |
| Enson Un, CERISSE1 | The same AFD--HLLC framework with WENO-Z5 or TENO5. A separate LLF split-flux branch uses WENO-Z with \(q=1\) and fixed \(\epsilon=10^{-40}\). | Source-level reference for the current AFD and LLF implementations. |
| Monal Patel | Characteristic global Lax--Friedrichs splitting and TENO6. | AMR vortex with fourth-order prolongation, global time stepping, no coarse--fine flux correction, CFL about 0.75, and 50 transits. |
| Omer Al-Hamdi | HLLC--TVD, hybrid DRP, and skew-symmetric methods. | Background for those methods. It is not direct LLF--WENO or LLF--TENO accuracy evidence. |

The AFD--HLLC results are the closest direct comparison with Tin and Enson.
The LLF methods require separate verification.

## LLF--WENO-Z5

The legacy LLF branch retains Enson's WENO-Z exponent \(q=1\). Changing the
exponent to \(q=2\) improved smooth MMS rates but added substantial damping in
the Shu--Osher test. The \(q=2\) diagnostic was therefore rejected.

The accepted candidate changes only the WENO regularisation. It forms one
scale from all characteristic split-flux values on a face stencil and uses
\(\epsilon=10^{-6}S_f^2\). The legacy fixed epsilon remains the default.

The fixed-time MMS uses \(N=16,32,64,128\) and
\(\Delta t\propto\Delta x^2\). The final \(L_2\) orders are:

| Variable | \(64\rightarrow128\) order |
|---|---:|
| Density | 5.019 |
| \(x\) momentum | 5.000 |
| \(y\) momentum | 4.986 |
| Total energy | 5.030 |

In the \(N=200\) Shu--Osher regression, the candidate changes density total
variation by 0.0526 per cent relative to the legacy default. The minimum
density is unchanged to the reported precision.

This is a new CERISSE all-fluid Cartesian option. It is not part of Tin's
principal AFD method and is not yet a general production default.

## LLF--TENO5

Reducing the LLF TENO cutoff from the robust default \(10^{-3}\) to
\(10^{-7}\) restores final MMS orders close to five. The corresponding
Shu--Osher density total variation changes by 7.58 per cent. The
\(10^{-7}\) value is therefore retained only for smooth verification. It is
not adopted as the shock-capturing default.

Tin's AFD--HLLC--TENO5 calculation uses a separate TENO path and a cutoff of
\(10^{-4}\). That value does not define the LLF--TENO5 default.

## AMR isentropic vortex

The earlier CERISSE interface test used limited linear prolongation,
subcycling, and refluxing. Its observed density rates were between about 1.5
and 1.8. This configuration did not match Patel's reported AMR coupling.

The final short-time sensitivity test uses:

- LLF--WENO-Z5 with \(q=1\) and the scaled epsilon candidate;
- AMReX conservative quartic prolongation;
- global time stepping;
- no refluxing;
- one transit at CFL 0.3;
- base grids \(40^2\), \(80^2\), and \(160^2\).

The observed orders are:

| Grid pair | Density \(L_2\) | Velocity \(L_2\) | Pressure \(L_2\) |
|---|---:|---:|---:|
| \(40\rightarrow80\) | 3.343 | 3.694 | 3.785 |
| \(80\rightarrow160\) | 3.923 | 4.727 | 4.566 |

The interpolation matrix shows that limited linear prolongation is an
important error source. The recovered rates also depend on the complete AMR
coupling and the scaled WENO regularisation. This is a Patel-aligned
short-time sensitivity test. It is not a reproduction of Patel's TENO6
calculation and does not establish a formal fourth-order AMR method.

## Result status

- AFD--HLLC--WENO-Z5 and AFD--HLLC--TENO5 provide the direct Tin and Enson
  comparison.
- Scaled-epsilon LLF--WENO-Z5 provides stable fifth-order smooth MMS evidence
  and a nearly neutral Shu--Osher regression. It remains opt-in.
- The \(10^{-7}\) LLF--TENO5 result is a parameter sensitivity result.
- The AMR result demonstrates recovery from the earlier low-order interface
  result under a reference-aligned coupling configuration.
- No IBM or EB production path is changed by these options.

The archived data are in:

- `mms_cartesian/results/llf_wenoz_q1_face_scaled_eps1e6_fixedt_h2_20260729_r1`;
- `mms_cartesian/results/llf_teno5_cutoff1e7_fixedt_h2_20260729_r1`;
- `shuosher/results/llf_weno_q1_face_epsilon_shuosher_20260729_r1`;
- `shuosher/results/llf_reference_parameter_regression_20260729_r1`;
- `isentropic_vortex/results/amr_reference_config_q1_gts_reflux0_20260729_r1`.
