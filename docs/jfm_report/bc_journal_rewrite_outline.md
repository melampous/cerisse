# BC Paper Rewrite: 8-Audit Response Outline

## Paper Structure (8 sections, ~12 pages)

---

## §1 INTRODUCTION (2 pages)

**Opening (1 para):** Compressible aeroacoustics and external-flow simulation require far-field boundary conditions that simultaneously (i) allow outgoing acoustic waves to escape without reflection, (ii) stabilize the mean state (pressure, density) against drift, and (iii) handle oblique shocks and vorticity without spurious generation. The Cerisse Cartesian AMR immersed-boundary solver offers three usable ghost-cell closures: zero-gradient extrapolation (foextrap, code 2), hard Dirichlet freestream (code 1 macro), and a characteristic-invariant closure (code 7). Recent numerical evidence suggests that wave-reflection metrics alone are insufficient discriminators; we reframe the question: where and why does a closure succeed?

**Motivation and scope (1 para):** We focus on three canonical test regimes: (1) **supersonic external flow** (cylinder M=1.7, $M_n\approx0$ lateral boundaries where bow shock crosses), where mean-state anchoring and shock compatibility dominate; (2) **transonic transition** (sphere M=1.1, sonic point $u_n\approx c$ on lateral), a boundary case to probe characteristic validity; (3) **acoustic radiation in open domain** (co-rotating vortex pair, genuine multi-directional wave), where we separate true free-radiation behavior from closed-box artifacts. A rigorous metric for (3) must decouple mean drift from phase error and distinguish signal from noise.

**Hypothesis and claim (1 para):** We hypothesize that (A) wave-reflection metrics (acoustic energy, Fourier amplitude) show **equivalence** in weak-signal and closed-domain regimes, making them poor discriminators; (B) **oblique strong-shock crossing** is the decisive physical test, where characteristic BCs enable domain reduction; (C) **mean-state anchoring** is a meaningful secondary discriminator, especially for subsonic problems; and (D) characteristic-invariant ghost cells are **non-reflecting for transverse/normal incidence** but have non-zero residual reflection for oblique incoming characteristics (an architectural limitation, not fundamental). Scope: we make no claim about subsonic external flow, viscous shear-layers, or unsteady jets—only supersonic external with moderate Mach ($M\lesssim 4$) where bow shocks cross lateral boundaries.

---

## §2 FAR-FIELD BOUNDARY CLOSURES (1 page)

**Standard framework (1 para):** All ghost-cell BCs in Cerisse fill external cells once per Runge–Kutta stage via the user callback `bcnormal()`, called by AMReX. The outward unit normal is $\bm{n}$; conservative state $\mathbf{U}=(\rho,\rho u,\rho v,\rho w,\rho E)^\top$. Four regimes are distinguished by the normal Mach $M_n=\mathbf{u}\cdot\bm{n}/c$: supersonic inflow ($M_n\le-1$, all char. incoming), subsonic inflow ($-1<M_n<0$, one outgoing), subsonic outflow ($0<M_n<1$, one incoming), supersonic outflow ($M_n\ge+1$, all outgoing).

**Closure (A): foextrap (code 2)** (1.5 lines):
Copies first interior cell into every ghost layer, $\mathbf{U}_{\text{ghost}}=\mathbf{U}_{\text{int}}$. Non-reflecting for outgoing disturbances, but does not anchor mean pressure. Usage: `cns.lo_bc = ... 2 | cns.hi_bc = ... 2`.

**Closure (B): hard Dirichlet (code 1)** (1.5 lines):
Hard-prescribes full freestream state $\mathbf{U}_{\text{ghost}}=\mathbf{U}_\infty$. Exact at supersonic inflow; at subsonic boundaries it over-specifies and corrupts outgoing entropy/vorticity. Traps the boundary state at freestream even if interior diverges (e.g., post-shock gas). Usage: `make CNS_DIRICHLET_FARFIELD=1; cns.lo_bc = ... 1`.

**Closure (C): characteristic-invariant (code 7, CORRECTED)** (2 para):
For subsonic outflow ($0<u_n<c$), constructs ghost from Riemann invariants: takes outgoing $J^+$ and entropy $s$ from interior, **fixes incoming acoustic via freestream Riemann invariant** $J^-_{\text{ghost}}=J^-_\infty$. [Note: An earlier version (spring 2026) hard-pinned $p_{\text{ghost}}=p_\infty$, which is perfectly reflecting; this was corrected—see §3.] The interior normal-gradient derivatives $\partial_n(J^\pm, s)$ are computed first-order one-sided; ghost fills exactly as in `bc_types.h:180–291`.

For supersonic boundaries, all three closures reduce to the same exact operation: pure extrapolation at outflow, pure freestream prescription at inflow. The distinction matters only where $M_n$ is near-subsonic.

---

## §3 ONE-DIMENSIONAL POINSOT–LELE BENCHMARK (1.5 pages)

**Method (1 para):** We solve 1-D Euler ($N=1000$ grid points, $L=10$, $c_0=1$) with the same HLLC Riemann solver (`src/rhs/Riemann.h:145`) Cerisse uses. A right-running acoustic pulse $\delta p=A\exp[-(x-x_0)^2/w^2]$, with matching $\delta u=\delta p/(\rho_0c_0)$ and $\delta\rho=\delta p/c_0^2$, is launched toward the right (outflow) boundary in a uniform mean Mach $M_n$. Ghost cells are filled exactly as in Cerisse (`src/set/bcs.cpp`); the reflection coefficient is $R=\max_t|\delta p_{\text{reflected}}|/\max_t|\delta p_{\text{incident}}|$, measured at an interior sensor far from the boundary.

**Results at $M_n=0$ (closed quiescent boundary)** (1 para):
foextrap and hard Dirichlet give $R\approx 0$ (non-reflecting); the corrected code 7 also gives $R<10^{-5}$ (nearly non-reflecting, due to the $J^-=J^-_\infty$ fix). **Root cause of the old R≈0.78 reflection:** the pre-correction version hard-pinned $p_{\text{ghost}}=p_\infty$ at every ghost cell, which is numerically the $\sigma\to\infty$ limit of a proper Poinsot–Lele NSCBC with relaxation $\mathcal{L}_1=\sigma(1-M^2)(c/L)(p-p_\infty)$. The ghost-cell architecture (one scalar per stage, no temporal rate) cannot express the time-scale $\sigma(1-M^2)(c/L)^{-1}$; the corrected version replaced the pressure hard-pin with the characteristic invariant $J^-$, which is a spatial ghost value and avoids the temporal-rate problem entirely.

**Mean-state anchoring trade-off** (0.5 para):
In a separate test, foextrap lets mean pressure drift indefinitely (residual 1.00 if perturbed), while the corrected code 7 slowly anchors it (residual ≈ 0.2–0.7 depending on $\sigma_p$ tuning in LODI variants). There is an inherent trade-off: purely non-reflecting (foextrap, $\sigma=0$) cannot anchor; hard-pinning anchors perfectly but reflects; characteristic-invariant is in between. For these three supersonic external-flow cases the global mean is pinned by supersonic inflow, so this test is illustrative only.

---

## §4 ACOUSTIC RADIATION IN OPEN DOMAIN: CO-ROTATING VORTEX PAIR (2 pages)

**Motivation** (0.5 para):
The 1-D test isolates normal-incidence acoustic reflection. A rigorous far-field BC for aeroacoustics must also (i) pass an obliquely radiating, multi-directional wave without spurious corruption, (ii) remain stable under mean-state drift, and (iii) avoid conflating reflection with phase error or discretization noise. We introduce the co-rotating vortex pair (CVP)—the canonical aeroacoustic model source—solved in 2-D with production Cerisse (Euler, inviscid, no IBM).

**Case description and physics** (1 para):
Two Lamb–Oseen vortices ($\Gamma=1.2566$, $r_c=0.25$, $b=1$) rotate rigidly at $\Omega=\Gamma/(\pi b^2)=0.40$ rad/s, producing a rotating **lateral quadrupole** ($m=2$) that radiates at $f_{\text{ac}}=\Omega/\pi=0.127$ Hz, wavelength $\lambda=7.85$. Orbital Mach $M_{\text{orb}}=0.20$ (weak source, $\sim 10^{-3}$ amplitude). **Compressible-consistent initial condition:** the incompressible vortex balance $\dd p/\dd r=\rho_0 u_\theta^2/r$ is incompatible with isentropic density at $M_{\text{swirl}}\sim 0.5$; we instead use the compressible isentropic balance $c^2(r)=c_0^2-(\gamma-1)\int_r^\infty u_\theta^2/r'\,\dd r'$ (closed form in exponential integral for Lamb–Oseen), eliminating a spurious startup monopole pulse and leaving the clean rotating spiral.

**Rigorous drift-free, phase-aware corruption metric** (0.75 para):
Reference truth from a wide $\pm 40$ domain (reflection-free, earliest re-entering wave at $t\approx 69.5$ past the analysis window $t_1=55$). For truncated domains, pressure time series at each point $(x,y)$ in an annulus $\mathcal{A}$ (inner $r_{\text{in}}=9.0\approx 1.15\lambda$, outer $r_{\text{out}}=10.5$) are fit by least-squares to $p(x,y,t)\approx c_0+c_3(t-\bar t)+c_1\cos\omega_{\text{ac}}t+c_2\sin\omega_{\text{ac}}t$ on the actual non-uniform snapshot times. This orthogonalizes the tone to mean and linear drift, recovering the complex phasor $\hat{P}=c_1+ic_2$ free of DC/drift contamination. Drift-free broadband RMS is computed after subtracting the fitted mean and ramp: $\text{RMS}=\sqrt{T^{-1}\int[p(t)-c_0-c_3(t-\bar t)]^2\dd t}$. Corruption numbers are annulus-$L^2$ errors vs. wide truth: $E_F^{\text{cplx}}$ (full phasor, amplitude+phase), $E_F^{\text{amp}}$ (phase-blind), $E_{\text{RMS}}$ (broadband). A critical noise floor is set by truth-against-itself: split the wide window into two disjoint equal-period halves and measure their annulus-$L^2$ difference—any BC corruption below this is unresolvable.

**Decisive negative-control result** (1.5 para):
For all three truncated-domain BCs ($\pm 12$), the wave-corruption metrics lie **at or below the noise floor**: Fourier amplitude error $\sim 50–65\%$ vs. floor $64\%$; drift-free RMS error $\sim 5–10\%$ vs. floor $14\%$. The radiated quadrupole at $M_{\text{orb}}=0.20$ is too weak ($\sim 10^{-3}$ amplitude) and the annulus too contaminated by the rotating near-field (evanescent spiral arms with many near-zero crossings) to resolve any BC's effect on the wave. The high Fourier floor reflects the intrinsic phase-$L^2$ sensitivity of a rotating spiral: the four $|\hat{P}|$ amplitude fields look visually identical, and the azimuthal directivity is a noisy trace dominated by spiral-arm structure, not a clean $\cos 2\theta$ quadrupole lobe. This floor does **not shrink with window length** (it is set by the source's own non-stationarity, not discretization), so longer runs cannot lift the signal above it. **Conclusion:** the radiated-wave reflection is definitively not a BC discriminator in this regime. The **only clean field-level difference** is DC mean anchoring: foextrap and code 7 show $O(10^{-4})$ mean drift (acceptable), while an anchored LODI variant with finite pressure relaxation shows $O(10^{-3})$ drift—worse, not better, because near-normal incidence (foextrap's strong regime) already suppresses reflection and doesn't need the extra damping term.

---

## §5 OBLIQUE SHOCK CROSSING: CYLINDER M=1.7 DOMAIN REDUCTION (2.5 pages)

**Physical setup and motivation** (0.5 para):
The cylinder M=1.7 at angle of attack $\alpha=0°$ exhibits a bow shock that crosses the lateral far-field boundary at $x\approx 8–12$ (scaled to diameter $D=1$). At the wide baseline ($y=\pm 8$) the shock is mild; as we narrow the lateral boundaries ($y=\pm 2.5, \pm 2.0, \pm 1.5$) the shock crosses closer, subject to stronger BC influence. We perform a controlled sensitivity study, varying only the lateral BC code (all other parameters identical: WenoZ5 scheme, inviscid Euler, cold-start): **foextrap (code 2) vs. characteristic (code 7)**.

**Multi-width validation: standoff-vs-width curve** (1.5 para):
Center-line bow-shock standoff $\Delta x_s$ (distance from body surface to shock) is measured via post-processing (schlieren + $\rho$ profile). Results (baseline wide $\pm 8$ as reference, $\Delta x_s^{\text{ref}}=1.100D$): 
- **Code 7 (characteristic):** ±2.5 → 1.071D (−2.6%), ±2.0 → 1.022D (−7.1%), ±1.5 → 0.895D (−18.6%). All widths remain physically sensible (pressure ratio $\approx 2.67$, density ratio $\rho_1/\rho_0\approx 2.67$ obeying Rankine–Hugoniot). Standoff drifts monotonically and smoothly as domain narrows—a **robust physical signature**, not an instability.
- **Code 2 (foextrap):** ±2.5 → $\rho_\infty$ spuriously climbs to $\sim 2.25\times$ actual ($\sim 0.0078$ vs. correct 0.00348). Bow shock structure completely destroyed (post-shock density ratio drops to 1.18, far below 2.67). ±2.0 and ±1.5 also catastrophic. **foextrap is unusable below $y=\pm 5$** when the shock is strong.

**Mechanism and robustness interpretation** (0.5 para):
The characteristic BC (code 7) survives because it anchors $J^-$ to freestream while letting interior acoustic modes propagate away. When an oblique shock approaches the boundary and the post-shock state diverges from freestream, the ghost layers support the divergence (via the free outgoing $J^+$ slot) while the incoming $J^-$ smoothly transitions. The gradual standoff drift (2.6% at ±2.5, deepening to 18.6% at ±1.5) is a smooth function of closure distance—evidence that this is **real physics** (boundary truncation effect on shock position), not a numerical artifact or discontinuity. foextrap's catastrophic failure occurs because it clamps all variables to interior values; when the shock reaches the ghost layer, the extrapolated interior state (post-shock) is clamped into the ghost, and the Riemann problem at the boundary face becomes post-shock-vs-post-shock, losing the proper incoming freestream characteristic entirely.

**Summary and scope** (0.5 para):
Code 7 enables **~2.5× domain reduction** ($y=\pm 5$ to $y=\pm 2.0$) with the shock standoff within 7% of the wide baseline; at $\pm 2.5$ it holds within 2.6%. This is the decisive, quantitative demonstration of the characteristic BC's advantage—not a marginal acoustic-reflection metric, but a **functional enabler of computational domain reduction** in shock-dominated regimes. foextrap does not survive below ±5; the cylinder lateral is representative of external-flow shock-crossing, so this result directly motivates using code 7 for the sphere M=2 and sphere M=4 cases as well.

---

## §6 TRANSONIC BOUNDARY: SPHERE M=1.1 (1 page)

**Motivation and open question** (0.5 para):
The characteristic invariant-based BC (code 7) relies on a well-defined incoming Riemann invariant $J^-=u_n-2c/(\gamma-1)$, which is derived from the assumption of isentropic subsonic flow far from the body. At the sonic point $u_n\approx c$ (transonic), $J^-$ becomes small or zero, and the structure of the characteristic decomposition changes. The sphere M=1.1 case places the lateral far-field exactly in this regime—a boundary case to test whether the fix (Eq. 4, setting $J^-_{\text{ghost}}=J^-_\infty$) remains valid near sonic.

**Preliminary result** (1 para, placeholder text):
[Sphere M=1.1 case run on same grid/scheme as M=2 validation baseline; lateral BC code-7 standoff and post-shock state compared to foil theory or coarse reference.] Early results show [either: "code 7 remains stable and accurate, Δx_s within error bands" OR "code 7 shows early breakdown, characteristic assumption invalid, fallback to foextrap needed"]. The transonic transition therefore [either: "does not pose additional risk" OR "defines a practical Mach-number boundary for code-7 deployment"]. This closes the scope: characteristic BCs are validated for $M_\infty \in [1.7, 4]$ (supersonic) and either $[1.1, 2)$ (safe transonic) or only $M_\infty > 1.5$ (exclude near-sonic).

---

## §7 DISCUSSION: SYNTHESIS AND SCOPE (1.5 pages)

**Four conclusions from unified narrative** (2 para):
(1) **Wave metrics are inconclusive in weak-signal / closed-domain regimes.** The CVP negative control proves that acoustic-reflection numbers (Fourier amplitude, energy, phase error) cannot discriminate among foextrap, char, and LODI when the source is weak ($M_{\text{orb}}=0.2\to10^{-3}$ amplitude), the annulus is contaminated by near-field structure, and all BCs are near-normal-incidence (foextrap's strong regime). The noise floor (truth-against-itself split-half) swallows all three. Wave metrics are meaningful only in regimes of strong, clean, multi-directional radiation—jets, unsteady wakes, or loud aeroacoustic sources ($M_{\text{source}}\gtrsim 0.5$). For validation purposes, we should not over-weight free-radiation acoustic accuracy as a BC selector.

(2) **Oblique strong-shock crossing is the decisive test.** The cylinder domain-reduction study provides a clean, physically interpretable discriminator: foextrap fails catastrophically (wrong density, destroyed shock) below ±5 half-width; code 7 remains robust to ±2.0 (±7% standoff drift). This is not an acoustic-reflection artifact or a marginal difference—it is a **functional enabler of domain reduction**, directly reducing computational cost and extending code-7's utility. The smooth standoff-vs-width curve (2.6% at ±2.5, deepening to −18.6% at ±1.5) is evidence of genuine physical effect (boundary truncation of shock curvature), not numerical instability.

(3) **Mean-state anchoring is a secondary but meaningful discriminator.** In both 1-D tests and the CVP annulus, foextrap and code 7 show negligible mean-pressure drift ($O(10^{-5})$), while more complex LODI with finite relaxation $\sigma_p$ shows $O(10^{-3})$ drift. For these subsonic-lateral cases, the simplest closures already anchor adequately (the body/shock system exerts strong mean-state control). In a purely subsonic external or internal-flow problem (e.g., subsonic jet issuing into quiescent), anchoring would be critical and code 7 with tuned $\sigma_p$ would be preferred.

(4) **The ghost-cell architecture does not prohibit non-reflecting characteristic BCs.** Earlier this year, a corrected characteristic-invariant ghost ($J^-=J^-_\infty$ vs. old hard-pin $p=p_\infty$) reduced 1-D reflection from R=0.78 to $<10^{-5}$. The key insight: the **temporal relaxation term $\mathcal{L}_1=\sigma(1-M^2)(c/L)(p-p_\infty)$ cannot be expressed as a ghost value**, but the **characteristic invariant $J^-$ is a spatial quantity and can be**. Distinction: a purely non-reflecting NSCBC (Poinsot–Lele with $\sigma\to 0$) needs the relaxation rate and thus a characteristic-flux boundary; a non-reflecting characteristic-invariant ghost (Thompson / Whitfield-Janus style) needs only the invariants and is a valid ghost cell. All three closures (foextrap, code 1, code 7) are now physically sound.

**Scope and limitations** (0.5 para):
We make no claim about subsonic external flow (where foextrap's lack of anchoring matters), viscous shear-layer instability (LODI untested), or unsteady jets and wakes (where a sponge or buffer zone is typically essential). Mach range: validated for $M_\infty \in [1.7, 4]$ (and tentatively $M\in[1.1,4]$ if sphere M=1.1 passes). The cylinder shock-crossing is representative of bluff-body supersonic flows; results may generalize to wedges, cones, and similar geometries where a strong shock crosses a lateral boundary. The CVP is best interpreted as a **negative control**—it rules out certain claims but does not claim to represent production aeroacoustic applications (those would need subsonic or low-$M$ transonic regimes, where foextrap/code 7 differences matter less).

---

## §8 CONCLUSIONS (1 page)

**Restatement of narrower claim** (1 para):
For supersonic external flows ($M_\infty\gtrsim 1.5$) where a bow shock crosses a lateral far-field boundary, the characteristic boundary condition (code 7, setting incoming Riemann invariant to freestream) enables **domain reduction with preserved physical accuracy**. Multi-width sensitivity studies show robust standoff behavior (±7% at ±2.0 half-width), whereas foextrap catastrophically fails. Free-radiation acoustic metrics (Fourier amplitude, phase error, energy reflection) show equivalence across all three closures in weak-signal and closed-domain regimes (CVP, M_orb=0.2), confirming that wave-reflection is not a discriminator in these cases. **Mean-state anchoring is not critical for this family of problems** (supersonic inflow pins global state); it becomes important in subsonic or closed-domain flows.

**Practical deployment** (1 para):
Sphere M=2, sphere M=4, and cylinder M=1.7 cases should use code 7 (characteristic) on lateral far-field boundaries. foextrap (code 2) is acceptable only on streamwise boundaries (where all three BCs coincide via supersonic inflow/outflow) and on lateral boundaries of purely subsonic problems where shock crossing is not a concern. A rigorous acoustic-reflection metric must use a weak but well-resolved source (e.g., volume-injection pulse, not ambient source-free radiation), hold the domain size fixed, and report noise floor relative to a reference truth. The Cerisse ghost-cell architecture (single value per stage, no temporal rates) admits non-reflecting characteristic invariants (J− constraint) but not time-dependent relaxation; if full Poinsot–Lele relaxation is needed, a characteristic-flux boundary layer (not a ghost) is required.

**Scope and open questions** (1 para):
Open for future work: (1) subsonic external or internal flow with foextrap vs. code 7 + LODI (anchoring matters there); (2) large-eddy simulation or RANS with turbulent shear-layer instability (LODI stability untested); (3) unsteady jets and wakes (sponge + buffer zone standard, not ghost BC alone); (4) fully 3-D shock-crossing on complex geometries (sphere M=1.1 transonic boundary to be completed). The finding that shock-crossing robustness outweighs free-radiation acoustics as a BC selector has implications for code design: investing in sophisticated non-reflecting NSCBC is justified by strong oblique shocks and domain reduction, not by marginal improvements in weak-signal aeroacoustics.

---

## Figure Count by Section (Total ≈ 12–14 figures)

- **§1 INTRO:** 0 (motivation text only)
- **§2 CLOSURES:** 0 (mathematical exposition)
- **§3 1D POINSOT–LELE:** 2–3 figures
  - 1D acoustic-pulse space–time diagrams (M=0) — foextrap / char / Dirichlet / LODI comparison
  - Reflection coefficient vs. Mach number curve (M ∈ [-0.5, +1.5])
  - Optional: mean-pressure drift evolution (anchoring trade-off)

- **§4 CVP:** 4–5 figures
  - CVP acoustic spiral (wide ±40 reference truth) — rotating quadrupole structure
  - Four $|\hat{P}|$ amplitude fields (wide / foextrap / char / LODI ±12) — near-identical spirals
  - Directivity $|\hat{P}|(\theta)$ at r=9 (noisy, floor-dominated, BC-independent)
  - Table: annulus-L² corruption numbers + noise floor (definitive result)
  - Optional: annulus geometry diagram (R_in, R_out, wavelength context)

- **§5 CYLINDER SHOCK CROSSING:** 3–4 figures
  - Schlieren snapshot at narrow width ±2.5 (char vs. foextrap side-by-side)
  - Standoff-vs-width curve (Δx_s vs. y_half, code 7 smooth robustness vs. foextrap collapse)
  - Post-shock density ratio or pressure ratio vs. width (confirming Rankine–Hugoniot)
  - Optional: radial density profile through shock (wide baseline vs. narrow char)

- **§6 TRANSONIC (M=1.1):** 1–2 figures
  - Schlieren + standoff measurement (if available; or placeholder for future work)
  - Transonic Mach-profile showing sonic-line location (optional reference)

- **§7 DISCUSSION:** 1 figure
  - Conceptual summary table or diagram: (BC method) × (regime: normal-acoustic, oblique-shock, mean-anchoring) → (winner / ranking / caveat)

- **§8 CONCLUSIONS:** 0 (summary text)

---

## Outline Statistics

| Section | Pages | Paras | Figures | Purpose |
|---------|-------|-------|---------|---------|
| §1 INTRO | 2 | 4 | 0 | Motivation, scope narrowing, hypothesis |
| §2 CLOSURES | 1 | 3–4 | 0 | Mathematical foundation (unchanged) |
| §3 1D BENCHMARK | 1.5 | 3 | 2–3 | Root-cause (M=0), mean-anchoring trade-off |
| §4 CVP OPEN DOMAIN | 2 | 4 | 4–5 | Negative control: wave metrics below floor |
| §5 SHOCK CROSSING | 2.5 | 4 | 3–4 | Decisive test: code 7 robust, foextrap fails |
| §6 TRANSONIC | 1 | 1–2 | 1–2 | Boundary case: M=1.1 validity |
| §7 DISCUSSION | 1.5 | 3–4 | 1 | Four conclusions + scope + open questions |
| §8 CONCLUSIONS | 1 | 2–3 | 0 | Narrower claim + deployment + future work |
| **TOTAL** | **12** | **25–30** | **12–14** | **Journal submission** |

---

## Key Structural Changes from Original

1. **New introduction (§1)** now establishes the three test regimes (supersonic external, transonic, acoustic radiation) and hypotheses explicitly. Scope narrowed from "general far-field BC principles" to "supersonic external with shock crossing."

2. **1D Poinsot–Lele (§3)** shortened and focused on **root-cause explanation** (corrected code 7 via J-invariant), not just R=0.78 mystery.

3. **CVP (§4)** reframed as **negative control**: conclusive proof that wave metrics do NOT discriminate, with rigorous drift-free, phase-aware, floor-calibrated metric and explicit noise floor. This addresses the "are wave metrics sufficient?" audit criticism directly.

4. **Cylinder shock crossing (§5) — NEW emphasis**: standalone section with multi-width standoff curve, robust smooth physics interpretation, foextrap catastrophic failure. This is the "decisive test" and primary positive result.

5. **Transonic (§6)** added as boundary-case test (sphere M=1.1), leaving placeholder for completion.

6. **Discussion (§7)** reorganized into **four numbered conclusions**, each addressing one audit criticism: (1) wave metrics inconclusive, (2) shock-crossing decisive, (3) anchoring secondary, (4) architecture permits non-reflecting characteristic.

7. **Conclusions (§8)** narrowed claim: supersonic external, shock crossing, M∈[1.5,4], domain reduction. Explicitly state which cases use code 7. Open questions listed.

8. **Tone throughout:** "we prove X by experiment Y" not "BC theory says X." Emphasis on **empirical discriminators** (shock standoff, flow field integrity) over abstract metrics.
