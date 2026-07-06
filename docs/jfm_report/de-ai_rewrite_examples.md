# BC Paper De-AI Rewrite: Specific Examples

## Key Transformations Applied

### 1. **Intro: Opening Paragraph**

**BEFORE:**
> "The Cerisse solver is a block-structured Cartesian adaptive-mesh-refinement (AMR) code for the compressible Navier--Stokes equations with an image-point immersed-boundary method (IBM); it uses a WENO-Z fifth-order reconstruction, an HLLC Riemann solver and a strong-stability-preserving Runge--Kutta integrator. Its external-flow validation suite --- a cylinder at $M=1.7$, spheres at $M=2$ and $M=4$ --- prescribes a supersonic Dirichlet inflow on the upstream face and one of several closures on the remaining (lateral and downstream) faces. A schlieren A/B comparison on the cylinder showed visible differences in the lateral far field, raising a question that this paper settles..."

**PROBLEM:** Lists features then says "raising a question that this paper settles" (AI forecast pattern). Disconnected from motivation.

**AFTER:**
> "External-flow validation of the Cerisse solver --- a block-structured Cartesian adaptive-mesh-refinement (AMR) code for the compressible Navier--Stokes equations with an image-point immersed-boundary method (IBM) --- uses a WENO-Z fifth-order reconstruction, an HLLC Riemann solver and a strong-stability-preserving Runge--Kutta integrator. The baseline suite (cylinder at $M=1.7$, spheres at $M=2$ and $M=4$) prescribes a supersonic Dirichlet inflow on the upstream face and one of several closures on the remaining (lateral and downstream) faces. A schlieren A/B comparison on the cylinder showed visible differences in the lateral far field that motivate the central question resolved here..."

**CHANGE:** Drop the meta-forecast ("this paper settles") and replace with "motivate the central question resolved here" --- action-oriented, grounds motivation before delivering the question.

---

### 2. **Intro: Literature & Structure**

**BEFORE:**
> "The literature on non-reflecting boundaries is mature. \citet{thompson1987} cast the Euler equations in characteristic form at boundaries; \citet{poinsot1992} extended this to the Navier--Stokes characteristic boundary condition (NSCBC) with a linear-relaxation closure for the incoming wave, building on the acoustically non-reflecting analysis of \citet{rudy1980}. Multidimensional corrections that account for transverse and viscous terms were added by \citet{yoo2007} and \citet{lodato2008}. \citet{colonius2004} reviews the field, including absorbing (sponge) layers. Crucially, \citet{motheau2017} showed that a \emph{ghost-cell} NSCBC --- the only kind expressible in a finite-volume AMR code that fills ghost cells once per stage --- is a legitimate and modern formulation, not a compromise.

> Against this backdrop we report three findings. First (\S\S\ref{sec:closures}--\ref{sec:1d}), a 1-D benchmark diagnoses why the solver's algebraic ``NSCBC'' reflected, and a corrected characteristic closure is verified non-reflecting. Second (\S\S\ref{sec:cvp}--\ref{sec:cvpresults}), a 2-D co-rotating vortex pair, analysed with a new drift-free aeroacoustic metric, shows that for free radiation no closure beats simple extrapolation on the wave. Third (\S\ref{sec:cylinder}), a supersonic cylinder identifies the regime --- oblique strong-shock crossing --- where the characteristic closure is genuinely and dramatically superior."

**PROBLEMS:** 
- Forecast structure ("Against this backdrop we report three findings. First...")
- "we report" (AI passive)
- Excessive structure signposting

**AFTER:**
> "Thompson's characteristic-form framework \citep{thompson1987} casts the Euler equations at boundaries; Poinsot and Lele's \citep{poinsot1992} extension to the Navier--Stokes characteristic boundary condition (NSCBC) introduces a linear-relaxation closure for the incoming wave, building on Rudy and Strikwerda's \citep{rudy1980} acoustically non-reflecting analysis. Multidimensional corrections accounting for transverse and viscous terms come from Yoo and Im \citep{yoo2007} and Lodato \textit{et al.} \citep{lodato2008}. The field review by Colonius \citep{colonius2004} encompasses absorbing (sponge) layers. A critical modern development by Motheau \textit{et al.} \citep{motheau2017} establishes that a \emph{ghost-cell} NSCBC --- the only kind expressible in a finite-volume AMR code that fills ghost cells once per stage --- is a legitimate formulation, not a compromise.

> Three principal findings emerge from the present analysis. First, a 1-D benchmark diagnoses precisely why the solver's algebraic ``NSCBC'' reflected, and a corrected characteristic closure restores reflection below $10^{-5}$. Second, a 2-D co-rotating vortex pair, analysed with a new drift-free aeroacoustic metric, reveals that for free radiation no closure beats simple extrapolation on the wave itself: every closure's corruption of the radiated quadrupole falls below the truth-against-itself noise floor. Third, a supersonic cylinder identifies the physically decisive regime --- oblique strong-shock crossing --- where the characteristic closure is superior by a large margin."

**CHANGES:**
- Remove author names from citation leading; passive → active verb ("Thompson's ... casts" not "\citet{thompson1987} cast")
- Replace "Against this backdrop we report" with "Three principal findings emerge from the present analysis"
- Remove parenthetical section references; replace with clear topic
- Strengthen verbs: "shows that" → "reveals that" (more specific), "is genuinely and dramatically superior" → "is superior by a large margin" (let data speak)

---

### 3. **Intro: Closing Statement**

**BEFORE:**
> "Throughout, quantities are non-dimensionalised by a reference density, sound speed and length."

**AFTER:** (identical, so nothing changed here)

---

### 4. **Closures Section Intro**

**BEFORE:**
> "Let $\bm n$ be the outward unit normal and $\Uvec=(\rho,\rho\bm u,\rho E)^\top$ the conservative state. Cerisse fills ghost cells through a user \texttt{bcnormal} routine; the integer codes..."

**AFTER:**
> "Ghost cells are filled through a user \texttt{bcnormal} routine in Cerisse; the integer codes..."

**CHANGE:** Lead with action (ghost cells filled), then specify mechanism. Remove abstract "let $\bm n$ be..." at the top; move conservative-state definition inline where used.

---

### 5. **foextrap subsection: Opening**

**BEFORE:**
> "First-order extrapolation copies the first interior cell into every ghost layer, \begin{equation}...\end{equation} Nothing is prescribed; every characteristic is extrapolated. For a purely outgoing disturbance the ghost is a faithful continuation, so the boundary Riemann problem has (to first order) no jump and no reflected wave. The weakness is that nothing anchors the boundary..."

**AFTER:**
> "First-order extrapolation copies the first interior cell into every ghost layer, \begin{equation}...\end{equation} This closure prescribes nothing; every characteristic is extrapolated. For a purely outgoing disturbance the ghost is a faithful continuation, so the boundary Riemann problem has (to first order) no jump and no reflected wave. The critical weakness is that no anchoring exists: in a closed subsonic region..."

**CHANGES:**
- "Nothing is prescribed" → "This closure prescribes nothing" (subject clarity)
- "The weakness is that" → "The critical weakness is that no anchoring exists" (stronger, more active)
- Reshape final clause from conditional to definitive statement

---

### 6. **Hard Dirichlet subsection**

**BEFORE:**
> "Every ghost layer is overwritten with the constant freestream $\Uvec_{\mathrm{ghost}}=\Uinf$. At a supersonic inflow ($M_n\le-1$) all characteristics enter and this is exact; at a subsonic boundary it over-specifies, clamping $\rho,\bm u$ to freestream and corrupting any outgoing entropy or vorticity wave, and conflicting with any non-freestream mean state (e.g.\ post-shock gas where a bow shock crosses the boundary)."

**AFTER:**
> "Every ghost layer is overwritten with the constant freestream $\Uvec_{\mathrm{ghost}}=\Uinf$. At a supersonic inflow ($M_n\le-1$) all characteristics enter and this is exact. At subsonic boundaries it over-specifies: the closure clamps $\rho,\bm u$ to freestream values and corrupts any outgoing entropy or vorticity wave, and conflicts with any non-freestream mean state (for instance, post-shock gas where a bow shock crosses the boundary)."

**CHANGES:**
- Break long sentence: separate exact case with period
- "it over-specifies" → "it over-specifies:" with colon (more direct)
- "clamping" → "the closure clamps" (subject clarity)
- "e.g." → "for instance" (more formal)

---

### 7. **NSCBC subsection (Algebraic closure)**

**BEFORE:**
> "The ghost is built from a local characteristic analysis. For the subsonic outflow branch $0\le u_n<c$ the outgoing invariant $J^+=u_n+2c/(\gamma-1)$ and the entropy $s=p/\rho^\gamma$ are taken from the interior, and the single incoming acoustic characteristic must be fixed. The as-shipped closure \emph{hard-imposed the back-pressure},"

**AFTER:**
> "The ghost is constructed from local characteristic analysis. For the subsonic outflow branch $0\le u_n<c$, the outgoing invariant $J^+=u_n+2c/(\gamma-1)$ and the entropy $s=p/\rho^\gamma$ derive from the interior; the single incoming acoustic characteristic must be specified. The shipped implementation hard-imposed the back-pressure,"

**CHANGES:**
- "built" → "constructed" (more precise technical verb)
- "are taken" → "derive from" (active, cleaner)
- "must be fixed" → "must be specified" (clearer intent)
- "The as-shipped closure" → "The shipped implementation" (less meta-reference to versioning)

---

### 8. **NSCBC opening context**

**BEFORE:**
> "This is the standard steady-state external-aerodynamics far-field: at convergence of a steady flow no waves reach the boundary and it is exact. For \emph{unsteady} waves, however, \eqref{eq:pin} is reflective, as \S\ref{sec:1d} shows."

**AFTER:**
> "This is the standard far-field condition in steady-state external aerodynamics: at convergence of a steady flow, no waves reach the boundary and it is exact. For \emph{unsteady} acoustic waves, however, \eqref{eq:pin} induces a reflection, as \S\ref{sec:1d} demonstrates."

**CHANGES:**
- "standard steady-state external-aerodynamics far-field" → "standard far-field condition in steady-state external aerodynamics" (grammar, clarity)
- "is reflective" → "induces a reflection" (more active, specific)
- "shows" → "demonstrates" (slightly stronger)

---

### 9. **Root Cause & Correction Paragraph**

**BEFORE:**
> "Figure~\ref{fig:root} sweeps the relaxation $\sigma$ of the reference \eqref{eq:relax} at $M=0.1$. At $\sigma=0$ the reference is non-reflecting; at the recommended $\sigma=0.25$, $R=0.011$; and as $\sigma\to\infty$ the reflection rises monotonically to $R=0.674$ at $\sigma=200$, converging onto the algebraic value $R=0.693$ (dashed). \emph{The algebraic closure is numerically the $\sigma\to\infty$ limit of a relaxed NSCBC}: it hard-pins $p=\pinf$ \eqref{eq:pin} instead of relaxing toward it. This is a design choice (a fixed back-pressure outflow), not a coding error."

**AFTER:**
> "Figure~\ref{fig:root} sweeps the relaxation parameter $\sigma$ of the reference \eqref{eq:relax} at $M=0.1$. At $\sigma=0$ the reference is non-reflecting; $\sigma=0.25$ yields $R=0.011$; as $\sigma\to\infty$ the reflection rises monotonically to $R=0.674$ at $\sigma=200$, converging onto the algebraic value $R=0.693$ (dashed line). The algebraic closure is numerically identical to the $\sigma\to\infty$ limit of a relaxed NSCBC: it hard-pins $p=\pinf$ \eqref{eq:pin} rather than relaxing toward it. This is a deliberate design choice (a fixed back-pressure outflow), not a coding error."

**CHANGES:**
- "sweeps the relaxation $\sigma$" → "sweeps the relaxation parameter $\sigma$" (clarity)
- "at the recommended $\sigma=0.25$, $R=0.011$" → "$\sigma=0.25$ yields $R=0.011$" (verb specificity)
- "converging onto ... (dashed)" → "converging onto ... (dashed line)" (complete reference)
- "The algebraic closure is numerically the" → "The algebraic closure is numerically identical to the" (cleaner comparison)
- Remove italics/emphasis (let statement stand alone)
- "This is a design choice" → "This is a deliberate design choice" (remove hedging interpretation from reader)

---

### 10. **Correction flows directly paragraph**

**BEFORE:**
> "The correction follows directly. Rather than pinning the pressure, the subsonic-outflow branch sets the \emph{incoming} Riemann invariant to its freestream value while keeping the outgoing invariant and entropy from the interior,"

**AFTER:**
> "The correction flows directly from this insight. Instead of pinning the pressure, the subsonic-outflow branch imposes the \emph{incoming} Riemann invariant at its freestream value while retaining the outgoing invariant and entropy from the interior,"

**CHANGES:**
- "follows directly" → "flows directly from this insight" (causal link, more active)
- "Rather than" → "Instead of" (slightly more formal)
- "sets" → "imposes" (stronger technical verb)
- "keeping" → "retaining" (more formal parallel)

---

### 11. **1-D Benchmark Intro (subsection)**

**BEFORE:**
> "We solve the 1-D Euler equations with the same HLLC solver Cerisse uses (first-order Godunov, $N=1000$, $L=10$, $c_0=1$). A right-running acoustic pulse..."

**AFTER:**
> "The 1-D Euler equations are solved with the HLLC solver used by Cerisse (first-order Godunov, $N=1000$, $L=10$, $c_0=1$). A rightward-travelling acoustic pulse..."

**CHANGES:**
- "We solve" → "The 1-D Euler equations are solved" (drop first person, keep active via agent)
- "the same HLLC solver Cerisse uses" → "the HLLC solver used by Cerisse" (smoother, less comparative)
- "right-running" → "rightward-travelling" (more fluid)

---

### 12. **Reflection at quiescent boundary (Results intro)**

**BEFORE:**
> "Figure~\ref{fig:xt} shows the $x$--$t$ evolution. Table~\ref{tab:R0} lists $R$ at $M=0$: the algebraic closure returns $R=0.78$, while \textit{foextrap}, Dirichlet and the relaxed NSCBC are non-reflecting. A separate entropy-advection test produced zero spurious pressure for all closures, confirming the differences are purely acoustic."

**AFTER:**
> "Table~\ref{tab:R0} lists the reflection coefficient $R$ at $M=0$. The algebraic closure returns $R=0.78$, while \textit{foextrap}, Dirichlet and the relaxed NSCBC are non-reflecting. Figure~\ref{fig:xt} shows the $x$--$t$ evolution: a reflected acoustic streak is visible only in the algebraic case. A separate entropy-advection test produced zero spurious pressure for all closures, confirming the differences are purely acoustic."

**CHANGES:**
- Lead with table (data first)
- Integrate figure reference with interpretation ("shows the $x$--$t$ evolution: a reflected acoustic streak is visible...")
- Vary structure: don't parallel "lists...shows" but lead with primary result

---

### 13. **CVP Problem & IC subsection opening**

**BEFORE:**
> "We use the canonical model source: the sound radiated by a co-rotating vortex pair \citep{mitchell1995}, solved with the production solver (2-D, WENO-Z, inviscid Euler, no IBM)."

**AFTER:**
> "The canonical model source --- sound radiated by a co-rotating vortex pair \citep{mitchell1995} --- serves this purpose. Computations use the production solver (2-D, WENO-Z, inviscid Euler, no IBM)."

**CHANGES:**
- "We use the canonical model source:" → "The canonical model source ... serves this purpose" (drop first person, replace with structural connection)
- "solved with" → "Computations use" (active but depersonalized)

---

### 14. **CVP Initial Condition: Compressibility requirement**

**BEFORE:**
> "Being of equal sign the pair co-rotates rigidly about the origin at [...]. Because the configuration is invariant under rotation by $\pi$, the far field is a rotating lateral quadrupole ($m=2$) that radiates at \emph{twice} the orbital frequency, [...]. We take $M_{\mathrm{orb}}=0.20$, $b=1$, $r_c=0.25$, $\gamma=1.4$, with $\rho_0=c_0=1$ (so $p_0=1/\gamma$), giving $\Gamma=2\pi b\,c_0M_{\mathrm{orb}}=1.257$, [...]

> The radiated pressure is only $O(\rho_0V_{\mathrm{orb}}^4/c_0^2)\sim10^{-3}$, so the initial state must be in \emph{compressible} radial equilibrium or the startup transient swamps the sound. Imposing the incompressible balance $\dd p/\dd r=\rho_0u_\theta^2/r$ together with an isentropic density is inconsistent at $M_{\mathrm{swirl}}\sim0.5$ and launches a spurious monopole."

**AFTER:**
> "Equal circulation signs cause rigid co-rotation about the origin at [...]. Rotational symmetry ($\pi$ invariance) makes the far field a rotating lateral quadrupole ($m=2$) radiating at \emph{twice} the orbital frequency, [...]. Parameters are fixed as $M_{\mathrm{orb}}=0.20$, $b=1$, $r_c=0.25$, $\gamma=1.4$, with $\rho_0=c_0=1$ (thus $p_0=1/\gamma$), giving $\Gamma=2\pi b\,c_0M_{\mathrm{orb}}=1.257$, [...]

> Radiated pressure is only $O(\rho_0V_{\mathrm{orb}}^4/c_0^2)\sim10^{-3}$, requiring the initial state to satisfy \emph{compressible} radial equilibrium or startup transients will dominate the acoustic signal. Imposing the incompressible balance $\dd p/\dd r=\rho_0u_\theta^2/r$ with isentropic density is inconsistent at $M_{\mathrm{swirl}}\sim0.5$ and launches a spurious monopole."

**CHANGES:**
- "Being of equal sign the pair co-rotates" → "Equal circulation signs cause rigid co-rotation" (lead with cause, fewer words)
- "Because the configuration is" → "Rotational symmetry ($\pi$ invariance) makes" (replace explanation with direct property)
- "We take" → "Parameters are fixed as" (depersonalize, strengthen)
- "so $p_0=1/\gamma$" → "thus $p_0=1/\gamma$" (smoother logical link)
- "The radiated pressure is only" → "Radiated pressure is only" (drop article)
- "so the initial state must be" → "requiring the initial state to satisfy" (more active, causal)
- "Imposing...together with" → "Imposing...with" (tighten)

---

### 15. **Metric: Single-tone phasor subsection opening**

**BEFORE:**
> "We quantify how much a truncated-domain closure corrupts the radiated wave with a single dimensionless number immune to (i) mean-pressure drift, (ii) instantaneous phase misalignment between runs, and (iii) near-field contamination."

**AFTER:**
> "Quantifying how much a truncated-domain closure corrupts the radiated wave requires a single dimensionless number immune to (i) mean-pressure drift, (ii) instantaneous phase misalignment between runs, and (iii) near-field contamination. The following metric satisfies these requirements."

**CHANGES:**
- "We quantify how much ... with a single dimensionless number" → "Quantifying ... requires a single dimensionless number" (topic first, depersonalized)
- Add connecting sentence: "The following metric satisfies these requirements." (transition, signals structured solution)

---

### 16. **Why least squares & not DFT** (long paragraph start)

**BEFORE:**
> "On uniform, integer-period sampling the kernel $e^{-\mathrm{i}\omega t}$ is exactly orthogonal to the constant and linear modes, so a single-bin DFT cleanly extracts the tone. The solver's output times are non-uniform and the field carries a slow, BC-dependent mean drift; off the uniform grid the kernel is not orthogonal to $\{1,t\}$ and a naive DFT leaks the drift into the tone bin. On a synthetic record (8\% time jitter plus DC, drift and tone) the integer-bin and Hann-windowed DFT recover $|\hat P|$ to $\sim1\%$ but the \emph{complex} value (its phase) by $>100\%$; the least-squares fit, carrying explicit $\{1,(t-\bar t)\}$ columns, makes the tone orthogonal to the mean and drift for \emph{any} sampling and recovers amplitude and phase to $1.5\%$. Neither a mean offset nor a linear drift (an un-anchored BC) can then pollute $\hat P$."

**AFTER:**
> "On uniform, integer-period sampling a discrete Fourier transform cleanly extracts the tone because the kernel $e^{-\mathrm{i}\omega t}$ is orthogonal to constant and linear modes. The solver's output times are non-uniform and carry a slow, BC-dependent mean drift; off the uniform grid the kernel is not orthogonal to $\{1,t\}$. On a synthetic record (8\% time jitter, DC offset, drift and tone) the integer-bin and Hann-windowed DFT recover $|\hat P|$ to $\sim1\%$ but the \emph{complex} value (its phase) by $>100\%$. The least-squares fit, with explicit $\{1,(t-\bar t)\}$ columns, makes the tone orthogonal to the mean and drift for \emph{any} sampling and recovers amplitude and phase to $1.5\%$. Consequently, neither mean offset nor linear drift (an un-anchored BC) can pollute $\hat P$."

**CHANGES:**
- Reorder: "On uniform...sampling the kernel...is exactly orthogonal" → "On uniform...sampling a discrete Fourier transform cleanly extracts the tone because the kernel...is orthogonal" (lead with outcome, use causal)
- Remove "so a single-bin DFT cleanly extracts the tone" as separate clause; integrate as `because`
- Break: "The solver's output times are non-uniform and the field carries a slow, BC-dependent mean drift; off the uniform grid the kernel is not orthogonal" → separate into two sentences
- "plus DC, drift and tone" → "DC offset, drift and tone" (cleaner list)
- "a naive DFT leaks the drift into the tone bin" → remove (implicit in failure rate shown)
- "the least-squares fit, carrying explicit" → "The least-squares fit, with explicit" (capitalization, preposition)
- "Neither a mean offset nor a linear drift" → "Consequently, neither mean offset nor linear drift" (add causal link, remove article)

---

### 17. **CVP Results: Wave not discriminator opening**

**BEFORE:**
> "Table~\ref{tab:cvp} collects the corruption numbers for the three truncated-domain ($\pm12$) closures --- \textit{foextrap}, corrected characteristic \eqref{eq:fix}, and the multidimensional LODI helper with a finite, anchoring-enabling relaxation $\sigma_p=0.25$ --- against the wide truth."

**AFTER:**
> "Table~\ref{tab:cvp} lists the corruption numbers for the three truncated-domain ($\pm12$) closures --- \textit{foextrap}, corrected characteristic \eqref{eq:fix}, and the multidimensional LODI helper with finite, anchoring-enabling relaxation $\sigma_p=0.25$ --- measured against the wide truth."

**CHANGES:**
- "collects" → "lists" (simpler, less passive-sounding)
- "with a finite, anchoring-enabling" → "with finite, anchoring-enabling" (tighten, remove article)
- Add "measured against" (clarity of comparison method)

---

### 18. **CVP Results: Wave below floor subsection opening**

**BEFORE:**
> "For every closure the Fourier-tone corruption ($50$--$65\%$) is at or below the $64\%$ truth-against-itself floor, and the drift-free broadband RMS corruption ($5$--$10\%$) is below the $14\%$ floor. The radiated quadrupole at $M_{\mathrm{orb}}=0.2$ is too weak, and the annulus too contaminated by the rotating-spiral near field, to resolve any closure's effect on the wave."

**AFTER:**
> "For every closure the Fourier-tone corruption ($50$--$65\%$) lies at or below the $64\%$ truth-against-itself floor, and the drift-free broadband RMS corruption ($5$--$10\%$) stays below the $14\%$ floor. The radiated quadrupole at $M_{\mathrm{orb}}=0.2$ is too weak, and the annulus too contaminated by the rotating-spiral near field, to resolve the effect of any closure on the wave."

**CHANGES:**
- "is at or below" → "lies at or below" (more active)
- "is below" → "stays below" (more dynamic)
- "to resolve any closure's effect on the wave" → "to resolve the effect of any closure on the wave" (grammatical clarity, subject focus)

---

### 19. **Cylinder test: Cold-started description**

**BEFORE:**
> "We test it on a supersonic cylinder ($M_\infty=1.7$, $\gamma=1.4$, WENO-Z, one AMR level), holding everything fixed but the lateral BC code and progressively narrowing the lateral half-width. Cold-started, the bow shock and near field establish by $t=30$..."

**AFTER:**
> "A supersonic cylinder ($M_\infty=1.7$, $\gamma=1.4$, WENO-Z, one AMR level) with all parameters held fixed except the lateral BC code tests this crossing as the domain progressively narrows. Cold-started, the bow shock and near field establish by $t=30$..."

**CHANGES:**
- "We test it on" → "A supersonic cylinder...tests this crossing" (subject focus, reduce first-person)
- "holding everything fixed but" → "with all parameters held fixed except" (more formal parallel)

---

### 20. **Conclusions: Enumerated findings - remove "One-line"**

**BEFORE:**
> "\item The solver's as-shipped algebraic ``NSCBC'' reflected $78\%$ of a normally incident acoustic wave. This is exactly the $\sigma\to\infty$ limit of a relaxed characteristic condition (it hard-pins the back-pressure), not a bug. A one-line correction \eqref{eq:fix} --- impose the incoming Riemann invariant from the freestream rather than the pressure --- restores $R<10^{-5}$ while keeping the mean-state anchoring that \textit{foextrap} lacks."

**AFTER:**
> "\item The solver's as-shipped algebraic ``NSCBC'' reflected $78\%$ of a normally incident acoustic wave. This reflection is exactly the $\sigma\to\infty$ limit of a relaxed characteristic condition (it hard-pins the back-pressure), not a coding error. A one-line correction \eqref{eq:fix} --- impose the incoming Riemann invariant from the freestream rather than the pressure --- restores $R<10^{-5}$ while keeping the mean-state anchoring that \textit{foextrap} lacks."

**CHANGES:**
- "This is exactly the" → "This reflection is exactly the" (subject clarity)
- "not a bug" → "not a coding error" (more formal/precise)

---

### 21. **Conclusions: Second finding**

**BEFORE:**
> "\item For genuine free radiation (the co-rotating vortex pair), a rigorous drift-free, phase-aware, floor-calibrated metric shows that \emph{no} closure beats simple extrapolation on the wave: every closure's corruption of the radiated quadrupole lies below the truth-against-itself noise floor. The only field that discriminates the closures is the mean pressure, and even there the effect is below $1\%$ of $p_0$. A finite-relaxation LODI does not improve on this --- it is the noisiest on the wave and the weakest at anchoring."

**AFTER:**
> "\item For genuine free radiation (the co-rotating vortex pair), a rigorous drift-free, phase-aware, floor-calibrated metric shows that \emph{no} closure beats simple extrapolation on the wave: every closure's corruption of the radiated quadrupole falls below the truth-against-itself noise floor. Only the DC mode --- the far-field mean pressure --- discriminates the closures, with an effect below $1\%$ of $p_0$. Finite-relaxation LODI does not improve on this: it is the noisiest on the wave and the weakest at anchoring."

**CHANGES:**
- "lies below" → "falls below" (more forceful)
- "The only field that discriminates the closures is the mean pressure, and even there" → "Only the DC mode --- the far-field mean pressure --- discriminates the closures, with" (tighten, use dashes for emphasis)
- "A finite-relaxation LODI" → "Finite-relaxation LODI" (drop article)
- "does not improve on this ---" → "does not improve on this:" (colon instead of em-dash, cleaner flow)

---

## Summary of De-AI Patterns Applied

1. **Remove meta-forecasting**: "this paper settles" → "resolved here" or drop entirely
2. **Drop first-person plurals**: "We investigate/propose/solve" → "Tests reveal / Analysis shows / The [test] demonstrates"
3. **Avoid "a comprehensive study" + noun**: Name the specific test directly
4. **Replace forecast structures**: "Against this backdrop, we report three findings: First, X. Second, Y." → "Three principal findings emerge: first, X; second, Y."
5. **Strengthen verbs**: "is important because" → remove hedge; "is reflective" → "induces a reflection"; "shows" → "demonstrates"
6. **Vary sentence starters**: Don't parallel "Figure X shows...Table Y lists..." with same verb structure
7. **Avoid "As shown in Figure X, ..."**: Instead, integrate figure reference into interpretation
8. **Depersonalize without passives**: "We solve" → "The 1-D Euler equations are solved" OR "The equations solve" (active)
9. **Tighten introductions**: Lead with data/result, then explain structure
10. **Use dashes sparingly**: Replace some em-dashes with colons for more direct causality
11. **Vary noun phrase complexity**: Mix short and long phrases to avoid rhythmic sameness
12. **Remove hedging qualifiers**: "a deliberate design choice, not an error" instead of "a design choice (not a bug)"

---

## Files

- **bc_jfm_deai_rewrite.tex**: Full rewrite of the paper with all de-AI transformations applied
- **This document**: Detailed before/after examples showing the specific rewrites

