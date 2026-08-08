# Screech Reproduction: Campaign Audit and Next-Case Requirements

Date: 2026-07-17

## Verdict

The campaign has resolved shock topology, mean shock-cell trends and several
unsteady near-field components. It has not yet reproduced screech under a
strict definition. A narrow spectral peak or an upstream phase slope is not
enough. A successful result must demonstrate one self-sustained global mode at
the experimental operating condition, with the correct azimuthal mode,
frequency, downstream instability wave, upstream return wave, source region,
far/near acoustic signature and saturated amplitude.

The present failures are not primarily caused by SSP-RK3, CFL=0.3 or an
insufficiently small integration time step. They are more plausibly caused by
insufficient loop gain and an incorrect loop phase in the nozzle/lip model,
with shear-layer resolution and acoustic propagation as the next controls.

Two objectives should be kept separate:

1. **Solver capability:** demonstrate that CERISSE can sustain a published
   screeching flow under a completely specified numerical configuration.
2. **Panda reproduction:** match the Panda nozzle-block experiment, frequency
   and mode.

The first objective should be passed before more expensive attempts at the
second are interpreted physically.

## The governing test

Screech is a global feedback oscillator. A useful bookkeeping model is

\[
G(f,m)=G_{KH}\,G_{shock}\,G_{up}\,R_{lip},
\]

where the factors represent downstream Kelvin-Helmholtz amplification,
shock/shear-layer sound generation, upstream propagation and lip/nozzle
receptivity. Self-excitation requires both

\[
|G(f,m)|\ge 1,
\qquad
\arg G(f,m)=2\pi n.
\]

Accurate first-cell geometry only constrains part of `G_shock`. It does not
show that the other three factors are present or that the loop phase closes.

## What the existing cases actually establish

| Case family | What it established | Why it is not a screech reproduction |
|---|---|---|
| Euler/top-hat and NPR sweeps | Shock cells and pressure-mismatch topology | No controlled viscous exit layer or defensible KH growth; discrete tones can be numerical |
| Old 3D NPR=3 `case2_powell_tam` | A candidate 19.7 kHz near-field upstream component | No matching experimental benchmark, global mode or saturated acoustic validation; the published internal coherence estimate used about 9.95 kHz Welch bins |
| Mj=1.42 RZ cases | Mean first-cell and shock motion | Panda's 5.4 kHz mode is helical (`m=1`) and cannot exist in RZ |
| Mj=1.42 3D pilots | Three-dimensional shock-cell evolution and modal capability | Simplified exit/flange geometry, marginal lip-layer resolution and short exploratory records preclude a negative or positive screech claim |
| Mj=1.19 old L6 cases | Fine near-lip resolution and axisymmetric dynamics | The former inlet/boundary model did not reproduce the experimental nozzle receptivity; no persistent 8.4 kHz global line appeared |
| Mj=1.19 L4 NSCBC case | Correct target mode, improved characteristic annulus and rigid flush flange | It contains a local 9.83 kHz component, but no persistent far-field/global lock, no complete phase branch and no kinematic closure |

The dedicated audit of the last row is reproducible with:

```bash
python3 analysis/panda_benchmark/m119_screech_audit.py
```

## Ranked campaign blockers

### 1. The nozzle-block geometry and loop phase are not the experiment

Panda used a 25.4 mm convergent nozzle inside a three-dimensional metal nozzle
block. The main flange was 305 mm (`12D`) in diameter and approximately
202 mm (`7.95D`) upstream of the exit. The exit body had an outer diameter of
38 mm (`1.5D`). The paper reports that covering the flange changed both
screech amplitude and frequency.

The current model instead places an ideal `12D` rigid plate in the exit plane.
At 8.4 kHz the ambient wavelength is about 41.2 mm (`1.62D`), so moving the
reflector by 202 mm changes the path by about 4.9 acoustic wavelengths before
accounting for the conical body and lip scattering. This is not a small
geometric approximation for a phase-closure problem.

The prescribed sonic patch also omits the internal nozzle and the physical lip
thickness. Its LODI annulus may be less reflective than the old Dirichlet patch,
but its jet-on receptivity has not been measured. A rigid-wall pulse test proves
the flange reflects sound; it does not prove that a returning wave excites a new
KH wavepacket at the aperture.

### 2. The Mj=1.19 L4 shear layer is marginal for the present flux

The L4 case has

\[
\Delta_{min}=D/256=99.2\ \mu\mathrm{m},
\qquad \delta=254\ \mu\mathrm{m}.
\]

For the prescribed tanh profile, the vorticity thickness is approximately
`2 delta`, or 5.12 L4 cells, while the momentum thickness is approximately
`delta/2`, or only 1.28 cells. A low-dissipation spectral-like solver may still
amplify a layer at this resolution; WENO-Z5 with local LLF cannot be assumed to
do so. A suppressed KH growth rate directly reduces `G_KH` below unity.

The finest grid must be chosen by a measured linear growth-rate error, not by
the ability to draw a sharp shock. A practical starting point is at least
10-12 cells across vorticity thickness in the lip and source corridors. For the
current profile this means at least D/512, with D/1024 retained as a resolution
audit.

### 3. The weak acoustic path has not been verified end to end

Gradient tagging follows shocks and vortices but not small-amplitude sound.
The source-to-lip/reflection path therefore needs a fixed acoustic corridor and
a frozen production BoxArray. Refluxing should remain enabled. Before a long
run, an 8.4 kHz oblique wave packet must traverse the actual AMR interfaces and
return through the nozzle-block geometry with measured complex transmission:
amplitude, phase and direction.

The existing wall-pulse tests are useful but incomplete. They do not include
mean jet flow, the sonic aperture, the source region or the complete return
path.

### 4. Several cases excluded the experimental azimuthal mode

The Mj=1.19 Panda condition is the correct first target because its measured
8.4 kHz mode is axisymmetric. A two-dimensional RZ calculation can host it.
The Mj=1.42, 5.4 kHz Panda condition is helical and requires full 3D. No amount
of RZ resolution or run time can reproduce that mode.

### 5. The diagnostics have repeatedly promoted candidates too early

The latest 6 ms record illustrates the problem. The original high-resolution
periodogram selected six near-field probes at 9.83 kHz, but a consistent Welch
estimate gives much weaker prominence in the outer probes. The apparent phase
speed fits only through about `x/D=1`; it changes from -241 to -539 m/s as two
more lip-line probes are added. The 9.83 kHz cross-spectral mode ranks seventh,
not first, in the 1-30 kHz band. The experimental 8.4 kHz line is absent.

This feature should be named a **local coherent near-field component of
unresolved origin**, not screech.

## Recommended next calculation

### Stage 0A: reproduce one published numerical screech benchmark

For the shortest unambiguous solver test, reproduce the Berland--Bogey--Bailly
planar LES before changing more Panda parameters. Its public specification
includes `Mj=1.55`, `Re_h=6e4`, `p_e/p_inf=2.09`, a lip thickness of `h/4`, a
short explicit nozzle, characteristic inflow, non-reflecting outer boundaries,
a downstream sponge, about seven points through the imposed boundary layer and
at least 100 screech periods. A static mesh removes AMR from the diagnosis.

If this benchmark screeches but Panda does not, the leading uncertainty is the
Panda geometry, mode and receptivity model. If it also fails, inspect CERISSE
acoustic dissipation, KH amplification and boundary transfer before any new
production campaign.

### Stage 0B: three inexpensive solver gates

Do not start another production run until all three pass on the intended mesh
and numerical scheme.

1. **Acoustic gate:** propagate 8.4 kHz and its first harmonic through the
   complete fixed AMR path for at least 10D. Require less than 5% amplitude loss
   and less than 10 degrees accumulated phase error, excluding the designed
   nozzle reflection.
2. **KH gate:** calculate the spatial growth of the prescribed compressible
   shear layer and compare with LSA or a uniformly refined reference. Require
   less than 10% error in the most amplified growth rate and frequency.
3. **Receptivity gate:** impose a weak upstream-returning wave with mean jet
   flow present, measure the induced downstream KH wave, and report the complex
   transfer coefficient at 8.4 kHz. Repeat with the source removed to quantify
   numerical/background excitation.

### Stage 1: Panda Mj=1.19 exact-geometry RZ case

This is the shortest credible route to the first successful screech result.

- Reproduce the nozzle block from Panda figure 1: the recessed `12D` flange,
  conical exterior, `1.5D` outer exit diameter and a finite internal convergent
  nozzle section.
- Prescribe stagnation conditions upstream of the nozzle through a
  characteristic boundary. Do not prescribe the complete sonic exit state in
  the aperture plane.
- Use adiabatic no-slip internal walls, or a documented wall model if the
  internal layer cannot be wall resolved. Treat exit momentum thickness as a
  controlled sensitivity parameter.
- Use Navier-Stokes, WENO/TENO and SSP-RK3 at CFL=0.3. These are not the current
  leading blockers; retain them while isolating geometry and resolution.
- Keep D/512 or finer through the lip layer and at least the first four shock
  interactions. Use a fixed D/32-or-better acoustic corridor after the
  propagation gate sets the actual requirement.
- Use non-reflecting outer boundaries plus a downstream sponge outside the
  physical analysis region. Freeze the grid after a tag-union startup stage.
- Seed with a very weak broadband axisymmetric impulse only if needed, switch
  it off, and require the tone to continue growing or remain saturated without
  forcing.
- Run until onset and saturation are observed, then record at least 100
  periods. At 8.4 kHz this means at least 11.9 ms of post-saturation data; a
  total simulation near 18-25 ms is a reasonable initial budget.

### Stage 2: only then attempt Mj=1.42 in 3D

Use the same validated nozzle/reflection model, a resolved lip annulus and a
weak broadband `m=1` seed. Ring probes must establish that `m=1` dominates
`m=0` and the Cartesian `m=4` imprint. An L5-only pilot can test whether a
large-scale helix appears, but it cannot validate growth rate, SPL or failure
to screech.

## Success criteria locked before running

A case counts as a screech reproduction only if all of the following pass:

1. One narrow peak persists in all time blocks after forcing is removed and
   its amplitude reaches a stationary distribution.
2. Frequency agrees with the experimental branch within a predeclared bound,
   initially 3-5%; the azimuthal mode is correct.
3. The same frequency appears in downstream KH, shock-source, upstream-return
   and upstream acoustic probes.
4. A `k-omega` or SPOD analysis separates the downstream and upstream branches
   and closes the phase around the full loop.
5. The source localizes near the experimentally relevant shock-cell region,
   rather than only at a boundary or AMR interface.
6. A far/near acoustic array shows the expected upstream radiation; probe
   spacing is at most one eighth of the shortest wavelength used for direction
   fitting.
7. The result survives grid, time-step, outer-boundary and weak-seed tests.

Until these gates pass, longer runs of the current flush-flange cases will
improve statistics on their existing dynamics but are unlikely to create the
missing physical feedback loop.

## Primary references

- J. Panda, *An experimental investigation of screech noise generation*, JFM
  378 (1999), DOI `10.1017/S0022112098003383`.
- G. Raman, *Cessation of screech in underexpanded jets*, JFM 336 (1997), DOI
  `10.1017/S002211209600451X`.
- J. Berland, C. Bogey and C. Bailly, *Large Eddy Simulation of Screech Tone
  Generation in a Planar Underexpanded Jet*, AIAA 2006-2496.
- X. Li and J. Gao, *Numerical simulation of the generation mechanism of
  axisymmetric supersonic jet screech tones*, Physics of Fluids 17, 085105
  (2005), DOI `10.1063/1.2033909`.
- J. Mackenzie et al., *A unifying theory of jet screech*, JFM 2022.
