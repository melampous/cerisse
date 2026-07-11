#include <AMReX_FluxRegister.H>
#include <AMReX_ParmParse.H>  // runtime knobs (RZ near-axis dissipation)
#include <AMReX_FillPatchUtil.H>  // coarse-fine ghosts for the overlap path (FPinfo patches, FillPatchInterp)
#include <CNS.h>
#include <prob.h>
#include <limits>  // for std::numeric_limits (RZ divergence floor)

// IBM and EB marker paths use different array types; combined mode is not yet supported.
#if defined(AMREX_USE_GPIBM) && defined(CNS_USE_EB)
#error "compute_rhs.cpp expects exactly one of AMREX_USE_GPIBM or CNS_USE_EB to be enabled"
#endif

#ifdef AMREX_USE_GPIBM
#include <ibm_solver.h>
#endif


using namespace amrex;

namespace {

// Detection idiom: does the flux scheme advertise rz_pressure_split_capable?
// Defaults to FALSE for any flux that does not define the trait (keep_euler_t,
// no_euler_t, and any future flux). Only weno_t opts in (=true). This guards the
// RZ pressure-split (which removes the radial pressure flux only in weno_t) from
// being silently enabled with an incompatible flux -> see compute_rhs guard.
template <typename T, typename = void>
struct rz_psplit_capable : std::false_type {};
template <typename T>
struct rz_psplit_capable<T, std::void_t<decltype(T::rz_pressure_split_capable)>>
    : std::bool_constant<T::rz_pressure_split_capable> {};

template <typename PrimArrayT, typename RhsArrayT, typename ClosuresT,
          typename ProbParmT>
auto try_rhs_nscbc(int, const Geometry& geom, const MFIter& mfi,
                   PrimArrayT const& prims, RhsArrayT const& rhs,
                   const ClosuresT* closures, const ProbParmT* pparm,
                   const Real dt, const Real time)
    -> decltype(rhs_nscbc(geom, mfi, prims, rhs, closures, pparm, dt, time),
                void())
{
  rhs_nscbc(geom, mfi, prims, rhs, closures, pparm, dt, time);
}

template <typename PrimArrayT, typename RhsArrayT, typename ClosuresT,
          typename ProbParmT>
void try_rhs_nscbc(long, const Geometry&, const MFIter&,
                   PrimArrayT const&, RhsArrayT const&,
                   const ClosuresT*, const ProbParmT*,
                   const Real, const Real)
{
}

// ---------------------------------------------------------------------------
// Detection idiom for the region-parameterized (Box instead of MFIter) flux
// overloads used by the comm/comp-overlap path. Flux schemes that do not
// provide them (e.g. Riemann/Rusanov/Skew) simply report false, and
// compute_rhs_overlap aborts at runtime instead of failing to compile.
template <typename T, typename = void>
struct has_region_eflux : std::false_type {};
template <typename T>
struct has_region_eflux<
    T, std::void_t<decltype(std::declval<T&>().eflux(
           std::declval<const Geometry&>(), std::declval<const Box&>(),
           std::declval<const Array4<const Real>&>(),
           std::declval<std::array<FArrayBox*, AMREX_SPACEDIM> const&>(),
           std::declval<const Array4<Real>&>(),
           std::declval<const PROB::ProbClosures*>()))>> : std::true_type {};

template <typename T, typename = void>
struct has_region_dflux : std::false_type {};
template <typename T>
struct has_region_dflux<
    T, std::void_t<decltype(std::declval<T&>().dflux(
           std::declval<const Geometry&>(), std::declval<const Box&>(),
           std::declval<const Array4<Real>&>(),
           std::declval<std::array<FArrayBox*, AMREX_SPACEDIM> const&>(),
           std::declval<const Array4<Real>&>(),
           std::declval<const PROB::ProbClosures*>()))>> : std::true_type {};

// Local equivalents of AMReX's FArrayBox patch makers used by
// FillPatchTwoLevels / FillPatcher (make_mf_crse_patch / make_mf_fine_patch).
// Reproduced here because their namespace differs across AMReX versions
// (anonymous namespace in older snapshots, amrex::detail in newer ones).
inline MultiFab ovl_make_mf_crse_patch(FabArrayBase::FPinfo const& fpc,
                                       int ncomp) {
  return MultiFab(fpc.ba_crse_patch, fpc.dm_patch, ncomp, 0, MFInfo(),
                  *fpc.fact_crse_patch);
}
inline MultiFab ovl_make_mf_fine_patch(FabArrayBase::FPinfo const& fpc,
                                       int ncomp) {
  return MultiFab(fpc.ba_fine_patch, fpc.dm_patch, ncomp, 0, MFInfo(),
                  *fpc.fact_fine_patch);
}

// Euler + diffusive fluxes and finite-volume flux divergence on an explicit
// cell region rbx (subset of a tilebox). Cells (and the faces interior to)
// the optional 'skip' box are not touched — used by the overlap shell pass,
// which covers the whole tilebox but skips the interior region already
// handled in Pass 1. Cartesian only (the overlap path aborts for RZ).
// Per-face and per-cell numerics are identical to the corresponding legacy
// compute_rhs code, so splitting a tilebox into interior + shell regions
// produces bitwise-identical RHS values (seam faces are recomputed from the
// same prims).
template <typename RhsT>
void region_flux_div(RhsT& prob_rhs, const Geometry& geom, const Box& rbx,
                     Array4<Real> const& prims, Array4<Real> const& rhs,
                     const PROB::ProbClosures* cls_d, const int ncons,
                     const Box& skip = Box())
{
  if constexpr (has_region_eflux<RhsT>::value && has_region_dflux<RhsT>::value) {
    const Box skipbox = skip;
    const bool skip_ok = skipbox.ok();

    // flux arrays sized on the faces of this region only; zero only the
    // faces this call will actually compute (masked like the flux kernels
    // themselves) — skipped faces are never read by the masked divergence
    std::array<FArrayBox, AMREX_SPACEDIM> fluxt;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const Box fbx = amrex::surroundingNodes(rbx, dir);
      fluxt[dir].resize(fbx, ncons, The_Async_Arena());
      auto const& f4 = fluxt[dir].array();
      const IntVect ivd = IntVect::TheDimensionVector(dir);
      ParallelFor(fbx, ncons,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    if (skip_ok) {
                      const IntVect iv{AMREX_D_DECL(i, j, k)};
                      if (skipbox.contains(iv) && skipbox.contains(iv - ivd)) {
                        return;
                      }
                    }
                    f4(i, j, k, n) = Real(0.0);
                  });
    }

    {
      BL_PROFILE_VAR("CNS::compute_rhs::eflux", prof_eflux);
      prob_rhs.eflux(geom, rbx, prims,
                     {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, rhs,
                     cls_d, skip);
      BL_PROFILE_VAR_STOP(prof_eflux);
    }
    {
      BL_PROFILE_VAR("CNS::compute_rhs::dflux", prof_dflux);
      prob_rhs.dflux(geom, rbx, prims,
                     {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, rhs,
                     cls_d, skip);
      BL_PROFILE_VAR_STOP(prof_dflux);
    }

    const auto dx = geom.CellSizeArray();
    auto const& fx = fluxt[0].array();
#if (AMREX_SPACEDIM >= 2)
    auto const& fy = fluxt[1].array();
#endif
#if (AMREX_SPACEDIM == 3)
    auto const& fz = fluxt[2].array();
#endif
    const int nc = ncons;

#if (AMREX_SPACEDIM == 1)
    const Real invdx = Real(1.0) / dx[0];
    ParallelFor(rbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (skip_ok && skipbox.contains(i, j, k)) { return; }
      for (int n = 0; n < nc; ++n) {
        rhs(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
      }
    });
#elif (AMREX_SPACEDIM == 2)
    const Real invdx = Real(1.0) / dx[0];
    const Real invdy = Real(1.0) / dx[1];
    ParallelFor(rbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (skip_ok && skipbox.contains(i, j, k)) { return; }
      for (int n = 0; n < nc; ++n) {
        rhs(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
        rhs(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
      }
    });
#else
    const Real invdx = Real(1.0) / dx[0];
    const Real invdy = Real(1.0) / dx[1];
    const Real invdz = Real(1.0) / dx[2];
    ParallelFor(rbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (skip_ok && skipbox.contains(i, j, k)) { return; }
      for (int n = 0; n < nc; ++n) {
        rhs(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
        rhs(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
        rhs(i, j, k, n) += (fz(i, j, k, n) - fz(i, j, k + 1, n)) * invdz;
      }
    });
#endif
  } else {
    amrex::ignore_unused(prob_rhs, geom, rbx, prims, rhs, cls_d, ncons, skip);
    amrex::Abort(
        "cns.overlap_comm=1: the configured flux scheme does not provide "
        "region-parameterized eflux/dflux overloads (only WENO/TENO eflux "
        "and viscous_t dflux do). Disable cns.overlap_comm.");
  }
}

}  // namespace

// Since we do not want to use expensive cudaMemCopy, we are storing all our
// data on the GPU to begin with. Concurrency on GPU using streams, parallel
// computation and data transfer, is not useful then. Therefore, we can have all
// grid point computations, per fab, in a single MFIter loop (single stream).

void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxRegister* fr_as_crse,
                      FluxRegister* fr_as_fine, Real stage_time) {
  BL_PROFILE("CNS::compute_rhs()");

  // Variables
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  const PROB::ProbParm* pparm_d = CNS::d_prob_parm;

  // time
  // RK stage time is explicit.  StateData::curTime() is not reliable for the
  // first stage after swapTimeLevels(): it already points at t^{n+1} while
  // FillPatch reads U^n at t^n.
  const Real cur_time = stage_time;

#ifdef AMREX_USE_GPIBM
  // Convert conserved variables to primitives level-wide, then apply IBM
  // ghost-point corrections in a single pass before the MFIter loop.
  // This avoids redundant per-fab conversions and ensures all ghost-point
  // data are consistent when each fab's flux kernel executes.
  MultiFab prims_mf(statemf.boxArray(), statemf.DistributionMap(),
                    cls_h.NPRIM, cls_h.NGHOST,
                    MFInfo().SetArena(The_Async_Arena()));

  // In 2D, prob_initdata and bcnormal do not write UMZ (z-momentum). That
  // leaves UMZ in valid cells and physical-BC ghosts as uninitialized memory
  // — occasionally NaN. cons2prims reads UMZ unconditionally, and a NaN
  // there cascades: uz = NaN → rhoke = NaN → E' = NaN → all prims NaN at
  // that cell → WENO stencils produce NaN flux → RHS = NaN → blow-up at
  // step 1 (typically exposed by AMR because extra fab allocations deplete
  // zero-pages and expose stale/NaN memory). Zero the UMZ component each
  // call to guarantee a clean 2D slice regardless of prob.h conventions.
#if (AMREX_SPACEDIM < 3)
  statemf.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, statemf.nGrow());
#endif
  {
    BL_PROFILE_VAR("CNS::compute_rhs::cons2prims", prof_cons2prims);
    for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
      cls_h.cons2prims(mfi, statemf.array(mfi), prims_mf.array(mfi));
    }
    BL_PROFILE_VAR_STOP(prof_cons2prims);
  }
  {
    // Sync wall-motion time before GP reconstruction so that
    // compute_surfIB() sees the correct wall velocity.
#ifdef CNS_USE_FSI
    PROB::Motion::sim_time = cur_time;
#endif

    BL_PROFILE_VAR("IBM::computeAllGPs", prof_gp);
    IBM::ib.computeAllGPs(prims_mf, cls_d, level);
    BL_PROFILE_VAR_STOP(prof_gp);
  }

  // Validation hook: poison only unreconstructed interior-solid primitives.
  // A marker-safe IBM flux must produce bitwise-identical fluid RHS values
  // with this enabled.  The default is off and adds no production kernel.
  static const bool poison_interior_solid = [] {
    int value = 0;
    ParmParse pp("ib");
    pp.query("poison_interior_solid", value);
    return value != 0;
  }();
  if (poison_interior_solid) {
    const Real poison = std::numeric_limits<Real>::quiet_NaN();
    auto& marker_mf = *IBM::ib.bmf_a[level];
    for (MFIter mfi(prims_mf, false); mfi.isValid(); ++mfi) {
      const Box bxg = mfi.growntilebox(cls_h.NGHOST);
      const auto prims = prims_mf.array(mfi);
      const auto marker = marker_mf.array(mfi);
      ParallelFor(bxg, cls_h.NPRIM,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    if (marker(i, j, k, 0) != 0 &&
                        marker(i, j, k, 1) == 0) {
                      prims(i, j, k, n) = poison;
                    }
                  });
    }
  }
#endif

  //...................................................................
  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& state = statemf.array(mfi);

    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls_h.NGHOST);
#ifdef CNS_USE_EB     
    const Box& bxflux = mfi.growntilebox(cls_h.NGHOST+1); // add 1 cell 
#else
    const Box& bxflux = mfi.growntilebox(cls_h.NGHOST); 
#endif    
    
    // primitives and fluxes arrays
#ifdef AMREX_USE_GPIBM
    // Prims already filled + GP-corrected in pre-loop; just alias
    Array4<Real> const& prims = prims_mf.array(mfi);
#else
    FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
    Array4<Real> const& prims = primf.array();
#endif

    
#ifdef CNS_USE_EB     
    // auxiliary arrays for redistribution 
    FArrayBox divcfab(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& divc = divcfab.array();    

    // store array cons 
    FArrayBox consfab(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& cons = consfab.array();    
    amrex::ParallelFor(bxg, cls_h.NCONS,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        cons(i,j,k,n) = state(i,j,k,n);
      });
#endif

    // flux arrays  
    std::array<FArrayBox ,AMREX_SPACEDIM> fluxt;
    for (int dir=0; dir < AMREX_SPACEDIM; ++dir)
    {
      fluxt[dir].resize(amrex::surroundingNodes(bxflux, dir),cls_h.NCONS, The_Async_Arena() );
      fluxt[dir].setVal<RunOn::Device>(0.);
    }
     
    // We want to minimise function calls. So, we call prims2cons, flux and
    // source term evaluations once per fab from CPU, to be run on GPU.
#ifndef AMREX_USE_GPIBM
    cls_h.cons2prims(mfi, state, prims);
#endif

    // Geometry markers (one of IBM/EB is expected to be enabled).
    // Alias geoMarkers to the underlying marker MultiFab to avoid
    // allocating an auxiliary fab + doing a device copy.
#if (AMREX_USE_GPIBM && !CNS_USE_EB)
    auto& marker_mf = *IBM::ib.bmf_a[level];
    auto const& geoMarkers = marker_mf.array(mfi);
#elif (CNS_USE_EB && !AMREX_USE_GPIBM)
    // EB markers (bool) -> convert to uint8_t for interface compatibility
    // with eflux_ibm/dflux_ibm which expect Array4<uint8_t>.
    // EBM core uses bool internally; we convert at the boundary here
    // to avoid modifying EBM code that is shared with other users.
    auto& eb_marker_mf = *EBM::eb.bmf_a[level];
    auto const& ebBoolMarkers = eb_marker_mf.array(mfi);
    BaseFab<uint8_t> geoMarkerFab(bxg, 2, The_Async_Arena());
    auto const& geoMarkers = geoMarkerFab.array();
    amrex::ParallelFor(bxg, 2,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        geoMarkers(i,j,k,n) = static_cast<uint8_t>(ebBoolMarkers(i,j,k,n));
      });
#endif
  
    // Euler/Diff Fluxes including boundary/discontinuity corrections
    // WARNING: state is the U array (cons)
    {
    BL_PROFILE_VAR("CNS::compute_rhs::eflux", prof_eflux);
#if (AMREX_USE_GPIBM || CNS_USE_EB)
    prob_rhs.eflux_ibm(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d, geoMarkers);
#else
    prob_rhs.eflux(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d);
#endif
    BL_PROFILE_VAR_STOP(prof_eflux);
    }
    {
    BL_PROFILE_VAR("CNS::compute_rhs::dflux", prof_dflux);
#if (AMREX_USE_GPIBM || CNS_USE_EB)
    prob_rhs.dflux_ibm(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d, geoMarkers);
#else
    prob_rhs.dflux(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d);
#endif
    BL_PROFILE_VAR_STOP(prof_dflux);
    }

    // compute rhs as finite-volume flux divergence, i.e.
    //   rhs += ((F·A)_lo - (F·A)_hi) / V
    // WARNING: state is now the RHS array
    // set RHS=0 (everywhere including ghost points)
    ParallelFor(bxg, cls_h.NCONS, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {state(i,j,k,n) = 0.0;});

    // Geometry-aware finite-volume flux divergence.  Metrics are evaluated
    // per-cell inside the kernel, avoiding temporary volume/area arrays.
    //   Cartesian : standard uniform-cell form, dU/dt += (F_lo - F_hi)/dx.
    //   RZ (2D)   : cylindrical (r,z) with r_i = prob_lo[0] + (i+1/2)*dr;
    //               the r-flux term uses the metric-consistent form
    //               -(1/r) d(rF_r)/dr ≈ 2(r_lo F_lo - r_hi F_hi)/(r_hi^2 - r_lo^2).
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    auto const& fx = fluxt[0].array();
#if (AMREX_SPACEDIM >= 2)
    auto const& fy = fluxt[1].array();
#endif
#if (AMREX_SPACEDIM == 3)
    auto const& fz = fluxt[2].array();
#endif

    const bool is_rz = geom.IsRZ();

    // Performance note: a 3D kernel with an inner loop over components reuses
    // per-cell metrics across all NCONS equations, avoiding redundant metric
    // evaluations that would arise from a 4D (i,j,k,n) ParallelFor.
    const int ncons = cls_h.NCONS;

#if (AMREX_SPACEDIM == 2)
    if (is_rz) {

        const Real dr = dx[0];
        const Real dz = dx[1];
        const Real r0 = prob_lo[0];
        const Box domain = geom.Domain();
        const int ilo = domain.smallEnd(0);
        const int ihi = domain.bigEnd(0);
        const int qpres   = PROB::ProbClosures::QPRES;
        const int qrho    = PROB::ProbClosures::QRHO;
        const int qu      = PROB::ProbClosures::QU;   // radial velocity (dir 0 = r)
        const int qc      = PROB::ProbClosures::QC;   // sound speed
        const int umom_r  = PROB::ProbClosures::UMX;
        const Real inv_dz = Real(1.0) / dz;

        // --- (Route B) near-axis radial-momentum dissipation -----------------
        // Damps the spurious odd-even / carbuncle-like u_r seeded at strong
        // on-axis shocks by the singular RZ axis discretisation (the p/r source
        // amplified by 1/r_c=2/dr at i=0).  OFF by default (coeff 0) so no other
        // case is affected; enable per-case with cns.rz_axis_diss=<O(0.5)>.
        static const Real s_rz_diss = []{ Real v=Real(0.0);
            amrex::ParmParse pp("cns"); pp.query("rz_axis_diss", v); return v; }();
        static const int  s_rz_nc   = []{ int n=3;
            amrex::ParmParse pp("cns"); pp.query("rz_axis_diss_ncell", n); return n; }();
        // Experimental WENO/TENO-only pressure split:
        //   radial momentum pressure force = -dp/dr, not
        //   -(1/r)d(rp)/dr + p/r.  The matching Weno.h path removes p from
        //   the r-momentum r-face flux.  Leave OFF unless that path is active.
        static const int s_rz_pressure_split = []{ int v=0;
            amrex::ParmParse pp("cns"); pp.query("rz_pressure_split", v); return v; }();
        // Near-axis well-balanced geometric pressure source (default OFF, 0 cells).
        // In the first N radial cells it replaces the cell-centre source p_c/r
        // with the flux-face-average 1/2(F_r(i)+F_r(i+1))/r.  Combined with the
        // metric divergence this collapses to the Cartesian flux difference
        // (F_r(i)-F_r(i+1))/dr, which is exactly well-balanced (uniform state ->
        // 0) and, since u_r->0 at the axis so F_r=rho*u_r^2+p -> p, reproduces
        // -dp/dr WITHOUT the 1/r-amplified [p_c - mean(p_face)] residual (the
        // eq.12 seed of the spurious on-axis radial momentum at the Mach disk).
        // Uses the SAME face fluxes as the divergence (exact consistency); bulk
        // cells (i-ilo >= N) keep the standard p_c/r source so the validated
        // shock structure / standoff is untouched.  Pairs with the baseline flux
        // (NOT rz_pressure_split, which removes p from the flux entirely).
        static const int s_rz_wb_ncell = []{ int n=0;
            amrex::ParmParse pp("cns"); pp.query("rz_wb_axis_ncell", n); return n; }();
        const Real rz_diss = s_rz_diss;   // local copies -> captured by value (GPU-safe)
        const int  rz_nc   = s_rz_nc;
        const int  rz_wb_ncell = s_rz_wb_ncell;
        const bool rz_pressure_split = (s_rz_pressure_split != 0);
        // HARD GUARD: the -dp/dr pressure-split form below MUST be paired with a
        // flux that removed the radial pressure flux. Only weno_t::eflux does so
        // (rz_pressure_split_capable=true). With any other flux scheme the radial
        // pressure would be counted twice (in the flux AND here) -> wrong equations.
        if (rz_pressure_split && !rz_psplit_capable<PROB::ProbRHS>::value) {
          amrex::Abort("cns.rz_pressure_split is experimental and WENO/TENO-only: "
                       "it requires weno_t::eflux to remove the radial pressure flux. "
                       "The active flux scheme does not, so the pressure would be "
                       "double-counted. Disable cns.rz_pressure_split or use weno_t.");
        }
        const Real dr_loc  = dr;
        const Real inv_2dr = Real(0.5) / dr;

        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                    // Axisymmetric RZ: dir=0 is r, dir=1 is z.
                    // In this branch:
                    //   fx = flux density through r-faces (i±1/2)
                    //   fy = flux density through z-faces (j±1/2)
                    const Real r_lo = r0 + Real(i) * dr;
                    const Real r_hi = r_lo + dr;
                    const Real r_c  = r_lo + Real(0.5) * dr;

                    // r2diff = r_hi^2 - r_lo^2 = (r_hi + r_lo)*(r_hi - r_lo)
                    //        = (r_hi + r_lo) * dr
                    const Real r2diff = (r_hi + r_lo) * dr;
                    // For well-posed RZ problems prob_lo[0]=0 and i>=0, so
                    // r2diff is strictly positive; the floor guards against
                    // pathological setups and is inactive in normal use.
                    const Real tiny = amrex::max(std::numeric_limits<Real>::min(), Real(1.0e-14) * dr * dr);
                    const Real inv_r2diff = Real(1.0) / amrex::max(r2diff, tiny);
                    const Real tiny_r = Real(1.0e-14) * dr;
                    const Real inv_r = Real(1.0) / amrex::max(r_c, tiny_r);

                    for (int n = 0; n < ncons; ++n) {
                        // FV form with metrics. After cancellation of common factors:
                        // - z term reduces to Cartesian difference / dz
                        // - r term is (2*r*F_r)|_lo - (2*r*F_r)|_hi over (r_hi^2 - r_lo^2)
                        
                        // r-direction: metric-consistent RZ divergence
                        //   -(1/r) d(rF_r)/dr ≈ 2(r_lo*F_lo - r_hi*F_hi) / (r_hi^2 - r_lo^2)
                        Real rhs_rz = (Real(2.0) * (r_lo * fx(i, j, k, n) - r_hi * fx(i + 1, j, k, n))) * inv_r2diff;
                        state(i, j, k, n) += rhs_rz;

                        // Axisymmetric Euler geometric source for radial momentum:
                        // +p/r term is not contained in -(1/r) d(r F_r)/dr - dF_z/dz
                        // when F_r uses the standard Cartesian-form momentum flux.
                        if (n == umom_r) {
                            if (rz_pressure_split) {
                                Real dpdr;
                                if (i <= ilo) {
                                    dpdr = (-Real(3.0) * prims(i, j, k, qpres)
                                            + Real(4.0) * prims(i + 1, j, k, qpres)
                                            - prims(i + 2, j, k, qpres)) * inv_2dr;
                                } else if (i >= ihi) {
                                    dpdr = ( Real(3.0) * prims(i, j, k, qpres)
                                            - Real(4.0) * prims(i - 1, j, k, qpres)
                                            + prims(i - 2, j, k, qpres)) * inv_2dr;
                                } else {
                                    // A2-a: monotone (minmod-limited) radial pressure
                                    // gradient. Well-balanced (uniform p -> 0) AND
                                    // non-oscillatory at shocks -- unlike the central
                                    // difference, which Gibbs-oscillates at the on-axis
                                    // Mach disk / bow shock and (via the 1/r metric)
                                    // seeds the spurious near-axis radial momentum.
                                    // Consistent with the cell-centre pressure that
                                    // weno_t removed from the advective r-momentum flux.
                                    const Real gL = prims(i,   j, k, qpres) - prims(i-1, j, k, qpres);
                                    const Real gR = prims(i+1, j, k, qpres) - prims(i,   j, k, qpres);
                                    const Real slope = (gL * gR <= Real(0.0)) ? Real(0.0)
                                        : (amrex::Math::abs(gL) < amrex::Math::abs(gR) ? gL : gR);
                                    dpdr = slope / dr_loc;
                                }
                                state(i, j, k, n) -= dpdr;
                            } else if (rz_wb_ncell > 0 && (i - ilo) < rz_wb_ncell) {
                                // Near-axis well-balanced source (see derivation
                                // above): flux-face-average in place of p_c/r.
                                // metric-div(line ~310) + this == (F_r(i)-F_r(i+1))/dr.
                                state(i, j, k, n) += Real(0.5)
                                    * (fx(i, j, k, n) + fx(i + 1, j, k, n)) * inv_r;
                            } else {
                                state(i, j, k, n) += prims(i, j, k, qpres) * inv_r;
                            }

                            // (Route B) near-axis radial-momentum dissipation:
                            // artificial radial viscosity nu = rz_diss*(|u_r|+a)*dr,
                            // RHS += nu/dr^2 * d2(rho*u_r)/dr2, ramped 1->0 over rz_nc cells.
                            if (rz_diss > Real(0.0) && i <= rz_nc) {
                                const Real mL = prims(i-1,j,k,qrho)*prims(i-1,j,k,qu);
                                const Real mC = prims(i,  j,k,qrho)*prims(i,  j,k,qu);
                                const Real mR = prims(i+1,j,k,qrho)*prims(i+1,j,k,qu);
                                const Real lap   = mL - Real(2.0)*mC + mR;
                                const Real speed = std::abs(prims(i,j,k,qu)) + prims(i,j,k,qc);
                                const Real wgt   = Real(1.0) - Real(i)/Real(rz_nc + 1);
                                state(i, j, k, n) += rz_diss * wgt * (speed / dr_loc) * lap;
                            }
                        }

                        state(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * inv_dz;
                    }
                });
    } else
#endif
    {
#if (AMREX_SPACEDIM == 1)
        const Real invdx = Real(1.0) / dx[0];
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < ncons; ++n) {
                        state(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
                    }
                });
#elif (AMREX_SPACEDIM == 2)
        const Real invdx = Real(1.0) / dx[0];
        const Real invdy = Real(1.0) / dx[1];
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < ncons; ++n) {
                        state(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
                        state(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
                    }
                });
#else
        const Real invdx = Real(1.0) / dx[0];
        const Real invdy = Real(1.0) / dx[1];
        const Real invdz = Real(1.0) / dx[2];
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < ncons; ++n) {
                        state(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
                        state(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
                        state(i, j, k, n) += (fz(i, j, k, n) - fz(i, j, k + 1, n)) * invdz;
                    }
                });
#endif
    }

    // RZ viscous geometric source (hoop-stress and related terms) not captured
    // by the metric FV divergence of the face fluxes.
    if (is_rz) {
        prob_rhs.rz_geometric_source(geom, mfi, prims, state, cls_d);
    }
                      
#if CNS_USE_EB    
    // internal geometry fluxes
    const Box&  ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const auto& flag = (*EBM::eb.ebflags_a[level])[mfi];
    FabType t = flag.getType(ebbox);

    const bool fab_with_eb = (FabType::singlevalued == t);  
    // EB flux     
    if (fab_with_eb) {
      EBM::eb.ebflux(geom,mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},state, cls_d,level);
    }

    // redistribution 
    // WARNING: state is  the RHS array, prims is the prims 
    // compute divc here
    amrex::ParallelFor(bxg, cls_h.NCONS,  
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      divc(i,j,k,n) = state(i,j,k,n);
    });  
    
    // do redistribution only in box with EB
    if (eb_redistribution && fab_with_eb){      

      EBM::eb.redist(geom,mfi,cons,divc, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},
                    state, cls_d,level,dt,h_phys_bc);
    }                    
#endif 

    // Optional RHS-level NSCBC hook.  This is intentionally after the finite-
    // volume flux divergence (and EB redistribution when present) so a problem
    // can directly replace boundary-cell dU/dt instead of only influencing
    // ghost-cell states before the Riemann solve.
    try_rhs_nscbc(0, geom, mfi, prims, state, cls_d, pparm_d, dt, cur_time);

    // Source terms (body forces, chemistry, etc.)
#if (AMREX_USE_GPIBM || CNS_USE_EB)
    prob_rhs.src(geom,mfi, prims, state, cls_d, dt, cur_time, geoMarkers);
#else
    prob_rhs.src(geom,mfi, prims, state, cls_d, dt, cur_time);
#endif

    // Zero the RHS inside solid cells (state holds RHS at this point)
#if (AMREX_USE_GPIBM || CNS_USE_EB)
        amrex::ParallelFor(bxg,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // IBM: geoMarkers(i,j,k,0) stores geometry_index in a uint8_t.
            // EB : geoMarkers(i,j,k,0) is a bool covered-cell marker.
            // In both cases, nonzero => solid.
            if (geoMarkers(i,j,k,0) != 0) {
                for (int n = 0; n < ncons; ++n) {
                    state(i,j,k,n) = Real(0.0);
                }
            }
        });
#endif
  }

}

// ============================================================================
// Communication/computation overlap variant of (FillPatch + compute_rhs) for
// one RK stage (cns.overlap_comm=1). Design:
//
//   (a)  valid region <- SRC (plain local copy; the RK3 driver arranges state
//        time levels so the legacy FillPatch's valid fill is exactly this)
//   (b1) coarse-fine boundary ghosts, phase 1 (level > 0): start non-blocking
//        fetches (ParallelCopy_nowait) of any coarse source data not in the
//        per-level cache. Split-phase re-implementation of
//        amrex::FillPatcher::fillCoarseFineBoundary (AMReX 25.12) — same
//        FPinfo patches, same time-interp formula, same coarse physbc, same
//        FillPatchInterp call, so the result is bitwise-identical; only the
//        communication is overlapped. Cache reset in CNS::post_timestep.
//   (c)  same-level ghost exchange started with FillBoundary_nowait
//   (d)  PASS 1 (overlapped with both exchanges): prims on the VALID box,
//        RHS zeroed on the interior, flux + divergence on the interior
//        region ibx = grow(bx, -NGHOST) whose stencils never reach ghosts
//   (e)  FillBoundary_finish  (the MPI wait now overlaps Pass-1 kernels)
//   (b2) coarse-fine boundary ghosts, phase 2: finish the fetches, time-
//        interpolate the coarse patch, coarse physbc, space-interpolate and
//        copy into statemf's coarse-fine ghosts
//   (f)  physical-BC ghosts via StateDataPhysBCFunct (same functor and time
//        as the legacy FillPatch applies last)
//   (g)  PASS 2: ghost-ring prims from the fresh ghost cons (one masked
//        kernel), ghost-ring + shell RHS zeroed, then ONE masked flux +
//        divergence call over the tilebox that skips the interior region
//        already handled in Pass 1; then the NSCBC hook and source terms
//        once per fab (legacy ordering)
//
// Seam faces between interior and shell are computed twice from identical
// prims, giving bitwise-identical fluxes => conservation preserved and the
// result is bitwise-identical to FillPatch + compute_rhs.
// ============================================================================
void CNS::compute_rhs_overlap(MultiFab& statemf, Real dt,
                              FluxRegister* fr_as_crse,
                              FluxRegister* fr_as_fine, Real t_fill,
                              MultiFab& SRC, MultiFab& prims_mf) {
  BL_PROFILE("CNS::compute_rhs_overlap()");

  amrex::ignore_unused(fr_as_crse, fr_as_fine);  // reflux not supported here

#if defined(AMREX_USE_GPIBM) || defined(CNS_USE_EB)
  amrex::ignore_unused(statemf, dt, t_fill, SRC, prims_mf);
  amrex::Abort("cns.overlap_comm=1 is not supported in IBM/EB builds");
#else
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  const PROB::ProbParm* pparm_d = CNS::d_prob_parm;

  // The overlap driver already supplies the exact RK abscissa.
  const Real cur_time = t_fill;

  const int ncons = cls_h.NCONS;
  const int ng = cls_h.NGHOST;

  if (geom.IsRZ()) {
    amrex::Abort("cns.overlap_comm=1 does not support RZ geometry");
  }
  AMREX_ALWAYS_ASSERT(statemf.nGrowVect().allGE(IntVect(ng)));
  AMREX_ALWAYS_ASSERT(prims_mf.nGrowVect().allGE(IntVect(ng)));

  // ---- (a) valid-region copy (local, no communication) --------------------
  MultiFab::Copy(statemf, SRC, 0, 0, ncons, 0);

  // ---- (c) same-level ghost exchange, non-blocking -------------------------
  // Posted FIRST so the sends leave as soon as the pack kernels finish; the
  // coarse-fine interpolation chain below would otherwise delay the pack
  // synchronization (and thus the neighbours' FillBoundary_finish).
  statemf.FillBoundary_nowait(geom.periodicity());

  // ---- (b1) coarse-fine ghosts, phase 1: start coarse fetches --------------
  const StateDescriptor& desc = AmrLevel::desc_lst[State_Type];
  InterpBase* cfb_mapper = desc.interp(0);
  FabArrayBase::FPinfo const* cfb_fpc = nullptr;
  Vector<MultiFab*> cfb_pending;  // fetches to finish in phase 2
  IntVect cfb_ratio;
  bool cfb_precomputed = false;  // chain done + ghost copy in flight from B1

  // Time-interp of the cached coarse patches + coarse physbc + space-interp
  // into m_ovl_cfb_fine. Replicates FillPatcher::fillCoarseFineBoundary
  // (AMReX 25.12) step for step, so the values are bitwise-identical to the
  // legacy FillPatchTwoLevels result. Requires all needed coarse source
  // patches present in m_ovl_cfb_data.
  auto run_cfb_chain = [&]() {
    if (m_ovl_cfb_tmp == nullptr) {
      m_ovl_cfb_tmp = std::make_unique<MultiFab>(
          ovl_make_mf_crse_patch(*cfb_fpc, ncons));
    }
    if (m_ovl_cfb_fine == nullptr) {
      m_ovl_cfb_fine = std::make_unique<MultiFab>(
          ovl_make_mf_fine_patch(*cfb_fpc, ncons));
    }

    AmrLevel& crse_level = parent->getLevel(level - 1);
    const Geometry& cgeom = crse_level.Geom();

    int const ng_space_interp = 8;  // must match FillPatcher
    Box domain = cgeom.growPeriodicDomain(ng_space_interp);
    domain.convert(statemf.ixType());

    int idata = -1;
    if (m_ovl_cfb_data.size() == 1) {
      idata = 0;
    } else if (m_ovl_cfb_data.size() == 2) {
      Real const teps =
          std::abs(m_ovl_cfb_data[1].first - m_ovl_cfb_data[0].first) *
          Real(1.e-3);
      if (t_fill > m_ovl_cfb_data[0].first - teps &&
          t_fill < m_ovl_cfb_data[0].first + teps) {
        idata = 0;
      } else if (t_fill > m_ovl_cfb_data[1].first - teps &&
                 t_fill < m_ovl_cfb_data[1].first + teps) {
        idata = 1;
      } else {
        idata = 2;
      }
    }

    if (idata == 0 || idata == 1) {
      auto const& dst = m_ovl_cfb_tmp->arrays();
      auto const& src = m_ovl_cfb_data[idata].second->const_arrays();
      amrex::ParallelFor(
          *m_ovl_cfb_tmp, IntVect(0), ncons,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (domain.contains(i, j, k)) {
              dst[bi](i, j, k, n) = src[bi](i, j, k, n);
            }
          });
    } else if (idata == 2) {
      Real t0 = m_ovl_cfb_data[0].first;
      Real t1 = m_ovl_cfb_data[1].first;
      Real alpha = (t1 - t_fill) / (t1 - t0);
      Real beta = (t_fill - t0) / (t1 - t0);
      auto const& a = m_ovl_cfb_tmp->arrays();
      auto const& a0 = m_ovl_cfb_data[0].second->const_arrays();
      auto const& a1 = m_ovl_cfb_data[1].second->const_arrays();
      amrex::ParallelFor(
          *m_ovl_cfb_tmp, IntVect(0), ncons,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (domain.contains(i, j, k)) {
              a[bi](i, j, k, n) =
                  alpha * a0[bi](i, j, k, n) + beta * a1[bi](i, j, k, n);
            }
          });
    } else {
      amrex::Abort(
          "compute_rhs_overlap: invalid coarse-fine cache (more than two "
          "coarse times — was resetOvlCFB() skipped?)");
    }
    Gpu::streamSynchronize();

    StateDataPhysBCFunct cbc(crse_level.get_state_data(State_Type), 0, cgeom);
    cbc(*m_ovl_cfb_tmp, 0, ncons, m_ovl_cfb_tmp->nGrowVect(), t_fill, 0);

    FillPatchInterp(*m_ovl_cfb_fine, 0, *m_ovl_cfb_tmp, 0, ncons, IntVect(0),
                    cgeom, geom,
                    amrex::grow(amrex::convert(geom.Domain(), statemf.ixType()),
                                IntVect(ng)),
                    cfb_ratio, cfb_mapper, desc.getBCs(), 0);
  };

  if (level > 0) {
    BL_PROFILE("CNS::compute_rhs_overlap::cfb_start");
    AmrLevel& crse_level = parent->getLevel(level - 1);
    const Geometry& cgeom = crse_level.Geom();

    if (level > 1 &&
        !amrex::ProperlyNested(crse_ratio, parent->blockingFactor(level), ng,
                               statemf.ixType(), cfb_mapper)) {
      amrex::Abort(
          "cns.overlap_comm=1: grids are not properly nested for the "
          "coarse-fine fill; increase amr.blocking_factor or disable "
          "cns.overlap_comm");
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      cfb_ratio[idim] =
          geom.Domain().length(idim) / cgeom.Domain().length(idim);
    }

    // Same FPinfo as FillPatcher::getFPinfo (cached inside amrex, keyed by
    // the BoxArray/DistributionMapping/nghost/coarsener)
    MultiFab sfine(grids, dmap, 1, IntVect(ng), MFInfo().SetAlloc(false));
    const InterpolaterBoxCoarsener& coarsener =
        cfb_mapper->BoxCoarsener(cfb_ratio);
    cfb_fpc = &FabArrayBase::TheFPinfo(sfine, sfine, IntVect(ng), coarsener,
                                       geom, cgeom, nullptr);

    if (!cfb_fpc->ba_crse_patch.empty()) {
      Vector<MultiFab*> smf_crse;
      Vector<Real> stime_crse;
      crse_level.get_state_data(State_Type).getData(smf_crse, stime_crse,
                                                    t_fill);
      for (int i = 0; i < smf_crse.size(); ++i) {
        const Real t = stime_crse[i];
        auto it = std::find_if(m_ovl_cfb_data.begin(), m_ovl_cfb_data.end(),
                               [=](auto const& x) {
                                 return amrex::almostEqual(x.first, t, 5);
                               });
        if (it == m_ovl_cfb_data.end()) {
          auto p = std::make_unique<MultiFab>(
              ovl_make_mf_crse_patch(*cfb_fpc, ncons));
          p->ParallelCopy_nowait(*smf_crse[i], cgeom.periodicity());
          cfb_pending.push_back(p.get());
          m_ovl_cfb_data.emplace_back(t, std::move(p));
        }
      }

      if (cfb_pending.empty()) {
        // All coarse source data already cached (typical for RK stages 2,3
        // and later substeps): run the whole interpolation chain NOW and
        // start the ghost copy non-blocking, so the FPinfo-patch exchange
        // (knapsack-distributed, i.e. real communication) also overlaps
        // Pass 1. A FabArray keeps FillBoundary and ParallelCopy requests
        // in separate handlers, so both can be in flight simultaneously.
        run_cfb_chain();
        statemf.ParallelCopy_nowait(*m_ovl_cfb_fine, 0, 0, ncons, IntVect{0},
                                    IntVect(ng));
        cfb_precomputed = true;
      }
    }
  }

  // No device synchronization is needed before Pass 1: the in-flight
  // exchanges read only (i) neighbouring fabs' valid data within NGHOST of
  // fab interfaces (shell cells — Pass 1 leaves the whole shell untouched
  // and only writes the interior region ibx), (ii) the coarse level's state,
  // and (iii) m_ovl_cfb_fine; their ghost writes target regions no Pass-1
  // kernel touches. All reads/writes are therefore disjoint from Pass 1.

  // ---- (d) PASS 1: interior (overlapped with the ghost exchange) ----------
  {
    BL_PROFILE("CNS::compute_rhs_overlap::pass1");
    for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& state4 = statemf.array(mfi);
      Array4<Real> const& prims4 = prims_mf.array(mfi);

      // prims for ALL valid cells (pointwise; needs no ghosts) — extracted
      // BEFORE any RHS zeroing wipes the conservative data
      cls_h.cons2prims(bx, state4, prims4);

      // interior region: stencils (WENO ng=3 <= NGHOST, viscous halfsten
      // <= NGHOST) never reach outside the valid box
      const Box ibx = amrex::grow(bx, -ng);
      if (ibx.ok()) {
        // zero RHS on the INTERIOR ONLY. The shell (and ghost ring) keep
        // their conservative data until Pass 2: the physical-BC fill in
        // step (f) extrapolates domain ghosts from near-boundary valid
        // cells (all within the shell), so they must stay intact here.
        ParallelFor(ibx, ncons,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                      state4(i, j, k, n) = 0.0;
                    });
        region_flux_div(prob_rhs, geom, ibx, prims4, state4, cls_d, ncons);
      }
    }
  }

  // ---- (e) finish the same-level exchange ----------------------------------
  statemf.FillBoundary_finish();

  // ---- (b2) coarse-fine ghosts, phase 2 ------------------------------------
  // Cache hit (typical): only finish the ghost copy started in (b1).
  // Cache miss (first stage on fresh coarse data): finish the overlapped
  // coarse fetches, run the interpolation chain, and copy the ghosts.
  if (level > 0 && cfb_fpc != nullptr && !cfb_fpc->ba_crse_patch.empty()) {
    BL_PROFILE("CNS::compute_rhs_overlap::cfb_finish");
    if (cfb_precomputed) {
      statemf.ParallelCopy_finish();
    } else {
      for (auto* p : cfb_pending) {
        p->ParallelCopy_finish();
      }
      run_cfb_chain();
      statemf.ParallelCopy(*m_ovl_cfb_fine, 0, 0, ncons, IntVect{0},
                           IntVect(ng));
    }
  }

  // ---- (f) physical-BC ghosts (same functor/time as legacy FillPatch) -----
  {
    BL_PROFILE("CNS::compute_rhs_overlap::physbc");
    StateDataPhysBCFunct physbcf(state[State_Type], 0, geom);
    physbcf(statemf, 0, ncons, statemf.nGrowVect(), t_fill, 0);
  }

  // ---- (g) PASS 2: shell + finalize ----------------------------------------
  {
    BL_PROFILE("CNS::compute_rhs_overlap::pass2");
    for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const Box& bxg = mfi.growntilebox(ng);
      Array4<Real> const& state4 = statemf.array(mfi);
      Array4<Real> const& prims4 = prims_mf.array(mfi);

      // interior region handled in Pass 1 (invalid box if the fab is too
      // small — then the whole tilebox is done here)
      const Box ibx0 = amrex::grow(bx, -ng);
      const Box ibx = ibx0.ok() ? ibx0 : Box();
      const bool have_ibx = ibx.ok();

      // ghost-ring prims from the freshly exchanged conservative data
      // (single masked kernel; valid-cell prims come from Pass 1) ...
      cls_h.cons2prims(bxg, state4, prims4, bx);

      // ... then zero the RHS on the ghost ring + shell (the interior was
      // zeroed in Pass 1; parity with the legacy bxg-wide zero)
      {
        const Box skipz = ibx;
        const bool skipz_ok = have_ibx;
        ParallelFor(bxg, ncons,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                      if (skipz_ok && skipz.contains(i, j, k)) { return; }
                      state4(i, j, k, n) = 0.0;
                    });
      }

      // fluxes + divergence for the shell in ONE masked call over the
      // tilebox: faces and cells interior to ibx are skipped (Pass 1 did
      // them); seam faces are recomputed from the same prims -> identical
      region_flux_div(prob_rhs, geom, bx, prims4, state4, cls_d, ncons, ibx);

      // NSCBC hook and source terms once per fab, after the divergence is
      // complete on the whole tilebox (same ordering as the legacy path)
      try_rhs_nscbc(0, geom, mfi, prims4, state4, cls_d, pparm_d, dt,
                    cur_time);
      prob_rhs.src(geom, mfi, prims4, state4, cls_d, dt, cur_time);
    }
  }
#endif  // !(AMREX_USE_GPIBM || CNS_USE_EB)
}
