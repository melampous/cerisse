#include <AMReX_FillPatchUtil.H>
#include <CNS.h>
#include <prob.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

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
[[maybe_unused]] inline MultiFab
ovl_make_mf_crse_patch(FabArrayBase::FPinfo const& fpc, int ncomp) {
  return MultiFab(fpc.ba_crse_patch, fpc.dm_patch, ncomp, 0, MFInfo(),
                  *fpc.fact_crse_patch);
}
[[maybe_unused]] inline MultiFab
ovl_make_mf_fine_patch(FabArrayBase::FPinfo const& fpc, int ncomp) {
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

  if (fr_as_crse != nullptr || fr_as_fine != nullptr) {
    amrex::Abort(
        "cns.overlap_comm=1 cannot be combined with cns.do_reflux=1: "
        "the split face-flux path does not expose a unique level flux MultiFab");
  }

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
