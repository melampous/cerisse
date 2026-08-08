#include <AMReX_FluxRegister.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_Reduce.H>
#include <CNS.h>
#include <CNSconstants.h>
#include <prob.h>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <type_traits>

#ifdef AMREX_USE_GPIBM
#include <ibm_positivity_limiter.h>
#include <ibm_solver.h>
#endif
using namespace amrex;

Real CNS::advance(Real time, Real dt, int iteration, int ncycle) {
  BL_PROFILE("CNS::advance()");

#ifdef AMREX_USE_GPIBM
  if (plot_surf && !ib_move && !geom.IsRZ() && level == 0) {
    if (!IBM::ib.conservativeSurfaceLoadWindowActive()) {
      IBM::ib.beginConservativeSurfaceLoadWindow(time);
    }
  }
  // Checkpoints do not serialize the transient RK-integrated IBM budget.
  // Initialize it from the synchronized accepted hierarchy before the first
  // restart advance swaps StateData time levels.  Deferring this to the first
  // CSV sample would discard the intervening crossing-flux impulse.
  if (ib_momentum_budget && level == 0 &&
      !ib_momentum_budget_initialized) {
    for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
      dynamic_cast<CNS&>(parent->getLevel(lev)).writeSurfFile(true);
    }
    initializeIBMMomentumBudget(time);
  }
#endif

  state[0].allocOldData();
  state[0].swapTimeLevels(dt);
  
  MultiFab& S1 = get_old_data(State_Type);
  MultiFab& S2 = get_new_data(State_Type);

  int ncons = d_prob_closures->NCONS;
  int nghost= d_prob_closures->NGHOST;
  MultiFab Stemp(grids,dmap,ncons,nghost,MFInfo(),Factory());
  std::unique_ptr<MultiFab> admissibility_fe_candidate;

#ifdef AMREX_USE_GPIBM
  std::unique_ptr<MultiFab> admissibility_stage_input;
  std::unique_ptr<MultiFab> admissibility_low_rhs;
  std::unique_ptr<MultiFab> admissibility_low_candidate;
  std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM>
      admissibility_high_face_flux;
  std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM>
      admissibility_low_face_flux;
  const bool use_production_positivity_limiter =
      CNS::ibm_positivity_flux_limiter;
  if (CNS::ibm_positivity_retry) {
    amrex::Abort(
        "cns.ibm_positivity_retry is not yet certified for the pure full-cell "
        "GP-IBM transaction path; disable it");
  }
  if (use_production_positivity_limiter) {
    const bool use_stage_local_gc_nscbc =
        use_nscbc &&
        nscbc_parm.ghost_update_model ==
            nscbc::GHOST_UPDATE_STAGE_LOCAL_GC;
    if (level != 0 || parent->finestLevel() != 0 || do_reflux != 0) {
      amrex::Abort(
          "IBM positivity limiting currently requires one uniform level "
          "with reflux disabled");
    }
    if (geom.IsRZ() || CNS::ib_move ||
        (use_nscbc && !use_stage_local_gc_nscbc)) {
      amrex::Abort(
          "IBM positivity limiting requires fixed Cartesian IBM; the only "
          "supported NSCBC coupling is the stage-local GC-NSCBC outflow");
    }
    if (geom.isAnyPeriodic() && !use_stage_local_gc_nscbc) {
      amrex::Abort(
          "IBM positivity limiting permits a periodic direction only for the "
          "restricted stage-local GC-NSCBC coupling");
    }
    constexpr bool supported_rhs =
        std::is_base_of_v<
            weno_t<ReconScheme::WenoZ5, PROB::ProbClosures>,
            PROB::ProbRHS> &&
        std::is_base_of_v<no_diffusive_t, PROB::ProbRHS> &&
        std::is_base_of_v<no_source_t, PROB::ProbRHS>;
    if constexpr (!supported_rhs) {
      amrex::Abort(
          "cns.ibm_positivity_flux_limiter=1 currently requires "
          "source-free LLF-WENO-Z5 Euler");
    }
    if constexpr (PROB::ProbClosures::NCONS != 5) {
      amrex::Abort(
          "cns.ibm_positivity_flux_limiter=1 currently requires the "
          "single-species five-equation ideal-gas system");
    }
    if constexpr (!PROB::ProbIB::stationary_slip_wall) {
      amrex::Abort(
          "cns.ibm_positivity_flux_limiter=1 currently requires a fixed "
          "Euler-slip wall");
    }
    admissibility_stage_input = std::make_unique<MultiFab>(
        grids, dmap, ncons, nghost, MFInfo(), Factory());
    admissibility_low_rhs = std::make_unique<MultiFab>(
        grids, dmap, ncons, 0, MFInfo(), Factory());
    admissibility_low_candidate = std::make_unique<MultiFab>(
        grids, dmap, ncons, 0, MFInfo(), Factory());
    admissibility_fe_candidate = std::make_unique<MultiFab>(
        grids, dmap, ncons, 0, MFInfo(), Factory());
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      BoxArray face_boxes =
          amrex::convert(grids, IntVect::TheDimensionVector(dir));
      admissibility_high_face_flux[dir] = std::make_unique<MultiFab>(
          face_boxes, dmap, ncons, 0, MFInfo().SetArena(The_Async_Arena()));
      admissibility_low_face_flux[dir] = std::make_unique<MultiFab>(
          face_boxes, dmap, ncons, 0, MFInfo().SetArena(The_Async_Arena()));
    }
  }
#endif

#ifdef AMREX_USE_GPIBM
  // Moving geometry must follow the RK abscissa, not merely the end-of-step
  // time.  Each AMR level receives its own (time,dt) from AMReX, so this also
  // makes fine-level subcycles use their actual substep times.
#ifdef CNS_USE_FSI
  std::unique_ptr<FabArray<BaseFab<uint8_t>>> old_markers;
  Real prepared_ibm_time = std::numeric_limits<Real>::quiet_NaN();
#endif
  if (CNS::ib_move) {
#ifdef CNS_USE_FSI
    auto& mfab = *IBM::ib.bmf_a[level];
    old_markers = std::make_unique<FabArray<BaseFab<uint8_t>>>(
        mfab.boxArray(), mfab.DistributionMap(), 1, mfab.nGrow(),
        MFInfo().SetArena(The_Async_Arena()));
#else
    amrex::Abort("ib.move=1 requires a USE_FSI=TRUE build");
#endif
  }

  auto prepare_ibm_stage = [&](Real stage_time, MultiFab& stage_state,
                               bool rebuild_surface) {
    if (!CNS::ib_move) return;

#ifdef CNS_USE_FSI
    const bool have_prepared_time = std::isfinite(prepared_ibm_time);
    const Real scale = have_prepared_time
        ? amrex::max(Real(1.0), amrex::max(std::abs(stage_time),
                                           std::abs(prepared_ibm_time)))
        : Real(1.0);
    const bool same_time = have_prepared_time
        && std::abs(stage_time - prepared_ibm_time)
               <= Real(64.0) * std::numeric_limits<Real>::epsilon() * scale;

    if (!same_time) {
      // Snapshot this level's previous topology before replacing its marker MF.
      auto& mfab_pre = *IBM::ib.bmf_a[level];
      for (MFIter mfi(mfab_pre, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.fabbox();
        auto const& dst = old_markers->array(mfi);
        auto const& src = mfab_pre.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          dst(i,j,k,0) = src(i,j,k,0);
        });
      }
      // rebuildIBM destroys mfab_pre; make the queued snapshot complete first.
      Gpu::streamSynchronize();

#ifdef CNS_FSI_DEFORMABLE
      PROB::update_geometry(stage_time, IBM::ib.geom_a, IBM::ib.ngeom);
      IBM::ib.rebuildGeometryData();
#else
      PROB::update_rigid_transforms(stage_time, IBM::ib.transform_a,
                                    IBM::ib.ngeom);
      for (int i = 0; i < IBM::ib.ngeom; ++i) {
        IBM::ib.updateRigidTransform(i, IBM::ib.transform_a[i]);
      }
#endif

      // Surface ownership/interpolation is only needed at the final synchronized
      // subcycle.  Marker and GP geometry, however, is required at every stage.
      rebuildIBM(false);
      IBM::ib.fixExposedCells(*old_markers, stage_state, level);
      prepared_ibm_time = stage_time;
    }

    PROB::Motion::sim_time = stage_time;

    if (rebuild_surface && level == parent->finestLevel()) {
      for (int lev = parent->finestLevel(); lev >= 0; --lev) {
        IBM::ib.computeSurfIndices(lev);
      }
    }
#else
    amrex::ignore_unused(stage_time, stage_state, rebuild_surface);
#endif
  };
#else
  auto prepare_ibm_stage = [](Real, MultiFab&, bool) {};
#endif

  FluxRegister* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    CNS& fine_level = getLevel(level + 1);
    if (!fine_level.flux_reg) {
      amrex::Abort(
          "advance: missing fine-level FluxRegister; restart/regrid did not "
          "rebuild the conservative AMR register");
    }
    fr_as_crse = fine_level.flux_reg.get();
  }

  FluxRegister* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    if (!flux_reg) {
      amrex::Abort(
          "advance: missing level FluxRegister; restart/regrid did not "
          "rebuild the conservative AMR register");
    }
    fr_as_fine = flux_reg.get();
  }

  if (fr_as_crse) {
    fr_as_crse->setVal(Real(0.0));
  }

  // Persistent NSCBC ghost cells are an ODE state coupled to the interior
  // solution.  Integrate them with exactly the same RK coefficients used for
  // S2; nscbc_ghost_state stores only the accepted end-of-step state.
  std::unique_ptr<MultiFab> G_old;
  std::unique_ptr<MultiFab> G_stage;
  std::unique_ptr<MultiFab> G_rhs;
  const bool use_persistent_nscbc =
    use_nscbc && nscbc_parm.ghost_update_model ==
      nscbc::GHOST_UPDATE_PERSISTENT_ODE;
  if (use_persistent_nscbc) {
    if (!nscbc_ghost_state) {
      amrex::Abort("NSCBC enabled but persistent ghost state is not allocated");
    }
    if (!nscbc_ghost_initialized) {
      initialize_nscbc_ghost_state(time);
    }
    G_old = std::make_unique<MultiFab>(
      grids, dmap, ncons, nghost, MFInfo(), Factory());
    G_stage = std::make_unique<MultiFab>(
      grids, dmap, ncons, nghost, MFInfo(), Factory());
    G_rhs = std::make_unique<MultiFab>(
      grids, dmap, ncons, nghost, MFInfo(), Factory());
    MultiFab::Copy(*G_old, *nscbc_ghost_state, 0, 0, ncons, nghost);
    MultiFab::Copy(*G_stage, *nscbc_ghost_state, 0, 0, ncons, nghost);
  }

  auto prepare_nscbc_stage = [&](MultiFab& stage_state) {
    if (!use_nscbc) return;
    if (use_persistent_nscbc) {
      copy_nscbc_ghost_to_state(stage_state, *G_stage);
      compute_nscbc_ghost_rhs(stage_state, *G_rhs);
    } else {
      fill_nscbc_gc_ghosts(stage_state);
    }
  };

  auto ghost_saxpy = [&](Real a) {
    if (use_persistent_nscbc) {
      MultiFab::Saxpy(*G_stage, a, *G_rhs, 0, 0, ncons, nghost);
    }
  };

#ifdef AMREX_USE_GPIBM
  auto capture_admissibility_stage_input = [&](MultiFab& stage_input) {
    if (!use_production_positivity_limiter) return;
#if (AMREX_SPACEDIM < 3)
    stage_input.setVal(
        Real(0.0), PROB::ProbClosures::UMZ, 1, stage_input.nGrow());
#endif
    MultiFab::Copy(*admissibility_stage_input, stage_input, 0, 0, ncons,
                   nghost);
  };

  auto captured_high_order_face_flux = [&]() {
    std::array<MultiFab*, AMREX_SPACEDIM> result{};
    const bool capture_flux = use_production_positivity_limiter;
    if (capture_flux) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = admissibility_high_face_flux[dir].get();
      }
    }
    return result;
  };

  auto accumulate_final_limited_flux = [&](Real stage_weight) {
    if (!use_production_positivity_limiter) return;
    if (!(stage_weight > Real(0.0)) || !std::isfinite(stage_weight)) {
      amrex::Abort(
          "post-limiter conservation accounting requires a positive finite "
          "SSPRK quadrature weight");
    }
    std::array<MultiFab*, AMREX_SPACEDIM> final_face_flux{};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      final_face_flux[dir] = admissibility_high_face_flux[dir].get();
    }
    if (ib_momentum_budget) {
      accumulateIBMBudgetFlux(final_face_flux, stage_weight);
    }
    if (!plot_surf) return;
    if (!IBM::ib.conservativeSurfaceLoadWindowActive()) {
      amrex::Abort(
          "post-limiter conservative IBM surface-load window is not active");
    }
    // This restricted production path is inviscid, so the finalized Euler and
    // combined update fluxes are identical. Record both channels to preserve
    // the existing ConservativeElementForce output contract.
    IBM::ib.accumulateConservativeSurfaceLoad(
        level, geom, final_face_flux, stage_weight, true);
    IBM::ib.accumulateConservativeSurfaceLoad(
        level, geom, final_face_flux, stage_weight, false);
  };

  auto apply_production_fe_limiter =
      [&](MultiFab& high_candidate, Real forward_euler_dt,
          const char* scheme, int stage_index, int stage_count,
          Real candidate_time) -> bool {
        if (!use_production_positivity_limiter) return true;

        const auto& closure = *CNS::h_prob_closures;
        const auto high_stats =
            IBM::positivity::collectStageAdmissibility(
                high_candidate, level, closure);
        if (high_stats.admissible()) {
          if (CNS::ibm_positivity_flux_limiter_verbose) {
            amrex::Print()
                << "[IBM-Positivity-Limiter] scheme=" << scheme
                << " FE_bracket=" << stage_index << '/' << stage_count
                << " t=" << std::setprecision(17) << candidate_time
                << " mode=HIGH_ORDER theta_min=1\n";
          }
          return true;
        }

        std::array<MultiFab*, AMREX_SPACEDIM> high_face_flux{};
        std::array<MultiFab*, AMREX_SPACEDIM> low_face_flux{};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          high_face_flux[dir] = admissibility_high_face_flux[dir].get();
          low_face_flux[dir] = admissibility_low_face_flux[dir].get();
        }
        const Real low_cfl =
            IBM::positivity::computePiecewiseConstantRusanovRhs<
                PROB::ProbRHS>(
                *admissibility_stage_input, *admissibility_low_rhs, geom,
                level, CNS::d_prob_closures, closure, forward_euler_dt,
                low_face_flux);
        MultiFab::Copy(*admissibility_low_candidate,
                       *admissibility_stage_input, 0, 0, ncons, 0);
        MultiFab::Saxpy(*admissibility_low_candidate, forward_euler_dt,
                        *admissibility_low_rhs, 0, 0, ncons, 0);

        const auto low_stats =
            IBM::positivity::collectStageAdmissibility(
                *admissibility_low_candidate, level, closure);
        if (!low_stats.admissible()) {
          amrex::Print()
              << "\n[IBM-Positivity-Limiter] scheme=" << scheme
              << " FE_bracket=" << stage_index << '/' << stage_count
              << " t=" << std::setprecision(17) << candidate_time
              << " dt_FE=" << forward_euler_dt
              << " multidimensional_cfl=" << low_cfl
              << " verdict=ALL_LOW_ORDER_INADMISSIBLE\n";
          IBM::positivity::printStageAdmissibility(
              "high-order", high_stats);
          IBM::positivity::printStageAdmissibility(
              "all-low-order", low_stats);
          if (CNS::ibm_positivity_retry) {
            return false;
          }
          amrex::Abort(
              "IBM positivity limiter found an inadmissible all-low update; "
              "enable cns.ibm_positivity_retry for transactional SSPRK43 "
              "time-interval subdivision (the step was terminated without "
              "clipping)");
        }


        const auto local =
            IBM::positivity::applyLocalSharedFaceConvexLimiter(
                high_candidate, *admissibility_low_candidate,
                high_face_flux, low_face_flux, geom, level, closure, scheme,
                stage_index, stage_count, candidate_time, forward_euler_dt,
                low_cfl, false, true);
        const char* fallback_mode = "LOCAL_SHARED_FACE";
        Real fallback_theta = local.minimum_theta;
        if (!local.limited_state_admissible) {
          const bool global_ok =
              IBM::positivity::tryGlobalThetaFluxFallback(
                  high_candidate, *admissibility_low_candidate,
                  high_face_flux, low_face_flux, geom, level, closure,
                  fallback_theta);
          if (global_ok) {
            fallback_mode = "GLOBAL_THETA";
          } else {
            IBM::positivity::replaceWithAllLowOrderFluxUpdate<
                PROB::ProbClosures>(
                high_candidate, *admissibility_low_candidate,
                high_face_flux, low_face_flux, geom);
            fallback_mode = "ALL_LOW_ORDER";
            fallback_theta = Real(0.0);
          }
        }

        const auto final_stats =
            IBM::positivity::collectStageAdmissibility(
                high_candidate, level, closure);
        amrex::Print()
            << "[IBM-Positivity-Limiter] scheme=" << scheme
            << " FE_bracket=" << stage_index << '/' << stage_count
            << " t=" << std::setprecision(17) << candidate_time
            << " dt_FE=" << forward_euler_dt
            << " multidimensional_cfl=" << low_cfl
            << " mode=" << fallback_mode
            << " theta_min=" << fallback_theta
            << " limited_faces=" << local.limited_faces
            << " limited_crossing_faces=" << local.limited_crossing_faces
            << " conservation_relative_linf="
            << local.conservation_relative_linf
            << " verdict="
            << (final_stats.admissible() ? "PASS" : "FAIL") << '\n';
        if (!final_stats.admissible()) {
          amrex::Abort(
              "IBM positivity limiter exhausted local, global-theta, and "
              "all-low-order fallbacks without an admissible FE bracket");
        }
        return true;
      };

  auto audit_production_stage =
      [&](const MultiFab& stage_state, const char* scheme, int stage_index,
          int stage_count, Real stage_time) {
        if (!use_production_positivity_limiter) return;
        const auto stats =
            IBM::positivity::collectStageAdmissibility(
                stage_state, level, *CNS::h_prob_closures);
        if (!stats.admissible()) {
          IBM::positivity::printStageAdmissibility(
              "post-SSP-convex-combination", stats);
          amrex::Abort(
              "IBM positivity limiter produced an inadmissible SSPRK state");
        }
        if (CNS::ibm_positivity_flux_limiter_verbose) {
          amrex::Print()
              << "[IBM-Positivity-Limiter] scheme=" << scheme
              << " stage=" << stage_index << '/' << stage_count
              << " t=" << std::setprecision(17) << stage_time
              << " post_SSP_audit=PASS\n";
        }
      };

#else
  auto capture_admissibility_stage_input = [](MultiFab&) {};
  auto captured_high_order_face_flux = []() {
    return std::array<MultiFab*, AMREX_SPACEDIM>{};
  };
  auto accumulate_final_limited_flux = [](Real) {};
  auto apply_production_fe_limiter =
      [](MultiFab&, Real, const char*, int, int, Real) { return true; };
  auto audit_production_stage =
      [](const MultiFab&, const char*, int, int, Real) {};
#endif

  if (order_rk == -2) {
    // Original time integration ///////////////////////////////////////////////
    // RK2 stage 1
    prepare_ibm_stage(time, S1, false);
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    prepare_nscbc_stage(Stemp);
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine, time,
                Real(0.5) * dt);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);
    ghost_saxpy(dt);
    // RK2 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n
    state[0].setNewTimeLevel(time + dt);
    prepare_ibm_stage(time + dt, S2, false);
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
    prepare_nscbc_stage(Stemp);
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine, time + dt,
                Real(0.5) * dt);
    // S_new = 0.5*(Sborder+S_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(S2, Real(0.5), S1, 0, Real(0.5), S2, 0, 0, ncons, 0);
    // S_new += 0.5*dt*dSdt
    MultiFab::Saxpy(S2, Real(0.5) * dt, Stemp, 0, 0, ncons, 0);
    if (use_persistent_nscbc) {
      MultiFab::LinComb(*G_stage, Real(0.5), *G_old, 0,
                        Real(0.5), *G_stage, 0, 0, ncons, nghost);
      ghost_saxpy(Real(0.5) * dt);
    }
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*


    ////////////////////////////////////////////////////////////////////////////
  } else if (order_rk == 0) {  // returns rhs
    if (CNS::ib_move) {
      amrex::Abort("order_rk=0 RHS-output mode does not support moving IBM");
    }
    prepare_ibm_stage(time, S1, false);
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    prepare_nscbc_stage(Stemp);
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine, time, Real(0.0));
    MultiFab::Copy(S2, Stemp, 0, 0, ncons, 0);
  } else if (order_rk == 1) {
    prepare_ibm_stage(time, S1, false);
    FillPatch(*this, Stemp, nghost, time, State_Type, 0,
              ncons);  // filled at t_n to evalulate f(t_n,y_n).
    prepare_nscbc_stage(Stemp);
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine, time, dt);
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);
    ghost_saxpy(dt);
  } else if (order_rk == 2) {
    // Low-storage SSP-RK(m,2): m stages, C=m-1, C_eff=1-1/m.
    // Ref: Gottlieb et al., "Strong Stability Preserving Runge-Kutta and
    // Multistep Time Discretizations", §4.2.
    int m = stages_rk;
    MultiFab::Copy(S2, S1, 0, 0, ncons, 0);
    state[0].setOldTimeLevel(time);
    state[0].setNewTimeLevel(time);
    // First m-1 forward-Euler increments
    for (int i = 1; i <= m - 1; i++) {
      const Real stage_time = time + dt * Real(i - 1) / (m - 1);
      prepare_ibm_stage(stage_time, S2, false);
      FillPatch(*this, Stemp, nghost, stage_time, State_Type, 0, ncons);
      prepare_nscbc_stage(Stemp);
      compute_rhs(Stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine,
                  stage_time, dt / Real(m));
      MultiFab::Saxpy(S2, dt / Real(m - 1), Stemp, 0, 0, ncons, 0);
      ghost_saxpy(dt / Real(m - 1));
      state[State_Type].setNewTimeLevel(
          time + dt * Real(i) /
                     (m - 1));  // important to do this for correct fillpatch
                                // interpolations for the proceeding stages
    }
    // final stage
    prepare_ibm_stage(time + dt, S2, false);
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
    prepare_nscbc_stage(Stemp);
    compute_rhs(Stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine, time + dt,
                dt / Real(m));
    MultiFab::LinComb(S2, Real(m - 1), S2, 0, dt, Stemp, 0, 0, ncons, 0);
    MultiFab::LinComb(S2, Real(1.0) / m, S1, 0, Real(1.0) / m, S2, 0, 0, ncons,
                      0);
    if (use_persistent_nscbc) {
      MultiFab::LinComb(*G_stage, Real(m - 1), *G_stage, 0,
                        dt, *G_rhs, 0, 0, ncons, nghost);
      MultiFab::LinComb(*G_stage, Real(1.0) / m, *G_old, 0,
                        Real(1.0) / m, *G_stage, 0, 0, ncons, nghost);
    }

    state[State_Type].setNewTimeLevel(time + dt);
  }

  // default
  else if (order_rk == 3) {
    if (stages_rk == 3) {
      // SSP-RK(3,3): http://ketch.github.io/numipedia/methods/SSPRK33.html

      // cns.overlap_comm=1: replace each (FillPatch + compute_rhs) pair by
      // compute_rhs_overlap, which overlaps the same-level ghost exchange
      // with the interior flux computation. The setOld/NewTimeLevel calls
      // are kept EXACTLY as in the legacy path (other machinery depends on
      // them). Both RHS paths receive the RK stage time explicitly.
#if defined(AMREX_USE_GPIBM) || defined(CNS_USE_EB)
      if (CNS::overlap_comm != 0) {
        amrex::Abort("cns.overlap_comm=1 is not supported in IBM/EB builds");
      }
      const bool use_ovl = false;
#else
      const bool use_ovl = (CNS::overlap_comm != 0);
#endif
      // level-wide prims workspace for the overlap path (allocated once per
      // advance, reused across the three stages)
      MultiFab prims_ovl;
      if (use_ovl) {
        prims_ovl.define(grids, dmap, PROB::ProbClosures::NPRIM, nghost,
                         MFInfo(), Factory());
      }

      state[0].setOldTimeLevel(time);
      prepare_ibm_stage(time, S1, false);
      if (use_ovl) {
        // stage-1 valid data = S1 (old), fill time = time
        compute_rhs_overlap(Stemp, dt, fr_as_crse, fr_as_fine, time, S1,
                            prims_ovl);
      } else {
        FillPatch(*this, Stemp, nghost, time, State_Type, 0,
                  ncons);  // filled at t_n to evalulate f(t_n,y_n).
        prepare_nscbc_stage(Stemp);
        compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine, time, dt / Real(6.0));
      }
      MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);
      ghost_saxpy(dt);

      state[0].setNewTimeLevel(
          time + dt);  // same time as upcoming FillPatch ensures we copy S2 to
                       // Sborder, without time interpolation
      prepare_ibm_stage(time + dt, S2, false);
      if (use_ovl) {
        // stage-2 valid data = S2 (new), fill time = time + dt
        compute_rhs_overlap(Stemp, dt / 4, fr_as_crse, fr_as_fine, time + dt,
                            S2, prims_ovl);
      } else {
        FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
        prepare_nscbc_stage(Stemp);
        compute_rhs(Stemp, dt / 4, fr_as_crse, fr_as_fine, time + dt,
                    dt / Real(6.0));
      }
      MultiFab::Xpay(Stemp, dt, S2, 0, 0, ncons, 0);
      MultiFab::LinComb(S2, Real(3.0) / 4, S1, 0, Real(1.0) / 4, Stemp, 0, 0,
                        ncons, 0);
      if (use_persistent_nscbc) {
        ghost_saxpy(dt);
        MultiFab::LinComb(*G_stage, Real(3.0) / 4, *G_old, 0,
                          Real(1.0) / 4, *G_stage, 0, 0, ncons, nghost);
      }

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      prepare_ibm_stage(time + dt / 2, S2, false);
      if (use_ovl) {
        // stage-3 valid data = S2 (new), fill time = time + dt/2
        compute_rhs_overlap(Stemp, dt * Real(2.0) / 3, fr_as_crse, fr_as_fine,
                            time + dt / 2, S2, prims_ovl);
      } else {
        FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
        prepare_nscbc_stage(Stemp);
        compute_rhs(Stemp, dt * Real(2.0) / 3, fr_as_crse, fr_as_fine,
                    time + dt / 2, dt * Real(2.0) / Real(3.0));
      }
      MultiFab::Xpay(Stemp, dt, S2, 0, 0, ncons, 0);
      MultiFab::LinComb(S2, Real(1.0) / 3, S1, 0, Real(2.0) / 3, Stemp, 0, 0,
                        ncons, 0);
      if (use_persistent_nscbc) {
        ghost_saxpy(dt);
        MultiFab::LinComb(*G_stage, Real(1.0) / 3, *G_old, 0,
                          Real(2.0) / 3, *G_stage, 0, 0, ncons, nghost);
      }

      state[State_Type].setNewTimeLevel(
          time + dt);  // important to do this for correct fillpatch
                       // interpolations for the proceeding stages
    }

    else if (stages_rk == 4) {
      if (CNS::ibm_positivity_retry) {
        amrex::Abort(
            "cns.ibm_positivity_retry is unavailable in the pure full-cell "
            "GP-IBM build until its state, NSCBC ghost, load, and budget "
            "rollback are certified");
      } else {
      // SSP-RK(4,3): http://ketch.github.io/numipedia/methods/SSPRK43.html
      // Ref: Gottlieb et al., §4.2, p. 85.

      state[0].setOldTimeLevel(time);
      prepare_ibm_stage(time, S1, false);
      FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
      prepare_nscbc_stage(Stemp);
      capture_admissibility_stage_input(Stemp);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine, time,
                  dt / Real(6.0), captured_high_order_face_flux());
      MultiFab::LinComb(S2, Real(1.0), S1, 0, dt / 2, Stemp, 0, 0, ncons, 0);
      ghost_saxpy(dt / 2);
      apply_production_fe_limiter(
          S2, dt / Real(2.0), "ssprk43", 1, 4,
          time + dt / Real(2.0));
      accumulate_final_limited_flux(dt / Real(6.0));
      audit_production_stage(
          S2, "ssprk43", 1, 4, time + dt / Real(2.0));

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      prepare_ibm_stage(time + dt / 2, S2, false);
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      prepare_nscbc_stage(Stemp);
      capture_admissibility_stage_input(Stemp);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine, time + dt / 2,
                  dt / Real(6.0), captured_high_order_face_flux());
      MultiFab::Saxpy(S2, dt / 2, Stemp, 0, 0, ncons, 0);
      ghost_saxpy(dt / 2);
      apply_production_fe_limiter(
          S2, dt / Real(2.0), "ssprk43", 2, 4, time + dt);
      accumulate_final_limited_flux(dt / Real(6.0));
      audit_production_stage(S2, "ssprk43", 2, 4, time + dt);

      state[0].setNewTimeLevel(
          time + dt);  // same time as upcoming FillPatch ensures we copy S2 to
                       // Sborder, without time interpolation
      prepare_ibm_stage(time + dt, S2, false);
      FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
      prepare_nscbc_stage(Stemp);
      capture_admissibility_stage_input(Stemp);
      if (CNS::ibm_positivity_flux_limiter) {
        // SSPRK(4,3) stage 3 is
        //   2/3 U^n + 1/3 [U^(2) + dt/2 L(U^(2))].
        // Limit the bracketed Forward-Euler state before taking the convex
        // combination. Limiting the already-combined stage would not establish
        // the SSP admissibility argument.
        compute_rhs(Stemp, dt / Real(2.0), fr_as_crse, fr_as_fine,
                    time + dt, dt / Real(6.0),
                    captured_high_order_face_flux());
        MultiFab::LinComb(*admissibility_fe_candidate, Real(1.0), S2, 0,
                          dt / Real(2.0), Stemp, 0, 0, ncons, 0);
        apply_production_fe_limiter(
            *admissibility_fe_candidate, dt / Real(2.0), "ssprk43", 3, 4,
            time + dt / Real(2.0));
        accumulate_final_limited_flux(dt / Real(6.0));
        MultiFab::LinComb(S2, Real(2.0) / Real(3.0), S1, 0,
                          Real(1.0) / Real(3.0),
                          *admissibility_fe_candidate, 0, 0, ncons, 0);
        audit_production_stage(
            S2, "ssprk43", 3, 4, time + dt / Real(2.0));
      } else {
        // Preserve the historical operation order bit for bit when the
        // production limiter is disabled.
        compute_rhs(Stemp, dt / Real(6.0), fr_as_crse, fr_as_fine,
                    time + dt, dt / Real(6.0),
                    captured_high_order_face_flux());
        MultiFab::LinComb(S2, Real(2.0) / 3, S1, 0,
                          Real(1.0) / 3, S2, 0, 0, ncons, 0);
        MultiFab::Saxpy(S2, dt / Real(6.0), Stemp, 0, 0, ncons, 0);
        if (use_persistent_nscbc) {
          MultiFab::LinComb(*G_stage, Real(2.0) / 3, *G_old, 0,
                            Real(1.0) / 3, *G_stage, 0, 0, ncons, nghost);
          ghost_saxpy(dt / Real(6.0));
        }
      }

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      prepare_ibm_stage(time + dt / 2, S2, false);
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      prepare_nscbc_stage(Stemp);
      capture_admissibility_stage_input(Stemp);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine, time + dt / 2,
                  dt / Real(2.0), captured_high_order_face_flux());
      MultiFab::Saxpy(S2, dt / 2, Stemp, 0, 0, ncons, 0);
      ghost_saxpy(dt / 2);
      apply_production_fe_limiter(
          S2, dt / Real(2.0), "ssprk43", 4, 4, time + dt);
      accumulate_final_limited_flux(dt / Real(2.0));
      audit_production_stage(S2, "ssprk43", 4, 4, time + dt);

      state[State_Type].setNewTimeLevel(
          time + dt);  // important to do this for correct fillpatch
                       // interpolations for the proceeding stages
      }
    }

    else {
      // General SSP-RK(n^2, 3), n>2: C=2, C_eff=1-1/n (not yet implemented).
      // Ref: Gottlieb et al., §4.2, p. 85.
      amrex::Abort("SSPRK(n^2,3) with n>2 is not yet implemented");
    }

  }


  if (use_persistent_nscbc) {
    MultiFab::Copy(*nscbc_ghost_state, *G_stage, 0, 0, ncons, nghost);
    nscbc_ghost_initialized = true;
  }

  // Leave markers, GP geometry, wall clock and surface ownership at the
  // accepted state time.  On subcycled AMR levels, surface indices are only
  // assembled after the final subcycle reaches the synchronized coarse time.
  prepare_ibm_stage(time + dt, S2, iteration >= ncycle);

#if ENSURE_MASSFRACSUM_ONE  
  clip_species_state(S2);
#endif

#ifdef AMREX_USE_GPIBM
  // ==========================================================================
  // End-of-step IBM correction (runs on S2 = state at t^{n+1})
  //
  // After the RK stages, the conservative state in IBM cells needs cleanup:
  //   [1] Ghost points: overwrite with wall-BC-reconstructed primitives
  //                     (always — also needed for static geometry).
  //   [2] Interior solid cells (FSI only): flood-fill from fluid/ghost
  //                     neighbors, then zero momentum. This keeps solid
  //                     cells carrying bounded, physically plausible data
  //                     so that AMR FillPatch/avgDown during regrid does
  //                     not interpolate garbage into fresh fluid cells.
  //
  // NO SAFETY NET: if any fluid cell becomes NaN/Inf/non-positive, that is
  // a numerical failure of the scheme (under-resolved shocks, wrong CFL,
  // missing positivity limiter, etc.), not something to silently patch.
  // We detect such cells and abort with a diagnostic — upstream policy.
  // ==========================================================================
  {
    const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
    const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
    auto& ib_mf = *IBM::ib.bmf_a[level];

    // ------------------------------------------------------------------------
    // Reconstruct ghost-point primitives from the current flow state + wall BC
    // ------------------------------------------------------------------------
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);

    // In 2D, prob_initdata / bcnormal may never write UMZ, so FillPatch can
    // hand back garbage in the z-momentum plane (valid cells at step 0 and
    // physical-BC ghosts every step). Zero it before cons2prims so Pass 1
    // does not embed garbage-derived GP values into S2. Mirrors the
    // per-stage sanitisation in compute_rhs.cpp.
#if (AMREX_SPACEDIM < 3)
    Stemp.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, Stemp.nGrow());
#endif

    MultiFab prims_mf(
                      Stemp.boxArray(), Stemp.DistributionMap(),
                      cls_h.NPRIM,
                      IBM::ib.volumeInterpolationNghost(level),
                      MFInfo().SetArena(The_Async_Arena()));
    for (MFIter mfi(Stemp, false); mfi.isValid(); ++mfi) {
      cls_h.cons2prims(mfi, Stemp.array(mfi), prims_mf.array(mfi));
    }

    // Sync wall-motion time so that compute_surfIB() sees the correct
    // instantaneous wall velocity for moving-wall BCs.
#ifdef CNS_USE_FSI
    PROB::Motion::sim_time = time + dt;
#endif

    const bool rz_annular_gp = IBM::ib.rzAnnularCellAverageEnabled(level);
    if (rz_annular_gp) {
      IBM::ib.computeAllGPsRZAnnular(prims_mf, Stemp, cls_d, level);
    } else {
      IBM::ib.computeAllGPs(prims_mf, cls_d, level);
    }
    // No streamSynchronize needed here: computeAllGPs launches its GP kernel
    // on the default AMReX gpuStream (its internal H2D pointer-table copy is
    // also stream-ordered, see ibm_solver.h), and Pass 1 below reads prims_mf
    // exclusively from device kernels on that same stream — stream order
    // already guarantees the GP writes are visible. No host code touches
    // prims_mf data in between (MFIter/array() are metadata-only).

    // ------------------------------------------------------------------------
    // Pass 1 (ALWAYS): Write GP-corrected primitives back to S2 as conservatives.
    //                   Only touches cells marked as ghost points (ibMarkers(,1)).
    // ------------------------------------------------------------------------
    for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& state = S2.array(mfi);
      Array4<Real> const& prims = prims_mf.array(mfi);
      Array4<Real const> const& annular_stage = Stemp.const_array(mfi);
      const auto& ibMarkers = ib_mf.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (ibMarkers(i, j, k, 1)) {
          if (rz_annular_gp) {
            for (int n = 0; n < PROB::ProbClosures::NCONS; ++n) {
              state(i, j, k, n) = annular_stage(i, j, k, n);
            }
          } else {
            IntVect iv(AMREX_D_DECL(i, j, k));
            Real cons[PROB::ProbClosures::NCONS];
            cls_h.prims2cons(iv, prims, cons);
            for (int n = 0; n < PROB::ProbClosures::NCONS; ++n) {
              state(i, j, k, n) = cons[n];
            }
          }
        }
      });
    }

    // ------------------------------------------------------------------------
    // Pass 2: Flood-fill interior solid cells + zero their momentum.
    //
    // Runs when:
    //   - ib_move = 1 (FSI): cells that were solid can become fluid as the
    //     body moves, so solid cells must carry bounded, physically plausible
    //     data that reflects something close to wall-BC conditions.
    //   - pass2_static = 1 (opt-in for static geometry): some shock–body
    //     interaction cases benefit from clean solid-cell data because
    //     WENO stencils reach across the surface and read solid values.
    //     Empirically, flood-filling helps simple geometries
    //     (2d_bvh_cpu, airfoil_static) but *destabilises* complex geometries
    //     (2DSphere, complex_geom) — the averaged post-shock state leaks
    //     into the body and the next step's WENO oscillates on the gradient.
    //     Hence opt-in, not default.
    //
    // Zeroing momentum prevents spurious velocity amplification when
    // averaging fluid neighbors from different acoustic phases.
    //
    // Sync S2 ghost cells from neighboring fabs first: the RK stages only
    // write the valid region of S2, so cross-fab ghost cells still hold
    // data from the previous step. The flood-fill reads 3×3 neighbors,
    // and at box boundaries those neighbors land in that stale ghost
    // region — without FillBoundary, Pass 2 averages current-step valid
    // cells with previous-step ghost values, silently polluting the solid
    // state and feeding garbage into the next step's WENO stencils.
    // ------------------------------------------------------------------------
    if (CNS::ib_move || CNS::pass2_static) {
    S2.FillBoundary(geom.periodicity());
    {
      // Fixed flood-fill iteration count, derived from the marker construction:
      //   - Seeds (tag 0) are fluid cells and the GP band. computeMarkers marks
      //     a solid cell as a ghost point iff fluid lies within +-ghost_layers
      //     (gl) along an axis, so seeds cover the first gl solid layers —
      //     unless initialiseGPs demoted a GP (imp_ninterp(0)==0 -> comp1=0).
      //   - Each iteration advances the fill front one cell (3^DIM stencil).
      //   - Downstream consumers (next step's stencils, FillPatch/avgDown)
      //     read at most NGHOST cells past the fluid interface, so
      //     gl + NGHOST iterations reach every consumer-visible cell even in
      //     the worst case where the entire GP band was demoted.
      // Deeper interior cells (thick bodies) were equally not guaranteed by
      // the previous 32-iteration cap: they keep bounded prior data and are
      // momentum-zeroed below. Running a fixed count removes the per-iteration
      // DeviceScalar alloc + streamSynchronize + 4 B readback (up to 32 per
      // FAB per step); convergence is verified ONCE per level after the loop.
      const int flood_iters = IBM::ib.ghost_layers + nghost;

      // Level-wide convergence counter: >0 after the loop means a
      // stencil-visible interior solid cell (see check below) was left
      // unfilled — abort loudly (no silent under-fill).
      Gpu::DeviceScalar<int> d_unconverged(0);
      int* p_unconv = d_unconverged.dataPtr();

      for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.tilebox();
        const Box& bxg = mfi.growntilebox(d_prob_closures->NGHOST);
        auto const& state     = S2.array(mfi);
        auto const& ibMarkers = ib_mf.array(mfi);
        const int nc = ncons;

        // Tag values: 0 = already valid (fluid or GP-corrected ghost point)
        //             1 = interior solid, needs fixing
        //             2 = solid, fixed in a previous iteration
        //
        // JACOBI DOUBLE BUFFER (F3 fix): the fill is a Jacobi iteration,
        // not Gauss-Seidel. tagfab has TWO components: comp iter%2 holds
        // the previous iterate's tags (read-only within an iteration),
        // comp 1-iter%2 receives the next iterate's; the buffers swap by
        // parity each iteration. The old in-place update let sibling GPU
        // threads in the SAME kernel launch observe (or miss) a neighbor's
        // just-set tag=2 — and then read its just-written state — depending
        // on warp scheduling: run-to-run nondeterministic fill values that
        // leak into fluid cells through next-step stencils (finding F3,
        // evidence: IBM/gpu_a100_test/fix135_ab/f3_flood_nondet/).
        // `state` itself needs no second buffer: within one iteration the
        // cells WRITTEN are exactly those with old-tag 1, while the cells
        // READ are those with old-tag 0 or 2 — disjoint sets, and each
        // writer touches only its own cell. Both backends execute this same
        // kernel; under Jacobi the CPU-sequential and GPU-parallel orders
        // compute identical fill values (cross-backend parity). NOTE: the
        // pre-fix sequential CPU sweep acted as Gauss-Seidel (neighbors
        // filled earlier in the same sweep fed the average), so fill values
        // for fresh solid cells change slightly vs the old CPU results;
        // this is an IC for newly-uncovered cells, physics-neutral — and
        // GPU results become deterministic. The one-cell-per-iteration
        // front advance assumed by the flood_iters derivation above is
        // exactly the Jacobi propagation rate, so the fixed count and the
        // single final convergence check are unchanged.
        //
        // Async arena (not Managed): all kernels below are queued without
        // any host sync before tagfab leaves scope, so its free must be
        // stream-ordered; the tag data is device-only anyway. One 2-comp
        // fab per FAB per step — no per-iteration allocations.
        BaseFab<int> tagfab(bxg, 2, The_Async_Arena());
        auto const& tag = tagfab.array();

        ParallelFor(bxg, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int t = (ibMarkers(i,j,k,0) == 0 || ibMarkers(i,j,k,1) != 0)
                        ? 0 : 1;
          // Init BOTH buffers: ghost-ring cells (bxg minus bx) are never
          // rewritten by the fill kernel (domain bx), so the two comps must
          // agree there from the start.
          tag(i,j,k,0) = t;
          tag(i,j,k,1) = t;
        });

        // Iterative flood-fill: each iteration propagates valid data one
        // cell deeper into the solid region. Fixed iteration count (see
        // derivation above); iterations past convergence are no-ops for
        // state (the front has stalled; tag comps just carry forward), so
        // no per-iteration convergence readback is needed.
        for (int iter = 0; iter < flood_iters; ++iter) {
          const int told = iter & 1;   // previous iterate (read)
          const int tnew = 1 - told;   // next iterate (write)
          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int t = tag(i,j,k,told);
            if (t != 1) {              // already valid or already fixed:
              tag(i,j,k,tnew) = t;     // carry tag into the next iterate
              return;
            }

            Real sum[PROB::ProbClosures::NCONS] = {};
            int count = 0;
            for (int dj = -1; dj <= 1; ++dj) {
              for (int di = -1; di <= 1; ++di) {
#if (AMREX_SPACEDIM == 3)
                for (int dk = -1; dk <= 1; ++dk) {
#else
                { int dk = 0;
#endif
                  if (di == 0 && dj == 0 && dk == 0) continue;
                  const int ii = i+di, jj = j+dj, kk = k+dk;
                  if (!bxg.contains(IntVect(AMREX_D_DECL(ii,jj,kk)))) continue;
                  const int tn = tag(ii,jj,kk,told);
                  if (tn == 0 || tn == 2) {
                    for (int n = 0; n < nc; ++n)
                      sum[n] += state(ii,jj,kk,n);
                    count++;
                  }
                }
              }
            }
            if (count > 0) {
              const Real inv = Real(1.0) / count;
              for (int n = 0; n < nc; ++n)
                state(i,j,k,n) = sum[n] * inv;
              tag(i,j,k,tnew) = 2;
            } else {
              tag(i,j,k,tnew) = 1;     // still unreached, try next iteration
            }
          });
          // No sync/readback: iterations chain on the stream; the tag comps
          // carry the fill front between kernels.
        }
        // Final tags live in the comp last written: flood_iters%2 (equals
        // the init comp 0 when flood_iters == 0).
        const int tfin = flood_iters & 1;

        // Convergence check, accumulated level-wide into d_unconverged.
        // Checked band: tag-1 cells with fluid within min(gl+1, NGHOST)
        // along an axis (clamped to NGHOST so marker reads stay inside bxg),
        // but ONLY when the straight path to that fluid stays inside the
        // valid box: the farthest intermediate cell (distance d-1) must be
        // in bx. Such a cell is provably fixable within check_band <=
        // flood_iters iterations — every in-box path cell either is a GP
        // seed or fills from the fluid seed by induction (fluid cells are
        // tag 0 anywhere in bxg via marker comp 0), so if it is still tag 1
        // here the fill kernel or the marker data is genuinely broken and
        // we abort rather than let stencils read the unfilled cell.
        // Paths that exit the valid box are NOT counted: ghost-ring cells
        // are never processed by the fill (kernel domain is bx) and their
        // GP mark (comp 1) is 0 in this FAB's private marker copy
        // (computeMarkers Step 2 writes comp 1 over the valid box only,
        // GP_BOX_EXTRA=0, and markers are never FillBoundary'd) — so
        // box-edge cells reachable only across the face are legitimately
        // unfillable and were silently left by the old per-iteration loop
        // as well; they stay bounded and are momentum-zeroed below.
        const int check_band = amrex::min(IBM::ib.ghost_layers + 1, nghost);
        const Dim3 bxlo = amrex::lbound(bx);
        const Dim3 bxhi = amrex::ubound(bx);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          if (tag(i,j,k,tfin) != 1) return;  // valid or fixed
          bool visible = false;
          for (int d = 1; d <= check_band; ++d) {
            visible = visible || (!ibMarkers(i-d, j,   k,   0) && (i-(d-1) >= bxlo.x))
                              || (!ibMarkers(i+d, j,   k,   0) && (i+(d-1) <= bxhi.x))
                              || (!ibMarkers(i,   j-d, k,   0) && (j-(d-1) >= bxlo.y))
                              || (!ibMarkers(i,   j+d, k,   0) && (j+(d-1) <= bxhi.y));
#if (AMREX_SPACEDIM == 3)
            visible = visible || (!ibMarkers(i,   j,   k-d, 0) && (k-(d-1) >= bxlo.z))
                              || (!ibMarkers(i,   j,   k+d, 0) && (k+(d-1) <= bxhi.z));
#endif
            if (visible) break;
          }
          if (visible) Gpu::Atomic::Add(p_unconv, 1);
        });

        // Zero momentum in interior solid cells.
        // Rationale: averaging fluid neighbors from different acoustic
        // phases can amplify velocity (observed solid |u| > fluid |u|
        // by 50-100%). During regrid, these spurious momenta feed into
        // FillPatch/avgDown and contaminate fresh fluid cells. Zeroing
        // is conservative and prevents this amplification loop.
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          if (ibMarkers(i,j,k,0) == 0) return;  // fluid
          if (ibMarkers(i,j,k,1) != 0) return;  // ghost point (handled in Pass 1)

          using PC = PROB::ProbClosures;
          const Real rho = state(i,j,k, PC::URHO);
          if (rho <= Real(0)) return;

          // Subtract kinetic energy from total energy (keep internal energy)
          const Real mx = state(i,j,k, PC::UMX);
          const Real my = state(i,j,k, PC::UMY);
#if (AMREX_SPACEDIM == 3)
          const Real mz = state(i,j,k, PC::UMZ);
          const Real ke = Real(0.5) * (mx*mx + my*my + mz*mz) / rho;
#else
          const Real ke = Real(0.5) * (mx*mx + my*my) / rho;
#endif
          state(i,j,k, PC::UMX) = Real(0.0);
          state(i,j,k, PC::UMY) = Real(0.0);
          state(i,j,k, PC::UMZ) = Real(0.0);  // always: index exists in 2D
          state(i,j,k, PC::UET) -= ke;
        });
      }

      // Single convergence readback for the whole level (replaces up to
      // 32 sync+readback pairs per FAB per step).
      Gpu::streamSynchronize();
      const int n_unconverged = d_unconverged.dataValue();
      if (n_unconverged > 0) {
        amrex::Print() << "\nadvance: Pass-2 flood-fill failed to converge on level "
                       << level << " at t=" << time + dt << ": "
                       << n_unconverged << " stencil-visible interior solid"
                       << " cell(s) left unfilled after " << flood_iters
                       << " iterations.\n";
        amrex::Abort("advance: Pass-2 flood-fill unconverged");
      }
    } // inner block
    } // end Pass 2 (ib_move || pass2_static)

    // ------------------------------------------------------------------------
    // Hard NaN/Inf/non-positive check on fluid cells. If any fluid cell is
    // broken, abort with a useful diagnostic. No silent repair — if this
    // fires, the scheme itself failed and needs fixing (resolution, CFL,
    // positivity limiter, viscosity, ...).
    //
    // When cns.strict_positivity = 1, we also abort if any fluid cell's
    // density or internal energy has collapsed to near the cons2prims
    // clipping floors. That catches
    // silent clipping — the scheme may have kept marching because prims
    // were clamped, but the underlying conservative state is unphysical.
    // The order_rk=0 diagnostic intentionally stores dU/dt in S2; an RHS is
    // not a conservative state and is not required to have positive entries.
    // ------------------------------------------------------------------------
    if (order_rk != 0) {
      const bool strict = CNS::strict_positivity;
      // Density uses 1e6 times its tiny hard floor. Internal energy uses the
      // same margin for a pressure-floor build, but a temperature-floor build
      // stays just above the configured temperature to avoid turning 10 K
      // into an unintended 10^7 K threshold.
      // Capture the floor constants as locals so CUDA device lambdas don't
      // reach back into CNSConstants:: namespace storage.
      const Real smallr_local = CNSConstants::smallr;
      const Real rho_floor = smallr_local * Real(1.0e6);
#if CLIP_TEMPERATURE_MIN
      // A 1e6 multiplier is appropriate for the tiny pressure floor but would
      // turn the configured 10 K temperature floor into 10^7 K.
      const Real floor_scale =
          Real(1.0) + Real(1024.0) * std::numeric_limits<Real>::epsilon();
#else
      const Real floor_scale = Real(1.0e6);
#endif

      ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpMin, ReduceOpMin> rop;
      ReduceData<int, int, Real, Real> rdata(rop);
      using RT = typename decltype(rdata)::Type;

      for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& state = S2.array(mfi);
        const auto& ibMarkers = ib_mf.array(mfi);

        rop.eval(bx, rdata, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT {
          using PC = PROB::ProbClosures;
          if (ibMarkers(i,j,k,0) != 0) return {0, 0, Real(1e300), Real(1e300)};
          const Real rho = state(i,j,k, PC::URHO);
          const Real E   = state(i,j,k, PC::UET);
          const Real mx  = state(i,j,k, PC::UMX);
          const Real my  = state(i,j,k, PC::UMY);
#if (AMREX_SPACEDIM == 3)
          const Real mz  = state(i,j,k, PC::UMZ);
#else
          const Real mz  = Real(0.0);  // UMZ not evolved in 2D
#endif
          // NaN momenta must abort too — they poison ke/eint downstream even
          // when rho and E are still finite (no silent repair, per policy).
          int nonfinite = (!std::isfinite(rho) || !std::isfinite(E) ||
                           !std::isfinite(mx)  || !std::isfinite(my) ||
                           !std::isfinite(mz)) ? 1 : 0;
          Real rhoe = Real(1e300);
          if (!nonfinite && rho > Real(0.0)) {
            rhoe = E - Real(0.5) * (mx*mx + my*my + mz*mz) / rho;
            if (!std::isfinite(rhoe)) nonfinite = 1;
          }

          // Physical admissibility requires positive internal-energy density,
          // not merely positive total energy. cons2prims() applies a primitive
          // floor, so omitting rho*e here silently allowed an inadmissible
          // conservative cell to feed the next RHS evaluation.
          int nonpos = (rho <= Real(0.0) || E <= Real(0.0) ||
                        !std::isfinite(rhoe) || rhoe <= Real(0.0)) ? 1 : 0;
          if (strict && !nonfinite && !nonpos) {
            const Real rhoe_floor =
                cls_d->get_rhoe_min(rho) * floor_scale;
            if (rho < rho_floor || rhoe < rhoe_floor) nonpos = 1;
          }
          return {nonfinite, nonpos, rho,
                  std::isfinite(rhoe) ? rhoe : Real(1e300)};
        });
      }
      auto hv = rdata.value(rop);
      int  icnt[2] = {amrex::get<0>(hv), amrex::get<1>(hv)};  // n_nonfinite, n_nonpos
      Real rmin[2] = {amrex::get<2>(hv), amrex::get<3>(hv)};  // rho_min, rhoe_min
      // Fused: one array allreduce per MPI op (was 4 separate scalar
      // allreduces every step per level — latency-bound on fragmented grids).
      ParallelDescriptor::ReduceIntSum(icnt, 2);
      ParallelDescriptor::ReduceRealMin(rmin, 2);
      int  n_nonfinite = icnt[0];
      int  n_nonpos    = icnt[1];
      Real rho_min     = rmin[0];
      Real rhoe_min    = rmin[1];

      // cns.soft_positivity = 1: clip recoverable bad cells (rho<=0,
      // rhoE<=0, or rho*e<=0, with all components finite) instead of aborting. NaN/Inf
      // remain fatal. The clip replaces each bad cell's (rho, E) with
      // (max(rho, rho_floor), rho*e_floor(rho) + 0.5 |m|^2/rho),
      // which preserves momentum and guarantees physical eint. Not
      // globally conservative; report the clip count every step.
      if (CNS::soft_positivity && n_nonfinite == 0 && n_nonpos > 0) {
        // Device ParallelFor (was a HOST triple loop dereferencing device
        // Array4 pointers — it only worked under amrex.the_arena_is_managed=1
        // and forced a full-FAB page migration whenever it fired). Each cell
        // reads/writes only itself, so results are bitwise identical to the
        // serial loop; the atomic counter is diagnostic-only.
        Gpu::DeviceScalar<int> d_clipped(0);
        int* p_clipped = d_clipped.dataPtr();
        for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          Array4<Real> const& state = S2.array(mfi);
          const auto& ibMarkers = ib_mf.array(mfi);
          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            using PC = PROB::ProbClosures;
            if (ibMarkers(i,j,k,0) != 0) return;
            Real rho_val = state(i,j,k, PC::URHO);
            Real E_val   = state(i,j,k, PC::UET);
            const Real mx = state(i,j,k, PC::UMX);
            const Real my = state(i,j,k, PC::UMY);
#if (AMREX_SPACEDIM == 3)
            const Real mz = state(i,j,k, PC::UMZ);
#else
            const Real mz = Real(0.0);
#endif
            bool need_clip = (rho_val <= Real(0.0) || E_val <= Real(0.0));
            if (!need_clip && rho_val > Real(0.0)) {
              const Real rhoe_val = E_val -
                  Real(0.5) * (mx*mx + my*my + mz*mz) / rho_val;
              need_clip = !std::isfinite(rhoe_val) || rhoe_val <= Real(0.0);
            }
            if (strict && !need_clip) {
              const Real rho_safe = amrex::max(rho_val, smallr_local);
              const Real eint_approx =
                E_val / rho_safe - Real(0.5) * (mx*mx + my*my + mz*mz) / (rho_safe*rho_safe);
              const Real ei_floor =
                  cls_d->get_ei_min(rho_safe) * floor_scale;
              if (rho_val < rho_floor || eint_approx < ei_floor) need_clip = true;
            }
            if (!need_clip) return;
            const Real rho_new = amrex::max(rho_val, rho_floor);
            const Real KE      = Real(0.5) * (mx*mx + my*my + mz*mz) / rho_new;
            const Real ei_floor_new =
                cls_d->get_ei_min(rho_new) * floor_scale;
            const Real ei_new =
                amrex::max((E_val - KE) / rho_new, ei_floor_new);
            state(i,j,k, PC::URHO) = rho_new;
            state(i,j,k, PC::UET)  = rho_new * ei_new + KE;
            Gpu::Atomic::Add(p_clipped, 1);
          });
        }
        Gpu::streamSynchronize();  // conditional path: one readback per triggered step
        int clipped = d_clipped.dataValue();
        ParallelDescriptor::ReduceIntSum(clipped);
        amrex::Print() << "[soft_positivity] level " << level
                       << " t=" << time + dt << " dt=" << dt
                       << " clipped " << clipped << " fluid cells"
                       << " (rho_min=" << rho_min
                       << ", rhoe_min=" << rhoe_min << ")\n";
        // continue — do not abort
      }
      else if (n_nonfinite > 0 || n_nonpos > 0) {
        const char* badlabel = strict
          ? "rho<=0, rhoE<=0, rho*e<=0, rho<rho_floor, or e<e_floor(rho) (strict_positivity=1)"
          : "rho<=0, rhoE<=0, or rho*e<=0";
        amrex::Print() << "\n========================================================\n"
                       << "NUMERICAL FAILURE at level " << level
                       << ", t = " << time + dt << ", dt = " << dt << "\n"
                       << "  fluid cells with NaN/Inf : " << n_nonfinite << "\n"
                       << "  fluid cells with " << badlabel << " : " << n_nonpos << "\n"
                       << "  min(rho) = " << rho_min
                       << "   (rho_floor = " << rho_floor << ")\n"
          << "  min(rho*e) = " << rhoe_min
          << "   (scaled rhoe_floor at min(rho) ~ "
          << cls_h.get_rhoe_min(amrex::max(rho_min, smallr_local)) * floor_scale
          << ")\n"
                       << "========================================================\n";

        // Failure reporting runs on the host.  With managed memory disabled,
        // S2 and the IBM marker FabArray are device allocations and must be
        // copied explicitly before they are inspected here.
        MultiFab state_host(
            S2.boxArray(), S2.DistributionMap(), S2.nComp(), S2.nGrowVect(),
            MFInfo().SetArena(The_Pinned_Arena()));
        IBMultiFab<uint8_t> marker_host(
            ib_mf.boxArray(), ib_mf.DistributionMap(), ib_mf.nComp(),
            ib_mf.nGrow(), MFInfo().SetArena(The_Pinned_Arena()));
        amrex::dtoh_memcpy(state_host, S2);
        amrex::dtoh_memcpy(marker_host, ib_mf);
        Gpu::streamSynchronize();

        // Second pass: locate bad cells and print (i,j,k, x,y,z, cons+neighbors)
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        int printed = 0;
        const int MAX_PRINT = 8;
        for (MFIter mfi(state_host, false); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          const Box& fabbox = mfi.fabbox();
          const auto state = state_host.const_array(mfi);
          const auto ibMarkers = marker_host.const_array(mfi);
          const int lo0 = bx.smallEnd(0), lo1 = bx.smallEnd(1);
          const int hi0 = bx.bigEnd(0),   hi1 = bx.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
          const int lo2 = bx.smallEnd(2), hi2 = bx.bigEnd(2);
#else
          const int lo2 = 0, hi2 = 0;
#endif
          for (int k = lo2; k <= hi2 && printed < MAX_PRINT; ++k) {
          for (int j = lo1; j <= hi1 && printed < MAX_PRINT; ++j) {
          for (int i = lo0; i <= hi0 && printed < MAX_PRINT; ++i) {
            using PC = PROB::ProbClosures;
            if (ibMarkers(i,j,k,0) != 0) continue;
            const Real rho = state(i,j,k, PC::URHO);
            const Real E   = state(i,j,k, PC::UET);
            const Real mxp = state(i,j,k, PC::UMX);
            const Real myp = state(i,j,k, PC::UMY);
#if (AMREX_SPACEDIM == 3)
            const Real mzp = state(i,j,k, PC::UMZ);
#else
            const Real mzp = Real(0.0);
#endif
            const bool finite = std::isfinite(rho) && std::isfinite(E) &&
                                std::isfinite(mxp) && std::isfinite(myp) &&
                                std::isfinite(mzp);
            bool bad = !finite || rho <= Real(0.0) || E <= Real(0.0);
            Real eint_approx = std::numeric_limits<Real>::quiet_NaN();
            if (finite && rho > Real(0.0)) {
              eint_approx = E / rho -
                  Real(0.5) * (mxp*mxp + myp*myp + mzp*mzp) / (rho*rho);
              if (!(eint_approx > Real(0.0)) ||
                  (strict &&
                   (rho < rho_floor ||
                    eint_approx < cls_h.get_ei_min(rho) * floor_scale))) {
                bad = true;
              }
            }
            if (!bad) continue;
            const Real x = plo[0] + (i + Real(0.5)) * dx[0];
            const Real y = plo[1] + (j + Real(0.5)) * dx[1];
            amrex::AllPrint() << "[BAD CELL #" << printed << "] level=" << level
                              << " (i,j,k)=(" << i << "," << j << "," << k << ")"
                              << " (x,y)=(" << x << "," << y << ")\n"
                              << "  center:  rho=" << rho
                              << "  mx=" << state(i,j,k,PC::UMX)
                              << "  my=" << state(i,j,k,PC::UMY)
                              << "  E=" << E
                              << "  ei=" << eint_approx
                              << "  ibm0=" << int(ibMarkers(i,j,k,0))
                              << "  ibm1=" << int(ibMarkers(i,j,k,1)) << "\n";
            // 3x3 neighborhood dump
            amrex::AllPrint() << "  rho 3x3 (j-1..j+1, i-1..i+1):\n";
            for (int dj = 1; dj >= -1; --dj) {
              amrex::AllPrint() << "    ";
              for (int di = -1; di <= 1; ++di) {
                const int ii = i+di, jj = j+dj;
                const IntVect niv(AMREX_D_DECL(ii, jj, k));
                if (fabbox.contains(niv)) {
                  amrex::AllPrint()
                      << std::setw(13) << std::setprecision(4)
                      << state(ii,jj,k,PC::URHO)
                      << "[" << int(ibMarkers(ii,jj,k,0)) << "] ";
                } else {
                  amrex::AllPrint() << std::setw(16) << "n/a" << " ";
                }
              }
              amrex::AllPrint() << "\n";
            }
            printed++;
          }}}
        }
        amrex::Print() << "\nFix the root cause — do not silently patch.\n"
                       << "========================================================\n";
        amrex::Abort("advance: non-finite / non-positive fluid state detected");
      }
    }
  }
#endif  // AMREX_USE_GPIBM

#ifndef AMREX_USE_GPIBM
  // Universal post-step state validity check for non-IBM builds (the detailed
  // IBM reporter above is compiled out here; long production runs previously
  // had NO NaN/positivity safety net). Cheap: two reductions per level-step.
  // Policy: NaN-abort with location diagnostics — never silently continue.
  {
    static const int s_check_state = [] {
      int v = 1;
      amrex::ParmParse pp("cns");
      pp.query("check_state", v);
      return v;
    }();
    // order_rk=0 stores dU/dt in the state MultiFab for verification.  The
    // resulting RHS is not a conservative state and need not be positive.
    if (s_check_state && order_rk != 0) {
      MultiFab& Schk = get_new_data(State_Type);
      const Real rho_min = Schk.min(PROB::ProbClosures::URHO, 0);
      const bool has_nan = Schk.contains_nan(0, PROB::ProbClosures::NCONS, 0) ||
                           Schk.contains_inf(0, PROB::ProbClosures::NCONS, 0);
      // Min internal-energy density rho*e = E - |m|^2/(2 rho): catches
      // negative-p/T/c^2 states that finite + rho>0 checks miss.
      ReduceOps<ReduceOpMin> ei_op;
      ReduceData<Real> ei_data(ei_op);
      using EiTuple = typename decltype(ei_data)::Type;
      for (MFIter mfi(Schk, false); mfi.isValid(); ++mfi) {
        const Box& bxc = mfi.tilebox();
        const Array4<Real>& a = Schk.array(mfi);
        ei_op.eval(bxc, ei_data,
                   [=] AMREX_GPU_DEVICE(int i, int j, int k) -> EiTuple {
          using icls = PROB::ProbClosures;
          const Real rhoc = a(i, j, k, icls::URHO);
          const Real ke = Real(0.5) *
              (a(i, j, k, icls::UMX) * a(i, j, k, icls::UMX) +
               a(i, j, k, icls::UMY) * a(i, j, k, icls::UMY) +
               a(i, j, k, icls::UMZ) * a(i, j, k, icls::UMZ)) /
              amrex::max(rhoc, Real(1e-300));
          return {a(i, j, k, icls::UET) - ke};
        });
      }
      const Real ei_min = amrex::get<0>(ei_data.value(ei_op));
      bool stats_bad = false;
      if (compute_stats) {
        stats_bad = get_new_data(Stats_Type)
                        .contains_nan(0, PROB::ProbClosures::NSTAT, 0);
      }
      if (has_nan || stats_bad || !(rho_min > Real(0.0)) ||
          !(ei_min > Real(0.0))) {
        amrex::Print() << "\n==================== STATE CHECK FAILED ====================\n"
                       << "  level = " << level
                       << "  step = " << parent->levelSteps(level)
                       << "  time = " << state[State_Type].curTime() << "\n"
                       << "  min(rho) = " << rho_min
                       << "  min(rho*e_int) = " << ei_min
                       << "  contains_nan = " << (has_nan ? "YES" : "no")
                       << "  stats_nan = " << (stats_bad ? "YES" : "no") << "\n"
                       << "  Fix the root cause — do not silently patch.\n"
                       << "============================================================\n";
        amrex::Abort("advance: non-finite / non-positive state (non-IBM check)");
      }
    }
  }
#endif

  return dt;
}
