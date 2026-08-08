#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>
#include <CNS.h>
#include <CNS_K.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <ibm_positivity_limiter.h>

#include <type_traits>
#include <utility>
#endif

#ifdef CNS_USE_FSI
#include <fsi/Kinematics.h>
#include <fsi/RigidBodyProperties.h>
#endif

#ifdef USE_PELEPHYSICS
#include "TransPele.h"

pele::physics::PeleParams<
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type
  >> trans_parms;
#endif


using namespace amrex;

#ifdef AMREX_USE_GPIBM
namespace {

template <typename ProbParm, typename = void>
struct HasIBMRefineSupportHook : std::false_type {};

template <typename ProbParm>
struct HasIBMRefineSupportHook<
    ProbParm,
    std::void_t<decltype(ibm_support_refine_allowed(
        std::declval<Real>(), std::declval<Real>(), std::declval<Real>(),
        std::declval<int>(), std::declval<ProbParm const&>()))>>
    : std::true_type {};

template <typename ProbParm>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
problem_ibm_support_refine_allowed(
    Real x, Real y, Real z, int parent_level,
    ProbParm const& prob_parm) noexcept
{
  if constexpr (HasIBMRefineSupportHook<ProbParm>::value) {
    return ibm_support_refine_allowed(
        x, y, z, parent_level, prob_parm);
  }
  return true;
}

} // namespace
#endif

bool CNS::verbose = true;
bool CNS::record_probe = false;
bool CNS::dt_dynamic = false;
bool CNS::ib_move = false;
bool CNS::plot_surf = false;

int CNS::surf_int = 10000000;
std::string CNS::surf_filename = "surfplot";


// utilities
bool CNS::use_utility = false; 
Utility CNS::utilidades;

Real CNS::eb_weight = 0.0;
bool CNS::eb_redistribution = false;
std::string CNS::eb_redistribution_type = "NoRedist";

int CNS::nstep_screen_output = 10;
int CNS::order_rk = 2;
int CNS::stages_rk = 2;
bool CNS::use_nscbc = false;
GpuArray<int, AMREX_SPACEDIM> CNS::nscbc_lo = {AMREX_D_DECL(0, 0, 0)};
GpuArray<int, AMREX_SPACEDIM> CNS::nscbc_hi = {AMREX_D_DECL(0, 0, 0)};
nscbc::Parm CNS::nscbc_parm{};
bool CNS::strict_positivity = false;
bool CNS::ibm_positivity_flux_limiter = false;
bool CNS::ibm_positivity_flux_limiter_verbose = false;
bool CNS::ibm_positivity_retry = false;
int CNS::ibm_positivity_max_step_halvings = 4;
bool CNS::soft_positivity = false;
bool CNS::pass2_static = false;
int CNS::do_reflux = 0; // default reflux is off
int CNS::amr_state_interp = 0;
bool CNS::regrid_prolongation_diagnostics = false;
int CNS::refine_max_dengrad_lev = -1;
Real CNS::cfl = 0.0_rt;
Real CNS::dt_constant = 0.0_rt;
Real CNS::dt_max = std::numeric_limits<Real>::max(); // disabled by default
bool CNS::restart_first_dt_from_cfl = false;
Real CNS::refine_dengrad = 1.0e10;
int  CNS::INDEX_THERM = 0;
bool CNS::compute_stats = false;
bool CNS::record_stats = false;
Real CNS::time_stats = 0.0;
Real CNS::time_stat_level[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Real CNS::stats_start_time = 0.0;

PROB::ProbClosures *CNS::h_prob_closures = nullptr;
PROB::ProbClosures *CNS::d_prob_closures = nullptr;
PROB::ProbParm *CNS::h_prob_parm = nullptr;
PROB::ProbParm *CNS::d_prob_parm = nullptr;
BCRec *CNS::h_phys_bc = nullptr;
BCRec *CNS::d_phys_bc = nullptr;

// needed for CNSBld - derived from LevelBld (abstract class, pure virtual
// functions must be implemented)

CNS::CNS() {}

CNS::CNS(Amr &papa, int lev, const Geometry &level_geom, const BoxArray &bl,
         const DistributionMapping &dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time) {
  if (do_reflux && level > 0) {
    flux_reg.reset(new FluxRegister(grids, dmap, crse_ratio, level,PROB::ProbClosures::NCONS));
  }

#ifdef AMREX_USE_GPIBM
  IBM::ib.build_mf(grids, dmap, level);
  IBM::ib.computeMarkers(level);
#endif

#ifdef CNS_USE_EB
  EBM::eb.build_mf(grids, dmap, level);
#endif

  buildMetrics();

  if (use_nscbc) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if ((nscbc_lo[dir] > 0 || nscbc_hi[dir] > 0) && Geom().isPeriodic(dir)) {
        amrex::Abort("NSCBC cannot be enabled on a periodic direction");
      }
    }
  }

  rz_sanity_check(Geom());
};

CNS::~CNS() {}
// -----------------------------------------------------------------------------

// init ------------------------------------------------------------------------

void CNS::read_params() {

  ParmParse pp("cns");

#ifdef AMREX_USE_GPIBM
  amrex::Print()
      << "[IBM-AMR-Support-Hook] detected="
      << int(HasIBMRefineSupportHook<PROB::ProbParm>::value)
      << " dispatch=automatic default=whole_interface\n";
#endif

  pp.query("nstep_screen_output", nstep_screen_output);
  pp.query("verbose", verbose);

  int restart_first_dt_from_cfl_value = 0;
  pp.query("restart_first_dt_from_cfl", restart_first_dt_from_cfl_value);
  if (restart_first_dt_from_cfl_value != 0 &&
      restart_first_dt_from_cfl_value != 1) {
    amrex::Abort("cns.restart_first_dt_from_cfl must be 0 or 1");
  }
  restart_first_dt_from_cfl = restart_first_dt_from_cfl_value == 1;
  if (restart_first_dt_from_cfl && !dt_dynamic) {
    amrex::Abort(
        "cns.restart_first_dt_from_cfl=1 requires CFL-based time stepping");
  }
  if (restart_first_dt_from_cfl) {
    amrex::Print()
        << "  cns.restart_first_dt_from_cfl = 1 (checkpoint restart only: "
           "the first normal computeNewDt uses the current CFL estimate)\n";
  }

  pp.query("record_probe", record_probe);
  pp.query("record_stats", record_stats);
  pp.query("stats_start_time", stats_start_time);

  Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    h_phys_bc->setLo(i, lo_bc[i]);
    h_phys_bc->setHi(i, hi_bc[i]);
  }

  Vector<int> nslo(AMREX_SPACEDIM, 0);
  Vector<int> nshi(AMREX_SPACEDIM, 0);
  pp.queryarr("nscbc_lo", nslo, 0, AMREX_SPACEDIM);
  pp.queryarr("nscbc_hi", nshi, 0, AMREX_SPACEDIM);
  use_nscbc = false;
  bool has_relaxed_inflow = false;
  bool has_pressure_outflow = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (nslo[dir] < 0 || nslo[dir] > 3 || nshi[dir] < 0 || nshi[dir] > 3) {
      amrex::Abort("cns.nscbc_lo/hi entries must be 0, 1, 2, or 3");
    }
    nscbc_lo[dir] = nslo[dir];
    nscbc_hi[dir] = nshi[dir];
    use_nscbc = use_nscbc || nslo[dir] > 0 || nshi[dir] > 0;
    has_relaxed_inflow = has_relaxed_inflow || nslo[dir] == 1 || nshi[dir] == 1;
    has_pressure_outflow =
      has_pressure_outflow || nslo[dir] == 3 || nshi[dir] == 3;
  }

  if (use_nscbc) {
    pp.query("nscbc_Lchar", nscbc_parm.Lchar);
    pp.query("nscbc_Mmax", nscbc_parm.Mmax);
    pp.query("nscbc_Ptarget", nscbc_parm.Ptarget);
    pp.query("nscbc_sigma", nscbc_parm.sigma);
    pp.query("nscbc_utarget", nscbc_parm.utarget);
    pp.query("nscbc_vtarget", nscbc_parm.vtarget);
    pp.query("nscbc_wtarget", nscbc_parm.wtarget);
    pp.query("nscbc_Ttarget", nscbc_parm.Ttarget);
    pp.query("nscbc_eta", nscbc_parm.eta);
    pp.query("nscbc_transverse_relax", nscbc_parm.transverse_relax);
    pp.query("nscbc_outflow_transverse_model",
             nscbc_parm.outflow_transverse_model);
    pp.query("nscbc_ghost_update_model", nscbc_parm.ghost_update_model);
    pp.query("nscbc_order", nscbc_parm.derivative_order);
    pp.query("nscbc_use_transverse", nscbc_parm.use_transverse);
    pp.query("nscbc_min_rho", nscbc_parm.min_rho);
    pp.query("nscbc_min_T", nscbc_parm.min_T);
    pp.query("nscbc_min_p", nscbc_parm.min_p);

    if (has_relaxed_inflow) {
      const bool targets_present =
        pp.contains("nscbc_utarget") && pp.contains("nscbc_Ttarget") &&
        pp.contains("nscbc_eta")
#if (AMREX_SPACEDIM >= 2)
        && pp.contains("nscbc_vtarget")
#endif
#if (AMREX_SPACEDIM == 3)
        && pp.contains("nscbc_wtarget")
#endif
        ;
      if (!targets_present) {
        amrex::Abort(
          "NSCBC type 1 requires explicit velocity, temperature, and eta targets "
          "in solver units");
      }
    }
    if (has_pressure_outflow &&
        !(pp.contains("nscbc_Ptarget") && pp.contains("nscbc_sigma") &&
          pp.contains("nscbc_Lchar") && pp.contains("nscbc_Mmax"))) {
      amrex::Abort(
        "NSCBC type 3 requires explicit Ptarget, sigma, Lchar, and Mmax");
    }

    if (nscbc_parm.Lchar <= Real(0.0)) {
      amrex::Abort("cns.nscbc_Lchar must be positive");
    }
    if (nscbc_parm.Mmax < Real(0.0) || nscbc_parm.Mmax >= Real(1.0)) {
      amrex::Abort("cns.nscbc_Mmax must satisfy 0 <= Mmax < 1");
    }
    if (nscbc_parm.sigma < Real(0.0) || nscbc_parm.eta < Real(0.0)) {
      amrex::Abort("cns.nscbc_sigma and cns.nscbc_eta must be non-negative");
    }
    if (nscbc_parm.use_transverse != 0 && nscbc_parm.use_transverse != 1) {
      amrex::Abort("cns.nscbc_use_transverse must be 0 or 1");
    }
    if (nscbc_parm.transverse_relax < Real(0.0)) {
      amrex::Abort("cns.nscbc_transverse_relax must be non-negative");
    }
    if (nscbc_parm.outflow_transverse_model <
          nscbc::OUTFLOW_TRANSVERSE_PROJECTED_NSCBC ||
        nscbc_parm.outflow_transverse_model >
          nscbc::OUTFLOW_TRANSVERSE_GILES2) {
      amrex::Abort(
        "cns.nscbc_outflow_transverse_model must be 0 (legacy projected "
        "NSCBC) or 1 (Giles second-order 2-D outflow)");
    }
    if (nscbc_parm.ghost_update_model <
          nscbc::GHOST_UPDATE_PERSISTENT_ODE ||
        nscbc_parm.ghost_update_model >
          nscbc::GHOST_UPDATE_STAGE_LOCAL_GC) {
      amrex::Abort(
        "cns.nscbc_ghost_update_model must be 0 (persistent ODE) or 1 "
        "(stage-local GC-NSCBC)");
    }
    if (nscbc_parm.ghost_update_model ==
          nscbc::GHOST_UPDATE_STAGE_LOCAL_GC &&
        nscbc_parm.outflow_transverse_model !=
          nscbc::OUTFLOW_TRANSVERSE_GILES2) {
      amrex::Abort(
        "stage-local GC-NSCBC currently requires the Giles second-order "
        "outflow model");
    }
    if (nscbc_parm.ghost_update_model ==
          nscbc::GHOST_UPDATE_STAGE_LOCAL_GC &&
        has_relaxed_inflow) {
      amrex::Abort(
        "stage-local GC-NSCBC currently supports subsonic outflow only; "
        "retain the accepted characteristic inlet implementation");
    }
    if (nscbc_parm.outflow_transverse_model ==
          nscbc::OUTFLOW_TRANSVERSE_GILES2 &&
        nscbc_parm.use_transverse == 0) {
      amrex::Abort(
        "Giles second-order outflow requires cns.nscbc_use_transverse=1");
    }
#if (AMREX_SPACEDIM != 2)
    if (nscbc_parm.outflow_transverse_model ==
          nscbc::OUTFLOW_TRANSVERSE_GILES2 ||
        nscbc_parm.ghost_update_model ==
          nscbc::GHOST_UPDATE_STAGE_LOCAL_GC) {
      amrex::Abort(
        "Giles/GC second-order outflow is currently implemented only in 2-D");
    }
#endif
    if (nscbc_parm.min_rho <= Real(0.0) || nscbc_parm.min_T <= Real(0.0) ||
        nscbc_parm.min_p <= Real(0.0)) {
      amrex::Abort("cns.nscbc_min_rho/min_T/min_p must be positive");
    }
    if (nscbc_parm.derivative_order != 1 &&
        nscbc_parm.derivative_order != 2) {
      amrex::Abort("cns.nscbc_order must be 1 or 2");
    }

    amrex::Print() << "  NSCBC boundary model: lo="
                   << AMREX_D_TERM(nscbc_lo[0], << " " << nscbc_lo[1],
                                   << " " << nscbc_lo[2])
                   << " hi="
                   << AMREX_D_TERM(nscbc_hi[0], << " " << nscbc_hi[1],
                                   << " " << nscbc_hi[2])
                   << " order=" << nscbc_parm.derivative_order
                   << " transverse=" << nscbc_parm.use_transverse
                   << " transverse_model="
                   << nscbc_parm.outflow_transverse_model
                   << " ghost_update_model="
                   << nscbc_parm.ghost_update_model
                   << " beta=" << nscbc_parm.transverse_relax << "\n";
  }

  pp.query("do_reflux", do_reflux);
  int regrid_prolongation_diagnostics_value = 0;
  pp.query("regrid_prolongation_diagnostics",
           regrid_prolongation_diagnostics_value);
  if (regrid_prolongation_diagnostics_value != 0 &&
      regrid_prolongation_diagnostics_value != 1) {
    amrex::Abort(
        "cns.regrid_prolongation_diagnostics must be 0 or 1");
  }
  regrid_prolongation_diagnostics =
      regrid_prolongation_diagnostics_value == 1;
#ifndef AMREX_USE_GPIBM
  if (regrid_prolongation_diagnostics) {
    amrex::Abort(
        "cns.regrid_prolongation_diagnostics=1 requires a GP-IBM build");
  }
#else
  if (regrid_prolongation_diagnostics) {
    amrex::Print()
        << "[AMR-Regrid-Prolongation-Diagnostic] enabled=1 read_only=1 "
           "retains_old_overlap_mask=1 retains_coarse_source_stencil=1\n";
  }
#endif
  std::string amr_state_interp_name = "linear";
  pp.query("amr_state_interp", amr_state_interp_name);
  if (amr_state_interp_name == "linear" ||
      amr_state_interp_name == "euler_admissible_linear") {
    amr_state_interp = 0;
  } else if (amr_state_interp_name == "conservative_quartic") {
    amr_state_interp = 1;
  } else if (amr_state_interp_name == "legacy_linear") {
    amr_state_interp = 2;
  } else {
    amrex::Abort(
        "cns.amr_state_interp must be linear, euler_admissible_linear, "
        "legacy_linear, or conservative_quartic");
  }
  if (amr_state_interp == 1) {
    ParmParse amr_pp("amr");
    Vector<int> refinement_ratios;
    if (amr_pp.queryarr("ref_ratio", refinement_ratios)) {
      for (const int ratio : refinement_ratios) {
        if (ratio != 2) {
          amrex::Abort(
              "cns.amr_state_interp=conservative_quartic requires "
              "amr.ref_ratio=2 on every refined level");
        }
      }
    }
  }

#if AMREX_USE_GPIBM || CNS_USE_EB
  if (amr_state_interp != 0) {
    amrex::Abort(
        "pure GP-IBM and EB builds require the closure-admissibility-preserving "
        "limited-linear coarse-fine interpolation selected by "
        "cns.amr_state_interp=linear");
  }
#endif

  if (!pp.query("order_rk", order_rk)) {
    amrex::Abort(
        "Need to specify SSPRK scheme order of accuracy, order_rk={-2, 1, 2, "
        "3}");
  }

  if (!pp.query("stages_rk", stages_rk)) {
    amrex::Abort("Need to specify SSPRK number of stages, stages_rk");
  } else {
    if (order_rk == 1 && stages_rk != 1) {
      amrex::Abort("Forward Euler number of stages must be 1");
    }

    if (order_rk == 2 && stages_rk < order_rk) {
      amrex::Abort(
          "SSPRK2 number of stages must equal or greater than order of "
          "accuracy");
    }
    if (order_rk == 3 && !(stages_rk == 4 || stages_rk == 3)) {
      amrex::Abort("SSPRK3 number of stages must equal 3 or 4");
    }
  }

  pp.query("strict_positivity", strict_positivity);
  if (strict_positivity) {
    amrex::Print() << "  cns.strict_positivity = 1 (abort if state approaches smallr/ei_min floors)\n";
  }

  int legacy_stage_positivity = 0;
  int legacy_stage_positivity_verbose = 0;
  pp.query("stage_positivity", legacy_stage_positivity);
  pp.query("stage_positivity_verbose", legacy_stage_positivity_verbose);
  if ((legacy_stage_positivity != 0 && legacy_stage_positivity != 1) ||
      (legacy_stage_positivity_verbose != 0 &&
       legacy_stage_positivity_verbose != 1)) {
    amrex::Abort(
        "cns.stage_positivity and cns.stage_positivity_verbose must be 0 or 1");
  }

  pp.query("ibm_positivity_flux_limiter", ibm_positivity_flux_limiter);
  pp.query("ibm_positivity_flux_limiter_verbose",
           ibm_positivity_flux_limiter_verbose);
  pp.query("ibm_positivity_retry", ibm_positivity_retry);
  pp.query("ibm_positivity_max_step_halvings",
           ibm_positivity_max_step_halvings);
  if (ibm_positivity_flux_limiter) {
#ifndef AMREX_USE_GPIBM
    amrex::Abort(
        "cns.ibm_positivity_flux_limiter=1 requires a GP-IBM build");
#endif
    if (order_rk != 3 || stages_rk != 4) {
      amrex::Abort(
          "cns.ibm_positivity_flux_limiter=1 requires SSPRK(4,3)");
    }
    amrex::Print()
        << "  cns.ibm_positivity_flux_limiter = 1 "
           "(per-face convex limiting of every SSPRK43 FE bracket; "
           "fixed pure-GP Cartesian-2D/3D or RZ-2D; multilevel reflux is fail-closed "
           "until an invariant-domain-preserving synchronization is "
           "qualified)\n";
  }
  if (legacy_stage_positivity != 0 &&
      !ibm_positivity_flux_limiter) {
    amrex::Abort(
        "cns.stage_positivity=1 was historically ignored and is not a "
        "standalone RK-stage check; enable "
        "cns.ibm_positivity_flux_limiter=1 for fail-closed SSPRK43 "
        "Forward-Euler-bracket protection");
  }
  if (legacy_stage_positivity_verbose != 0 &&
      legacy_stage_positivity == 0) {
    amrex::Abort(
        "cns.stage_positivity_verbose=1 requires cns.stage_positivity=1");
  }
  if (legacy_stage_positivity != 0) {
    amrex::Print()
        << "  cns.stage_positivity = 1 (legacy request validated; actual "
           "stage protection is cns.ibm_positivity_flux_limiter=1)\n";
  }
  if (ibm_positivity_flux_limiter_verbose &&
      !ibm_positivity_flux_limiter) {
    amrex::Abort(
        "cns.ibm_positivity_flux_limiter_verbose=1 requires "
        "cns.ibm_positivity_flux_limiter=1");
  }
  if (ibm_positivity_retry && !ibm_positivity_flux_limiter) {
    amrex::Abort(
        "cns.ibm_positivity_retry=1 requires "
        "cns.ibm_positivity_flux_limiter=1");
  }
  if (ibm_positivity_max_step_halvings < 0 ||
      ibm_positivity_max_step_halvings > 20) {
    amrex::Abort(
        "cns.ibm_positivity_max_step_halvings must be in [0,20]");
  }
  if (ibm_positivity_retry) {
    amrex::Abort(
        "cns.ibm_positivity_retry is reserved but not yet certified for the "
        "pure full-cell GP-IBM SSPRK43 transaction path; disable it");
  }


  pp.query("soft_positivity", soft_positivity);
  if (soft_positivity) {
    amrex::Print() << "  cns.soft_positivity = 1 (clip bad fluid cells to floors instead of aborting; NaN/Inf still aborts)\n";
  }
  if (ibm_positivity_flux_limiter && soft_positivity) {
    amrex::Abort(
        "cns.ibm_positivity_flux_limiter is conservative and incompatible "
        "with cns.soft_positivity=1 clipping");
  }

  pp.query("pass2_static", pass2_static);
  if (pass2_static) {
    amrex::Print() << "  cns.pass2_static = 1 (Pass 2 flood-fill runs for static geometry too)\n";
  }

  //  Utilities options ----------------------------------------------------
  // specific keywords for Utilities
  pp.query("use_utility",use_utility);
  if (use_utility)
  {
    amrex::Print() << " Using Utilities " << std::endl;      
    ParmParse pp_util("util");

    // PMF
    bool use_PMF=false;
    pp_util.query("use_PMF",use_PMF); 
    if (use_PMF){
#ifdef USE_PELEPHYSICS      
      amrex::Print() << " Reading PMF from file.. " << std::endl;      
      CNS::utilidades.initPMF();
#else      
      amrex::Abort("using PMF files need PelePhysics");
#endif
    }  

    // Read from file
    bool use_turb_file  = false;
    pp_util.query("use_turb_file",use_turb_file); 
    if (use_turb_file){
      std::string turbfilename;
      pp_util.query("turb_file",turbfilename); 
      amrex::Print() << " Reading turbulence from file: " << turbfilename << std::endl; 
      CNS::utilidades.initTurbulenceFile(turbfilename);
    }

  }

#ifdef AMREX_USE_GPIBM
  // IBM-specific input keywords (ib.* namespace)
  ParmParse ppib("ib");
  if (!ppib.query("move", ib_move)) {
    amrex::Abort("ib.move not specified (0=false, 1=true)");
  }
#ifndef CNS_USE_FSI
  if (ib_move) {
    amrex::Abort("ib.move=1 requires a USE_FSI=TRUE build");
  }
#endif
  if (!ppib.query("plot_surf", plot_surf)) {
    amrex::Abort("ib.plot_surf not specified (0=false, 1=true)");
  }
  if (plot_surf) {
    ppib.query("surf_int", surf_int);
    ppib.query("surf_file", surf_filename);
    if (surf_int <= 0) {
      amrex::Abort("ib.surf_int must be positive when ib.plot_surf=1");
    }
  }
#endif
  
#if CNS_USE_EB 
  // specific keywords for EB boundaries
  ParmParse ppeb2("eb2");  
  ppeb2.query("eb_weight",eb_weight); 
  ppeb2.query("redistribution_type", eb_redistribution_type);
  if (eb_redistribution_type != "StateRedist" && eb_redistribution_type != "FluxRedist" &&
      eb_redistribution_type != "NoRedist"    && eb_redistribution_type != "NewRedist") {
    amrex::Abort( " input file: redistribution_type must be StateRedist/FluxRedist/NewRedist/NoRedist");
  }
  if (eb_redistribution_type != "NoRedist") eb_redistribution =true;
  // This communicates to the class (not very elegant)
  EBM::eb.eb_weight = eb_weight;
  EBM::eb.redistribution_type = eb_redistribution_type; 

#endif


#ifdef USE_PELEPHYSICS
  // One-time transport parameter initialization (host->device)
  static bool trans_inited = false;
  if (!trans_inited) {
    trans_parms.initialize();
    trans_inited = true;
  }
#endif


#if AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(d_prob_closures, h_prob_closures,
                          sizeof(PROB::ProbClosures));
  amrex::Gpu::htod_memcpy(d_prob_parm, h_prob_parm, sizeof(PROB::ProbParm));
  amrex::Gpu::htod_memcpy(d_phys_bc, h_phys_bc, sizeof(BCRec));
#endif
}

#ifdef AMREX_USE_GPIBM
void CNS::captureRegridProlongationContext(const CNS* old_level,
                                           const Real time) {
  regrid_new_from_coarse_mask.reset();
  regrid_coarse_source_state.reset();
  if (!regrid_prolongation_diagnostics || level == 0) return;

  regrid_new_from_coarse_mask =
      std::make_unique<iMultiFab>(grids, dmap, 1, 0);
  regrid_new_from_coarse_mask->setVal(1);

  if (old_level != nullptr) {
    iMultiFab old_fine_coverage(old_level->boxArray(),
                                old_level->DistributionMap(), 1, 0);
    old_fine_coverage.setVal(0);
    regrid_new_from_coarse_mask->ParallelCopy(
        old_fine_coverage, 0, 0, 1, 0, 0, geom.periodicity());
  }

  const IntVect ref_ratio = parent->refRatio(level - 1);
  BoxArray coarse_source_boxes = grids;
  coarse_source_boxes.coarsen(ref_ratio);
  regrid_coarse_source_state = std::make_unique<MultiFab>(
      coarse_source_boxes, dmap, PROB::ProbClosures::NCONS, 1);
  regrid_coarse_source_state->setVal(
      std::numeric_limits<Real>::quiet_NaN());
  CNS& coarse_level = getLevel(level - 1);
  FillPatch(coarse_level, *regrid_coarse_source_state, 1, time, State_Type, 0,
            PROB::ProbClosures::NCONS);

  const Long new_from_coarse =
      regrid_new_from_coarse_mask->sum(0, 0, false);
  amrex::Print()
      << "[AMR-Regrid-Prolongation-Context] level=" << level
      << " t=" << std::setprecision(17) << time
      << " fine_valid=" << grids.numPts()
      << " old_overlap=" << grids.numPts() - new_from_coarse
      << " new_from_coarse=" << new_from_coarse
      << " coarse_source_ngrow=1 ref_ratio=" << ref_ratio << '\n';
}
#endif

void CNS::init(AmrLevel &old) {
  auto &oldlev = dynamic_cast<CNS &>(old);

  // Preserve a pending restart-only dt exception if level 0 is rebuilt for
  // load balancing/regridding before it is consumed.
  restart_first_dt_from_cfl_control =
      oldlev.restart_first_dt_from_cfl_control;

  // amrex::Print( ) << " oo CNS::init (AMR recast) -----  " << std::endl;  

  Real dt_new = parent->dtLevel(level);
  Real cur_time = oldlev.state[State_Type].curTime();
  Real prev_time = oldlev.state[State_Type].prevTime();
  Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

#ifdef AMREX_USE_GPIBM
  captureRegridProlongationContext(&oldlev, cur_time);
#endif

  MultiFab &S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0,PROB::ProbClosures::NCONS);

  if (compute_stats){
    MultiFab &Sstat_new = get_new_data(Stats_Type);
    FillPatch(old, Sstat_new, 0, cur_time, Stats_Type, 0,PROB::ProbClosures::NSTAT);
  }

}

void CNS::init() {

  // amrex::Print( ) << " oo CNS::init -----  " << std::endl;  

  Real dt = parent->dtLevel(level);
  Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  Real prev_time = getLevel(level - 1).state[State_Type].prevTime();
  Real dt_old = (cur_time - prev_time) /
                static_cast<Real>(parent->MaxRefRatio(level - 1));
  setTimeLevel(cur_time, dt_old, dt);

#ifdef AMREX_USE_GPIBM
  captureRegridProlongationContext(nullptr, cur_time);
#endif

  MultiFab &S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0,PROB::ProbClosures::NCONS);

  if (compute_stats) {
    // A brand-new level has no accumulation history: start its statistics
    // fresh. Interpolating running means from the coarse level would import
    // the coarse level's accumulation clock; StateData::define leaves the
    // MultiFab uninitialized otherwise (garbage propagates via avgDown).
    get_new_data(Stats_Type).setVal(0.0);
    time_stat_level[level] = 0.0;
  }
};

//-----
void CNS::initData() {
  BL_PROFILE("CNS::initData()");

  const auto geomdata = geom.data();
  MultiFab &S_new = get_new_data(State_Type);
  auto const &sma = S_new.arrays();

  PROB::ProbClosures const *lclosures = d_prob_closures;
  PROB::ProbParm const *lprobparm = d_prob_parm;

  //amrex::Print( ) << "  calling  prob_init in prob.h ...  " << std::endl; 

  // Initialise problem by calling user-given prob.h
#if USE_UTILITY
  amrex::ParallelFor(
    S_new, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      prob_initdata(i, j, k, sma[box_no], geomdata, *lclosures, *lprobparm, use_utility ? &CNS::utilidades : nullptr);
  });
#else
  amrex::ParallelFor(
      S_new, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        prob_initdata(i, j, k, sma[box_no], geomdata, *lclosures, *lprobparm);        
  });
#endif

            
  // Initialise stats 
  if (compute_stats) {
    setupStats();
  }


   prob_rhs.init_coeffs(); // CGPT dixit


}

void CNS::buildMetrics() {

  if (verbose) {
    const Real *dx = geom.CellSize();
    amrex::Print() << "Mesh size (dx,dy,dz) = ";
    amrex::Print() << AMREX_D_TERM(dx[0], << "  " << dx[1], << "  " << dx[2]) << "  \n";
  }  

}

void CNS::post_init(Real stop_time) {

  //amrex::Print() << " oo CNS::post_init level= "  << level << std::endl;

  if (use_nscbc) {
    initialize_nscbc_ghost_state(state[State_Type].curTime());
  }

  if (level > 0) {
    return;
  };

  for (int k = parent->finestLevel() - 1; k >= 0; --k) {
    getLevel(k).avgDown();
  }

  if (verbose) {
    printTotal();
  }
  
  // Set up diagnostics
  if (record_probe) {
    setupTimeProbe();
  }


  
#if CNS_USE_EB
  EBM::eb.check_geometry(level);
#endif
  

}
// -----------------------------------------------------------------------------

// Time-stepping ---------------------------------------------------------------
void CNS::computeInitialDt(int finest_level, int sub_cycle,
                           Vector<int> &n_cycle,  // no. of subcycling steps
                           const Vector<IntVect> &ref_ratio,
                           Vector<Real> &dt_level, Real stop_time) {
  BL_PROFILE("CNS::computeInitialDt()");
  //amrex::Print() << " oo CNS::computeInitialDt  " << std::endl;

  Real dt0 = std::numeric_limits<Real>::max();
  Vector<GpuArray<Real,AMREX_SPACEDIM>> eigenvals_level;
  Vector<Real> CFL_level;
  eigenvals_level.resize(finest_level + 1);
  CFL_level.resize(finest_level + 1);

  // Compute max eigenvalues in all directions on each level
  for (int i = 0; i <= finest_level; i++) {
    eigenvals_level[i] = getLevel(i).maxEigen();
  }

  if (dt_dynamic) {
    // Dynamic dt
    for (int i = 0; i <= finest_level; i++) {
      const GpuArray<Real, AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
#if (AMREX_SPACEDIM == 1)
      dt_level[i] = cfl * dx[0] / eigenvals_level[i][0];
#elif (AMREX_SPACEDIM == 2)
      dt_level[i] =
          cfl * amrex::min( dx[0] / eigenvals_level[i][0],
                            dx[1] / eigenvals_level[i][1] );                                        
#else
      dt_level[i] =
          cfl * amrex::min(AMREX_D_DECL(dx[0] / eigenvals_level[i][0],
                                        dx[1] / eigenvals_level[i][1],
                                        dx[2] / eigenvals_level[i][2]));
#endif
    }
    // Absolute cap (cns.dt_max — disabled by default)
    for (int i = 0; i <= finest_level; i++) {
      dt_level[i] = std::min(dt_level[i], dt_max);
    }
    // Find min dt across all levels
    int nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt0 = std::min(dt0, nfactor * dt_level[i]);
    }
  }
  else {
    // If constant dt
    dt0 = dt_constant;
  }
  // Set dt for all levels
  int nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    const GpuArray<Real,AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
    nfactor *= n_cycle[i];
    dt_level[i] = dt0 / nfactor;
#if (AMREX_SPACEDIM == 1)
    CFL_level[i] = dt_level[i] * eigenvals_level[i][0] / dx[0];
#elif (AMREX_SPACEDIM == 2)
     CFL_level[i] =
        dt_level[i] * amrex::max( eigenvals_level[i][0] / dx[0],
                                  eigenvals_level[i][1] / dx[1]);
#else
    CFL_level[i] =
        dt_level[i] * amrex::max(AMREX_D_DECL(eigenvals_level[i][0] / dx[0],
                                              eigenvals_level[i][1] / dx[1],
                                              eigenvals_level[i][2] / dx[2]));
#endif                                          
  }
  // Print
  if (ParallelDescriptor::IOProcessor()) {
    for (int i = 0; i <= finest_level; i++) {
      printf("[computeInitialDt] Level %d, Max CFL= %f   Max eigenvalues =( ",i,CFL_level[i]);
      for (int dim=0; dim < AMREX_SPACEDIM; dim++)
      { printf(" %f ",eigenvals_level[i][dim]);}
      printf(" ) \n");
      
    }
  }
}

// Called at the end of a coarse grid timecycle or after regrid, to compute the
// dt (time step) for all levels, for the next step.
// Output dt_level
void CNS::computeNewDt(int finest_level, int sub_cycle, Vector<int> &n_cycle,
                       const Vector<IntVect> &ref_ratio, Vector<Real> &dt_min,
                       Vector<Real> &dt_level, Real stop_time,
                       int post_regrid_flag) {
  BL_PROFILE("CNS::computeNewDt()");

  //amrex::Print() << " oo CNS::computeNewDt " << std::endl;

  Real dt0 = std::numeric_limits<Real>::max();
  Vector<GpuArray<Real,AMREX_SPACEDIM>> eigenvals_level;
  Vector<Real> CFL_level;
  eigenvals_level.resize(finest_level + 1);
  CFL_level.resize(finest_level + 1);

  // Compute max eigenvalues in all directions on each level
  for (int i = 0; i <= finest_level; i++) {
    eigenvals_level[i] = getLevel(i).maxEigen();
  }

  if (dt_dynamic) {
    const bool post_regrid = post_regrid_flag == 1;
    const bool bypass_restart_growth_limit =
        restart_first_dt_from_cfl_control.consume(post_regrid);
    const Real restart_old_dt0 = dt_level[0];

    // Estimate timestep across all points, levels, and procs
    for (int i = 0; i <= finest_level; i++) {
      const GpuArray<Real,AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
#if (AMREX_SPACEDIM == 1)
      dt_min[i] = cfl * dx[0] / eigenvals_level[i][0];
#elif (AMREX_SPACEDIM ==2) 
      dt_min[i] = cfl * amrex::min( dx[0] / eigenvals_level[i][0],
                                    dx[1] / eigenvals_level[i][1] );   
#else     
      dt_min[i] = cfl * amrex::min(AMREX_D_DECL(dx[0] / eigenvals_level[i][0],
                                                dx[1] / eigenvals_level[i][1],
                                                dx[2] / eigenvals_level[i][2]));
#endif                                          
    };

    // Post-regrid estimates never grow.  Normal estimates retain the legacy
    // 1.1x cap except for the one opt-in call armed by level-0 post_restart().
    static constexpr Real change_max = Real(1.1);
    for (int i = 0; i <= finest_level; i++) {
      dt_min[i] = cerisse::time_step::limit_cfl_dt(
          dt_min[i], dt_level[i], change_max, post_regrid,
          bypass_restart_growth_limit);
    }
    // Absolute cap (cns.dt_max — disabled by default)
    for (int i = 0; i <= finest_level; i++) {
      dt_min[i] = std::min(dt_min[i], dt_max);
    }
    // Find the minimum over all levels
    int nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt0 = std::min(dt0, nfactor * dt_min[i]);
    }
    const Real cfl_selected_dt0 = dt0;
    // Limit dt0 by the value of stop_time.
    const Real eps = 0.001_rt * dt0; ////////// PARAMETER /////////////
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0_rt) {
      if ((cur_time + dt0) > (stop_time - eps)) {
        dt0 = stop_time - cur_time;
      }
    }
    // Set dt at all levels
    nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt_level[i] = dt0 / nfactor;
    }
    if (bypass_restart_growth_limit) {
      amrex::Print()
          << "[CNS-Dt] restart-first CFL reset consumed: old_coarse_dt="
          << restart_old_dt0 << " cfl_candidate_dt=" << cfl_selected_dt0
          << " selected_dt=" << dt0 << "\n";
    }
  } else {
    // If constant dt
    dt0 = dt_constant;
  }
  // Set dt for all levels
  int nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    const GpuArray<Real,AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
    nfactor *= n_cycle[i];
    dt_level[i] = dt0 / nfactor;
#if (AMREX_SPACEDIM == 1)
    CFL_level[i] = dt_level[i] * eigenvals_level[i][0] / dx[0];
#elif (AMREX_SPACEDIM ==2) 
    CFL_level[i] =
        dt_level[i] * amrex::max(eigenvals_level[i][0] / dx[0],
                                 eigenvals_level[i][1] / dx[1]);    
#else    
    CFL_level[i] =
        dt_level[i] * amrex::max(AMREX_D_DECL(eigenvals_level[i][0] / dx[0],
                                              eigenvals_level[i][1] / dx[1],
                                              eigenvals_level[i][2] / dx[2]));
#endif                                          
  }
}

// Returns maximum eigenvalue in each direction
 [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> CNS::maxEigen() {
  BL_PROFILE("CNS::maxEigen()");

  PROB::ProbClosures const *d_cls = d_prob_closures;

  // Get multifabs
  MultiFab& consmf = get_new_data(State_Type);

  GpuArray<Real,AMREX_SPACEDIM> h_max_eigenvals;

  // Use ReduceOps for proper GPU-parallel reduction (replaces the original
  // AsyncArray approach which had a data race on the device max-reduction).
#if (AMREX_SPACEDIM == 1)
  ReduceOps<ReduceOpMax> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (MFIter mfi(consmf, false); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    const Array4<Real>& cons = consmf.array(mfi);
    reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      GpuArray<int, 3> vdir0 = {1, 0, 0};
      auto temp0 = d_cls->cons2eigenvals(i, j, k, cons, vdir0);
      Real maxe0 = Real(0.0);
      for (int iw = 0; iw < PROB::ProbClosures::NWAVES; iw++)
        maxe0 = amrex::max(maxe0, std::abs(temp0[iw]));
      return {maxe0};
    });
  }
  auto hv = reduce_data.value(reduce_op);
  h_max_eigenvals[0] = amrex::get<0>(hv);

#elif (AMREX_SPACEDIM == 2)
  ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_GPIBM
  // IBM-aware CFL: skip solid cells when computing max eigenvalue.
  // Solid cells don't participate in time integration (their RHS is zeroed
  // in compute_rhs), so their density/pressure values — which come from
  // extrapolation, not from the flow — should not constrain the timestep.
  // Without this, a solid cell with extrapolated ρ near zero produces a
  // huge sound speed → dt → 0 → simulation stalls.
  auto& ib_mf = *IBM::ib.bmf_a[level];
#endif

  for (MFIter mfi(consmf, false); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    const Array4<Real>& cons = consmf.array(mfi);
#ifdef AMREX_USE_GPIBM
    const auto& ibm = ib_mf.const_array(mfi);
#endif
    reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
#ifdef AMREX_USE_GPIBM
      // Skip solid cells — they don't evolve and shouldn't constrain CFL
      if (ibm(i,j,k,0) != 0) return {Real(0.0), Real(0.0)};
#endif
      GpuArray<int, 3> vdir0 = {1, 0, 0};
      GpuArray<int, 3> vdir1 = {0, 1, 0};
      auto temp0 = d_cls->cons2eigenvals(i, j, k, cons, vdir0);
      auto temp1 = d_cls->cons2eigenvals(i, j, k, cons, vdir1);
      Real maxe0 = Real(0.0), maxe1 = Real(0.0);
      for (int iw = 0; iw < PROB::ProbClosures::NWAVES; iw++) {
        maxe0 = amrex::max(maxe0, std::abs(temp0[iw]));
        maxe1 = amrex::max(maxe1, std::abs(temp1[iw]));
      }
      return {maxe0, maxe1};
    });
  }
  auto hv = reduce_data.value(reduce_op);
  h_max_eigenvals[0] = amrex::get<0>(hv);
  h_max_eigenvals[1] = amrex::get<1>(hv);

#else // 3D
  ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_GPIBM
  auto& ib_mf = *IBM::ib.bmf_a[level];
#endif

  for (MFIter mfi(consmf, false); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    const Array4<Real>& cons = consmf.array(mfi);
#ifdef AMREX_USE_GPIBM
    const auto& ibm = ib_mf.const_array(mfi);
#endif
    reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
#ifdef AMREX_USE_GPIBM
      if (ibm(i,j,k,0) != 0) return {Real(0.0), Real(0.0), Real(0.0)};
#endif
      GpuArray<int, 3> vdir0 = {1, 0, 0};
      GpuArray<int, 3> vdir1 = {0, 1, 0};
      GpuArray<int, 3> vdir2 = {0, 0, 1};
      auto temp0 = d_cls->cons2eigenvals(i, j, k, cons, vdir0);
      auto temp1 = d_cls->cons2eigenvals(i, j, k, cons, vdir1);
      auto temp2 = d_cls->cons2eigenvals(i, j, k, cons, vdir2);
      Real maxe0 = Real(0.0), maxe1 = Real(0.0), maxe2 = Real(0.0);
      for (int iw = 0; iw < PROB::ProbClosures::NWAVES; iw++) {
        maxe0 = amrex::max(maxe0, std::abs(temp0[iw]));
        maxe1 = amrex::max(maxe1, std::abs(temp1[iw]));
        maxe2 = amrex::max(maxe2, std::abs(temp2[iw]));
      }
      return {maxe0, maxe1, maxe2};
    });
  }
  auto hv = reduce_data.value(reduce_op);
  h_max_eigenvals[0] = amrex::get<0>(hv);
  h_max_eigenvals[1] = amrex::get<1>(hv);
  h_max_eigenvals[2] = amrex::get<2>(hv);
#endif

  // Communicate across processors
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    ParallelDescriptor::ReduceRealMax(h_max_eigenvals[idir]);
  }

  return h_max_eigenvals;
}

#ifdef AMREX_USE_GPIBM
void CNS::refreshAcceptedSharedGPState(Real time)
{
  BL_PROFILE("CNS::refreshAcceptedSharedGPState");

  MultiFab& accepted = get_new_data(State_Type);
  const int ncons = d_prob_closures->NCONS;
  const int conservative_nghost = d_prob_closures->NGHOST;
  const int primitive_nghost = IBM::ib.volumeInterpolationNghost(level);

  MultiFab conservative(
      grids, dmap, ncons, conservative_nghost,
      MFInfo().SetArena(The_Async_Arena()), Factory());
  FillPatch(
      *this, conservative, conservative_nghost, time, State_Type, 0, ncons);

#if (AMREX_SPACEDIM < 3)
  conservative.setVal(
      Real(0.0), PROB::ProbClosures::UMZ, 1, conservative.nGrow());
#endif

  const PROB::ProbClosures& cls_h = *h_prob_closures;
  const PROB::ProbClosures* cls_d = d_prob_closures;
  MultiFab primitives(
      grids, dmap, cls_h.NPRIM, primitive_nghost,
      MFInfo().SetArena(The_Async_Arena()));
  for (MFIter mfi(conservative, false); mfi.isValid(); ++mfi) {
    cls_h.cons2prims(
        mfi, conservative.array(mfi), primitives.array(mfi));
  }

#ifdef CNS_USE_FSI
  PROB::Motion::sim_time = time;
#endif

  const bool rz_annular_gp = IBM::ib.rzAnnularCellAverageEnabled(level);
  if (rz_annular_gp) {
    IBM::ib.computeAllGPsRZAnnular(
        primitives, conservative, cls_d, level);
  } else {
    IBM::ib.computeAllGPs(primitives, cls_d, level);
  }

  auto& markers = *IBM::ib.bmf_a[level];
  for (MFIter mfi(accepted, false); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();
    const auto state = accepted.array(mfi);
    const auto prims = primitives.const_array(mfi);
    const auto annular_state = conservative.const_array(mfi);
    const auto marker = markers.const_array(mfi);
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (marker(i, j, k, 1) != 0) {
        if (rz_annular_gp) {
          for (int n = 0; n < PROB::ProbClosures::NCONS; ++n) {
            state(i, j, k, n) = annular_state(i, j, k, n);
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
}
#endif

void CNS::post_timestep(int /* iteration*/) {
  BL_PROFILE("post_timestep");
  //amrex::Print() << " oo CNS::post_timestep " << std::endl;

  // The register contains time-integrated, area-weighted Euler+diffusive face
  // fluxes accumulated with the final RK quadrature weights.  Apply the
  // coarse-fine correction before average-down so the synchronized composite
  // state retains the conservative interface balance.
  if (do_reflux && level < parent->finestLevel()) {
    MultiFab& S = get_new_data(State_Type);
    CNS& fine_level = getLevel(level + 1);
    // Use the metric cell volumes: the constant-volume Reflux overload divides
    // by dx*dy, which mis-scales every RZ correction by the 2*pi*r volume
    // factor. GetVolume is exact in Cartesian too (dx*dy), so this is
    // behavior-preserving outside RZ.
    MultiFab volume(grids, dmap, 1, 0);
    geom.GetVolume(volume);
    fine_level.flux_reg->Reflux(
        S, volume, Real(1.0), 0, 0, PROB::ProbClosures::NCONS, geom);
  }

  if (level < parent->finestLevel()) {
    avgDown();
#ifdef AMREX_USE_GPIBM
    // Reflux and average-down occur after advance() has published the
    // end-of-step GP state. Reconstruct only the derived shared GP cells so
    // the synchronized accepted state is a valid future coarse FillPatch
    // source. Active fluid cells retain the conservative AMR correction.
    refreshAcceptedSharedGPState(state[State_Type].curTime());
#endif
  }

#ifdef AMREX_USE_GPIBM
  if (ibm_positivity_flux_limiter) {
    const MultiFab& synchronized_state = get_new_data(State_Type);
    const auto synchronized_stats =
        IBM::positivity::collectStageAdmissibility(
            synchronized_state, level, *h_prob_closures);
    if (!synchronized_stats.admissible()) {
      IBM::positivity::printStageAdmissibility(
          "post-reflux-average-down", synchronized_stats);
      amrex::Abort(
          "IBM positivity limiter: AMR reflux/average-down produced an "
          "inadmissible synchronized state");
    }
    if (ibm_positivity_flux_limiter_verbose) {
      amrex::Print()
          << "[IBM-Positivity-Limiter] level=" << level
          << " post_reflux_average_down_audit=PASS\n";
    }
  }
#endif

  // Record time statistics
  if (record_probe) {
    recordTimeProbe();    
  }

  // Record statistics
  if (record_stats) {
    time_stat_level[level] += parent->dtLevel(level);
    computeStats();
  }


  
    
}

void CNS::postCoarseTimeStep(Real time) {

  // amrex::Print() << " oo CNS::postCoarseTimeStep " << std::endl;
  amrex::ignore_unused(time);

#ifdef AMREX_USE_GPIBM
  // Surface fields are also required for FSI loads when file output is off.
  const int istep = parent->levelSteps(0);
  const bool surface_output_due = plot_surf && (istep % surf_int == 0);
#ifdef CNS_PROBLEM_IBM_POST_COARSE_DIAGNOSTIC
  const bool problem_surface_diagnostic_due =
      PROB::ibmPostCoarseDiagnosticDue(istep);
#else
  const bool problem_surface_diagnostic_due = false;
#endif
#ifdef CNS_USE_FSI
  const bool surface_compute_due = true;
#else
  const bool surface_compute_due =
      surface_output_due || problem_surface_diagnostic_due;
#endif
  Real surface_reconstruction_wall_s = Real(0.0);
  if (surface_compute_due) {
      const Real surface_start = amrex::second();
      for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
          dynamic_cast<CNS&>(parent->getLevel(lev)).writeSurfFile(
              !surface_output_due);
      }
      surface_reconstruction_wall_s = amrex::second() - surface_start;
  }

#ifdef CNS_PROBLEM_IBM_POST_COARSE_DIAGNOSTIC
  if (problem_surface_diagnostic_due) {
      PROB::ibmPostCoarseDiagnostic(
          IBM::ib, istep, time, parent->dtLevel(0),
          surface_reconstruction_wall_s);
  }
#endif

#ifdef CNS_USE_FSI
  // Calculate and print FSI loads and properties
  {
    auto& ib = IBM::ib;
    if (ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "\n=== FSI Loads (Step " << istep
                       << ", Time " << time << ") ===\n";
    }
    for (int i = 0; i < ib.ngeom; ++i) {
        auto props = FSI::RigidBodyProperties::readOrCompute(ib.geom_a[i], i);
        auto loads = FSI::Kinematics::computeLoads(i, props.xcenter);
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Geometry " << i << ":\n"
                           << "  Mass: " << props.mass << "\n"
                           << "  Center of Mass: " << props.xcenter << "\n"
                           << "  Inertia Tensor:\n";
            for (int r = 0; r < 3; ++r) {
                amrex::Print() << "    [ " << props.inertia[r][0] << ", "
                               << props.inertia[r][1] << ", "
                               << props.inertia[r][2] << " ]\n";
            }
            amrex::Print() << "  Fluid Force: " << loads.force << "\n"
                           << "  Fluid Moment (about CM): " << loads.moment << "\n"
                           << "----------------------------------------\n";
        }
    }
  }
#endif  // CNS_USE_FSI
#endif  // AMREX_USE_GPIBM

  if (verbose && ((this->nStep() % nstep_screen_output) == 0)) {
    printTotal();
  }
}
// -----------------------------------------------------------------------------

// Gridding -------------------------------------------------------------------
// Called for each level from 0,1...nlevs-1

void CNS::post_regrid(int lbase, int new_finest) {

#ifdef AMREX_USE_GPIBM
  rebuildIBM();

  // Interior-solid conservative values are opaque storage.  The IBM flux and
  // viscous operators mask them, and moving-body solid-to-fluid transitions
  // are repaired by the level-wide Jacobi solve in fixExposedCells().  Do not
  // flood-fill them here: this hook runs before state ghosts are guaranteed to
  // be current, so a per-FAB fill is decomposition-dependent and can import
  // uninitialised ghost values into valid cells.

#endif  // AMREX_USE_GPIBM

#ifdef CNS_USE_EB
  EBM::eb.destroy_mf(level);
  EBM::eb.build_mf(grids, dmap, level);

  // update volfrac and relevant EB data
  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  EBM::eb.volmf_a[level]  = &(ebfactory.getVolFrac()); 
  EBM::eb.normmcf_a[level] = &(ebfactory.getBndryNormal());
  EBM::eb.areamcf_a[level] = ebfactory.getAreaFrac();  
  EBM::eb.ebflags_a[level] = &(ebfactory.getMultiEBCellFlagFab());
  EBM::eb.bcareamcf_a[level] = &(ebfactory.getBndryArea());
  EBM::eb.bndrycent_a[level] = &(ebfactory.getBndryCent());
  EBM::eb.volcent_a[level] = &(ebfactory.getCentroid());
  

  // Level mask for redistribution (stored as object not pointer)
  EBM::eb.level_mask_a[level].clear();
  EBM::eb.level_mask_a[level].define(grids, dmap, 1, 3);
  EBM::eb.level_mask_a[level].BuildMask(
        geom.Domain(), geom.periodicity(), CNSConstants::level_mask_covered,
        CNSConstants::level_mask_notcovered, CNSConstants::level_mask_physbnd,
        CNSConstants::level_mask_interior);

  // EBM::eb.facecent  = ebfactory.getFaceCent();

  // Calculate markers  
  EBM::eb.computeMarkers(level);
  
  EBM::eb.check_geometry(level);


#endif

  if (use_nscbc) {
    initialize_nscbc_ghost_state(state[State_Type].curTime());
  }

}

void CNS::errorEst(TagBoxArray &tags, int /*clearval*/, int tagval,
                   Real time, int /*n_error_buf*/, int /*ngrow*/) {

 // amrex::Print() << " oo CNS::errorEst " << std::endl;

  // MF without ghost points filled (why?)
  MultiFab sdata(get_new_data(State_Type).boxArray(),
                 get_new_data(State_Type).DistributionMap(), PROB::ProbClosures::NCONS, PROB::ProbClosures::NGHOST,
                 MFInfo(), Factory());

  // filling ghost points (copied from PeleC)
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this, sdata, PROB::ProbClosures::NGHOST, cur_time, State_Type, 0, PROB::ProbClosures::NCONS);
  const auto geomdata = geom.data();

  
  // fill ghost points of stats (SNM needed? )
  // if (compute_stats) {
  //   MultiFab sdata_stat(get_new_data(Stats_Type).boxArray(),
  //                get_new_data(Stats_Type).DistributionMap(), PROB::ProbClosures::NSTAT, PROB::ProbClosures::NGHOST,
  //                MFInfo(), Factory());
  //   FillPatch(*this, sdata_stat, PROB::ProbClosures::NGHOST, cur_time, Stats_Type, 0, PROB::ProbClosures::NSTAT);
  
  // }

#ifdef AMREX_USE_GPIBM
  // call function from cns_prob
  auto &ibdata = (*IBM::ib.bmf_a[level]);
#endif
  for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    auto const &tagfab = tags.array(mfi);
    auto const &sdatafab = sdata.array(mfi);
#ifdef AMREX_USE_GPIBM
    auto const &ibfab = ibdata.array(mfi);
#endif
    int lev = level;
    int nt_lev = nStep();
    PROB::ProbParm const *lprobparm = d_prob_parm;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_GPIBM
      // call function from cns_prob
      user_tagging(i, j, k, nt_lev, tagfab, sdatafab, ibfab, geomdata,
                   *lprobparm, lev);
#else
      user_tagging(i, j, k, nt_lev, tagfab, sdatafab, geomdata ,*lprobparm, lev);
#endif
    });
  }

#ifdef AMREX_USE_GPIBM
  // Keep the next level's complete image-point support on same-level valid
  // cells. Build a one-cell interface mask and dilate it separably instead of
  // scanning an O(radius^D) neighbourhood for every coarse cell. Tags are only
  // added, so problem-specific refinement remains authoritative elsewhere.
  if (level < parent->maxLevel() && IBM::ib.amrSupportBufferEnabled()) {
    BL_PROFILE("CNS::errorEst::IBMAMRSupportBuffer");
    const IntVect support_radius = IBM::ib.amrSupportTagRadius(level);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        ibdata.nGrowVect().allGE(support_radius),
        "IBM marker halo is smaller than the AMR support-tag radius");

    ibdata.FillBoundary(geom.periodicity());
    iMultiFab band_a(tags.boxArray(), tags.DistributionMap(), 1,
                     support_radius);
    iMultiFab band_b(tags.boxArray(), tags.DistributionMap(), 1,
                     support_radius);
    band_a.setVal(0);
    band_b.setVal(0);

    const Box domain = geom.Domain();
    GpuArray<int, AMREX_SPACEDIM> periodic;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      periodic[d] = geom.isPeriodic(d) ? 1 : 0;
    }
    const int support_lev = level;
    PROB::ProbParm const* support_probparm = d_prob_parm;

    // Mark both sides of every marker transition. A radius-one cube preserves
    // the old Chebyshev-neighbourhood semantics, including diagonal cuts.
    for (MFIter mfi(band_a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& ibfab = ibdata.const_array(mfi);
      auto const& interface = band_a.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bool has_fluid = false;
        bool has_solid = false;
#if (AMREX_SPACEDIM == 3)
        for (int kk = -1; kk <= 1; ++kk) {
#endif
          for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
              const int ni = i + ii;
              const int nj = j + jj;
#if (AMREX_SPACEDIM == 3)
              const int nk = k + kk;
#else
              const int nk = k;
#endif
              bool available =
                  (periodic[0] ||
                   (ni >= domain.smallEnd(0) && ni <= domain.bigEnd(0))) &&
                  (periodic[1] ||
                   (nj >= domain.smallEnd(1) && nj <= domain.bigEnd(1)));
#if (AMREX_SPACEDIM == 3)
              available = available &&
                  (periodic[2] ||
                   (nk >= domain.smallEnd(2) && nk <= domain.bigEnd(2)));
#endif
              if (!available) continue;
              const bool solid = ibfab(ni, nj, nk, 0) != 0;
              has_solid = has_solid || solid;
              has_fluid = has_fluid || !solid;
            }
          }
#if (AMREX_SPACEDIM == 3)
        }
#endif
        // Some geometries need a spatially varying finest IBM level. Detect
        // the problem hook directly so a copied build file cannot silently
        // disable the case's GP-support contract. The hook limits only the
        // baseline interface seed; the interpolation-order-derived dilation
        // and strict support audit remain active. Cases without the hook retain
        // the original whole-interface behaviour.
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
#if (AMREX_SPACEDIM == 3)
        const Real z = prob_lo[2] + (Real(k) + Real(0.5)) * dx[2];
#else
        const Real z = Real(0.0);
#endif
        const bool support_refine_allowed =
            problem_ibm_support_refine_allowed(
            x, y, z, support_lev, *support_probparm);
        interface(i, j, k) =
            (has_solid && has_fluid && support_refine_allowed) ? 1 : 0;
      });
    }
    band_a.FillBoundary(geom.periodicity());

    iMultiFab* source = &band_a;
    iMultiFab* destination = &band_b;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const int dilation = amrex::max(support_radius[dir] - 1, 0);
      destination->setVal(0);
      for (MFIter mfi(*destination, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& src = source->const_array(mfi);
        auto const& dst = destination->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          int marked = 0;
          for (int offset = -dilation; offset <= dilation; ++offset) {
            int ii = i;
            int jj = j;
            int kk = k;
            if (dir == 0) ii += offset;
            if (dir == 1) jj += offset;
#if (AMREX_SPACEDIM == 3)
            if (dir == 2) kk += offset;
#endif
            marked = amrex::max(marked, src(ii, jj, kk, 0));
          }
          dst(i, j, k, 0) = marked;
        });
      }
      destination->FillBoundary(geom.periodicity());
      std::swap(source, destination);
    }

    const auto set_tag = static_cast<TagBox::TagType>(tagval);
    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& tagfab = tags.array(mfi);
      auto const& band = source->const_array(mfi);
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (band(i, j, k, 0) != 0) tagfab(i, j, k) = set_tag;
      });
    }
  }
#endif
}

void CNS::post_restart() {

  if (level == 0) {
    restart_first_dt_from_cfl_control.arm(restart_first_dt_from_cfl);
    if (restart_first_dt_from_cfl) {
      amrex::Print()
          << "[CNS-Dt] armed restart-first CFL reset for the first normal "
             "computeNewDt call\n";
    }
  }

// recreate markers
amrex::Print() << " recreate markers " << std::endl;

  // Restart constructs levels through the default constructor and restores
  // AmrLevel state afterward, so the parameterized constructor never defines
  // this level's coarse-fine register.  Reflux was historically disabled and
  // hid the null/empty register; with conservative AMR flux accumulation it
  // must be rebuilt from the restored hierarchy before the first advance.
  if (do_reflux && level > 0) {
    flux_reg = std::make_unique<FluxRegister>(
        grids, dmap, crse_ratio, level, PROB::ProbClosures::NCONS);
  } else {
    flux_reg.reset();
  }


#ifdef AMREX_USE_GPIBM
  rebuildIBM();
#endif

#ifdef CNS_USE_EB
  EBM::eb.destroy_mf(level);
  EBM::eb.build_mf(grids, dmap, level);

  // update volfrac and relevant EB data
  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  EBM::eb.volmf_a[level]  = &(ebfactory.getVolFrac()); 
  EBM::eb.normmcf_a[level] = &(ebfactory.getBndryNormal());
  EBM::eb.areamcf_a[level] = ebfactory.getAreaFrac();  
  EBM::eb.ebflags_a[level] = &(ebfactory.getMultiEBCellFlagFab());
  EBM::eb.bcareamcf_a[level] = &(ebfactory.getBndryArea());
  EBM::eb.bndrycent_a[level] = &(ebfactory.getBndryCent());
  EBM::eb.volcent_a[level] = &(ebfactory.getCentroid());
  

  // Level mask for redistribution (stored as object not pointer)
  EBM::eb.level_mask_a[level].clear();
  EBM::eb.level_mask_a[level].define(grids, dmap, 1, 3);
  EBM::eb.level_mask_a[level].BuildMask(
        geom.Domain(), geom.periodicity(), CNSConstants::level_mask_covered,
        CNSConstants::level_mask_notcovered, CNSConstants::level_mask_physbnd,
        CNSConstants::level_mask_interior);

  // EBM::eb.facecent  = ebfactory.getFaceCent();

  // Calculate markers  
  EBM::eb.computeMarkers(level);

  EBM::eb.check_geometry(level);

#endif

  // A legacy checkpoint has no Stats_Type, whereas a statistics-enabled
  // checkpoint contains running means that must survive continuation.
  if (compute_stats) {
    const bool checkpoint_has_stats = state[Stats_Type].hasNewData();
    const Real cur_time = state[State_Type].curTime();
    if (!checkpoint_has_stats) {
      const Real dt_old   = cur_time - state[State_Type].prevTime();
      state[Stats_Type].define(geom.Domain(), grids, dmap,
                               desc_lst[Stats_Type], cur_time, dt_old,
                               Factory());
      setupStats();
    }
    time_stat_level[level] = checkpoint_has_stats
        ? amrex::max(cur_time - stats_start_time, Real(0.0))
        : Real(0.0);
  }

  // Set up diagnostics after restart
  if (record_probe) {
    setupTimeProbe();
  }

  if (use_nscbc) {
    initialize_nscbc_ghost_state(state[State_Type].curTime());
  }

}

void CNS::checkPointPost(const std::string& dir, std::ostream& os)
{
  AmrLevel::checkPointPost(dir, os);
#ifdef AMREX_USE_GPIBM
#else
  amrex::ignore_unused(dir);
#endif
}

void CNS::set_state_in_checkpoint(Vector<int>& state_in_checkpoint) {
  // This is only called when the checkpoint has fewer state types than
  // the current code.  Mark Stats_Type as absent so AMReX skips reading it.
  if (compute_stats) {
    state_in_checkpoint[Stats_Type] = 0;
  }
}

// 
void CNS::avgDown() {
  BL_PROFILE("CNS::avgDown()");

  if (level == parent->finestLevel()) return;

  auto &fine_lev = getLevel(level + 1);

  MultiFab &S_crse = get_new_data(State_Type);
  MultiFab &S_fine = fine_lev.get_new_data(State_Type);

  amrex::average_down(
      S_fine, S_crse, fine_lev.geom, geom, 0, S_fine.nComp(),
      parent->refRatio(level));
  Gpu::streamSynchronize();  // ensure GPU average_down is complete before CPU read

#ifdef AMREX_USE_GPIBM
  // IBM-aware avgDown correction.
  // Standard average_down has already run above. Now correct coarse cells
  // that overlap MIXED fine regions (some solid + some fluid fine sub-cells)
  // by re-averaging using only fluid fine cells. Pure-fluid coarse cells
  // keep the standard average; pure-solid coarse cells keep their own value.
  if (IBM::ib.bmf_a[level + 1] != nullptr)
  {
    auto& fine_ibmf = *IBM::ib.bmf_a[level + 1];
    const int nc = S_fine.nComp();
    const IntVect rr = parent->refRatio(level);

    // Build coarsened box array on fine's DistributionMap
    BoxArray cba = S_fine.boxArray();
    cba.coarsen(rr);

    // Use ncons+1 components: first ncons are data, last is validity flag.
    // The geometry-aware AMReX average_down above uses metric cell volumes in
    // 2-D.  Preserve that RZ semantics when the IBM correction excludes solid
    // children: equal dr cells are unequal annular volumes (most visibly the
    // 1:3 pair adjacent to the axis).  Cartesian children remain on the legacy
    // arithmetic path so this correction is behavior-preserving there.
    const int ncp1 = nc + 1;
    MultiFab S_corr(cba, S_fine.DistributionMap(), ncp1, 0);
    S_corr.setVal(Real(0.0));

    const bool metric_weighted = fine_lev.geom.IsRZ();
    std::unique_ptr<MultiFab> fine_volume;
    if (metric_weighted) {
      fine_volume = std::make_unique<MultiFab>(
          S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
      fine_lev.geom.GetVolume(*fine_volume, S_fine.boxArray(),
                              S_fine.DistributionMap(), 0);
    }

    // Populate S_corr with fluid-only averages for MIXED cells
    // Note: use GPU ParallelFor because in CUDA builds, MultiFab data
    // lives in device memory and can't be accessed via raw CPU loops.
    for (MFIter fmfi(S_fine, false); fmfi.isValid(); ++fmfi) {
      const Box fbx = fmfi.tilebox();
      const Box cbx = amrex::coarsen(fbx, rr);
      auto const fine = S_fine.const_array(fmfi);
      auto const fmk  = fine_ibmf.const_array(fmfi);
      auto const corr = S_corr.array(fmfi);
      Array4<const Real> fvol;
      if (metric_weighted) fvol = fine_volume->const_array(fmfi);
      const int ncomp = nc;
      const IntVect ratio = rr;

      amrex::ParallelFor(cbx,
      [=] AMREX_GPU_DEVICE (int ci, int cj, int ck) noexcept
      {
        int n_fluid = 0, n_solid = 0;
        Real sum[PROB::ProbClosures::NCONS] = {};
        Real fluid_volume = Real(0.0);
#if (AMREX_SPACEDIM == 2)
        for (int fj = cj*ratio[1]; fj < (cj+1)*ratio[1]; ++fj)
        for (int fi = ci*ratio[0]; fi < (ci+1)*ratio[0]; ++fi) {
          if (fmk(fi,fj,0,0) == 0) {
            const Real weight = metric_weighted
                ? fvol(fi,fj,0,0) : Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
              sum[n] += weight * fine(fi,fj,0,n);
            }
            fluid_volume += weight;
            ++n_fluid;
          } else ++n_solid;
        }
#else
        for (int fk = ck*ratio[2]; fk < (ck+1)*ratio[2]; ++fk)
        for (int fj = cj*ratio[1]; fj < (cj+1)*ratio[1]; ++fj)
        for (int fi = ci*ratio[0]; fi < (ci+1)*ratio[0]; ++fi) {
          if (fmk(fi,fj,fk,0) == 0) {
            for (int n = 0; n < ncomp; ++n) sum[n] += fine(fi,fj,fk,n);
            fluid_volume += Real(1.0);
            ++n_fluid;
          } else ++n_solid;
        }
#endif
        // Only overwrite MIXED cells
        if (n_fluid > 0 && n_solid > 0) {
          const Real inv = Real(1.0) / fluid_volume;
          for (int n = 0; n < ncomp; ++n) corr(ci,cj,ck,n) = sum[n] * inv;
          corr(ci,cj,ck,ncomp) = Real(1.0);  // mark valid
        }
      });
    }

    // ParallelCopy S_corr to a MultiFab aligned with S_crse
    MultiFab S_corr_aligned(S_crse.boxArray(), S_crse.DistributionMap(), ncp1, 0);
    S_corr_aligned.setVal(Real(0.0));
    S_corr_aligned.ParallelCopy(S_corr, 0, 0, ncp1);

    // Apply corrections: where validity flag == 1, overwrite S_crse
    for (MFIter mfi(S_crse, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& crse = S_crse.array(mfi);
      auto const& corr = S_corr_aligned.const_array(mfi);
      const int ncomp = nc;

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (corr(i,j,k,ncomp) > Real(0.5)) {
          for (int n = 0; n < ncomp; ++n)
            crse(i,j,k,n) = corr(i,j,k,n);
        }
      });
    }
  }
#endif

  if (compute_stats) {
    MultiFab &Sstat_crse = get_new_data(Stats_Type);
    MultiFab &Sstat_fine = fine_lev.get_new_data(Stats_Type);
    amrex::average_down(Sstat_fine, Sstat_crse, fine_lev.geom, geom, 0, Sstat_fine.nComp(),
                      parent->refRatio(level));
  }

}

void CNS::printTotal() const {
  // Get conservatives multifab
  const MultiFab& consmf = get_new_data(State_Type);

  // Volume-weighted integral of conserved variables (works for Cartesian and RZ)
  MultiFab volume(consmf.boxArray(), consmf.DistributionMap(), 1, 0);
  geom.GetVolume(volume, consmf.boxArray(), consmf.DistributionMap(), 0);

  std::array<Real, PROB::ProbClosures::NCONS> tot{};
  for (int comp = 0; comp < PROB::ProbClosures::NCONS; ++comp) {
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    auto const& a = consmf.const_arrays();
    auto const& v = volume.const_arrays();
    reduce_op.eval(consmf, IntVect(0), reduce_data,
                   [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept -> Real {
                     return a[box_no](i, j, k, comp) * v[box_no](i, j, k);
                   });
    Gpu::streamSynchronize();
    auto const& hv = reduce_data.value(reduce_op);
    tot[comp] = amrex::get<0>(hv);
  }

  // Communicate across processors
  ParallelDescriptor::ReduceRealSum(tot.data(), PROB::ProbClosures::NCONS,
                                    ParallelDescriptor::IOProcessorNumber());
  // Print
  Vector<std::string> names= PROB::ProbClosures::get_cons_vars_names();
  for (int comp = 0; comp < PROB::ProbClosures::NCONS; ++comp) {
    amrex::Print().SetPrecision(17) << "   Total " << names[comp] << " = " << tot[comp] << "\n";
  }
}

void CNS::variableCleanUp() {
  delete h_prob_closures;
  delete h_phys_bc;

#ifdef AMREX_USE_GPU
  The_Arena()->free(d_prob_closures);
  The_Arena()->free(d_phys_bc);
#endif
  desc_lst.clear();
  derive_lst.clear();
}

// Plotting
//------------------------------------------------------------------------------
void CNS::writePlotFile(const std::string &dir, std::ostream &os,
                        VisMF::How how) {
  int i, n;
  //
  // The list of indices of State to write to plotfile.
  // first component of pair is state_type,
  // second component of pair is component # within the state_type
  //
  std::vector<std::pair<int, int>> plot_var_map;
  for (int typ = 0; typ < desc_lst.size(); typ++) {
    for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
      if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
          desc_lst[typ].getType() == IndexType::TheCellType()) {
        plot_var_map.push_back(std::pair<int, int>(typ, comp));
      }
    }
  }

  int num_derive = 0;
  std::vector<std::string> derive_names;
  const std::list<DeriveRec> &dlist = derive_lst.dlist();
  for (auto const &d : dlist) {
    if (parent->isDerivePlotVar(d.name())) {
      derive_names.push_back(d.name());
      num_derive += d.numDerive();
    }
  }

  int n_data_items = plot_var_map.size() + num_derive;

//----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
  n_data_items += 2;
#endif
#ifdef CNS_USE_EB
  n_data_items += 1;
#endif
  //------------------------------------------------------------------------------

  // get the time from the first State_Type
  // if the State_Type is ::Interval, this will get t^{n+1/2} instead of t^n
  Real cur_time = state[0].curTime();

  if (level == 0 && ParallelDescriptor::IOProcessor()) {
    //
    // The first thing we write out is the plotfile type.
    //
    os << thePlotFileType() << '\n';

    if (n_data_items == 0)
      amrex::Error("Must specify at least one valid data item to plot");

    os << n_data_items << '\n';

    //
    // Names of variables
    //
    for (i = 0; i < static_cast<int>(plot_var_map.size()); i++) {
      int typ = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      os << desc_lst[typ].name(comp) << '\n';
    }

    // derived
    for (auto const &dname : derive_names) {
      const DeriveRec *rec = derive_lst.get(dname);
      for (i = 0; i < rec->numDerive(); ++i) {
        os << rec->variableName(i) << '\n';
      }
    }

    //----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
    os << "sld\n";
    os << "ghs\n";
#endif
#ifdef CNS_USE_EB
    os << "vfrac\n";
#endif

    //------------------------------------------------------------------------------

    os << AMREX_SPACEDIM << '\n';
    os << parent->cumTime() << '\n';
    int f_lev = parent->finestLevel();
    os << f_lev << '\n';
    for (i = 0; i < AMREX_SPACEDIM; i++) os << Geom().ProbLo(i) << ' ';
    os << '\n';
    for (i = 0; i < AMREX_SPACEDIM; i++) os << Geom().ProbHi(i) << ' ';
    os << '\n';
    for (i = 0; i < f_lev; i++) os << parent->refRatio(i)[0] << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) os << parent->Geom(i).Domain() << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) os << parent->levelSteps(i) << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) {
      for (int k = 0; k < AMREX_SPACEDIM; k++)
        os << parent->Geom(i).CellSize()[k] << ' ';
      os << '\n';
    }
    os << (int)Geom().Coord() << '\n';
    os << "0\n";  // Write bndry data.
  }
  // Build the directory to hold the MultiFab at this level.
  // The name is relative to the directory containing the Header file.
  //
  static const std::string BaseName = "/Cell";
  char buf[64];
  snprintf(buf, sizeof buf, "Level_%d", level);
  std::string sLevel = buf;
  //
  // Now for the full pathname of that directory.
  //
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.size() - 1] != '/') {
    FullPath += '/';
  }
  FullPath += sLevel;
  //
  // Only the I/O processor makes the directory if it doesn't already exist.
  //
  if (!levelDirectoryCreated) {
    if (ParallelDescriptor::IOProcessor()) {
      if (!amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier();
  }

  if (ParallelDescriptor::IOProcessor()) {
    os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
    os << parent->levelSteps(level) << '\n';

    for (i = 0; i < grids.size(); ++i) {
      RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
      for (n = 0; n < AMREX_SPACEDIM; n++)
        os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
    }
    //
    // The full relative pathname of the MultiFabs at this level.
    // The name is relative to the Header file containing this name.
    // It's the name that gets written into the Header.
    //
    if (n_data_items > 0) {
      std::string PathNameInHeader = sLevel;
      PathNameInHeader += BaseName;
      os << PathNameInHeader << '\n';
    }
  }
  //
  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  int cnt = 0;
  const int nGrow = 0;
  MultiFab plotMF(grids, dmap, n_data_items, nGrow, MFInfo(), Factory());
  MultiFab *this_dat = 0;
  //
  // Cull data from state variables -- use no ghost cells.
  //
  for (i = 0; i < static_cast<int>(plot_var_map.size()); i++) {
    int typ = plot_var_map[i].first;
    int comp = plot_var_map[i].second;
    this_dat = &state[typ].newData();
    MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
    cnt++;
  }

  // derived
  if (derive_names.size() > 0) {
    for (auto const &dname : derive_names) {
      derive(dname, cur_time, plotMF, cnt);
      cnt += derive_lst.get(dname)->numDerive();
    }  //exit(1);
  }

  //------------------------------------------------------------------------------
  // additional plotting ...
  //------------------------------------------------------------------------------
  
#ifdef AMREX_USE_GPIBM
  plotMF.setVal(0.0_rt, cnt, 2, nGrow);
  IBM::ib.bmf_a[level]->copytoRealMF(plotMF, 0, cnt);
  cnt+=2;
#endif

#ifdef CNS_USE_EB
  plotMF.setVal(0.0_rt, cnt, 0, nGrow); 
 // EBM::eb.bmf_a[level]->copytoRealMF(plotMF, 0, cnt);  // boolean 
  const MultiFab *vfrac = EBM::eb.volmf_a[level];
  MultiFab::Copy(plotMF, *vfrac, 0, cnt, 1, 0);
  cnt++;
#endif

  //------------------------------------------------------------------------------

  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  if (AsyncOut::UseAsyncOut()) {
    VisMF::AsyncWrite(plotMF, TheFullPath);
  } else {
    VisMF::Write(plotMF, TheFullPath, how, true);
  }

  levelDirectoryCreated = false;  // ---- now that the plotfile is finished
}

// This is called once per level on write timestep.
// surf_int must be the same as plot_int
void CNS::writePlotFilePost(const std::string &dir, std::ostream &os) {

#if AMREX_USE_GPIBM
  // writeSurfFile();
#endif

}

// this subroutine is called from the main loop 
// should be called per level
#if AMREX_USE_GPIBM

void CNS::rebuildIBM(bool rebuild_surface) {
  const Real ibm_rebuild_timing_start = IBM::ib.beginPerformanceTiming();
  IBM::ib.destroy_mf(level);
  IBM::ib.build_mf(grids, dmap, level);
  IBM::ib.computeMarkers(level);
  IBM::ib.initialiseGPs(level);
  // Surface indices depend on all levels having valid bmf_a, so we can
  // only rebuild them once the entire regrid cascade is complete (i.e.,
  // when the finest level calls rebuildIBM).
  if (rebuild_surface && level == parent->finestLevel()) {
     // --- Runaway-regrid guard (benefits ALL IBM+AMR cases) -------------------
     // A refinement criterion that is NOT grid-independent (e.g. a pure flow-
     // gradient tag whose stencil straddles the artificial IBM ghost/jet jump)
     // can make the AMR regrid cascade never converge: every regrid produces a
     // different fine-grid layout, so this re-surfacing fires endlessly within a
     // single coarse step with no time advance -- the run wedges silently (CPU
     // pinned, only plt00000, no STEP). Rather than spin forever, fail loud with
     // an actionable diagnostic (consistent with the solver's no-silent-fix /
     // abort-with-diagnostics policy).
     static int s_last_step = -1;
     static int s_rebuilds_this_step = 0;
     const int cur_step = parent->levelSteps(0);
     if (cur_step != s_last_step) { s_last_step = cur_step; s_rebuilds_this_step = 0; }
     const int max_rebuilds = 32 * (parent->maxLevel() + 1);
     if (++s_rebuilds_this_step > max_rebuilds) {
        amrex::Abort(
          "IBM+AMR: surface re-build fired " + std::to_string(s_rebuilds_this_step) +
          " times within coarse step " + std::to_string(cur_step) +
          " (cap " + std::to_string(max_rebuilds) + ") -- the AMR regrid cascade is "
          "NOT converging. Cause: a non-grid-independent refinement tag near the IBM "
          "surface (a flow-gradient stencil that straddles the IBM ghost/jet jump "
          "oscillates as FillPatch re-interpolates each regrid). Fix: add a "
          "GEOMETRY-ANCHORED (grid-independent) refinement region near the IBM "
          "surface in user_tagging, and skip the gradient tag where the stencil "
          "touches solid/ghost cells.");
     }
     for (int lev = parent->finestLevel(); lev >= 0; --lev) {
        IBM::ib.computeSurfIndices(lev);
     }
  }
  IBM::ib.finishPerformanceTiming("rebuild_ibm_total", level,
                                  ibm_rebuild_timing_start);
}

void CNS::writeSurfFile(bool force_compute, bool force_output) {
      
  // calculate and  write surface data  
  int istep = parent->levelSteps(0);

  const bool output_due =
      plot_surf && (force_output || istep % surf_int == 0);


  if (output_due || force_compute)  {
     
    MultiFab& state_data = get_new_data(State_Type);

    int ncons = CNS::d_prob_closures->NCONS;
    int nghost= CNS::d_prob_closures->NGHOST;

    Real time = parent->cumTime();

    if (output_due && this->level == parent->finestLevel()) {
      Print() << "Computing surface properties ";
      Print() << " at time= " << time << " and step= " << istep << std::endl;
    }
    
    FillPatch(*this, state_data, nghost, time, State_Type, 0, ncons);

    MultiFab& Sdata = state_data;

    // In 2D, prob_initdata / bcnormal may never write UMZ (the index exists
    // even in 2D), so valid cells and physical-BC ghosts can hold garbage.
    // compute_rhs() sanitises its own working copy each stage, but this
    // export path reads Sdata directly — zero UMZ here too so cons2prims
    // does not turn garbage z-kinetic-energy into floor-clipped P/T on the
    // exported surface. Mirrors compute_rhs.cpp.
#if (AMREX_SPACEDIM < 3)
    Sdata.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, Sdata.nGrow());
#endif

    // Convert conservative to primitive variables for surface interpolation
    int nprim = PROB::ProbClosures::NPRIM;
    MultiFab prims_mf(Sdata.boxArray(), Sdata.DistributionMap(),
                      nprim, IBM::ib.interpolationMarkerNghost(this->level),
                      MFInfo().SetArena(The_Async_Arena()));
    for (MFIter mfi(Sdata, false); mfi.isValid(); ++mfi) {
      CNS::h_prob_closures->cons2prims(mfi, Sdata.array(mfi), prims_mf.array(mfi));
    }

    const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
    const PROB::ProbClosures* cls_h = CNS::h_prob_closures; 
    amrex::ignore_unused(cls_h);

    if (IBM::ib.rzAnnularCellAverageEnabled(this->level)) {
      IBM::ib.recoverRZCentrePrimitives(
          Sdata, prims_mf, cls_d, this->level,
          /*active_fluid_only=*/true);
      if constexpr (PROB::ProbIB::rz_bic_support_centre_recovery) {
        IBM::ib.recoverRZBICSupportCentrePrimitives(
            Sdata, prims_mf, cls_d, this->level);
      }
      IBM::ib.fillRZAxisPrimitiveParity(prims_mf, this->level);
    }

    MultiFab* surface_prims = &prims_mf;

    bool use_surface_cell_average_recovery =
        IBM::ib.nsSurfaceCellAverageRecoveryEnabled();
    MultiFab* surface_conservative_averages = nullptr;
    std::unique_ptr<MultiFab> surface_conservative_storage;
    if (use_surface_cell_average_recovery) {
      const int surface_ngrow =
          IBM::ib.surfaceInterpolationNghost(this->level);
      if (Sdata.nGrow() >= surface_ngrow) {
        surface_conservative_averages = &Sdata;
      } else {
        surface_conservative_storage = std::make_unique<MultiFab>(
            Sdata.boxArray(), Sdata.DistributionMap(), ncons, surface_ngrow,
            MFInfo().SetArena(The_Async_Arena()));
        FillPatch(
            *this, *surface_conservative_storage, surface_ngrow, time,
            State_Type, 0, ncons);
#if (AMREX_SPACEDIM < 3)
        surface_conservative_storage->setVal(
            Real(0.0), PROB::ProbClosures::UMZ, 1, surface_ngrow);
#endif
        surface_conservative_averages =
            surface_conservative_storage.get();
      }
    }

    IBM::ib.computeSURFs(
        prims_mf, surface_conservative_averages, cls_d,
        this->level); // computed from low to high

    // Only gather and write on the finest level to ensure all levels are processed
    if (output_due && this->level == parent->finestLevel()){
      // collect data to rank 0
      IBM::ib.gatherSurfData(); 

      if (amrex::ParallelDescriptor::IOProcessor()){
        IBM::ib.plotSURF(time, istep, surf_filename); 
      } 
    }

  }
}
#endif
