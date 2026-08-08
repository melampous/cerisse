#include <CNS.h>
#include <CNS_K.h>
#include <EulerAdmissibleProlongation.H>
#include <prob.h>

#include <limits>
#include <memory>
#include <type_traits>

using namespace amrex;

int CNS::num_state_data_types = 0;

namespace {

template <typename T, typename = void>
struct inviscid_scheme_manifest {
  static constexpr const char* value = "unreported";
};

template <typename T>
struct inviscid_scheme_manifest<
    T, std::void_t<decltype(T::scheme_name)>> {
  static constexpr const char* value = T::scheme_name;
};

template <typename T, typename = void>
struct rz_paired_pressure_capability : std::false_type {};

template <typename T>
struct rz_paired_pressure_capability<
    T, std::void_t<decltype(T::rz_paired_pressure_flux_capable)>>
    : std::bool_constant<T::rz_paired_pressure_flux_capable> {};

template <typename T, typename = void>
struct local_llf_fallback_manifest {
  static void print() {}
};

template <typename T>
struct local_llf_fallback_manifest<
    T, std::void_t<decltype(T::print_local_llf_fallback_manifest())>> {
  static void print() { T::print_local_llf_fallback_manifest(); }
};

}  // namespace

static Box the_same_box(const Box& b) { return b; }
// static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

using BndryFunc = StateDescriptor::BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall, User, FarField
//  (0)       (1)     (2)      (3)       (4)       (5)         (6)   (7)
static int scalar_bc[] = {BCType::int_dir,      BCType::ext_dir,
                          BCType::foextrap,     BCType::reflect_even,
                          BCType::reflect_even, BCType::reflect_even,
                          BCType::ext_dir,      BCType::ext_dir};

static int norm_vel_bc[] = {BCType::int_dir,     BCType::ext_dir,
                            BCType::foextrap,    BCType::reflect_odd,
                            BCType::reflect_odd, BCType::reflect_odd,
                            BCType::ext_dir,     BCType::ext_dir};

static int tang_vel_bc[] = {BCType::int_dir,      BCType::ext_dir,
                            BCType::foextrap,     BCType::reflect_even,
                            BCType::reflect_even, BCType::reflect_odd,
                            BCType::ext_dir,      BCType::ext_dir};

// If AMReX returns a negative BC (e.g. -1), treat as Interior (0) (this is to avoid weird errors while using periodic BCs)
inline int safe_bc_index(int bc) {    
  if (bc < 0) return 0;
  if (bc > 7) amrex::Abort("Unsupported cns.lo_bc/hi_bc value");
  return bc;
}

// BC type for a STATISTICS component in one direction.
// Stats ghosts are only consumed by interpolation at regrid; they must never
// take the ext_dir/bcnormal path (bcnormal writes conservative-variable
// values, meaningless for stats). Mirror parity applies at Symmetry/SlipWall:
// odd for components with an odd number of u_d factors in direction d.
static int stat_bc_type(int physbc, bool odd) {
  switch (safe_bc_index(physbc)) {
    case 0: return BCType::int_dir;                                   // Interior
    case 3:                                                            // Symmetry
    case 4: return odd ? BCType::reflect_odd : BCType::reflect_even;   // SlipWall
    default: return BCType::foextrap;  // Inflow/Outflow/NoSlipWall/User/FarField
  }
}

static void set_scalar_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    bc.setLo(i, scalar_bc[safe_bc_index(lo_bc[i])]);
    bc.setHi(i, scalar_bc[safe_bc_index(hi_bc[i])]);
  }
}

static void set_x_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();

  bc.setLo(0, norm_vel_bc[safe_bc_index(lo_bc[0])]);
  bc.setHi(0, norm_vel_bc[safe_bc_index(hi_bc[0])]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_vel_bc[safe_bc_index(lo_bc[1])]);
  bc.setHi(1, tang_vel_bc[safe_bc_index(hi_bc[1])]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_vel_bc[safe_bc_index(lo_bc[2])]);
  bc.setHi(2, tang_vel_bc[safe_bc_index(hi_bc[2])]);
#endif
}

static void set_y_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();

  bc.setLo(0, tang_vel_bc[safe_bc_index(lo_bc[0])]);
  bc.setHi(0, tang_vel_bc[safe_bc_index(hi_bc[0])]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, norm_vel_bc[safe_bc_index(lo_bc[1])]);
  bc.setHi(1, norm_vel_bc[safe_bc_index(hi_bc[1])]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_vel_bc[safe_bc_index(lo_bc[2])]);
  bc.setHi(2, tang_vel_bc[safe_bc_index(hi_bc[2])]);
#endif
}

static void set_z_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();

  bc.setLo(0, tang_vel_bc[safe_bc_index(lo_bc[0])]);
  bc.setHi(0, tang_vel_bc[safe_bc_index(hi_bc[0])]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_vel_bc[safe_bc_index(lo_bc[1])]);
  bc.setHi(1, tang_vel_bc[safe_bc_index(hi_bc[1])]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, norm_vel_bc[safe_bc_index(lo_bc[2])]);
  bc.setHi(2, norm_vel_bc[safe_bc_index(hi_bc[2])]);
#endif

}

void CNS::variableSetUp() {

  //amrex::Print( ) << " oo CNS::variableSetUp " << std::endl; 

  // Closures and Problem structures (available on both CPU and GPU)
  CNS::h_prob_closures = new PROB::ProbClosures{};
  CNS::h_prob_parm = new PROB::ProbParm{};
  CNS::h_phys_bc = new BCRec{};
#ifdef AMREX_USE_GPU
  CNS::d_prob_closures =
      (PROB::ProbClosures*)The_Arena()->alloc(sizeof(PROB::ProbClosures));
  CNS::d_prob_parm =
      (PROB::ProbParm*)The_Arena()->alloc(sizeof(PROB::ProbParm));
  CNS::d_phys_bc = (BCRec*)The_Arena()->alloc(sizeof(BCRec));
#else
  CNS::d_prob_closures = h_prob_closures;
  CNS::d_prob_parm = h_prob_parm;
  CNS::d_phys_bc = h_phys_bc;
#endif

  // Read input parameters
  read_params();
  amrex::Print() << "[Numerics] inviscid_scheme="
                 << inviscid_scheme_manifest<PROB::ProbRHS>::value << '\n';
  local_llf_fallback_manifest<PROB::ProbRHS>::print();
  int removed_rz_companion_option = 0;
  if (ParmParse("cns").query(
          "rz_companion_pressure_flux", removed_rz_companion_option)) {
    amrex::Abort(
        "cns.rz_companion_pressure_flux has been removed; the paired "
        "radial-pressure operator is automatic for every compatible R-Z "
        "inviscid scheme, including pure shared-GP builds");
  }
#if !defined(CNS_USE_EB)
  if constexpr (rz_paired_pressure_capability<PROB::ProbRHS>::value) {
    amrex::Print()
        << "[Numerics] rz_paired_pressure_operator="
        << "paired_metric_advection_and_pressure_gradient_when_RZ\n"
        << "[Numerics] rz_axis_pressure_auxiliary="
        << "order_matched_parity_reconstruction"
        << " selector="
        << inviscid_scheme_manifest<PROB::ProbRHS>::value
#ifdef AMREX_USE_GPIBM
        << " stencil=marker_aware_shared_gp"
#else
        << " stencil=all_fluid"
#endif
        << '\n';
#ifdef AMREX_USE_GPIBM
    if (CNS::ibm_positivity_flux_limiter) {
      amrex::Print()
          << "[Numerics] ibm_positivity_rz_operator="
          << "shared_face_dual_channel_metric_flux_and_pressure_gradient\n"
          << "[Numerics] ibm_positivity_rz_pressure_blend="
          << "same_theta_as_complete_flux\n"
          << "[Numerics] ibm_positivity_amr_reflux="
          << "multilevel_fail_closed_pending_invariant_domain_sync\n"
          << "[Numerics] ibm_positivity_active_rhs_hooks="
          << "forbidden_face_flux_operator_only\n";
    }
#endif
  }
#endif

  // Independent (solved) variables and their boundary condition types
  bool state_data_extrap = false;
  bool store_in_checkpoint = true;
  static std::unique_ptr<cerisse::amr::EulerAdmissibleConservativeLinear>
      admissible_linear_interp;
  Interpolater* state_interp = nullptr;
  const char* state_interp_name = nullptr;
  if (amr_state_interp == 1) {
    state_interp = &quartic_interp;
    state_interp_name =
        "AMReX conservative quartic (smooth-flow verification only)";
  } else if (amr_state_interp == 2) {
    state_interp = &lincc_interp;
    state_interp_name = "legacy limited conservative linear";
  } else {
    Real specific_internal_energy_floor = Real(0.0);
    Real internal_energy_density_floor =
        CNSConstants::min_press() / h_prob_closures->gamma_m1;
#if CLIP_TEMPERATURE_MIN
    // Match the strict post-step audit while leaving a small representable
    // interior margin. This prevents regrid interpolation from creating a
    // conservative state that cons2prims would silently temperature-clip.
    const Real thermodynamic_floor_scale =
        Real(1.0) + Real(1024.0) * std::numeric_limits<Real>::epsilon();
    specific_internal_energy_floor =
        h_prob_closures->cv * CNSConstants::min_temp() *
        thermodynamic_floor_scale;
    internal_energy_density_floor *= thermodynamic_floor_scale;
#endif
    if (!admissible_linear_interp) {
      admissible_linear_interp = std::make_unique<
          cerisse::amr::EulerAdmissibleConservativeLinear>(
          cerisse::amr::EulerStateLayout{
              PROB::ProbClosures::UMX,
              PROB::ProbClosures::UMY,
              PROB::ProbClosures::UMZ,
              PROB::ProbClosures::UET,
              PROB::ProbClosures::URHO,
              CNSConstants::small_rho(),
              internal_energy_density_floor,
              specific_internal_energy_floor});
    }
    state_interp = admissible_linear_interp.get();
    state_interp_name =
        "closure-admissibility-preserving limited conservative linear "
        "(parent-block conservative theta; configured rho-e floor)";
    amrex::Print()
        << "  AMR state interpolation floors: rho="
        << CNSConstants::small_rho()
        << " rhoe=" << internal_energy_density_floor
        << " ei=" << specific_internal_energy_floor << '\n';
  }
  amrex::Print() << "  AMR state interpolation = "
                 << state_interp_name << "\n";
  desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                         StateDescriptor::Point, h_prob_closures->NGHOST,
                         h_prob_closures->NCONS, state_interp,
                         state_data_extrap, store_in_checkpoint);
  // https://github.com/AMReX-Codes/amrex/issues/396

  // Statistical variables
  if (h_prob_closures->NSTAT > 0 ) {
    compute_stats=true;
    desc_lst.addDescriptor(Stats_Type, IndexType::TheCellType(),
                          StateDescriptor::Point, 0, h_prob_closures->NSTAT, &lincc_interp,
                          state_data_extrap, store_in_checkpoint);
  }                        
  
  // Physical boundary conditions ////////////////////////////////////////////
  Vector<BCRec> bcs(PROB::ProbClosures::NCONS);
  Vector<int> cons_vars_type = indicies_t::get_cons_vars_type();

  for (int cnt=0;cnt<h_prob_closures->NCONS;cnt++) {

    switch (cons_vars_type[cnt])
    {
      case 0:
        set_scalar_bc(bcs[cnt], h_phys_bc);
        break;
      case 1:
        set_x_vel_bc (bcs[cnt], h_phys_bc);
        break;
      case 2:  
        set_y_vel_bc (bcs[cnt], h_phys_bc); 
        break;
      case 3:  
        set_z_vel_bc (bcs[cnt], h_phys_bc);    
        break;
      default:
        std::cout << " error ... variableSetUp" << std::endl;
        exit(1); 
    }

  }

  // Boundary conditions
  StateDescriptor::BndryFunc bndryfunc(cns_bcfill);
  StateDescriptor::setBndryFuncThreadSafety(true);
  bndryfunc.setRunOnGPU(true);
  // applies bndry func to all variables in desc_lst starting from from 0.
  desc_lst.setComponent(State_Type, 0, PROB::ProbClosures::get_cons_vars_names(), bcs, bndryfunc);

  num_state_data_types = desc_lst.size();

  // SET-UP Stats Type 
  ////////////////////////////////////////////////////////////////////////////
  if (h_prob_closures->NSTAT > 0) { 


    int NSTAT        = h_prob_closures->NSTAT;         
    Vector<BCRec> stats_bcs(NSTAT);
    Vector<std::string>  stats_name(NSTAT);

    amrex::Print() << " Storing  Flow Statistics " << "\n";
    amrex::Print() << " NSTAT=" << NSTAT  << " record_stats=" << record_stats << "\n";

    int statv = -1;
    // names velocity
    if (h_prob_closures->record_velocity > 0) {      
      statv++; stats_name[statv] = "x_velocityMEAN";
#if AMREX_SPACEDIM >1    
      statv++; stats_name[statv] = "y_velocityMEAN";
#endif
#if AMREX_SPACEDIM == 3    
      statv++; stats_name[statv] = "z_velocityMEAN";
#endif
      statv++; stats_name[statv] = "x_velocitySQR";
#if AMREX_SPACEDIM >1    
      statv++; stats_name[statv] = "y_velocitySQR";
#endif
#if AMREX_SPACEDIM == 3    
      statv++; stats_name[statv] = "z_velocitySQR";
#endif
#if AMREX_SPACEDIM >1    
      statv++; stats_name[statv] = "xy_velocityMEAN";
#endif
#if AMREX_SPACEDIM == 3    
      statv++; stats_name[statv] = "xz_velocityMEAN";
      statv++; stats_name[statv] = "yz_velocityMEAN";    
#endif    
    }
    // names P,T,rho
    if (h_prob_closures->record_PTrho > 0) {
      statv++; stats_name[statv] = "pressureMEAN"; INDEX_THERM=statv;
      statv++; stats_name[statv] = "temperatureMEAN";
      statv++; stats_name[statv] = "DensityMEAN";
      statv++; stats_name[statv] = "pressureSQR";
      statv++; stats_name[statv] = "temperatureSQR";
      statv++; stats_name[statv] = "DensitySQR";            
    }
    // names species (TODO)
  
    // bc
    // Per-component reflection parity per direction: odd iff the statistic
    // contains an odd number of u_d factors (u_d MEAN and u_d*u_e cross
    // moments are odd; squares and P/T/rho stats are even).
    {
      constexpr int NS = PROB::ProbClosures::NSTAT;
      int stat_odd[NS][AMREX_SPACEDIM] = {};
      if (h_prob_closures->record_velocity > 0) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) stat_odd[d][d] = 1;   // u_d MEAN
#if (AMREX_SPACEDIM >= 2)
        stat_odd[2*AMREX_SPACEDIM][0] = 1;                              // xy cross
        stat_odd[2*AMREX_SPACEDIM][1] = 1;
#endif
#if (AMREX_SPACEDIM == 3)
        stat_odd[2*AMREX_SPACEDIM+1][0] = 1;                            // xz cross
        stat_odd[2*AMREX_SPACEDIM+1][2] = 1;
        stat_odd[2*AMREX_SPACEDIM+2][1] = 1;                            // yz cross
        stat_odd[2*AMREX_SPACEDIM+2][2] = 1;
#endif
      }
      for (statv = 0; statv < h_prob_closures->NSTAT; statv++) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          stats_bcs[statv].setLo(d, stat_bc_type(h_phys_bc->lo(d), stat_odd[statv][d] != 0));
          stats_bcs[statv].setHi(d, stat_bc_type(h_phys_bc->hi(d), stat_odd[statv][d] != 0));
        }
      }
    }
    StateDescriptor::BndryFunc bndryfuncstats( cns_bcfill);
    bndryfuncstats.setRunOnGPU(true);
    desc_lst.setComponent(Stats_Type, 0, stats_name, stats_bcs, bndryfuncstats);
  }
  //printf("num_state_data_types(2) %d \n",num_state_data_types);

  ////////////////////////////////////////////////////////////////////////////

  // Define derived quantities ///////////////////////////////////////////////
  // Pressure
  derive_lst.add("pressure", IndexType::TheCellType(), 1, derpres,
                 the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Temperature
  derive_lst.add("temperature", IndexType::TheCellType(), 1, dertemp,
                 the_same_box);
  derive_lst.addComponent("temperature", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Velocities
  derive_lst.add("velocity", IndexType::TheCellType(), amrex::SpaceDim,
                 {AMREX_D_DECL("x_velocity", "y_velocity", "z_velocity")},
                 dervel, the_same_box);
  derive_lst.addComponent("velocity", desc_lst, State_Type, 0,
                          h_prob_closures->NCONS);

  // Density
  derive_lst.add("density", IndexType::TheCellType(), 1, derdensity,
                 the_same_box);
  derive_lst.addComponent("density", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Kinetic energy
  derive_lst.add("kinetic_energy", IndexType::TheCellType(), 1, derkineticenergy,
                 the_same_box);
  derive_lst.addComponent("kinetic_energy", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Vorticity magnitude
  derive_lst.add("magvort", IndexType::TheCellType(), 1, dermagvort,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("magvort", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Enstrophy
  derive_lst.add("enstrophy", IndexType::TheCellType(), 1, derenstrophy,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("enstrophy", desc_lst, State_Type, 0, h_prob_closures->NCONS);
}
