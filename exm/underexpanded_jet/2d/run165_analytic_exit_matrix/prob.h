#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <Closures.h>
#include <Constants.h>
#include <RHS.h>
#include <bc_types.h>

#include <cmath>

using namespace amrex;
using namespace universal_constants;

#if (AMREX_SPACEDIM != 2)
#error "Run165 analytic-exit matrix requires DIM=2"
#endif

#ifndef RUN165_JET_ENABLED
#define RUN165_JET_ENABLED 1
#endif

#ifndef RUN165_JET_SCALE
#define RUN165_JET_SCALE 1.0
#endif

#ifndef RUN165_USE_NS
#define RUN165_USE_NS 0
#endif

#ifndef RUN165_USE_LES
#define RUN165_USE_LES 0
#endif

#ifndef RUN165_LES_CS
#define RUN165_LES_CS 0.16
#endif

#ifndef RUN165_EULER_SCHEME
#define RUN165_EULER_SCHEME 0
#endif

#ifndef RUN165_SKEW_C2
#define RUN165_SKEW_C2 1.5
#endif

#ifndef RUN165_SKEW_C4
#define RUN165_SKEW_C4 0.016
#endif

static_assert(RUN165_JET_ENABLED == 0 || RUN165_JET_ENABLED == 1,
              "RUN165_JET_ENABLED must be zero or one");
static_assert((RUN165_JET_SCALE) > 0.0 && (RUN165_JET_SCALE) <= 1.0,
              "RUN165_JET_SCALE must lie in (0,1]");
static_assert(RUN165_USE_NS == 0 || RUN165_USE_NS == 1,
              "RUN165_USE_NS must be zero or one");
static_assert(RUN165_USE_LES == 0 || RUN165_USE_LES == 1,
              "RUN165_USE_LES must be zero or one");
static_assert(!(RUN165_USE_NS && RUN165_USE_LES),
              "select either laminar NS or LES, not both mode switches");
static_assert((RUN165_LES_CS) > 0.0 && (RUN165_LES_CS) <= 0.3,
              "RUN165_LES_CS must lie in (0,0.3]");
static_assert(RUN165_EULER_SCHEME >= 0 && RUN165_EULER_SCHEME <= 5,
              "RUN165_EULER_SCHEME must be 0 (LLF-WENO-Z5), 1 (Skew-JST), "
              "2 (HLLC-MUSCL), 3 (AFD-HLLC-WENO-Z5), "
              "4 (LLF-TENO5), or 5 (AFD-HLLC-TENO5)");
static_assert((RUN165_SKEW_C2) >= 0.0 && (RUN165_SKEW_C4) >= 0.0,
              "Run165 Skew-JST coefficients must be non-negative");

namespace PROB {

// Exact external state used by the IBM Run165 case.
static constexpr bool jet_enabled = RUN165_JET_ENABLED != 0;
static constexpr bool use_ns = RUN165_USE_NS != 0;
static constexpr bool use_les = RUN165_USE_LES != 0;
static constexpr bool use_skew = RUN165_EULER_SCHEME == 1;
static constexpr bool use_hllc_muscl = RUN165_EULER_SCHEME == 2;
static constexpr bool use_afd_hllc_weno = RUN165_EULER_SCHEME == 3;
static constexpr bool use_llf_teno = RUN165_EULER_SCHEME == 4;
static constexpr bool use_afd_hllc_teno = RUN165_EULER_SCHEME == 5;
static constexpr Real jet_scale = Real(RUN165_JET_SCALE);
static constexpr Real gamma_air = 1.4_rt;
static constexpr Real molecular_weight_air = 0.02896_rt;
static constexpr Real R_air = gas_constant / molecular_weight_air;
static constexpr Real p_inf = 534.5806702256042_rt;
static constexpr Real T_inf = 64.72986748216108_rt;
static constexpr Real rho_inf = 0.028765564804356037_rt;
static constexpr Real u_inf = 741.9796978999816_rt;
static constexpr Real e_inf = 46460.12288142623_rt;

// Isentropic virtual-exit state implied by the exact reconstructed N1 area
// ratio and the same reservoir state used by the IBM throat homotopy.
static constexpr Real x_exit = 0.002441268805570889_rt;
static constexpr Real r_exit = 0.00636778_rt;
static constexpr Real exit_mach = 2.953619236_rt;
static constexpr Real exit_temperature = 125.9120371_rt;
static constexpr Real exit_pressure_full = 120537.5495_rt;
static constexpr Real exit_density_full = 3.334413689_rt;
static constexpr Real exit_velocity = -664.4617178_rt;
static constexpr Real exit_pressure = jet_scale * exit_pressure_full;
static constexpr Real exit_density = jet_scale * exit_density_full;

struct ProbParm {
  // All cells are covered through this AMR level.  This is intentionally a
  // uniform-grid convergence block: it has no coarse/fine interface at all.
  int target_level = 1;
  // Optional startup regularisation.  A positive value ramps the imposed
  // aperture state from the exact freestream conservative state to the exact
  // Run165 exit conservative state with a C1 smoothstep.  The final boundary
  // state and mass flow are unchanged.
  Real jet_ramp_time = 0.0_rt;

  ProbParm()
  {
    ParmParse pp("prob");
    pp.query("target_level", target_level);
    pp.query("jet_ramp_time", jet_ramp_time);
    if (target_level < 1 || target_level > 3) {
      amrex::Abort("prob.target_level must lie in [1,3]");
    }
    if (!(jet_ramp_time >= 0.0_rt) || !std::isfinite(jet_ramp_time)) {
      amrex::Abort("prob.jet_ramp_time must be finite and non-negative");
    }
  }
};

struct methodparm_t {
  static constexpr int order = 2;
  static constexpr Real viscosity = 0.0_rt;
  static constexpr Real conductivity = 0.0_rt;
  static constexpr bool use_LES = RUN165_USE_LES != 0;
};

struct lesparm_t {
  static constexpr int order = 2;
  static constexpr Real Pr_o_Prsgs = 0.8_rt;
  static constexpr Real Scsgs = 0.7_rt;
  static constexpr Real Cs = Real(RUN165_LES_CS);
  static constexpr Real CI = 0.08_rt;
  static constexpr bool fixDelta = false;
  static constexpr Real Delta = 0.0_rt;
};

struct skewparm_t {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr int sensor_power = 1;
  static constexpr Real C2skew = Real(RUN165_SKEW_C2);
  static constexpr Real C4skew = Real(RUN165_SKEW_C4);
};

struct gasparm_t {
  static constexpr Real gamma = gamma_air;
  static constexpr Real molecular_weight = molecular_weight_air;
};

#if RUN165_USE_LES
using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                perfect_gas_t<gasparm_t, indicies_t>,
                Smagorinsky_t<lesparm_t, indicies_t>>;
#elif RUN165_USE_NS
using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                perfect_gas_t<gasparm_t, indicies_t>>;
#else
using ProbClosures =
    closures_dt<indicies_t, transport_const_t<methodparm_t>,
                perfect_gas_t<gasparm_t, indicies_t>>;
#endif
using GlobalBC = manual_bc_t<ProbClosures>;
#if RUN165_EULER_SCHEME == 1
using ProbEuler = skew_t<skewparm_t, ProbClosures>;
#elif RUN165_EULER_SCHEME == 2
using ProbEuler = riemann_t<false, ProbClosures>;
#elif RUN165_EULER_SCHEME == 3
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
#elif RUN165_EULER_SCHEME == 4
using ProbEuler = weno_t<ReconScheme::Teno5, ProbClosures>;
#elif RUN165_EULER_SCHEME == 5
using ProbEuler = afd_hllc_teno5_t<ProbClosures>;
#else
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
#endif
#if RUN165_USE_LES
using ProbRHS =
    rhs_dt<ProbEuler,
           viscousLES_t<methodparm_t, ProbClosures>, no_source_t>;
#elif RUN165_USE_NS
using ProbRHS =
    rhs_dt<ProbEuler,
           viscous_t<methodparm_t, ProbClosures>, no_source_t>;
#else
using ProbRHS =
    rhs_dt<ProbEuler, no_diffusive_t,
           no_source_t>;
#endif

inline void inputs()
{
  const ProbParm pparm;
  const Real mdot = exit_density * (-exit_velocity) *
                    3.14159265358979323846_rt * r_exit * r_exit;
  amrex::Print() << " ****** Run165 RZ analytic virtual-exit matrix: jet="
                 << (jet_enabled ? "ON" : "OFF")
                 << "; Euler operator="
                 << (use_skew
                         ? "Skew-JST-O4"
                         : (use_hllc_muscl
                                ? "HLLC-MUSCL-O2"
                                : (use_afd_hllc_weno
                                       ? "AFD-HLLC-WENO-Z5"
                                       : (use_llf_teno
                                              ? "LLF-TENO5"
                                              : (use_afd_hllc_teno
                                                     ? "AFD-HLLC-TENO5"
                                                     : "LLF-WENO-Z5")))))
                 << "; equations="
                 << (use_les ? "Navier-Stokes+Smagorinsky-SGS"
                             : (use_ns ? "laminar Navier-Stokes" : "Euler"))
                 << "; LES=" << (use_les ? "ON" : "OFF");
  if constexpr (use_les) {
    amrex::Print() << "; Cs=" << lesparm_t::Cs;
  }
  if constexpr (use_skew) {
    amrex::Print() << "; C2=" << skewparm_t::C2skew
                   << "; C4=" << skewparm_t::C4skew;
  }
  amrex::Print()
      << " ******\n"
      << " target uniform level=L" << pparm.target_level
      << "; coordinates=(r,x_TSP); no IBM/solid geometry; jet_ramp_time="
      << pparm.jet_ramp_time << " s\n"
      << " freestream: p=" << p_inf << " Pa T=" << T_inf
      << " K rho=" << rho_inf << " kg/m3 u_x=" << u_inf << " m/s\n"
      << " exit: x=" << x_exit << " m r=" << r_exit
      << " m M=" << exit_mach << " T=" << exit_temperature
      << " K scale=" << jet_scale << " p=" << exit_pressure
      << " Pa rho=" << exit_density << " kg/m3 u_x="
      << exit_velocity << " m/s mdot=" << mdot << " kg/s\n";
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const&, ProbClosures const& cls,
                   ProbParm const&)
{
  amrex::ignore_unused(i, j, k);
  state(i, j, k, cls.URHO) = rho_inf;
  state(i, j, k, cls.UMX) = 0.0_rt;
  state(i, j, k, cls.UMY) = rho_inf * u_inf;
  state(i, j, k, cls.UMZ) = 0.0_rt;
  state(i, j, k, cls.UET) =
      rho_inf * e_inf + 0.5_rt * rho_inf * u_inf * u_inf;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real x[AMREX_SPACEDIM], Real,
              const Real s_int[ProbClosures::NCONS],
              const Real[ProbClosures::NCONS],
              Real s_ext[ProbClosures::NCONS], int idir, int sgn, Real time,
              GeometryData const&, ProbClosures const& cls,
              ProbParm const& pparm)
{
  for (int n = 0; n < ProbClosures::NCONS; ++n) s_ext[n] = s_int[n];
  if (idir != 1) return;

  // Axial-low is the exact M=4.6 Run165 inflow, directed toward +x.
  if (sgn == 1) {
    s_ext[cls.URHO] = rho_inf;
    s_ext[cls.UMX] = 0.0_rt;
    s_ext[cls.UMY] = rho_inf * u_inf;
    s_ext[cls.UMZ] = 0.0_rt;
    s_ext[cls.UET] =
        rho_inf * e_inf + 0.5_rt * rho_inf * u_inf * u_inf;
    return;
  }

  // Axial-high is outflow except for the face-aligned supersonic jet disk.
  // With dx0=r_exit/3, the aperture contains exactly 6/12/24 complete
  // annuli on L1/L2/L3; no cell-centre area alias contaminates H effects.
  if constexpr (jet_enabled) {
    if (x[0] < r_exit) {
      Real ramp = 1.0_rt;
      if (pparm.jet_ramp_time > 0.0_rt) {
        const Real fraction = amrex::min(
            1.0_rt, amrex::max(0.0_rt, time / pparm.jet_ramp_time));
        ramp = fraction * fraction * (3.0_rt - 2.0_rt * fraction);
      }
      const Real rho_jet = exit_density;
      const Real mom_jet = exit_density * exit_velocity;
      const Real energy_jet =
          exit_pressure / (gamma_air - 1.0_rt) +
          0.5_rt * exit_density * exit_velocity * exit_velocity;
      const Real mom_inf = rho_inf * u_inf;
      const Real energy_inf =
          rho_inf * e_inf + 0.5_rt * rho_inf * u_inf * u_inf;
      s_ext[cls.URHO] = (1.0_rt - ramp) * rho_inf + ramp * rho_jet;
      s_ext[cls.UMX] = 0.0_rt;
      s_ext[cls.UMY] = (1.0_rt - ramp) * mom_inf + ramp * mom_jet;
      s_ext[cls.UMZ] = 0.0_rt;
      s_ext[cls.UET] =
          (1.0_rt - ramp) * energy_inf + ramp * energy_jet;
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int, auto& tagfab, const auto&,
                  const auto&, const ProbParm& pparm, int level)
{
  // Full-domain coverage through target_level is cheap in this deliberately
  // small matched-aperture block and removes AMR as a causal variable.
  tagfab(i, j, k) = level < pparm.target_level
                        ? amrex::TagBox::SET
                        : amrex::TagBox::CLEAR;
}

}  // namespace PROB

#endif
