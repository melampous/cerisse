#ifndef CNS_ISENTROPIC_VORTEX_PROB_H_
#define CNS_ISENTROPIC_VORTEX_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

// C1.6 convecting isentropic vortex parameters.  The default Mach number and
// vortex strength are the "fast vortex" values used for the compact thesis
// test; both can be overridden from the input file.
struct ProbParm {
  Real mach = Real(0.5);
  Real beta = Real(0.2);
  Real p_inf = Real(1.0e5);
  Real T_inf = Real(300.0);
  Real vortex_radius = Real(0.005);
  Real x_centre = Real(0.05);
  Real y_centre = Real(0.05);
  Real flow_angle_deg = Real(0.0);

  int static_refinement = 0;
  Real refine_x_lo = Real(0.05);
  Real refine_x_hi = Real(0.10);

  ProbParm() {
    ParmParse pp("prob");
    pp.query("mach", mach);
    pp.query("beta", beta);
    pp.query("p_inf", p_inf);
    pp.query("T_inf", T_inf);
    pp.query("vortex_radius", vortex_radius);
    pp.query("x_centre", x_centre);
    pp.query("y_centre", y_centre);
    pp.query("flow_angle_deg", flow_angle_deg);
    pp.query("static_refinement", static_refinement);
    pp.query("refine_x_lo", refine_x_lo);
    pp.query("refine_x_hi", refine_x_hi);
  }
};

using ProbClosures =
    closures_dt<indicies_t, calorifically_perfect_gas_t<indicies_t>>;

// The two Skew executables use the same order-four central operator.  Separate
// compile-time parameter types keep the undamped operator and the full JST
// method distinct in archived executables and source snapshots.
struct SkewCentralO4Parm {
  static constexpr bool dissipation = false;
  static constexpr int order = 4;
  static constexpr Real C2skew = Real(0.0);
  static constexpr Real C4skew = Real(0.0);
};

struct SkewJstO4Parm {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = Real(1.5);
  static constexpr Real C4skew = Real(0.016);
};

#ifndef IV_EULER_SCHEME_ID
#define IV_EULER_SCHEME_ID 0
#endif

#if (IV_EULER_SCHEME_ID == 0)
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
#elif (IV_EULER_SCHEME_ID == 1)
using ProbEuler = weno_t<ReconScheme::Teno5, ProbClosures>;
#elif (IV_EULER_SCHEME_ID == 2)
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
#elif (IV_EULER_SCHEME_ID == 3)
using ProbEuler = afd_hllc_teno5_t<ProbClosures>;
#elif (IV_EULER_SCHEME_ID == 4)
using ProbEuler = skew_t<SkewCentralO4Parm, ProbClosures>;
#elif (IV_EULER_SCHEME_ID == 5)
using ProbEuler = skew_t<SkewJstO4Parm, ProbClosures>;
#else
#error "Unknown IV_EULER_SCHEME_ID"
#endif

using ProbRHS = rhs_dt<ProbEuler, no_diffusive_t, no_source_t>;

// Time integration is selected explicitly in the archived input files.
inline void inputs() {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_euler_state(Real *state, const Real rho, const Real u, const Real v,
                const Real p, ProbClosures const &cls) {
  state[cls.URHO] = rho;
  state[cls.UMX] = rho * u;
  state[cls.UMY] = rho * v;
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = p / (cls.gamma - Real(1.0))
                 + Real(0.5) * rho * (u * u + v * v);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &pparm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];

  const Real length_x = prob_hi[0] - prob_lo[0];
  const Real length_y = prob_hi[1] - prob_lo[1];
  Real x_distance = x - pparm.x_centre;
  Real y_distance = y - pparm.y_centre;
  if (x_distance > Real(0.5) * length_x) x_distance -= length_x;
  if (x_distance < -Real(0.5) * length_x) x_distance += length_x;
  if (y_distance > Real(0.5) * length_y) y_distance -= length_y;
  if (y_distance < -Real(0.5) * length_y) y_distance += length_y;

  const Real xr = x_distance / pparm.vortex_radius;
  const Real yr = y_distance / pparm.vortex_radius;
  const Real radius_sq = xr * xr + yr * yr;
  const Real envelope = std::exp(-Real(0.5) * radius_sq);

  const Real sound_speed =
      std::sqrt(cls.gamma * cls.Rspec * pparm.T_inf);
  const Real mean_speed = pparm.mach * sound_speed;
  const Real flow_angle =
      pparm.flow_angle_deg * std::acos(Real(-1.0)) / Real(180.0);
  const Real u_inf = mean_speed * std::cos(flow_angle);
  const Real v_inf = mean_speed * std::sin(flow_angle);
  const Real circulation_speed = mean_speed * pparm.beta;

  const Real u = u_inf - circulation_speed * yr * envelope;
  const Real v = v_inf + circulation_speed * xr * envelope;
  const Real temperature = pparm.T_inf
      - Real(0.5) * circulation_speed * circulation_speed / cls.cp
            * std::exp(-radius_sq);

  const Real rho_inf = pparm.p_inf / (cls.Rspec * pparm.T_inf);
  const Real rho = rho_inf
      * std::pow(temperature / pparm.T_inf,
                 Real(1.0) / (cls.gamma - Real(1.0)));
  const Real p = rho * cls.Rspec * temperature;

  Real local_state[ProbClosures::NCONS] = {Real(0.0)};
  set_euler_state(local_state, rho, u, v, p, cls);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = local_state[n];
  }
}

// All physical directions are periodic for this problem.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS], Real[ProbClosures::NCONS],
         int, int, Real, GeometryData const &, ProbClosures const &,
         ProbParm const &) {
  amrex::Abort("isentropic_vortex: bcnormal called on a periodic problem");
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int, int, int, const auto &, const auto &, const ProbParm &,
            ProbClosures const &, auto const) {}

// Optional stationary refined band.  It is used only by inputs.amr to expose
// the vortex to two coarse--fine crossings; the uniform-grid verification does
// not use this path.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int, auto &tagfab, const auto &,
             const auto &geomdata, const ProbParm &pparm, int level) {
  if (pparm.static_refinement == 0 || level != 0) {
    return;
  }
  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  if (x >= pparm.refine_x_lo && x < pparm.refine_x_hi) {
    tagfab(i, j, k) = true;
  }
}

} // namespace PROB

#endif
