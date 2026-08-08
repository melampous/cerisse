#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <Constants.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <ibm_walltypes.h>

#include <cmath>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

struct ProbParm
{
  Real pressure = Real(101325.0);
  Real temperature = Real(300.0);
  Real interface_z = Real(0.25);
  Real forced_velocity = Real(2500.0);
  Real pressure_trench_ratio = Real(1.0e-4);
  int initial_mode = 0;
  int refine_patch = 0;

  ProbParm()
  {
    ParmParse pp("prob");
    pp.query("initial_mode", initial_mode);
    pp.query("interface_z", interface_z);
    pp.query("forced_velocity", forced_velocity);
    pp.query("pressure_trench_ratio", pressure_trench_ratio);
    pp.query("refine_patch", refine_patch);
    AMREX_ALWAYS_ASSERT(initial_mode >= 0 && initial_mode <= 3);
    AMREX_ALWAYS_ASSERT(refine_patch == 0 || refine_patch == 1);
    AMREX_ALWAYS_ASSERT(pressure_trench_ratio > Real(0.0));
    if (initial_mode != 0) {
      amrex::Print()
          << "[TestManifest] forced_limiter_initial_mode=" << initial_mode
          << " interface_z=" << interface_z
          << " forced_velocity=" << forced_velocity
          << " pressure_trench_ratio=" << pressure_trench_ratio << '\n';
    }
  }
};

struct methodparm_t
{
  static constexpr int order = 2;
  static constexpr Real conductivity = Real(1.0);
  static constexpr Real viscosity = Real(1.0);
  static constexpr bool use_LES = false;
};

struct ibmparm_t
{
  static constexpr int interp_order = 2;
  static constexpr int extrap_order = 1;
  static constexpr Real alpha = Real(0.6);
  static constexpr int interp_order_surf = 2;
  static constexpr int extrap_order_surf = 1;
  static constexpr Real alpha_surf = Real(0.6);
  static constexpr int ghost_layers = 1;
  static constexpr bool interior_is_solid = true;

  static constexpr int support_visibility_mode = 2;
  static constexpr Real support_visibility_condition_max = Real(1.0e10);
  static constexpr bool expanded_bic_support = true;
  static constexpr ibm_shared_gp_reconstruction_t shared_gp_reconstruction =
      ibm_shared_gp_reconstruction_t::boundary_intercept_constrained;
  static constexpr bool cell_average_bic_normal_momentum = true;
  static constexpr bool one_sided_entropy_jet = true;
  static constexpr bool rz_annular_bic_cell_average = true;
  static constexpr bool rz_bic_support_centre_recovery = true;
  static constexpr bool rz_axis_regular_entropy_jet = true;
  static constexpr bool smooth_wall_geometry = true;
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_t::euler_slip_analytic_curvature;

  static constexpr Real circle_r = Real(0.55);
  static constexpr Real circle_z = Real(0.50);
  static constexpr Real circle_radius = Real(0.12);

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
  smooth_wall_boundary_intercept(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& raw_point,
      Array1D<Real, 0, AMREX_SPACEDIM - 1>& boundary_intercept) noexcept
  {
    const Real dr = raw_point(0) - circle_r;
    const Real dz = raw_point(1) - circle_z;
    const Real distance = std::sqrt(dr * dr + dz * dz);
    if (!(distance > Real(1.0e-14)) ||
        !amrex::Math::isfinite(distance)) {
      return false;
    }
    const Real scale = circle_radius / distance;
    boundary_intercept(0) = circle_r + scale * dr;
    boundary_intercept(1) = circle_z + scale * dz;
    return true;
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
  smooth_wall_local_frame(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& boundary_intercept,
      LocalFrame& frame) noexcept
  {
    const Real dr = boundary_intercept(0) - circle_r;
    const Real dz = boundary_intercept(1) - circle_z;
    const Real distance = std::sqrt(dr * dr + dz * dz);
    if (!(distance > Real(1.0e-14)) ||
        !amrex::Math::isfinite(distance)) {
      return false;
    }
    const Real inverse_distance = Real(1.0) / distance;
    frame.normal[0] = dr * inverse_distance;
    frame.normal[1] = dz * inverse_distance;
    frame.tangent1[0] = -frame.normal[1];
    frame.tangent1[1] = frame.normal[0];
    return true;
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ibm_shape_operator_t
  wall_shape_operator(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& xyz,
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& normal,
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tangent1,
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tangent2) noexcept
  {
    amrex::ignore_unused(xyz, normal, tangent1, tangent2);
    return {Real(1.0) / circle_radius, Real(0.0), Real(0.0), 1};
  }
};

using ProbClosures =
    closures_dt<indicies_t, transport_const_t<methodparm_t>,
                calorifically_perfect_gas_t<indicies_t>>;

struct SkewJstO4Parm
{
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = Real(1.5);
  static constexpr Real C4skew = Real(0.016);
};

using ProbRHS =
    rhs_dt<skew_t<SkewJstO4Parm, ProbClosures>, no_diffusive_t,
           no_source_t>;
using TypeWall = ibm_euler_slip_wall_t<ibmparm_t, ProbClosures>;
using ProbIB = ibm_solver_t<TypeWall, ibmparm_t, ProbClosures>;

inline void update_geometry(Real, Vector<GeomType>&, int) {}

inline void inputs()
{
  amrex::Print()
      << "[TestManifest] case=rz_axis_shared_face_pure_gp_integration"
      << " state=uniform_quiescent solid=off_axis_analytic_circle"
      << " rhs=source_free_skew_jst_o4 rk=ssprk43"
      << " grid=full_background_cells\n";
}

#ifdef CERISSE_TEST_ACTIVE_RHS_HOOK
// Negative qualification fixture: merely defining an active-cell RHS hook
// must make the production face-flux positivity path fail closed.  The body
// is intentionally empty; the incompatibility is the unpaired operator
// channel itself, not the magnitude of a particular case correction.
inline void rhs_nscbc(
    const Geometry&, const MFIter&, const auto&, const auto&,
    const ProbClosures*, const ProbParm*, Real, Real)
{}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata, ProbClosures const& closure,
                   ProbParm const& problem)
{
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx = geomdata.CellSize();
  const Real z = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
  const Real density =
      problem.pressure / (closure.Rspec * problem.temperature);
  Real pressure = problem.pressure;
  Real axial_velocity = Real(0.0);
  if (problem.initial_mode == 1) {
    // Grid-aligned double rarefaction. Density is continuous, so the
    // independent density-only WENO shock fallback does not pre-empt this
    // positivity-limiter integration test.
    axial_velocity = z < problem.interface_z
                           ? -problem.forced_velocity
                           : problem.forced_velocity;
  } else if (problem.initial_mode == 2) {
    // One-cell kinetic-energy spike: a compact unresolved stencil designed
    // to exercise the low Rusanov candidate without changing density.
    if (std::abs(z - problem.interface_z) <= Real(0.51) * dx[1]) {
      axial_velocity = problem.forced_velocity;
    }
  } else if (problem.initial_mode == 3) {
    // One-cell pressure trench, also at constant density. This remains a
    // physical ideal-gas state at t=0 and is useful if the velocity tests are
    // kept admissible by a future high-order reconstruction improvement.
    if (std::abs(z - problem.interface_z) <= Real(0.51) * dx[1]) {
      pressure *= problem.pressure_trench_ratio;
    }
  }
  state(i, j, k, ProbClosures::URHO) = density;
  state(i, j, k, ProbClosures::UMX) = Real(0.0);
  state(i, j, k, ProbClosures::UMY) = density * axial_velocity;
  state(i, j, k, ProbClosures::UMZ) = Real(0.0);
  state(i, j, k, ProbClosures::UET) =
      pressure / (closure.gamma - Real(1.0)) +
      Real(0.5) * density * axial_velocity * axial_velocity;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int, auto& tagfab,
                  const auto& state, const auto& markers,
                  const auto& geomdata, const ProbParm& problem, int level)
{
  amrex::ignore_unused(k, state, markers);
  if (problem.refine_patch == 0 || level != 0) {
    return;
  }
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx = geomdata.CellSize();
  const Real r = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real z = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
  // One rectangular fine patch intersects r=0, the forced axial interface,
  // and the off-axis circle.  It exists only to exercise same-level masks,
  // axis faces, coarse-fine boundaries, and shared GP reconstruction together.
  if (r < Real(0.82) && z > Real(0.10) && z < Real(0.72)) {
    tagfab(i, j, k) = true;
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real x[AMREX_SPACEDIM], Real,
              const Real s_int[ProbClosures::NCONS],
              const Real s_refl[ProbClosures::NCONS],
              Real s_ext[ProbClosures::NCONS], int idir, int sgn, Real time,
              GeometryData const& geomdata, ProbClosures const& closure,
              ProbParm const& problem)
{
  amrex::ignore_unused(x, s_refl, idir, sgn, time, geomdata, closure, problem);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    s_ext[n] = s_int[n];
  }
}

} // namespace PROB

#endif
