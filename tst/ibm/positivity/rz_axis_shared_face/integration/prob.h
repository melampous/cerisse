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

#include <array>
#include <cmath>
#include <limits>
#include <vector>

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
#ifdef CERISSE_TEST_RZ_GP_PARITY_PROBE
  // Deliberately select the generic point-state GP publication path.  This is
  // a test-only A/B against the annular production candidate; it does not
  // enable a cut-control, EB, aggregate, or partial-face update.
  static constexpr bool rz_annular_bic_cell_average = false;
  static constexpr bool rz_bic_support_centre_recovery = false;
  static constexpr bool rz_axis_regular_entropy_jet = false;
#else
  static constexpr bool rz_annular_bic_cell_average = true;
  static constexpr bool rz_bic_support_centre_recovery = true;
  static constexpr bool rz_axis_regular_entropy_jet = true;
#endif
  static constexpr bool smooth_wall_geometry = true;
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_t::euler_slip_analytic_curvature;

#ifdef CERISSE_TEST_RZ_GP_PARITY_PROBE
  static constexpr Real circle_r = Real(0.0);
  static constexpr Real circle_radius = Real(0.15);
#else
  static constexpr Real circle_r = Real(0.55);
  static constexpr Real circle_radius = Real(0.12);
#endif
  static constexpr Real circle_z = Real(0.50);

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
using ProbRHS =
    rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,
           no_source_t>;
using TypeWall = ibm_euler_slip_wall_t<ibmparm_t, ProbClosures>;
#ifdef CERISSE_TEST_RZ_GP_PARITY_PROBE
class ProbIB : public ibm_solver_t<TypeWall, ibmparm_t, ProbClosures>
{
  using Base = ibm_solver_t<TypeWall, ibmparm_t, ProbClosures>;

  struct AxisPair
  {
    int depth = 0;
    int axial_index = 0;
    int mirror_is_gp = 0;
    std::array<Real, 4> negative{};
    std::array<Real, 4> positive{};
  };

  using Snapshot = std::vector<AxisPair>;

  Snapshot snapshot_axis_pairs(const MultiFab& prims_mf, int lev) const
  {
    Gpu::streamSynchronize();
    const auto& geom = this->amr_p->Geom(lev);
    AMREX_ALWAYS_ASSERT(geom.IsRZ());
    AMREX_ALWAYS_ASSERT(geom.Domain().smallEnd(0) == 0);
    AMREX_ALWAYS_ASSERT(prims_mf.nGrow() >= 3);

    Snapshot result;
    const auto& markers = *this->bmf_a[lev];
    for (MFIter mfi(prims_mf, false); mfi.isValid(); ++mfi) {
      const Box& valid = mfi.validbox();
      if (valid.smallEnd(0) != 0) continue;
      const auto prims = prims_mf.const_array(mfi);
      const auto marker = markers.const_array(mfi);
      for (int depth = 1; depth <= 3; ++depth) {
        const int negative_i = -depth;
        const int positive_i = depth - 1;
        AMREX_ALWAYS_ASSERT(prims_mf[mfi].box().contains(
            IntVect(AMREX_D_DECL(negative_i, valid.smallEnd(1), 0))));
        for (int j = valid.smallEnd(1); j <= valid.bigEnd(1); ++j) {
          AxisPair pair;
          pair.depth = depth;
          pair.axial_index = j;
          pair.mirror_is_gp = marker(positive_i, j, 0, 1) != 0 ? 1 : 0;
          constexpr std::array<int, 4> components{
              ProbClosures::QRHO, ProbClosures::QPRES,
              ProbClosures::QV, ProbClosures::QU};
          for (int n = 0; n < 4; ++n) {
            pair.negative[n] = prims(negative_i, j, 0, components[n]);
            pair.positive[n] = prims(positive_i, j, 0, components[n]);
          }
          result.push_back(pair);
        }
      }
    }
    return result;
  }

  static Real parity_error(const AxisPair& pair, int component)
  {
    // rho, p, and axial velocity are even; radial velocity is odd.
    return component == 3
               ? std::abs(pair.negative[component] + pair.positive[component])
               : std::abs(pair.negative[component] - pair.positive[component]);
  }

  static Real parity_tolerance()
  {
    static const Real value = [] {
      Real tolerance = Real(1.0e-12);
      ParmParse pp("prob");
      pp.query("parity_tolerance", tolerance);
      AMREX_ALWAYS_ASSERT(tolerance >= Real(0.0));
      return tolerance;
    }();
    return value;
  }

  static void report_snapshots(const Snapshot& before,
                               const Snapshot& after, int call)
  {
    AMREX_ALWAYS_ASSERT(before.size() == after.size());
    constexpr std::array<const char*, 4> names{"rho", "p", "uz", "ur"};
    Real worst_post_gp_error = Real(0.0);
    Print() << "[RZ-GP-PARITY] call=" << call
            << " method=pure_shared_gp_full_cartesian"
            << " gp_path=generic_nonannular_point_state"
            << " cut=0 eb=0 aggregate=0\n";

    for (int depth = 1; depth <= 3; ++depth) {
      std::array<Real, 4> pre_max{};
      std::array<Real, 4> post_max{};
      std::array<Real, 4> pre_gp_max{};
      std::array<Real, 4> post_gp_max{};
      int gp_count = 0;
      int worst_index = -1;
      Real worst_error = Real(-1.0);
      for (int row = 0; row < static_cast<int>(before.size()); ++row) {
        if (before[row].depth != depth) continue;
        if (before[row].mirror_is_gp != after[row].mirror_is_gp) {
          Abort("GP marker changed during primitive publication");
        }
        if (after[row].mirror_is_gp) ++gp_count;
        for (int n = 0; n < 4; ++n) {
          const Real pre_error = parity_error(before[row], n);
          const Real post_error = parity_error(after[row], n);
          pre_max[n] = amrex::max(pre_max[n], pre_error);
          post_max[n] = amrex::max(post_max[n], post_error);
          if (after[row].mirror_is_gp) {
            pre_gp_max[n] = amrex::max(pre_gp_max[n], pre_error);
            post_gp_max[n] = amrex::max(post_gp_max[n], post_error);
            worst_post_gp_error =
                amrex::max(worst_post_gp_error, post_error);
            if (post_error > worst_error) {
              worst_error = post_error;
              worst_index = row;
            }
          }
        }
      }

      Print() << "[RZ-GP-PARITY-MAX] depth=" << depth
              << " mirror_i=" << depth - 1 << " negative_i=" << -depth
              << " gp_mirror_count=" << gp_count;
      for (int n = 0; n < 4; ++n) {
        Print() << " pre_" << names[n] << '=' << pre_max[n]
                << " post_" << names[n] << '=' << post_max[n]
                << " pre_gp_" << names[n] << '=' << pre_gp_max[n]
                << " post_gp_" << names[n] << '=' << post_gp_max[n];
      }
      Print() << '\n';

      if (worst_index >= 0) {
        const AxisPair& pre = before[worst_index];
        const AxisPair& post = after[worst_index];
        Print() << "[RZ-GP-PARITY-SAMPLE] depth=" << depth
                << " j=" << post.axial_index
                << " mirror_is_gp=" << post.mirror_is_gp;
        for (int n = 0; n < 4; ++n) {
          Print() << " pre_neg_" << names[n] << '=' << pre.negative[n]
                  << " pre_pos_" << names[n] << '=' << pre.positive[n]
                  << " post_neg_" << names[n] << '=' << post.negative[n]
                  << " post_pos_" << names[n] << '=' << post.positive[n];
        }
        Print() << '\n';
      }
    }
    if (worst_post_gp_error > parity_tolerance()) {
      Abort("generic point-GP physical-axis parity gate failed");
    }
  }

  static bool apply_post_gp_axis_parity()
  {
    static const bool enabled = [] {
      int value = 0;
      ParmParse pp("prob");
      pp.query("apply_post_gp_axis_parity", value);
      AMREX_ALWAYS_ASSERT(value == 0 || value == 1);
      return value != 0;
    }();
    return enabled;
  }

  static void report_refresh(const Snapshot& stale,
                             const Snapshot& refreshed)
  {
    AMREX_ALWAYS_ASSERT(stale.size() == refreshed.size());
    constexpr std::array<const char*, 4> names{"rho", "p", "uz", "ur"};
    for (int depth = 1; depth <= 3; ++depth) {
      std::array<Real, 4> stale_gp_max{};
      std::array<Real, 4> refreshed_gp_max{};
      for (int row = 0; row < static_cast<int>(stale.size()); ++row) {
        if (stale[row].depth != depth || !stale[row].mirror_is_gp) continue;
        for (int n = 0; n < 4; ++n) {
          stale_gp_max[n] = amrex::max(
              stale_gp_max[n], parity_error(stale[row], n));
          refreshed_gp_max[n] = amrex::max(
              refreshed_gp_max[n], parity_error(refreshed[row], n));
        }
      }
      Print() << "[RZ-GP-PARITY-REFRESH] depth=" << depth;
      for (int n = 0; n < 4; ++n) {
        Print() << " stale_gp_" << names[n] << '=' << stale_gp_max[n]
                << " refreshed_gp_" << names[n] << '='
                << refreshed_gp_max[n];
      }
      Print() << '\n';
    }
  }

public:
  void computeAllGPs(MultiFab& prims_mf, const ProbClosures* cls,
                     int lev, int rz_point_target = -1)
  {
    static int call_count = 0;
    const int call = call_count++;
    if (call != 0) {
      Base::computeAllGPs(prims_mf, cls, lev, rz_point_target);
      if (apply_post_gp_axis_parity()) {
        this->fillRZAxisPrimitiveParity(prims_mf, lev);
      }
      return;
    }

    // Match the first operation of Base::computeAllGPs, then retain an exact
    // host snapshot immediately before the GP kernel publishes its states.
    prims_mf.FillBoundary(this->amr_p->Geom(lev).periodicity());
    const Snapshot before = snapshot_axis_pairs(prims_mf, lev);
    Base::computeAllGPs(prims_mf, cls, lev, rz_point_target);
    const Snapshot after = snapshot_axis_pairs(prims_mf, lev);
    report_snapshots(before, after, call);
    if (apply_post_gp_axis_parity()) {
      this->fillRZAxisPrimitiveParity(prims_mf, lev);
      const Snapshot refreshed = snapshot_axis_pairs(prims_mf, lev);
      report_refresh(after, refreshed);
    }
  }
};
#else
using ProbIB = ibm_solver_t<TypeWall, ibmparm_t, ProbClosures>;
#endif

inline void update_geometry(Real, Vector<GeomType>&, int) {}

inline void inputs()
{
#ifdef CERISSE_TEST_RZ_GP_PARITY_PROBE
  amrex::Print()
      << "[TestManifest] case=rz_generic_point_gp_axis_parity_microgate"
      << " method=pure_shared_gp_full_cartesian"
      << " gp_path=generic_nonannular_point_state"
      << " geometry=axis_intersecting_semicircle"
      << " cut=0 eb=0 aggregate=0\n";
#else
  amrex::Print()
      << "[TestManifest] case=rz_axis_shared_face_pure_gp_integration"
      << " state=uniform_quiescent solid=off_axis_analytic_circle"
      << " rhs=source_free_llf_wenoz5 rk=ssprk43"
      << " grid=full_background_cells\n";
#endif
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
  const Real r = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real z = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
  Real pressure = problem.pressure;
  Real temperature = problem.temperature;
  Real radial_velocity = Real(0.0);
  Real axial_velocity = Real(0.0);
#ifdef CERISSE_TEST_RZ_GP_PARITY_PROBE
  // Axis-regular, deliberately non-wall-compatible point field.  FillPatch
  // therefore gives an exact even/even/even/odd axis mirror before GP
  // publication, while the slip-wall GP reconstruction measurably changes
  // solid-side mirror cells adjacent to the axis-intersecting body.
  pressure *= Real(1.0) + Real(0.2) * r * r;
  temperature *= Real(1.0) + Real(0.1) * r * r;
  radial_velocity = Real(80.0) * r;
  axial_velocity = Real(40.0) + Real(5.0) * r * r;
#else
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
#endif
  const Real density = pressure / (closure.Rspec * temperature);
  state(i, j, k, ProbClosures::URHO) = density;
  state(i, j, k, ProbClosures::UMX) = density * radial_velocity;
  state(i, j, k, ProbClosures::UMY) = density * axial_velocity;
  state(i, j, k, ProbClosures::UMZ) = Real(0.0);
  state(i, j, k, ProbClosures::UET) =
      pressure / (closure.gamma - Real(1.0)) +
      Real(0.5) * density *
          (radial_velocity * radial_velocity +
           axial_velocity * axial_velocity);
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
