#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

#include <cmath>

using namespace amrex;

namespace PROB {

// Pressure profiles used to isolate the automatic all-fluid R-Z paired
// pressure operator.  The state is otherwise static: rho=const and u=0.
struct ProbParm {
  int profile = 0;
  int cell_average = 0;
  Real gamma = 1.4;
  Real rho0 = 1.0;
  Real p0 = 2.0;
  Real uz0 = 0.0;
  Real ur_slope = 0.0;
  Real amplitude = 0.25;
  Real quartic = 0.10;
  Real layer_z = 0.5;
  Real curvature = 0.20;
  Real width = 0.04;
  // Optional constant-pressure density contact.  It is disabled for the
  // pressure-operator convergence cases above.  orientation=0 places a
  // stationary cylindrical contact at r=contact_location; orientation=1
  // places a stationary plane contact at z=contact_location.
  int density_contact = 0;
  int contact_orientation = 0;
  Real rho_inner = 4.0;
  Real rho_outer = 1.0;
  Real contact_location = 0.0625;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("profile", profile);
    pp.query("cell_average", cell_average);
    pp.query("gamma", gamma);
    pp.query("rho0", rho0);
    pp.query("p0", p0);
    pp.query("uz0", uz0);
    pp.query("ur_slope", ur_slope);
    pp.query("amplitude", amplitude);
    pp.query("quartic", quartic);
    pp.query("layer_z", layer_z);
    pp.query("curvature", curvature);
    pp.query("width", width);
    pp.query("density_contact", density_contact);
    pp.query("contact_orientation", contact_orientation);
    pp.query("rho_inner", rho_inner);
    pp.query("rho_outer", rho_outer);
    pp.query("contact_location", contact_location);
    if (cell_average != 0) {
      amrex::Abort(
          "rz_paired_pressure_axis is a POINT-state regression: "
          "prob.cell_average must be 0");
    }
    if (std::abs(gamma - Real(1.4)) > Real(1.0e-12)) {
      amrex::Abort(
          "rz_paired_pressure_axis uses the gamma=1.4 production closure: "
          "prob.gamma must be 1.4");
    }
  }
};

using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                calorifically_perfect_gas_t<indicies_t>>;

// Match the current Run165 Skew4 definition.  The central variant keeps the
// identical fourth-order split form but disables both JST terms; the JST
// variant uses the Run165 coefficients and the historical first-power sensor.
struct SkewCentralO4Parm {
  static constexpr bool dissipation = false;
  static constexpr int order = 4;
  static constexpr int sensor_power = 1;
  static constexpr Real C2skew = Real(0.0);
  static constexpr Real C4skew = Real(0.0);
};

struct SkewJstO4Parm {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr int sensor_power = 1;
  static constexpr Real C2skew = Real(1.5);
  static constexpr Real C4skew = Real(0.016);
};

#ifndef RZ_AXIS_SCHEME_ID
#define RZ_AXIS_SCHEME_ID 0
#endif

#if RZ_AXIS_SCHEME_ID == 0
using ProbConvectiveFlux = weno_t<ReconScheme::WenoZ5, ProbClosures>;
static constexpr const char* convective_scheme_name = "llf-wenoz5";
#elif RZ_AXIS_SCHEME_ID == 1
using ProbConvectiveFlux = weno_t<ReconScheme::Teno5, ProbClosures>;
static constexpr const char* convective_scheme_name = "llf-teno5";
#elif RZ_AXIS_SCHEME_ID == 2
using ProbConvectiveFlux = weno_t<ReconScheme::Teno6, ProbClosures>;
static constexpr const char* convective_scheme_name = "llf-teno6";
#elif RZ_AXIS_SCHEME_ID == 3
using ProbConvectiveFlux = afd_hllc_wenoz5_t<ProbClosures>;
static constexpr const char* convective_scheme_name = "afd-hllc-wenoz5";
#elif RZ_AXIS_SCHEME_ID == 4
using ProbConvectiveFlux = afd_hllc_teno5_t<ProbClosures>;
static constexpr const char* convective_scheme_name = "afd-hllc-teno5";
#elif RZ_AXIS_SCHEME_ID == 5
using ProbConvectiveFlux = skew_t<SkewCentralO4Parm, ProbClosures>;
static constexpr const char* convective_scheme_name = "skew4-central";
#elif RZ_AXIS_SCHEME_ID == 6
using ProbConvectiveFlux = skew_t<SkewJstO4Parm, ProbClosures>;
static constexpr const char* convective_scheme_name = "skew4-jst";
#else
#error "Unsupported RZ_AXIS_SCHEME_ID"
#endif

using ProbRHS = rhs_dt<ProbConvectiveFlux, no_diffusive_t, no_source_t>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
pressure(const Real r, const Real z, const ProbParm& p) noexcept {
  if (p.profile == 0) {
    return p.p0 + p.amplitude * r * r;
  }
  if (p.profile == 2) {
    const Real r2 = r * r;
    return p.p0 + p.amplitude * r2 + p.quartic * r2 * r2;
  }
  const Real arg =
      (z - p.layer_z - p.curvature * r * r) / p.width;
  return p.p0 + p.amplitude * std::tanh(arg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
density(const Real r, const Real z, const ProbParm& p) noexcept {
  if (p.density_contact == 0) return p.rho0;
  const Real coordinate = p.contact_orientation == 0 ? r : z;
  return coordinate < p.contact_location ? p.rho_inner : p.rho_outer;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
cell_average_pressure(const Real r, const Real z, const Real dr,
                      const Real dz, const ProbParm& p) noexcept {
  const Real rlo = r - Real(0.5) * dr;
  const Real rhi = r + Real(0.5) * dr;
  if (p.profile == 0 || p.profile == 2) {
    const Real rlo2 = rlo * rlo;
    const Real rhi2 = rhi * rhi;
    const Real average_r2 = Real(0.5) * (rlo2 + rhi2);
    const Real average_r4 =
        (rhi2 * rhi2 * rhi2 - rlo2 * rlo2 * rlo2) /
        (Real(3.0) * (rhi2 - rlo2));
    return p.p0 + p.amplitude * average_r2 +
           (p.profile == 2 ? p.quartic * average_r4 : Real(0.0));
  }

  constexpr Real nodes[4] = {
      Real(-0.8611363115940525752), Real(-0.3399810435848562648),
      Real(0.3399810435848562648), Real(0.8611363115940525752)};
  constexpr Real weights[4] = {
      Real(0.3478548451374538574), Real(0.6521451548625461426),
      Real(0.6521451548625461426), Real(0.3478548451374538574)};
  Real numerator = Real(0.0);
  Real denominator = Real(0.0);
  for (int ir = 0; ir < 4; ++ir) {
    const Real rq = r + Real(0.5) * dr * nodes[ir];
    for (int iz = 0; iz < 4; ++iz) {
      const Real zq = z + Real(0.5) * dz * nodes[iz];
      const Real weight = weights[ir] * weights[iz] * rq;
      numerator += weight * pressure(rq, zq, p);
      denominator += weight;
    }
  }
  return numerator / denominator;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_state(const Real pressure_value, const Real density_value,
          const Real radial_coordinate, Real* state, const ProbClosures& cls,
          const ProbParm& p) noexcept {
  const Real ur = p.ur_slope * radial_coordinate;
  state[cls.URHO] = density_value;
  state[cls.UMX] = density_value * ur;
  state[cls.UMY] = density_value * p.uz0;
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = pressure_value / (p.gamma - Real(1.0)) +
                   Real(0.5) * density_value *
                       (ur * ur + p.uz0 * p.uz0);
}

inline void inputs() {
  ProbParm p;
  amrex::Print() << "[TestManifest] case=rz_paired_pressure_axis scheme="
                 << convective_scheme_name << '\n';
  amrex::Print() << "Current R-Z pressure-axis regression: profile="
                 << p.profile << " cell_average=" << p.cell_average
                 << " density_contact=" << p.density_contact
                 << " contact_orientation=" << p.contact_orientation
                 << " uz0=" << p.uz0
                 << " ur_slope=" << p.ur_slope
                 << " scheme=" << convective_scheme_name << '\n';
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const& cls,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real r = geom.ProbLo(0) +
                 (Real(i) + Real(0.5)) * geom.CellSize(0);
  const Real z = geom.ProbLo(1) +
                 (Real(j) + Real(0.5)) * geom.CellSize(1);
  const Real p_state = p.cell_average
      ? cell_average_pressure(
            r, z, geom.CellSize(0), geom.CellSize(1), p)
      : pressure(r, z, p);
  Real conserved[ProbClosures::NCONS];
  set_state(p_state, density(r, z, p), r, conserved, cls, p);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = conserved[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int, const int, Real,
         GeometryData const& geom, ProbClosures const& cls,
         ProbParm const& p) {
  const Real radial = std::abs(x[0]);
  const Real p_state = p.cell_average
      ? cell_average_pressure(
            radial, x[1], geom.CellSize(0), geom.CellSize(1), p)
      : pressure(radial, x[1], p);
  set_state(
      p_state, density(radial, x[1], p), radial, s_ext, cls, p);
}

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, TagFab&, const StateFab&, const GeomData&,
             const ProbParm&, int) {}

}  // namespace PROB

#endif
