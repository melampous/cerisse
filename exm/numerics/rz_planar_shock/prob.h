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

// A stationary normal shock whose exact solution depends only on z.  In R-Z,
// radial uniformity, u_r=0, and omega_theta=0 must be preserved independently
// of the radial metric and of any radial coarse-fine interfaces.
struct ProbParm {
  Real gamma = 1.4;
  Real rho_up = 1.0;
  Real p_up = 1.0;
  Real mach_up = 3.0;
  Real sound_up = 0.0;
  Real uz_up = 0.0;

  Real density_ratio = 0.0;
  Real pressure_ratio = 0.0;
  Real rho_down = 0.0;
  Real p_down = 0.0;
  Real uz_down = 0.0;

  Real shock_z = 0.5;
  Real fine_radius_l1 = 0.16;
  Real fine_radius_l2 = 0.10;
  Real fine_halfwidth_l1 = 0.18;
  Real fine_halfwidth_l2 = 0.12;
  int amr_full_radial = 0;
  Real shock_checkerboard_cells = 0.0;
  int shock_checkerboard_radial_cells = 0;
  int flow_sign = 1;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("mach_up", mach_up);
    pp.query("flow_sign", flow_sign);
    pp.query("amr_full_radial", amr_full_radial);
    pp.query("shock_checkerboard_cells", shock_checkerboard_cells);
    pp.query("shock_checkerboard_radial_cells",
             shock_checkerboard_radial_cells);

    flow_sign = flow_sign < 0 ? -1 : 1;
    sound_up = std::sqrt(gamma * p_up / rho_up);
    uz_up = Real(flow_sign) * mach_up * sound_up;
    density_ratio =
        (gamma + Real(1.0)) * mach_up * mach_up /
        ((gamma - Real(1.0)) * mach_up * mach_up + Real(2.0));
    pressure_ratio =
        Real(1.0) + Real(2.0) * gamma / (gamma + Real(1.0)) *
                        (mach_up * mach_up - Real(1.0));
    rho_down = rho_up * density_ratio;
    p_down = p_up * pressure_ratio;
    uz_down = uz_up / density_ratio;
  }
};

using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                calorifically_perfect_gas_t<indicies_t>>;

// These are the same fourth-order split-form/JST choices used by the current
// Run165 all-fluid control.  They live in this regression case only; selecting
// either option does not add a case branch to the production operator.
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

#ifndef RZ_PLANAR_SHOCK_SCHEME_ID
#define RZ_PLANAR_SHOCK_SCHEME_ID 0
#endif

#if RZ_PLANAR_SHOCK_SCHEME_ID == 0
using ProbConvectiveFlux =
    weno_t<ReconScheme::WenoZ5, ProbClosures>;
static constexpr const char* convective_scheme_name = "llf-wenoz5";
#elif RZ_PLANAR_SHOCK_SCHEME_ID == 1
using ProbConvectiveFlux =
    weno_t<ReconScheme::Teno5, ProbClosures>;
static constexpr const char* convective_scheme_name = "llf-teno5";
#elif RZ_PLANAR_SHOCK_SCHEME_ID == 2
using ProbConvectiveFlux =
    weno_t<ReconScheme::Teno6, ProbClosures>;
static constexpr const char* convective_scheme_name = "llf-teno6";
#elif RZ_PLANAR_SHOCK_SCHEME_ID == 3
using ProbConvectiveFlux = afd_hllc_wenoz5_t<ProbClosures>;
static constexpr const char* convective_scheme_name = "afd-hllc-wenoz5";
#elif RZ_PLANAR_SHOCK_SCHEME_ID == 4
using ProbConvectiveFlux = afd_hllc_teno5_t<ProbClosures>;
static constexpr const char* convective_scheme_name = "afd-hllc-teno5";
#elif RZ_PLANAR_SHOCK_SCHEME_ID == 5
using ProbConvectiveFlux = skew_t<SkewCentralO4Parm, ProbClosures>;
static constexpr const char* convective_scheme_name = "skew4-central";
#elif RZ_PLANAR_SHOCK_SCHEME_ID == 6
using ProbConvectiveFlux = skew_t<SkewJstO4Parm, ProbClosures>;
static constexpr const char* convective_scheme_name = "skew4-jst";
#else
#error "Unsupported RZ_PLANAR_SHOCK_SCHEME_ID"
#endif

using ProbRHS =
    rhs_dt<ProbConvectiveFlux, no_diffusive_t, no_source_t>;

inline void inputs() {
  ProbParm p;
  amrex::Print() << "[TestManifest] case=rz_planar_shock scheme="
                 << convective_scheme_name << '\n';
  amrex::Print() << "R-Z radially uniform stationary normal shock: M1="
                 << p.mach_up << " rho2/rho1=" << p.density_ratio
                 << " p2/p1=" << p.pressure_ratio
                 << " flow_sign=" << p.flow_sign
                 << " scheme=" << convective_scheme_name << '\n';
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
set_normal_shock_state(const bool upstream, Real* state,
                       const ProbClosures& cls, const ProbParm& p) {
  const Real rho = upstream ? p.rho_up : p.rho_down;
  const Real pressure = upstream ? p.p_up : p.p_down;
  const Real uz = upstream ? p.uz_up : p.uz_down;

  state[cls.URHO] = rho;
  state[cls.UMX] = Real(0.0);
  state[cls.UMY] = rho * uz;
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = pressure / (p.gamma - Real(1.0)) +
                   Real(0.5) * rho * uz * uz;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geomdata, ProbClosures const& cls,
              ProbParm const& p) {
  amrex::ignore_unused(i, k);
  const Real z = geomdata.ProbLo(1) +
                 (Real(j) + Real(0.5)) * geomdata.CellSize(1);
  Real local_shock_z = p.shock_z;
  if (p.shock_checkerboard_radial_cells > 0 &&
      i < p.shock_checkerboard_radial_cells) {
    const Real sign = (i & 1) == 0 ? Real(-1.0) : Real(1.0);
    local_shock_z += sign * p.shock_checkerboard_cells *
                     geomdata.CellSize(1);
  }
  Real conserved[ProbClosures::NCONS];
  const bool upstream = p.flow_sign > 0 ? z < local_shock_z
                                        : z > local_shock_z;
  set_normal_shock_state(upstream, conserved, cls, p);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = conserved[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real[AMREX_SPACEDIM], Real,
         const Real s_int[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int idir, const int sgn, Real,
         GeometryData const&, ProbClosures const& cls, ProbParm const& p) {
  for (int n = 0; n < ProbClosures::NCONS; ++n) s_ext[n] = s_int[n];
  if (idir == 1) {
    // The inflow side is upstream; reversing the flow also reverses which
    // physical boundary carries the supersonic state.
    set_normal_shock_state(sgn == p.flow_sign, s_ext, cls, p);
  }
}

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int, TagFab& tagfab,
             const StateFab&, const GeomData& geomdata, const ProbParm& p,
             int level) {
  const Real r = geomdata.ProbLo(0) +
                 (Real(i) + Real(0.5)) * geomdata.CellSize(0);
  const Real z = geomdata.ProbLo(1) +
                 (Real(j) + Real(0.5)) * geomdata.CellSize(1);
  const Real radius = level == 0 ? p.fine_radius_l1 : p.fine_radius_l2;
  const Real halfwidth =
      level == 0 ? p.fine_halfwidth_l1 : p.fine_halfwidth_l2;
  if (level < 2 && (p.amr_full_radial != 0 || r < radius) &&
      amrex::Math::abs(z - p.shock_z) < halfwidth) {
    tagfab(i, j, k) = true;
  }
}

}  // namespace PROB

#endif
