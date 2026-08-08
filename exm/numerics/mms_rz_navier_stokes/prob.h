#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

#ifndef MMS_VISCOUS
#define MMS_VISCOUS 1
#endif

#ifndef RZ_MMS_EULER_SCHEME_ID
#define RZ_MMS_EULER_SCHEME_ID 0
#endif

#ifndef RZ_MMS_VISCOUS_ORDER
#define RZ_MMS_VISCOUS_ORDER 2
#endif

#ifndef RZ_MMS_VARIABLE_TRANSPORT
#define RZ_MMS_VARIABLE_TRANSPORT 0
#endif

namespace PROB {

constexpr Real mms_gamma = Real(1.4);
constexpr Real mms_rspec = Real(1.0);
constexpr Real mms_mu = Real(0.01);
constexpr Real mms_conductivity = Real(0.02);
// Nondimensional Sutherland calibration matching the Cartesian MMS: the laws
// return mms_mu and mms_conductivity at T = mms_temperature_ref, with
// moderate constants so the coefficient variation is large enough to expose
// variable-coefficient truncation errors.
constexpr Real mms_temperature_ref = Real(1.0);
constexpr Real mms_viscosity_sutherland = Real(0.5);
constexpr Real mms_conductivity_sutherland = Real(0.8);

static_assert(RZ_MMS_VARIABLE_TRANSPORT == 0 || RZ_MMS_VARIABLE_TRANSPORT == 1,
              "RZ_MMS_VARIABLE_TRANSPORT must be 0 or 1");

struct ProbParm {
  // Optional fixed Level-1 radial band for smooth AMR interface tests.
  // The explicit enable flag keeps every established uniform run unchanged.
  int static_radial_refinement = 0;
  Real refine_radial_lo = 0.25;
  Real refine_radial_hi = 0.75;
  Real rho0 = 1.0;
  Real temperature0 = 1.0;
  Real uz0 = 0.30;
  Real rho_amp = 0.050;
  Real temperature_amp = 0.060;
  Real ur_amp = 0.080;
  Real uz_amp = 0.050;
  Real rho_r2 = 0.20;
  Real temperature_r2 = 0.25;
  Real ur_r2 = 0.15;
  Real uz_r2 = 0.10;
  Real temperature_phase = 0.21;
  Real ur_phase = 0.13;
  Real uz_phase = 0.35;
  Real omega = 1.0;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("static_radial_refinement", static_radial_refinement);
    pp.query("refine_radial_lo", refine_radial_lo);
    pp.query("refine_radial_hi", refine_radial_hi);
    pp.query("rho0", rho0);
    pp.query("temperature0", temperature0);
    pp.query("uz0", uz0);
    pp.query("rho_amp", rho_amp);
    pp.query("temperature_amp", temperature_amp);
    pp.query("ur_amp", ur_amp);
    pp.query("uz_amp", uz_amp);
    pp.query("rho_r2", rho_r2);
    pp.query("temperature_r2", temperature_r2);
    pp.query("ur_r2", ur_r2);
    pp.query("uz_r2", uz_r2);
    pp.query("temperature_phase", temperature_phase);
    pp.query("ur_phase", ur_phase);
    pp.query("uz_phase", uz_phase);
    pp.query("omega", omega);
  }
};

struct TransportParm {
  static_assert(RZ_MMS_VISCOUS_ORDER == 2 || RZ_MMS_VISCOUS_ORDER == 4 ||
                    RZ_MMS_VISCOUS_ORDER == 6,
                "RZ_MMS_VISCOUS_ORDER must be 2, 4, or 6");
  static constexpr int order = RZ_MMS_VISCOUS_ORDER;
  static constexpr bool use_LES = false;
  static constexpr Real viscosity = mms_mu;
  static constexpr Real conductivity = mms_conductivity;
};

struct GasParm {
  static constexpr Real gamma = mms_gamma;
  // This gives the nondimensional gas constant R=Ru/M=1.
  static constexpr Real molecular_weight =
      universal_constants::gas_constant;
};

// Nondimensional Sutherland transport for the solver side, mirroring the
// Cartesian MMS implementation.
class MmsSutherlandTransport {
 public:
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
  visc(const Real& temperature) const noexcept {
    const Real ratio = temperature / mms_temperature_ref;
    return mms_mu * ratio * ::sqrt(ratio) *
           (mms_temperature_ref + mms_viscosity_sutherland) /
           (temperature + mms_viscosity_sutherland);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
  cond(const Real& temperature) const noexcept {
    const Real ratio = temperature / mms_temperature_ref;
    return mms_conductivity * ratio * ::sqrt(ratio) *
           (mms_temperature_ref + mms_conductivity_sutherland) /
           (temperature + mms_conductivity_sutherland);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
  xi(const Real&) const noexcept {
    return Real(0.0);
  }
};

#if RZ_MMS_VARIABLE_TRANSPORT
using MmsTransport = MmsSutherlandTransport;
#else
using MmsTransport = transport_const_t<TransportParm>;
#endif

using ProbClosures =
    closures_dt<indicies_t, MmsTransport,
                perfect_gas_t<GasParm, indicies_t>>;

template <typename cls_t>
class user_source_t;

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

struct SkewCentralO6Parm {
  static constexpr bool dissipation = false;
  static constexpr int order = 6;
  static constexpr Real C2skew = Real(0.0);
  static constexpr Real C4skew = Real(0.0);
};

#if (RZ_MMS_EULER_SCHEME_ID == 0)
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "llf-wenoz5";
#elif (RZ_MMS_EULER_SCHEME_ID == 1)
using ProbEuler = skew_t<SkewCentralO4Parm, ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "skew-central-o4";
#elif (RZ_MMS_EULER_SCHEME_ID == 2)
using ProbEuler = skew_t<SkewJstO4Parm, ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "skew-jst-o4";
#elif (RZ_MMS_EULER_SCHEME_ID == 3)
using ProbEuler = weno_t<ReconScheme::Teno5, ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "llf-teno5";
#elif (RZ_MMS_EULER_SCHEME_ID == 4)
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "afd-hllc-wenoz5";
#elif (RZ_MMS_EULER_SCHEME_ID == 5)
using ProbEuler = afd_hllc_teno5_t<ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "afd-hllc-teno5";
#elif (RZ_MMS_EULER_SCHEME_ID == 6)
using ProbEuler = skew_t<SkewCentralO6Parm, ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "skew-central-o6";
#elif (RZ_MMS_EULER_SCHEME_ID == 7)
using ProbEuler = weno_t<ReconScheme::Teno6, ProbClosures>;
static constexpr const char* mms_convective_scheme_name = "llf-teno6";
#else
#error "Unknown RZ_MMS_EULER_SCHEME_ID"
#endif

#if MMS_VISCOUS
using ProbRHS =
    rhs_dt<ProbEuler, viscous_t<TransportParm, ProbClosures>,
           user_source_t<ProbClosures>>;
#else
using ProbRHS =
    rhs_dt<ProbEuler, no_diffusive_t, user_source_t<ProbClosures>>;
#endif

// Value, first derivatives in r,z,t, and the spatial Hessian.
struct Jet {
  Real v = 0.0;
  Real dr = 0.0;
  Real dz = 0.0;
  Real dt = 0.0;
  Real drr = 0.0;
  Real drz = 0.0;
  Real dzz = 0.0;

  AMREX_GPU_HOST_DEVICE Jet() = default;
  AMREX_GPU_HOST_DEVICE explicit Jet(const Real value) : v(value) {}
  AMREX_GPU_HOST_DEVICE
  Jet(const Real value, const Real radial, const Real axial,
      const Real temporal, const Real radial2 = 0.0,
      const Real radial_axial = 0.0, const Real axial2 = 0.0)
      : v(value),
        dr(radial),
        dz(axial),
        dt(temporal),
        drr(radial2),
        drz(radial_axial),
        dzz(axial2) {}
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator+(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v + b.v, a.dr + b.dr, a.dz + b.dz, a.dt + b.dt,
             a.drr + b.drr, a.drz + b.drz, a.dzz + b.dzz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v - b.v, a.dr - b.dr, a.dz - b.dz, a.dt - b.dt,
             a.drr - b.drr, a.drz - b.drz, a.dzz - b.dzz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a) noexcept {
  return Jet(-a.v, -a.dr, -a.dz, -a.dt, -a.drr, -a.drz, -a.dzz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator*(const Jet& a, const Jet& b) noexcept {
  return Jet(
      a.v * b.v,
      a.dr * b.v + a.v * b.dr,
      a.dz * b.v + a.v * b.dz,
      a.dt * b.v + a.v * b.dt,
      a.drr * b.v + Real(2.0) * a.dr * b.dr + a.v * b.drr,
      a.drz * b.v + a.dr * b.dz + a.dz * b.dr + a.v * b.drz,
      a.dzz * b.v + Real(2.0) * a.dz * b.dz + a.v * b.dzz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_inverse(const Jet& a) noexcept {
  const Real inv = Real(1.0) / a.v;
  const Real inv2 = inv * inv;
  const Real inv3 = inv2 * inv;
  return Jet(inv, -a.dr * inv2, -a.dz * inv2, -a.dt * inv2,
             Real(2.0) * a.dr * a.dr * inv3 - a.drr * inv2,
             Real(2.0) * a.dr * a.dz * inv3 - a.drz * inv2,
             Real(2.0) * a.dz * a.dz * inv3 - a.dzz * inv2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator/(const Jet& a, const Jet& b) noexcept {
  return a * jet_inverse(b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator+(const Jet& a, const Real b) noexcept {
  return a + Jet(b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator+(const Real a, const Jet& b) noexcept {
  return Jet(a) + b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a, const Real b) noexcept {
  return a - Jet(b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Real a, const Jet& b) noexcept {
  return Jet(a) - b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator*(const Jet& a, const Real b) noexcept {
  return Jet(a.v * b, a.dr * b, a.dz * b, a.dt * b,
             a.drr * b, a.drz * b, a.dzz * b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator*(const Real a, const Jet& b) noexcept {
  return b * a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator/(const Jet& a, const Real b) noexcept {
  return a * (Real(1.0) / b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_sin(const Jet& a) noexcept {
  const Real sine = ::sin(a.v);
  const Real cosine = ::cos(a.v);
  return Jet(
      sine, cosine * a.dr, cosine * a.dz, cosine * a.dt,
      cosine * a.drr - sine * a.dr * a.dr,
      cosine * a.drz - sine * a.dr * a.dz,
      cosine * a.dzz - sine * a.dz * a.dz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_cos(const Jet& a) noexcept {
  const Real sine = ::sin(a.v);
  const Real cosine = ::cos(a.v);
  return Jet(
      cosine, -sine * a.dr, -sine * a.dz, -sine * a.dt,
      -sine * a.drr - cosine * a.dr * a.dr,
      -sine * a.drz - cosine * a.dr * a.dz,
      -sine * a.dzz - cosine * a.dz * a.dz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_sqrt(const Jet& a) noexcept {
  const Real root = ::sqrt(a.v);
  const Real first = Real(0.5) / root;
  const Real second = -Real(0.25) / (a.v * root);
  return Jet(
      root, first * a.dr, first * a.dz, first * a.dt,
      first * a.drr + second * a.dr * a.dr,
      first * a.drz + second * a.dr * a.dz,
      first * a.dzz + second * a.dz * a.dz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
mms_sutherland_law(const Jet& temperature, const Real reference_value,
                   const Real sutherland_constant) noexcept {
  const Jet ratio = temperature / mms_temperature_ref;
  return reference_value * ratio * jet_sqrt(ratio) *
         (mms_temperature_ref + sutherland_constant) /
         (temperature + sutherland_constant);
}

// Only the value and first spatial derivatives of these derivative fields
// are used when the viscous flux divergence is formed.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
radial_derivative(const Jet& a) noexcept {
  return Jet(a.dr, a.drr, a.drz, Real(0.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
axial_derivative(const Jet& a) noexcept {
  return Jet(a.dz, a.drz, a.dzz, Real(0.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
manufactured_point(const Real radial, const Real axial, const Real time,
                   const ProbParm& p, Jet* conserved, Jet* convective_r,
                   Jet* convective_z, Jet* viscous_r, Jet* viscous_z,
                   Jet& pressure, Jet& tau_theta_theta) noexcept {
  constexpr Real pi = 3.141592653589793238462643383279502884;
  const Jet r(radial, Real(1.0), Real(0.0), Real(0.0));
  const Jet z(axial, Real(0.0), Real(1.0), Real(0.0));
  const Jet t(time, Real(0.0), Real(0.0), Real(1.0));
  const Jet r2 = r * r;
  const Jet phase = Real(2.0) * pi * z - p.omega * t;

  const Jet rho =
      p.rho0 +
      p.rho_amp * (Real(1.0) + p.rho_r2 * r2) * jet_sin(phase);
  const Jet temperature =
      p.temperature0 +
      p.temperature_amp *
          (Real(1.0) + p.temperature_r2 * r2) *
          jet_cos(phase - p.temperature_phase);
  const Jet ur =
      p.ur_amp * r * (Real(1.0) + p.ur_r2 * r2) *
      jet_cos(phase + p.ur_phase);
  const Jet uz =
      p.uz0 +
      p.uz_amp * (Real(1.0) + p.uz_r2 * r2) *
          jet_sin(phase + p.uz_phase);
  pressure = rho * mms_rspec * temperature;

  const Jet energy =
      pressure / (mms_gamma - Real(1.0)) +
      Real(0.5) * rho * (ur * ur + uz * uz);

  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    conserved[n] = Jet();
    convective_r[n] = Jet();
    convective_z[n] = Jet();
    viscous_r[n] = Jet();
    viscous_z[n] = Jet();
  }

  conserved[ProbClosures::URHO] = rho;
  conserved[ProbClosures::UMX] = rho * ur;
  conserved[ProbClosures::UMY] = rho * uz;
  conserved[ProbClosures::UMZ] = Jet();
  conserved[ProbClosures::UET] = energy;

  convective_r[ProbClosures::URHO] = rho * ur;
  convective_r[ProbClosures::UMX] = rho * ur * ur + pressure;
  convective_r[ProbClosures::UMY] = rho * ur * uz;
  convective_r[ProbClosures::UET] = ur * (energy + pressure);

  convective_z[ProbClosures::URHO] = rho * uz;
  convective_z[ProbClosures::UMX] = rho * ur * uz;
  convective_z[ProbClosures::UMY] = rho * uz * uz + pressure;
  convective_z[ProbClosures::UET] = uz * (energy + pressure);

  const Jet dur_dr = radial_derivative(ur);
  const Jet dur_dz = axial_derivative(ur);
  const Jet duz_dr = radial_derivative(uz);
  const Jet duz_dz = axial_derivative(uz);
  const Jet ur_over_r = ur / r;
  const Jet theta = dur_dr + ur_over_r + duz_dz;

#if RZ_MMS_VARIABLE_TRANSPORT
  const Jet viscosity =
      mms_sutherland_law(temperature, mms_mu, mms_viscosity_sutherland);
  const Jet conductivity_jet = mms_sutherland_law(
      temperature, mms_conductivity, mms_conductivity_sutherland);
#else
  const Jet viscosity(mms_mu);
  const Jet conductivity_jet(mms_conductivity);
#endif

  const Jet tau_rr =
      viscosity * (Real(2.0) * dur_dr - Real(2.0 / 3.0) * theta);
  const Jet tau_rz = viscosity * (dur_dz + duz_dr);
  const Jet tau_zz =
      viscosity * (Real(2.0) * duz_dz - Real(2.0 / 3.0) * theta);
  tau_theta_theta =
      viscosity * (Real(2.0) * ur_over_r - Real(2.0 / 3.0) * theta);

  const Jet dtemperature_dr = radial_derivative(temperature);
  const Jet dtemperature_dz = axial_derivative(temperature);

  viscous_r[ProbClosures::UMX] = tau_rr;
  viscous_r[ProbClosures::UMY] = tau_rz;
  viscous_r[ProbClosures::UET] =
      ur * tau_rr + uz * tau_rz +
      conductivity_jet * dtemperature_dr;

  viscous_z[ProbClosures::UMX] = tau_rz;
  viscous_z[ProbClosures::UMY] = tau_zz;
  viscous_z[ProbClosures::UET] =
      ur * tau_rz + uz * tau_zz +
      conductivity_jet * dtemperature_dz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
point_state(const Real r, const Real z, const Real time,
            const ProbParm& p, Real* state) noexcept {
  Jet conserved[ProbClosures::NCONS];
  Jet convective_r[ProbClosures::NCONS];
  Jet convective_z[ProbClosures::NCONS];
  Jet viscous_r[ProbClosures::NCONS];
  Jet viscous_z[ProbClosures::NCONS];
  Jet pressure;
  Jet tau_theta_theta;
  manufactured_point(r, z, time, p, conserved, convective_r, convective_z,
                     viscous_r, viscous_z, pressure, tau_theta_theta);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state[n] = conserved[n].v;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
point_source(const Real r, const Real z, const Real time,
             const ProbParm& p, Real* source) noexcept {
  Jet conserved[ProbClosures::NCONS];
  Jet convective_r[ProbClosures::NCONS];
  Jet convective_z[ProbClosures::NCONS];
  Jet viscous_r[ProbClosures::NCONS];
  Jet viscous_z[ProbClosures::NCONS];
  Jet pressure;
  Jet tau_theta_theta;
  manufactured_point(r, z, time, p, conserved, convective_r, convective_z,
                     viscous_r, viscous_z, pressure, tau_theta_theta);

  const Real inv_r = Real(1.0) / r;
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
#if MMS_VISCOUS
    source[n] =
        conserved[n].dt + convective_r[n].dr + convective_z[n].dz -
        viscous_r[n].dr - viscous_z[n].dz;

    if (n == ProbClosures::UMX) {
      // The solver adds +p/r and -tau_theta_theta/r separately. Combine
      // the corresponding analytic terms before division. This avoids
      // cancellation of individually singular terms near r=0.
      source[n] +=
          (convective_r[n].v - pressure.v -
           viscous_r[n].v + tau_theta_theta.v) *
          inv_r;
    } else {
      source[n] +=
          (convective_r[n].v - viscous_r[n].v) * inv_r;
    }
#else
    source[n] =
        conserved[n].dt + convective_r[n].dr + convective_z[n].dz;
    if (n == ProbClosures::UMX) {
      source[n] += (convective_r[n].v - pressure.v) * inv_r;
    } else {
      source[n] += convective_r[n].v * inv_r;
    }
#endif
  }
}

inline void inputs() {
  const ProbParm p;
  if (p.static_radial_refinement != 0 &&
      p.static_radial_refinement != 1) {
    amrex::Abort("prob.static_radial_refinement must be zero or one");
  }
  if (p.static_radial_refinement != 0 &&
      !(p.refine_radial_lo < p.refine_radial_hi)) {
    amrex::Abort(
        "prob.refine_radial_lo must be below prob.refine_radial_hi");
  }
  const Real radial_extent = Real(1.0);
  const Real rho_bound =
      amrex::Math::abs(p.rho_amp) *
      (Real(1.0) + amrex::Math::abs(p.rho_r2) * radial_extent);
  const Real temperature_bound =
      amrex::Math::abs(p.temperature_amp) *
      (Real(1.0) +
       amrex::Math::abs(p.temperature_r2) * radial_extent);
  if (p.rho0 <= rho_bound || p.temperature0 <= temperature_bound) {
    amrex::Abort(
        "R-Z Navier-Stokes MMS parameters do not guarantee positive rho and T");
  }
  amrex::Print()
      << "[TestManifest] case=mms_rz_navier_stokes physics="
#if MMS_VISCOUS
      << "navier-stokes"
#else
      << "euler"
#endif
      << " scheme=" << mms_convective_scheme_name << '\n';
  amrex::Print()
      << "R-Z time-dependent point-sample MMS, physics="
#if MMS_VISCOUS
      << "Navier-Stokes\n"
      << "  mu=" << mms_mu
      << "  conductivity=" << mms_conductivity
      << "  R=" << mms_rspec << '\n';
#else
      << "Euler\n"
      << "  R=" << mms_rspec << '\n';
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const&,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real r =
      geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
  const Real z =
      geom.ProbLo(1) + (Real(j) + Real(0.5)) * geom.CellSize(1);
  Real exact[ProbClosures::NCONS];
  point_state(r, z, Real(0.0), p, exact);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = exact[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int, const int,
         const Real time, GeometryData const&, ProbClosures const&,
         ProbParm const& p) {
  point_state(x[0], x[1], time, p, s_ext);
}

template <typename cls_t>
class user_source_t {
 public:
  user_source_t() = default;

  void src(const Geometry& geometry, const MFIter& mfi,
           const Array4<const Real>& prims, const Array4<Real>& rhs,
           const cls_t* closures, const Real dt, const Real time) {
    amrex::ignore_unused(prims, closures, dt);
    const Box box = mfi.tilebox();
    const GeometryData geom = geometry.data();
    const ProbParm p = parameters_;
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real r =
          geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
      const Real z =
          geom.ProbLo(1) + (Real(j) + Real(0.5)) * geom.CellSize(1);
      Real source[ProbClosures::NCONS];
      point_source(r, z, time, p, source);
      for (int n = 0; n < ProbClosures::NCONS; ++n) {
        rhs(i, j, k, n) += source[n];
      }
    });
  }

 private:
  ProbParm parameters_{};
};

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int, TagFab& tagfab, const StateFab&,
             const GeomData& geom, const ProbParm& p, int level) {
  if (p.static_radial_refinement == 0 || level != 0) {
    return;
  }

  const Real radial_coordinate =
      geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
  if (radial_coordinate >= p.refine_radial_lo &&
      radial_coordinate < p.refine_radial_hi) {
    tagfab(i, j, k) = true;
  }
}

}  // namespace PROB

#endif
