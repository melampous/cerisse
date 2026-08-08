#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>
#include <Weno_old.h>

using namespace amrex;

#ifndef MMS_VISCOUS
#define MMS_VISCOUS 0
#endif

#ifndef MMS_VISCOUS_ORDER
#define MMS_VISCOUS_ORDER 2
#endif

#ifndef MMS_EULER_SCHEME_ID
#define MMS_EULER_SCHEME_ID 0
#endif

#ifndef MMS_VARIABLE_TRANSPORT
#define MMS_VARIABLE_TRANSPORT 0
#endif

#ifndef MMS_SKEW_ORDER
#define MMS_SKEW_ORDER 4
#endif

#ifndef MMS_SKEW_DISSIPATION
#define MMS_SKEW_DISSIPATION 0
#endif

#ifndef MMS_SKEW_DAMP_ORDER
#define MMS_SKEW_DAMP_ORDER 0
#endif

#ifndef MMS_SKEW_SENSOR_POWER
#define MMS_SKEW_SENSOR_POWER 1
#endif

namespace PROB {

constexpr Real mms_gamma = Real(1.4);
constexpr Real mms_rspec = Real(1.0);
constexpr Real mms_mu = Real(0.01);
constexpr Real mms_conductivity = Real(0.02);
constexpr Real mms_temperature_ref = Real(1.0);
constexpr Real mms_viscosity_sutherland = Real(0.5);
constexpr Real mms_conductivity_sutherland = Real(0.8);

static_assert(MMS_VARIABLE_TRANSPORT == 0 || MMS_VARIABLE_TRANSPORT == 1,
              "MMS_VARIABLE_TRANSPORT must be 0 or 1");

struct ProbParm {
  Real rho0 = 1.0;
  Real p0 = 1.0;
  Real u0 = 0.45;
  Real v0 = 0.30;
  Real rho_x_amp = 0.08;
  Real rho_y_amp = 0.05;
  Real u_x_amp = 0.07;
  Real u_y_amp = 0.06;
  Real v_x_amp = 0.05;
  Real v_y_amp = 0.04;
  Real p_x_amp = 0.07;
  Real p_y_amp = 0.06;
  Real omega = 1.0;
  // Static refinement box for AMR convergence studies.  The default
  // (lo > hi) tags nothing, so single-level runs are unaffected.
  Real refine_xlo = 1.0;
  Real refine_xhi = -1.0;
  Real refine_ylo = 1.0;
  Real refine_yhi = -1.0;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("rho0", rho0);
    pp.query("p0", p0);
    pp.query("u0", u0);
    pp.query("v0", v0);
    pp.query("rho_x_amp", rho_x_amp);
    pp.query("rho_y_amp", rho_y_amp);
    pp.query("u_x_amp", u_x_amp);
    pp.query("u_y_amp", u_y_amp);
    pp.query("v_x_amp", v_x_amp);
    pp.query("v_y_amp", v_y_amp);
    pp.query("p_x_amp", p_x_amp);
    pp.query("p_y_amp", p_y_amp);
    pp.query("omega", omega);
    pp.query("refine_xlo", refine_xlo);
    pp.query("refine_xhi", refine_xhi);
    pp.query("refine_ylo", refine_ylo);
    pp.query("refine_yhi", refine_yhi);
  }
};

struct TransportParm {
  static_assert(MMS_VISCOUS_ORDER == 2 || MMS_VISCOUS_ORDER == 4 ||
                    MMS_VISCOUS_ORDER == 6,
                "MMS_VISCOUS_ORDER must be 2, 4, or 6");
  static constexpr int order = MMS_VISCOUS_ORDER;
  static constexpr bool use_LES = false;
  static constexpr Real viscosity = mms_mu;
  static constexpr Real conductivity = mms_conductivity;
};

struct GasParm {
  static constexpr Real gamma = mms_gamma;
  // The choice makes R=Ru/M=1 in the nondimensional verification problem.
  static constexpr Real molecular_weight = gas_constant;
};

// Nondimensional Sutherland laws calibrated so that mu=mms_mu and
// k=mms_conductivity at T=mms_temperature_ref.  The moderate Sutherland
// constants keep both coefficients positive and make their variation large
// enough to expose variable-coefficient truncation errors in this MMS.
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

#if MMS_VARIABLE_TRANSPORT
using MmsTransport = MmsSutherlandTransport;
#else
using MmsTransport = transport_const_t<TransportParm>;
#endif

using ProbClosures =
    closures_dt<indicies_t, MmsTransport,
                perfect_gas_t<GasParm, indicies_t>>;

template <typename cls_t>
class user_source_t;

#if (MMS_EULER_SCHEME_ID == 0)
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 1)
using ProbEuler = weno_t<ReconScheme::Teno5, ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 2)
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 3)
using ProbEuler = afd_hllc_teno5_t<ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 4)
using ProbEuler = weno_t<ReconScheme::Teno6, ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 5)
using ProbEuler = weno_old_t<OldReconScheme::WenoZ5, ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 6)
using ProbEuler = weno_old_t<OldReconScheme::Teno5, ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 7)
using ProbEuler = weno_old_t<OldReconScheme::Teno6, ProbClosures>;
#elif (MMS_EULER_SCHEME_ID == 8)
struct SkewParm {
  static constexpr bool dissipation = (MMS_SKEW_DISSIPATION != 0);
  static constexpr int order = MMS_SKEW_ORDER;
  // Documented defaults for the JST shock and background damping constants.
  static constexpr Real C2skew = 0.1;
  static constexpr Real C4skew = 0.016;
  // 0 selects the historical damping stencil (difference order = order-1).
  static constexpr int damp_difference_order = MMS_SKEW_DAMP_ORDER;
  static constexpr int sensor_power = MMS_SKEW_SENSOR_POWER;
};
using ProbEuler = skew_t<SkewParm, ProbClosures>;
#else
#error "Unknown MMS_EULER_SCHEME_ID"
#endif

#if MMS_VISCOUS
using ProbRHS = rhs_dt<ProbEuler, viscous_t<TransportParm, ProbClosures>,
                       user_source_t<ProbClosures>>;
#else
using ProbRHS =
    rhs_dt<ProbEuler, no_diffusive_t, user_source_t<ProbClosures>>;
#endif

// Value, first derivatives in x,y,t, and the spatial Hessian.
struct Jet {
  Real v = 0.0;
  Real dx = 0.0;
  Real dy = 0.0;
  Real dt = 0.0;
  Real dxx = 0.0;
  Real dxy = 0.0;
  Real dyy = 0.0;

  AMREX_GPU_HOST_DEVICE Jet() = default;
  AMREX_GPU_HOST_DEVICE explicit Jet(const Real value) : v(value) {}
  AMREX_GPU_HOST_DEVICE
  Jet(const Real value, const Real x1, const Real y1, const Real t1,
      const Real x2 = 0.0, const Real xy = 0.0, const Real y2 = 0.0)
      : v(value), dx(x1), dy(y1), dt(t1), dxx(x2), dxy(xy), dyy(y2) {}
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator+(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v + b.v, a.dx + b.dx, a.dy + b.dy, a.dt + b.dt,
             a.dxx + b.dxx, a.dxy + b.dxy, a.dyy + b.dyy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v - b.v, a.dx - b.dx, a.dy - b.dy, a.dt - b.dt,
             a.dxx - b.dxx, a.dxy - b.dxy, a.dyy - b.dyy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a) noexcept {
  return Jet(-a.v, -a.dx, -a.dy, -a.dt, -a.dxx, -a.dxy, -a.dyy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator*(const Jet& a, const Jet& b) noexcept {
  return Jet(
      a.v * b.v,
      a.dx * b.v + a.v * b.dx,
      a.dy * b.v + a.v * b.dy,
      a.dt * b.v + a.v * b.dt,
      a.dxx * b.v + Real(2.0) * a.dx * b.dx + a.v * b.dxx,
      a.dxy * b.v + a.dx * b.dy + a.dy * b.dx + a.v * b.dxy,
      a.dyy * b.v + Real(2.0) * a.dy * b.dy + a.v * b.dyy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_inverse(const Jet& a) noexcept {
  const Real inv = Real(1.0) / a.v;
  const Real inv2 = inv * inv;
  const Real inv3 = inv2 * inv;
  return Jet(inv, -a.dx * inv2, -a.dy * inv2, -a.dt * inv2,
             Real(2.0) * a.dx * a.dx * inv3 - a.dxx * inv2,
             Real(2.0) * a.dx * a.dy * inv3 - a.dxy * inv2,
             Real(2.0) * a.dy * a.dy * inv3 - a.dyy * inv2);
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
  return Jet(a.v * b, a.dx * b, a.dy * b, a.dt * b,
             a.dxx * b, a.dxy * b, a.dyy * b);
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
      sine, cosine * a.dx, cosine * a.dy, cosine * a.dt,
      cosine * a.dxx - sine * a.dx * a.dx,
      cosine * a.dxy - sine * a.dx * a.dy,
      cosine * a.dyy - sine * a.dy * a.dy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_cos(const Jet& a) noexcept {
  const Real sine = ::sin(a.v);
  const Real cosine = ::cos(a.v);
  return Jet(
      cosine, -sine * a.dx, -sine * a.dy, -sine * a.dt,
      -sine * a.dxx - cosine * a.dx * a.dx,
      -sine * a.dxy - cosine * a.dx * a.dy,
      -sine * a.dyy - cosine * a.dy * a.dy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_sqrt(const Jet& a) noexcept {
  const Real root = ::sqrt(a.v);
  const Real first = Real(0.5) / root;
  const Real second = -Real(0.25) / (a.v * root);
  return Jet(
      root, first * a.dx, first * a.dy, first * a.dt,
      first * a.dxx + second * a.dx * a.dx,
      first * a.dxy + second * a.dx * a.dy,
      first * a.dyy + second * a.dy * a.dy);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
mms_sutherland_law(const Jet& temperature, const Real reference_value,
                   const Real sutherland_constant) noexcept {
  const Jet ratio = temperature / mms_temperature_ref;
  return reference_value * ratio * jet_sqrt(ratio) *
         (mms_temperature_ref + sutherland_constant) /
         (temperature + sutherland_constant);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
x_derivative(const Jet& a) noexcept {
  return Jet(a.dx, a.dxx, a.dxy, Real(0.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
y_derivative(const Jet& a) noexcept {
  return Jet(a.dy, a.dxy, a.dyy, Real(0.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
manufactured_point(const Real x_value, const Real y_value,
                   const Real time, const ProbParm& p,
                   Jet* conserved, Jet* flux_x, Jet* flux_y,
                   Jet* viscous_x, Jet* viscous_y) noexcept {
  constexpr Real pi = 3.141592653589793238462643383279502884;
  const Jet x(x_value, Real(1.0), Real(0.0), Real(0.0));
  const Jet y(y_value, Real(0.0), Real(1.0), Real(0.0));
  const Jet t(time, Real(0.0), Real(0.0), Real(1.0));
  const Real k = Real(2.0) * pi;

  const Jet rho =
      p.rho0 +
      p.rho_x_amp * jet_sin(k * x - Real(0.70) * p.omega * t) +
      p.rho_y_amp * jet_cos(k * y + Real(0.30) * p.omega * t);
  const Jet u =
      p.u0 +
      p.u_x_amp * jet_sin(k * x + Real(0.20) * p.omega * t) +
      p.u_y_amp * jet_cos(k * y - Real(0.40) * p.omega * t);
  const Jet v =
      p.v0 +
      p.v_x_amp * jet_cos(k * x - Real(0.60) * p.omega * t) +
      p.v_y_amp * jet_sin(k * y + Real(0.50) * p.omega * t);
  const Jet pressure =
      p.p0 +
      p.p_x_amp * jet_cos(k * x + Real(0.30) * p.omega * t) +
      p.p_y_amp * jet_sin(k * y - Real(0.20) * p.omega * t);

  const Jet energy =
      pressure / (mms_gamma - Real(1.0)) +
      Real(0.5) * rho * (u * u + v * v);

  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    conserved[n] = Jet();
    flux_x[n] = Jet();
    flux_y[n] = Jet();
    viscous_x[n] = Jet();
    viscous_y[n] = Jet();
  }

  conserved[ProbClosures::URHO] = rho;
  conserved[ProbClosures::UMX] = rho * u;
  conserved[ProbClosures::UMY] = rho * v;
  conserved[ProbClosures::UMZ] = Jet();
  conserved[ProbClosures::UET] = energy;

  flux_x[ProbClosures::URHO] = rho * u;
  flux_x[ProbClosures::UMX] = rho * u * u + pressure;
  flux_x[ProbClosures::UMY] = rho * u * v;
  flux_x[ProbClosures::UET] = u * (energy + pressure);

  flux_y[ProbClosures::URHO] = rho * v;
  flux_y[ProbClosures::UMX] = rho * u * v;
  flux_y[ProbClosures::UMY] = rho * v * v + pressure;
  flux_y[ProbClosures::UET] = v * (energy + pressure);

  const Jet ux = x_derivative(u);
  const Jet uy = y_derivative(u);
  const Jet vx = x_derivative(v);
  const Jet vy = y_derivative(v);
  const Jet theta = ux + vy;
  const Jet temperature = pressure / (rho * mms_rspec);
#if MMS_VARIABLE_TRANSPORT
  const Jet viscosity =
      mms_sutherland_law(temperature, mms_mu,
                         mms_viscosity_sutherland);
  const Jet conductivity =
      mms_sutherland_law(temperature, mms_conductivity,
                         mms_conductivity_sutherland);
  const Jet tau_xx =
      viscosity * (Real(2.0) * ux - Real(2.0 / 3.0) * theta);
  const Jet tau_xy = viscosity * (uy + vx);
  const Jet tau_yy =
      viscosity * (Real(2.0) * vy - Real(2.0 / 3.0) * theta);
#else
  const Jet tau_xx =
      mms_mu * (Real(2.0) * ux - Real(2.0 / 3.0) * theta);
  const Jet tau_xy = mms_mu * (uy + vx);
  const Jet tau_yy =
      mms_mu * (Real(2.0) * vy - Real(2.0 / 3.0) * theta);
#endif
  const Jet tx = x_derivative(temperature);
  const Jet ty = y_derivative(temperature);

  viscous_x[ProbClosures::UMX] = tau_xx;
  viscous_x[ProbClosures::UMY] = tau_xy;
  viscous_x[ProbClosures::UET] =
#if MMS_VARIABLE_TRANSPORT
      u * tau_xx + v * tau_xy + conductivity * tx;
#else
      u * tau_xx + v * tau_xy + mms_conductivity * tx;
#endif

  viscous_y[ProbClosures::UMX] = tau_xy;
  viscous_y[ProbClosures::UMY] = tau_yy;
  viscous_y[ProbClosures::UET] =
#if MMS_VARIABLE_TRANSPORT
      u * tau_xy + v * tau_yy + conductivity * ty;
#else
      u * tau_xy + v * tau_yy + mms_conductivity * ty;
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
point_state(const Real x, const Real y, const Real time,
            const ProbParm& p, Real* state) noexcept {
  Jet conserved[ProbClosures::NCONS];
  Jet flux_x[ProbClosures::NCONS];
  Jet flux_y[ProbClosures::NCONS];
  Jet viscous_x[ProbClosures::NCONS];
  Jet viscous_y[ProbClosures::NCONS];
  manufactured_point(x, y, time, p, conserved, flux_x, flux_y,
                     viscous_x, viscous_y);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state[n] = conserved[n].v;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
point_source(const Real x, const Real y, const Real time,
             const ProbParm& p, Real* source) noexcept {
  Jet conserved[ProbClosures::NCONS];
  Jet flux_x[ProbClosures::NCONS];
  Jet flux_y[ProbClosures::NCONS];
  Jet viscous_x[ProbClosures::NCONS];
  Jet viscous_y[ProbClosures::NCONS];
  manufactured_point(x, y, time, p, conserved, flux_x, flux_y,
                     viscous_x, viscous_y);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    source[n] =
        conserved[n].dt + flux_x[n].dx + flux_y[n].dy;
#if MMS_VISCOUS
    source[n] -= viscous_x[n].dx + viscous_y[n].dy;
#endif
  }
}

inline void inputs() {
  const ProbParm p;
  if (p.rho0 <= amrex::Math::abs(p.rho_x_amp) +
                    amrex::Math::abs(p.rho_y_amp) ||
      p.p0 <= amrex::Math::abs(p.p_x_amp) +
                  amrex::Math::abs(p.p_y_amp)) {
    amrex::Abort("Cartesian MMS parameters do not guarantee positive rho and p");
  }
  amrex::Print()
      << "Cartesian time-dependent point-sample MMS, physics="
#if MMS_VISCOUS
      << "Navier-Stokes, transport="
#if MMS_VARIABLE_TRANSPORT
      << "nondimensional Sutherland\n";
#else
      << "constant\n";
#endif
#else
      << "Euler\n";
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const&,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real x =
      geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
  const Real y =
      geom.ProbLo(1) + (Real(j) + Real(0.5)) * geom.CellSize(1);
  Real exact[ProbClosures::NCONS];
  point_state(x, y, Real(0.0), p, exact);
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
      const Real x =
          geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
      const Real y =
          geom.ProbLo(1) + (Real(j) + Real(0.5)) * geom.CellSize(1);
      Real source[ProbClosures::NCONS];
      point_source(x, y, time, p, source);
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
             const GeomData& geomdata, const ProbParm& p, int) {
  const Real x =
      geomdata.ProbLo(0) + (Real(i) + Real(0.5)) * geomdata.CellSize(0);
  const Real y =
      geomdata.ProbLo(1) + (Real(j) + Real(0.5)) * geomdata.CellSize(1);
  if (x > p.refine_xlo && x < p.refine_xhi && y > p.refine_ylo &&
      y < p.refine_yhi) {
    tagfab(i, j, k) = true;
  }
}

}  // namespace PROB

#endif
