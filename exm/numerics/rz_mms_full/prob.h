#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>
#include <Weno_old.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  Real rho0 = 1.0;
  Real p0 = 1.0;
  Real uz0 = 0.30;
  Real rho_amp = 0.050;
  Real ur_amp = 0.080;
  Real uz_amp = 0.050;
  Real p_amp = 0.080;
  Real rho_r2 = 0.20;
  Real ur_r2 = 0.15;
  Real uz_r2 = 0.10;
  Real p_r2 = 0.25;
  Real omega = 1.0;
  Real uz_phase = 0.35;
  Real p_phase = 0.21;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("rho0", rho0);
    pp.query("p0", p0);
    pp.query("uz0", uz0);
    pp.query("rho_amp", rho_amp);
    pp.query("ur_amp", ur_amp);
    pp.query("uz_amp", uz_amp);
    pp.query("p_amp", p_amp);
    pp.query("rho_r2", rho_r2);
    pp.query("ur_r2", ur_r2);
    pp.query("uz_r2", uz_r2);
    pp.query("p_r2", p_r2);
    pp.query("omega", omega);
    pp.query("uz_phase", uz_phase);
    pp.query("p_phase", p_phase);
  }
};

using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                calorifically_perfect_gas_t<indicies_t>>;

template <typename cls_t>
class user_source_t;

using ProbRHS =
    rhs_dt<weno_old_t<OldReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,
           user_source_t<ProbClosures>>;

// A value and its exact first derivatives with respect to (r,z,t).
struct Jet {
  Real v = 0.0;
  Real dr = 0.0;
  Real dz = 0.0;
  Real dt = 0.0;

  AMREX_GPU_HOST_DEVICE Jet() = default;
  AMREX_GPU_HOST_DEVICE explicit Jet(const Real value) : v(value) {}
  AMREX_GPU_HOST_DEVICE Jet(const Real value, const Real radial,
                            const Real axial, const Real temporal)
      : v(value), dr(radial), dz(axial), dt(temporal) {}
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator+(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v + b.v, a.dr + b.dr, a.dz + b.dz, a.dt + b.dt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v - b.v, a.dr - b.dr, a.dz - b.dz, a.dt - b.dt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator-(const Jet& a) noexcept {
  return Jet(-a.v, -a.dr, -a.dz, -a.dt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
operator*(const Jet& a, const Jet& b) noexcept {
  return Jet(a.v * b.v, a.dr * b.v + a.v * b.dr,
             a.dz * b.v + a.v * b.dz,
             a.dt * b.v + a.v * b.dt);
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
  return Jet(a.v * b, a.dr * b, a.dz * b, a.dt * b);
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
  return Jet(sine, cosine * a.dr, cosine * a.dz, cosine * a.dt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Jet
jet_cos(const Jet& a) noexcept {
  const Real sine = ::sin(a.v);
  const Real cosine = ::cos(a.v);
  return Jet(cosine, -sine * a.dr, -sine * a.dz, -sine * a.dt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
gauss_node(const int q) noexcept {
  constexpr Real x0 = 0.86113631159405257522;
  constexpr Real x1 = 0.33998104358485626480;
  return q == 0 ? -x0 : (q == 1 ? -x1 : (q == 2 ? x1 : x0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
gauss_weight(const int q) noexcept {
  constexpr Real w0 = 0.34785484513745385737;
  constexpr Real w1 = 0.65214515486254614263;
  return (q == 0 || q == 3) ? w0 : w1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
manufactured_point(const Real radial, const Real axial, const Real time,
                   const ProbParm& p, Jet* conserved, Jet* radial_flux,
                   Jet* axial_flux, Jet& pressure) noexcept {
  constexpr Real pi = 3.141592653589793238462643383279502884;
  const Jet r(radial, Real(1.0), Real(0.0), Real(0.0));
  const Jet z(axial, Real(0.0), Real(1.0), Real(0.0));
  const Jet t(time, Real(0.0), Real(0.0), Real(1.0));
  const Jet r2 = r * r;
  const Jet theta = Real(2.0) * pi * z - p.omega * t;

  const Jet rho = p.rho0 +
                  p.rho_amp * (Real(1.0) + p.rho_r2 * r2) * jet_sin(theta);
  const Jet ur = p.ur_amp * r * (Real(1.0) + p.ur_r2 * r2) *
                 jet_cos(theta);
  const Jet uz = p.uz0 +
                 p.uz_amp * (Real(1.0) + p.uz_r2 * r2) *
                     jet_sin(theta + p.uz_phase);
  pressure = p.p0 +
             p.p_amp * (Real(1.0) + p.p_r2 * r2) *
                 jet_cos(theta - p.p_phase);

  const Jet energy = pressure / (ProbClosures::gamma - Real(1.0)) +
                     Real(0.5) * rho * (ur * ur + uz * uz);

  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    conserved[n] = Jet();
    radial_flux[n] = Jet();
    axial_flux[n] = Jet();
  }

  conserved[ProbClosures::URHO] = rho;
  conserved[ProbClosures::UMX] = rho * ur;
  conserved[ProbClosures::UMY] = rho * uz;
  conserved[ProbClosures::UMZ] = Jet();
  conserved[ProbClosures::UET] = energy;

  radial_flux[ProbClosures::URHO] = rho * ur;
  radial_flux[ProbClosures::UMX] = rho * ur * ur + pressure;
  radial_flux[ProbClosures::UMY] = rho * ur * uz;
  radial_flux[ProbClosures::UMZ] = Jet();
  radial_flux[ProbClosures::UET] = ur * (energy + pressure);

  axial_flux[ProbClosures::URHO] = rho * uz;
  axial_flux[ProbClosures::UMX] = rho * ur * uz;
  axial_flux[ProbClosures::UMY] = rho * uz * uz + pressure;
  axial_flux[ProbClosures::UMZ] = Jet();
  axial_flux[ProbClosures::UET] = uz * (energy + pressure);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
point_state(const Real r, const Real z, const Real time, const ProbParm& p,
            Real* state) noexcept {
  Jet conserved[ProbClosures::NCONS];
  Jet radial_flux[ProbClosures::NCONS];
  Jet axial_flux[ProbClosures::NCONS];
  Jet pressure;
  manufactured_point(r, z, time, p, conserved, radial_flux, axial_flux,
                     pressure);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state[n] = conserved[n].v;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
point_source(const Real r, const Real z, const Real time, const ProbParm& p,
             Real* source) noexcept {
  Jet conserved[ProbClosures::NCONS];
  Jet radial_flux[ProbClosures::NCONS];
  Jet axial_flux[ProbClosures::NCONS];
  Jet pressure;
  manufactured_point(r, z, time, p, conserved, radial_flux, axial_flux,
                     pressure);

  const Real inv_r = Real(1.0) / r;
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    // The solver already adds +p/r to radial momentum. Subtracting that
    // geometric term analytically leaves (rho*u_r^2)/r and avoids cancellation
    // of two individually singular pressure terms near the axis.
    const Real metric_flux =
        n == ProbClosures::UMX
            ? (radial_flux[n].v - pressure.v) * inv_r
            : radial_flux[n].v * inv_r;
    source[n] = conserved[n].dt + radial_flux[n].dr + metric_flux +
                axial_flux[n].dz;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cylindrical_average(const Real r_center, const Real z_center, const Real dr,
                    const Real dz, const Real time, const ProbParm& p,
                    Real* average, const bool source_average) noexcept {
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    average[n] = Real(0.0);
  }
  Real normalization = Real(0.0);

  for (int qr = 0; qr < 4; ++qr) {
    const Real r = r_center + Real(0.5) * dr * gauss_node(qr);
    const Real wr = gauss_weight(qr);
    for (int qz = 0; qz < 4; ++qz) {
      const Real z = z_center + Real(0.5) * dz * gauss_node(qz);
      const Real weight = wr * gauss_weight(qz) * r;
      Real point[ProbClosures::NCONS];
      if (source_average) {
        point_source(r, z, time, p, point);
      } else {
        point_state(r, z, time, p, point);
      }
      normalization += weight;
      for (int n = 0; n < ProbClosures::NCONS; ++n) {
        average[n] += weight * point[n];
      }
    }
  }

  const Real inverse_normalization = Real(1.0) / normalization;
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    average[n] *= inverse_normalization;
  }
}

inline void inputs() {
  const ProbParm p;
  if (p.rho0 <= amrex::Math::abs(p.rho_amp) *
                    (Real(1.0) + amrex::Math::abs(p.rho_r2)) ||
      p.p0 <= amrex::Math::abs(p.p_amp) *
                  (Real(1.0) + amrex::Math::abs(p.p_r2))) {
    amrex::Abort("rz_mms_full parameters do not guarantee positive rho and p");
  }
  amrex::Print()
      << "R-Z full Euler time-dependent MMS (cylindrical averages)\n";
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const&,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real r = geom.ProbLo(0) +
                 (Real(i) + Real(0.5)) * geom.CellSize(0);
  const Real z = geom.ProbLo(1) +
                 (Real(j) + Real(0.5)) * geom.CellSize(1);
  Real average[ProbClosures::NCONS];
  cylindrical_average(r, z, geom.CellSize(0), geom.CellSize(1), Real(0.0),
                      p, average, false);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = average[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int, const int,
         const Real time, GeometryData const& geom, ProbClosures const&,
         ProbParm const& p) {
  cylindrical_average(x[0], x[1], geom.CellSize(0), geom.CellSize(1), time, p,
                      s_ext, false);
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
      const Real r = geom.ProbLo(0) +
                     (Real(i) + Real(0.5)) * geom.CellSize(0);
      const Real z = geom.ProbLo(1) +
                     (Real(j) + Real(0.5)) * geom.CellSize(1);
      Real average[ProbClosures::NCONS];
      cylindrical_average(r, z, geom.CellSize(0), geom.CellSize(1), time, p,
                          average, true);
      for (int n = 0; n < ProbClosures::NCONS; ++n) {
        rhs(i, j, k, n) += average[n];
      }
    });
  }

 private:
  ProbParm parameters_{};
};

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, TagFab&, const StateFab&, const GeomData&,
             const ProbParm&, int) {}

}  // namespace PROB

#endif
