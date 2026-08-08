#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>
#include <Weno_old.h>

#include <cmath>

using namespace amrex;

namespace PROB {

struct ProbParm {
  int profile = 0;
  int cell_average = 1;
  Real gamma = 1.4;
  Real rho0 = 1.0;
  Real p0 = 2.0;
  Real amplitude = 0.25;
  Real quartic = 0.10;
  Real shock_z = 0.5;
  Real curvature = 0.20;
  Real width = 0.04;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("profile", profile);
    pp.query("cell_average", cell_average);
    pp.query("p0", p0);
    pp.query("amplitude", amplitude);
    pp.query("quartic", quartic);
    pp.query("shock_z", shock_z);
    pp.query("curvature", curvature);
    pp.query("width", width);
  }
};

using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                calorifically_perfect_gas_t<indicies_t>>;
using ProbRHS =
    rhs_dt<weno_old_t<OldReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,
           no_source_t>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
pressure(const Real r, const Real z, const ProbParm& p) noexcept {
  if (p.profile == 0) return p.p0 + p.amplitude * r * r;
  if (p.profile == 2) {
    const Real r2 = r * r;
    return p.p0 + p.amplitude * r2 + p.quartic * r2 * r2;
  }
  const Real arg =
      (z - p.shock_z - p.curvature * r * r) / p.width;
  return p.p0 + p.amplitude * std::tanh(arg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_state(const Real pressure_average, Real* state,
          const ProbClosures& cls, const ProbParm& p) noexcept {
  state[cls.URHO] = p.rho0;
  state[cls.UMX] = Real(0.0);
  state[cls.UMY] = Real(0.0);
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = pressure_average / (p.gamma - Real(1.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
cell_average_pressure(const Real r, const Real z, const Real dr,
                      const Real dz, const ProbParm& p) noexcept {
  const Real rlo = r - Real(0.5) * dr;
  const Real rhi = r + Real(0.5) * dr;
  if (p.profile == 0 || p.profile == 2) {
    // Cylindrical volume average: integral(p*r*dr)/integral(r*dr).
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

inline void inputs() {
  ProbParm p;
  amrex::Print() << "R-Z pressure RHS regression, profile=" << p.profile
                 << '\n';
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const& cls,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real r = geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
  const Real z = geom.ProbLo(1) + (Real(j) + Real(0.5)) * geom.CellSize(1);
  Real conserved[ProbClosures::NCONS];
  const Real p_state = p.cell_average
      ? cell_average_pressure(r, z, geom.CellSize(0), geom.CellSize(1), p)
      : pressure(r, z, p);
  set_state(p_state, conserved, cls, p);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = conserved[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int, const int, Real,
         GeometryData const& geom, ProbClosures const& cls, ProbParm const& p) {
  const Real p_state = p.cell_average
      ? cell_average_pressure(
            x[0], x[1], geom.CellSize(0), geom.CellSize(1), p)
      : pressure(x[0], x[1], p);
  set_state(p_state, s_ext, cls, p);
}

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, TagFab&, const StateFab&, const GeomData&,
             const ProbParm&, int) {}

}  // namespace PROB

#endif
