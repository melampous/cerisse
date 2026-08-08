#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  int mode = 0;
  int cell_average = 1;
  Real gamma = 1.4;
  Real rho0 = 1.0;
  Real p0 = 10.0;
  Real a = 0.20;
  Real b = 0.10;
  Real c = 0.20;
  Real d = 0.20;
  Real z0 = 0.5;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("mode", mode);
    pp.query("cell_average", cell_average);
    pp.query("a", a);
    pp.query("b", b);
    pp.query("c", c);
    pp.query("d", d);
  }
};

struct methodparm_t {
  static constexpr int order = 2;
#if defined(RZ_VISCOUS_TEST_LES) && RZ_VISCOUS_TEST_LES
  static constexpr bool use_LES = true;
#else
  static constexpr bool use_LES = false;
#endif
  static constexpr Real viscosity = 0.01;
  static constexpr Real conductivity = 0.0;
};

struct lesparm_t {
  static constexpr int order = 2;
  static constexpr Real Pr_o_Prsgs = 0.8;
  static constexpr Real Scsgs = 0.7;
  static constexpr Real Cs = 0.16;
  static constexpr Real CI = 0.08;
  static constexpr bool fixDelta = false;
  static constexpr Real Delta = 0.02;
};

#if defined(RZ_VISCOUS_TEST_LES) && RZ_VISCOUS_TEST_LES
using ProbClosures =
    closures_dt<indicies_t, visc_const_t<methodparm_t>,
                cond_const_t<methodparm_t>,
                calorifically_perfect_gas_t<indicies_t>,
                Smagorinsky_t<lesparm_t, indicies_t>>;
#else
using ProbClosures =
    closures_dt<indicies_t, visc_const_t<methodparm_t>,
                cond_const_t<methodparm_t>,
                calorifically_perfect_gas_t<indicies_t>>;
#endif
using ProbRHS =
    rhs_dt<no_euler_t, viscous_t<methodparm_t, ProbClosures>, no_source_t>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
velocity(const Real r, const Real z, const ProbParm& p,
         Real& ur, Real& uz) noexcept {
  ur = Real(0.0);
  uz = Real(0.0);
  if (p.mode == 0) ur = p.a * r;
  if (p.mode == 1) ur = p.a * r + p.b * r * r * r;
  if (p.mode == 2) uz = p.c * r * r;
  if (p.mode == 3) uz = p.d * (z - p.z0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
cylindrical_moment(const Real r, const Real dr, const int power) noexcept {
  const Real rlo = r - Real(0.5) * dr;
  const Real rhi = r + Real(0.5) * dr;
  const Real denominator = Real(0.5) * (rhi * rhi - rlo * rlo);
  if (power == 1) {
    return (rhi * rhi * rhi - rlo * rlo * rlo) /
           (Real(3.0) * denominator);
  }
  if (power == 2) {
    return (rhi * rhi * rhi * rhi - rlo * rlo * rlo * rlo) /
           (Real(4.0) * denominator);
  }
  return (rhi * rhi * rhi * rhi * rhi -
          rlo * rlo * rlo * rlo * rlo) /
         (Real(5.0) * denominator);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cell_average_velocity(const Real r, const Real z, const Real dr,
                      const ProbParm& p, Real& ur, Real& uz) noexcept {
  if (!p.cell_average) {
    velocity(r, z, p, ur, uz);
    return;
  }

  ur = Real(0.0);
  uz = Real(0.0);
  if (p.mode == 0) ur = p.a * cylindrical_moment(r, dr, 1);
  if (p.mode == 1) {
    ur = p.a * cylindrical_moment(r, dr, 1) +
         p.b * cylindrical_moment(r, dr, 3);
  }
  if (p.mode == 2) uz = p.c * cylindrical_moment(r, dr, 2);
  if (p.mode == 3) uz = p.d * (z - p.z0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_state(const Real r, const Real z, const Real dr, Real* state,
          const ProbClosures& cls, const ProbParm& p) noexcept {
  Real ur;
  Real uz;
  cell_average_velocity(r, z, dr, p, ur, uz);
  state[cls.URHO] = p.rho0;
  state[cls.UMX] = p.rho0 * ur;
  state[cls.UMY] = p.rho0 * uz;
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = p.p0 / (p.gamma - Real(1.0)) +
                   Real(0.5) * p.rho0 * (ur * ur + uz * uz);
}

inline void inputs() {
  ProbParm p;
  amrex::Print() << "R-Z viscous RHS regression, mode=" << p.mode << '\n';
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const& cls,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real r = geom.ProbLo(0) + (Real(i) + Real(0.5)) * geom.CellSize(0);
  const Real z = geom.ProbLo(1) + (Real(j) + Real(0.5)) * geom.CellSize(1);
  Real conserved[ProbClosures::NCONS];
  set_state(r, z, geom.CellSize(0), conserved, cls, p);
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
  set_state(x[0], x[1], geom.CellSize(0), s_ext, cls, p);
}

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, TagFab&, const StateFab&, const GeomData&,
             const ProbParm&, int) {}

}  // namespace PROB

#endif
