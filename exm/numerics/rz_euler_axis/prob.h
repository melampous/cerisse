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

  ProbParm() {
    ParmParse pp("prob");
    pp.query("mode", mode);
    pp.query("cell_average", cell_average);
    pp.query("a", a);
    pp.query("b", b);
  }
};

using ProbClosures =
    closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                calorifically_perfect_gas_t<indicies_t>>;
using ProbRHS =
    rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,
           no_source_t>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
cylindrical_moment(const Real r, const Real dr, const int power) noexcept {
  const Real rlo = r - Real(0.5) * dr;
  const Real rhi = r + Real(0.5) * dr;
  const Real denominator = Real(0.5) * (rhi * rhi - rlo * rlo);
  if (power == 1) {
    return (rhi * rhi * rhi - rlo * rlo * rlo) /
           (Real(3.0) * denominator);
  }
  if (power == 3) {
    return (rhi * rhi * rhi * rhi * rhi -
            rlo * rlo * rlo * rlo * rlo) /
           (Real(5.0) * denominator);
  }
  return Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
radial_velocity(const Real r, const Real dr, const ProbParm& p) noexcept {
  if (!p.cell_average) {
    return p.a * r + (p.mode == 1 ? p.b * r * r * r : Real(0.0));
  }
  return p.a * cylindrical_moment(r, dr, 1) +
         (p.mode == 1
              ? p.b * cylindrical_moment(r, dr, 3)
              : Real(0.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_state(const Real r, const Real dr, Real* state,
          const ProbClosures& cls, const ProbParm& p) noexcept {
  const Real ur = radial_velocity(r, dr, p);
  state[cls.URHO] = p.rho0;
  state[cls.UMX] = p.rho0 * ur;
  state[cls.UMY] = Real(0.0);
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = p.p0 / (p.gamma - Real(1.0)) +
                   Real(0.5) * p.rho0 * ur * ur;
}

inline void inputs() {
  ProbParm p;
  amrex::Print() << "R-Z Euler axis regression, mode=" << p.mode << '\n';
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geom, ProbClosures const& cls,
              ProbParm const& p) {
  amrex::ignore_unused(j, k);
  const Real r = geom.ProbLo(0) +
                 (Real(i) + Real(0.5)) * geom.CellSize(0);
  Real conserved[ProbClosures::NCONS];
  set_state(r, geom.CellSize(0), conserved, cls, p);
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
  set_state(x[0], geom.CellSize(0), s_ext, cls, p);
}

template <typename TagFab, typename StateFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, TagFab&, const StateFab&, const GeomData&,
             const ProbParm&, int) {}

}  // namespace PROB

#endif
