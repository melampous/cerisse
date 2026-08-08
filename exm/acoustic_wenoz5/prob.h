#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  int mode = 0;
  Real rho0 = 1.0;
  Real p0 = 1.0;
  Real gamma = 1.4;
  Real amplitude = 1.0e-5;
  int n_waves = 2;
  Real mean_mach = 0.0;
  Real packet_x0 = 0.20;
  Real packet_sigma = 0.05;
  int packet_n_waves = 8;
  Real refine_lo = 0.40;
  Real refine_hi = 0.70;

  ProbParm()
  {
    ParmParse pp("prob");
    pp.query("mode", mode);
    pp.query("rho0", rho0);
    pp.query("p0", p0);
    pp.query("gamma", gamma);
    pp.query("amplitude", amplitude);
    pp.query("n_waves", n_waves);
    pp.query("mean_mach", mean_mach);
    pp.query("packet_x0", packet_x0);
    pp.query("packet_sigma", packet_sigma);
    pp.query("packet_n_waves", packet_n_waves);
    pp.query("refine_lo", refine_lo);
    pp.query("refine_hi", refine_hi);
  }
};

using ProbClosures = closures_dt<
    indicies_t, visc_suth_t, cond_suth_t,
    calorifically_perfect_gas_t<indicies_t>>;

using ProbRHS = rhs_dt<
    weno_t<ReconScheme::WenoZ5, ProbClosures>,
    no_diffusive_t,
    no_source_t>;

inline void inputs()
{
  ParmParse pp;
  pp.add("cns.order_rk", 3);
  pp.add("cns.stages_rk", 3);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata,
                   ProbClosures const& cls,
                   ProbParm const& pparm)
{
  amrex::ignore_unused(j, k);
  const Real* plo = geomdata.ProbLo();
  const Real* phi = geomdata.ProbHi();
  const Real* dx = geomdata.CellSize();
  const Real x = plo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real lx = phi[0] - plo[0];
  const Real c0 = std::sqrt(pparm.gamma * pparm.p0 / pparm.rho0);
  const Real u0 = pparm.mean_mach * c0;

  Real shape = Real(0.0);
  if (pparm.mode == 0) {
    const Real twopi = Real(2.0) * Real(3.14159265358979323846);
    shape = std::sin(twopi * Real(pparm.n_waves) * (x - plo[0]) / lx);
  } else {
    Real distance = x - pparm.packet_x0;
    distance -= lx * std::round(distance / lx);
    const Real xi = distance / pparm.packet_sigma;
    shape = std::exp(-Real(0.5) * xi * xi);
    if (pparm.mode == 2) {
      const Real twopi = Real(2.0) * Real(3.14159265358979323846);
      shape *= std::cos(twopi * Real(pparm.packet_n_waves) *
                        distance / lx);
    }
  }

  // Linear right-running acoustic eigenmode:
  // rho'/rho0 = u'/c0 = p'/(gamma*p0) = amplitude*shape.
  const Real rho = pparm.rho0 * (Real(1.0) + pparm.amplitude * shape);
  const Real u = u0 + pparm.amplitude * c0 * shape;
  const Real p = pparm.p0 *
    (Real(1.0) + pparm.gamma * pparm.amplitude * shape);

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX) = rho * u;
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  state(i, j, k, cls.UET) =
    p / (pparm.gamma - Real(1.0)) + Real(0.5) * rho * u * u;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real[AMREX_SPACEDIM], Real,
              const Real[ProbClosures::NCONS],
              const Real[ProbClosures::NCONS],
              Real[ProbClosures::NCONS], int, int, Real,
              GeometryData const&, ProbClosures const&, ProbParm const&)
{}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int, int, int, const auto&, const auto&,
                 const ProbParm&, ProbClosures const&, auto const)
{}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int, auto& tagfab,
                  const auto&, const auto& geomdata,
                  const ProbParm& pparm, int level)
{
  amrex::ignore_unused(j, k);
  if ((pparm.mode != 1 && pparm.mode != 2) || level != 0) return;
  const Real x = geomdata.ProbLo()[0] +
    (Real(i) + Real(0.5)) * geomdata.CellSize()[0];
  if (x >= pparm.refine_lo && x <= pparm.refine_hi) {
    tagfab(i, j, k) = true;
  }
}

} // namespace PROB

#endif
