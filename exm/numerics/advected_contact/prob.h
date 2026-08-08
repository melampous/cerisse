#ifndef CNS_ADVECTED_CONTACT_PROB_H_
#define CNS_ADVECTED_CONTACT_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

// A square density pulse is transported at constant velocity and pressure.
// Its two edges are exact contact discontinuities of the Euler equations.
struct ProbParm {
  Real rho_low = Real(1.0);
  Real rho_high = Real(2.0);
  Real p0 = Real(1.0);
  Real mach_reference = Real(0.5);
  Real pulse_lo = Real(0.25);
  Real pulse_hi = Real(0.75);

  ProbParm() {
    ParmParse pp("prob");
    pp.query("rho_low", rho_low);
    pp.query("rho_high", rho_high);
    pp.query("p0", p0);
    pp.query("mach_reference", mach_reference);
    pp.query("pulse_lo", pulse_lo);
    pp.query("pulse_hi", pulse_hi);
  }
};

using ProbClosures =
    closures_dt<indicies_t, calorifically_perfect_gas_t<indicies_t>>;
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
using ProbRHS = rhs_dt<ProbEuler, no_diffusive_t, no_source_t>;

inline void inputs() {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_euler_state(Real *state, const Real rho, const Real ux, const Real p,
                ProbClosures const &cls) {
  state[cls.URHO] = rho;
  state[cls.UMX] = rho * ux;
  state[cls.UMY] = Real(0.0);
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = p / (cls.gamma - Real(1.0))
                 + Real(0.5) * rho * ux * ux;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &pparm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real rho =
      (x >= pparm.pulse_lo && x < pparm.pulse_hi)
          ? pparm.rho_high
          : pparm.rho_low;
  const Real c_ref = std::sqrt(cls.gamma * pparm.p0 / pparm.rho_low);
  const Real ux = pparm.mach_reference * c_ref;

  Real local_state[ProbClosures::NCONS] = {Real(0.0)};
  set_euler_state(local_state, rho, ux, pparm.p0, cls);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = local_state[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS], Real[ProbClosures::NCONS],
         int, int, Real, GeometryData const &, ProbClosures const &,
         ProbParm const &) {
  amrex::Abort("advected_contact: bcnormal called for a periodic problem");
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int, int, int, const auto &, const auto &, const ProbParm &,
            ProbClosures const &, auto const) {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, auto &, const auto &, const auto &,
             const ProbParm &, int) {}

} // namespace PROB

#endif
