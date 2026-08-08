#ifndef CNS_ADVECTED_SHEAR_WAVE_PROB_H_
#define CNS_ADVECTED_SHEAR_WAVE_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

// A small transverse velocity wave is advected by a uniform x velocity.  The
// density and pressure are uniform.  This is an exact solution of the
// compressible Euler equations.
struct ProbParm {
  Real rho0 = Real(1.0);
  Real p0 = Real(1.0);
  Real mach = Real(0.5);
  Real transverse_amplitude_over_c0 = Real(1.0e-5);
  int mode = 1;

  ProbParm() {
    ParmParse pp("prob");
    pp.query("rho0", rho0);
    pp.query("p0", p0);
    pp.query("mach", mach);
    pp.query("transverse_amplitude_over_c0",
             transverse_amplitude_over_c0);
    pp.query("mode", mode);
  }
};

// These values are deliberately local to this verification executable.
// They do not change the production solver defaults.
struct SkewJstParm {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = Real(1.5);
  static constexpr Real C4skew = Real(0.016);
};

using ProbClosures =
    closures_dt<indicies_t, calorifically_perfect_gas_t<indicies_t>>;

#ifndef SHEAR_EULER_SCHEME_ID
#define SHEAR_EULER_SCHEME_ID 0
#endif

#if (SHEAR_EULER_SCHEME_ID == 0)
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
#elif (SHEAR_EULER_SCHEME_ID == 1)
using ProbEuler = weno_t<ReconScheme::Teno5, ProbClosures>;
#elif (SHEAR_EULER_SCHEME_ID == 2)
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
#elif (SHEAR_EULER_SCHEME_ID == 3)
using ProbEuler = skew_t<SkewJstParm, ProbClosures>;
#else
#error "Unknown SHEAR_EULER_SCHEME_ID"
#endif

using ProbRHS = rhs_dt<ProbEuler, no_diffusive_t, no_source_t>;

// Time integration is stated in the archived input file.
inline void inputs() {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_euler_state(Real *state, const Real rho, const Real ux, const Real uy,
                const Real p, ProbClosures const &cls) {
  state[cls.URHO] = rho;
  state[cls.UMX] = rho * ux;
  state[cls.UMY] = rho * uy;
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = p / (cls.gamma - Real(1.0))
                 + Real(0.5) * rho * (ux * ux + uy * uy);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &pparm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real length = prob_hi[0] - prob_lo[0];
  const Real pi = std::acos(Real(-1.0));
  const Real wave_number = Real(2.0) * pi * Real(pparm.mode) / length;
  const Real c0 = std::sqrt(cls.gamma * pparm.p0 / pparm.rho0);
  const Real ux = pparm.mach * c0;
  const Real amplitude = pparm.transverse_amplitude_over_c0 * c0;
  const Real uy = amplitude * std::sin(wave_number * (x - prob_lo[0]));

  Real local_state[ProbClosures::NCONS] = {Real(0.0)};
  set_euler_state(local_state, pparm.rho0, ux, uy, pparm.p0, cls);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = local_state[n];
  }
}

// Every physical direction is periodic.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS], Real[ProbClosures::NCONS],
         int, int, Real, GeometryData const &, ProbClosures const &,
         ProbParm const &) {
  amrex::Abort("advected_shear_wave: bcnormal called for a periodic problem");
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int, int, int, const auto &, const auto &, const ProbParm &,
            ProbClosures const &, auto const) {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, auto &, const auto &, const auto &,
             const ProbParm &, int) {}

} // namespace PROB

#endif
