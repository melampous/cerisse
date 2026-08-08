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

// Standard Shu--Osher shock--entropy-wave interaction.  Defaults preserve the
// historical shifted [0,10] setup; inputs.thesis selects the canonical
// [-5,5] coordinates and x_s=-4.
struct ProbParm {
  Real shock_x = Real(1.0);
  Real rho_left = Real(3.857142857142857);
  Real u_left = Real(2.629369);
  Real p_left = Real(10.333333333333333);
  Real rho_wave_amplitude = Real(0.2);
  Real wave_number = Real(5.0);
  Real p_right = Real(1.0);

  ProbParm() {
    ParmParse pp("prob");
    pp.query("shock_x", shock_x);
    pp.query("rho_left", rho_left);
    pp.query("u_left", u_left);
    pp.query("p_left", p_left);
    pp.query("rho_wave_amplitude", rho_wave_amplitude);
    pp.query("wave_number", wave_number);
    pp.query("p_right", p_right);
  }
};

// numerical method parameters
#ifndef SHU_SKEW_C2_VALUE
#define SHU_SKEW_C2_VALUE 1.5
#endif

#ifndef SHU_SKEW_C4_VALUE
#define SHU_SKEW_C4_VALUE 0.016
#endif

struct methodparm_t {

  public:

  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = Real(SHU_SKEW_C2_VALUE);
  static constexpr Real C4skew = Real(SHU_SKEW_C4_VALUE);

};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

using ProbClosures =
    closures_dt<indicies_t, calorifically_perfect_gas_t<indicies_t>>;

#ifndef SHU_EULER_SCHEME_ID
#define SHU_EULER_SCHEME_ID 0
#endif

#if (SHU_EULER_SCHEME_ID == 0)
using ProbEuler = skew_t<methodparm_t, ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 1)
using ProbEuler = afd_hllc_teno5_t<ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 2)
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 3)
using ProbEuler = weno_t<ReconScheme::Teno5, ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 4)
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 5)
using ProbEuler = weno_t<ReconScheme::Teno6, ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 6)
using ProbEuler = weno_old_t<OldReconScheme::Teno5, ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 7)
using ProbEuler = weno_old_t<OldReconScheme::WenoZ5, ProbClosures>;
#elif (SHU_EULER_SCHEME_ID == 8)
using ProbEuler = weno_old_t<OldReconScheme::Teno6, ProbClosures>;
#else
#error "Unknown SHU_EULER_SCHEME_ID"
#endif

typedef rhs_dt<ProbEuler, no_diffusive_t, no_source_t> ProbRHS;

// The RK scheme is stated explicitly in each input file so that the archived
// input remains a complete and auditable record of the run.
inline void inputs() {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_euler_state(Real *state, const Real rho, const Real u, const Real p,
                ProbClosures const &cls) {
  state[cls.URHO] = rho;
  state[cls.UMX] = rho * u;
  state[cls.UMY] = Real(0.0);
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = p / (cls.gamma - Real(1.0))
                 + Real(0.5) * rho * u * u;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];

  const bool left_of_shock = x <= prob_parm.shock_x;
  const Real rho = left_of_shock
      ? prob_parm.rho_left
      : Real(1.0) + prob_parm.rho_wave_amplitude
                        * std::sin(prob_parm.wave_number * x);
  const Real u = left_of_shock ? prob_parm.u_left : Real(0.0);
  const Real p = left_of_shock ? prob_parm.p_left : prob_parm.p_right;

  Real local_state[ProbClosures::NCONS] = {Real(0.0)};
  set_euler_state(local_state, rho, u, p, cls);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) = local_state[n];
  }
}

// The left boundary is held at the uniform post-shock state.  At the right
// boundary, the undisturbed entropy wave is continued analytically.  The shock
// does not reach x=5 by t=1.8, so this removes the legacy foextrap artefact
// without influencing the interaction region.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real /*dratio*/,
         const Real /*s_int*/[ProbClosures::NCONS],
         const Real /*s_refl*/[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int idir, const int sgn,
         const Real /*time*/, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {
  if (idir != 0) {
    return;
  }

  if (sgn > 0) {
    set_euler_state(s_ext, prob_parm.rho_left, prob_parm.u_left,
                    prob_parm.p_left, closures);
  } else {
    const Real rho = Real(1.0) + prob_parm.rho_wave_amplitude
                                  * std::sin(prob_parm.wave_number * x[0]);
    set_euler_state(s_ext, rho, Real(0.0), prob_parm.p_right, closures);
  }
}

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int i, int j, int k, const auto &state, const auto &rhs,
            const ProbParm &lprobparm, ProbClosures const &closures,
            auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {

}
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
