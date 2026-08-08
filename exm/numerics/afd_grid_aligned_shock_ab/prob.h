#ifndef CNS_AFD_GRID_ALIGNED_SHOCK_PROB_H_
#define CNS_AFD_GRID_ALIGNED_SHOCK_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

#include <cmath>

using namespace amrex;

namespace PROB {

// A stationary normal shock in a two-dimensional Cartesian channel.  The
// upstream and downstream states satisfy the exact Rankine-Hugoniot relations.
// A small odd-even displacement of the initial shock can be prescribed as a
// fraction of dx.  The displaced discontinuity is integrated exactly within
// the one intersected cell, so subcell perturbations remain measurable.
struct ProbParm {
  Real gamma = Real(1.4);
  Real rho_up = Real(1.0);
  Real p_up = Real(1.0);
  Real mach_up = Real(10.0);
  Real shock_x = Real(0.5);
  Real shock_odd_even_cells = Real(0.01);

  Real sound_up = Real(0.0);
  Real u_up = Real(0.0);
  Real density_ratio = Real(0.0);
  Real pressure_ratio = Real(0.0);
  Real rho_down = Real(0.0);
  Real p_down = Real(0.0);
  Real u_down = Real(0.0);

  ProbParm() {
    ParmParse pp("prob");
    pp.query("mach_up", mach_up);
    pp.query("shock_x", shock_x);
    pp.query("shock_odd_even_cells", shock_odd_even_cells);

    sound_up = std::sqrt(gamma * p_up / rho_up);
    u_up = mach_up * sound_up;
    density_ratio =
        (gamma + Real(1.0)) * mach_up * mach_up /
        ((gamma - Real(1.0)) * mach_up * mach_up + Real(2.0));
    pressure_ratio =
        Real(1.0) + Real(2.0) * gamma / (gamma + Real(1.0)) *
                        (mach_up * mach_up - Real(1.0));
    rho_down = rho_up * density_ratio;
    p_down = p_up * pressure_ratio;
    u_down = u_up / density_ratio;
  }
};

using ProbClosures =
    closures_dt<indicies_t, calorifically_perfect_gas_t<indicies_t>>;
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
using ProbRHS = rhs_dt<ProbEuler, no_diffusive_t, no_source_t>;

inline void inputs() {
  const ProbParm p;
  amrex::Print() << "Cartesian stationary normal shock: M1=" << p.mach_up
                 << " rho2/rho1=" << p.density_ratio
                 << " p2/p1=" << p.pressure_ratio
                 << " odd_even_amplitude_dx="
                 << p.shock_odd_even_cells << '\n';
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
set_normal_shock_state(const bool upstream, Real* state,
                       const ProbClosures& cls, const ProbParm& p) {
  const Real rho = upstream ? p.rho_up : p.rho_down;
  const Real pressure = upstream ? p.p_up : p.p_down;
  const Real u = upstream ? p.u_up : p.u_down;

  state[cls.URHO] = rho;
  state[cls.UMX] = rho * u;
  state[cls.UMY] = Real(0.0);
  state[cls.UMZ] = Real(0.0);
  state[cls.UET] = pressure / (p.gamma - Real(1.0)) +
                   Real(0.5) * rho * u * u;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geomdata, ProbClosures const& cls,
              ProbParm const& p) {
  amrex::ignore_unused(k);
  const Real dx = geomdata.CellSize(0);
  const Real cell_lo = geomdata.ProbLo(0) + Real(i) * dx;
  const Real cell_hi = cell_lo + dx;
  const Real row_sign = (j & 1) == 0 ? Real(-1.0) : Real(1.0);
  const Real local_shock =
      p.shock_x + row_sign * p.shock_odd_even_cells * dx;
  const Real upstream_fraction =
      amrex::max(Real(0.0),
                 amrex::min(Real(1.0), (local_shock - cell_lo) / dx));

  Real upstream[ProbClosures::NCONS] = {Real(0.0)};
  Real downstream[ProbClosures::NCONS] = {Real(0.0)};
  set_normal_shock_state(true, upstream, cls, p);
  set_normal_shock_state(false, downstream, cls, p);
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    state(i, j, k, n) =
        upstream_fraction * upstream[n] +
        (Real(1.0) - upstream_fraction) * downstream[n];
  }
}

// Both streamwise boundaries are held at their exact stationary-shock state.
// The transverse direction is periodic and therefore never enters this
// function.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real[AMREX_SPACEDIM], Real,
         const Real[ProbClosures::NCONS],
         const Real[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int idir, const int sgn, Real,
         GeometryData const&, ProbClosures const& cls, ProbParm const& p) {
  if (idir == 0) {
    set_normal_shock_state(sgn > 0, s_ext, cls, p);
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int, int, int, const auto&, const auto&, const ProbParm&,
            ProbClosures const&, auto const) {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int, int, int, int, auto&, const auto&, const auto&,
             const ProbParm&, int) {}

}  // namespace PROB

#endif
