#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>
#include <bc_types.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  int mode = 0;     // 0 normal, 1 acoustic blob, 2 acoustic train, 3 entropy,
                    // 4 vortex, 5 mixed, 6 weak shock front, 7 wake packet,
                    // 8 broad low-frequency pressure pulse
  // 2 is the legacy plane-wave acoustic extrapolation prototype; 3 calls the
  // generic ghost-fill LODI helper; 4 applies RHS-level LODI/NSCBC.
  int bc_mode = 1;  // 0 Dirichlet, 1 characteristic, 2 plane, 3 ghost LODI, 4 RHS LODI

  Real gamma = 1.4;
  Real rho0 = 1.0;
  Real p0 = 1.0;
  Real mach = 0.30;
  Real amp = 1.0e-4;  // velocity perturbation amplitude

  Real theta_deg = 30.0;
  Real x0 = 0.30;
  Real y0 = 0.50;
  Real sigma_x = 0.045;
  Real sigma_y = 0.090;
  int carrier_y = 4;
  Real lodi_relax = 0.0;
  Real lodi_transverse = 1.0;
  Real lodi_shock_pressure_jump = 1.0e-2;
  Real lodi_shock_compression = 1.0e-3;
  Real lodi_acoustic_coherence = 0.35;
  Real lodi_convective_velocity_jump = 1.0e-5;
  Real lodi_convective_pressure_fraction = 0.25;
  Real lodi_pressure_relax_min_weight = 0.10;
  int lodi_use_shock_sensor = 1;
  int lodi_use_convective_sensor = 1;
  int lodi_use_convective_extrapolation = 0;
  int lodi_use_pressure_relax_sensor = 1;

  Real sponge_strength = 0.0;
  Real sponge_width = 0.20;
  Real sponge_power = 2.0;
  int sponge_x_lo = 0;
  int sponge_x_hi = 1;
  int sponge_y_lo = 0;
  int sponge_y_hi = 0;

  Real c0 = 0.0;
  Real u0 = 0.0;
  Real theta = 0.0;

  ProbParm()
  {
    ParmParse pp("prob");
    pp.query("mode", mode);
    pp.query("bc_mode", bc_mode);
    pp.query("mach", mach);
    pp.query("amp", amp);
    pp.query("theta_deg", theta_deg);
    pp.query("x0", x0);
    pp.query("y0", y0);
    pp.query("sigma_x", sigma_x);
    pp.query("sigma_y", sigma_y);
    pp.query("carrier_y", carrier_y);
    pp.query("lodi_relax", lodi_relax);
    pp.query("lodi_transverse", lodi_transverse);
    pp.query("lodi_shock_pressure_jump", lodi_shock_pressure_jump);
    pp.query("lodi_shock_compression", lodi_shock_compression);
    pp.query("lodi_acoustic_coherence", lodi_acoustic_coherence);
    pp.query("lodi_convective_velocity_jump", lodi_convective_velocity_jump);
    pp.query("lodi_convective_pressure_fraction", lodi_convective_pressure_fraction);
    pp.query("lodi_pressure_relax_min_weight", lodi_pressure_relax_min_weight);
    pp.query("lodi_use_shock_sensor", lodi_use_shock_sensor);
    pp.query("lodi_use_convective_sensor", lodi_use_convective_sensor);
    pp.query("lodi_use_convective_extrapolation", lodi_use_convective_extrapolation);
    pp.query("lodi_use_pressure_relax_sensor", lodi_use_pressure_relax_sensor);
    pp.query("sponge_strength", sponge_strength);
    pp.query("sponge_width", sponge_width);
    pp.query("sponge_power", sponge_power);
    pp.query("sponge_x_lo", sponge_x_lo);
    pp.query("sponge_x_hi", sponge_x_hi);
    pp.query("sponge_y_lo", sponge_y_lo);
    pp.query("sponge_y_hi", sponge_y_hi);

    c0 = std::sqrt(gamma * p0 / rho0);
    u0 = mach * c0;
    theta = theta_deg * Real(3.14159265358979323846264338327950288) / Real(180.0);
  }
};

struct methodparm_t {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = 0.5;
  static constexpr Real C4skew = 0.016;
};

using ProbClosures = closures_dt<
    indicies_t,
    calorifically_perfect_gas_t<indicies_t>>;
template <typename cls_t> class sponge_source_t;
using ProbRHS = rhs_dt<
    weno_t<ReconScheme::WenoZ5, ProbClosures>,
    no_diffusive_t,
    sponge_source_t<ProbClosures>>;
using GlobalBC = manual_bc_t<ProbClosures>;

inline void inputs() {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_state(const ProbClosures& cls, const ProbParm& pparm,
               const Real rho, const Real u, const Real v, const Real p,
               Real* s)
{
  Real y[NUM_SPECIES] = {Real(1.0)};
  Real eint = Real(0.0);
  cls.RYP2E(rho, y, p, eint);

  s[ProbClosures::URHO] = rho;
  s[ProbClosures::UMX] = rho * u;
  s[ProbClosures::UMY] = rho * v;
  s[ProbClosures::UMZ] = Real(0.0);
  s[ProbClosures::UET] = rho * eint + Real(0.5) * rho * (u * u + v * v);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void load_prims(const IntVect& iv, Array4<Real> const& state,
                const ProbClosures& cls, Real* q)
{
  Real u[ProbClosures::NCONS] = {Real(0.0)};
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    u[n] = state(iv, n);
  }
  cls.cons2prims_point(u, q);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real clamp_real(const Real x, const Real lo, const Real hi)
{
  return amrex::min(hi, amrex::max(lo, x));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real gaussian_train(const Real x, const Real y, const Real lx, const Real ly,
                    const ProbParm& pparm)
{
  const Real ky = Real(2.0) * Real(3.14159265358979323846264338327950288) *
                  Real(pparm.carrier_y) / ly;
  const Real tan_theta = amrex::max(std::tan(pparm.theta), Real(1.e-12));
  const Real kx = ky / tan_theta;
  const Real rx = (x - pparm.x0) / pparm.sigma_x;
  const Real phase = kx * (x - pparm.x0) + ky * (y - pparm.y0);
  return std::exp(-Real(0.5) * rx * rx) * std::cos(phase);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata,
                   ProbClosures const& cls,
                   ProbParm const& pparm)
{
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx = geomdata.CellSize();

  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];

  Real nx = Real(1.0);
  Real ny = Real(0.0);
  Real envelope = Real(0.0);
  const Real lx = prob_hi[0] - prob_lo[0];
  const Real ly = prob_hi[1] - prob_lo[1];

  if (pparm.mode == 0) {
    const Real rx = (x - pparm.x0) / pparm.sigma_x;
    envelope = std::exp(-Real(0.5) * rx * rx);
  } else if (pparm.mode == 1) {
    nx = std::cos(pparm.theta);
    ny = std::sin(pparm.theta);

    Real dy = y - pparm.y0;
    if (dy > Real(0.5) * ly) dy -= ly;
    if (dy < -Real(0.5) * ly) dy += ly;

    const Real dx0 = x - pparm.x0;
    const Real xi = dx0 * nx + dy * ny;
    const Real eta = -dx0 * ny + dy * nx;
    const Real rxi = xi / pparm.sigma_x;
    const Real reta = eta / pparm.sigma_y;
    envelope = std::exp(-Real(0.5) * (rxi * rxi + reta * reta));
  } else if (pparm.mode == 2) {
    nx = std::cos(pparm.theta);
    ny = std::sin(pparm.theta);
    envelope = gaussian_train(x, y, lx, ly, pparm);
  } else {
    nx = std::cos(pparm.theta);
    ny = std::sin(pparm.theta);
    envelope = gaussian_train(x, y, lx, ly, pparm);
  }

  const Real eps = pparm.amp / amrex::max(pparm.c0, Real(1.e-12));
  Real drho = Real(0.0);
  Real dp = Real(0.0);
  Real du = Real(0.0);
  Real dv = Real(0.0);

  if (pparm.mode <= 2) {
    du = pparm.amp * envelope;
    dp = pparm.rho0 * pparm.c0 * du;
    drho = dp / (pparm.c0 * pparm.c0);
    dv = du * ny;
    du *= nx;
  } else if (pparm.mode == 3) {
    drho = pparm.rho0 * eps * envelope;
  } else if (pparm.mode == 4) {
    du = -pparm.amp * ny * envelope;
    dv =  pparm.amp * nx * envelope;
  } else if (pparm.mode == 5) {
    const Real acoustic = Real(0.50) * pparm.amp * envelope;
    dp = pparm.rho0 * pparm.c0 * acoustic;
    drho = dp / (pparm.c0 * pparm.c0) + Real(0.35) * pparm.rho0 * eps * envelope;
    du = acoustic * nx - Real(0.35) * pparm.amp * ny * envelope;
    dv = acoustic * ny + Real(0.35) * pparm.amp * nx * envelope;
  } else if (pparm.mode == 6) {
    const Real width = amrex::max(pparm.sigma_x, Real(2.e-3));
    const Real front = Real(0.5) * (Real(1.0) - std::tanh((x - pparm.x0) / width));
    dp = pparm.p0 * eps * front;
    drho = dp / (pparm.c0 * pparm.c0);
    du = dp / (pparm.rho0 * pparm.c0);
  } else if (pparm.mode == 7) {
    const Real dx0 = (x - pparm.x0) / pparm.sigma_x;
    Real dy = y - pparm.y0;
    if (dy > Real(0.5) * ly) dy -= ly;
    if (dy < -Real(0.5) * ly) dy += ly;
    const Real dy0 = dy / pparm.sigma_y;
    const Real wake = std::exp(-Real(0.5) * (dx0 * dx0 + dy0 * dy0));
    du = -pparm.amp * wake;
    drho = Real(0.20) * pparm.rho0 * eps * wake;
  } else if (pparm.mode == 8) {
    const Real dx0 = (x - pparm.x0) / pparm.sigma_x;
    Real dy = y - pparm.y0;
    if (dy > Real(0.5) * ly) dy -= ly;
    if (dy < -Real(0.5) * ly) dy += ly;
    const Real dy0 = dy / pparm.sigma_y;
    const Real pulse = std::exp(-Real(0.5) * (dx0 * dx0 + dy0 * dy0));
    dp = pparm.rho0 * pparm.c0 * pparm.amp * pulse;
    drho = dp / (pparm.c0 * pparm.c0);
    du = Real(0.25) * dp / (pparm.rho0 * pparm.c0);
  }

  const Real rho = pparm.rho0 + drho;
  const Real p = pparm.p0 + dp;
  const Real u = pparm.u0 + du;
  const Real v = dv;

  state(i, j, k, ProbClosures::URHO) = rho;
  state(i, j, k, ProbClosures::UMX) = rho * u;
  state(i, j, k, ProbClosures::UMY) = rho * v;
  state(i, j, k, ProbClosures::UMZ) = Real(0.0);
  state(i, j, k, ProbClosures::UET) =
      p / (pparm.gamma - Real(1.0)) + Real(0.5) * rho * (u * u + v * v);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real /*x*/[AMREX_SPACEDIM], Real /*dratio*/,
              const Real s_int[ProbClosures::NCONS],
              const Real /*s_refl*/[ProbClosures::NCONS],
              Real s_ext[ProbClosures::NCONS],
              const int idir, const int sgn, const Real /*time*/,
              GeometryData const& /*geomdata*/,
              ProbClosures const& closures,
              ProbParm const& pparm)
{
  if (idir == 0 && sgn == 1) {
    set_state(closures, pparm, pparm.rho0, pparm.u0, Real(0.0), pparm.p0, s_ext);
    return;
  }

  if (idir == 0 && sgn == -1 &&
      (pparm.bc_mode == 1 || pparm.bc_mode == 2 ||
       pparm.bc_mode == 3 || pparm.bc_mode == 4)) {
    GlobalBC::bc_nscbc_farfield(Real(1.0), Real(0.0), Real(0.0), &closures,
                                pparm.rho0, pparm.u0, Real(0.0), Real(0.0),
                                pparm.p0, s_int, s_ext);
    return;
  }

  set_state(closures, pparm, pparm.rho0, pparm.u0, Real(0.0), pparm.p0, s_ext);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool bcnormal_lodi(const IntVect& iv, Array4<Real> const& state,
                   const Real* /*x*/, Real /*dratio*/,
                   const Real* /*s_int*/, const Real* /*s_refl*/, Real* s_ext,
                   const int idir, const int sgn, const Real /*time*/,
                   GeometryData const& geomdata,
                   ProbClosures const& closures,
                   ProbParm const& pparm)
{
#if (AMREX_SPACEDIM < 2)
  return false;
#else
  if ((pparm.bc_mode != 2 && pparm.bc_mode != 3 && pparm.bc_mode != 4) ||
      idir != 0 || sgn != -1) {
    return false;
  }

  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();

  if (pparm.bc_mode == 3 || pparm.bc_mode == 4) {
    GlobalBC::nscbc_lodi_bc_parm_t bp;
    bp.rho_inf = pparm.rho0;
    bp.u_inf = pparm.u0;
    bp.v_inf = Real(0.0);
    bp.w_inf = Real(0.0);
    bp.p_inf = pparm.p0;
    bp.pressure_relax = pparm.lodi_relax;
    bp.transverse_relax = pparm.lodi_transverse;
    bp.l_ref = prob_hi[0] - prob_lo[0];
    bp.shock_pressure_jump = pparm.lodi_shock_pressure_jump;
    bp.shock_compression = pparm.lodi_shock_compression;
    bp.acoustic_coherence_min = pparm.lodi_acoustic_coherence;
    bp.convective_velocity_jump = pparm.lodi_convective_velocity_jump;
    bp.convective_pressure_fraction = pparm.lodi_convective_pressure_fraction;
    bp.pressure_relax_min_weight = pparm.lodi_pressure_relax_min_weight;
    bp.use_transverse = 1;
    bp.use_acoustic_gate = 1;
    bp.use_ghost_projection = 1;
    bp.use_shock_sensor = pparm.lodi_use_shock_sensor;
    bp.use_convective_sensor = pparm.lodi_use_convective_sensor;
    bp.use_convective_extrapolation = pparm.lodi_use_convective_extrapolation;
    bp.use_pressure_relax_sensor = pparm.lodi_use_pressure_relax_sensor;
    bp.max_rel_thermo_change = Real(0.50);
    bp.max_vel_change_c = Real(2.0);
    return GlobalBC::bc_nscbc_lodi_farfield(
        iv, state, s_ext, idir, sgn, geomdata, &closures, bp);
  }

  // Validation-only plane-wave model: infer one local acoustic angle from dp and
  // tangential velocity, then extrapolate primitives along that direction.
  // Mixed acoustic/vortical/entropy/shock content should use bc_mode=3 or a
  // sponge/buffer strategy, not this old matched-wave prototype.
  const int* domhi = geomdata.Domain().hiVect();
  const Real* dx = geomdata.CellSize();

  IntVect first(iv);
  first[0] = domhi[0];

  IntVect jp(first);
  IntVect jm(first);
  jp[1] += 1;
  jm[1] -= 1;

  Real qi[ProbClosures::NPRIM] = {Real(0.0)};
  Real qp[ProbClosures::NPRIM] = {Real(0.0)};
  Real qm[ProbClosures::NPRIM] = {Real(0.0)};
  load_prims(first, state, closures, qi);
  load_prims(jp, state, closures, qp);
  load_prims(jm, state, closures, qm);

  const Real dp = qi[ProbClosures::QPRES] - pparm.p0;
  const Real dv = qi[ProbClosures::QV];
  const Real acoustic_scale = pparm.rho0 * pparm.c0 * pparm.amp;

  if (amrex::Math::abs(dp) < Real(1.e-4) * acoustic_scale) {
    return false;
  }

  Real sin_theta = pparm.rho0 * pparm.c0 * dv / dp;
  sin_theta = clamp_real(sin_theta, Real(-0.98), Real(0.98));
  if (amrex::Math::abs(sin_theta) < Real(0.08)) {
    return false;
  }

  const Real cos_theta =
      std::sqrt(amrex::max(Real(1.0) - sin_theta * sin_theta, Real(1.e-8)));
  const Real cot_theta = cos_theta / sin_theta;
  const Real inv_2dy = Real(0.5) / dx[1];
  const Real dist = Real(iv[0] - first[0]) * dx[0];

  const Real dpdx = cot_theta *
      (qp[ProbClosures::QPRES] - qm[ProbClosures::QPRES]) * inv_2dy;
  const Real dudx = cot_theta *
      (qp[ProbClosures::QU] - qm[ProbClosures::QU]) * inv_2dy;
  const Real dvdx = cot_theta *
      (qp[ProbClosures::QV] - qm[ProbClosures::QV]) * inv_2dy;

  const Real max_dp = Real(20.0) * acoustic_scale;
  const Real p = pparm.p0 + clamp_real(
      qi[ProbClosures::QPRES] + dist * dpdx - pparm.p0, -max_dp, max_dp);
  const Real u = pparm.u0 + clamp_real(
      qi[ProbClosures::QU] + dist * dudx - pparm.u0,
      -Real(20.0) * pparm.amp, Real(20.0) * pparm.amp);
  const Real v = clamp_real(
      qi[ProbClosures::QV] + dist * dvdx,
      -Real(20.0) * pparm.amp, Real(20.0) * pparm.amp);
  const Real rho = pparm.rho0 + (p - pparm.p0) / (pparm.c0 * pparm.c0);

  set_state(closures, pparm,
            amrex::max(rho, Real(1.e-12)), u, v,
            amrex::max(p, Real(1.e-12)), s_ext);
  return true;
#endif
}

inline void rhs_nscbc(const Geometry& geom, const amrex::MFIter& mfi,
                      const auto& prims, const auto& rhs,
                      const ProbClosures* closures, const ProbParm* /*pparm_d*/,
                      amrex::Real /*dt*/, amrex::Real /*time*/)
{
#if (AMREX_SPACEDIM < 2)
  return;
#else
  const ProbParm pparm;
  if (pparm.bc_mode != 4) return;

  const Real* prob_lo = geom.ProbLo();
  const Real* prob_hi = geom.ProbHi();

  GlobalBC::nscbc_lodi_bc_parm_t bp;
  bp.rho_inf = pparm.rho0;
  bp.u_inf = pparm.u0;
  bp.v_inf = Real(0.0);
  bp.w_inf = Real(0.0);
  bp.p_inf = pparm.p0;
  bp.pressure_relax = pparm.lodi_relax;
  bp.transverse_relax = pparm.lodi_transverse;
  bp.l_ref = prob_hi[0] - prob_lo[0];
  bp.shock_pressure_jump = pparm.lodi_shock_pressure_jump;
  bp.shock_compression = pparm.lodi_shock_compression;
  bp.acoustic_coherence_min = pparm.lodi_acoustic_coherence;
  bp.convective_velocity_jump = pparm.lodi_convective_velocity_jump;
  bp.convective_pressure_fraction = pparm.lodi_convective_pressure_fraction;
  bp.pressure_relax_min_weight = pparm.lodi_pressure_relax_min_weight;
  bp.use_transverse = 1;
  bp.use_acoustic_gate = 1;
  bp.use_ghost_projection = 1;
  bp.use_shock_sensor = pparm.lodi_use_shock_sensor;
  bp.use_convective_sensor = pparm.lodi_use_convective_sensor;
  bp.use_convective_extrapolation = 0;
  bp.use_pressure_relax_sensor = pparm.lodi_use_pressure_relax_sensor;
  bp.max_rel_thermo_change = Real(0.50);
  bp.max_vel_change_c = Real(2.0);

  GlobalBC::rhs_nscbc_lodi_farfield(
      mfi.tilebox(), prims, rhs, 0, -1, geom, closures, bp);
#endif
}

template <typename cls_t>
class sponge_source_t {
public:
  void inline src(const Geometry& geomdata, const amrex::MFIter& mfi,
                  const amrex::Array4<const amrex::Real>& prims,
                  const amrex::Array4<amrex::Real>& rhs,
                  const cls_t* /*cls_d*/, amrex::Real /*dt*/,
                  amrex::Real /*time*/)
  {
    const ProbParm pparm;
    if (pparm.sponge_strength <= Real(0.0) ||
        pparm.sponge_width <= Real(0.0)) {
      return;
    }

    const Box bx = mfi.tilebox();
    const auto prob_lo = geomdata.ProbLoArray();
    const auto prob_hi = geomdata.ProbHiArray();
    const auto dx = geomdata.CellSizeArray();
    const Real strength = pparm.sponge_strength;
    const Real width = pparm.sponge_width;
    const Real power = pparm.sponge_power;
    const Real rho_inf = pparm.rho0;
    const Real u_inf = pparm.u0;
    const Real p_inf = pparm.p0;
    const Real e_inf = p_inf / (pparm.gamma - Real(1.0)) +
      Real(0.5) * rho_inf * u_inf * u_inf;
    const int use_x_lo = pparm.sponge_x_lo;
    const int use_x_hi = pparm.sponge_x_hi;
    const int use_y_lo = pparm.sponge_y_lo;
    const int use_y_hi = pparm.sponge_y_hi;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
      Real sigma = Real(0.0);

      if (use_x_lo) {
        const Real xi = (prob_lo[0] + width - x) / width;
        if (xi > Real(0.0)) {
          sigma += strength * std::pow(amrex::min(xi, Real(1.0)), power);
        }
      }
      if (use_x_hi) {
        const Real xi = (x - (prob_hi[0] - width)) / width;
        if (xi > Real(0.0)) {
          sigma += strength * std::pow(amrex::min(xi, Real(1.0)), power);
        }
      }

#if (AMREX_SPACEDIM >= 2)
      const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
      if (use_y_lo) {
        const Real eta = (prob_lo[1] + width - y) / width;
        if (eta > Real(0.0)) {
          sigma += strength * std::pow(amrex::min(eta, Real(1.0)), power);
        }
      }
      if (use_y_hi) {
        const Real eta = (y - (prob_hi[1] - width)) / width;
        if (eta > Real(0.0)) {
          sigma += strength * std::pow(amrex::min(eta, Real(1.0)), power);
        }
      }
#endif

      if (sigma <= Real(0.0)) return;

      const Real rho = prims(i, j, k, cls_t::QRHO);
      const Real u = prims(i, j, k, cls_t::QU);
      const Real v = prims(i, j, k, cls_t::QV);
      const Real w = prims(i, j, k, cls_t::QW);
      const Real p = prims(i, j, k, cls_t::QPRES);
      const Real e = p / (pparm.gamma - Real(1.0)) +
        Real(0.5) * rho * (u * u + v * v + w * w);

      rhs(i, j, k, cls_t::URHO) += -sigma * (rho - rho_inf);
      rhs(i, j, k, cls_t::UMX) += -sigma * (rho * u - rho_inf * u_inf);
      rhs(i, j, k, cls_t::UMY) += -sigma * (rho * v);
      rhs(i, j, k, cls_t::UMZ) += -sigma * (rho * w);
      rhs(i, j, k, cls_t::UET) += -sigma * (e - e_inf);
    });
  }
};

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int, int, int, const auto&, const auto&,
                 const ProbParm&, ProbClosures const&, auto const) {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int, int, int, int, auto&,
                  const auto&, const auto&,
                  const ProbParm&, int) {}

} // namespace PROB

#endif
