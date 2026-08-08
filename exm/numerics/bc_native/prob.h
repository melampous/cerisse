#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <limits>

#include <Closures.h>
#include <RHS.h>
#include <bc_types.h>
#include <nscbc.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  int mode = 0;     // 0 normal, 1 acoustic blob, 2 acoustic train, 3 entropy,
                    // 4 vortex, 5 mixed, 6 weak shock front, 7 wake packet,
                    // 8 broad low-frequency pressure pulse, 9 tangential shear,
                    // 10 analytically divergence-free vortical packet
  // 2 is the legacy plane-wave acoustic extrapolation prototype; 3 calls the
  // generic ghost-fill LODI helper; 4 applies RHS-level LODI/NSCBC.
  int bc_mode = 1;  // 0 Dirichlet, 1 characteristic, 2 plane, 3 ghost LODI, 4 RHS LODI
  // Stage-3 algebra gates.  Zero preserves the historical bc_native setup.
  // inlet_model:  0=fixed static state, 1=subsonic (pt,Tt,direction),
  //               2=supersonic prescribed full state, 3=far-field invariants
  // outlet_model: 0=historical bc_mode selection, 1=subsonic p_back with
  //               automatic supersonic extrapolation, 2=far-field invariants
  int inlet_model = 0;
  int outlet_model = 0;
  int all_farfield = 0;
  Real back_pressure_ratio = Real(1.0);

  Real gamma = 1.4;
  Real rho0 = 1.0;
  Real p0 = 1.0;
  Real mach = 0.30;
  Real amp = 1.0e-4;  // velocity perturbation amplitude

  Real theta_deg = 30.0;
  Real flow_angle_deg = 0.0;
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
  Real v0 = 0.0;
  Real theta = 0.0;

  ProbParm()
  {
    ParmParse pp("prob");
    pp.query("mode", mode);
    pp.query("bc_mode", bc_mode);
    pp.query("inlet_model", inlet_model);
    pp.query("outlet_model", outlet_model);
    pp.query("all_farfield", all_farfield);
    pp.query("back_pressure_ratio", back_pressure_ratio);
    pp.query("mach", mach);
    pp.query("amp", amp);
    pp.query("theta_deg", theta_deg);
    pp.query("flow_angle_deg", flow_angle_deg);
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
    const Real flow_angle = flow_angle_deg *
        Real(3.14159265358979323846264338327950288) / Real(180.0);
    u0 = mach * c0 * std::cos(flow_angle);
    v0 = mach * c0 * std::sin(flow_angle);
    theta = theta_deg * Real(3.14159265358979323846264338327950288) / Real(180.0);

    if (inlet_model < 0 || inlet_model > 3) {
      amrex::Abort("prob.inlet_model must be 0, 1, 2, or 3");
    }
    if (outlet_model < 0 || outlet_model > 2) {
      amrex::Abort("prob.outlet_model must be 0, 1, or 2");
    }
    if (all_farfield < 0 || all_farfield > 1 ||
        !std::isfinite(flow_angle_deg)) {
      amrex::Abort("prob.all_farfield must be 0 or 1 and flow_angle_deg must be finite");
    }
    if (!(back_pressure_ratio > Real(0.0)) ||
        !std::isfinite(back_pressure_ratio)) {
      amrex::Abort("prob.back_pressure_ratio must be finite and positive");
    }
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

inline void inputs()
{
#ifdef AMREX_USE_GPU
  // The production boundary callback executes on the device.  Cerisse's
  // pointwise conservative-to-primitive closure is device-only in CUDA
  // builds, so the host algebra oracle below is exercised by the CPU build;
  // CUDA coverage comes from the stage-local ghost-fill runtime tests.
  amrex::Print()
      << "[Stage3-Characteristic-BC] host algebra self-test: CPU build only\n";
#else
  ProbClosures cls;
  constexpr Real tolerance = Real(2.e-12);
  const Real gamma = calorifically_perfect_gas_t<indicies_t>::gamma;
  const Real sound = std::sqrt(gamma);

  auto make_state = [&](const Real rho, const Real u, const Real v,
                        const Real p, Real* state) {
    Real Y[NUM_SPECIES] = {Real(1.0)};
    Real internal_energy = Real(0.0);
    cls.RYP2E(rho, Y, p, internal_energy);
    state[ProbClosures::URHO] = rho;
    state[ProbClosures::UMX] = rho * u;
    state[ProbClosures::UMY] = rho * v;
    state[ProbClosures::UMZ] = Real(0.0);
    state[ProbClosures::UET] =
        rho * internal_energy + Real(0.5) * rho * (u * u + v * v);
  };
  auto relative_error = [](const Real a, const Real b) {
    return amrex::Math::abs(a - b) /
           amrex::max(Real(1.0), amrex::max(amrex::Math::abs(a),
                                            amrex::Math::abs(b)));
  };
  auto require = [](const bool condition, const char* message) {
    if (!condition) amrex::Abort(message);
  };

  const Real angle = Real(20.0) *
      Real(3.14159265358979323846264338327950288) / Real(180.0);
  const Real mach = Real(0.3);
  const Real speed = mach * sound;
  const Real u = speed * std::cos(angle);
  const Real v = speed * std::sin(angle);
  Real interior[ProbClosures::NCONS] = {Real(0.0)};
  make_state(Real(1.0), u, v, Real(1.0), interior);

  const Real static_temperature = Real(1.0) / cls.Rspec;
  const Real total_ratio =
      Real(1.0) + Real(0.5) * (gamma - Real(1.0)) * mach * mach;
  const Real total_temperature = static_temperature * total_ratio;
  const Real total_pressure =
      std::pow(total_ratio, gamma / (gamma - Real(1.0)));
  Real reconstructed[ProbClosures::NCONS] = {Real(0.0)};
  auto status = GlobalBC::bc_characteristic_total_inlet(
      Real(-1.0), Real(0.0), Real(0.0), &cls,
      total_pressure, total_temperature, std::cos(angle), std::sin(angle),
      Real(0.0), interior, reconstructed);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic total-inlet self-test returned a failure status");
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    require(relative_error(reconstructed[n], interior[n]) < tolerance,
            "characteristic total-inlet self-test failed exact-state recovery");
  }

  Real outlet[ProbClosures::NCONS] = {Real(0.0)};
  const Real back_pressure = Real(0.98);
  status = GlobalBC::bc_characteristic_pressure_outlet(
      Real(1.0), Real(0.0), Real(0.0), &cls, back_pressure,
      interior, outlet);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic pressure-outlet self-test returned a failure status");
  auto decode_ideal_gas = [&](const Real* state, Real& rho, Real& u,
                              Real& v, Real& p, Real& c) {
    rho = state[ProbClosures::URHO];
    u = state[ProbClosures::UMX] / rho;
    v = state[ProbClosures::UMY] / rho;
    const Real w = state[ProbClosures::UMZ] / rho;
    const Real kinetic = Real(0.5) * rho * (u * u + v * v + w * w);
    p = (gamma - Real(1.0)) * (state[ProbClosures::UET] - kinetic);
    c = std::sqrt(gamma * p / rho);
  };
  Real rho_i, u_i, v_i, p_i, c_i;
  Real rho_o, u_o, v_o, p_o, c_o;
  decode_ideal_gas(interior, rho_i, u_i, v_i, p_i, c_i);
  decode_ideal_gas(outlet, rho_o, u_o, v_o, p_o, c_o);
  const Real jplus_i = u_i + Real(2.0) * c_i / (gamma - Real(1.0));
  const Real jplus_o = u_o + Real(2.0) * c_o / (gamma - Real(1.0));
  const Real entropy_i = p_i / std::pow(rho_i, gamma);
  const Real entropy_o = p_o / std::pow(rho_o, gamma);
  require(relative_error(p_o, back_pressure) < tolerance,
          "characteristic pressure-outlet self-test failed p_back");
  require(relative_error(jplus_i, jplus_o) < tolerance,
          "characteristic pressure-outlet self-test failed J+ preservation");
  require(relative_error(entropy_i, entropy_o) < tolerance,
          "characteristic pressure-outlet self-test failed entropy preservation");
  require(relative_error(v_i, v_o) < tolerance,
          "characteristic pressure-outlet self-test failed tangential velocity preservation");

  Real supersonic[ProbClosures::NCONS] = {Real(0.0)};
  make_state(Real(1.0), Real(2.0) * sound, Real(0.0), Real(1.0), supersonic);
  Real supersonic_out[ProbClosures::NCONS] = {Real(0.0)};
  status = GlobalBC::bc_characteristic_pressure_outlet(
      Real(1.0), Real(0.0), Real(0.0), &cls, Real(0.1),
      supersonic, supersonic_out);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic supersonic-outlet self-test returned a failure status");
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    require(supersonic_out[n] == supersonic[n],
            "characteristic supersonic outlet depends on p_back");
  }

  Real reverse[ProbClosures::NCONS] = {Real(0.0)};
  make_state(Real(1.0), -Real(0.1) * sound, Real(0.0), Real(1.0), reverse);
  status = GlobalBC::bc_characteristic_pressure_outlet(
      Real(1.0), Real(0.0), Real(0.0), &cls, Real(1.0), reverse, outlet);
  require(status == GlobalBC::characteristic_bc_status_t::reverse_flow,
          "characteristic pressure outlet did not fail closed on reverse flow");

  Real farfield[ProbClosures::NCONS] = {Real(0.0)};
  status = GlobalBC::bc_characteristic_farfield(
      Real(-1.0), Real(0.0), Real(0.0), &cls,
      Real(1.0), u, v, Real(0.0), Real(1.0), interior, farfield);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic far-field inlet self-test returned a failure status");
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    require(relative_error(farfield[n], interior[n]) < tolerance,
            "characteristic far-field inlet failed uniform-state recovery");
  }
  status = GlobalBC::bc_characteristic_farfield(
      Real(1.0), Real(0.0), Real(0.0), &cls,
      Real(1.0), u, v, Real(0.0), Real(1.0), interior, farfield);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic far-field outlet self-test returned a failure status");
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    require(relative_error(farfield[n], interior[n]) < tolerance,
            "characteristic far-field outlet failed uniform-state recovery");
  }

  // A roundoff-scale sign on a nominally stagnant normal velocity must not
  // switch the boundary into a contradictory inflow/outflow regime.
  for (const Real signed_dust :
       {-std::numeric_limits<Real>::denorm_min(),
         std::numeric_limits<Real>::denorm_min()}) {
    Real stagnant[ProbClosures::NCONS] = {Real(0.0)};
    make_state(Real(1.0), Real(0.0), signed_dust, Real(1.0), stagnant);
    status = GlobalBC::bc_characteristic_farfield(
        Real(0.0), Real(1.0), Real(0.0), &cls,
        Real(1.0), Real(0.0), Real(0.0), Real(0.0), Real(1.0),
        stagnant, farfield);
    require(status == GlobalBC::characteristic_bc_status_t::success,
            "characteristic far-field rejected a stagnant state with "
            "roundoff-scale normal velocity");
    for (int n = 0; n < ProbClosures::NCONS; ++n) {
      require(relative_error(farfield[n], stagnant[n]) < tolerance,
              "characteristic far-field failed stagnant-state recovery");
    }
  }

  // A weak outgoing pressure wave can reverse the reconstructed boundary
  // velocity even when the adjacent cell still has a small inward velocity.
  // The entropy and tangential modes must follow the reconstructed boundary
  // regime rather than make this regular zero crossing fail closed.
  Real reversing[ProbClosures::NCONS] = {Real(0.0)};
  make_state(Real(0.997), Real(-0.19), Real(1.e-5),
             Real(1.001), reversing);
  status = GlobalBC::bc_characteristic_farfield(
      Real(0.0), Real(-1.0), Real(0.0), &cls,
      Real(1.0), Real(0.0), Real(0.0), Real(0.0), Real(1.0),
      reversing, farfield);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic far-field rejected a regular subsonic "
          "flow-direction crossing");
  Real rho_cross, u_cross, v_cross, p_cross, c_cross;
  decode_ideal_gas(farfield, rho_cross, u_cross, v_cross, p_cross, c_cross);
  require(-v_cross > Real(0.0) && -v_cross < c_cross,
          "characteristic far-field selected an inconsistent boundary regime");

  // A finite outgoing wave may carry the reconstructed ghost state through
  // M_n=1 even while the adjacent interior state remains marginally
  // subsonic.  This is a regular transition to an all-outgoing boundary, not
  // a failure of the characteristic closure.
  Real sonic_transition[ProbClosures::NCONS] = {Real(0.0)};
  make_state(
      Real(4.7284795935263375), Real(0.8895434519546643),
      Real(-0.7466090664326998), Real(1.8841832224188872),
      sonic_transition);
  status = GlobalBC::bc_characteristic_farfield(
      Real(0.0), Real(-1.0), Real(0.0), &cls,
      Real(1.0), Real(0.0), Real(0.0), Real(0.0), Real(1.0),
      sonic_transition, farfield);
  require(status == GlobalBC::characteristic_bc_status_t::success,
          "characteristic far-field rejected a regular sonic-outflow "
          "transition");
  Real rho_sonic, u_sonic, v_sonic, p_sonic, c_sonic;
  decode_ideal_gas(
      farfield, rho_sonic, u_sonic, v_sonic, p_sonic, c_sonic);
  require(-v_sonic >= c_sonic,
          "characteristic far-field sonic-transition oracle did not enter "
          "the all-outgoing regime");

  // GC-NSCBC recurrence must exactly reproduce a linear profile through all
  // three ghost layers required by Cerisse.
  const Real q_boundary = Real(1.25);
  const Real hqn = Real(0.04);
  const Real q_in = q_boundary - hqn;
  Real qg[3] = {Real(0.0), Real(0.0), Real(0.0)};
  for (int layer = 1; layer <= 3; ++layer) {
    qg[layer - 1] = nscbc::gc_nscbc_value(
      layer, q_boundary, q_in, hqn, qg[0], qg[1]);
    require(relative_error(qg[layer - 1],
                           q_boundary + Real(layer) * hqn) < tolerance,
            "GC-NSCBC recurrence failed linear-profile reproduction");
  }

  // The same three-layer recurrence is exact for a quadratic normal profile
  // when supplied with the exact boundary value, first interior value, and
  // boundary-normal derivative.  This is the polynomial consistency needed
  // by the second-order stage-local boundary closure.
  const Real quadratic_curvature = Real(0.007);
  const Real q_in_quadratic =
    q_boundary - hqn + quadratic_curvature;
  Real qg_quadratic[3] = {Real(0.0), Real(0.0), Real(0.0)};
  for (int layer = 1; layer <= 3; ++layer) {
    qg_quadratic[layer - 1] = nscbc::gc_nscbc_value(
      layer, q_boundary, q_in_quadratic, hqn,
      qg_quadratic[0], qg_quadratic[1]);
    const Real layer_r = Real(layer);
    const Real exact = q_boundary + layer_r * hqn +
      layer_r * layer_r * quadratic_curvature;
    require(relative_error(qg_quadratic[layer - 1], exact) < tolerance,
            "GC-NSCBC recurrence failed quadratic-profile reproduction");
  }

  // For a divergence-free convected vorticity mode, Giles Eq. (174) must
  // recover the interior pressure and normal-velocity derivatives exactly.
  Real rho_n, un_n, ut_n, p_n;
  nscbc::giles_second_order_normal_derivatives(
    Real(1.0), Real(2.0), Real(0.5),
    Real(0.11), Real(0.3), Real(0.2), Real(0.0), Real(-0.3),
    Real(0.0), rho_n, un_n, ut_n, p_n);
  require(relative_error(rho_n, Real(0.11)) < tolerance &&
          relative_error(un_n, Real(0.3)) < tolerance &&
          relative_error(ut_n, Real(0.2)) < tolerance &&
          relative_error(p_n, Real(0.0)) < tolerance,
          "Giles second-order outflow failed the vorticity-mode oracle");

  amrex::Print() << "[Stage3-Characteristic-BC] algebra self-test PASS\n";
#endif
}

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
void require_characteristic_status(
    const GlobalBC::characteristic_bc_status_t status,
    Real* state)
{
  if (status == GlobalBC::characteristic_bc_status_t::success) return;
  const Real invalid = std::numeric_limits<Real>::quiet_NaN();
  for (int n = 0; n < ProbClosures::NCONS; ++n) state[n] = invalid;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real gaussian_train(const Real x, const Real y, const Real lx, const Real ly,
                    const ProbParm& pparm)
{
  const Real ky = Real(2.0) * Real(3.14159265358979323846264338327950288) *
                  Real(pparm.carrier_y) / ly;
  const Real tan_raw = std::tan(pparm.theta);
  const Real tan_theta = amrex::Math::abs(tan_raw) > Real(1.e-12)
      ? tan_raw
      : ((tan_raw < Real(0.0)) ? -Real(1.e-12) : Real(1.e-12));
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
  } else if (pparm.mode == 9) {
    const Real rx = (x - pparm.x0) / pparm.sigma_x;
    envelope = std::exp(-Real(0.5) * rx * rx);
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
  } else if (pparm.mode == 9) {
    dv = pparm.amp * envelope;
  } else if (pparm.mode == 10) {
    // Construct the velocity from a streamfunction so the initial
    // perturbation is exactly divergence-free.  The older mode 4 omitted the
    // derivative of the x-envelope and therefore contained an acoustic
    // component that invalidated a pure vorticity-to-acoustic conversion test.
    const Real ky = Real(2.0) *
        Real(3.14159265358979323846264338327950288) *
        Real(pparm.carrier_y) / ly;
    const Real tan_raw = std::tan(pparm.theta);
    const Real tan_theta = amrex::Math::abs(tan_raw) > Real(1.e-12)
        ? tan_raw
        : ((tan_raw < Real(0.0)) ? -Real(1.e-12) : Real(1.e-12));
    const Real kx = ky / tan_theta;
    const Real wave_number = std::sqrt(kx * kx + ky * ky);
    const Real dx0 = x - pparm.x0;
    const Real gaussian =
        std::exp(-Real(0.5) * dx0 * dx0 /
                 (pparm.sigma_x * pparm.sigma_x));
    const Real phase = kx * dx0 + ky * (y - pparm.y0);
    const Real psi_scale = pparm.amp / wave_number;
    du = -psi_scale * ky * gaussian * std::sin(phase);
    dv = psi_scale * gaussian *
        (dx0 * std::cos(phase) /
             (pparm.sigma_x * pparm.sigma_x) +
         kx * std::sin(phase));
  }

  const Real rho = pparm.rho0 + drho;
  const Real p = pparm.p0 + dp;
  const Real u = pparm.u0 + du;
  const Real v = pparm.v0 + dv;

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
  if (pparm.all_farfield != 0) {
    Real nx = Real(0.0);
    Real ny = Real(0.0);
    Real nz = Real(0.0);
    const Real outward = (sgn == 1) ? Real(-1.0) : Real(1.0);
    if (idir == 0) nx = outward;
#if (AMREX_SPACEDIM >= 2)
    if (idir == 1) ny = outward;
#endif
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) nz = outward;
#endif
    const auto status = GlobalBC::bc_characteristic_farfield(
        nx, ny, nz, &closures, pparm.rho0, pparm.u0, pparm.v0,
        Real(0.0), pparm.p0, s_int, s_ext);
    require_characteristic_status(status, s_ext);
    return;
  }

  if (idir == 0 && sgn == 1) {
    if (pparm.inlet_model == 1) {
      const Real static_temperature =
          pparm.p0 / (pparm.rho0 * closures.Rspec);
      const Real total_ratio = Real(1.0) +
          Real(0.5) * (pparm.gamma - Real(1.0)) *
              pparm.mach * pparm.mach;
      const Real total_temperature = static_temperature * total_ratio;
      const Real total_pressure = pparm.p0 *
          std::pow(total_ratio, pparm.gamma / (pparm.gamma - Real(1.0)));
      const auto status = GlobalBC::bc_characteristic_total_inlet(
          Real(-1.0), Real(0.0), Real(0.0), &closures,
          total_pressure, total_temperature,
          pparm.u0, pparm.v0, Real(0.0), s_int, s_ext);
      require_characteristic_status(status, s_ext);
      return;
    }
    if (pparm.inlet_model == 2) {
      const auto status = GlobalBC::bc_characteristic_supersonic_inlet(
          Real(-1.0), Real(0.0), Real(0.0), &closures,
          pparm.rho0, pparm.u0, pparm.v0, Real(0.0), pparm.p0,
          s_int, s_ext);
      require_characteristic_status(status, s_ext);
      return;
    }
    if (pparm.inlet_model == 3) {
      const auto status = GlobalBC::bc_characteristic_farfield(
          Real(-1.0), Real(0.0), Real(0.0), &closures,
          pparm.rho0, pparm.u0, pparm.v0, Real(0.0), pparm.p0,
          s_int, s_ext);
      require_characteristic_status(status, s_ext);
      return;
    }
    set_state(closures, pparm, pparm.rho0, pparm.u0, pparm.v0, pparm.p0, s_ext);
    return;
  }

  if (idir == 0 && sgn == -1 && pparm.outlet_model == 1) {
    const auto status = GlobalBC::bc_characteristic_pressure_outlet(
        Real(1.0), Real(0.0), Real(0.0), &closures,
        pparm.back_pressure_ratio * pparm.p0, s_int, s_ext);
    require_characteristic_status(status, s_ext);
    return;
  }

  if (idir == 0 && sgn == -1 && pparm.outlet_model == 2) {
    const auto status = GlobalBC::bc_characteristic_farfield(
        Real(1.0), Real(0.0), Real(0.0), &closures,
        pparm.rho0, pparm.u0, pparm.v0, Real(0.0), pparm.p0,
        s_int, s_ext);
    require_characteristic_status(status, s_ext);
    return;
  }

  if (idir == 0 && sgn == -1 &&
      (pparm.bc_mode == 1 || pparm.bc_mode == 2 ||
       pparm.bc_mode == 3 || pparm.bc_mode == 4)) {
    GlobalBC::bc_nscbc_farfield(Real(1.0), Real(0.0), Real(0.0), &closures,
                                pparm.rho0, pparm.u0, pparm.v0, Real(0.0),
                                pparm.p0, s_int, s_ext);
    return;
  }

  set_state(closures, pparm, pparm.rho0, pparm.u0, pparm.v0, pparm.p0, s_ext);
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
    bp.v_inf = pparm.v0;
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
  bp.v_inf = pparm.v0;
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
    const Real v_inf = pparm.v0;
    const Real p_inf = pparm.p0;
    const Real e_inf = p_inf / (pparm.gamma - Real(1.0)) +
      Real(0.5) * rho_inf * (u_inf * u_inf + v_inf * v_inf);
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
      rhs(i, j, k, cls_t::UMY) += -sigma * (rho * v - rho_inf * v_inf);
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
