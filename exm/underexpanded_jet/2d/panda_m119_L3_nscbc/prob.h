#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>
#include <Constants.h>
#include <bc_types.h>

// 2D Axisymmetric (RZ) underexpanded jet — Panda & Seasholtz 1999, Mj=1.19.
// NSCBC/LODI variant (L3 pilot, 2026-07-15): unlock lip receptivity + open faces.
// Dimension mapping: i -> r (radial, [0, r_max]), j -> z (streamwise, [0, z_max]).
//
// z-low (j=0) boundary, FOUR radial zones (prob.use_lodi = 1):
//   M(r) >= lodi_mach_gate       : sonic/supersonic jet core — Dirichlet (all
//                                  characteristics enter; profile imposed)
//   subsonic tanh rim, r < r_jet : CHARACTERISTIC inflow toward the LOCAL
//                                  profile targets (receptivity unlocked: the
//                                  outgoing u-c invariant is read from the
//                                  interior instead of being overwritten)
//   r in [r_jet, 6D]             : FLANGE — adiabatic slip wall (mirror ghosts)
//   r >  6D                      : ambient via characteristic closure
// r-high / z-high (was foextrap): subsonic-outflow LODI, ambient targets,
//   pressure relax sigma_out (1D-validated R ~ 1e-2; foextrap measured ~3e-2).
// prob.use_lodi = 0 reproduces the legacy hard-Dirichlet/foextrap case exactly
// (baseline twin for A/B receptivity tests).
// Domain 12D x 32D; flange radius 6D = Panda's 305 mm plate.

using namespace amrex;
using namespace universal_constants;

namespace PROB {

static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant / Mw;
static constexpr Real Cv       = Rgas / (gam - 1.0);
static constexpr Real Cp       = gam * Cv;

struct ProbParm
{
  Real p_amb    = 99780.0;  // Pa; paper-exact: rho_inf=1.17 @ T_inf=297.15K (derived)
  Real T_amb    = 297.15;   // K = T0_jet (unheated, paper 24C)
  Real Mach_amb = 0.0;      // Panda: quiescent ambient, NO coflow

  Real p_jet    = 99780.0 * 2.393613 / 1.89293;  // NPR=2.3936 (Mj=1.19) on exp-day p_amb
  Real T_jet    = 247.625;  // = 297.15/1.2 (sonic exit static T)
  Real Mach_jet = 1.0;

  Real r_jet     = 0.0127;    // 12.7 mm, D_e = 25.4 mm (Panda nozzle)
  Real delta_jet = 254.0e-6;  // delta/De = 0.01, tanh lip layer

  Real r_flange  = 0.1524;    // 6D — Panda 305 mm flange radius (adiabatic slip wall)

  // ---- NSCBC/LODI controls (runtime, prob.*) ----
  int  use_lodi       = 1;     // 0 = legacy Dirichlet/foextrap everywhere (baseline twin)
  Real sigma_in       = 0.25;  // lip-annulus inflow pressure relax (bp.pressure_relax)
  Real sigma_out      = 0.25;  // open-face outflow pressure relax
  Real lodi_mach_gate = 0.99;  // zone-1 cells with local M >= gate stay Dirichlet

  ProbParm() {
    amrex::ParmParse pp("prob");
    pp.query("use_lodi", use_lodi);
    pp.query("sigma_in", sigma_in);
    pp.query("sigma_out", sigma_out);
    pp.query("lodi_mach_gate", lodi_mach_gate);
  }
};

struct methodparm_t {
  static constexpr int  order = 2;
  static constexpr bool use_LES = false;
  static constexpr Real conductivity = 0.0262;
  static constexpr Real viscosity   = 1.85e-5;
};

using ProbClosures =
    closures_dt<indicies_stat_t, visc_suth_t, cond_suth_t, calorifically_perfect_gas_t<indicies_t>>;

typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t> ProbRHS;

using GlobalBC = manual_bc_t<ProbClosures>;


void inline inputs() {
  ProbParm p;
  const Real rho_amb  = p.p_amb / (Rgas * p.T_amb);
  const Real eint_amb = Cv * p.T_amb;
  const Real rho_e = p.p_jet / (Rgas * p.T_jet);
  const Real c_e   = std::sqrt(gam * Rgas * p.T_jet);
  const Real U_e   = p.Mach_jet * c_e;
  const Real NPR0  = (p.p_jet / p.p_amb) * std::pow(1.0 + 0.5*(gam-1.0)*p.Mach_jet*p.Mach_jet, gam/(gam-1.0));

  // Reynolds numbers — TWO conventions, do not mix:
  //   Re_e: sonic-exit state (rho_e, U_e). Re_j: fully-expanded state — the
  //   convention Panda reports (0.93e6 / 1.26e6).
  auto mu_suth = [](Real T) { return 1.458e-6 * T * std::sqrt(T) / (T + 110.4); };
  const Real Mj   = std::sqrt(5.0 * (std::pow(NPR0, 2.0/7.0) - 1.0));
  const Real Tj   = p.T_amb / (1.0 + 0.2*Mj*Mj);            // T0 = T_amb (unheated)
  const Real Uj   = Mj * std::sqrt(gam * Rgas * Tj);
  const Real rhoj = p.p_amb / (Rgas * Tj);
  const Real Re_e = rho_e * U_e * 2.0*p.r_jet / mu_suth(p.T_jet);
  const Real Re_j = rhoj  * Uj  * 2.0*p.r_jet / mu_suth(Tj);

  // ACTUAL inlet fluxes: integrate the tanh profile (centered at R-3delta so
  // u_z(R)~0; the flange inner edge is the physical lip r=R).
  const Real r0 = p.r_jet - 3.0*p.delta_jet;
  const int  N  = 20000;
  Real mdot = 0.0, Jz = 0.0, Hdot = 0.0;
  for (int i = 0; i < N; ++i) {
    const Real r   = (Real(i)+0.5) * p.r_jet / N;
    const Real f   = 0.5*(1.0 - std::tanh((r - r0)/p.delta_jet));
    const Real uz  = U_e * f;
    const Real T   = p.T_amb - 0.5*uz*uz/Cp;   // constant T0, adiabatic BL
    const Real rho = p.p_jet / (Rgas * T);
    const Real pr  = p.p_jet;                  // uniform static p = p_e
    const Real dA  = 2.0*M_PI*r * (p.r_jet/N);
    mdot += rho*uz*dA;
    Jz   += (rho*uz*uz + pr)*dA;
    Hdot += rho*uz*(Cp*T + 0.5*uz*uz)*dA;
  }
  const Real mdot_id = rho_e*U_e*M_PI*p.r_jet*p.r_jet;
  const Real x_s1    = 0.67*std::sqrt(NPR0)*2.0*p.r_jet;  // Crist scaling

  amrex::Print() << "\n====== 2D Axisym Panda jet m119 v2 (flanged, 12Dx32D) ======\n";
  amrex::Print() << "  D_e = " << 2.0*p.r_jet*1000.0 << " mm   p_e/p_a = " << p.p_jet/p.p_amb
                 << "   NPR_0 = " << NPR0 << "   Mj = " << Mj << "\n";
  amrex::Print() << "  U_e = " << U_e << " m/s   rho_e = " << rho_e << " kg/m3\n";
  amrex::Print() << "  Re_e = " << Re_e << "   Re_j(fully-expanded, Panda convention) = " << Re_j << "\n";
  amrex::Print() << "  Inlet fluxes (adiabatic BL: p=p_e uniform, T0 const, u_z(R)~0):  mdot = " << mdot
                 << " kg/s  (ideal " << mdot_id << ", Cd_eff = " << mdot/mdot_id << ")\n";
  amrex::Print() << "    Jz = " << Jz << " N   Hdot0 = " << Hdot << " W\n";
  amrex::Print() << "  Flange: adiabatic slip wall,  r in [" << p.r_jet*1000.0
                 << ", " << p.r_flange*1000.0 << "] mm on z-low (inner edge = physical lip)\n";
  amrex::Print() << "  x_shock1_est(Crist) = " << x_s1/(2.0*p.r_jet)
                 << " D — GRID-EXTENT ESTIMATE ONLY (Mj<1.5: conical closure, no Mach disk; not a validation quantity)\n";
  amrex::Print() << "=======================================================\n\n";
}


////////////////////////////// Initial conditions //////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, Array4<Real> const& state,
                    GeometryData const& geomdata, ProbClosures const& cls,
                    ProbParm const& pparm)
{
  amrex::ignore_unused(i, j, k, geomdata);

  // Quiescent ambient
  const Real rho  = pparm.p_amb / (Rgas * pparm.T_amb);
  const Real c    = std::sqrt(gam * Rgas * pparm.T_amb);
  const Real uz   = pparm.Mach_amb * c;              // axial (z) velocity
  const Real eint = Cv * pparm.T_amb;
  const Real rhoE = rho * (eint + 0.5_rt * uz * uz);

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = 0.0_rt;                  // radial (r) momentum
  state(i, j, k, cls.UMY)  = rho * uz;                // axial  (z) momentum
  state(i, j, k, cls.UMZ)  = 0.0_rt;                  // theta momentum (identically 0 in RZ no-swirl)
  state(i, j, k, cls.UET)  = rhoE;
}


////////////////////////////// Boundary conditions /////////////////////////////
// z-low (idir=1, sgn=+1), three radial zones (see header). Other faces:
// r=0 axis symmetry (type 3), r-high/z-high foextrap (type 2).
// Slip wall uses per-row EXACT mirror ghosts (s_refl): u_z odd, u_r even,
// rho/E even -> u_n=0, zero shear, dT/dn=0 (adiabatic), pressure reflection R=+1.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio,
         const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS],
         amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& geomdata, ProbClosures const& cls,
         ProbParm const& pparm)
{
  amrex::ignore_unused(dratio, time, geomdata);

  for (int n = 0; n < ProbClosures::NCONS; ++n)
    s_ext[n] = s_int[n];

  // Only specify z-low
  if (!(idir == 1 && sgn == 1)) return;

  const Real r = x[0];   // radial coordinate

  if (r < pparm.r_jet) {
    // ---- zone 1: choked inlet, r < R (physical lip). The tanh lip profile
    // sits fully INSIDE the orifice, centered at R-3*delta, so u_z(R)~0 and
    // the flange inner edge is the physical lip r=R (acoustic aperture = D).
    // Cost: mdot deficit ~ (1-3delta/R)^2 vs ideal — an effective discharge
    // coefficient, printed by the banner integral.
    // Adiabatic exit-BL profile (review 2026-07-15): STATIC PRESSURE is
    // UNIFORM p = p_e over the WHOLE orifice — the underexpansion (pressure)
    // aperture is the full D. Only the velocity carries the lip layer; T
    // follows constant total temperature T0 = T_amb. (The previous full-blend
    // ramp of rho/T/p shrank the effective pressure aperture to ~0.94D and
    // weakened the first expansion fan.)
    const Real c_jet = std::sqrt(gam * Rgas * pparm.T_jet);
    const Real uz_e  = pparm.Mach_jet * c_jet;
    const Real r0 = pparm.r_jet - Real(3.0) * pparm.delta_jet;
    Real fjet;
    if (pparm.delta_jet > Real(0.0)) {
      fjet = Real(0.5) * (Real(1.0) - std::tanh((r - r0) / pparm.delta_jet));
    } else {
      fjet = (r <= r0) ? Real(1.0) : Real(0.0);
    }
    const Real uz   = uz_e * fjet;
    const Real T    = pparm.T_amb - Real(0.5) * uz * uz / Cp;  // T0 = T_amb
    const Real rho  = pparm.p_jet / (Rgas * T);
    const Real eint = Cv * T;

    s_ext[cls.URHO] = rho;
    s_ext[cls.UMX]  = Real(0.0);
    s_ext[cls.UMY]  = rho * uz;
    s_ext[cls.UET]  = rho * (eint + Real(0.5) * uz * uz);
  } else if (r <= pparm.r_flange) {
    // ---- zone 2: FLANGE — adiabatic slip wall (exact mirror of the reflected cell)
    s_ext[cls.URHO] =  s_refl[cls.URHO];
    s_ext[cls.UMX]  =  s_refl[cls.UMX];   // tangential (radial) momentum: even
    s_ext[cls.UMY]  = -s_refl[cls.UMY];   // wall-normal (axial) momentum: odd -> u_n = 0
    s_ext[cls.UET]  =  s_refl[cls.UET];   // |u| unchanged -> same KE; T even -> adiabatic
  } else {
    // ---- zone 3: quiescent ambient Dirichlet beyond the flange edge
    const Real rho_amb  = pparm.p_amb / (Rgas * pparm.T_amb);
    const Real eint_amb = Cv * pparm.T_amb;
    s_ext[cls.URHO] = rho_amb;
    s_ext[cls.UMX]  = Real(0.0);
    s_ext[cls.UMY]  = Real(0.0);
    s_ext[cls.UET]  = rho_amb * eint_amb;
  }

  s_ext[cls.UMZ] = Real(0.0);  // theta momentum: identically 0 (RZ, no swirl)
}


///////////////////////// NSCBC/LODI boundary hook //////////////////////////
// Called by bcs.cpp BEFORE bcnormal on every ext_dir ghost point (SFINAE
// hook). Return true = handled here; false = fall through to bcnormal above.
// Zones handled here:
//   z-low subsonic lip annulus  -> characteristic inflow, per-point profile
//                                  targets (receptivity fix: kernel keeps the
//                                  outgoing interior invariant)
//   z-low ambient beyond flange -> characteristic closure, ambient targets
//   r-high / z-high open faces  -> subsonic-outflow LODI, ambient targets
// Left to bcnormal: sonic/supersonic jet core (Dirichlet), flange mirror wall.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
bcnormal_lodi(const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& state,
              const amrex::Real* x, amrex::Real /*dratio*/,
              const amrex::Real* /*s_int*/, const amrex::Real* /*s_refl*/,
              amrex::Real* s_ext,
              const int idir, const int sgn, const amrex::Real /*time*/,
              amrex::GeometryData const& geomdata, ProbClosures const& cls,
              ProbParm const& pparm)
{
#if (AMREX_SPACEDIM < 2)
  return false;
#else
  if (!pparm.use_lodi) return false;

  const Real rho_amb = pparm.p_amb / (Rgas * pparm.T_amb);

  // ---- open faces: r-high (idir 0, sgn -1) and z-high (idir 1, sgn -1).
  // Subsonic outflow takes the full LODI path (pressure-relaxed, transverse
  // terms); locally supersonic outflow (jet core crossing z-high) and any
  // transient inflow are delegated by the kernel to the algebraic
  // characteristic far-field closure.
  if (sgn == -1 && (idir == 0 || idir == 1)) {
    GlobalBC::nscbc_lodi_bc_parm_t bp;
    bp.rho_inf = rho_amb;
    bp.u_inf = Real(0.0);
    bp.v_inf = Real(0.0);
    bp.w_inf = Real(0.0);
    bp.p_inf = pparm.p_amb;
    bp.pressure_relax = pparm.sigma_out;
    // bp.l_ref <= 0 -> kernel uses the domain extent along idir
    return GlobalBC::bc_nscbc_lodi_farfield(iv, state, s_ext, idir, sgn,
                                            geomdata, &cls, bp);
  }

  if (idir == 1 && sgn == 1) {  // z-low plane
    const Real r = x[0];

    if (r < pparm.r_jet) {
      // Local adiabatic-BL exit profile (identical math to bcnormal zone 1).
      const Real c_jet = std::sqrt(gam * Rgas * pparm.T_jet);
      const Real uz_e  = pparm.Mach_jet * c_jet;
      const Real r0 = pparm.r_jet - Real(3.0) * pparm.delta_jet;
      Real fjet;
      if (pparm.delta_jet > Real(0.0)) {
        fjet = Real(0.5) * (Real(1.0) - std::tanh((r - r0) / pparm.delta_jet));
      } else {
        fjet = (r <= r0) ? Real(1.0) : Real(0.0);
      }
      const Real uz    = uz_e * fjet;
      const Real T     = pparm.T_amb - Real(0.5) * uz * uz / Cp;  // T0 = T_amb
      const Real c_loc = std::sqrt(gam * Rgas * T);

      // Sonic/supersonic core: every characteristic enters the domain, the
      // profile Dirichlet in bcnormal is the correct closure. Keep it.
      if (uz >= pparm.lodi_mach_gate * c_loc) return false;

      // Subsonic lip annulus: characteristic inflow relaxing toward the LOCAL
      // profile state. The kernel reads the interior state, so the outgoing
      // (u-c) invariant survives — this is the receptivity unlock.
      GlobalBC::nscbc_lodi_bc_parm_t bp;
      bp.rho_inf = pparm.p_jet / (Rgas * T);
      bp.u_inf = Real(0.0);   // radial target
      bp.v_inf = uz;          // axial target = local profile
      bp.w_inf = Real(0.0);
      bp.p_inf = pparm.p_jet;
      bp.pressure_relax = pparm.sigma_in;
      bp.l_ref = Real(2.0) * pparm.r_jet;  // relaxation length = D
      return GlobalBC::bc_nscbc_lodi_farfield(iv, state, s_ext, idir, sgn,
                                              geomdata, &cls, bp);
    }

    if (r <= pparm.r_flange) return false;  // flange mirror wall in bcnormal

    // Ambient segment beyond the flange edge: characteristic closure instead
    // of a hard Dirichlet (quiescent targets; un~0 -> algebraic branch).
    GlobalBC::nscbc_lodi_bc_parm_t bp;
    bp.rho_inf = rho_amb;
    bp.u_inf = Real(0.0);
    bp.v_inf = Real(0.0);
    bp.w_inf = Real(0.0);
    bp.p_inf = pparm.p_amb;
    bp.pressure_relax = pparm.sigma_out;
    return GlobalBC::bc_nscbc_lodi_farfield(iv, state, s_ext, idir, sgn,
                                            geomdata, &cls, bp);
  }

  return false;  // r-low axis is reflect-type, never reaches ext_dir fill
#endif
}


////////////////////////////////// tagging /////////////////////////////////////
// v2 "balanced" AMR (2026-07-14):
//   - dimensionless density-gradient sensor kept (0.03 at lev0, 0.05 above),
//     gated by per-level spatial caps (AND of r,z) so far field stays coarse;
//   - three forced zones guarantee resolution independent of the sensor:
//       A: jet core      L4  (z<=8D,  r<=1.5D)
//       B: shear corridor L5 (z<=6D,  0.35+0.04*zeta <= r/D <= 0.65+0.12*zeta)
//       C: lip corridor   L6 (z<=1D,  same corridor)   zeta = z/D
template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level)
{
  if (nt <= 0) return;

  using idx = indicies_t;
  const Real rho_c = sdatafab(i, j, k, idx::URHO);
  const Real drhor = amrex::Math::abs(sdatafab(i+1,j,k,idx::URHO) - sdatafab(i-1,j,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real drhoz = amrex::Math::abs(sdatafab(i,j+1,k,idx::URHO) - sdatafab(i,j-1,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real grad_rho = std::sqrt(drhor*drhor + drhoz*drhoz);

  const Real xr = geomdata.ProbLo()[0] + (Real(i)+Real(0.5))*geomdata.CellSize()[0];
  const Real xz = geomdata.ProbLo()[1] + (Real(j)+Real(0.5))*geomdata.CellSize()[1];

  const Real D    = Real(2.0) * pparm.r_jet;   // 25.4 mm
  const Real zeta = xz / D;

  // ---- per-level spatial caps: tagging at lev l creates level l+1 (AND, not OR)
  //                 lev:      0        1        2        3        4        5
  const Real rcap[6] = { Real(0.1778), Real(0.1270), Real(0.1016), Real(0.0762), Real(0.0508), Real(0.03175) }; // 7D 5D 4D 3D 2D 1.25D
  const Real zcap[6] = { Real(0.6096), Real(0.5080), Real(0.4064), Real(0.3048), Real(0.2032), Real(0.1016)  }; // 24D 20D 16D 12D 8D 4D
  const int  lc = (level < 6) ? level : 5;
  const bool allowed = (xr <= rcap[lc]) && (xz <= zcap[lc]);

  // ---- dynamic sensor (dimensionless). Sensor value scales ~dx for smooth
  // gradients, so FINE levels get the LOWER threshold: coarse (lev<3) 0.05,
  // fine (lev>=3) 0.03 — lets shock-vortex structures cascade to L5/L6
  // inside the caps instead of self-limiting.
  const Real thresh = (level < 3) ? Real(0.05) : Real(0.03);
  const bool sensor = (grad_rho > thresh);

  // ---- forced zones
  const Real rlo_cor = (Real(0.35) + Real(0.04) * zeta) * D;
  const Real rhi_cor = (Real(0.65) + Real(0.12) * zeta) * D;
  const bool in_corridor = (xr >= rlo_cor) && (xr <= rhi_cor);
  const bool forceA = (level < 4) && (xz <= Real(8.0)*D) && (xr <= Real(1.5)*D);  // L4 jet core
  const bool forceB = (level < 5) && (zeta <= Real(6.0)) && in_corridor;          // L5 shear corridor
  const bool forceC = (level < 6) && (zeta <= Real(1.0)) && in_corridor;          // L6 lip corridor

  // caps gate EVERYTHING (forced zones lie inside the caps by construction,
  // but the AND makes that a hard invariant, not a coincidence)
  const bool tag = (sensor || forceA || forceB || forceC) && allowed;

  tagfab(i, j, k) = tagfab(i, j, k) || tag;
}

} // namespace PROB

#endif
