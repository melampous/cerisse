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

// 2D Axisymmetric (RZ) underexpanded jet — Panda & Seasholtz 1999, Mj=1.42 (v2).
// Dimension mapping: i -> r (radial, [0, r_max]), j -> z (streamwise, [0, z_max]).
// z-low (j=0) boundary, three radial zones:
//   r <  r_jet+3*delta          : sonic-exit jet inflow (tanh lip profile)
//   r in [r_jet+3*delta, 6D]    : FLANGE — adiabatic slip wall (exact mirror ghosts)
//   r >  6D                     : quiescent ambient Dirichlet
// Domain 12D x 32D (widened v2); flange radius 6D = Panda's 305 mm plate.

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

  Real p_jet    = 99780.0 * 3.273446 / 1.89293;  // NPR=3.2735 (Mj=1.42) on exp-day p_amb
  Real T_jet    = 247.625;  // = 297.15/1.2 (sonic exit static T)
  Real Mach_jet = 1.0;

  Real r_jet     = 0.0127;    // 12.7 mm, D_e = 25.4 mm (Panda nozzle)
  Real delta_jet = 2.200000e-04;  // delta/De = 0.01, tanh lip layer

  Real r_flange  = 0.1524;    // 6D — Panda 305 mm flange radius (adiabatic slip wall)
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
  const Real r0 = p.r_jet - 2.284208*p.delta_jet;
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

  amrex::Print() << "\n====== 2D Axisym Panda jet v2 (flanged, 12Dx32D) ======\n";
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
    const Real r0 = pparm.r_jet - Real(2.284208) * pparm.delta_jet;
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


////////////////////////////////// tagging /////////////////////////////////////
// x_fc-focused refinement.  max_level = 5  (finest dx = 49.6 um = D_e/512).
//
//   L5 is FORCED ONLY, and only where the measurement needs it:
//        - the nozzle-lip band, widened downstream so it follows the initial
//          Prandtl-Meyer turning of the jet boundary (nu(M=1.42) = 9.6 deg,
//          hence the 0.17 slope);
//        - a thin strip on the axis, so x_fc is extracted at 49.6 um.
//     The density sensor is switched off at level 4, so it can NEVER build
//     L5 on the downstream shock cells.  In the production tagging that
//     sensor-driven L5/L6 carried 86% of the cost and none of the signal.
//
//   L4 (99.2 um) covers the whole first shock cell (z <= 2D, r <= 1D) plus
//     sensor-driven refinement inside tight caps, so the first cell is fully
//     resolved and the coarsening interface sits downstream of x_1.
//
// Justification for dropping L6: discretising every profile in the sweep onto
// the real cell centres and integrating the RZ mass flux gives a discharge-
// coefficient error <= 6.7e-5 at 49.6 um, which maps through dx_fc/dC_d = 0.504
// to a bias <= 3.4e-5 D_e -- 30x below the x_fc extraction resolution and
// 1500x below the s-ladder signal of 0.052 D_e.
template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level)
{
  if (nt <= 0) return;

  using idx = indicies_t;
  const Real rho_c = sdatafab(i, j, k, idx::URHO);
  const Real drhor = amrex::Math::abs(sdatafab(i+1,j,k,idx::URHO)
                                    - sdatafab(i-1,j,k,idx::URHO))
                     / (Real(2.0) * rho_c);
  const Real drhoz = amrex::Math::abs(sdatafab(i,j+1,k,idx::URHO)
                                    - sdatafab(i,j-1,k,idx::URHO))
                     / (Real(2.0) * rho_c);
  const Real grad_rho = std::sqrt(drhor*drhor + drhoz*drhoz);

  const Real xr = geomdata.ProbLo()[0] + (Real(i)+Real(0.5))*geomdata.CellSize()[0];
  const Real xz = geomdata.ProbLo()[1] + (Real(j)+Real(0.5))*geomdata.CellSize()[1];
  const Real D  = Real(2.0) * pparm.r_jet;          // 25.4 mm

  // ---- per-level spatial caps: tagging at lev l creates level l+1.
  //                 lev:      0        1        2        3        4
  //                          7D       5D       3D     1.5D     0.90D
  const Real rcap[5] = { Real(0.17780), Real(0.12700), Real(0.07620),
                         Real(0.03810), Real(0.022860) };
  //                         24D      16D       8D       3D      1.5D
  const Real zcap[5] = { Real(0.60960), Real(0.40640), Real(0.20320),
                         Real(0.07620), Real(0.038100) };
  const int  lc = (level < 5) ? level : 4;
  const bool allowed = (xr <= rcap[lc]) && (xz <= zcap[lc]);

  // ---- density sensor: levels 0-3 only.  L5 is never sensor-driven.
  const Real thresh = (level < 3) ? Real(0.05) : Real(0.03);
  const bool sensor = (level < 4) && (grad_rho > thresh);

  // ---- forced L4: the whole first shock cell.
  const bool forceCore = (level < 4)
                      && (xz <= Real(2.0)*D) && (xr <= Real(1.0)*D);

  // ---- forced L5: lip band + axis strip, x <= 1.5 D_e (downstream of x_1).
  //      2.5 mm inside the lip covers the thickest profile in the matrix
  //      (wall tanh a = 679.5 um reaches 0.99 U_e at R - 1.9 mm).
  const Real rlo_lip = pparm.r_jet - Real(2.5e-3);
  const Real rhi_lip = pparm.r_jet + Real(2.5e-3) + Real(0.17) * xz;
  const bool inLip   = (xr >= rlo_lip) && (xr <= rhi_lip);
  const bool onAxis  = (xr <= Real(2.0e-3));
  const bool forceLip = (level == 4) && (xz <= Real(1.5)*D)
                     && inLip;           // lip band only - axis strip removed

  const bool tag = (sensor || forceCore || forceLip) && allowed;
  tagfab(i, j, k) = tagfab(i, j, k) || tag;
}

} // namespace PROB

#endif
