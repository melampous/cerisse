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

// 2D Axisymmetric (RZ) underexpanded jet — FIXED NPR = 10.0 standalone case.
// Base: Panda v2 flanged case (12D x 32D, prescribed sonic exit, tanh lip).
// Dimension mapping: i -> r (radial, [0, r_max]), j -> z (streamwise, [0, z_max]).
// z-low (j=0) boundary, three radial zones:
//   r <  r_jet+3*delta          : sonic-exit jet inflow (tanh lip profile)
//   r in [r_jet+3*delta, 6D]    : FLANGE — adiabatic slip wall (exact mirror ghosts)
//   r >  6D                     : quiescent ambient Dirichlet
//
// SWEEP: NPR(t) staircase with tanh transitions. Cold sonic exit => T_e and
// u_e are NPR-independent; only p_e and rho_e scale linearly with NPR(t).
// Plateaus (8):  1.893(perfectly expanded, p_e=p_amb)  2.394(m119)  2.800
//                3.273(m142)  3.671(Mj=1.5)  4.200  4.900  5.746(m180)
// Timeline: P0 hold [0, 3.3] ms (2.0 develop + 1.3 stats); P1..P7 hold 2.5 ms;
// tanh transitions 0.3 ms wide (tau = 75 us), centers tc_k = 3.45 + 2.8k ms.
// End of schedule: 22.9 ms. Quasi-steady margin: dNPR/dt within a transition
// ~1.6/ms over 0.3 ms << plateau spacing; near-field flushing time ~0.55 ms
// is absorbed by the 1.2 ms per-plateau flush segment before statistics.

using namespace amrex;
using namespace universal_constants;

namespace PROB {

static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant / Mw;
static constexpr Real Cv       = Rgas / (gam - 1.0);
static constexpr Real Cp       = gam * Cv;

// ---- NPR sweep schedule ----------------------------------------------------
static constexpr int  NPR_NPLAT = 8;
static constexpr Real NPR_TC0 = 3.45e-3;   // s, center of first transition
static constexpr Real NPR_DTC = 2.80e-3;   // s, transition-center spacing
static constexpr Real NPR_TAU = 7.5e-5;    // s, tanh width (full ramp ~4*tau = 0.3 ms)
static constexpr Real NPR_CRIT = 1.89293;  // (1.2)^3.5, p0/p_e at sonic exit

// Plateau table lives INSIDE the function: a namespace-scope constexpr array
// is not usable from CUDA device code (only scalars constant-fold), which
// broke the first build (2026-07-20). Host-device so the banner can call it.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real npr_of_t (Real t)
{
  amrex::ignore_unused(t);
  return Real(10.0);   // fixed NPR, no sweep
}

// hold-interval midpoint of plateau k (host helper for the banner)
inline Real npr_plateau_mid_t (int k)
{
  return (k == 0) ? 1.65e-3 : 4.85e-3 + Real(k-1) * NPR_DTC;
}

struct ProbParm
{
  Real p_amb    = 99780.0;  // Pa; paper-exact: rho_inf=1.17 @ T_inf=297.15K (derived)
  Real T_amb    = 297.15;   // K = T0_jet (unheated, paper 24C)
  Real Mach_amb = 0.0;      // quiescent ambient, NO coflow

  // NOTE (sweep): p_jet is NOT used by bcnormal here — the inflow static
  // pressure is p_e(t) = p_amb * npr_of_t(t) / NPR_CRIT. Kept for reference.
  Real p_jet    = 99780.0 * 3.273446 / 1.89293;
  Real T_jet    = 247.625;  // = 297.15/1.2 (sonic exit static T, NPR-independent)
  Real Mach_jet = 1.0;

  Real r_jet     = 0.0127;    // 12.7 mm, D_e = 25.4 mm (Panda nozzle)
  Real delta_jet = 254.0e-6;  // delta/De = 0.01, tanh lip layer

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
  const Real c_e = std::sqrt(gam * Rgas * p.T_jet);
  const Real U_e = p.Mach_jet * c_e;

  amrex::Print() << "\n====== 2D Axisym NPR sweep (staircase+tanh, flanged 12Dx32D) ======\n";
  amrex::Print() << "  D_e = " << 2.0*p.r_jet*1000.0 << " mm   U_e = " << U_e
                 << " m/s (NPR-independent, cold sonic exit)\n";
  amrex::Print() << "  plateaus (NPR):";
  for (int kk = 0; kk < NPR_NPLAT; ++kk) amrex::Print() << "  " << npr_of_t(npr_plateau_mid_t(kk));
  amrex::Print() << "\n  schedule: P0 hold [0,3.3]ms; P1..P7 hold 2.5 ms; tanh ramps 0.3 ms\n";
  amrex::Print() << "  transition centers tc_k = 3.45 + 2.8k ms (k=0..6); end 22.9 ms\n";
  for (int kk = 0; kk < NPR_NPLAT; ++kk) {
    const Real npr = npr_of_t(npr_plateau_mid_t(kk));
    const Real Mj  = std::sqrt(5.0 * (std::pow(npr, 2.0/7.0) - 1.0));
    const Real pe  = p.p_amb * npr / NPR_CRIT;
    amrex::Print() << "    P" << kk << ": NPR = " << npr << "   Mj = " << Mj
                   << "   p_e = " << pe << " Pa   rho_e = " << pe/(Rgas*p.T_jet) << " kg/m3\n";
  }
  amrex::Print() << "===================================================================\n\n";
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
  amrex::ignore_unused(dratio, geomdata);

  for (int n = 0; n < ProbClosures::NCONS; ++n)
    s_ext[n] = s_int[n];

  // Only specify z-low
  if (!(idir == 1 && sgn == 1)) return;

  const Real r = x[0];   // radial coordinate

  if (r < pparm.r_jet) {
    // ---- zone 1: choked inlet, r < R (physical lip). Same adiabatic exit-BL
    // profile as the Panda v2 case (uniform static p over the whole orifice,
    // constant T0, tanh lip in velocity only), with the SWEEP modification:
    //   p_e(t) = p_amb * NPR(t) / NPR_CRIT   (sonic exit static pressure)
    // T_e and u_e are NPR-independent for the cold jet, so only p and rho ramp.
    const Real p_e   = pparm.p_amb * npr_of_t(time) / NPR_CRIT;
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
    const Real rho  = p_e / (Rgas * T);
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
// v2 "balanced" AMR (unchanged from the Panda v2 case): dimensionless
// density-gradient sensor gated by per-level spatial caps, plus forced zones.
// At max_level = 4 (this sweep) only forceA is active; the sensor tracks the
// wave pattern as it moves with NPR, which is exactly what a sweep needs.
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
  // fine (lev>=3) 0.03.
  const Real thresh = (level < 3) ? Real(0.05) : Real(0.03);
  const bool sensor = (grad_rho > thresh);

  // ---- forced zones
  const Real rlo_cor = (Real(0.35) + Real(0.04) * zeta) * D;
  const Real rhi_cor = (Real(0.65) + Real(0.12) * zeta) * D;
  const bool in_corridor = (xr >= rlo_cor) && (xr <= rhi_cor);
  const bool forceA = (level < 4) && (xz <= Real(8.0)*D) && (xr <= Real(1.5)*D);  // L4 jet core
  const bool forceB = (level < 5) && (zeta <= Real(6.0)) && in_corridor;          // L5 shear corridor (inactive at max_level=4)
  const bool forceC = (level < 6) && (zeta <= Real(1.0)) && in_corridor;          // L6 lip corridor  (inactive at max_level=4)

  const bool tag = (sensor || forceA || forceB || forceC) && allowed;

  tagfab(i, j, k) = tagfab(i, j, k) || tag;
}

} // namespace PROB

#endif
