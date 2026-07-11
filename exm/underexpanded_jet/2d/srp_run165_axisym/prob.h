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

// 2D Axisymmetric (RZ) Supersonic Retropropulsion (SRP) — Run 165
// ------------------------------------------------------------------
//   Reference: Berry et al. Phase-I UPWT SRP test — Run 165 flight
//   conditions (most common CFD validation case).
//   Vehicle flies in +z into M=4.6 freestream; retro jet fires in +z
//   (counter to freestream) from a 4:1 area-ratio conical nozzle.
//   In this simulation the vehicle frame is fixed:
//     - freestream  enters from z-high, flowing in -z  (axial)
//     - jet         enters from z-low  at r < r_jet, flowing in +z
//     - r=0 is the axis of symmetry (forebody centreline)
// ------------------------------------------------------------------
//   Freestream (Run 165):
//     M_inf = 4.6,  T_inf = 65 K,  Re/ft = 1.5e6  (UPWT)
//     With Sutherland-air: mu(65K) ~ 4.36e-6 Pa*s
//       => rho_inf ~ 0.0289 kg/m^3,  p_inf ~ 538 Pa,  V_inf ~ 743 m/s
//
//   Jet stagnation (Run 165):
//     p0_jet / p_inf = 7724,  T0_jet / T_inf = 5.34
//     => p0_jet ~ 4.156 MPa,  T0_jet ~ 347 K
//     Nozzle A_e/A* = 4  =>  M_e ~ 2.94 (isentropic)
//     => p_e ~ 122.4 kPa,  T_e ~ 127 K,  V_e ~ 665 m/s,  rho_e ~ 3.36 kg/m^3
//     mdot = 0.62 lbm/s = 0.281 kg/s  =>  r_jet ~ 6.3 mm (matches hardware)
// ------------------------------------------------------------------

using namespace amrex;
using namespace universal_constants;

namespace PROB {

static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant / Mw;       // ~287.06 J/(kg K)
static constexpr Real Cv       = Rgas / (gam - 1.0);
static constexpr Real Cp       = gam * Cv;

struct ProbParm
{
  //------------------- Freestream (Run 165) -------------------
  Real p_inf    = 538.0;           // Pa (from Re/ft=1.5e6, T=65K)
  Real T_inf    = 65.0;            // K
  Real Mach_inf = 4.6;             // flows in -z (into the vehicle)

  //----------------- Jet exit (after 4:1 nozzle) --------------
  // Isentropic from stagnation (p0,T0) with M_e = 2.94
  Real p_jet    = 122400.0;        // Pa
  Real T_jet    = 127.2;           // K
  Real Mach_jet = 2.94;            // flows in +z (out of vehicle)
  Real r_jet    = 0.0063;          // 6.3 mm  (mdot-consistent)
  Real delta_jet = 0.0;            // sharp top-hat
};

struct methodparm_t {
  static constexpr int  order = 2;
  static constexpr bool use_LES = false;
  static constexpr Real conductivity = 0.0262;
  static constexpr Real viscosity   = 1.85e-5;
};

using ProbClosures =
    closures_dt<indicies_stat_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>, calorifically_perfect_gas_t<indicies_stat_t>>;

typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t, no_source_t> ProbRHS;


void inline inputs() {
  ProbParm p;
  const Real rho_inf = p.p_inf / (Rgas * p.T_inf);
  const Real c_inf   = std::sqrt(gam * Rgas * p.T_inf);
  const Real V_inf   = p.Mach_inf * c_inf;

  const Real rho_jet = p.p_jet / (Rgas * p.T_jet);
  const Real c_jet   = std::sqrt(gam * Rgas * p.T_jet);
  const Real V_jet   = p.Mach_jet * c_jet;

  const Real pe_pinf   = p.p_jet / p.p_inf;
  const Real mdot      = rho_jet * V_jet * M_PI * p.r_jet * p.r_jet;

  // Reference thrust = (rho_e V_e^2 + (p_e - p_inf)) * A_e
  const Real A_e       = M_PI * p.r_jet * p.r_jet;
  const Real Thrust    = rho_jet * V_jet * V_jet * A_e + (p.p_jet - p.p_inf) * A_e;

  amrex::Print() << "\n========== 2D RZ SRP Run-165 ==========\n";
  amrex::Print() << "  FREESTREAM:  M_inf = " << p.Mach_inf
                 << "  T_inf = " << p.T_inf << " K"
                 << "  p_inf = " << p.p_inf << " Pa\n";
  amrex::Print() << "               rho_inf = " << rho_inf
                 << " kg/m^3,  V_inf = " << V_inf << " m/s (-z)\n";
  amrex::Print() << "  JET EXIT:    M_e = " << p.Mach_jet
                 << "  T_e = " << p.T_jet << " K"
                 << "  p_e = " << p.p_jet << " Pa\n";
  amrex::Print() << "               rho_e = " << rho_jet
                 << " kg/m^3,  V_e = " << V_jet << " m/s (+z)\n";
  amrex::Print() << "               r_jet = " << p.r_jet*1000.0 << " mm,  "
                 << "A_e = " << A_e*1e6 << " mm^2\n";
  amrex::Print() << "  RATIOS:      p_e/p_inf = " << pe_pinf << "\n";
  amrex::Print() << "               mdot     = " << mdot << " kg/s\n";
  amrex::Print() << "               Thrust   = " << Thrust << " N\n";
  amrex::Print() << "=======================================\n\n";
}


////////////////////////////// Initial conditions //////////////////////////////
// Initialise the whole domain with the freestream state (supersonic in -z).
// This avoids a long transient establishing the bow shock layer.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, Array4<Real> const& state,
                    GeometryData const& geomdata, ProbClosures const& cls,
                    ProbParm const& pparm)
{
  amrex::ignore_unused(i, j, k, geomdata);

  const Real rho  = pparm.p_inf / (Rgas * pparm.T_inf);
  const Real c    = std::sqrt(gam * Rgas * pparm.T_inf);
  const Real uz   = -pparm.Mach_inf * c;      // freestream flows in -z
  const Real eint = Cv * pparm.T_inf;
  const Real rhoE = rho * (eint + 0.5_rt * uz * uz);

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = 0.0_rt;           // radial momentum
  state(i, j, k, cls.UMY)  = rho * uz;         // axial  momentum (-z)
#if (AMREX_SPACEDIM == 3)
  state(i, j, k, cls.UMZ)  = 0.0_rt;
#endif
  state(i, j, k, cls.UET)  = rhoE;
}


////////////////////////////// Boundary conditions /////////////////////////////
// z-lo (idir=1, sgn=+1):  jet for r < r_jet (supersonic inflow in +z),
//                         outflow (foextrap) elsewhere.
// z-hi (idir=1, sgn=-1):  supersonic freestream inflow in -z everywhere.
// r=0  (x-lo):  symmetry  (handled automatically by BC type 3)
// r-hi (x-hi):  outflow   (handled automatically by BC type 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio,
         const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS],
         amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& geomdata, ProbClosures const& cls,
         ProbParm const& pparm)
{
  amrex::ignore_unused(dratio, s_refl, time, geomdata);

  // Default: foextrap / outflow (copy interior state outward)
  for (int n = 0; n < ProbClosures::NCONS; ++n)
    s_ext[n] = s_int[n];

  // We only specialise the two axial faces (idir=1)
  if (idir != 1) return;

  // Pre-compute freestream and jet primitive states (cheap, inlined)
  const Real rho_inf  = pparm.p_inf / (Rgas * pparm.T_inf);
  const Real c_inf    = std::sqrt(gam * Rgas * pparm.T_inf);
  const Real uz_inf   = -pparm.Mach_inf * c_inf;             // -z direction
  const Real eint_inf = Cv * pparm.T_inf;

  const Real rho_jet  = pparm.p_jet / (Rgas * pparm.T_jet);
  const Real c_jet    = std::sqrt(gam * Rgas * pparm.T_jet);
  const Real uz_jet   = pparm.Mach_jet * c_jet;              // +z direction
  const Real eint_jet = Cv * pparm.T_jet;

  if (sgn == 1) {
    // ---------- z-lo face : jet + outflow ----------
    const Real r = x[0];

    Real fjet;
    if (pparm.delta_jet > Real(0.0)) {
      fjet = Real(0.5) * (Real(1.0) - std::tanh((r - pparm.r_jet) / pparm.delta_jet));
    } else {
      fjet = (r <= pparm.r_jet) ? Real(1.0) : Real(0.0);
    }

    if (fjet > Real(0.0)) {
      // Jet inflow state (blend with interior through fjet when delta_jet>0)
      const Real rho  = fjet * rho_jet  + (Real(1.0) - fjet) * s_int[cls.URHO];
      const Real uz   = fjet * uz_jet;                        // zero outside jet core
      const Real eint = fjet * eint_jet + (Real(1.0) - fjet) * eint_inf;
      const Real ke   = Real(0.5) * uz * uz;

      s_ext[cls.URHO] = rho;
      s_ext[cls.UMX]  = Real(0.0);
      s_ext[cls.UMY]  = rho * uz;
      s_ext[cls.UET]  = rho * (eint + ke);
    }
    // else: keep foextrap default (outflow for plume wake outside jet core)

  } else {
    // ---------- z-hi face : supersonic freestream inflow ----------
    const Real ke = Real(0.5) * uz_inf * uz_inf;
    s_ext[cls.URHO] = rho_inf;
    s_ext[cls.UMX]  = Real(0.0);
    s_ext[cls.UMY]  = rho_inf * uz_inf;
    s_ext[cls.UET]  = rho_inf * (eint_inf + ke);
  }
}


////////////////////////////////// tagging /////////////////////////////////////
// Refine on density gradient: picks up bow shock, jet plume, Mach disk, and
// shear layers. SRP has strong gradients so thresholds are moderate.
template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level)
{
  amrex::ignore_unused(geomdata, pparm);
  if (nt <= 0) return;

  using idx = indicies_stat_t;
  const Real rho_c = sdatafab(i, j, k, idx::URHO);
  const Real drhor = amrex::Math::abs(sdatafab(i+1,j,k,idx::URHO) - sdatafab(i-1,j,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real drhoz = amrex::Math::abs(sdatafab(i,j+1,k,idx::URHO) - sdatafab(i,j-1,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real grad_rho = std::sqrt(drhor*drhor + drhoz*drhoz);

  if (level == 0) {
    tagfab(i, j, k) = tagfab(i, j, k) || (grad_rho > Real(0.05));
  } else {
    tagfab(i, j, k) = tagfab(i, j, k) || (grad_rho > Real(0.10));
  }
}

} // namespace PROB

#endif
