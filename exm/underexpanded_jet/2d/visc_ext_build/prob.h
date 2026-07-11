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

// 2D Axisymmetric (RZ) underexpanded jet.
// Dimension mapping: i → r (radial, [0, r_max]), j → z (streamwise, [0, z_max])
// Nozzle on z-low (j=0) boundary for r < r_jet.
// NPR_0 = 3.0  (p_0/p_a = 3)  =>  pe/pa ~ 1.585

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
  Real p_amb    = 101325.0;
  Real T_amb    = 300.0;
  Real Mach_amb = 0.05;          // weak coflow in z

  // NPR_0 = 3.0 => pe/pa = NPR_0 / (1+0.2)^3.5 = 3.0 / 1.89293
  Real p_jet    = 101325.0 * 20 / 1.89293;
  Real T_jet    = 300.0;
  Real Mach_jet = 1.0;

  Real r_jet     = 0.005;         // 5 mm, D_e = 10 mm
  Real delta_jet = 0.0;           // sharp top-hat
};

struct methodparm_t {
  static constexpr int  order = 2;
  static constexpr bool use_LES = false;
  static constexpr Real conductivity = 0.0262;
  static constexpr Real viscosity   = 1.85e-5;
};

using ProbClosures =
    closures_dt<indicies_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>, calorifically_perfect_gas_t<indicies_t>>;

typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t> ProbRHS;


void inline inputs() {
  ProbParm p;
  const Real rho_amb = p.p_amb / (Rgas * p.T_amb);
  const Real c_jet   = std::sqrt(gam * Rgas * p.T_jet);
  const Real V_jet   = p.Mach_jet * c_jet;
  const Real NPR0    = (p.p_jet / p.p_amb) * std::pow(1.0 + 0.5*(gam-1.0)*p.Mach_jet*p.Mach_jet, gam/(gam-1.0));
  const Real x_MD    = 0.67 * std::sqrt(NPR0) * 2.0 * p.r_jet;

  amrex::Print() << "\n========== 2D Axisymmetric Jet (NPR_0=3) ==========\n";
  amrex::Print() << "  D_e   = " << 2.0*p.r_jet*1000.0 << " mm\n";
  amrex::Print() << "  p_e/p_a = " << p.p_jet/p.p_amb << ",  NPR_0 = " << NPR0 << "\n";
  amrex::Print() << "  V_jet = " << V_jet << " m/s\n";
  amrex::Print() << "  Predicted Mach disk x_MD = " << x_MD*1000.0 << " mm  (x_MD/D_e = " << x_MD/(2.0*p.r_jet) << ")\n";
  amrex::Print() << "===================================================\n\n";
}


////////////////////////////// Initial conditions //////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, Array4<Real> const& state,
                    GeometryData const& geomdata, ProbClosures const& cls,
                    ProbParm const& pparm)
{
  amrex::ignore_unused(i, j, k, geomdata);

  // Quiescent ambient (weak coflow in z direction = UMY)
  const Real rho  = pparm.p_amb / (Rgas * pparm.T_amb);
  const Real c    = std::sqrt(gam * Rgas * pparm.T_amb);
  const Real uz   = pparm.Mach_amb * c;              // axial (z) velocity
  const Real eint = Cv * pparm.T_amb;
  const Real rhoE = rho * (eint + 0.5_rt * uz * uz);

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = 0.0_rt;                  // radial (r) momentum
  state(i, j, k, cls.UMY)  = rho * uz;                // axial  (z) momentum
#if (AMREX_SPACEDIM == 3)
  state(i, j, k, cls.UMZ)  = 0.0_rt;
#endif
  state(i, j, k, cls.UET)  = rhoE;
}


////////////////////////////// Boundary conditions /////////////////////////////
// Jet inflow on z-low (idir=1, sgn=+1):
//   - radius r = x[0] (the first coord in 2D RZ)
//   - if r < r_jet: jet state; else: ambient coflow
// All other faces default: external extrapolation / foextrap (set via BC type 2)
// r=0 axis (x-lo): symmetry BC (type 3) — reflect_even scalars + reflect_odd radial velocity.
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

  for (int n = 0; n < ProbClosures::NCONS; ++n)
    s_ext[n] = s_int[n];

  // Only specify z-low inflow
  if (!(idir == 1 && sgn == 1)) return;

  const Real r = x[0];   // radial coordinate

  const Real rho_amb  = pparm.p_amb / (Rgas * pparm.T_amb);
  const Real c_amb    = std::sqrt(gam * Rgas * pparm.T_amb);
  const Real uz_amb   = pparm.Mach_amb * c_amb;
  const Real eint_amb = Cv * pparm.T_amb;

  const Real rho_jet  = pparm.p_jet / (Rgas * pparm.T_jet);
  const Real c_jet    = std::sqrt(gam * Rgas * pparm.T_jet);
  const Real uz_jet   = pparm.Mach_jet * c_jet;
  const Real eint_jet = Cv * pparm.T_jet;

  Real fjet;
  if (pparm.delta_jet > Real(0.0)) {
    fjet = Real(0.5) * (Real(1.0) - std::tanh((r - pparm.r_jet) / pparm.delta_jet));
  } else {
    fjet = (r <= pparm.r_jet) ? Real(1.0) : Real(0.0);
  }

  const Real rho  = fjet * rho_jet  + (Real(1.0) - fjet) * rho_amb;
  const Real uz   = fjet * uz_jet   + (Real(1.0) - fjet) * uz_amb;
  const Real eint = fjet * eint_jet + (Real(1.0) - fjet) * eint_amb;
  const Real ke   = Real(0.5) * uz * uz;

  s_ext[cls.URHO] = rho;
  s_ext[cls.UMX]  = Real(0.0);       // no radial momentum at inflow
  s_ext[cls.UMY]  = rho * uz;
  s_ext[cls.UET]  = rho * (eint + ke);
}


////////////////////////////////// tagging /////////////////////////////////////
template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level)
{
  amrex::ignore_unused(geomdata, pparm);
  if (nt <= 0) return;

  using idx = indicies_t;
  const Real rho_c = sdatafab(i, j, k, idx::URHO);
  const Real drhor = amrex::Math::abs(sdatafab(i+1,j,k,idx::URHO) - sdatafab(i-1,j,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real drhoz = amrex::Math::abs(sdatafab(i,j+1,k,idx::URHO) - sdatafab(i,j-1,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real grad_rho = std::sqrt(drhor*drhor + drhoz*drhoz);

  if (level == 0) {
    tagfab(i, j, k) = tagfab(i, j, k) || (grad_rho > Real(0.03));
  } else {
    tagfab(i, j, k) = tagfab(i, j, k) || (grad_rho > Real(0.06));
  }
}

} // namespace PROB

#endif
