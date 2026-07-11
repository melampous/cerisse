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

// 3D Cartesian underexpanded jet — NPR_0 = 3.0
// VISCOUS (Navier–Stokes) with tanh inflow profile at nozzle lip.
// Extra AMR levels (L3, L4) forced near nozzle exit for shear-layer resolution.
// L4 dx ~ 39 um  =>  256 pts/D_e at lip

using namespace amrex;
using namespace universal_constants;

namespace PROB {

static constexpr Real Mw   = 28.96e-3;
static constexpr Real gam  = 1.4;
static constexpr Real Rgas = gas_constant / Mw;
static constexpr Real Cv   = Rgas / (gam - 1.0);
static constexpr Real Cp   = gam * Cv;

struct ProbParm
{
  Real p_amb    = 101325.0;
  Real T_amb    = 300.0;
  Real Mach_amb = 0.05;

  // NPR_0 = 3.0  =>  pe/pa = 3.0 / 1.89293 = 1.5848
  Real p_jet    = 101325.0 * 3.0 / 1.89293;
  Real T_jet    = 300.0;
  Real Mach_jet = 1.0;

  Real y_jet_center = 0.0;
  Real z_jet_center = 0.0;
  Real r_jet        = 0.005;      // 5 mm radius

  // Tanh smoothing width at nozzle lip  (Version A profile)
  // fjet = 0.5*(1 - tanh((r - R) / a)),   a = delta_jet
  // 100 um  ~  2.6 L4 cells ;  visible shear-layer width ~7a ~ 700 um
  Real delta_jet    = 1.0e-4;     // [m]
};

struct skewparm_t {
  static constexpr bool dissipation = true;
  static constexpr int  order = 4;
  static constexpr Real C2skew = 1.6, C4skew = 0.008;
};

struct ibmparm_t {
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha = 0.6;
  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.6;
  static constexpr int  ghost_layers = 3;
  static constexpr int  ghost_switch = 1;
};

struct methodparm_t {
  static constexpr int  order = 2;
  static constexpr bool use_LES = false;
  static constexpr Real conductivity = 0.0262;    // W/(m·K)  air @ 300 K
  static constexpr Real viscosity    = 1.85e-5;   // Pa·s     air @ 300 K
};

using ProbClosures =
    closures_dt<indicies_t,
                visc_const_t<methodparm_t>,
                cond_const_t<methodparm_t>,
                calorifically_perfect_gas_t<indicies_t>>;

#ifndef JET_EULER_SCHEME_ID
#define JET_EULER_SCHEME_ID 0
#endif

#if (JET_EULER_SCHEME_ID == 0)
using ProbEuler = weno_t<ReconScheme::WenoZ5, ProbClosures>;
#elif (JET_EULER_SCHEME_ID == 1)
using ProbEuler = afd_hllc_teno5_t<ProbClosures>;
#elif (JET_EULER_SCHEME_ID == 2)
using ProbEuler = afd_hllc_wenoz5_t<ProbClosures>;
#else
#error "Unknown JET_EULER_SCHEME_ID"
#endif

// Navier-Stokes RHS: selectable convective flux + viscous diffusion (no source).
typedef rhs_dt<ProbEuler, viscous_t<methodparm_t, ProbClosures>, no_source_t>
    ProbRHS;


void inline inputs() {
  ProbParm p;
  const Real rho_jet = p.p_jet / (Rgas * p.T_jet);
  const Real c_jet   = std::sqrt(gam * Rgas * p.T_jet);
  const Real V_jet   = p.Mach_jet * c_jet;
  const Real D_e     = 2.0 * p.r_jet;
  const Real NPR0    = (p.p_jet / p.p_amb)
                       * std::pow(1.0 + 0.5*(gam-1.0)*p.Mach_jet*p.Mach_jet,
                                  gam/(gam-1.0));
  const Real x_MD    = 0.67 * std::sqrt(NPR0) * D_e;
  const Real Re_D    = rho_jet * V_jet * D_e / methodparm_t::viscosity;

  amrex::Print() << "\n====== 3D Underexpanded Jet  NPR_0=3  VISCOUS+tanh ======\n";
  amrex::Print() << "  D_e   = " << D_e*1e3 << " mm\n";
  amrex::Print() << "  pe/pa = " << p.p_jet/p.p_amb << "   NPR_0 = " << NPR0 << "\n";
  amrex::Print() << "  M_e   = " << p.Mach_jet << "   V_jet = " << V_jet << " m/s\n";
  amrex::Print() << "  Re_D  = " << Re_D << "\n";
  amrex::Print() << "  tanh a = " << p.delta_jet*1e6 << " um\n";
  amrex::Print() << "  Mach disk  x_MD ~ " << x_MD*1e3 << " mm\n";
  amrex::Print() << "==========================================================\n\n";
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, Array4<Real> const& state,
                    GeometryData const& geomdata, ProbClosures const& cls,
                    ProbParm const& pparm)
{
  amrex::ignore_unused(i, j, k, geomdata);

  const Real rho  = pparm.p_amb / (Rgas * pparm.T_amb);
  const Real c    = std::sqrt(gam * Rgas * pparm.T_amb);
  const Real u    = pparm.Mach_amb * c;
  const Real eint = Cv * pparm.T_amb;
  const Real rhoE = rho * (eint + 0.5_rt * u * u);

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * u;
  state(i, j, k, cls.UMY)  = 0.0_rt;
#if (AMREX_SPACEDIM == 3)
  state(i, j, k, cls.UMZ)  = 0.0_rt;
#endif
  state(i, j, k, cls.UET)  = rhoE;
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio,
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

  if (!(idir == 0 && sgn == 1)) return;

  const Real rho_amb  = pparm.p_amb / (Rgas * pparm.T_amb);
  const Real c_amb    = std::sqrt(gam * Rgas * pparm.T_amb);
  const Real u_amb    = pparm.Mach_amb * c_amb;
  const Real eint_amb = Cv * pparm.T_amb;

  const Real rho_jet  = pparm.p_jet / (Rgas * pparm.T_jet);
  const Real c_jet    = std::sqrt(gam * Rgas * pparm.T_jet);
  const Real u_jet    = pparm.Mach_jet * c_jet;
  const Real eint_jet = Cv * pparm.T_jet;

  const Real dy = x[1] - pparm.y_jet_center;
#if (AMREX_SPACEDIM == 3)
  const Real dz = x[2] - pparm.z_jet_center;
  const Real r2 = dy * dy + dz * dz;
#else
  const Real r2 = dy * dy;
#endif

  Real fjet;
  if (pparm.delta_jet > Real(0.0)) {
    const Real r = std::sqrt(r2);
    fjet = Real(0.5) * (Real(1.0) - std::tanh((r - pparm.r_jet) / pparm.delta_jet));
  } else {
    fjet = (r2 <= pparm.r_jet * pparm.r_jet) ? Real(1.0) : Real(0.0);
  }

  const Real rho  = fjet * rho_jet  + (Real(1.0) - fjet) * rho_amb;
  const Real u    = fjet * u_jet    + (Real(1.0) - fjet) * u_amb;
  const Real eint = fjet * eint_jet + (Real(1.0) - fjet) * eint_amb;
  const Real ke   = Real(0.5) * u * u;
  const Real rhoE = rho * (eint + ke);

  s_ext[cls.URHO] = rho;
  s_ext[cls.UMX]  = rho * u;
  s_ext[cls.UMY]  = Real(0.0);
#if (AMREX_SPACEDIM == 3)
  s_ext[cls.UMZ]  = Real(0.0);
#endif
  s_ext[cls.UET]  = rhoE;
}


// Tagging: levels 0-1 gradient-based (same as inviscid case),
//          levels 2-3 add forced refinement near nozzle lip.
template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level)
{
  if (nt <= 0) return;

  using idx = indicies_t;

  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();

  const Real xc = prob_lo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
  const Real yc = prob_lo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#if (AMREX_SPACEDIM == 3)
  const Real zc = prob_lo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
  const Real zc = Real(0.0);
#endif

  // ---------- Density-gradient sensor ----------
  const Real rho_c = sdatafab(i, j, k, idx::URHO);
  const Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,idx::URHO) - sdatafab(i-1,j,k,idx::URHO)) / (Real(2.0) * rho_c);
  const Real drhoy = amrex::Math::abs(sdatafab(i,j+1,k,idx::URHO) - sdatafab(i,j-1,k,idx::URHO)) / (Real(2.0) * rho_c);
#if (AMREX_SPACEDIM == 3)
  const Real drhoz = amrex::Math::abs(sdatafab(i,j,k+1,idx::URHO) - sdatafab(i,j,k-1,idx::URHO)) / (Real(2.0) * rho_c);
  const Real grad_rho = std::sqrt(drhox*drhox + drhoy*drhoy + drhoz*drhoz);
#else
  const Real grad_rho = std::sqrt(drhox*drhox + drhoy*drhoy);
#endif

  // ---------- Radial distance from jet axis & lip proximity ----------
  const Real dy_lip = yc - pparm.y_jet_center;
  const Real dz_lip = zc - pparm.z_jet_center;
  const Real r_lip  = std::sqrt(dy_lip*dy_lip + dz_lip*dz_lip);
  const Real dr_lip = amrex::Math::abs(r_lip - pparm.r_jet);  // distance from lip circle

  // Level 0 → creates L1:  gradient-based
  if (level == 0) {
    tagfab(i,j,k) = tagfab(i,j,k) || (grad_rho > Real(0.03));
  }

  // Level 1 → creates L2:  gradient-based
  if (level == 1) {
    tagfab(i,j,k) = tagfab(i,j,k) || (grad_rho > Real(0.06));
  }

  // Level 2 → creates L3:  gradient + forced lip annulus
  //   lip region:  x < 2 D_e (20 mm)  and  |r - r_jet| < 2 mm
  if (level == 2) {
    tagfab(i,j,k) = tagfab(i,j,k) || (grad_rho > Real(0.06));
    tagfab(i,j,k) = tagfab(i,j,k) || (xc < Real(0.020) && dr_lip < Real(0.002));
  }

  // Level 3 → creates L4:  tighter lip focus only
  //   lip region:  x < 1 D_e (10 mm)  and  |r - r_jet| < 1 mm
  if (level == 3) {
    tagfab(i,j,k) = tagfab(i,j,k) || (xc < Real(0.010) && dr_lip < Real(0.001));
  }
}


} // namespace PROB

#endif
