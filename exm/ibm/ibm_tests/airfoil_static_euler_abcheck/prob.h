#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <Constants.h>
#include <ibm_walltypes.h>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// A/B equivalence check of the static diamond-wedge benchmark.
//
// The baseline in exm/ibm/ibm_tests/airfoil_static_euler/ enforces a
// nonzero angle of attack by *rotating the geometry* (pre-rotated .dat
// files) in a uniform +x freestream.
//
// This companion case enforces the same relative attack angle by the
// opposite convention: it uses the *un-rotated* geometry
// (diamond_wedge_mid_a00.dat) and *tilts the freestream* by ALPHA_FS
// degrees CCW, so that u_oo = U cos(ALPHA_FS), v_oo = U sin(ALPHA_FS).
//
// Purpose: confirm that the integrated forces, surface Cp, and
// extracted shock angle are nearly identical to the rotate-geometry
// case, which is a direct check on the IBM treatment of obliquely
// oriented embedded boundaries.
static constexpr Real Mach     = 2.0;
static constexpr Real ALPHA_FS = 12.0;  // deg, freestream rotated CCW
static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant/Mw;
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;
static constexpr Real PI_      = 3.14159265358979323846;

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{
  Real p_oo    = 101325.0;
  Real T_oo    = 300.0;
  Real rho_oo  = p_oo/(Rgas*T_oo);
  Real c_oo    = sqrt(gam*Rgas*T_oo);
  Real U_oo    = c_oo*Mach;
  Real u_oo    = U_oo * cos(ALPHA_FS * PI_ / 180.0);
  Real v_oo    = U_oo * sin(ALPHA_FS * PI_ / 180.0);
  Real eint_oo = rho_oo*Cv*T_oo;
  Real kin_oo  = 0.5*rho_oo*(u_oo*u_oo + v_oo*v_oo);
};

// no viscosity / conductivity for the Euler canonical benchmark
struct methodparm_t {
  public:
  static constexpr int  order = 2;
  static constexpr Real conductivity = 0.0;
  static constexpr Real viscosity    = 0.0;
  static constexpr bool use_LES = false;
};

struct skewparm_t {
  public:
  static constexpr bool dissipation = true;
  static constexpr int  order = 4;
  static constexpr Real C2skew = 1.5, C4skew = 0.016;
};

struct ibmparm_t {
  public:
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha = 0.4;

  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.4;

  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true;
};

// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// Euler: WENO-Z5 inviscid flux + no diffusive flux + no source
typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>,
               no_diffusive_t, no_source_t> ProbRHS;

// IBM wall: adiabatic slip wall — matches the canonical inviscid benchmark
typedef ibm_adiabatic_slip_wall_t<ibmparm_t, ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall, ibmparm_t, ProbClosures> ProbIB;

// Static geometry: no update needed
inline void update_geometry(Real /*time*/,
                            Vector<GeomType>& /*geom_a*/,
                            int /*ngeom*/) {}

void inline inputs() {
  amrex::Print() << " ****** Static Diamond Wedge  A/B check "
                    "(un-rotated body, freestream tilted "
                 << ALPHA_FS << " deg CCW) *******\n";
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& /*geomdata*/, ProbClosures const& cls, ProbParm const& pp) {

  // Uniform freestream tilted by ALPHA_FS (see top of file).
  state(i, j, k, cls.URHO) = pp.rho_oo;
  state(i, j, k, cls.UMX)  = pp.rho_oo * pp.u_oo;
  state(i, j, k, cls.UMY)  = pp.rho_oo * pp.v_oo;
  state(i, j, k, cls.UMZ)  = 0.0;
  state(i, j, k, cls.UET)  = pp.rho_oo * pp.eint_oo + pp.kin_oo;
}

//////////////////////////// AMR tagging ///////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int /*nt*/, auto& tagfab, const auto& sdatafab,
                  const Array4<const unsigned char>& ibfab, const auto& /*geomdata*/,
                  const ProbParm& /*pp*/, int level) {

  if (level >= 3) return;  // no refinement beyond L3

  // (1) body-adjacency tagging, level-aware ring radius:
  //      L=0 -> 5x5  (dx0 ~ 2.5e-3 m -> ring ~ 1.2e-2 m ~ chord/8)
  //      L=1 -> 3x3  (dx1 ~ 1.3e-3 m -> ring ~ 3.9e-3 m ~ chord/26)
  //      L=2 -> center cell only (L2->L3 stays tight to the surface)
  const int ring = (level == 0) ? 2 : (level == 1) ? 1 : 0;
  bool body_tag = false;
  for (int ii = -ring; ii <= ring; ii++) {
    for (int jj = -ring; jj <= ring; jj++) {
#if (AMREX_SPACEDIM == 3)
      for (int kk = -ring; kk <= ring; kk++) {
        body_tag = ibfab(i+ii, j+jj, k+kk, 1) || body_tag;
      }
#else
      body_tag = ibfab(i+ii, j+jj, k, 1) || body_tag;
#endif
    }
  }

  // (2) density-gradient sensor (normalised). Thresholds are a decade lower
  //     than multi_body (which handles moving Mach-4 shocks); here the only
  //     gradients are the wedge-generated oblique shocks and expansions.
  using PC = ProbClosures;
  const Real rhop  = sdatafab(i, j, k, PC::URHO);
  const Real drhox = std::abs(sdatafab(i+1, j, k, PC::URHO)
                              - sdatafab(i-1, j, k, PC::URHO));
  const Real drhoy = std::abs(sdatafab(i, j+1, k, PC::URHO)
                              - sdatafab(i, j-1, k, PC::URHO));
  const Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy)
                       / amrex::max(rhop, Real(1.0e-30));
  constexpr Real thr[3] = {0.03_rt, 0.08_rt, 1.0e30_rt};
  const bool grad_tag = rhofluc > thr[level];

  tagfab(i, j, k) = (body_tag || grad_tag);
}

//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real /*dratio*/,
         const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real /*s_refl*/[ProbClosures::NCONS],
         amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real /*time*/,
         amrex::GeometryData const& /*geomdata*/,
         ProbClosures const& /*closures*/, ProbParm const& pp)
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir + 1) * sgn;

  // face map: idir=0 sgn=+1 -> x-low (face=1), idir=1 sgn=+1 -> y-low
  // (face=2). Both carry supersonic inflow because the tilted freestream
  // has u>0 AND v>0, i.e. fluid enters through x-low AND y-low.
  // x-high and y-high remain outflow (zero-gradient).
  switch (face) {
    case 1:  // x-low: supersonic inflow
    case 2:  // y-low: supersonic inflow (needed for ALPHA_FS > 0)
      s_ext[URHO] = pp.rho_oo;
      s_ext[UMX]  = pp.rho_oo * pp.u_oo;
      s_ext[UMY]  = pp.rho_oo * pp.v_oo;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = pp.rho_oo * pp.eint_oo + pp.kin_oo;
      break;
    default:
      // x-hi, y-hi: outflow.
      break;
  }
}

} // namespace PROB
#endif
