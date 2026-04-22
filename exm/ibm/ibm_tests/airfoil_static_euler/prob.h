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

// Canonical static diamond-wedge benchmark (Euler).
//
// Freestream: M = 2, fixed along +x.
// Angle of attack is realised by rotating the geometry clockwise about
// the origin (see gen_rotated_airfoils.py); the solver sees a static
// body in a uniform +x inflow at every AoA, so ib.move = 0.
static constexpr Real Mach = 2.0;
static constexpr Real Mw   = 28.96e-3;
static constexpr Real gam  = 1.4;
static constexpr Real Rgas = gas_constant/Mw;
static constexpr Real Cv   = Rgas/(gam - 1.0);
static constexpr Real Cp   = gam*Cv;

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{
  Real p_oo    = 101325.0;
  Real T_oo    = 300.0;
  Real rho_oo  = p_oo/(Rgas*T_oo);
  Real c_oo    = sqrt(gam*Rgas*T_oo);
  Real u_oo    = c_oo*Mach;
  Real eint_oo = rho_oo*Cv*T_oo;
  Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;
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
  amrex::Print() << " ****** Static Diamond Wedge (Euler, M=2) *******\n";
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& /*geomdata*/, ProbClosures const& cls, ProbParm const& pp) {

  // Uniform +x freestream everywhere; the wedge is embedded via IBM.
  state(i, j, k, cls.URHO) = pp.rho_oo;
  state(i, j, k, cls.UMX)  = pp.rho_oo * pp.u_oo;
  state(i, j, k, cls.UMY)  = 0.0;
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

  switch (face) {
    case 1:  // x-low: supersonic inflow (Dirichlet)
      s_ext[URHO] = pp.rho_oo;
      s_ext[UMX]  = pp.rho_oo * pp.u_oo;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = pp.rho_oo * pp.eint_oo + pp.kin_oo;
      break;
    default:
      // x-hi, y-lo, y-hi: outflow (zero-gradient). s_ext already filled with s_int.
      break;
  }
}

} // namespace PROB
#endif
