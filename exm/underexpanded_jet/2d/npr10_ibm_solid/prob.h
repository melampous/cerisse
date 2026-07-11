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

// =============================================================================
//  2D Axisymmetric (RZ) underexpanded SONIC jet emitted from an IBM SOLID plate.
//  Control case for the IBM jet-emission modelling: identical NPR=3 conditions
//  to the boundary-inflow case npr3_lev3_euler_axisym, but the jet now issues
//  from the top face of a flat solid plate (plate.dat) via the IBM wall model.
//  Validation target: Ashkenas-Sherman Mach-disk x_MD = 0.67*sqrt(NPR0)*De
//                     = 0.67*sqrt(3)*10mm = 11.6 mm.
//  x=r (idir0, r=0 axis), y=z (idir1); jet fires +z (plate face normal +z).
// =============================================================================

namespace PROB {

static constexpr Real Mw   = 28.96e-3;
static constexpr Real gam  = 1.4;
static constexpr Real Rgas = gas_constant / Mw;
static constexpr Real Cv   = Rgas / (gam - 1.0);
static constexpr Real Cp   = gam * Cv;

static constexpr Real p_amb    = 101325.0;
static constexpr Real T_amb    = 300.0;
static constexpr Real Mach_amb = 0.05;            // weak +z coflow
static const     Real rho_amb  = p_amb / (Rgas * T_amb);
static const     Real c_amb    = std::sqrt(gam * Rgas * T_amb);
static const     Real u_amb    = Mach_amb * c_amb;

// NPR_0 = p0/pa = 3  ->  sonic (M=1) exit pe/pa = 3/1.89293 = 1.585 (underexpanded)
static constexpr Real p_jet = 101325.0 * 10.0 / 1.89293;  // 1.607e5 Pa
static constexpr Real T_jet = 300.0;
static constexpr Real V_jet = 347.19;                    // M=1: sqrt(gam*R*300)

// jet-hole on the plate top face (plate.dat): z = z_jet, r < r_jet
static constexpr Real z_jet  = 0.0;
static constexpr Real r_jet  = 0.005;     // 5 mm, D_e = 10 mm
static constexpr Real ztol   = 0.0015;

struct ProbParm
{
  Real p_oo    = p_amb;
  Real T_oo    = T_amb;
  Real rho_oo  = rho_amb;
  Real c_oo    = c_amb;
  Real u_oo    = 0.0;          // radial
  Real v_oo    = u_amb;        // axial +z (weak coflow)
  Real eint_oo = rho_amb * Cv * T_amb;
  Real kin_oo  = 0.5 * rho_amb * u_amb * u_amb;

  static constexpr Real zj_loc = z_jet;
  static constexpr Real rj_val = r_jet;
  static constexpr Real ztol_v = ztol;
  static constexpr Real pj_val = p_jet;
  static constexpr Real Tj_val = T_jet;
  static constexpr Real Vj_val = V_jet;
};

struct methodparm_t {
  public:
  static constexpr int  order        = 2;
  static constexpr Real conductivity = 1.0;
  static constexpr Real viscosity    = 1.0;
  static constexpr bool use_LES      = false;
};

struct ibmparm_t {
  public:
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha        = 0.6;
  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf        = 0.6;
  static constexpr int  ghost_layers      = 1;
  static constexpr bool interior_is_solid = true;
};

typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
// inviscid Euler (no_diffusive) + WENO-Z5 shock capture
typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,
               no_source_t> ProbRHS;

template < typename param, typename cls_t > class ibm_user_t;
typedef ibm_user_t<ProbParm, ProbClosures>             TypeWall;
typedef ibm_solver_t<TypeWall, ibmparm_t, ProbClosures> ProbIB;

inline void update_geometry(Real, Vector<GeomType>&, int) {}

void inline inputs() {
  const Real NPR0 = (p_jet / p_amb) * std::pow(1.0 + 0.2 * 1.0, 3.5);
  const Real x_MD = 0.67 * std::sqrt(NPR0) * 2.0 * r_jet;
  amrex::Print() << " ****** IBM-solid underexpanded jet (RZ), NPR0=" << NPR0 << " ******\n";
  amrex::Print() << "  jet(M=1): p_e=" << p_jet << " T_e=" << T_jet << " V=" << V_jet
                 << "  pe/pa=" << p_jet/p_amb << "  De=" << 2*r_jet*1000 << "mm\n";
  amrex::Print() << "  Ashkenas-Sherman x_MD = " << x_MD*1000 << " mm (x_MD/De="
                 << x_MD/(2*r_jet) << ")\n";
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const&, ProbClosures const& cls,
                   ProbParm const& pparm)
{
  state(i, j, k, cls.URHO) = pparm.rho_oo;
  state(i, j, k, cls.UMX)  = 0.0;
  state(i, j, k, cls.UMY)  = pparm.rho_oo * pparm.v_oo;
  state(i, j, k, cls.UMZ)  = 0.0;
  state(i, j, k, cls.UET)  = pparm.eint_oo + pparm.kin_oo;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int, auto& tagfab, const auto& sdatafab,
                  const Array4<const unsigned char>& ibfab, const auto& geomdata,
                  const ProbParm&, int level)
{
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();
  const Real rc = prob_lo[0] + (i + Real(0.5)) * dx[0];   // r
  const Real zc = prob_lo[1] + (j + Real(0.5)) * dx[1];   // z

  // (1) GEOMETRY-ANCHORED jet-column refinement -- grid-INDEPENDENT, so the AMR
  //     regrid cascade CONVERGES.  A pure-gradient criterion (below) oscillates
  //     at t=0: the only density gradient is the artificial IBM jet-jump (ghost
  //     cells carrying the injected sonic state), which FillPatch re-interpolates
  //     differently on every regrid pass -> tagged cell set keeps changing ->
  //     grids never stabilise -> endless post_regrid/computeSurfIndices cascade
  //     (the IBM+AMR "init/first-step loop").  Nested physical boxes per level
  //     resolve the barrel shock + shock-cell train down to the finest level
  //     deterministically and cheaply.
  const Real r_ref = (level == 0) ? Real(0.024)
                   : (level == 1) ? Real(0.018)
                                  : Real(0.014);
  const Real z_top = (level == 0) ? Real(0.10)
                   : (level == 1) ? Real(0.075)
                                  : Real(0.055);
  if (rc < r_ref && zc > Real(-0.001) && zc < z_top) tagfab(i, j, k) = true;

  // (2) Density-gradient tag for shock structure -- but SKIP any stencil that
  //     straddles the IBM solid/ghost band (comp 0 = solid, comp 1 = ghost).
  //     That artificial jump is what drives the tag oscillation; refining the
  //     real shocks away from the wall is deterministic and safe.
  const bool stencil_clean =
        !ibfab(i,   j,   k, 0) && !ibfab(i,   j,   k, 1)
     && !ibfab(i+1, j,   k, 0) && !ibfab(i-1, j,   k, 0)
     && !ibfab(i,   j+1, k, 0) && !ibfab(i,   j-1, k, 0)
     && !ibfab(i+1, j,   k, 1) && !ibfab(i-1, j,   k, 1)
     && !ibfab(i,   j+1, k, 1) && !ibfab(i,   j-1, k, 1);
  if (stencil_clean) {
    const Real rho  = sdatafab(i, j, k, 0);
    const Real drr  = sdatafab(i + 1, j, k, 0) - sdatafab(i - 1, j, k, 0);
    const Real drz  = sdatafab(i, j + 1, k, 0) - sdatafab(i, j - 1, k, 0);
    const Real g    = std::sqrt(drr*drr + drz*drz) / (rho + Real(1e-10));
    const Real thr  = (level == 0) ? Real(0.03) : Real(0.06);
    if (g > thr) tagfab(i, j, k) = true;
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real,
         const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const Real,
         GeometryData const&, ProbClosures const&, ProbParm const& pparm)
{
  // all open boundaries: ambient / outflow (foextrap-like copy)
  for (int n = 0; n < ProbClosures::NCONS; ++n) s_ext[n] = s_int[n];
}

// Wall model: slip plate; sonic (M=1) jet on the top-face hole
// (|z - z_jet| < ztol, r < r_jet, outward normal +z).
template < typename param, typename cls_t >
class ibm_user_t
{
  public:
  static constexpr int eorder_tparm = ibmparm_t::extrap_order;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
                             const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
                             const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t1*/,
                             const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
                             Array2D<Real,0,eorder_tparm+1,
                                     0,cls_t::NPRIM-1>& q,
                             int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    // --- Default: inviscid SLIP wall (plate) ---
    q(1,cls_t::QU)    = 0.0;                 // un  = 0
    q(1,cls_t::QV)    = q(2,cls_t::QV);      // tangential slip
    q(1,cls_t::QW)    = q(2,cls_t::QW);
    q(1,cls_t::QPRES) = q(2,cls_t::QPRES);
    q(1,cls_t::QT)    = q(2,cls_t::QT);

    // --- Jet: sonic (M=1) on the plate top-face hole, blowing +z ---
    const Real r = xyz(0);
    const Real z = xyz(1);
    const bool on_jet = (amrex::Math::abs(z - param::zj_loc) < param::ztol_v)
                        && (r < param::rj_val);
    const bool fwd    = (norm(1) > Real(0.3));   // outward normal +z
    if (on_jet && fwd) {
      q(1,cls_t::QU)    = param::Vj_val;     // sonic, blowing along normal (+z)
      q(1,cls_t::QV)    = 0.0;
      q(1,cls_t::QW)    = 0.0;
      q(1,cls_t::QPRES) = param::pj_val;
      q(1,cls_t::QT)    = param::Tj_val;
    }
    ibm_copy_species<cls_t, eorder_tparm>(q);
  }
};

} // namespace PROB

#endif
