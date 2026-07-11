#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <cmath>

#include <Closures.h>
#include <RHS.h>
#include <bc_types.h>

// =====================================================================
//  2-D co-rotating Lamb-Oseen vortex pair (CVP) aeroacoustic showcase.
//
//  Two identical SAME-sign vortices orbit their midpoint and radiate a
//  rotating m=2 quadrupole (spiral) acoustic field. The radiated sound
//  is O(rho0*Vorb^4/c^2) ~ 1e-4, i.e. 2-3 orders weaker than the
//  hydrodynamic core dip, so it is extremely sensitive to boundary
//  reflection. Purpose: showcase non-reflecting BC superiority on a
//  NARROW domain (foextrap=code 2 vs characteristic code-7 vs LODI)
//  against a WIDE-domain "truth".
//
//  Non-dim: rho0=1, c0=1, gamma=1.4 => p0=1/gamma=0.7142857.
//  Background is QUIESCENT (u_inf=v_inf=0), so far-field target on all
//  4 faces is (rho0, 0, 0, p0).  NO IBM, Euler (no_diffusive_t).
// =====================================================================

using namespace amrex;

namespace PROB {

struct ProbParm {
  // --- core physical knobs (runtime-overridable via ParmParse "prob") ---
  Real gamma = 1.4;
  Real rho0  = 1.0;
  Real c0    = 1.0;        // reference sound speed (derived from rho0,p0)
  Real p0    = 1.0 / 1.4;  // = 1/gamma = 0.7142857 (quiescent background)

  Real Gamma = 1.2566370614359172;  // circulation, = M_rot*c0*2*pi*b, M_rot=0.20
  Real b     = 1.0;                  // vortex separation (length unit)
  Real rc    = 0.25;                 // Lamb-Oseen core radius (b/rc = 4)

  // --- BC control ---
  // bc_mode: 1 = characteristic far-field (code-7 -> bc_nscbc_farfield)
  //          3 = multidim LODI far-field (code-7 -> bc_nscbc_lodi_farfield)
  // foextrap is selected purely by setting cns.lo_bc/hi_bc = 2 2 in inputs
  // (AMReX handles it; bcnormal is then never called).
  int  bc_mode  = 1;
  int  use_lodi = 0;       // 1 -> enable bcnormal_lodi (set together with bc_mode=3)

  // LODI tuning (Lodato transverse term; only used when use_lodi=1)
  Real lodi_pressure_relax   = 0.0;   // 0 -> fully non-reflecting normal char
  Real lodi_transverse_relax = 0.25;  // Lodato beta for oblique spiral fronts

  // --- AMR tagging knobs ---
  // Tag where |omega_z| exceeds a fraction of the peak core vorticity, and
  // (optionally) restrict refinement to a disk around the orbiting pair.
  Real tag_vort_frac = 0.15;   // fraction of (Gamma/(pi*rc^2)) to tag on
  Real tag_radius    = 1.5;    // only refine within this radius of origin (b units)
  Real tag_graddens  = 0.05;   // |grad rho|/rho threshold (backup tag)

  // --- derived (filled in constructor) ---
  Real omega   = 0.0;   // orbital angular velocity Gamma/(pi*b^2)
  Real Trot    = 0.0;   // rotation period 2*pi/omega
  Real f_ac    = 0.0;   // acoustic frequency omega/pi
  Real lambda  = 0.0;   // acoustic wavelength c0/f_ac
  Real D0      = 0.0;   // core pressure deficit p0-p(r=0)

  ProbParm()
  {
    ParmParse pp("prob");
    pp.query("gamma", gamma);
    pp.query("rho0",  rho0);
    pp.query("p0",    p0);
    pp.query("Gamma", Gamma);
    pp.query("b",     b);
    pp.query("rc",    rc);
    pp.query("bc_mode",  bc_mode);
    pp.query("use_lodi", use_lodi);
    pp.query("lodi_pressure_relax",   lodi_pressure_relax);
    pp.query("lodi_transverse_relax", lodi_transverse_relax);
    pp.query("tag_vort_frac", tag_vort_frac);
    pp.query("tag_radius",    tag_radius);
    pp.query("tag_graddens",  tag_graddens);

    c0 = std::sqrt(gamma * p0 / rho0);

    const Real pi = Real(3.14159265358979323846);
    omega  = Gamma / (pi * b * b);
    Trot   = (omega > Real(0.0)) ? Real(2.0) * pi / omega : Real(0.0);
    f_ac   = omega / pi;
    lambda = (f_ac > Real(0.0)) ? c0 / f_ac : Real(0.0);
    // Core (r->0) limit of the superposed-deficit formula for one vortex:
    //   D0_single = rho0*Gamma^2*ln2/(4*pi^2*rc^2)
    D0 = rho0 * Gamma * Gamma * std::log(Real(2.0)) /
         (Real(4.0) * pi * pi * rc * rc);
  }
};

// --------------------------------------------------------------------
// Numerical method: inviscid (Euler), WenoZ5, no source.  Matches the
// bc_native template's closures_dt / weno_t<WenoZ5> / rhs_dt pattern.
// --------------------------------------------------------------------
using ProbClosures = closures_dt<
    indicies_t,
    calorifically_perfect_gas_t<indicies_t>>;

using ProbRHS = rhs_dt<
    weno_t<ReconScheme::WenoZ5, ProbClosures>,
    no_diffusive_t,
    no_source_t>;

using GlobalBC = manual_bc_t<ProbClosures>;

inline void inputs() {}

// --------------------------------------------------------------------
// Helper: write conserved state from (rho,u,v,p).  Mirrors bc_native.
// --------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_state(const ProbClosures& cls, const ProbParm& pparm,
               const Real rho, const Real u, const Real v, const Real p,
               Real* s)
{
  Real y[NUM_SPECIES] = {Real(1.0)};
  Real eint = Real(0.0);
  cls.RYP2E(rho, y, p, eint);

  s[ProbClosures::URHO] = rho;
  s[ProbClosures::UMX]  = rho * u;
  s[ProbClosures::UMY]  = rho * v;
  s[ProbClosures::UMZ]  = Real(0.0);
  s[ProbClosures::UET]  = rho * eint + Real(0.5) * rho * (u * u + v * v);
}

// --------------------------------------------------------------------
// Closed-form single-vortex pressure deficit D(r) = p0 - p(r) from the
// incompressible radial balance dp/dr = rho0*u_theta^2/r for the
// Lamb-Oseen profile u_theta(r)=(Gamma/(2*pi*r))*(1-exp(-r^2/rc^2)):
//
//   x = r^2/rc^2
//   boundary = (1 - e^{-x})^2 / x
//   g(x) = E1(x) - E1(2x) = integral_x^{2x} e^{-t}/t dt   (GL4 quadrature)
//   J = boundary + 2*g
//   D(r) = rho0 * (Gamma^2/(4*pi^2)) * J / (2*rc^2)
//
// GL4 on [x,2x]: t = 1.5x + 0.5x*xi, g = 0.5x*sum(wi*exp(-t)/t).
// Verified GPU-safe (only exp), rel err <2e-6 for x<=1, negligible far field.
// --------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real vortex_deficit(const Real r, const ProbParm& pparm)
{
  const Real pi   = Real(3.14159265358979323846);
  const Real rc2  = pparm.rc * pparm.rc;
  Real x = r * r / rc2;
  if (x < Real(1.0e-12)) x = Real(1.0e-12);  // guard 1/x and 1/t

  const Real ex   = std::exp(-x);
  const Real omex = Real(1.0) - ex;
  const Real boundary = omex * omex / x;

  // 4-point Gauss-Legendre on [x, 2x] of e^{-t}/t  =>  g = E1(x)-E1(2x).
  const Real xi[4] = { Real(-0.8611363115940526), Real(-0.3399810435848563),
                       Real( 0.3399810435848563), Real( 0.8611363115940526) };
  const Real wi[4] = { Real( 0.3478548451374538), Real( 0.6521451548625461),
                       Real( 0.6521451548625461), Real( 0.3478548451374538) };
  Real g = Real(0.0);
  for (int q = 0; q < 4; ++q) {
    const Real t = Real(1.5) * x + Real(0.5) * x * xi[q];
    g += wi[q] * std::exp(-t) / t;
  }
  g *= Real(0.5) * x;

  const Real J = boundary + Real(2.0) * g;
  return pparm.rho0 * (pparm.Gamma * pparm.Gamma / (Real(4.0) * pi * pi)) * J /
         (Real(2.0) * rc2);
}

// --------------------------------------------------------------------
// Initial condition: superpose two SAME-sign (+Gamma) Lamb-Oseen
// vortices at (+b/2,0) and (-b/2,0).  Velocity = vector superposition;
// pressure = p0 - D(r1) - D(r2) (scalar deficit superposition); density
// isentropic.  Energy via RYP2E in set_state.
// --------------------------------------------------------------------
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata,
                   ProbClosures const& cls,
                   ProbParm const& pparm)
{
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();

  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];

  const Real pi   = Real(3.14159265358979323846);
  const Real rc2  = pparm.rc * pparm.rc;
  const Real xc[2] = { Real(0.5) * pparm.b, -Real(0.5) * pparm.b };
  const Real yc[2] = { Real(0.0), Real(0.0) };

  Real u = Real(0.0);
  Real v = Real(0.0);
  Real Dsum = Real(0.0);   // sum over vortices of  rho0 * INT_r^inf u_theta^2/r' dr'

  for (int m = 0; m < 2; ++m) {
    const Real dxk = x - xc[m];
    const Real dyk = y - yc[m];
    const Real rk  = std::sqrt(dxk * dxk + dyk * dyk);
    if (rk > Real(1.0e-9)) {
      // single-vortex tangential speed
      const Real utheta = (pparm.Gamma / (Real(2.0) * pi * rk)) *
                          (Real(1.0) - std::exp(-rk * rk / rc2));
      // counter-clockwise (+Gamma) Cartesian components
      u += -utheta * dyk / rk;
      v +=  utheta * dxk / rk;
    }
    // accumulate this vortex's radial-balance integral (= rho0 * INT u^2/r dr)
    Dsum += vortex_deficit(rk, pparm);
  }

  // CONSISTENT COMPRESSIBLE isentropic radial equilibrium.  For p/rho^gamma=const
  //   dp/dr = rho u^2/r   <=>   d(c^2)/dr = (gamma-1) u^2/r
  //   => c^2(r) = c0^2 - (gamma-1) * INT_r^inf u^2/r' dr'
  // then p = p0 (c^2/c0^2)^{gamma/(gamma-1)},  rho = rho0 (c^2/c0^2)^{1/(gamma-1)}.
  // (The old incompressible deficit p=p0-D was inconsistent at M_swirl~0.5 and
  //  launched a strong spurious startup pulse that swamps the radiated sound.)
  const Real c0sq = pparm.gamma * pparm.p0 / pparm.rho0;
  Real csq = c0sq - (pparm.gamma - Real(1.0)) * Dsum / pparm.rho0;
  csq = amrex::max(csq, Real(1.0e-6) * c0sq);            // positivity floor
  const Real ratio = csq / c0sq;                         // (c/c0)^2
  const Real p   = pparm.p0   *
      std::pow(ratio, pparm.gamma / (pparm.gamma - Real(1.0)));
  const Real rho = pparm.rho0 *
      std::pow(ratio, Real(1.0) / (pparm.gamma - Real(1.0)));

  // Write conserved state component-by-component. NOTE: an Array4's
  // components are strided, NOT contiguous after (i,j,k), so we cannot
  // pass &state(i,j,k,0) to the flat-pointer set_state helper.
  Real yspec[NUM_SPECIES] = {Real(1.0)};
  Real eint = Real(0.0);
  cls.RYP2E(rho, yspec, p, eint);

  state(i, j, k, ProbClosures::URHO) = rho;
  state(i, j, k, ProbClosures::UMX)  = rho * u;
  state(i, j, k, ProbClosures::UMY)  = rho * v;
  state(i, j, k, ProbClosures::UMZ)  = Real(0.0);
  state(i, j, k, ProbClosures::UET)  =
      rho * eint + Real(0.5) * rho * (u * u + v * v);
}

// --------------------------------------------------------------------
// bcnormal: QUIESCENT far-field on every ext_dir (code-7) face.
// Outward-normal convention (matches manual_bc / bc_native): n[idir] =
// -sgn (lo face sgn=+1 -> n=-1; hi face sgn=-1 -> n=+1), tangential 0.
// Target state is the quiescent background (rho0, 0, 0, p0) on ALL faces.
// foextrap faces (code 2) never reach here.
// --------------------------------------------------------------------
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
  // outward normal of this face
  const Real nx = (idir == 0) ? -Real(sgn) : Real(0.0);
  const Real ny = (idir == 1) ? -Real(sgn) : Real(0.0);
  const Real nz = (idir == 2) ? -Real(sgn) : Real(0.0);

  GlobalBC::bc_nscbc_farfield(nx, ny, nz, &closures,
                              pparm.rho0, Real(0.0), Real(0.0), Real(0.0),
                              pparm.p0, s_int, s_ext);
}

// --------------------------------------------------------------------
// bcnormal_lodi: multidim LODI far-field on every code-7 face when
// use_lodi=1.  Quiescent target (rho0,0,0,0,p0), l_ref = domain size.
// Returns false (-> fall back to bc_nscbc_farfield) when use_lodi=0.
// --------------------------------------------------------------------
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
  if (pparm.use_lodi == 0) return false;

  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();

  GlobalBC::nscbc_lodi_bc_parm_t bp;
  bp.rho_inf = pparm.rho0;
  bp.u_inf   = Real(0.0);
  bp.v_inf   = Real(0.0);
  bp.w_inf   = Real(0.0);
  bp.p_inf   = pparm.p0;
  bp.pressure_relax   = pparm.lodi_pressure_relax;    // 0 = non-reflecting
  bp.transverse_relax = pparm.lodi_transverse_relax;  // Lodato beta
  bp.l_ref   = prob_hi[idir] - prob_lo[idir];
  bp.use_transverse        = 1;
  bp.use_acoustic_gate     = 1;
  bp.use_ghost_projection  = 1;
  bp.use_shock_sensor      = 0;   // smooth acoustics, no shock
  bp.use_convective_sensor = 0;
  bp.use_convective_extrapolation = 0;
  bp.use_pressure_relax_sensor    = 0;

  return GlobalBC::bc_nscbc_lodi_farfield(
      iv, state, s_ext, idir, sgn, geomdata, &closures, bp);
#endif
}

// --------------------------------------------------------------------
// AMR tagging: refine on the two vortex cores.  Tag where |omega_z|
// exceeds tag_vort_frac * (Gamma/(pi*rc^2)) (cores) OR the density
// gradient is large, restricted to a disk of radius tag_radius around
// the orbit center so the far acoustic field is NOT over-refined.
// Signature matches the non-IBM convention (covo / sod).
// NOTE: conserved index ordering here is URHO=4, UMX=0, UMY=1 (Index.h);
// use the named ProbClosures constants, not raw 0.
// --------------------------------------------------------------------
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int /*nt_level*/, auto& tagfab,
                  const auto& sdatafab, const auto& geomdata,
                  const ProbParm& pparm, int /*level*/)
{
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();

  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];

  // restrict refinement to the orbiting region (vortex centers at r=b/2)
  const Real rr = std::sqrt(x * x + y * y);
  if (rr > pparm.tag_radius) return;

  // local velocity from conserved state (URHO is index 4, momenta 0/1)
  auto vel = [&](int ii, int jj, Real& uu, Real& vv) {
    const Real rho = sdatafab(ii, jj, k, ProbClosures::URHO);
    const Real rinv = Real(1.0) / amrex::max(rho, Real(1.0e-12));
    uu = sdatafab(ii, jj, k, ProbClosures::UMX) * rinv;
    vv = sdatafab(ii, jj, k, ProbClosures::UMY) * rinv;
  };

  Real up, um, vp, vm;
  vel(i + 1, j, up, vp);
  vel(i - 1, j, um, vm);
  const Real dvdx = (vp - vm) / (Real(2.0) * dx[0]);

  Real up2, um2, vp2, vm2;
  vel(i, j + 1, up2, vp2);
  vel(i, j - 1, um2, vm2);
  const Real dudy = (up2 - um2) / (Real(2.0) * dx[1]);

  const Real omega_z = dvdx - dudy;

  const Real pi = Real(3.14159265358979323846);
  const Real omega_ref = pparm.Gamma / (pi * pparm.rc * pparm.rc);
  if (amrex::Math::abs(omega_z) > pparm.tag_vort_frac * omega_ref) {
    tagfab(i, j, k) = true;
    return;
  }

  // backup: density gradient (catches the core dip rim)
  const Real rho = amrex::max(sdatafab(i, j, k, ProbClosures::URHO),
                              Real(1.0e-12));
  const Real drx = sdatafab(i + 1, j, k, ProbClosures::URHO) -
                   sdatafab(i - 1, j, k, ProbClosures::URHO);
  const Real dry = sdatafab(i, j + 1, k, ProbClosures::URHO) -
                   sdatafab(i, j - 1, k, ProbClosures::URHO);
  const Real gmag = std::sqrt(drx * drx + dry * dry);
  if (gmag / rho > pparm.tag_graddens) {
    tagfab(i, j, k) = true;
  }
}

// no extra cell-wise source term (Euler, no_source_t handles the RHS)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int, int, int, const auto&, const auto&,
                 const ProbParm&, ProbClosures const&, auto const) {}

} // namespace PROB

#endif
