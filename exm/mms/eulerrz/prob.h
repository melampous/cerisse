#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>

#include <Closures.h>
#include <RHS.h>

// 2D axisymmetric (RZ) Euler MMS, axis-regular manufactured solution.
// env(r) = cos^2(pi r / (2R)), R = 0.5; z-period 1; U0 = 0.5.
// rho,u_z,p = const + A*env*cos(2 pi z);  u_r = 0.05*(r/R)*env*sin(2 pi z).
// All perturbations AND their radial slopes vanish at r = R, so the outer
// ghost state is the constant (1, 0, U0, 1) to second order in curvature.
// Sources = sympy residual of the continuous axisymmetric Euler equations
// (gen_rz_mms.py); independent of the discrete RZ flux splitting.

using namespace amrex;

namespace PROB {

static constexpr Real gam_mms = 1.4;
static constexpr Real U0_mms  = 0.5;

struct ProbParm {};

struct methodparm_t {
  static constexpr int  order = 2;
  static constexpr bool use_LES = false;
  static constexpr Real conductivity = 0.0262;
  static constexpr Real viscosity   = 1.85e-5;
};

using ProbClosures = closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                                 calorifically_perfect_gas_t<indicies_t>>;

template <typename cls_t> class user_source_t;

typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,
               user_source_t<ProbClosures>> ProbRHS;

void inline inputs() {
  amrex::Print() << " MMS Euler RZ (axis-regular manufactured solution) \n";
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mms_exact(const Real r, const Real z, Real& rho_o, Real& ur_o, Real& uz_o,
          Real& p_o) {
  rho_o = (3.0/20.0)*pow(cos(M_PI*r), 2)*cos(2*M_PI*z) + 1;
  ur_o  = (1.0/10.0)*r*sin(2*M_PI*z)*pow(cos(M_PI*r), 2);
  uz_o  = (1.0/10.0)*pow(cos(M_PI*r), 2)*cos(2*M_PI*z) + 1.0/2.0;
  p_o   = (1.0/10.0)*pow(cos(M_PI*r), 2)*cos(2*M_PI*z) + 1;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const& state,
              GeometryData const& geomdata, ProbClosures const& cls,
              ProbParm const& pparm) {
  amrex::ignore_unused(k, pparm);
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx = geomdata.CellSize();
  const Real r = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real z = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
  Real rho, ur, uz, p;
  mms_exact(r, z, rho, ur, uz, p);
  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * ur;
  state(i, j, k, cls.UMY)  = rho * uz;
  state(i, j, k, cls.UMZ)  = Real(0.0);
  state(i, j, k, cls.UET)  = p / (gam_mms - Real(1.0))
                           + Real(0.5) * rho * (ur * ur + uz * uz);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio,
         const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS],
         amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& geomdata, ProbClosures const& cls,
         ProbParm const& pparm) {
  amrex::ignore_unused(dratio, s_refl, time, geomdata, pparm);
  for (int n = 0; n < ProbClosures::NCONS; ++n) s_ext[n] = s_int[n];
  if (idir == 0 && sgn == -1) {
    // Exact manufactured state at the ghost location (the manufactured
    // fields are smooth for all r, so this is exact Dirichlet data).
    Real rho, ur, uz, p;
    mms_exact(x[0], x[1], rho, ur, uz, p);
    s_ext[cls.URHO] = rho;
    s_ext[cls.UMX]  = rho * ur;
    s_ext[cls.UMY]  = rho * uz;
    s_ext[cls.UMZ]  = Real(0.0);
    s_ext[cls.UET]  = p / (gam_mms - Real(1.0))
                    + Real(0.5) * rho * (ur * ur + uz * uz);
  }
}

template <typename cls_t>
class user_source_t {
 public:
  void inline src(const amrex::Geometry& geomdata, const amrex::MFIter& mfi,
                  const amrex::Array4<const amrex::Real>& prims,
                  const amrex::Array4<amrex::Real>& rhs, const cls_t* cls_d,
                  amrex::Real dt, amrex::Real time) {
    amrex::ignore_unused(prims, cls_d, dt, time);
    const Box bx = mfi.tilebox();
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    using idx = indicies_t;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real r = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
      const Real z = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
      const Real src_rho = (1.0/100.0)*(-6*M_PI*r*sin(M_PI*r)*pow(cos(M_PI*r), 2)*cos(2*M_PI*z) - 20*M_PI*r*sin(M_PI*r) - 6*M_PI*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) + 3*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) - 35*M_PI*cos(M_PI*r) + 20*cos(M_PI*r))*sin(2*M_PI*z)*cos(M_PI*r);
      const Real src_mr  = (1.0/2000.0)*(-18*M_PI*pow(r, 2)*sin(M_PI*r)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 4)*cos(2*M_PI*z) - 80*M_PI*pow(r, 2)*sin(M_PI*r)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 2) - 12*M_PI*r*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 5)*cos(2*M_PI*z) + 9*r*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 5)*cos(2*M_PI*z) - 70*M_PI*r*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 3) + 60*r*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 3) + 6*M_PI*r*pow(cos(M_PI*r), 5)*pow(cos(2*M_PI*z), 3) + 70*M_PI*r*pow(cos(M_PI*r), 3)*pow(cos(2*M_PI*z), 2) + 200*M_PI*r*cos(M_PI*r)*cos(2*M_PI*z) - 400*M_PI*sin(M_PI*r)*cos(2*M_PI*z))*cos(M_PI*r);
      const Real src_mz  = (1.0/1000.0)*(-9*M_PI*r*sin(M_PI*r)*pow(cos(M_PI*r), 4)*pow(cos(2*M_PI*z), 2) - 70*M_PI*r*sin(M_PI*r)*pow(cos(M_PI*r), 2)*cos(2*M_PI*z) - 100*M_PI*r*sin(M_PI*r) - 9*M_PI*pow(cos(M_PI*r), 5)*pow(cos(2*M_PI*z), 2) + 3*pow(cos(M_PI*r), 5)*pow(cos(2*M_PI*z), 2) - 100*M_PI*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) + 35*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) - 475*M_PI*cos(M_PI*r) + 100*cos(M_PI*r))*sin(2*M_PI*z)*cos(M_PI*r);
      const Real src_E   = (1.0/20000.0)*(-12*M_PI*pow(r, 3)*sin(M_PI*r)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 6)*cos(2*M_PI*z) - 60*M_PI*pow(r, 3)*sin(M_PI*r)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 4) - 6*M_PI*pow(r, 2)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 7)*cos(2*M_PI*z) + 6*pow(r, 2)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 7)*cos(2*M_PI*z) - 35*M_PI*pow(r, 2)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 5) + 40*pow(r, 2)*pow(sin(2*M_PI*z), 2)*pow(cos(M_PI*r), 5) + 6*M_PI*pow(r, 2)*pow(cos(M_PI*r), 7)*pow(cos(2*M_PI*z), 3) + 70*M_PI*pow(r, 2)*pow(cos(M_PI*r), 5)*pow(cos(2*M_PI*z), 2) + 200*M_PI*pow(r, 2)*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) - 12*M_PI*r*sin(M_PI*r)*pow(cos(M_PI*r), 6)*pow(cos(2*M_PI*z), 3) - 150*M_PI*r*sin(M_PI*r)*pow(cos(M_PI*r), 4)*pow(cos(2*M_PI*z), 2) - 3350*M_PI*r*sin(M_PI*r)*pow(cos(M_PI*r), 2)*cos(2*M_PI*z) - 14500*M_PI*r*sin(M_PI*r) - 12*M_PI*pow(cos(M_PI*r), 7)*pow(cos(2*M_PI*z), 3) + 3*pow(cos(M_PI*r), 7)*pow(cos(2*M_PI*z), 3) - 195*M_PI*pow(cos(M_PI*r), 5)*pow(cos(2*M_PI*z), 2) + 50*pow(cos(M_PI*r), 5)*pow(cos(2*M_PI*z), 2) - 3850*M_PI*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) + 1675*pow(cos(M_PI*r), 3)*cos(2*M_PI*z) - 22875*M_PI*cos(M_PI*r) + 14500*cos(M_PI*r))*sin(2*M_PI*z)*cos(M_PI*r);
      rhs(i, j, k, idx::URHO) += src_rho;
      rhs(i, j, k, idx::UMX)  += src_mr;
      rhs(i, j, k, idx::UMY)  += src_mz;
      rhs(i, j, k, idx::UET)  += src_E;
    });
  }
};

template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level) {
  amrex::ignore_unused(i, j, k, nt, tagfab, sdatafab, geomdata, pparm, level);
}

}  // namespace PROB

#endif
