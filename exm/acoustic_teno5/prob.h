#ifndef CNS_PROB_H_
#define CNS_PROB_H_
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
using namespace amrex;
namespace PROB {
struct ProbParm {
  Real rho0=1.0, p0=1.0, gamma=1.4;
  Real A=1.0e-3;            // acoustic amplitude (linear)
  Real k=2.0*3.14159265358979*2.0;  // 2 waves in 16 cells = 8 PPW  // 1 wavelength in unit domain
};
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
typedef rhs_dt<weno_t<ReconScheme::Teno5,ProbClosures>, no_diffusive_t, no_source_t> ProbRHS;
void inline inputs() {
  ParmParse pp; pp.add("cns.order_rk",3); pp.add("cns.stages_rk",3);
}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i,int j,int k,Array4<Real> const &state,GeometryData const &geomdata,
              ProbClosures const &cls, ProbParm const &pp){
  const Real *plo=geomdata.ProbLo(); const Real *dx=geomdata.CellSize();
  Real x=plo[0]+(i+Real(0.5))*dx[0];
  Real c0=std::sqrt(pp.gamma*pp.p0/pp.rho0);
  Real s=std::sin(pp.k*x);
  // right-traveling linear acoustic wave: rho'/rho0 = u'/c0 = p'/(gamma p0)
  Real rho=pp.rho0*(1.0+pp.A*s);
  Real u  =pp.A*c0*s;
  Real p  =pp.p0*(1.0+pp.gamma*pp.A*s);
  state(i,j,k,cls.URHO)=rho;
  state(i,j,k,cls.UMX)=rho*u; state(i,j,k,cls.UMY)=0.0; state(i,j,k,cls.UMZ)=0.0;
  state(i,j,k,cls.UET)=p/(cls.gamma-1.0)+0.5*rho*u*u;
}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM],Real dr,const Real s_int[5],const Real s_refl[ProbClosures::NCONS],
         Real s_ext[5],const int idir,const int sgn,const Real time,GeometryData const&,
         ProbClosures const&,ProbParm const&){}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int i,int j,int k,const auto &state,const auto &rhs,const ProbParm &pp,
            ProbClosures const &cls,auto const dx){}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i,int j,int k,int nt,auto &tagfab,const auto &sdatafab,const auto &geomdata,
             const ProbParm &pp,int level){}
} // namespace PROB
#endif
