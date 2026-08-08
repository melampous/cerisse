#ifndef CentralDif_H_
#define CentralDif_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

#include "diff_ops.H"
#include "IBMSharedGPFluxUtils.h"
template <bool isAD, bool isIB, int order, typename cls_t>
class centraldif_t {
  public:
  AMREX_GPU_HOST
  centraldif_t() {
    // initialize coefficients for flux interpolation based on order
    calc_CDcoeffs<order>(INTcoef,CDcoef);
  }

  AMREX_GPU_HOST_DEVICE
  ~centraldif_t() {}

  // vars accessed by functions 
  int order_sch=order;  
  int halfsten = order / 2;

  typedef Array1D<Real, 0, order> arrayNumCoef;
  arrayNumCoef CDcoef,INTcoef;

  
  //////////////////////////////////////////////////////
#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void inline eflux_ibm(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& cons, const cls_t* cls,
                    const Array4<uint8_t>& ibMarkers) {
#else
  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& cons, const cls_t* cls) {
#endif
                      
    // const Box& bx  = mfi.growntilebox(0);
    // const Box& bxg = mfi.growntilebox(cls->NGHOST);
    // ---------------------------------------------------------------------  //
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      auto const& flx = flxt[dir]->array(); 
      const Box bxface = mfi.grownnodaltilebox(dir, 0);

      // compute interface fluxes at i-1/2, j-1/2, k-1/2
      ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if (AMREX_USE_GPIBM || CNS_USE_EB)
                    // Make the capture explicit before the constexpr branch;
                    // NVCC rejects a variable first captured from inside an
                    // extended-lambda if-constexpr context.
                    const cls_t* local_cls = cls;
                    const auto local_prims = prims;
                    const IntVect iv(AMREX_D_DECL(i, j, k));
                    const IntVect ivd = IntVect::TheDimensionVector(dir);
                    if (ibm_flux::is_solid_solid_face(iv, ivd, ibMarkers)) {
                      ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
                      return;
                    }
                    const bool full_stencil = ibm_flux::stencil_all_fluid(
                        iv, ivd, -halfsten, order, ibMarkers);
                    if (!full_stencil) {
                      if constexpr (order >= 4) {
                        if (ibm_flux::one_sided_polynomial_flux(
                                iv, dir, local_prims, flx, ibMarkers,
                                *local_cls)) {
                          return;
                        }
                      }
                      const bool adjacent_valid =
                          ibm_flux::adjacent_states_valid<
                              decltype(prims), decltype(ibMarkers), cls_t>(
                              iv, ivd, local_prims, ibMarkers);
                      if (!adjacent_valid) {
                        ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
                        return;
                      }
                      // A symmetric two-point face flux is the matched
                      // second-order closure for all higher-order central
                      // variants when only one reconstructed GP layer exists.
                      ibm_flux::two_point_central_flux(
                          iv, dir, local_prims, flx, *local_cls);
                      return;
                    }
#endif
                    this->flux_dir(i, j, k,dir, vdir, cons, prims, flx, cls);
                  });
    }
  }  

  // compute flux in each direction at f[i-1/2]   stored in i,j,k
  // central formulation following f[i-1/2] = 1/2 (fi + fi-1)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int dir,const GpuArray<int, 3>& vdir, const Array4<Real>& /*cons*/, const Array4<Real>& prims, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    Real flux_l[cls_t::NCONS];
    
    // prepare flux arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   
        
    for (int l = 0; l < order; l++) {  

      IntVect iv(AMREX_D_DECL(il, jl, kl));

      // evaluate flux from primitive
      cls->prims2flux(iv,dir,prims,flux_l);

      // compute flux
      for (int n=0;n< cls_t::NCONS;n++){
        flx(i, j, k,n) += flux_l[n]*INTcoef(l);
      }  

      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];

    }

  }
  ////////////////////////////////////////////////////////////////////////////////////////


  };




#endif
