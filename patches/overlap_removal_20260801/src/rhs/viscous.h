#ifndef VISCOUS_H_
#define VISCOUS_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <Constants.h>
#include <TransPele.h>
#include <LES.h>

#include <string>


#include "diff_ops.H"
#include "IBMSharedGPFluxUtils.h"

// param
//      :: order     spatial order of central derivatives
//      :: useLES    use LES by modifying viscosity  (default false)

template <typename param, typename cls_t>
class viscous_t {

  public:
  // trivial ctor/dtor
  AMREX_GPU_HOST_DEVICE
  constexpr viscous_t() = default;

  AMREX_GPU_HOST_DEVICE
  ~viscous_t() = default;

  // half stencil size 
  //int halfsten = param::order / 2;
  static constexpr int halfsten = param::order / 2;

#if NUM_SPECIES > 1
  typedef Array1D<Real, 0, param::order> arrayNumCoef;
  arrayNumCoef CDcoef,INTcoef;
  // Host-only initialization of arrays
  AMREX_GPU_HOST
  void init_coeffs()
  {
    calc_CDcoeffs<param::order>(INTcoef, CDcoef);
  }
#else
  // No-op init, so ProbRHS::init_coeffs() is always valid
  AMREX_GPU_HOST
  void init_coeffs() {}
#endif

  static bool rz_viscous_annular_average_enabled()
  {
    static const int enabled = [] {
      int value = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_viscous_annular_average", value);
      return value;
    }();
    return enabled != 0;
  }

#if (AMREX_SPACEDIM == 2)
  // A finite-volume R-Z state stores an annular average.  For the regular
  // radial velocity u_r = r w, recover the even quantity w rather than using
  // the annular average of u_r as a point value.  This makes the radial and
  // hoop stresses exactly compatible for u_r proportional to r at the axis.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real rz_annular_centroid(
      const amrex::IntVect& iv,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo) const noexcept {
    const Real rlo = prob_lo[0] + Real(iv[0]) * dx[0];
    const Real rhi = rlo + dx[0];
    const Real denominator = rhi * rhi - rlo * rlo;
    const Real tiny = Real(1.0e-14) * dx[0] * dx[0];
    if (amrex::Math::abs(denominator) <= tiny) {
      return Real(0.5) * (rlo + rhi);
    }
    return Real(2.0 / 3.0) *
           (rhi * rhi * rhi - rlo * rlo * rlo) / denominator;
  }

  template <typename T>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real rz_regular_w(
      const amrex::IntVect& iv,
      amrex::Array4<T> const& q,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo) const noexcept {
    const Real centroid = rz_annular_centroid(iv, dx, prob_lo);
    const Real tiny = Real(1.0e-14) * dx[0];
    return (amrex::Math::abs(centroid) > tiny)
        ? q(iv, cls_t::QU) / centroid
        : Real(0.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real rz_regular_ur_face_value(
      const amrex::IntVect& iv, const int face_dir,
      amrex::Array4<const Real> const& q,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo) const noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(face_dir);
    const Real wface = Real(0.5) *
        (rz_regular_w(iv, q, dx, prob_lo) +
         rz_regular_w(iv - ivn, q, dx, prob_lo));
    const Real rface = (face_dir == 0)
        ? prob_lo[0] + Real(iv[0]) * dx[0]
        : prob_lo[0] + (Real(iv[0]) + Real(0.5)) * dx[0];
    return rface * wface;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real rz_regular_ur_face_derivative(
      const amrex::IntVect& iv, const int face_dir, const int derivative_dir,
      amrex::Array4<const Real> const& q,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo) const noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(face_dir);
    const auto ivd = amrex::IntVect::TheDimensionVector(derivative_dir);
    const Real wface = Real(0.5) *
        (rz_regular_w(iv, q, dx, prob_lo) +
         rz_regular_w(iv - ivn, q, dx, prob_lo));

    Real dwdx;
    if (derivative_dir == face_dir) {
      dwdx = (rz_regular_w(iv, q, dx, prob_lo) -
              rz_regular_w(iv - ivn, q, dx, prob_lo)) * dxinv[derivative_dir];
    } else {
      dwdx = (rz_regular_w(iv + ivd, q, dx, prob_lo) -
              rz_regular_w(iv - ivd, q, dx, prob_lo) +
              rz_regular_w(iv - ivn + ivd, q, dx, prob_lo) -
              rz_regular_w(iv - ivn - ivd, q, dx, prob_lo)) *
             Real(0.25) * dxinv[derivative_dir];
    }

    const Real rface = (face_dir == 0)
        ? prob_lo[0] + Real(iv[0]) * dx[0]
        : prob_lo[0] + (Real(iv[0]) + Real(0.5)) * dx[0];
    return (derivative_dir == 0 ? wface : Real(0.0)) + rface * dwdx;
  }
#endif
   


#if (AMREX_USE_GPIBM || CNS_USE_EB )  
  void inline dflux_ibm(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,            
            const Array4<Real>& /*cons*/, const cls_t* cls,
            const Array4<uint8_t>& ibMarkers) {
#else
  void inline dflux(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt, 
            const Array4<Real>& /*cons*/, const cls_t* cls) {
#endif

    // LES options 
    constexpr bool useLES = []{
    if constexpr (requires { param::use_LES; })
        return param::use_LES;
    else
        return false;
    }();

    // mesh sizes
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    const bool is_rz = geom.IsRZ();
    const bool rz_annular_average =
        is_rz && param::order == 2 && rz_viscous_annular_average_enabled();
    if constexpr (useLES) {
      if (rz_annular_average) {
        amrex::Abort("cns.rz_viscous_annular_average is not yet compatible with LES");
      }
    }


    // grid
   // const Box& bx = mfi.tilebox();        
    const Box& bxg = mfi.growntilebox(cls->NGHOST);     // to handle high-order 

    // allocate arrays for transport properties  
    FArrayBox coeffs(bxg, cls_t::NCOEF, The_Async_Arena());
    const int CMU    = cls_t::CMU;
    const int CLAM   = cls_t::CLAM;
    const int CXI    = cls_t::CXI;
    const int CRHOD  = cls_t::CRHOD;
        
    const auto& mu_arr   = coeffs.array(CMU);     // dynamic viscosity
    const auto& lam_arr  = coeffs.array(CLAM);    // thermal conductivity 
    const auto& xi_arr   = coeffs.array(CXI);     // bulk viscosity
    const auto& rhoD_arr = coeffs.array(CRHOD);   // species diffusivity (times rho)

    // pointer to array of transport coefficients    
    const amrex::Array4<const amrex::Real>& coeftrans = coeffs.array();
    
    // calculate all transport properties and store in array (up to ghost points)
#ifdef USE_PELEPHYSICS

    FArrayBox qfab(bxg, cls_t::NPRIM, The_Async_Arena());     // prep space q-arrays
    auto const& q = qfab.array();
    // fill it with prims data
    amrex::ParallelFor( bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
      for (int n=0;n<cls_t::NPRIM; n++){ q(i,j,k,n) = prims(i,j,k,n);}   

      q(i,j,k,cls_t::QRHO) *=rho_si2cgs; // convert to cgs for PelePhysics   

    });  

    // pointers/arrays to communicate with pelephysics 
    auto const& q_y   = qfab.const_array(cls_t::QFS);             // species mass fraction
    auto const& q_T   = qfab.const_array(cls_t::QT);              // temperature 
    auto const& q_rho = qfab.const_array(cls_t::QRHO);            // density (change units below)
    
     
    BL_PROFILE("PelePhysics::get_transport_coeffs()");
    //Array4<Real> chi; // dummy Soret effect coef (not ready yet)
    // Soret effect (not used yet)
    const auto& chi_arr = coeffs.array(cls_t::CSORET); 
    
#if (PELEPVERSION==23)   
    trans_parms.allocate(); 
    auto const* ltransparm = trans_parms.device_trans_parm();
#else
    auto const* ltransparm = trans_parms.device_parm();
#endif    
    
    amrex::launch(bxg, [=] AMREX_GPU_DEVICE(Box const& tbx) {

            auto trans = pele::physics::PhysicsType::transport();                      
            trans.get_transport_coeffs(tbx, q_y, q_T, q_rho, 
                rhoD_arr, chi_arr, mu_arr,xi_arr, lam_arr, ltransparm);
          });

    // change units
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
        mu_arr(i,j,k) *= visc_cgs2si;
        lam_arr(i,j,k)*= cond_cgs2si;    
        for (int n=0;n<NUM_SPECIES; n++){        
          rhoD_arr(i,j,k,n) *= rhodiff_cgs2si;
        }   
        xi_arr(i,j,k) *= visc_cgs2si;     
        });        
    //    
#else
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
#if (AMREX_USE_GPIBM || CNS_USE_EB)
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        if (!ibm_flux::is_usable(iv, ibMarkers)) {
          mu_arr(i,j,k) = 0.0;
          lam_arr(i,j,k) = 0.0;
          xi_arr(i,j,k) = 0.0;
          return;
        }
#endif
        mu_arr(i,j,k)  = cls->visc(prims(i,j,k,cls_t::QT));
        lam_arr(i,j,k) = cls->cond(prims(i,j,k,cls_t::QT));       
        xi_arr(i,j,k)  = 0.0;
        });
#endif     

    // -------  LES Options  ----------- //
    if constexpr(useLES)
    {
      Real Delta = cls->calc_delta(dx); // compute filter width
      // loop over cells (including enough ghost to build stencil)
      const Box& bxgs = mfi.growntilebox(halfsten);
      // BEWARE cannot go over all the ghost cell !!
      // for viscous order 2, LES order can be  2,4
      // for viscous order 4, LES order can be  2
      // for viscous order 6, LES cannot be used
      amrex::ParallelFor( bxgs, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

        Real mu_sgs, cond_sgs, diff_sgs;   // per-cell SGS outputs (must be lambda-local, not captured-by-value)
        Real Cp_o_Pr = lam_arr(i,j,k)/(mu_arr(i,j,k)+1.e-15);
        // axisymmetric (r-z) hoop strain S_thetatheta = u_r/r passed to the SGS model
        Real hoop = Real(0.0);
        if (is_rz) {
          const Real rr = prob_lo[0] + (Real(i)+Real(0.5))*dx[0];
          hoop = (rr > Real(1.0e-12)) ? prims(i,j,k,cls_t::QU)/rr : Real(0.0);
        }
        cls-> compute_sgsterms(i,j,k,prims, dxinv, Delta, Cp_o_Pr,  mu_sgs, cond_sgs, diff_sgs, hoop, is_rz);
        mu_arr(i,j,k) += mu_sgs;
        lam_arr(i,j,k)+= cond_sgs;   
        for (int n=0;n<NUM_SPECIES; n++){        
          rhoD_arr(i,j,k,n) += diff_sgs;
        }   
      }); 
    }          

  
    // loop over directions -----------------------------------------------
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
     // GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
      auto const& flx = flxt[dir]->array(); 
      const Box bxface = mfi.grownnodaltilebox(dir, 0);

      // Yosihizawa model  tau_kk
      // if constexpr(useLES)
      // {
      //   Real Delta = cls->calc_delta(dx); // compute filter width
      //   amrex::ParallelFor(bxgnodal,
      //             [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                     
      //               flx(i,j,k,cls_t::UMX+dir) += cls->compute_xisgs(i,j,k,dir,prims, dxinv, Delta);
      //             });        
      // }

      // compute diffusion fluxes
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
      amrex::ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                                       
                    this->cns_diff_ibm(i, j, k,dir, prims,flx,coeftrans,
                                       dxinv, dx, prob_lo, is_rz,
                                       rz_annular_average, cls,ibMarkers
                                       );
                  });                      
#else
      amrex::ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->cns_diff(i, j, k,dir, prims,flx,coeftrans,
                                   dxinv, dx, prob_lo, is_rz,
                                   rz_annular_average, cls);
                  });
#endif

    }
    // end loop  ------------------------------------------------------

  }

#if !(AMREX_USE_GPIBM || CNS_USE_EB)
  // ----------------------------------------------------------------------------------------------
  // Region-parameterized overload (comm/comp overlap). Computes diffusion
  // fluxes only on the faces of surroundingNodes(rbx, dir) for each dir,
  // with transport coefficients evaluated on grow(rbx, halfsten) — exactly
  // the cells read by the face interpolations. Faces whose BOTH adjacent
  // cells lie inside the optional 'skip_cells' box are skipped (overlap
  // shell pass; those faces were computed in Pass 1), and coefficient cells
  // read only by skipped faces (grow(skip_cells, -halfsten)) are skipped
  // too. Per-face numerics are identical to the MFIter version above, so
  // splitting a tilebox into interior + shell regions yields
  // bitwise-identical fluxes at every face.
  // Requires prims valid on grow(rbx, max(halfsten (+1 for LES), NGHOST-safe)).
  void inline dflux(const Geometry& geom, const Box& rbx,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
            const Array4<Real>& /*cons*/, const cls_t* cls,
            const Box& skip_cells = Box()) {

    // LES options
    constexpr bool useLES = []{
    if constexpr (requires { param::use_LES; })
        return param::use_LES;
    else
        return false;
    }();

    // mesh sizes
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    const bool is_rz = geom.IsRZ();
    const bool rz_annular_average =
        is_rz && param::order == 2 && rz_viscous_annular_average_enabled();
    if constexpr (useLES) {
      if (rz_annular_average) {
        amrex::Abort("cns.rz_viscous_annular_average is not yet compatible with LES");
      }
    }

    // face-skip box and coefficient-skip box (cells read ONLY by skipped
    // faces; see proof in the header comment: a cell c in
    // grow(skip_cells,-halfsten) is only read by faces whose two adjacent
    // cells are both inside skip_cells, which are all skipped)
    const Box skipbox = skip_cells;
    const bool skip_ok = skipbox.ok();
    const Box cskipbox = skip_ok ? amrex::grow(skipbox, -halfsten) : Box();
    const bool cskip_ok = cskipbox.ok();

    // coefficient cells needed by the face interpolations of this region
    const Box bxg = amrex::grow(rbx, halfsten);

    // allocate arrays for transport properties
    FArrayBox coeffs(bxg, cls_t::NCOEF, The_Async_Arena());
    const int CMU    = cls_t::CMU;
    const int CLAM   = cls_t::CLAM;
    const int CXI    = cls_t::CXI;
    const int CRHOD  = cls_t::CRHOD;

    const auto& mu_arr   = coeffs.array(CMU);     // dynamic viscosity
    const auto& lam_arr  = coeffs.array(CLAM);    // thermal conductivity
    const auto& xi_arr   = coeffs.array(CXI);     // bulk viscosity
    const auto& rhoD_arr = coeffs.array(CRHOD);   // species diffusivity (times rho)

    // pointer to array of transport coefficients
    const amrex::Array4<const amrex::Real>& coeftrans = coeffs.array();

    // calculate all transport properties and store in array
#ifdef USE_PELEPHYSICS

    FArrayBox qfab(bxg, cls_t::NPRIM, The_Async_Arena());     // prep space q-arrays
    auto const& q = qfab.array();
    // fill it with prims data
    amrex::ParallelFor( bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      for (int n=0;n<cls_t::NPRIM; n++){ q(i,j,k,n) = prims(i,j,k,n);}

      q(i,j,k,cls_t::QRHO) *=rho_si2cgs; // convert to cgs for PelePhysics

    });

    // pointers/arrays to communicate with pelephysics
    auto const& q_y   = qfab.const_array(cls_t::QFS);             // species mass fraction
    auto const& q_T   = qfab.const_array(cls_t::QT);              // temperature
    auto const& q_rho = qfab.const_array(cls_t::QRHO);            // density (change units below)


    BL_PROFILE("PelePhysics::get_transport_coeffs()");
    // Soret effect (not used yet)
    const auto& chi_arr = coeffs.array(cls_t::CSORET);

#if (PELEPVERSION==23)
    trans_parms.allocate();
    auto const* ltransparm = trans_parms.device_trans_parm();
#else
    auto const* ltransparm = trans_parms.device_parm();
#endif

    amrex::launch(bxg, [=] AMREX_GPU_DEVICE(Box const& tbx) {

            auto trans = pele::physics::PhysicsType::transport();
            trans.get_transport_coeffs(tbx, q_y, q_T, q_rho,
                rhoD_arr, chi_arr, mu_arr,xi_arr, lam_arr, ltransparm);
          });

    // change units
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        mu_arr(i,j,k) *= visc_cgs2si;
        lam_arr(i,j,k)*= cond_cgs2si;
        for (int n=0;n<NUM_SPECIES; n++){
          rhoD_arr(i,j,k,n) *= rhodiff_cgs2si;
        }
        xi_arr(i,j,k) *= visc_cgs2si;
        });
    //
#else
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (cskip_ok && cskipbox.contains(i,j,k)) { return; }
        mu_arr(i,j,k)  = cls->visc(prims(i,j,k,cls_t::QT));
        lam_arr(i,j,k) = cls->cond(prims(i,j,k,cls_t::QT));
        xi_arr(i,j,k)  = 0.0;
        });
#endif

    // -------  LES Options  ----------- //
    if constexpr(useLES)
    {
      Real Delta = cls->calc_delta(dx); // compute filter width
      // SGS terms are added on the same cells the face interpolations read.
      // compute_sgsterms reads a +-1 stencil around each cell, so prims must
      // be valid on grow(bxg, 1) as well (guaranteed: halfsten+1 <= NGHOST).
      amrex::ParallelFor( bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (cskip_ok && cskipbox.contains(i,j,k)) { return; }

        Real mu_sgs, cond_sgs, diff_sgs;   // per-cell SGS outputs (must be lambda-local, not captured-by-value)
        Real Cp_o_Pr = lam_arr(i,j,k)/(mu_arr(i,j,k)+1.e-15);
        // axisymmetric (r-z) hoop strain S_thetatheta = u_r/r passed to the SGS model
        Real hoop = Real(0.0);
        if (is_rz) {
          const Real rr = prob_lo[0] + (Real(i)+Real(0.5))*dx[0];
          hoop = (rr > Real(1.0e-12)) ? prims(i,j,k,cls_t::QU)/rr : Real(0.0);
        }
        cls-> compute_sgsterms(i,j,k,prims, dxinv, Delta, Cp_o_Pr,  mu_sgs, cond_sgs, diff_sgs, hoop, is_rz);
        mu_arr(i,j,k) += mu_sgs;
        lam_arr(i,j,k)+= cond_sgs;
        for (int n=0;n<NUM_SPECIES; n++){
          rhoD_arr(i,j,k,n) += diff_sgs;
        }
      });
    }


    // loop over directions -----------------------------------------------
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      auto const& flx = flxt[dir]->array();

      const Box& rbxnodal = amrex::surroundingNodes(rbx, dir);
      const IntVect ivd = IntVect::TheDimensionVector(dir);

      // compute diffusion fluxes
      amrex::ParallelFor(rbxnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (skip_ok) {
                      const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
                      if (skipbox.contains(iv) && skipbox.contains(iv - ivd)) {
                        return;  // face computed in Pass 1
                      }
                    }
                    this->cns_diff(i, j, k,dir, prims,flx,coeftrans,
                                   dxinv, dx, prob_lo, is_rz,
                                   rz_annular_average, cls);
                  });

    }
    // end loop  ------------------------------------------------------
  }
#endif  // !(AMREX_USE_GPIBM || CNS_USE_EB)

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real velocity_normal_diff(
      const amrex::IntVect& iv, const int face_dir, const int comp,
      amrex::Array4<const Real> const& q,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
      const bool rz_annular_average) const noexcept {
#if (AMREX_SPACEDIM == 2)
    if (rz_annular_average && comp == cls_t::QU) {
      return rz_regular_ur_face_derivative(
          iv, face_dir, face_dir, q, dxinv, dx, prob_lo);
    }
#else
    amrex::ignore_unused(dx, prob_lo, rz_annular_average);
#endif
    return normal_diff<param::order>(iv, face_dir, comp, q, dxinv);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real velocity_tangent_diff(
      const amrex::IntVect& iv, const int face_dir, const int tangent_dir,
      const int comp, amrex::Array4<const Real> const& q,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
      const bool rz_annular_average) const noexcept {
#if (AMREX_SPACEDIM == 2)
    if (rz_annular_average && comp == cls_t::QU) {
      return rz_regular_ur_face_derivative(
          iv, face_dir, tangent_dir, q, dxinv, dx, prob_lo);
    }
#else
    amrex::ignore_unused(dx, prob_lo, rz_annular_average);
#endif
    return tangent_diff<param::order>(
        iv, face_dir, tangent_dir, comp, q, dxinv);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real velocity_face_interp(
      const amrex::IntVect& iv, const int face_dir, const int comp,
      amrex::Array4<const Real> const& q,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& dx,
      amrex::GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
      const bool rz_annular_average) const noexcept {
#if (AMREX_SPACEDIM == 2)
    if (rz_annular_average && comp == cls_t::QU) {
      return rz_regular_ur_face_value(iv, face_dir, q, dx, prob_lo);
    }
#else
    amrex::ignore_unused(dx, prob_lo, rz_annular_average);
#endif
    return interp<param::order>(iv, face_dir, comp, q);
  }


  // ----------------------------------------------------------------------------------------------
  /**
  * @brief Compute diffusion fluxes (viscosity + heat + diffusion).
  *        Calculates flux[i] which correspond to flux(i-1/2) between i and i-1
  *
  * @param i,j,k  x, y, z index.cls_t::CLAM
  * @param d1    direction, 0:x, 1:y, 2:z (dir)
  * @param q      primitive variables.
  * @param[out] flx  output diffusion fluxes.
  * @param coeffs transport coefficients.
  * @param dxinv  1/dx
  * @param cls_t  ProbClosures 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff(
      const int i, const int j, const int k, const int d1,
      amrex::Array4<const amrex::Real> const& q,
      amrex::Array4<amrex::Real> const& flx,
	      amrex::Array4<const amrex::Real> const& coeffs,
	      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
	      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
	      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& prob_lo,
	      const bool is_rz,
	      const bool rz_annular_average,
	      const cls_t* /*cls*/) const {
    
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(d1);

    const int d2 = d1 == 0 ? 1 : 0;
#if (AMREX_SPACEDIM == 3)    
    const int d3 = d1 == 2 ? 1 : 2;
#endif    

    AMREX_D_TERM(const int QU1 = cls_t::QU + d1;,  const int QU2 = cls_t::QU + d2;
               , const int QU3 = cls_t::QU + d3;)
    AMREX_D_TERM(const int UM1 = cls_t::UMX + d1;, const int UM2 = cls_t::UMX + d2;
               , const int UM3 = cls_t::UMX + d3;)

    // Aij = dA_i/dx_j
    const Real dTdn = normal_diff<param::order>(iv, d1, cls_t::QT, q, dxinv);
    const Real u11  = velocity_normal_diff(
        iv, d1, QU1, q, dxinv, dx, prob_lo, rz_annular_average);
#if (AMREX_SPACEDIM >= 2)
    const Real u21  = velocity_normal_diff(
        iv, d1, QU2, q, dxinv, dx, prob_lo, rz_annular_average);
    const Real u12  = velocity_tangent_diff(
        iv, d1, d2, QU1, q, dxinv, dx, prob_lo, rz_annular_average);
    const Real u22  = velocity_tangent_diff(
        iv, d1, d2, QU2, q, dxinv, dx, prob_lo, rz_annular_average);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real u31  = normal_diff<param::order>(iv, d1, QU3, q, dxinv);
    const Real u13  = tangent_diff<param::order>(iv, d1, d3, QU1, q, dxinv);
    const Real u33  = tangent_diff<param::order>(iv, d1, d3, QU3, q, dxinv);
#endif
	    Real divu = AMREX_D_TERM(u11, +u22, +u33);
#if (AMREX_SPACEDIM == 2)
	    if (is_rz) {
	      const Real r_face = (d1 == 0)
	          ? prob_lo[0] + Real(i) * dx[0]
	          : prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
	      const Real dudr_face = (d1 == 0) ? u11 : u22;
	      const Real ur_face = velocity_face_interp(
	          iv, d1, cls_t::QU, q, dx, prob_lo, rz_annular_average);
	      const Real tiny_r = Real(1.0e-14) * dx[0];
	      divu += (r_face > tiny_r) ? (ur_face / r_face) : dudr_face;
	    }
#endif
    
    const Real muf    = interp<param::order>(iv, d1, cls_t::CMU, coeffs);
    const Real xif    = interp<param::order>(iv, d1, cls_t::CXI, coeffs);    
    const Real lamf   = interp<param::order>(iv, d1, cls_t::CLAM, coeffs);
     
    AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - 2.0 / 3.0 * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

    // momentum
    AMREX_D_TERM(flx(iv, UM1) -= tau11;, flx(iv, UM2) -= tau12;, flx(iv, UM3) -= tau13;)

    // interpolate velocity
    // note: this is the velocity at the face, not at the cell center
    const Real u1    = velocity_face_interp(
        iv, d1, QU1, q, dx, prob_lo, rz_annular_average);
#if (AMREX_SPACEDIM >= 2)    
    const Real u2    = velocity_face_interp(
        iv, d1, QU2, q, dx, prob_lo, rz_annular_average);
#endif    
#if (AMREX_SPACEDIM == 3)
    const Real u3    = interp<param::order>(iv, d1, QU3, q);
#endif

    // energy heat flux
    flx(iv, cls_t::UET) -= AMREX_D_TERM(u1 * tau11, +u2 * tau12, +u3 * tau13) + lamf* dTdn;
    
#if NUM_SPECIES > 1    
    // --------------------------------------------------------------------------
    // diffusion species  
    Real ymass[param::order][NUM_SPECIES],xmole[param::order][NUM_SPECIES];
    Real hi[param::order][NUM_SPECIES];
    Real yaux[NUM_SPECIES],xaux[NUM_SPECIES],haux[NUM_SPECIES];

    auto thermo = typename cls_t::multispecies_pele_gas_t();
    
    // Declare clipping arrays unconditionally
    Real maxy[NUM_SPECIES], maxx[NUM_SPECIES], miny[NUM_SPECIES], minx[NUM_SPECIES];

    // Initialize only if needed
    if constexpr(param::order > 2)
    {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        maxy[n] = 0.0;
        maxx[n] = 0.0;
        miny[n] = 1.0;
        minx[n] = 1.0;
      }
    }

    amrex::IntVect ivp(iv -halfsten*amrex::IntVect::TheDimensionVector(d1));
    for (int l = 0; l < param::order; l++) {
      
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        yaux[n] = q(ivp,  cls_t::QFS + n); 
      }
      // calculate xmol and specific enthalpies
      thermo.Y2X(yaux, xaux);    
      thermo.RTY2Hi(q(ivp,  cls_t::QRHO), q(ivp,   cls_t::QT), yaux, haux);

      // store in temp array 
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        xmole[l][n] = xaux[n];
        ymass[l][n] = yaux[n]; 
        hi[l][n]    = haux[n];
      }
      //
      ivp +=  amrex::IntVect::TheDimensionVector(d1);  
      
      // calculate max and min to clip reconstruction (high-order)
      if constexpr(param::order > 2)
      {
        for (int n = 0; n < NUM_SPECIES; ++n) { 
          maxx[n] = max(maxx[n],xaux[n]);
          maxy[n] = max(maxy[n],yaux[n]);
          minx[n] = min(minx[n],xaux[n]);
          miny[n] = min(miny[n],yaux[n]); 
        }
      }

    }
    
    const Real dpdx  = normal_diff<param::order>(iv, d1, cls_t::QPRES, q, dxinv); 
    const Real pface = interp<param::order>(iv, d1, cls_t::QPRES, q);
    const Real dlnp = dpdx/pface; //Real dlnp = 0.0;  
     
    Real Vc = 0.0;    
    Real Yf[NUM_SPECIES],hf[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      Real Xface = 0.0, Yface = 0.0, hface = 0.0, dXdx = 0.0;
      for (int l = 0; l < param::order; l++) {
        Xface += xmole[l][n]*INTcoef(l);
        Yface += ymass[l][n]*INTcoef(l);
        hface += hi[l][n]*INTcoef(l);
        dXdx  += xmole[l][n]*CDcoef(l);
      }      
      dXdx  *= dxinv[d1];      
      // prevent extrema in high-order
      if constexpr(param::order > 2)
      {      
        Xface = min( max(Xface,minx[n]) ,maxx[n]);
        Yface = min( max(Yface,miny[n]) ,maxy[n]);
      }
      Yf[n] = Yface; hf[n] = hface;

      const Real rhoD_f = interp<param::order>(iv, d1, cls_t::CRHOD + n, coeffs); 
      const Real Vd = -rhoD_f * (dXdx + (Xface - Yface) * dlnp);
      Vc += Vd;
      flx(iv, cls_t::UFS + n) += Vd; 
      flx(iv, cls_t::UET)     += Vd * hface; 
     }
    // Add correction velocity to fluxes so sum(Vd) = 0
    for (int n = 0; n < NUM_SPECIES; ++n) {       
      flx(iv, cls_t::UFS + n)-= Yf[n] * Vc;
      flx(iv, cls_t::UET)    -= Yf[n] * hf[n] * Vc; 
    }

    // --------------------------------------------------------------------------
#endif                                           
      
  }
  // ---------------------------------------------------------------------------------------------
#if (AMREX_USE_GPIBM || CNS_USE_EB)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool ibm_diff_cell_usable(
      const amrex::IntVect& iv, const Array4<uint8_t>& marker) const noexcept {
    return (marker(iv, 0) == 0) || (marker(iv, 1) != 0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool ibm_diff_cell_fluid(
      const amrex::IntVect& iv, const Array4<uint8_t>& marker) const noexcept {
    return marker(iv, 0) == 0;
  }

  template <int order>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool ibm_standard_diff_stencil_clean(
      const amrex::IntVect& iv, int idir, const Array4<uint8_t>& marker) const noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);

    // Normal face derivative/interpolation footprint:
    // order 2: iv-1, iv; order 4: iv-2..iv+1; order 6: iv-3..iv+2.
    constexpr int half = order / 2;
    for (int dn = -half; dn <= half - 1; ++dn) {
      if (!ibm_diff_cell_fluid(iv + dn * ivn, marker)) return false;
    }

    // Tangential derivative footprint used by tangent_diff<order>. If any
    // point is a GP or unreconstructed solid cell, the formal high-order
    // Cartesian stencil is no longer clean relative to the immersed wall.
    for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
      if (tdir == idir) continue;
      const auto ivt = amrex::IntVect::TheDimensionVector(tdir);
      for (int dn = -half; dn <= half - 1; ++dn) {
        const amrex::IntVect base = iv + dn * ivn;
        for (int dt = 1; dt <= half; ++dt) {
          if (!ibm_diff_cell_fluid(base + dt * ivt, marker)) return false;
          if (!ibm_diff_cell_fluid(base - dt * ivt, marker)) return false;
        }
      }
    }
    return true;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real normal_diff_ibm2(
      const amrex::IntVect& iv, int idir, int comp,
      amrex::Array4<const amrex::Real> const& q,
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
      const Array4<uint8_t>& marker) const noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    const amrex::IntVect ivm = iv - ivn;
    if (ibm_diff_cell_usable(iv, marker) && ibm_diff_cell_usable(ivm, marker)) {
      return (q(iv, comp) - q(ivm, comp)) * dxinv[idir];
    }
    return Real(0.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real interp_ibm2(
      const amrex::IntVect& iv, int idir, int comp,
      amrex::Array4<const amrex::Real> const& q,
      const Array4<uint8_t>& marker) const noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    const amrex::IntVect ivm = iv - ivn;
    const bool okp = ibm_diff_cell_usable(iv, marker);
    const bool okm = ibm_diff_cell_usable(ivm, marker);
    if (okp && okm) return Real(0.5) * (q(iv, comp) + q(ivm, comp));
    if (okp) return q(iv, comp);
    if (okm) return q(ivm, comp);
    return Real(0.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real tangent_diff_cell_ibm2(
      const amrex::IntVect& iv, int tdir, int comp,
      amrex::Array4<const amrex::Real> const& q,
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
      const Array4<uint8_t>& marker,
      bool& second_order) const noexcept {
    const auto ivt = amrex::IntVect::TheDimensionVector(tdir);
    const amrex::IntVect ivp = iv + ivt;
    const amrex::IntVect ivm = iv - ivt;
    const amrex::IntVect ivpp = iv + 2 * ivt;
    const amrex::IntVect ivmm = iv - 2 * ivt;
    const bool okp = ibm_diff_cell_usable(ivp, marker);
    const bool ok0 = ibm_diff_cell_usable(iv, marker);
    const bool okm = ibm_diff_cell_usable(ivm, marker);
    const bool okpp = ibm_diff_cell_usable(ivpp, marker);
    const bool okmm = ibm_diff_cell_usable(ivmm, marker);

    second_order = false;
    if (!ok0) return Real(0.0);

    if (okp && okm) {
      second_order = true;
      return Real(0.5) * (q(ivp, comp) - q(ivm, comp)) * dxinv[tdir];
    }
    if (okp && okpp) {
      second_order = true;
      return Real(0.5) * (-Real(3.0) * q(iv, comp) +
                          Real(4.0) * q(ivp, comp) - q(ivpp, comp)) * dxinv[tdir];
    }
    if (okm && okmm) {
      second_order = true;
      return Real(0.5) * (Real(3.0) * q(iv, comp) -
                          Real(4.0) * q(ivm, comp) + q(ivmm, comp)) * dxinv[tdir];
    }

    // Non-contiguous three-point alternatives for sharp/grid-aligned corners.
    // These are the derivatives at x=0 of the quadratic through offsets
    // {0,+1,-2} and {0,-1,+2}, respectively.
    if (okp && okmm) {
      second_order = true;
      return (-Real(0.5) * q(iv, comp) + Real(2.0 / 3.0) * q(ivp, comp) -
              Real(1.0 / 6.0) * q(ivmm, comp)) * dxinv[tdir];
    }
    if (okm && okpp) {
      second_order = true;
      return (Real(0.5) * q(iv, comp) - Real(2.0 / 3.0) * q(ivm, comp) +
              Real(1.0 / 6.0) * q(ivpp, comp)) * dxinv[tdir];
    }

    // Geometry can leave only two usable points at an unresolved corner.
    // Keep a bounded first-order fallback; tangent_diff_ibm2 propagates the
    // quality flag so strict quality checks can reject this case explicitly.
    if (okp && ok0) return (q(ivp, comp) - q(iv, comp)) * dxinv[tdir];
    if (ok0 && okm) return (q(iv, comp) - q(ivm, comp)) * dxinv[tdir];
    return Real(0.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real tangent_diff_ibm2(
      const amrex::IntVect& iv, int idir, int tdir, int comp,
      amrex::Array4<const amrex::Real> const& q,
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
      const Array4<uint8_t>& marker,
      bool& second_order) const noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    bool second_order_p = false;
    bool second_order_m = false;
    const Real dp = tangent_diff_cell_ibm2(
        iv, tdir, comp, q, dxinv, marker, second_order_p);
    const Real dm = tangent_diff_cell_ibm2(
        iv - ivn, tdir, comp, q, dxinv, marker, second_order_m);
    second_order = second_order_p && second_order_m;
    return Real(0.5) * (dp + dm);
  }

#endif

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  // ---------------------------------------------------------------------------------------------
  /**
  * @brief Compute diffusion fluxes in IB/EB.
  *
  * @param i,j,k  x, y, z index.cls_t::CLAM
  * @param dir    direction, 0:x, 1:y, 2:z.
  * @param q      primitive variables.
  * @param[out] flx  output diffusion fluxes.
  * @param coeffs transport coefficients.
  * @param dxinv  1/dx
  * @param cls_t  ProbClosures 
  * @param marker  geometry markers  (sld and cutcells) 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_ibm(
      const int i, const int j, const int k, const int d1,
      amrex::Array4<const amrex::Real> const& q,
      amrex::Array4<amrex::Real> const& flx,
	      amrex::Array4<const amrex::Real> const& coeffs,
	      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
	      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
	      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& prob_lo,
	      const bool is_rz,
	      const bool rz_annular_average,
	      const cls_t* cls, const Array4<uint8_t>& marker
        ) const {
    
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
	    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(d1);

    const int d2 = d1 == 0 ? 1 : 0;
#if (AMREX_SPACEDIM == 3)
    const int d3 = d1 == 2 ? 1 : 2;
#endif
    AMREX_D_TERM(const int QU1 = cls_t::QU + d1;, const int QU2 = cls_t::QU + d2;
               , const int QU3 = cls_t::QU + d3;)
    AMREX_D_TERM(const int UM1 = cls_t::UMX + d1;, const int UM2 = cls_t::UMX + d2;
               , const int UM3 = cls_t::UMX + d3;)

    // ivm is iv -1.  Use the high-order viscous stencil only when its whole
    // footprint is pure fluid.  Otherwise, use IBM-aware second-order
    // operators that never read unreconstructed solid cells.
    const bool close_to_wall = !ibm_standard_diff_stencil_clean<param::order>(iv, d1, marker);
    const bool intersolid_flx = marker(iv,0) &&  marker(ivm,0);    // inter-flux

    if (intersolid_flx) return;  // flux =0  inside solid   

		    Real u11,dTdn,u21,u12,u22,muf,xif,lamf;
		    Real u1f,u2f;
#if (AMREX_SPACEDIM == 3)
		    Real u31,u13,u33,u3f;
#endif
#if (AMREX_SPACEDIM == 2)
	    Real ur_face = Real(0.0);
#endif
    
#if NUM_SPECIES > 1
    Real rhoD_f[NUM_SPECIES];
#endif
    // reduce interpolation and differentiation to second order across the wall
    if (close_to_wall)
    {        
      bool tangent_quality = ibm_diff_cell_usable(iv, marker) &&
                             ibm_diff_cell_usable(ivm, marker);
      bool deriv_quality = true;
      dTdn = normal_diff_ibm2(iv, d1, cls_t::QT, q, dxinv, marker);
      u11  = normal_diff_ibm2(iv, d1, QU1, q, dxinv, marker);
#if (AMREX_SPACEDIM >= 2)
      u21  = normal_diff_ibm2(iv, d1, QU2, q, dxinv, marker);
      u12  = tangent_diff_ibm2(iv, d1, d2, QU1, q, dxinv, marker, deriv_quality);
      tangent_quality = tangent_quality && deriv_quality;
      u22  = tangent_diff_ibm2(iv, d1, d2, QU2, q, dxinv, marker, deriv_quality);
      tangent_quality = tangent_quality && deriv_quality;
#endif
#if (AMREX_SPACEDIM == 3)
      u31  = normal_diff_ibm2(iv, d1, QU3, q, dxinv, marker);
      u13  = tangent_diff_ibm2(iv, d1, d3, QU1, q, dxinv, marker, deriv_quality);
      tangent_quality = tangent_quality && deriv_quality;
      u33  = tangent_diff_ibm2(iv, d1, d3, QU3, q, dxinv, marker, deriv_quality);
      tangent_quality = tangent_quality && deriv_quality;
#endif  
      amrex::ignore_unused(tangent_quality);
      // properties
      muf  = interp_ibm2(iv, d1, cls_t::CMU, coeffs, marker);
      xif  = interp_ibm2(iv, d1, cls_t::CXI, coeffs, marker);
	      lamf = interp_ibm2(iv, d1, cls_t::CLAM, coeffs, marker);
      u1f = interp_ibm2(iv, d1, QU1, q, marker);
#if (AMREX_SPACEDIM >= 2)
      u2f = interp_ibm2(iv, d1, QU2, q, marker);
#else
      u2f = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
      u3f = interp_ibm2(iv, d1, QU3, q, marker);
#endif
#if (AMREX_SPACEDIM == 2)
	      ur_face = interp_ibm2(iv, d1, cls_t::QU, q, marker);
#endif
#if NUM_SPECIES > 1          
      for (int n = 0; n < NUM_SPECIES; ++n) {  
        rhoD_f[n] = interp_ibm2(iv, d1, cls_t::CRHOD + n, coeffs, marker);
      }
#endif      
    }
    else
    {
      dTdn = normal_diff<param::order>(iv, d1, cls_t::QT, q, dxinv);
      u11  = velocity_normal_diff(
          iv, d1, QU1, q, dxinv, dx, prob_lo, rz_annular_average);
#if (AMREX_SPACEDIM >= 2)
      u21  = velocity_normal_diff(
          iv, d1, QU2, q, dxinv, dx, prob_lo, rz_annular_average);
      u12  = velocity_tangent_diff(
          iv, d1, d2, QU1, q, dxinv, dx, prob_lo, rz_annular_average);
      u22  = velocity_tangent_diff(
          iv, d1, d2, QU2, q, dxinv, dx, prob_lo, rz_annular_average);
#endif
#if (AMREX_SPACEDIM == 3)
      u31  = normal_diff<param::order>(iv, d1, QU3, q, dxinv);
      u13  = tangent_diff<param::order>(iv, d1, d3, QU1, q, dxinv);
      u33  = tangent_diff<param::order>(iv, d1, d3, QU3, q, dxinv);
#endif  
      // properties
      muf  = interp<param::order>(iv, d1, cls_t::CMU, coeffs);
      xif  = interp<param::order>(iv, d1, cls_t::CXI, coeffs);
	      lamf = interp<param::order>(iv, d1, cls_t::CLAM, coeffs);
      u1f = velocity_face_interp(
          iv, d1, QU1, q, dx, prob_lo, rz_annular_average);
#if (AMREX_SPACEDIM >= 2)
      u2f = velocity_face_interp(
          iv, d1, QU2, q, dx, prob_lo, rz_annular_average);
#else
      u2f = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
      u3f = interp<param::order>(iv, d1, QU3, q);
#endif
#if (AMREX_SPACEDIM == 2)
	      ur_face = velocity_face_interp(
	          iv, d1, cls_t::QU, q, dx, prob_lo, rz_annular_average);
#endif
#if NUM_SPECIES > 1          
      for (int n = 0; n < NUM_SPECIES; ++n) {  
        rhoD_f[n] = interp<param::order>(iv, d1, cls_t::CRHOD + n, coeffs);             
      }
#endif      

    }  
    
	    Real divu = AMREX_D_TERM(u11, +u22, +u33);
#if (AMREX_SPACEDIM == 2)
	    if (is_rz) {
	      const Real r_face = (d1 == 0)
	          ? prob_lo[0] + Real(i) * dx[0]
	          : prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
	      const Real dudr_face = (d1 == 0) ? u11 : u22;
	      const Real tiny_r = Real(1.0e-14) * dx[0];
	      divu += (r_face > tiny_r) ? (ur_face / r_face) : dudr_face;
	    }
#endif
    
    AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - (2.0 / 3.0) * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

    // momentum
    AMREX_D_TERM(flx(iv, UM1) -= tau11;, flx(iv, UM2) -= tau12;, flx(iv, UM3) -= tau13;)
   
    // Energy flux uses the same face interpolation order and IBM mask as the
    // stress/transport terms.  The old unconditional two-cell average both
    // capped pure-fluid energy flux at second order and bypassed the IBM mask.
    flx(iv, cls_t::UET) -= AMREX_D_TERM(u1f * tau11,
                                      +u2f * tau12,
                                      +u3f * tau13) + lamf * dTdn;
#if NUM_SPECIES > 1    
    // --------------------------------------------------------------------------
    // diffusion species  (this array should be order_local)
    Real ymass[param::order][NUM_SPECIES],xmole[param::order][NUM_SPECIES];
    Real hi[param::order][NUM_SPECIES];
    Real yaux[NUM_SPECIES],xaux[NUM_SPECIES],haux[NUM_SPECIES];

    auto thermo = typename cls_t::multispecies_pele_gas_t();


    amrex::IntVect ivp(iv -halfsten*amrex::IntVect::TheDimensionVector(d1));

    // Declare clipping arrays unconditionally
    Real maxy[NUM_SPECIES], maxx[NUM_SPECIES], miny[NUM_SPECIES], minx[NUM_SPECIES];

    // Initialize only if needed
    if constexpr(param::order > 2)
    {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        maxy[n] = 0.0;
        maxx[n] = 0.0;
        miny[n] = 1.0;
        minx[n] = 1.0;
      }
    }

    for (int l = 0; l < param::order; l++) {
 
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        yaux[n] = q(ivp,  cls_t::QFS + n); 
      }
      // calculate xmol and specific enthalpies
      thermo.Y2X(yaux, xaux);    
      thermo.RTY2Hi(q(ivp,  cls_t::QRHO), q(ivp,   cls_t::QT), yaux, haux);

      // store in temp array 
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        xmole[l][n] = xaux[n];
        ymass[l][n] = yaux[n]; 
        hi[l][n]    = haux[n];
      }

      ivp +=  amrex::IntVect::TheDimensionVector(d1);
       
      // calculate max and min to clip reconstruction (high-order)
      if constexpr(param::order > 2)
      {
        for (int n = 0; n < NUM_SPECIES; ++n) { 
          maxx[n] = max(maxx[n],xaux[n]);
          maxy[n] = max(maxy[n],yaux[n]);
          minx[n] = min(minx[n],xaux[n]);
          miny[n] = min(miny[n],yaux[n]); 
        }
      }

    } // end l arrays

    const Real dpdx  = normal_diff<param::order>(iv, d1, cls_t::QPRES, q, dxinv);     
    const Real pface = interp<param::order>(iv, d1, cls_t::QPRES, q);
    const Real dlnp = dpdx/pface; 
     
    Real Vc = 0.0;    
    Real Yf[NUM_SPECIES],hf[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      Real Xface = 0.0, Yface = 0.0, hface = 0.0, dXdx = 0.0;
      for (int l = 0; l < param::order; l++) {
        Xface += xmole[l][n]*INTcoef(l);
        Yface += ymass[l][n]*INTcoef(l);
        hface += hi[l][n]*INTcoef(l);
        dXdx  += xmole[l][n]*CDcoef(l);
      }      
      dXdx  *= dxinv[d1];  //<<<<<<<<<<<-  BUG Fixed (should be dXdx*dxinv)     

      // prevent extrema in high-order
      if constexpr(param::order > 2)
      {      
        Xface = min( max(Xface,minx[n]) ,maxx[n]);
        Yface = min( max(Yface,miny[n]) ,maxy[n]);
      }
      
      Yf[n] = Yface; hf[n] = hface;
    
      const Real Vd = -rhoD_f[n] * (dXdx + (Xface - Yface) * dlnp);
      Vc += Vd;
      flx(iv, cls_t::UFS + n) += Vd; 
      flx(iv, cls_t::UET)     += Vd * hface;    
     }

    // Add correction velocity to fluxes so sum(Vd) = 0
    for (int n = 0; n < NUM_SPECIES; ++n) {       
      flx(iv, cls_t::UFS + n)-= Yf[n] * Vc;
      flx(iv, cls_t::UET)    -= Yf[n] * hf[n] * Vc;
    } 

  // --------------------------------------------------------------------------
#endif 


  }
#endif


	  // RZ viscous hoop-stress source for radial momentum:
	  //   RHS(rho u_r) += -tau_theta_theta / r,
	  //   tau_theta_theta = 2 mu u_r/r - 2/3 mu Theta + xi Theta,
	  //   Theta = du_r/dr + u_r/r + du_z/dz.
	  // The remaining viscous terms are already handled by the metric FV
	  // divergence of the face fluxes.
	  void inline rz_geometric_source(const amrex::Geometry& geom,
	                                  const amrex::MFIter& mfi,
	                                  const amrex::Array4<amrex::Real>& prims,
	                                  const amrex::Array4<amrex::Real>& state,
	                                  const cls_t* cls) {
#if (AMREX_SPACEDIM == 2)
	    if (!geom.IsRZ()) return;

	    // Validation toggle (default ON): cns.rz_visc_hoop=0 disables the hoop
	    // source for A/B testing its effect on near-axis vorticity. Default
	    // preserves the as-implemented behaviour.
	    static const int s_hoop = []{ int v = 1;
	        amrex::ParmParse pp("cns"); pp.query("rz_visc_hoop", v); return v; }();
	    if (!s_hoop) return;

	    const auto dx = geom.CellSizeArray();
	    const auto prob_lo = geom.ProbLoArray();
	    const auto dxinv = geom.InvCellSizeArray();
	    const Box& bx = mfi.tilebox();
	    const bool rz_annular_average =
	        param::order == 2 && rz_viscous_annular_average_enabled();

	    // Rebuild the cell transport coefficients used by the face stresses.
	    // In particular, PelePhysics does not expose its mixture coefficients
	    // through cls->visc(), and LES adds mu_sgs after the molecular lookup.
	    // The hoop source must use that same effective mu and bulk viscosity xi.
	    FArrayBox coeffs(bx, cls_t::NCOEF, The_Async_Arena());
	    const auto& mu_arr = coeffs.array(cls_t::CMU);
	    const auto& xi_arr = coeffs.array(cls_t::CXI);
	    const auto& lam_arr = coeffs.array(cls_t::CLAM);
#ifdef USE_PELEPHYSICS
	    const auto& rhoD_arr = coeffs.array(cls_t::CRHOD);
	    const auto& chi_arr = coeffs.array(cls_t::CSORET);
	    FArrayBox qfab(bx, cls_t::NPRIM, The_Async_Arena());
	    const auto& q = qfab.array();
	    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      for (int n = 0; n < cls_t::NPRIM; ++n) {
	        q(i,j,k,n) = prims(i,j,k,n);
	      }
	      q(i,j,k,cls_t::QRHO) *= rho_si2cgs;
	    });
#if (PELEPVERSION==23)
	    trans_parms.allocate();
	    auto const* ltransparm = trans_parms.device_trans_parm();
#else
	    auto const* ltransparm = trans_parms.device_parm();
#endif
	    const auto& q_y = qfab.const_array(cls_t::QFS);
	    const auto& q_T = qfab.const_array(cls_t::QT);
	    const auto& q_rho = qfab.const_array(cls_t::QRHO);
	    amrex::launch(bx, [=] AMREX_GPU_DEVICE(Box const& tbx) {
	      auto trans = pele::physics::PhysicsType::transport();
	      trans.get_transport_coeffs(tbx, q_y, q_T, q_rho, rhoD_arr,
	                                 chi_arr, mu_arr, xi_arr, lam_arr,
	                                 ltransparm);
	    });
	    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      mu_arr(i,j,k) *= visc_cgs2si;
	      xi_arr(i,j,k) *= visc_cgs2si;
	      lam_arr(i,j,k) *= cond_cgs2si;
	    });
#else
	    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      Real temperature = prims(i,j,k,cls_t::QT);
	      mu_arr(i,j,k) = cls->visc(temperature);
	      xi_arr(i,j,k) = Real(0.0);
	      lam_arr(i,j,k) = cls->cond(temperature);
	    });
#endif

	    constexpr bool useLES = [] {
	      if constexpr (requires { param::use_LES; }) {
	        return bool(param::use_LES);
	      }
	      return false;
	    }();
	    if constexpr (useLES) {
	      if (rz_annular_average) {
	        amrex::Abort("cns.rz_viscous_annular_average is not yet compatible with LES");
	      }
	    }
	    if constexpr (useLES) {
	      const Real Delta = cls->calc_delta(dx);
	      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	        Real mu_sgs = Real(0.0);
	        Real cond_sgs = Real(0.0);
	        Real diff_sgs = Real(0.0);
	        const Real Cp_o_Pr =
	            lam_arr(i,j,k) / (mu_arr(i,j,k) + Real(1.0e-15));
	        const Real r = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
	        const Real hoop = (r > Real(1.0e-12))
	            ? prims(i,j,k,cls_t::QU) / r
	            : Real(0.0);
	        cls->compute_sgsterms(i, j, k, prims, dxinv, Delta, Cp_o_Pr,
	                              mu_sgs, cond_sgs, diff_sgs, hoop, true);
	        mu_arr(i,j,k) += mu_sgs;
	      });
	    }

	    amrex::ParallelFor(bx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      const Real r = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
	      const Real tiny_r = Real(1.0e-14) * dx[0];
	      const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
	      Real ur_over_r;
	      Real durdr;
	      if (rz_annular_average) {
	        const auto ivr = amrex::IntVect::TheDimensionVector(0);
	        const Real w = rz_regular_w(iv, prims, dx, prob_lo);
	        const Real dwdr = Real(0.5) *
	            (rz_regular_w(iv + ivr, prims, dx, prob_lo) -
	             rz_regular_w(iv - ivr, prims, dx, prob_lo)) * dxinv[0];
	        ur_over_r = w;
	        durdr = w + r * dwdr;
	      } else {
	        const Real ur = prims(i,j,k,cls_t::QU);
	        ur_over_r = (r > tiny_r)
	            ? (ur / r)
	            : normal_diff_cc<param::order>(iv, 0, cls_t::QU, prims, dxinv);
	        durdr = normal_diff_cc<param::order>(
	            iv, 0, cls_t::QU, prims, dxinv);
	      }
	      const Real duzdz = normal_diff_cc<param::order>(iv, 1, cls_t::QV, prims, dxinv);
	      const Real theta = durdr + ur_over_r + duzdz;

	      const Real mu = mu_arr(i,j,k);
	      const Real xi = xi_arr(i,j,k);
	      const Real tau_tt = mu * (Real(2.0) * ur_over_r - Real(2.0/3.0) * theta)
	                        + xi * theta;
	      state(i,j,k,cls_t::UMX) -= tau_tt / amrex::max(r, tiny_r);
	    });
#else
	    amrex::ignore_unused(geom, mfi, prims, state, cls);
#endif
	  }

  };

//---------------------------------------------
#endif
