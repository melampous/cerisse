#ifndef IBM_WALLMODEL_H_
#define IBM_WALLMODEL_H_

#include <ibm_containers.h>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_Derive.H>

//--------------------------------------------------------------------------//
// \brief Templates for different wall types for IB
// param is a struct with the following options:
// \param Twall : wall temperature (required for isothermal)
// \param alpha : array of coefficients (required for generic bc)
// \param beta  : array of coefficients (required for generic bc)
//
// compute_surfIB calculates surface primitive variables:
//   normal and tangential velocities, P, T and mass fractions
//   as a function of x, y, z, normal and interpolated vars
//--------------------------------------------------------------------------//

/// \brief Copy and renormalise species mass fractions from the image point
///        to the ghost point (row 1 from row 2).
///
/// If the sum of mass fractions is positive, each species is divided by the
/// sum so that the array sums to 1. Otherwise the raw values are copied
/// (fallback for non-reacting or single-species runs).
template <typename cls_t, int eorder>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ibm_copy_species(Array2D<Real,0,eorder+1,0,cls_t::NPRIM-1>& q)
{
#if NUM_SPECIES > 1
    Real sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      sumY += q(2,cls_t::QFS+n);
    }
    if (sumY > 0.0_rt) {
      Real inv = 1.0_rt / sumY;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        q(1,cls_t::QFS+n) = q(2,cls_t::QFS+n) * inv;
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        q(1,cls_t::QFS+n) = q(2,cls_t::QFS+n);
      }
    }
#else
    amrex::ignore_unused(q);
#endif
}

/// \brief disIM-aware species copy: impose zero normal gradient on the mass
///        fractions to 2nd order (general spacing, surface + IP1 + IP2) and
///        renormalise to sum to 1.  Degrades to the 1st-order single-image-point
///        form (identical to the overload above) when EO < 2 or n_valid < 2.
template <typename cls_t, int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ibm_copy_species(QArr& q, const DisArr& disIM, int n_valid)
{
#if NUM_SPECIES > 1
    Real sumY = 0.0_rt;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      q(1,cls_t::QFS+n) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QFS+n, disIM, n_valid);
      sumY += q(1,cls_t::QFS+n);
    }
    if (sumY > 0.0_rt) {
      Real inv = 1.0_rt / sumY;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        q(1,cls_t::QFS+n) *= inv;
      }
    }
#else
    amrex::ignore_unused(q, disIM, n_valid);
#endif
}

//--------------------------------------------------------------------------//
// Isothermal slip wall
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_isothermal_slip_wall_t
{
public:
  static constexpr Real Twall = param::Twall;
  static constexpr int eorder_tparm = param::extrap_order;

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t1*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;                                                       // un  = 0 (no penetration)
    q(1,cls_t::QV) = ibm_zero_grad_surface<EO>(q, cls_t::QV,    disIM, n_valid);   // ut1 = slip (zero normal grad)
    q(1,cls_t::QW) = ibm_zero_grad_surface<EO>(q, cls_t::QW,    disIM, n_valid);   // ut2 = slip
    q(1,cls_t::QPRES) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QPRES, disIM, n_valid);// zero-gradient pressure
    q(1,cls_t::QT)    = param::Twall;                                              // prescribed wall temperature
    ibm_copy_species<cls_t, EO>(q, disIM, n_valid);
  }
};

//--------------------------------------------------------------------------//
// Isothermal no-slip wall
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_isothermal_noslip_wall_t
{
public:
  static constexpr Real Twall = param::Twall;
  static constexpr int eorder_tparm = param::extrap_order;

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t1*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;   // un  = 0
    q(1,cls_t::QV) = 0.0_rt;   // ut1 = 0
    q(1,cls_t::QW) = 0.0_rt;   // ut2 = 0
    q(1,cls_t::QPRES) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QPRES, disIM, n_valid); // zero-gradient pressure
    q(1,cls_t::QT)    = param::Twall;                                              // prescribed wall temperature
    ibm_copy_species<cls_t, EO>(q, disIM, n_valid);
  }
};

//--------------------------------------------------------------------------//
// Adiabatic slip wall
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_adiabatic_slip_wall_t
{
public:
  static constexpr int eorder_tparm = param::extrap_order;

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t1*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;                                                       // un  = 0
    q(1,cls_t::QV) = ibm_zero_grad_surface<EO>(q, cls_t::QV,    disIM, n_valid);   // ut1 = slip (zero normal grad)
    q(1,cls_t::QW) = ibm_zero_grad_surface<EO>(q, cls_t::QW,    disIM, n_valid);   // ut2 = slip
    q(1,cls_t::QPRES) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QPRES, disIM, n_valid);// zero-gradient pressure
    q(1,cls_t::QT)    = ibm_zero_grad_surface_pos<EO>(q, cls_t::QT,    disIM, n_valid);// zero-gradient temperature (adiabatic)
    ibm_copy_species<cls_t, EO>(q, disIM, n_valid);
  }
};

//--------------------------------------------------------------------------//
// Adiabatic no-slip wall
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_adiabatic_noslip_wall_t
{
public:
  static constexpr int eorder_tparm = param::extrap_order;

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t1*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;   // un  = 0
    q(1,cls_t::QV) = 0.0_rt;   // ut1 = 0
    q(1,cls_t::QW) = 0.0_rt;   // ut2 = 0
    q(1,cls_t::QPRES) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QPRES, disIM, n_valid); // zero-gradient pressure
    q(1,cls_t::QT)    = ibm_zero_grad_surface_pos<EO>(q, cls_t::QT,    disIM, n_valid); // zero-gradient temperature (adiabatic)
    ibm_copy_species<cls_t, EO>(q, disIM, n_valid);
  }
};

//--------------------------------------------------------------------------//
// General boundary condition
// Imposes BC of the form: phi(1) = alpha * phi_zg + beta
//   alpha=1, beta=0     → dphi/dn = 0  (Neumann)
//   alpha=0, beta=PHIBC → phi = PHIBC   (Dirichlet)
// where phi_zg is the zero-normal-gradient surface value reconstructed to 2nd
// order from IP1, IP2 (general spacing).  At EO=1 / n_valid<2 phi_zg = phi(2),
// so this reduces *exactly* to the previous "alpha*phi(2)+beta" form.
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_general_wall_t
{
public:
  static constexpr int eorder_tparm = param::extrap_order;

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t1*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    for (int n = 0; n <= cls_t::QLS; ++n) {
      const Real phi_zg = ibm_zero_grad_surface<EO>(q, n, disIM, n_valid);
      q(1,n) = param::alpha[n] * phi_zg + param::beta[n];
    }
  }
};

//--------------------------------------------------------------------------//
#endif

