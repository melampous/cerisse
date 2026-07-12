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

template <typename param>
AMREX_GPU_HOST_DEVICE constexpr ibm_pressure_closure_t
ibm_pressure_closure_value()
{
  if constexpr (requires { param::pressure_closure; }) {
    return param::pressure_closure;
  }
  return ibm_pressure_closure_t::zero_gradient;
}

template <typename param>
AMREX_GPU_HOST_DEVICE constexpr Real ibm_pressure_sensor_low()
{
  if constexpr (requires { param::pressure_sensor_low; }) {
    return Real(param::pressure_sensor_low);
  }
  return Real(0.03);
}

template <typename param>
AMREX_GPU_HOST_DEVICE constexpr Real ibm_pressure_sensor_high()
{
  if constexpr (requires { param::pressure_sensor_high; }) {
    return Real(param::pressure_sensor_high);
  }
  return Real(0.10);
}

/// Select the pressure closure for the built-in walls.  Existing parameter
/// structs that do not declare pressure_closure retain the historical
/// zero-normal-gradient condition.  fluid_extrapolation leaves pressure to the
/// fluid-side solution; shock_aware_fluid_extrapolation limits that extrapolate
/// when the first two image points span an irreversible compression;
/// prescribed_gradient accepts a problem-supplied dp/dn; and
/// euler_slip_analytic_curvature evaluates the stationary Euler compatibility
/// condition
///
///   dp/dn = rho * u_t^T (grad n) u_t
///
/// from an analytic shape operator supplied by the problem.  The latter is
/// deliberately restricted to slip walls: selecting it for a no-slip wall is
/// a compile-time error, and no curvature is inferred from STL facets.
template <typename param, typename cls_t, int EO, bool IsSlip,
          typename XYZArr, typename NormArr, typename TanArr,
          typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_wall_pressure_surface(
    const XYZArr& xyz, const NormArr& norm,
    const TanArr& tangent1, const TanArr& tangent2, const QArr& q,
    const DisArr& disIM, int n_valid)
{
    constexpr auto closure = ibm_pressure_closure_value<param>();
    if constexpr (closure == ibm_pressure_closure_t::zero_gradient) {
      amrex::ignore_unused(xyz, norm, tangent1, tangent2);
      return ibm_zero_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid);
    } else if constexpr (closure == ibm_pressure_closure_t::fluid_extrapolation) {
      amrex::ignore_unused(xyz, norm, tangent1, tangent2);
      return ibm_fluid_extrap_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid);
    } else if constexpr (
        closure == ibm_pressure_closure_t::shock_aware_fluid_extrapolation) {
      amrex::ignore_unused(xyz, norm, tangent1, tangent2);
      constexpr Real sensor_low = ibm_pressure_sensor_low<param>();
      constexpr Real sensor_high = ibm_pressure_sensor_high<param>();
      static_assert(sensor_low >= Real(0.0) && sensor_low < sensor_high &&
                    sensor_high <= Real(2.0),
                    "IBM pressure sensor thresholds require "
                    "0 <= low < high <= 2");
      return ibm_shock_aware_pressure_surface<EO>(
          q, cls_t::QPRES, cls_t::QRHO, cls_t::QC,
          disIM, n_valid, sensor_low, sensor_high).value;
    } else if constexpr (closure == ibm_pressure_closure_t::prescribed_gradient) {
      amrex::ignore_unused(tangent1, tangent2);
      static_assert(
          requires { param::pressure_normal_derivative(xyz, norm); },
          "prescribed IBM pressure closure requires "
          "param::pressure_normal_derivative(xyz, normal)");
      const Real dpdn = param::pressure_normal_derivative(xyz, norm);
      return ibm_normal_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid, dpdn);
    } else if constexpr (
        closure == ibm_pressure_closure_t::euler_slip_analytic_curvature) {
      static_assert(
          IsSlip,
          "euler_slip_analytic_curvature is valid only for an IBM slip wall");
      static_assert(
          requires {
            param::wall_shape_operator(xyz, norm, tangent1, tangent2);
          },
          "euler_slip_analytic_curvature requires "
          "param::wall_shape_operator(xyz, normal, tangent1, tangent2)");

      const auto shape =
          param::wall_shape_operator(xyz, norm, tangent1, tangent2);
      const bool shape_valid =
          shape.valid != 0 &&
          amrex::Math::isfinite(shape.k11) &&
          amrex::Math::isfinite(shape.k12) &&
          amrex::Math::isfinite(shape.k22);
      if (!shape_valid) {
        return ibm_fluid_extrap_surface_pos<EO>(
            q, cls_t::QPRES, disIM, n_valid);
      }

      // Density is not an independently prescribed Euler wall variable.  A
      // fluid-side extrapolate supplies rho_s without coupling this explicit
      // compatibility condition back through p_s/T_s.  Its O(h^2) error is
      // multiplied by the O(h) inverse derivative coefficient, so it does not
      // reduce the second-order wall pressure reconstruction.
      const Real rho_s = ibm_fluid_extrap_surface_pos<EO>(
          q, cls_t::QRHO, disIM, n_valid);
      const Real ut1 = q(1, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
      const Real ut2 = q(1, cls_t::QW);
#else
      const Real ut2 = Real(0.0);
#endif
      const Real curvature_acceleration =
          shape.k11 * ut1 * ut1 + Real(2.0) * shape.k12 * ut1 * ut2 +
          shape.k22 * ut2 * ut2;
      const Real dpdn = rho_s * curvature_acceleration;
      if (!(rho_s > Real(0.0)) || !amrex::Math::isfinite(rho_s) ||
          !amrex::Math::isfinite(dpdn)) {
        return ibm_fluid_extrap_surface_pos<EO>(
            q, cls_t::QPRES, disIM, n_valid);
      }
      return ibm_normal_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid, dpdn);
    } else {
      static_assert(
          closure == ibm_pressure_closure_t::zero_gradient,
          "unknown IBM pressure closure");
      amrex::ignore_unused(xyz, norm, tangent1, tangent2);
      return ibm_zero_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid);
    }
}

//--------------------------------------------------------------------------//
// Isothermal slip wall
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_isothermal_slip_wall_t
{
public:
  using ibm_param_t = param;
  static constexpr Real Twall = param::Twall;
  static constexpr int eorder_tparm = param::extrap_order;
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_value<param>();
  static constexpr Real pressure_sensor_low = ibm_pressure_sensor_low<param>();
  static constexpr Real pressure_sensor_high = ibm_pressure_sensor_high<param>();

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t2,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;                                                       // un  = 0 (no penetration)
    q(1,cls_t::QV) = ibm_zero_grad_surface<EO>(q, cls_t::QV,    disIM, n_valid);   // ut1 = slip (zero normal grad)
    q(1,cls_t::QW) = ibm_zero_grad_surface<EO>(q, cls_t::QW,    disIM, n_valid);   // ut2 = slip
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, true>(
        xyz, norm, t1, t2, q, disIM, n_valid);
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
  using ibm_param_t = param;
  static constexpr Real Twall = param::Twall;
  static constexpr int eorder_tparm = param::extrap_order;
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_value<param>();
  static constexpr Real pressure_sensor_low = ibm_pressure_sensor_low<param>();
  static constexpr Real pressure_sensor_high = ibm_pressure_sensor_high<param>();

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t2,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;   // un  = 0
    q(1,cls_t::QV) = 0.0_rt;   // ut1 = 0
    q(1,cls_t::QW) = 0.0_rt;   // ut2 = 0
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, false>(
        xyz, norm, t1, t2, q, disIM, n_valid);
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
  using ibm_param_t = param;
  static constexpr int eorder_tparm = param::extrap_order;
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_value<param>();
  static constexpr Real pressure_sensor_low = ibm_pressure_sensor_low<param>();
  static constexpr Real pressure_sensor_high = ibm_pressure_sensor_high<param>();

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t2,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;                                                       // un  = 0
    q(1,cls_t::QV) = ibm_zero_grad_surface<EO>(q, cls_t::QV,    disIM, n_valid);   // ut1 = slip (zero normal grad)
    q(1,cls_t::QW) = ibm_zero_grad_surface<EO>(q, cls_t::QW,    disIM, n_valid);   // ut2 = slip
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, true>(
        xyz, norm, t1, t2, q, disIM, n_valid);
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
  using ibm_param_t = param;
  static constexpr int eorder_tparm = param::extrap_order;
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_value<param>();
  static constexpr Real pressure_sensor_low = ibm_pressure_sensor_low<param>();
  static constexpr Real pressure_sensor_high = ibm_pressure_sensor_high<param>();

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int,EO> /*eo_tag*/,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& t2,
    const Array1D<Real,0,EO-1>& disIM,
    int n_valid,
    Array2D<Real,0,EO+1,0,cls_t::NPRIM-1>& q,
    int /*type_solid_bc*/, const cls_t* /*cls*/)
  {
    q(1,cls_t::QU) = 0.0_rt;   // un  = 0
    q(1,cls_t::QV) = 0.0_rt;   // ut1 = 0
    q(1,cls_t::QW) = 0.0_rt;   // ut2 = 0
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, false>(
        xyz, norm, t1, t2, q, disIM, n_valid);
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
