#ifndef IBM_WALLMODEL_H_
#define IBM_WALLMODEL_H_

#include <ibm_containers.h>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_Derive.H>
#include <CNSconstants.h>

#include <cmath>

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
/// a compile-time error, and no curvature is inferred from STL facets.  The
/// Navier--Stokes no-slip policy consumes a solver-side Cartesian WLS estimate
/// of dp/dn = n dot div(tau); an invalid estimate falls back to fluid-side
/// extrapolation.
template <typename param, typename cls_t, int EO, bool IsSlip,
          typename XYZArr, typename NormArr, typename TanArr,
          typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_wall_pressure_surface(
    const XYZArr& xyz, const NormArr& norm,
    const TanArr& tangent1, const TanArr& tangent2, const QArr& q,
    const DisArr& disIM, int n_valid,
    const ibm_pressure_compatibility_t& pressure_compat)
{
    constexpr auto closure = ibm_pressure_closure_value<param>();
    if constexpr (closure == ibm_pressure_closure_t::zero_gradient) {
      amrex::ignore_unused(xyz, norm, tangent1, tangent2, pressure_compat);
      return ibm_zero_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid);
    } else if constexpr (closure == ibm_pressure_closure_t::fluid_extrapolation) {
      amrex::ignore_unused(xyz, norm, tangent1, tangent2, pressure_compat);
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
          disIM, n_valid, sensor_low, sensor_high,
          pressure_compat.feature_pressure_fallback).value;
    } else if constexpr (closure == ibm_pressure_closure_t::prescribed_gradient) {
      amrex::ignore_unused(tangent1, tangent2, pressure_compat);
      static_assert(
          requires { param::pressure_normal_derivative(xyz, norm); },
          "prescribed IBM pressure closure requires "
          "param::pressure_normal_derivative(xyz, normal)");
      const Real dpdn = param::pressure_normal_derivative(xyz, norm);
      return ibm_normal_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid, dpdn);
    } else if constexpr (
        closure ==
        ibm_pressure_closure_t::navier_stokes_noslip_viscous_compatibility) {
      static_assert(
          !IsSlip,
          "navier_stokes_noslip_viscous_compatibility requires an IBM no-slip wall");
      amrex::ignore_unused(xyz, norm, tangent1, tangent2);
      if (pressure_compat.valid == 0 ||
          !amrex::Math::isfinite(pressure_compat.dpdn)) {
        return ibm_fluid_extrap_surface_pos<EO>(
            q, cls_t::QPRES, disIM, n_valid);
      }
      const Real candidate = ibm_normal_grad_surface<EO>(
          q, cls_t::QPRES, disIM, n_valid, pressure_compat.dpdn);
      return (candidate > Real(0.0) && amrex::Math::isfinite(candidate))
                 ? candidate
                 : ibm_fluid_extrap_surface_pos<EO>(
                       q, cls_t::QPRES, disIM, n_valid);
    } else if constexpr (
        closure == ibm_pressure_closure_t::euler_slip_analytic_curvature) {
      amrex::ignore_unused(pressure_compat);
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
      amrex::ignore_unused(xyz, norm, tangent1, tangent2, pressure_compat);
      return ibm_zero_grad_surface_pos<EO>(
          q, cls_t::QPRES, disIM, n_valid);
    }
}

/// Reconstruct the ideal-gas entropy proxy
///
///   sigma = log(p) - gamma log(rho) = log(p / rho^gamma)
///
/// from the fluid-side image states, then recover rho and T from the
/// independently reconstructed pressure. The legacy path uses a homogeneous
/// auxiliary normal extension. The production one-sided-jet path instead
/// extrapolates sigma to the BI, because Euler prescribes no d_n(sigma)=0 wall
/// condition and a smooth solution may carry a non-zero normal entropy
/// gradient.
template <typename cls_t, int EO, bool OneSidedEntropy = false,
          typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_entropy_surface_candidate(
    const QArr& q, const DisArr& disIM, int n_valid, const cls_t* cls,
    Real pressure, Real& density, Real& temperature)
{
  static_assert(
      requires(const cls_t* closure) {
        closure->gamma;
        closure->Rspec;
      },
      "Entropy-based IBM wall extension requires ideal-gas gamma and Rspec");

  Array2D<Real, 0, EO + 1, 0, 0> entropy_proxy{};
  bool support_valid = n_valid >= 1 && cls->gamma > Real(1.0) &&
                       cls->Rspec > Real(0.0);
  for (int image = 0; image < EO; ++image) {
    if (image >= n_valid) break;
    const Real pressure = q(2 + image, cls_t::QPRES);
    const Real density = q(2 + image, cls_t::QRHO);
    const bool state_valid =
        pressure > CNSConstants::min_press() && density > Real(0.0) &&
        amrex::Math::isfinite(pressure + density);
    support_valid = support_valid && state_valid;
    if (state_valid) {
      entropy_proxy(2 + image, 0) =
          std::log(pressure) - cls->gamma * std::log(density);
    }
  }

  density = Real(-1.0);
  temperature = Real(-1.0);
  if (support_valid && pressure > CNSConstants::min_press()) {
    const Real sigma = OneSidedEntropy
        ? ibm_fluid_extrap_surface<EO>(
              entropy_proxy, 0, disIM, n_valid)
        : ibm_zero_grad_surface<EO>(
              entropy_proxy, 0, disIM, n_valid);
    density =
        std::exp((std::log(pressure) - sigma) / cls->gamma);
    temperature = pressure / (density * cls->Rspec);
    const bool candidate_valid =
        pressure > CNSConstants::min_press() && density > Real(0.0) &&
        temperature > Real(0.0) &&
        amrex::Math::isfinite(sigma + density + temperature);
    return candidate_valid;
  }
  return false;
}

template <typename cls_t, int EO, bool OneSidedEntropy = false,
          typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_entropy_surface_state(
    QArr& q, const DisArr& disIM, int n_valid, const cls_t* cls)
{
  const Real pressure = q(1, cls_t::QPRES);
  Real density = Real(-1.0);
  Real temperature = Real(-1.0);
  if (ibm_entropy_surface_candidate<cls_t, EO, OneSidedEntropy>(
          q, disIM, n_valid, cls, pressure, density, temperature)) {
    q(1, cls_t::QRHO) = density;
    q(1, cls_t::QT) = temperature;
    return true;
  }

  // Complete thermodynamic fallback to the nearest fluid image state. The
  // pressure closure is reduced together with rho and T, avoiding a mixed
  // state assembled from incompatible reconstruction orders.
  q(1, cls_t::QPRES) = q(2, cls_t::QPRES);
  q(1, cls_t::QRHO) = q(2, cls_t::QRHO);
  q(1, cls_t::QT) = q(2, cls_t::QT);
  return false;
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
  static constexpr bool stationary_no_slip = false;
  static constexpr bool stationary_slip = true;
  static constexpr bool isothermal_wall = true;
  static constexpr bool adiabatic_wall = false;
  static constexpr bool entropy_extension = false;
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
    int /*type_solid_bc*/, const cls_t* /*cls*/,
    const ibm_pressure_compatibility_t& pressure_compat)
  {
    q(1,cls_t::QU) = 0.0_rt;                                                       // un  = 0 (no penetration)
    q(1,cls_t::QV) = ibm_zero_grad_surface<EO>(q, cls_t::QV, disIM, n_valid);       // ut1 = slip
    q(1,cls_t::QW) = ibm_zero_grad_surface<EO>(q, cls_t::QW, disIM, n_valid);       // ut2 = slip
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, true>(
        xyz, norm, t1, t2, q, disIM, n_valid, pressure_compat);
    q(1,cls_t::QT)    = param::Twall;                                              // prescribed wall temperature
    ibm_copy_species<cls_t, EO>(q, disIM, n_valid);
  }
};

//--------------------------------------------------------------------------//
// Stationary Euler slip wall
//
// The wall imposes no penetration and the selected pressure compatibility
// relation. Tangential velocity and the entropy proxy receive homogeneous
// auxiliary normal extensions. Unlike the isothermal/adiabatic wall types,
// this class does not impose a viscous thermal boundary condition on Euler.
//--------------------------------------------------------------------------//
template <typename param, typename cls_t>
class ibm_euler_slip_wall_t
{
public:
  using ibm_param_t = param;
  static constexpr int eorder_tparm = param::extrap_order;
  static constexpr bool stationary_no_slip = false;
  static constexpr bool stationary_slip = true;
  static constexpr bool isothermal_wall = false;
  static constexpr bool adiabatic_wall = false;
  static constexpr bool entropy_extension = true;
  static constexpr bool one_sided_entropy_surface = [] {
    if constexpr (requires { param::one_sided_entropy_jet; }) {
      return bool(param::one_sided_entropy_jet);
    }
    return false;
  }();
  static constexpr ibm_pressure_closure_t pressure_closure =
      ibm_pressure_closure_value<param>();
  static constexpr Real pressure_sensor_low = ibm_pressure_sensor_low<param>();
  static constexpr Real pressure_sensor_high = ibm_pressure_sensor_high<param>();

  template <int EO>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(std::integral_constant<int, EO> /*eo_tag*/,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& xyz,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& norm,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& t1,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& t2,
    const Array1D<Real, 0, EO - 1>& disIM,
    int n_valid,
    Array2D<Real, 0, EO + 1, 0, cls_t::NPRIM - 1>& q,
    int /*type_solid_bc*/, const cls_t* cls,
    const ibm_pressure_compatibility_t& pressure_compat)
  {
    q(1, cls_t::QU) = Real(0.0);
    q(1, cls_t::QV) =
        ibm_zero_grad_surface<EO>(q, cls_t::QV, disIM, n_valid);
    q(1, cls_t::QW) =
        ibm_zero_grad_surface<EO>(q, cls_t::QW, disIM, n_valid);
    q(1, cls_t::QPRES) =
        ibm_wall_pressure_surface<param, cls_t, EO, true>(
            xyz, norm, t1, t2, q, disIM, n_valid, pressure_compat);
    ibm_entropy_surface_state<
        cls_t, EO, one_sided_entropy_surface>(
            q, disIM, n_valid, cls);
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
  static constexpr bool stationary_no_slip = true;
  static constexpr bool stationary_slip = false;
  static constexpr bool isothermal_wall = true;
  static constexpr bool adiabatic_wall = false;
  static constexpr bool entropy_extension = false;
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
    int /*type_solid_bc*/, const cls_t* /*cls*/,
    const ibm_pressure_compatibility_t& pressure_compat)
  {
    q(1,cls_t::QU) = 0.0_rt;   // un  = 0
    q(1,cls_t::QV) = 0.0_rt;   // ut1 = 0
    q(1,cls_t::QW) = 0.0_rt;   // ut2 = 0
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, false>(
        xyz, norm, t1, t2, q, disIM, n_valid, pressure_compat);
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
  static constexpr bool stationary_no_slip = false;
  static constexpr bool stationary_slip = true;
  static constexpr bool isothermal_wall = false;
  static constexpr bool adiabatic_wall = true;
  static constexpr bool entropy_extension = false;
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
    int /*type_solid_bc*/, const cls_t* /*cls*/,
    const ibm_pressure_compatibility_t& pressure_compat)
  {
    q(1,cls_t::QU) = 0.0_rt;                                                       // un  = 0
    q(1,cls_t::QV) = ibm_zero_grad_surface<EO>(q, cls_t::QV, disIM, n_valid);       // ut1 = slip
    q(1,cls_t::QW) = ibm_zero_grad_surface<EO>(q, cls_t::QW, disIM, n_valid);       // ut2 = slip
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, true>(
        xyz, norm, t1, t2, q, disIM, n_valid, pressure_compat);
    q(1,cls_t::QT) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QT, disIM, n_valid);  // adiabatic T
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
  static constexpr bool stationary_no_slip = true;
  static constexpr bool stationary_slip = false;
  static constexpr bool isothermal_wall = false;
  static constexpr bool adiabatic_wall = true;
  static constexpr bool entropy_extension = false;
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
    int /*type_solid_bc*/, const cls_t* /*cls*/,
    const ibm_pressure_compatibility_t& pressure_compat)
  {
    q(1,cls_t::QU) = 0.0_rt;   // un  = 0
    q(1,cls_t::QV) = 0.0_rt;   // ut1 = 0
    q(1,cls_t::QW) = 0.0_rt;   // ut2 = 0
    q(1,cls_t::QPRES) = ibm_wall_pressure_surface<param, cls_t, EO, false>(
        xyz, norm, t1, t2, q, disIM, n_valid, pressure_compat);
    q(1,cls_t::QT) = ibm_zero_grad_surface_pos<EO>(q, cls_t::QT, disIM, n_valid);  // adiabatic T
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
