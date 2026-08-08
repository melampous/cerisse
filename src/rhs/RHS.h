#ifndef RHS_H_
#define RHS_H_

#include <type_traits>

//#include <Index.h>
#include <CNS.h>

// Euler numerical methods
#include <Weno.h>
#include <CentralKEEP.h>
#include <CentralDif.h>
#include <Riemann.h>
#include <Afd.h>
#include <Rusanov.h>
#include <Skew.h>

// viscous templates
#include <DiffusionCD.h>
#include <viscous.h>
#include <viscousLES.h>


#ifdef USE_PELEPHYSICS
#include "react.h"
#include "react_source.h"
#include "react_sourceLES.h"
#endif

namespace rhs_detail {

template <typename T, typename = void>
struct rz_axis_face_flux_is_metric : std::false_type {};

template <typename T>
struct rz_axis_face_flux_is_metric<
    T, std::void_t<decltype(T::rz_radial_axis_face_flux_is_metric)>>
    : std::bool_constant<T::rz_radial_axis_face_flux_is_metric> {};

template <typename T, typename = void>
struct rz_axis_face_addition_is_zero : std::false_type {};

template <typename T>
struct rz_axis_face_addition_is_zero<
    T, std::void_t<decltype(T::rz_radial_axis_face_addition_is_zero)>>
    : std::bool_constant<T::rz_radial_axis_face_addition_is_zero> {};

template <typename T, typename = void>
struct rz_paired_pressure_flux_capable : std::false_type {};

template <typename T>
struct rz_paired_pressure_flux_capable<
    T, std::void_t<decltype(T::rz_paired_pressure_flux_capable)>>
    : std::bool_constant<T::rz_paired_pressure_flux_capable> {};

template <typename T, typename = void>
struct rz_paired_pressure_flux_required : std::false_type {};

template <typename T>
struct rz_paired_pressure_flux_required<
    T, std::void_t<decltype(T::rz_paired_pressure_flux_required)>>
    : std::bool_constant<T::rz_paired_pressure_flux_required> {};

}  // namespace rhs_detail

// _dt stands for derived type
template <typename euler, typename diffusive, typename source>
class rhs_dt : public euler, public diffusive, public source
{
private:
public:
  static constexpr bool rz_radial_axis_face_flux_is_metric =
      rhs_detail::rz_axis_face_flux_is_metric<euler>::value &&
      rhs_detail::rz_axis_face_addition_is_zero<diffusive>::value;
  static constexpr bool rz_paired_pressure_flux_capable =
      rhs_detail::rz_paired_pressure_flux_capable<euler>::value &&
      rhs_detail::rz_axis_face_addition_is_zero<diffusive>::value;
  static constexpr bool rz_paired_pressure_flux_required =
      rhs_detail::rz_paired_pressure_flux_required<euler>::value;
};

// no euler flux
class no_euler_t
{
public:
  // A diffusion-only RHS must not receive the Euler +p/r geometric source.
  static constexpr bool rz_euler_geometric_source_active = false;
  template<typename... Args>
#if (AMREX_USE_GPIBM || CNS_USE_EB )  
  // void eflux_ibm(Args&&... args){}
  void eflux_ibm(const Geometry& /*geom*/, const MFIter& /*mfi*/,
                    const Array4<Real>& /*prims*/, std::array<FArrayBox*, AMREX_SPACEDIM> const /*&flxt*/,
                    const Array4<Real>& /*cons*/, Args&&... args ) { }
#else
  //void eflux(Args&&... args){}
  void eflux(const Geometry& /*geom*/, const MFIter& /*mfi*/,
            const Array4<Real>& /*prims*/, std::array<FArrayBox*, AMREX_SPACEDIM> const /*&flxt*/,
            const Array4<Real>& /*rhs*/, Args&&... args) { }
#endif
};

// no diffusive flux
class no_diffusive_t
{
public:
  static constexpr bool rz_radial_axis_face_addition_is_zero = true;
  static constexpr bool ibm_wall_heat_flux_capable = false;
  static constexpr bool ibm_wall_viscous_flux_capable = false;

  // No-op init, so ProbRHS::init_coeffs() is always valid
  AMREX_GPU_HOST
  void init_coeffs() {}

  template<typename... Args>
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
  //void dflux_ibm(Args&&... args) {}
  void dflux_ibm(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,            
            const Array4<Real>& rhs, Args&&... args) { }

#else
  //void dflux(Args&&... args) {}
  void dflux(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
            const Array4<Real>& rhs, Args&&... args) { }

#endif

  // no-op RZ geometric viscous source (no diffusion => no hoop stress)
  void inline rz_geometric_source(const Geometry& /*geom*/, const MFIter& /*mfi*/,
            const Array4<Real>& /*prims*/, const Array4<Real>& /*state*/,
            const auto* /*cls*/) { }
};

// no source
class no_source_t
{
public:
  template<typename... Args>
  void src(Args&&... args) {}
};

#endif
