#ifndef Skew_H_
#define Skew_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>
#include <type_traits>

#include "IBMSharedGPFluxUtils.h"

namespace skew_detail {

// Optional parameter-struct members with backward-compatible defaults.
// sensor_power: exponent applied to the discontinuity sensor in the
// second-difference (shock) dissipation term.  The default 1 reproduces the
// original Jameson-Schmidt-Turkel scaling.
template <typename P, typename = void>
struct sensor_power : std::integral_constant<int, 1> {};
template <typename P>
struct sensor_power<P, std::void_t<decltype(P::sensor_power)>>
    : std::integral_constant<int, P::sensor_power> {};

// damp_difference_order: difference order of the background damping flux
// (3 = third difference, 5 = fifth difference).  The sentinel 0 selects the
// historical stencil bound to the base order (order-1).
template <typename P, typename = void>
struct damp_difference_order : std::integral_constant<int, 0> {};
template <typename P>
struct damp_difference_order<
    P, std::void_t<decltype(P::damp_difference_order)>>
    : std::integral_constant<int, P::damp_difference_order> {};

// The shared face arrays normally store a per-unit-area flux.  In the radial
// R-Z direction Skew/JST is constructed for the metric variables rF and rU.
// Away from the axis it is divided by the face radius before storage.  The
// axis slot instead stores the auxiliary metric flux directly because the
// finite-difference h-flux of rF need not vanish at r=0.
struct radial_metric_t {
  bool active;
  amrex::Real origin;
  amrex::Real spacing;
  int axis_face_index;
  int domain_high_index;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real face_radius(const int face_index) const noexcept
  {
    return origin +
           amrex::Real(face_index - axis_face_index) * spacing;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool is_axis_face(const int face_index) const noexcept
  {
    return active && origin == amrex::Real(0.0) &&
           face_index == axis_face_index;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool is_innermost_axis_interior_face(
      const int face_index) const noexcept
  {
    return active && origin == amrex::Real(0.0) &&
           face_index == axis_face_index + 1;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real cell_to_stored_flux_scale(
      const int cell_index, const int face_index) const noexcept
  {
    if (!active) {
      return amrex::Real(1.0);
    }
    const amrex::Real cell_radius =
        origin +
        (amrex::Real(cell_index - axis_face_index) + amrex::Real(0.5)) *
            spacing;
    if (is_axis_face(face_index)) {
      return cell_radius;
    }
    return cell_radius / face_radius(face_index);
  }
};

}  // namespace skew_detail


//-------------------
// discontinuity sensor function
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real disconSensor(Real pp, Real pl,
                                                      Real pr) {
  Real pjst = pr + 2.0_rt * pp + pl;
  Real ptvd = std::abs(pr - pp) + std::abs(pp - pl);
  return std::abs(2.0_rt * (pr - 2.0_rt * pp + pl) /
                  (pjst + ptvd + Real(1.0e-40)));
  }

// int imask
template <typename param, typename cls_t>
class skew_t {
  public:
  using DirectionalFaceFluxArrays =
      std::array<FArrayBox*, AMREX_SPACEDIM>;

  static constexpr const char* scheme_name =
      param::dissipation
          ? (param::order == 4 ? "skew4-jst"
                               : (param::order == 6 ? "skew6-jst"
                                                    : "skew-jst"))
          : (param::order == 4 ? "skew4-central"
                               : (param::order == 6 ? "skew6-central"
                                                    : "skew-central"));

  // At the radial axis the face array contains the auxiliary metric flux
  // h_(rF), not a per-unit-area flux.  The R-Z divergence assembly already
  // supports this storage convention.
  static constexpr bool rz_radial_axis_face_flux_is_metric = true;
  // In a selected first-ring rescue the UMX axis slot contains the auxiliary
  // metric advective flux required by the odd radial-momentum parity closure.
  static constexpr bool rz_paired_axis_advective_flux_is_metric = true;
  static constexpr bool rz_paired_pressure_flux_capable = true;
  static constexpr bool rz_paired_pressure_flux_required = true;
  static constexpr bool rz_annular_pressure_consistency_capable = false;
  static constexpr bool rz_euler_geometric_source_active = true;

  AMREX_GPU_HOST_DEVICE
  skew_t() {
    
    // init coefficients for skew-symmetric
    for (int l = 0; l < order; l++) {       
      for (int m = 0; m < order; m++) { 
        coefskew(l,m) = Real(0.0); 
      }
    }
    // init interpolation, damping and shock coeffcients      
    for (int l = 0; l < order; l++) {       
      coefdamp(l)  = Real(0.0);
      coefP(l)     = Real(0.0);
      coefshock(l) = Real(0.0);
    }    
    
    switch (order)
    {
    case 2:
      coefskew(0,0) =  Real(0.25);
      coefskew(1,0) =  Real(0.25);
      coefskew(0,1) =  Real(0.25);
      coefskew(1,1) =  Real(0.25);
      // pressure interpolation
      coefP(0)      =  Real(0.5);
      coefP(1)      =  Real(0.5);
      // damping coefficients
      coefdamp(0)   = Real(0.0);
      coefdamp(1)   = Real(0.0);
      // shock
      coefshock(0) =  -Real(1.0);
      coefshock(1) =  +Real(1.0);            
      break;
    case 4:
      coefskew(0,0) =  - Real(1.0/24.0);
      coefskew(0,1) =  Real(0.0);
      coefskew(0,2) =  - Real(1.0/24.0);
      coefskew(0,3) =   Real(0.0);      

      coefskew(1,0) =  Real(0.0);
      coefskew(1,1) =  Real(1.0/3.0) - Real(1.0/24.0);  
      coefskew(1,2) =  Real(1.0/3.0);
      coefskew(1,3) =  - Real(1.0/24.0);
     
      coefskew(2,0) =  - Real(1.0/24.0);      
      coefskew(2,1) =  Real(1.0/3.0);        
      coefskew(2,2) =  Real(1.0/3.0) - Real(1.0/24.0); 
      coefskew(2,3) =  - Real(0.0);
     
      coefskew(3,0) =  - Real(0.0);
      coefskew(3,1) =  - Real(1.0/24.0);
      coefskew(3,2) =  - Real(0.0);
      coefskew(3,3) =  - Real(1.0/24.0);
      
      // pressure interpolation
      coefP(0) = -Real(1.0/12.0);
      coefP(1) =  Real(7.0/12.0);
      coefP(2) =  Real(7.0/12.0);
      coefP(3) = -Real(1.0/12.0);   

      // damping coefficients     Ui+2 - 3 Ui+1 + 3 Ui - Ui-1
      coefdamp(0) = -Real(1.0);
      coefdamp(1) = +Real(3.0);
      coefdamp(2) = -Real(3.0);
      coefdamp(3) = +Real(1.0);

      // shock coefficients
      coefshock(0) = -Real(1.0/6.0);
      coefshock(1) = -Real(5.0/6.0);
      coefshock(2) = +Real(7.0/6.0);
      coefshock(3) = -Real(1.0/6.0);

      break; 
    case 6:       
      //   F(i+1/2) = beta(k,p) U(i+p) V(i+k)

      // coefskew(-2+ k, -2+p) = beta(k,p) 

      // beta(k= -2 ,p)   U(i+p)V(i-2)              
      coefskew(0,0) = Real(1.0/120.0);  //  p = -2
      coefskew(0,1) = Real(0.0);        //  -1
      coefskew(0,2) = Real(0.0);        //  0
      coefskew(0,3) = Real(1.0/120.0);  //  1
      coefskew(0,4) = Real(0.0);        //  2 
      coefskew(0,5) = Real(0.0);        //  3 
      
      // beta(k= -1 ,p)  U(i+p)V(i-1)
      coefskew(1,0) =  Real(0.0);         // p = -2
      coefskew(1,1) = -Real(8.0/120.0);   // -1
      coefskew(1,2) =  Real(0.0);         // 0
      coefskew(1,3) = -Real(9.0/120.0);   // 1 
      coefskew(1,4) =  Real(1.0/120.0);   // 2 
      coefskew(1,5) =  Real(0.0);         // 3

      // beta(k= 0 ,p)  U(i+p)V(i) 
      coefskew(2,0) =  Real(0.0);         // -2
      coefskew(2,1) =  Real(0.0);         // -1
      coefskew(2,2) =  Real(37.0/120.0);  // 0
      coefskew(2,3) =  Real(45.0/120.0);  // 1
      coefskew(2,4) = -Real(9.0/120.0);   // 2 
      coefskew(2,5) =  Real(1.0/120.0);   // 3

      // beta(k= 1 ,p) U(i+p)V(i+1)
      coefskew(3,0) =  Real(1.0/120.0);   // -2
      coefskew(3,1) = -Real(9.0/120.0);   // -1 
      coefskew(3,2) =  Real(45.0/120.0);  // 0  
      coefskew(3,3) =  Real(37.0/120.0);  // 1  
      coefskew(3,4) =  Real(0.0);         // 2
      coefskew(3,5) =  Real(0.0);         // 3

      // beta(k= 2 ,p) U(i+p)V(i+2)
      coefskew(4,0) =  Real(0.0);         // -2
      coefskew(4,1) =  Real(1.0/120.0);   // -1
      coefskew(4,2) = -Real(9.0/120.0);   // 0 
      coefskew(4,3) =  Real(0.0);         // 1
      coefskew(4,4) = -Real(8.0/120.0);   // 2
      coefskew(4,5) =  Real(0.0);         // 3

      // beta(k= 3 ,p) U(i+p)V(i+3)
      coefskew(5,0) = Real(0.0);          // -2
      coefskew(5,1) = Real(0.0);          // -1
      coefskew(5,2) = Real(1.0/120.0);    // 0
      coefskew(5,3) = Real(0.0);          // 1      
      coefskew(5,4) = Real(0.0);          // 2
      coefskew(5,5) = Real(1.0/120.0);    // 3

      // no FV correction
        
      // pressure interpolation 
      coefP(0) = Real(1.0/60.0);
      coefP(1) = -Real(8.0/60.0);
      coefP(2) = Real(37.0/60.0);
      coefP(3) = Real(37.0/60.0);        
      coefP(4) = -Real(8.0/60.0);
      coefP(5) = Real(1.0/60.0);

      // damping coefficients  -(Ui+3 -5 U i+2 + 10 Ui+1 -10 Ui  + 5 Ui+1 - Ui-2)
      coefdamp(0) = Real(1.0);
      coefdamp(1) = -Real(5.0);
      coefdamp(2) = +Real(10.0);
      coefdamp(3) = -Real(10.0);
      coefdamp(4) = +Real(5.0);
      coefdamp(5) = -Real(1.0);
      
      // shock coefficients
      coefshock(0) =  Real(0.0);      
      coefshock(1) = -Real(1.0/6.0);
      coefshock(2) = -Real(5.0/6.0);
      coefshock(3) = +Real(7.0/6.0);
      coefshock(4) = -Real(1.0/6.0);
      coefshock(5) =  Real(0.0);                             
      break;

    default: 

      amrex::Abort("Skew-symmetric order available 2/4/6 MUST specify one of those ");      

      break;
    }
    
    // Decoupled background-damping coefficients: minus the third or fifth
    // difference evaluated at the face, over damp_width points starting at
    // face_index - damp_halfsten.
    for (int l = 0; l <= 6; ++l) { coefdamp_ext(l) = Real(0.0); }
    if (damp_difference_order == 3) {
      coefdamp_ext(0) = -Real(1.0);
      coefdamp_ext(1) = +Real(3.0);
      coefdamp_ext(2) = -Real(3.0);
      coefdamp_ext(3) = +Real(1.0);
    } else {
      coefdamp_ext(0) = +Real(1.0);
      coefdamp_ext(1) = -Real(5.0);
      coefdamp_ext(2) = +Real(10.0);
      coefdamp_ext(3) = -Real(10.0);
      coefdamp_ext(4) = +Real(5.0);
      coefdamp_ext(5) = -Real(1.0);
    }

    // inside skew_t() constructor
    NSEN[0] = cls_t::QRHO;
    NSEN[1] = cls_t::QPRES;
// #if NUM_SPECIES > 1    
//     for (int nv = 2; nv < NVARSEN; nv++) {
//       NSEN[nv] = cls_t::QFS + (nv - 2);
//     }
// #endif

    // (manual sensor) param::manualsensor
    // const int NVARSEN=param::NVARSEN;
    // int NSEN[NVARSEN];
    //  for (int nv=0; nv<NVARSEN; nv++) {
    //  NSEN[nv] = param::NSEN[nv];    
    // }
    // 

    // no masking
    for (int l = 0; l < AMREX_SPACEDIM; l++) {  
      mask_sen(l,1) = -1;     
      mask_sen(l,2) = 9999999;     
    }  
    
  }

  AMREX_GPU_HOST_DEVICE
  ~skew_t() {}

  struct LocalLlfFallbackOptions {
    int enabled;
    Real pressure_sensor_threshold;
    Real density_sensor_threshold;
    Real density_jump_threshold;
    Real compression_threshold;
    int axis_collar_faces;
  };

  static const LocalLlfFallbackOptions& local_llf_fallback_options()
  {
    static const LocalLlfFallbackOptions options = [] {
      LocalLlfFallbackOptions value{
          1, Real(1.0e-1), Real(1.0e-1), Real(3.0e-3), Real(1.0e-3),
          0};
      amrex::ParmParse parameters("cns");
      parameters.query("skew_first_order_fallback", value.enabled);
      parameters.query(
          "skew_fallback_pressure_sensor",
          value.pressure_sensor_threshold);
      parameters.query(
          "skew_fallback_density_sensor",
          value.density_sensor_threshold);
      parameters.query(
          "skew_fallback_density_jump", value.density_jump_threshold);
      parameters.query(
          "skew_fallback_compression", value.compression_threshold);
      parameters.query(
          "skew_axis_collar_faces", value.axis_collar_faces);

      if (value.enabled != 0 && value.enabled != 1) {
        amrex::Abort("cns.skew_first_order_fallback must be 0 or 1");
      }
      const auto valid_unit_threshold = [](const Real threshold) {
        return std::isfinite(threshold) && threshold >= Real(0.0) &&
               threshold <= Real(1.0);
      };
      if (!valid_unit_threshold(value.pressure_sensor_threshold)) {
        amrex::Abort(
            "cns.skew_fallback_pressure_sensor must be in [0,1]");
      }
      if (!valid_unit_threshold(value.density_sensor_threshold)) {
        amrex::Abort(
            "cns.skew_fallback_density_sensor must be in [0,1]");
      }
      if (!std::isfinite(value.density_jump_threshold) ||
          value.density_jump_threshold < Real(0.0)) {
        amrex::Abort(
            "cns.skew_fallback_density_jump must be finite and non-negative");
      }
      if (!valid_unit_threshold(value.compression_threshold)) {
        amrex::Abort("cns.skew_fallback_compression must be in [0,1]");
      }
      if (value.axis_collar_faces < 0 || value.axis_collar_faces > 16) {
        amrex::Abort("cns.skew_axis_collar_faces must be in [0,16]");
      }
      return value;
    }();
    return options;
  }

  static void print_local_llf_fallback_manifest()
  {
    const auto& options = local_llf_fallback_options();
    amrex::Print()
        << "[Numerics] skew_local_first_order_llf_fallback="
        << (options.enabled != 0 ? "enabled" : "disabled")
        << " pressure_sensor=" << options.pressure_sensor_threshold
        << " density_sensor=" << options.density_sensor_threshold
        << " density_jump=" << options.density_jump_threshold
        << " compression=" << options.compression_threshold
        << " rz_axis_collar_faces=" << options.axis_collar_faces
        << " rz_axis_collar_sensor=density_or_paired_pressure_jump"
        << " rz_axis_collar_closure=metric_llf_plus_parity_axis_aux_and_h2"
        << " scope=cartesian_rz_all_fluid_shared_gp"
        << " action=replace_complete_face_flux_and_selected_axis_rescue\n";
  }

  // parameters for damping and shock capturing   
  static const int order = param::order;
  // Dissipation configuration (backward-compatible defaults; see the
  // skew_detail traits).  damp_difference_order 3 or 5 selects the third- or
  // fifth-difference background damping flux.  The historical default couples
  // it to the base order for orders four and six.  Skew2 has no background
  // damping term and uses the third-difference width only as a valid sentinel.
  static constexpr int sensor_power =
      skew_detail::sensor_power<param>::value;
  static constexpr int damp_difference_order =
      skew_detail::damp_difference_order<param>::value == 0
          ? (order == 2 ? 3 : order - 1)
          : skew_detail::damp_difference_order<param>::value;
  static constexpr int damp_width = damp_difference_order + 1;
  static constexpr int damp_halfsten = damp_width / 2;
  static_assert(damp_difference_order == 3 || damp_difference_order == 5,
                "skew_t: damp_difference_order must be 3 or 5");
  static_assert(sensor_power >= 1 && sensor_power <= 3,
                "skew_t: sensor_power must be 1, 2, or 3");
  static_assert(cls_t::NGHOST >= damp_width / 2,
                "skew_t: extended damping stencil exceeds NGHOST");
#if (AMREX_USE_GPIBM || CNS_USE_EB)
  static_assert(order == 2 || damp_difference_order == order - 1,
                "skew_t: extended damping is not implemented for the "
                "IBM/EB dissipation path");
#endif
  // default values Cshock = 0.1 Cdamp = 0.016; 
  Real Cshock = param::C2skew;
  Real Cdamp  = param::C4skew; 


#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void eflux_ibm(
      const Geometry& geom, const MFIter& mfi,
      const Array4<const Real>& prims,
      const DirectionalFaceFluxArrays& flxt,
      const Array4<Real>& cons, const cls_t* cls,
      const Array4<uint8_t>& ibMarkers)
  {
    Array4<Real> unused_pressure_face_flux;
    compute_face_fluxes<false>(
        geom, mfi, prims, flxt, cons, cls, ibMarkers,
        prims, cls_t::QPRES, unused_pressure_face_flux);
  }

#ifdef AMREX_USE_GPIBM
  void eflux_ibm_with_rz_paired_pressure(
      const Geometry& geom, const MFIter& mfi,
      const Array4<const Real>& prims,
      const DirectionalFaceFluxArrays& flxt,
      const Array4<Real>& cons, const cls_t* cls,
      const Array4<uint8_t>& ibMarkers,
      const Array4<const Real>& pressure_states,
      const int pressure_component,
      const Array4<Real>& radial_pressure_face_flux)
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.IsRZ(),
        "Skew paired pressure flux requires cylindrical R-Z geometry");
    compute_face_fluxes<true>(
        geom, mfi, prims, flxt, cons, cls, ibMarkers,
        pressure_states, pressure_component, radial_pressure_face_flux);
  }
#endif

  template <bool ComputeRzPairedPressure>
  void compute_face_fluxes(
      const Geometry& geom, const MFIter& mfi,
      const Array4<const Real>& prims,
      const DirectionalFaceFluxArrays& flxt,
      const Array4<Real>& cons, const cls_t* cls,
      const Array4<uint8_t>& ibMarkers,
      const Array4<const Real>& pressure_states,
      const int pressure_component,
      const Array4<Real>& radial_pressure_face_flux)
  {
#else
  void eflux(
      const Geometry& geom, const MFIter& mfi,
      const Array4<const Real>& prims,
      const DirectionalFaceFluxArrays& flxt,
      const Array4<Real>& cons, const cls_t* cls)
  {
    Array4<Real> unused_pressure_face_flux;
    compute_face_fluxes<false>(
        geom, mfi, prims, flxt, cons, cls,
        prims, cls_t::QPRES, unused_pressure_face_flux);
  }

  void eflux_with_rz_paired_pressure(
      const Geometry& geom, const MFIter& mfi,
      const Array4<const Real>& prims,
      const DirectionalFaceFluxArrays& flxt,
      const Array4<Real>& cons, const cls_t* cls,
      const Array4<const Real>& pressure_states,
      const int pressure_component,
      const Array4<Real>& radial_pressure_face_flux)
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.IsRZ(),
        "Skew paired pressure flux requires cylindrical R-Z geometry");
    compute_face_fluxes<true>(
        geom, mfi, prims, flxt, cons, cls,
        pressure_states, pressure_component, radial_pressure_face_flux);
  }

  template <bool ComputeRzPairedPressure>
  void compute_face_fluxes(
      const Geometry& geom, const MFIter& mfi,
      const Array4<const Real>& prims,
      const DirectionalFaceFluxArrays& flxt,
      const Array4<Real>& cons, const cls_t* cls,
      const Array4<const Real>& pressure_states,
      const int pressure_component,
      const Array4<Real>& radial_pressure_face_flux)
  {
#endif

    const auto cell_size = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const bool is_rz = geom.IsRZ();
    if constexpr (!ComputeRzPairedPressure) {
      amrex::ignore_unused(
          pressure_states, pressure_component,
          radial_pressure_face_flux);
    }
    const auto fallback_options = local_llf_fallback_options();
    const int fallback_enabled = fallback_options.enabled;
    const Real fallback_pressure_sensor =
        fallback_options.pressure_sensor_threshold;
    const Real fallback_density_sensor =
        fallback_options.density_sensor_threshold;
    const Real fallback_density_jump =
        fallback_options.density_jump_threshold;
    const Real fallback_compression =
        fallback_options.compression_threshold;
    const int fallback_axis_collar_faces =
        fallback_options.axis_collar_faces;

    //const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls->NGHOST);
    //FArrayBox consf(bxg, cls_t::NCONS, The_Async_Arena());
    FArrayBox lambda_maxf(bxg, 1, The_Async_Arena());

    // create lambda(0) and cons array
    // Array4<Real> cons   = consf.array();
    Array4<Real> lambda_max = lambda_maxf.array();

    // copy conservative variables from rhs to cons and clear rhs
    // ParallelFor(bxg, cls_t::NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
    //   cons(i,j,k,n) = rhs(i,j,k,n);
    //   rhs(i, j, k, n)=0.0;              
    //   });

    // int imask = 3; // reduce order 
    // const int* domlo = geom.Domain().loVect();
    // const int* domhi = geom.Domain().hiVect();

    // masking BC here function of global geometry  (only non-periodic dirs)       
    // for (int l = 0; l < AMREX_SPACEDIM; l++) {  
    //   if (geom.isPeriodic(l)==0) {
    //     mask_sen(l,1) = domlo[l] + imask;     
    //     mask_sen(l,2) = domhi[l] - imask;     
    //   }
    // }  

    // ---------------------------------------------------------------------  //
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
      const skew_detail::radial_metric_t radial_metric{
          is_rz && dir == 0, prob_lo[0], cell_size[0],
          geom.Domain().smallEnd(0), geom.Domain().bigEnd(0)};
      const bool paired_radial_flux =
          ComputeRzPairedPressure && radial_metric.active;

      int Qdir =  cls_t::QRHO + dir + 1; 

      auto const& flx = flxt[dir]->array(); 
      const Box bxface = mfi.grownnodaltilebox(dir, 0);
  

#if (AMREX_USE_GPIBM || CNS_USE_EB )  
      ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->template flux_dir_ibm<ComputeRzPairedPressure>(
                        i, j, k, Qdir, vdir, cons, prims, lambda_max,
                        flx, cls, ibMarkers, radial_metric,
                        pressure_states, pressure_component,
                        radial_pressure_face_flux);
                  });
#else    
      ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->template flux_dir<ComputeRzPairedPressure>(
                        i, j, k, Qdir, vdir, cons, prims, lambda_max,
                        flx, cls, radial_metric, pressure_states,
                        pressure_component, radial_pressure_face_flux);
                  });                  
#endif

      
      // dissipative fluxes 
      if constexpr (param::dissipation) {

#if (AMREX_USE_GPIBM || CNS_USE_EB )  

        ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->fluxdissip_dir_ibm(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls, ibMarkers, radial_metric);
                  });             
#else
        ParallelFor(bxface,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->fluxdissip_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls, radial_metric);
                  });  
#endif
        
      }

      // A troubled face must contain either the completed Skew/JST flux or the
      // complete two-state LLF flux.  The LLF contribution is therefore a
      // final overwrite, not another term added to the central/JST result.
      if (fallback_enabled != 0) {
#if AMREX_USE_GPIBM
        ParallelFor(
            bxface,
            [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              this->template apply_local_llf_fallback_ibm<
                  ComputeRzPairedPressure>(
                  i, j, k, Qdir, vdir, prims, flx, cls, ibMarkers,
                  radial_metric, fallback_pressure_sensor,
                  fallback_density_sensor, fallback_density_jump,
                  fallback_compression, fallback_axis_collar_faces,
                  pressure_states,
                  pressure_component, radial_pressure_face_flux);
            });
#elif !CNS_USE_EB
        ParallelFor(
            bxface,
            [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              this->template apply_local_llf_fallback<
                  ComputeRzPairedPressure>(
                  i, j, k, Qdir, vdir, prims, flx, cls, radial_metric,
                  fallback_pressure_sensor, fallback_density_sensor,
                  fallback_density_jump, fallback_compression,
                  fallback_axis_collar_faces,
                  pressure_states, pressure_component,
                  radial_pressure_face_flux);
            });
#endif
      }

      if (paired_radial_flux) {
        ParallelFor(
            bxface,
            [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              if (radial_metric.is_axis_face(i)) {
                const IntVect axis_face(AMREX_D_DECL(i, j, k));
                const IntVect radial_direction =
                    IntVect::TheDimensionVector(0);
                const IntVect second_ring = axis_face + radial_direction;
                const IntVect first_interior_face = second_ring;
                bool collar_states_usable =
                    fallback_enabled != 0 &&
                    second_ring[0] <= radial_metric.domain_high_index;
                bool four_point_sensor_available = true;
#if AMREX_USE_GPIBM
                collar_states_usable =
                    collar_states_usable &&
                    !ibm_flux::is_solid_solid_face(
                        first_interior_face, radial_direction, ibMarkers) &&
                    this->gp_state_usable(
                        axis_face, ibMarkers, radial_metric) &&
                    this->gp_state_usable(
                        second_ring, ibMarkers, radial_metric);
                four_point_sensor_available =
                    collar_states_usable &&
                    this->gp_state_usable(
                        axis_face - radial_direction,
                        ibMarkers, radial_metric) &&
                    this->gp_state_usable(
                        second_ring + radial_direction,
                        ibMarkers, radial_metric);
#elif CNS_USE_EB
                // The legacy EB-only path has no local fallback pass.  Keep
                // its existing axis closure unchanged.
                collar_states_usable = false;
#endif
                bool axis_collar_selected = false;
                if (collar_states_usable) {
                  const auto indicators = this->troubled_face_indicators(
                      first_interior_face, radial_direction, Qdir, prims,
                      four_point_sensor_available);
                  axis_collar_selected =
                      (fallback_axis_collar_faces > 0 &&
                       this->axis_collar_fallback_needed(
                           first_interior_face, radial_direction, prims,
                           pressure_states, pressure_component,
                           fallback_density_jump)) ||
                      this->local_llf_fallback_needed(
                          indicators, fallback_pressure_sensor,
                          fallback_density_sensor, fallback_density_jump,
                          fallback_compression);
                }

                if (axis_collar_selected) {
                  // Reproduce the actual first-interior-face selection and
                  // close it with a parity-aware axis auxiliary flux.  For
                  // n!=UMX the shifted axis value makes the first-ring
                  // difference algebraically identical to a standard
                  // physical two-state LLF flux with zero physical axis
                  // area.  UMX additionally cancels the leading LLF jump of
                  // the odd radial momentum.  This is a nonlinear axis rescue,
                  // not an invariant-domain proof for the full Skew scheme.
                  const Real axis_pressure =
                      pressure_states(axis_face, pressure_component);
                  const Real face_radius =
                      radial_metric.face_radius(first_interior_face[0]);
                  const Real ring1_radius =
                      radial_metric.origin + Real(0.5) * radial_metric.spacing;
                  const Real ring2_radius =
                      ring1_radius + radial_metric.spacing;
                  const Real alpha = amrex::max(
                      std::abs(prims(axis_face, cls_t::QU)) +
                          prims(axis_face, cls_t::QC),
                      std::abs(prims(second_ring, cls_t::QU)) +
                          prims(second_ring, cls_t::QC));
                  Real conservative_ring1[cls_t::NCONS];
                  Real physical_flux_ring1[cls_t::NCONS];
                  Real physical_flux_ring2[cls_t::NCONS];
                  cls->prims2cons(
                      axis_face, prims, conservative_ring1);
                  cls->prims2flux(
                      axis_face, 0, prims, physical_flux_ring1);
                  cls->prims2flux(
                      second_ring, 0, prims, physical_flux_ring2);
                  physical_flux_ring1[cls_t::UMX] -=
                      prims(axis_face, cls_t::QPRES);
                  physical_flux_ring2[cls_t::UMX] -=
                      prims(second_ring, cls_t::QPRES);

                  bool valid = std::isfinite(axis_pressure) &&
                               axis_pressure > Real(0.0) &&
                               std::isfinite(alpha) && alpha > Real(0.0) &&
                               face_radius > Real(0.0);
                  for (int n = 0; n < cls_t::NCONS; ++n) {
                    valid = valid &&
                            std::isfinite(conservative_ring1[n]) &&
                            std::isfinite(physical_flux_ring1[n]) &&
                            std::isfinite(physical_flux_ring2[n]);
                  }
                  if (!valid) {
                    this->template mark_flux_pair_invalid<true>(
                        axis_face, flx, radial_metric,
                        radial_pressure_face_flux);
                  } else {
                    for (int n = 0; n < cls_t::NCONS; ++n) {
                      flx(axis_face, n) =
                          Real(0.5) *
                          ((ring1_radius - face_radius) *
                               physical_flux_ring1[n] +
                           (ring2_radius - face_radius) *
                               physical_flux_ring2[n]);
                    }
                    // m_r is odd.  Replacing the ordinary jump m2-m1 by the
                    // parity-consistent m2-3m1 removes the O(1) first-ring
                    // LLF defect for m_r=a*r+O(r^3).
                    flx(axis_face, cls_t::UMX) -=
                        face_radius * alpha *
                        conservative_ring1[cls_t::UMX];
                    radial_pressure_face_flux(axis_face, 0) = axis_pressure;
                  }
                } else {
                  flx(axis_face, cls_t::UMX) = Real(0.0);
                }
              }
            });
      }

    }

  }

  struct TroubledFaceIndicators {
    Real pressure_sensor;
    Real density_sensor;
    Real relative_density_jump;
    Real compression;
  };

  template <typename PrimitiveArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE TroubledFaceIndicators
  troubled_face_indicators(
      const IntVect& face, const IntVect& direction,
      const int normal_velocity_component,
      const PrimitiveArray& prims,
      const bool four_point_sensor_available) const noexcept
  {
    const IntVect left = face - direction;
    const Real density_left = prims(left, cls_t::QRHO);
    const Real density_right = prims(face, cls_t::QRHO);
    const Real pressure_left = prims(left, cls_t::QPRES);
    const Real pressure_right = prims(face, cls_t::QPRES);
    const Real tiny = std::numeric_limits<Real>::min();

    const Real density_floor =
        amrex::max(amrex::min(density_left, density_right), tiny);
    const Real relative_density_jump =
        std::abs(density_right - density_left) / density_floor;

    Real density_sensor;
    Real pressure_sensor;
    if (four_point_sensor_available) {
      const IntVect left_outer = left - direction;
      const IntVect right_outer = face + direction;
      density_sensor = amrex::max(
          disconSensor(
              density_left, prims(left_outer, cls_t::QRHO), density_right),
          disconSensor(
              density_right, density_left,
              prims(right_outer, cls_t::QRHO)));
      pressure_sensor = amrex::max(
          disconSensor(
              pressure_left, prims(left_outer, cls_t::QPRES),
              pressure_right),
          disconSensor(
              pressure_right, pressure_left,
              prims(right_outer, cls_t::QPRES)));
    } else {
      // Near a GP stencil boundary only the two states sharing the face are
      // guaranteed to exist.  Their relative jumps provide a conservative
      // detector without reading an unreconstructed deep-solid state.
      density_sensor = amrex::min(relative_density_jump, Real(1.0));
      pressure_sensor = amrex::min(
          std::abs(pressure_right - pressure_left) /
              amrex::max(
                  amrex::max(std::abs(pressure_left),
                             std::abs(pressure_right)),
                  tiny),
          Real(1.0));
    }

    const Real velocity_left = prims(left, normal_velocity_component);
    const Real velocity_right = prims(face, normal_velocity_component);
    const Real acoustic_scale = amrex::max(
        amrex::max(prims(left, cls_t::QC), prims(face, cls_t::QC)), tiny);
    const Real compression =
        amrex::max(velocity_left - velocity_right, Real(0.0)) /
        acoustic_scale;

    return {pressure_sensor, density_sensor, relative_density_jump,
            compression};
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool local_llf_fallback_needed(
      const TroubledFaceIndicators& indicators,
      const Real pressure_sensor_threshold,
      const Real density_sensor_threshold,
      const Real density_jump_threshold,
      const Real compression_threshold) const noexcept
  {
    const bool shock_face =
        indicators.pressure_sensor > pressure_sensor_threshold &&
        indicators.compression > compression_threshold;
    const bool contact_face =
        indicators.density_sensor > density_sensor_threshold &&
        indicators.relative_density_jump > density_jump_threshold;
    return shock_face || contact_face;
  }

  template <typename PrimitiveArray, typename PressureArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool axis_collar_fallback_needed(
      const IntVect& face, const IntVect& radial_direction,
      const PrimitiveArray& prims, const PressureArray& pressure_states,
      const int pressure_component,
      const Real jump_threshold) const noexcept
  {
    const IntVect left = face - radial_direction;
    const Real density_left = prims(left, cls_t::QRHO);
    const Real density_right = prims(face, cls_t::QRHO);
    const Real pressure_left = pressure_states(left, pressure_component);
    const Real pressure_right = pressure_states(face, pressure_component);

    // A non-finite/non-positive input is already outside the admissible
    // state.  Select the fail-closed low-order path; its finite checks poison
    // the face instead of silently retaining a high-order flux.
    if (!std::isfinite(density_left) || !std::isfinite(density_right) ||
        !std::isfinite(pressure_left) || !std::isfinite(pressure_right) ||
        !(density_left > Real(0.0)) || !(density_right > Real(0.0)) ||
        !(pressure_left > Real(0.0)) || !(pressure_right > Real(0.0))) {
      return true;
    }

    const Real density_jump =
        std::abs(density_right - density_left) /
        amrex::min(density_left, density_right);
    const Real pressure_jump =
        std::abs(pressure_right - pressure_left) /
        amrex::min(pressure_left, pressure_right);
    // Density catches contacts; paired pressure catches acoustic/pressure
    // jumps (including equal-density ones).  Compression is intentionally not
    // required at the symmetry axis.
    return density_jump > jump_threshold ||
           pressure_jump > jump_threshold;
  }

  template <bool ComputeRzPairedPressure, typename PrimitiveArray,
            typename FaceFluxArray, typename PressureArray,
            typename PressureFluxArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void write_first_order_llf_flux(
      const IntVect& face, const int direction_index,
      const PrimitiveArray& prims, const FaceFluxArray& flux,
      const cls_t* closure,
      const skew_detail::radial_metric_t radial_metric,
      const PressureArray& pressure_states,
      const int pressure_component,
      const PressureFluxArray& radial_pressure_face_flux) const noexcept
  {
    const IntVect direction =
        IntVect::TheDimensionVector(direction_index);
    const IntVect left = face - direction;

    Real conservative_left[cls_t::NCONS];
    Real conservative_right[cls_t::NCONS];
    Real physical_flux_left[cls_t::NCONS];
    Real physical_flux_right[cls_t::NCONS];
    closure->prims2cons(left, prims, conservative_left);
    closure->prims2cons(face, prims, conservative_right);
    closure->prims2flux(
        left, direction_index, prims, physical_flux_left);
    closure->prims2flux(
        face, direction_index, prims, physical_flux_right);

    const Real alpha = amrex::max(
        std::abs(prims(left, cls_t::QU + direction_index)) +
            prims(left, cls_t::QC),
        std::abs(prims(face, cls_t::QU + direction_index)) +
            prims(face, cls_t::QC));
    if (!(alpha > Real(0.0)) || !std::isfinite(alpha)) {
      const Real invalid = std::numeric_limits<Real>::quiet_NaN();
      for (int n = 0; n < cls_t::NCONS; ++n) {
        flux(face, n) = invalid;
      }
      if constexpr (ComputeRzPairedPressure) {
        if (radial_metric.active) {
          radial_pressure_face_flux(face, 0) = invalid;
        }
      }
      return;
    }

    if (radial_metric.active) {
      const Real face_radius = radial_metric.face_radius(face[0]);
      if (!(face_radius > Real(0.0))) {
        // The generic two-state overwrite is intentionally excluded from the
        // symmetry axis.  Its high-order auxiliary metric flux is supplied by
        // the parity closure already used by the Skew operator.
        return;
      }
      const Real left_scale =
          radial_metric.cell_to_stored_flux_scale(left[0], face[0]);
      const Real right_scale =
          radial_metric.cell_to_stored_flux_scale(face[0], face[0]);
      Real pressure_face = Real(0.0);
      if constexpr (ComputeRzPairedPressure) {
        physical_flux_left[cls_t::UMX] -=
            prims(left, cls_t::QPRES);
        physical_flux_right[cls_t::UMX] -=
            prims(face, cls_t::QPRES);
        pressure_face = Real(0.5) *
                        (pressure_states(left, pressure_component) +
                         pressure_states(face, pressure_component));
        if (!std::isfinite(pressure_face)) {
          const Real invalid = std::numeric_limits<Real>::quiet_NaN();
          for (int n = 0; n < cls_t::NCONS; ++n) {
            flux(face, n) = invalid;
          }
          radial_pressure_face_flux(face, 0) = invalid;
          return;
        }
      }
      for (int n = 0; n < cls_t::NCONS; ++n) {
        flux(face, n) =
            Real(0.5) *
                (left_scale * physical_flux_left[n] +
                 right_scale * physical_flux_right[n]) -
            Real(0.5) * alpha *
                (conservative_right[n] - conservative_left[n]);
      }
      if constexpr (ComputeRzPairedPressure) {
        flux(face, cls_t::UMX) += pressure_face;
        radial_pressure_face_flux(face, 0) = pressure_face;
      }
      return;
    }

    for (int n = 0; n < cls_t::NCONS; ++n) {
      flux(face, n) =
          Real(0.5) *
          (physical_flux_left[n] + physical_flux_right[n] -
           alpha * (conservative_right[n] - conservative_left[n]));
    }
  }

  template <bool ComputeRzPairedPressure,
            typename PressureArray, typename PressureFluxArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void apply_local_llf_fallback(
      const int i, const int j, const int k,
      const int normal_velocity_component,
      const GpuArray<int, 3>& direction_components,
      const Array4<const Real>& prims, const Array4<Real>& flux,
      const cls_t* closure,
      const skew_detail::radial_metric_t radial_metric,
      const Real pressure_sensor_threshold,
      const Real density_sensor_threshold,
      const Real density_jump_threshold,
      const Real compression_threshold,
      const int axis_collar_faces,
      const PressureArray& pressure_states,
      const int pressure_component,
      const PressureFluxArray& radial_pressure_face_flux) const noexcept
  {
    const IntVect face(AMREX_D_DECL(i, j, k));
    if (radial_metric.is_axis_face(i)) {
      return;
    }
    const IntVect direction(AMREX_D_DECL(
        direction_components[0], direction_components[1],
        direction_components[2]));
    const auto indicators = troubled_face_indicators(
        face, direction, normal_velocity_component, prims, true);
    // An optional fixed-width axis collar can add a direct adjacent-state
    // jump detector to the generic shock/contact selector.  A width of zero
    // disables this extra collar; the ordinary troubled-face fallback and
    // its parity-synchronised axis closure remain active.  This avoids
    // imposing a permanent low/high-order interface in otherwise resolved
    // near-axis flow.
    const int radial_face_offset = i - radial_metric.axis_face_index;
    const bool axis_collar_troubled =
        ComputeRzPairedPressure && radial_metric.active &&
        radial_metric.origin == Real(0.0) &&
        radial_face_offset >= 1 &&
        radial_face_offset <= axis_collar_faces &&
        axis_collar_fallback_needed(
            face, direction, prims, pressure_states,
            pressure_component, density_jump_threshold);
    if (!axis_collar_troubled &&
        !local_llf_fallback_needed(
            indicators, pressure_sensor_threshold,
            density_sensor_threshold, density_jump_threshold,
            compression_threshold)) {
      return;
    }
    write_first_order_llf_flux<ComputeRzPairedPressure>(
        face, normal_velocity_component - cls_t::QRHO - 1,
        prims, flux, closure, radial_metric, pressure_states,
        pressure_component, radial_pressure_face_flux);
  }

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  template <typename MarkerArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE IntVect marker_index_with_axis_parity(
      IntVect point, const MarkerArray&,
      const skew_detail::radial_metric_t radial_metric) const noexcept
  {
    if (radial_metric.active && radial_metric.origin == Real(0.0) &&
        point[0] < radial_metric.axis_face_index) {
      point[0] = 2 * radial_metric.axis_face_index - 1 - point[0];
    }
    return point;
  }

  template <typename MarkerArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool gp_state_usable(
      const IntVect& point, const MarkerArray& marker,
      const skew_detail::radial_metric_t radial_metric) const noexcept
  {
    const IntVect marker_point =
        marker_index_with_axis_parity(point, marker, radial_metric);
    return ibm_flux::is_usable(marker_point, marker);
  }

  template <typename MarkerArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool gp_cell_is_solid(
      const IntVect& point, const MarkerArray& marker,
      const skew_detail::radial_metric_t radial_metric) const noexcept
  {
    const IntVect marker_point =
        marker_index_with_axis_parity(point, marker, radial_metric);
    return ibm_flux::is_solid(marker_point, marker);
  }

  template <typename MarkerArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool gp_stencil_usable(
      const IntVect& face, const IntVect& direction,
      const int first_offset, const int count,
      const MarkerArray& marker,
      const skew_detail::radial_metric_t radial_metric) const noexcept
  {
    for (int n = 0; n < count; ++n) {
      if (!gp_state_usable(
              face + (first_offset + n) * direction,
              marker, radial_metric)) {
        return false;
      }
    }
    return true;
  }

  template <typename MarkerArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool gp_stencil_all_fluid(
      const IntVect& face, const IntVect& direction,
      const int first_offset, const int count,
      const MarkerArray& marker,
      const skew_detail::radial_metric_t radial_metric) const noexcept
  {
    for (int n = 0; n < count; ++n) {
      if (gp_cell_is_solid(
              face + (first_offset + n) * direction,
              marker, radial_metric)) {
        return false;
      }
    }
    return true;
  }

  template <bool ComputeRzPairedPressure,
            typename PressureArray, typename PressureFluxArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void apply_local_llf_fallback_ibm(
      const int i, const int j, const int k,
      const int normal_velocity_component,
      const GpuArray<int, 3>& direction_components,
      const Array4<Real const>& prims, const Array4<Real>& flux,
      const cls_t* closure, const Array4<uint8_t>& marker,
      const skew_detail::radial_metric_t radial_metric,
      const Real pressure_sensor_threshold,
      const Real density_sensor_threshold,
      const Real density_jump_threshold,
      const Real compression_threshold,
      const int axis_collar_faces,
      const PressureArray& pressure_states,
      const int pressure_component,
      const PressureFluxArray& radial_pressure_face_flux) const noexcept
  {
    const IntVect face(AMREX_D_DECL(i, j, k));
    if (radial_metric.is_axis_face(i)) {
      return;
    }
    const IntVect direction(AMREX_D_DECL(
        direction_components[0], direction_components[1],
        direction_components[2]));
    const IntVect left = face - direction;
    if (ibm_flux::is_solid_solid_face(face, direction, marker)) {
      return;
    }
    if (!gp_state_usable(left, marker, radial_metric) ||
        !gp_state_usable(face, marker, radial_metric)) {
      return;
    }

    const bool four_point_sensor_available =
        gp_state_usable(left - direction, marker, radial_metric) &&
        gp_state_usable(face + direction, marker, radial_metric);
    const auto indicators = troubled_face_indicators(
        face, direction, normal_velocity_component, prims,
        four_point_sensor_available);
    // Use the same optional R-Z collar for shared-GP states.  Both adjacent
    // states have already passed the marker-aware usability check above; no
    // deep-solid state, cut-cell geometry, or face fraction is consulted.
    const int radial_face_offset = i - radial_metric.axis_face_index;
    const bool axis_collar_troubled =
        ComputeRzPairedPressure && radial_metric.active &&
        radial_metric.origin == Real(0.0) &&
        radial_face_offset >= 1 &&
        radial_face_offset <= axis_collar_faces &&
        axis_collar_fallback_needed(
            face, direction, prims, pressure_states,
            pressure_component, density_jump_threshold);
    if (!axis_collar_troubled &&
        !local_llf_fallback_needed(
            indicators, pressure_sensor_threshold,
            density_sensor_threshold, density_jump_threshold,
            compression_threshold)) {
      return;
    }
    write_first_order_llf_flux<ComputeRzPairedPressure>(
        face, normal_velocity_component - cls_t::QRHO - 1,
        prims, flux, closure, radial_metric, pressure_states,
        pressure_component, radial_pressure_face_flux);
  }
#endif

  template <bool ComputeRzPairedPressure,
            typename FaceFluxArray, typename PressureFluxArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void mark_flux_pair_invalid(
      const IntVect& face, const FaceFluxArray& flux,
      const skew_detail::radial_metric_t radial_metric,
      const PressureFluxArray& radial_pressure_face_flux) const noexcept
  {
    const Real invalid = std::numeric_limits<Real>::quiet_NaN();
    for (int n = 0; n < cls_t::NCONS; ++n) {
      flux(face, n) = invalid;
    }
    if constexpr (ComputeRzPairedPressure) {
      if (radial_metric.active) {
        radial_pressure_face_flux(face, 0) = invalid;
      }
    }
  }

  template <typename PressureArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real axis_pressure_auxiliary_flux(
      const IntVect& face, const IntVect& radial_direction,
      const PressureArray& pressure_states,
      const int pressure_component,
      const int radial_domain_high_index) const noexcept
  {
    const Real pressure_ring1 =
        pressure_states(face, pressure_component);
    if constexpr (order == 2) {
      return pressure_ring1;
    }
    if (face[0] + 1 > radial_domain_high_index) {
      return pressure_ring1;
    }
    const Real pressure_ring2 =
        pressure_states(face + radial_direction, pressure_component);
    if constexpr (order == 4) {
      return (Real(7.0) * pressure_ring1 - pressure_ring2) /
             Real(6.0);
    }
    if (face[0] + 2 > radial_domain_high_index) {
      return (Real(7.0) * pressure_ring1 - pressure_ring2) /
             Real(6.0);
    }
    const Real pressure_ring3 =
        pressure_states(face + 2 * radial_direction, pressure_component);
    return (Real(37.0) * pressure_ring1 -
            Real(8.0) * pressure_ring2 + pressure_ring3) /
           Real(30.0);
  }

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  template <typename PressureArray, typename MarkerArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
  axis_pressure_auxiliary_flux_gp(
      const IntVect& face, const IntVect& radial_direction,
      const PressureArray& pressure_states,
      const int pressure_component, const MarkerArray& marker,
      const skew_detail::radial_metric_t radial_metric,
      bool& valid) const noexcept
  {
    const bool ring1_available =
        face[0] <= radial_metric.domain_high_index &&
        gp_state_usable(face, marker, radial_metric);
    if (!ring1_available) {
      valid = false;
      return Real(0.0);
    }
    const Real pressure_ring1 =
        pressure_states(face, pressure_component);
    if (!std::isfinite(pressure_ring1)) {
      valid = false;
      return Real(0.0);
    }

    Real pressure_face = pressure_ring1;
    if constexpr (order >= 4) {
      const IntVect ring2 = face + radial_direction;
      const bool ring2_available =
          ring2[0] <= radial_metric.domain_high_index &&
          gp_state_usable(ring2, marker, radial_metric);
      if (ring2_available) {
        const Real pressure_ring2 =
            pressure_states(ring2, pressure_component);
        if (!std::isfinite(pressure_ring2)) {
          valid = false;
          return Real(0.0);
        }
        pressure_face =
            (Real(7.0) * pressure_ring1 - pressure_ring2) /
            Real(6.0);
        if constexpr (order >= 6) {
          const IntVect ring3 = face + 2 * radial_direction;
          const bool ring3_available =
              ring3[0] <= radial_metric.domain_high_index &&
              gp_state_usable(ring3, marker, radial_metric);
          if (ring3_available) {
            const Real pressure_ring3 =
                pressure_states(ring3, pressure_component);
            if (!std::isfinite(pressure_ring3)) {
              valid = false;
              return Real(0.0);
            }
            pressure_face =
                (Real(37.0) * pressure_ring1 -
                 Real(8.0) * pressure_ring2 + pressure_ring3) /
                Real(30.0);
          }
        }
      }
    }
    return pressure_face;
  }

  template <typename PrimitiveArray, typename FaceFluxArray,
            typename MarkerArray, typename PressureArray,
            typename PressureFluxArray>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  one_sided_paired_pressure_flux(
      const IntVect& face, const IntVect& direction,
      const PrimitiveArray& prims, const FaceFluxArray& flux,
      const MarkerArray& marker, const cls_t* closure,
      const skew_detail::radial_metric_t radial_metric,
      const PressureArray& pressure_states,
      const int pressure_component,
      const PressureFluxArray& radial_pressure_face_flux) const noexcept
  {
    const Real face_radius = radial_metric.face_radius(face[0]);
    if (!(face_radius > Real(0.0)) || !std::isfinite(face_radius)) {
      return false;
    }

    int start = 99;
    int count = 0;
    if (gp_stencil_usable(face, direction, -2, 4, marker, radial_metric)) {
      start = -2;
      count = 4;
    } else if (gp_stencil_usable(
                   face, direction, -1, 4, marker, radial_metric)) {
      start = -1;
      count = 4;
    } else if (gp_stencil_usable(
                   face, direction, -3, 4, marker, radial_metric)) {
      start = -3;
      count = 4;
    } else if (gp_stencil_usable(
                   face, direction, -1, 3, marker, radial_metric)) {
      start = -1;
      count = 3;
    } else if (gp_stencil_usable(
                   face, direction, -2, 3, marker, radial_metric)) {
      start = -2;
      count = 3;
    } else {
      return false;
    }

    Real weight[4] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};
    if (count == 4 && start == -2) {
      weight[0] = Real(-1.0 / 16.0);
      weight[1] = Real(9.0 / 16.0);
      weight[2] = Real(9.0 / 16.0);
      weight[3] = Real(-1.0 / 16.0);
    } else if (count == 4 && start == -1) {
      weight[0] = Real(5.0 / 16.0);
      weight[1] = Real(15.0 / 16.0);
      weight[2] = Real(-5.0 / 16.0);
      weight[3] = Real(1.0 / 16.0);
    } else if (count == 4) {
      weight[0] = Real(1.0 / 16.0);
      weight[1] = Real(-5.0 / 16.0);
      weight[2] = Real(15.0 / 16.0);
      weight[3] = Real(5.0 / 16.0);
    } else if (start == -1) {
      weight[0] = Real(3.0 / 8.0);
      weight[1] = Real(3.0 / 4.0);
      weight[2] = Real(-1.0 / 8.0);
    } else {
      weight[0] = Real(-1.0 / 8.0);
      weight[1] = Real(3.0 / 4.0);
      weight[2] = Real(3.0 / 8.0);
    }

    Real face_value[cls_t::NCONS] = {};
    Real cell_flux[cls_t::NCONS];
    Real pressure_face = Real(0.0);
    for (int m = 0; m < count; ++m) {
      const IntVect cell = face + (start + m) * direction;
      closure->prims2flux(cell, 0, prims, cell_flux);
      cell_flux[cls_t::UMX] -= prims(cell, cls_t::QPRES);
      const Real metric_scale =
          radial_metric.cell_to_stored_flux_scale(cell[0], face[0]);
      for (int n = 0; n < cls_t::NCONS; ++n) {
        face_value[n] += weight[m] * metric_scale * cell_flux[n];
      }
      pressure_face +=
          weight[m] * pressure_states(cell, pressure_component);
    }
    if (!std::isfinite(pressure_face)) {
      mark_flux_pair_invalid<true>(
          face, flux, radial_metric, radial_pressure_face_flux);
      return true;
    }
    for (int n = 0; n < cls_t::NCONS; ++n) {
      flux(face, n) = face_value[n];
    }
    flux(face, cls_t::UMX) += pressure_face;
    radial_pressure_face_flux(face, 0) = pressure_face;
    return true;
  }
#endif


  // ............................................................. 
  // compute flux in each direction at i-1/2 //  
  // skew-symmetric formulation following f = U*V
  // U vector of conservative  vars (rho, rho ux, rho uy, rho uz rho e) //
  // V velocity vector               (ux,uy,uz)    //
  // fi = 1/2 ( U + Ui-1) * 1/2 *(V + Vi-1)  (example of second order)
  template <bool ComputeRzPairedPressure>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<const Real>& prims,
    const Array4<Real>& /*lambda_max*/, const Array4<Real>& flx, const cls_t* /*cls*/,
    const skew_detail::radial_metric_t radial_metric,
    const Array4<const Real>& pressure_states,
    const int pressure_component,
    const Array4<Real>& radial_pressure_face_flux) const {
    Real V[order],P[order];
    Real U[order][cls_t::NCONS];

    // prepare arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   

    for (int l = 0; l < order; l++) {  
      V[l] = prims(il,jl,kl,Qdir);       
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {U[l][nvar] = cons(il,jl,kl,nvar);}
      if constexpr (ComputeRzPairedPressure) {
        P[l] = radial_metric.active
                   ? pressure_states(il, jl, kl, pressure_component)
                   : prims(il, jl, kl, cls_t::QPRES);
      } else {
        P[l] = prims(il, jl, kl, cls_t::QPRES);
      }
      U[l][cls_t::UET] += P[l];
      if constexpr (ComputeRzPairedPressure) {
        if (radial_metric.active) {
          U[l][cls_t::UET] +=
              prims(il, jl, kl, cls_t::QPRES) - P[l];
        }
      }
      if (radial_metric.active) {
        const Real metric_scale =
            radial_metric.cell_to_stored_flux_scale(il, i);
        for (int nvar = 0; nvar < cls_t::NCONS; ++nvar) {
          U[l][nvar] *= metric_scale;
        }
        if constexpr (!ComputeRzPairedPressure) {
          P[l] *= metric_scale;
        }
      }
      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
    }

    // compute fluxes
    for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
      flx(i, j, k, nvar) =  0.0_rt;     
      for (int l = 0; l < order; l++) { 
        for (int m = 0; m < order; m++) { 
          flx(i, j, k, nvar) += coefskew(l,m)*U[l][nvar]*V[m];  
        } 
      }
    }  

    Real pressure_face = Real(0.0);
    for (int l = 0; l < order; l++) {
      pressure_face += P[l] * coefP(l);
    }
    if constexpr (ComputeRzPairedPressure) {
      if (radial_metric.active) {
        const IntVect face(AMREX_D_DECL(i, j, k));
        const IntVect radial_direction(
            AMREX_D_DECL(vdir[0], vdir[1], vdir[2]));
        if (radial_metric.is_axis_face(i)) {
          pressure_face = axis_pressure_auxiliary_flux(
              face, radial_direction, pressure_states,
              pressure_component, radial_metric.domain_high_index);
        }
        if (!std::isfinite(pressure_face)) {
          mark_flux_pair_invalid<ComputeRzPairedPressure>(
              face, flx, radial_metric, radial_pressure_face_flux);
          return;
        }
        radial_pressure_face_flux(face, 0) = pressure_face;
        flx(face, cls_t::UMX) = radial_metric.is_axis_face(i)
                                    ? Real(0.0)
                                    : flx(face, cls_t::UMX) + pressure_face;
        return;
      }
    }
    flx(i, j, k, Qdir - 1) += pressure_face;
 
  }  
  // .............................................................
#if (AMREX_USE_GPIBM || CNS_USE_EB)
  template <bool ComputeRzPairedPressure>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir_ibm(
      int i, int j, int k, int Qdir,
      const GpuArray<int, 3>& vdir, const Array4<Real>& /*cons*/,
      const Array4<const Real>& prims,
      const Array4<Real>& /*lambda_max*/, const Array4<Real>& flx,
      const cls_t* cls, const Array4<uint8_t>& marker,
      const skew_detail::radial_metric_t radial_metric,
      const Array4<const Real>& pressure_states,
      const int pressure_component,
      const Array4<Real>& radial_pressure_face_flux) const
  {
    Real velocity[order];
    Real pressure[order];
    Real state[order][cls_t::NCONS];

    const IntVect face(AMREX_D_DECL(i, j, k));
    const IntVect direction(
        AMREX_D_DECL(vdir[0], vdir[1], vdir[2]));
    const IntVect left = face - direction;
    const bool paired_radial =
        ComputeRzPairedPressure && radial_metric.active;

    if (gp_cell_is_solid(face, marker, radial_metric) &&
        gp_cell_is_solid(left, marker, radial_metric)) {
      ibm_flux::zero_flux<decltype(flx), cls_t>(face, flx);
      if constexpr (ComputeRzPairedPressure) {
        if (radial_metric.active) {
          radial_pressure_face_flux(face, 0) = Real(0.0);
        }
      }
      return;
    }

    const bool full_stencil = gp_stencil_usable(
        face, direction, -halfsten, order, marker, radial_metric);
    if (!full_stencil) {
      if constexpr (order >= 4) {
        if constexpr (ComputeRzPairedPressure) {
          if (radial_metric.active &&
              one_sided_paired_pressure_flux(
                  face, direction, prims, flx, marker, cls,
                  radial_metric, pressure_states, pressure_component,
                  radial_pressure_face_flux)) {
            return;
          }
        }
        if (!paired_radial &&
            ibm_flux::one_sided_polynomial_flux(
                face, Qdir - 1, prims, flx, marker, *cls,
                radial_metric.active, radial_metric.origin,
                radial_metric.spacing)) {
          return;
        }
      }

      if (!gp_state_usable(left, marker, radial_metric) ||
          !gp_state_usable(face, marker, radial_metric)) {
        mark_flux_pair_invalid<ComputeRzPairedPressure>(
            face, flx, radial_metric, radial_pressure_face_flux);
        return;
      }

      IntVect point = left;
      for (int l = 0; l < 2; ++l) {
        velocity[l] = prims(point, Qdir);
        cls->prims2cons(point, prims, state[l]);
        state[l][cls_t::UET] += prims(point, cls_t::QPRES);
        pressure[l] = paired_radial
                          ? pressure_states(point, pressure_component)
                          : prims(point, cls_t::QPRES);
        if (radial_metric.active) {
          const Real metric_scale =
              radial_metric.cell_to_stored_flux_scale(point[0], i);
          for (int n = 0; n < cls_t::NCONS; ++n) {
            state[l][n] *= metric_scale;
          }
          if (!paired_radial) {
            pressure[l] *= metric_scale;
          }
        }
        point += direction;
      }

      for (int n = 0; n < cls_t::NCONS; ++n) {
        flx(face, n) = Real(0.25) *
                       (state[0][n] + state[1][n]) *
                       (velocity[0] + velocity[1]);
      }
      Real pressure_face = Real(0.5) * (pressure[0] + pressure[1]);
      if constexpr (ComputeRzPairedPressure) {
        if (radial_metric.active) {
          if (radial_metric.is_axis_face(i)) {
            bool valid = true;
            pressure_face = axis_pressure_auxiliary_flux_gp(
                face, direction, pressure_states, pressure_component,
                marker, radial_metric, valid);
            if (!valid) {
              mark_flux_pair_invalid<ComputeRzPairedPressure>(
                  face, flx, radial_metric, radial_pressure_face_flux);
              return;
            }
          }
          radial_pressure_face_flux(face, 0) = pressure_face;
          flx(face, cls_t::UMX) = radial_metric.is_axis_face(i)
                                      ? Real(0.0)
                                      : flx(face, cls_t::UMX) + pressure_face;
          return;
        }
      }
      flx(face, Qdir - 1) += pressure_face;
      return;
    }

    IntVect point = face - halfsten * direction;
    for (int l = 0; l < order; ++l) {
      velocity[l] = prims(point, Qdir);
      cls->prims2cons(point, prims, state[l]);
      state[l][cls_t::UET] += prims(point, cls_t::QPRES);
      pressure[l] = paired_radial
                        ? pressure_states(point, pressure_component)
                        : prims(point, cls_t::QPRES);
      if (radial_metric.active) {
        const Real metric_scale =
            radial_metric.cell_to_stored_flux_scale(point[0], i);
        for (int n = 0; n < cls_t::NCONS; ++n) {
          state[l][n] *= metric_scale;
        }
        if (!paired_radial) {
          pressure[l] *= metric_scale;
        }
      }
      point += direction;
    }

    for (int n = 0; n < cls_t::NCONS; ++n) {
      flx(face, n) = Real(0.0);
      for (int l = 0; l < order; ++l) {
        for (int m = 0; m < order; ++m) {
          flx(face, n) +=
              coefskew(l, m) * state[l][n] * velocity[m];
        }
      }
    }
    Real pressure_face = Real(0.0);
    for (int l = 0; l < order; ++l) {
      pressure_face += pressure[l] * coefP(l);
    }
    if constexpr (ComputeRzPairedPressure) {
      if (radial_metric.active) {
        if (radial_metric.is_axis_face(i)) {
          bool valid = true;
          pressure_face = axis_pressure_auxiliary_flux_gp(
              face, direction, pressure_states, pressure_component,
              marker, radial_metric, valid);
          if (!valid) {
            mark_flux_pair_invalid<ComputeRzPairedPressure>(
                face, flx, radial_metric, radial_pressure_face_flux);
            return;
          }
        }
        radial_pressure_face_flux(face, 0) = pressure_face;
        flx(face, cls_t::UMX) = radial_metric.is_axis_face(i)
                                    ? Real(0.0)
                                    : flx(face, cls_t::UMX) + pressure_face;
        return;
      }
    }
    flx(face, Qdir - 1) += pressure_face;
  }
#endif
  // .............................................................
#if (AMREX_USE_GPIBM || CNS_USE_EB)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void fluxdissip_dir_ibm(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& /*cons*/,
    const Array4<Real const>& prims, const Array4<Real>& /*lambda*/, const Array4<Real>& flx,
    const cls_t* cls,const Array4<uint8_t>& marker,
    const skew_detail::radial_metric_t radial_metric) const {

    const IntVect iv(AMREX_D_DECL(i, j, k));
    const IntVect ivd(AMREX_D_DECL(vdir[0], vdir[1], vdir[2]));
    const IntVect ivl = iv - ivd;

    // Away from the axis the R-Z divergence supplies the face-radius metric,
    // so JST must dissipate jumps/differences of the physical conservative
    // state U, not of rU.  At r=0 the face slot stores an auxiliary metric
    // central flux; a physical JST face contribution has zero metric area and
    // must not be mixed into that slot.
    if (radial_metric.is_axis_face(i)) return;

    if (gp_cell_is_solid(iv, marker, radial_metric) &&
        gp_cell_is_solid(ivl, marker, radial_metric)) return;
    if (!gp_state_usable(ivl, marker, radial_metric) ||
        !gp_state_usable(iv, marker, radial_metric)) return;

    // The sensor reads offsets [-2,+1].  For Skew2 that is wider than the
    // convective footprint, so both footprints must pass the marker audit.
    const bool full_diss_stencil =
        gp_stencil_usable(
            iv, ivd, -halfsten, order, marker, radial_metric) &&
        gp_stencil_usable(iv, ivd, -2, 4, marker, radial_metric);

    Real sen_num = Real(0.0);
    Real sen_denom = Real(1.0e-16);
    Real sen = Real(0.0);
    if (full_diss_stencil) {
      for (int l = 0; l < NVARSEN; ++l) {
        const int nv = NSEN[l];
        const Real p0 = prims(iv - 2 * ivd, nv);
        const Real p1 = prims(ivl, nv);
        const Real p2 = prims(iv, nv);
        const Real p3 = prims(iv + ivd, nv);
        // disconSensor(pp, pl, pr) centres the second difference on its FIRST
        // argument, so the two face-adjacent sensors are centred on p1 and p2.
        // Passing (p0,p1,p2) evaluates 2(p2 - 2p0 + p1), which is non-zero for
        // a linear profile and fires C2 on smooth gradients.
        const Real local = amrex::max(disconSensor(p1, p0, p2),
                                      disconSensor(p2, p1, p3));
        sen_num += local * local;
        sen_denom += local;
      }
      sen = sen_num / sen_denom;
    } else {
      const bool left_triplet =
          gp_state_usable(iv - 2 * ivd, marker, radial_metric) &&
          gp_state_usable(ivl, marker, radial_metric) &&
          gp_state_usable(iv, marker, radial_metric);
      const bool right_triplet =
          gp_state_usable(ivl, marker, radial_metric) &&
          gp_state_usable(iv, marker, radial_metric) &&
          gp_state_usable(iv + ivd, marker, radial_metric);
      if (left_triplet || right_triplet) {
        const IntVect first = left_triplet ? iv - 2 * ivd : ivl;
        for (int l = 0; l < NVARSEN; ++l) {
          const int nv = NSEN[l];
          // centre of the triplet is first + ivd, so it goes first
          const Real local = disconSensor(prims(first + ivd, nv),
                                          prims(first, nv),
                                          prims(first + 2 * ivd, nv));
          sen_num += local * local;
          sen_denom += local;
        }
        sen = sen_num / sen_denom;
      } else {
        // No safe three-point sensor exists in a one-cell fluid gap.
        sen = Real(1.0);
      }
    }

    const Real rr = amrex::max(
        std::abs(prims(ivl, Qdir)) + prims(ivl, cls_t::QC),
        std::abs(prims(iv, Qdir)) + prims(iv, cls_t::QC));
    Real sen_eff = sen;
    for (int q = 1; q < sensor_power; ++q) { sen_eff *= sen; }
    const Real eps2 = Cshock * rr * sen_eff;
    const Real eps4 = amrex::max(Real(0.0), Cdamp * rr - eps2);

    if (!full_diss_stencil) {
      Real right_state[cls_t::NCONS];
      Real left_state[cls_t::NCONS];
      cls->prims2cons(iv, prims, right_state);
      cls->prims2cons(ivl, prims, left_state);
      for (int nvar = 0; nvar < cls_t::NCONS; ++nvar) {
        flx(iv, nvar) -= eps2 * (right_state[nvar] - left_state[nvar]);
      }
      return;
    }

    IntVect point = iv - halfsten * ivd;
    for (int l = 0; l < order; ++l) {
      Real point_state[cls_t::NCONS];
      cls->prims2cons(point, prims, point_state);
      for (int nvar = 0; nvar < cls_t::NCONS; ++nvar) {
        flx(iv, nvar) -= eps2 * coefshock(l) * point_state[nvar];
        flx(iv, nvar) += eps4 * coefdamp(l) * point_state[nvar];
      }
      point += ivd;
    }
  }
#endif


  // .............................................................
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void fluxdissip_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<const Real>& prims, const Array4<Real>& /* lambda */, const Array4<Real>& flx,
    const cls_t* /*cls*/,
    const skew_detail::radial_metric_t radial_metric) const {

    // The radial face metric is applied by the R-Z divergence.  Dissipating
    // rU here would leave a deterministic pseudo-source whenever the JST
    // coefficient changes between neighbouring faces.  The axis face itself
    // has zero physical metric area, while its array slot is reserved for the
    // auxiliary metric central/paired-pressure closure, so JST adds nothing
    // there.
    if (radial_metric.is_axis_face(i)) return;

    int ir = i+vdir[0];   int jr = j+vdir[1];   int kr = k+vdir[2];   
    int il = i-vdir[0];   int jl = j-vdir[1];   int kl = k-vdir[2];   
    int ill= il-vdir[0]; int jll = jl-vdir[1]; int kll= kl-vdir[2];   
 
    // const int idir = Qdir -1;

    // calculate sensor    
    Real p0,p1,p2,p3;
    Real sen_num= Real(0.0),sen_denom=Real(1.0e-16);
    // loop over sensor variables
    int nv = 0;
    Real sen = Real(0.0);
    for (int l=0;l<NVARSEN;l++)
    { 
      nv = NSEN[l];       
      p0 =  prims(ill,jll,kll,nv);
      p1 =  prims(il,jl,kl,nv);
      p2 =  prims(i,j,k,nv);
      p3 =  prims(ir,jr,kr,nv);    
      // centred second differences: sensor at p1 and at p2 (see the note above)
      sen  = std::max(disconSensor(p1,p0,p2), disconSensor(p2,p1,p3) );
      sen_num += sen*sen;sen_denom +=sen;
      sen = sen_num/sen_denom;
    }

    // reduce order close to BC by making sensor  = 1   
    // sen = (i < mask_sen(idir,1)) ? 1.0 : sen;  
    // sen = (i > mask_sen(idir,2)) ? 1.0 : sen;  
    
    // spectral radius Jacobian matrix (u + c)
    //Real rr = std::max(lambda(il, jl, kl, 0), lambda(i, j, k, 0));
    const Real rr = amrex::max(
        std::abs(prims(il, jl, kl, Qdir)) +
            prims(il, jl, kl, cls_t::QC),
        std::abs(prims(i, j, k, Qdir)) +
            prims(i, j, k, cls_t::QC));
    // Raising the sensor to sensor_power sharpens its smooth-flow decay
    // (sen = O(dx^2) per power) so the shock term does not limit the formal
    // order of the base scheme.
    Real sen_eff = sen;
    for (int q = 1; q < sensor_power; ++q) { sen_eff *= sen; }
    Real eps2 = Cshock*rr*sen_eff;
    Real eps4 = std::max(0.0, Cdamp*rr - eps2);

    // shock capturing (sensor-scaled second-difference family)
    int ii= i-halfsten*vdir[0]; int jj= j-halfsten*vdir[1]; int kk= k-halfsten*vdir[2];
    for (int l = 0; l < order; l++) {
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
        flx(i, j, k, nvar) -= eps2*coefshock(l)*cons(ii,jj,kk,nvar);
      }
      ii +=  vdir[0];jj +=  vdir[1];kk +=  vdir[2];
    }

    // Skew2 historically has no background damping term.  Higher-order
    // variants use the decoupled third- or fifth-difference stencil.
    if constexpr (order > 2) {
      ii= i-damp_halfsten*vdir[0]; jj= j-damp_halfsten*vdir[1]; kk= k-damp_halfsten*vdir[2];
      for (int l = 0; l < damp_width; l++) {
        for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
          flx(i, j, k, nvar) += eps4*coefdamp_ext(l)*cons(ii,jj,kk,nvar);
        }
        ii +=  vdir[0];jj +=  vdir[1];kk +=  vdir[2];
      }
    }
  }
  // .............................................................  

  
  // coefficents
  //static const int ordermax = 6;
  typedef Array2D<Real, 0, order, 0, order> arrCoeff_t;
  arrCoeff_t coefskew;
  typedef Array1D<Real, 0, order> arrayNumCoef;
  arrayNumCoef coefdamp,coefshock,coefP;
  // Damping coefficients on the decoupled stencil (width damp_width, gathered
  // from face_index - damp_halfsten).  Sign convention matches coefdamp:
  // the stored stencil evaluates minus the (damp_difference_order)-th
  // difference at the face.
  Array1D<Real, 0, 6> coefdamp_ext;

  int halfsten = order / 2;

// sensor variables (density and pressure by default)
// #if NUM_SPECIES > 1    
//     static constexpr int NVARSEN = 2 + NUM_SPECIES;
// #else
    static  constexpr int NVARSEN=2;
//#endif    
    int NSEN[NVARSEN];

  // masking sensor  
  typedef Array2D<int, 0, AMREX_SPACEDIM, 0, 2> arrIntCoeff_t;
  arrIntCoeff_t mask_sen;

  //-----------------------------------------------------------------------------------
  };


#endif
