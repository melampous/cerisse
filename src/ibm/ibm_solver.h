#ifndef IBM_SOLVER_H_
#define IBM_SOLVER_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Scan.H>
#include <AMReX_Reduce.H>
#include <AMReX_iMultiFab.H>
#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#endif

#include "ibm_containers.h"
#include "ibm_method_config.h"
#include "ibm_rz_moments.h"
#include "IBMSharedGPFluxUtils.h"
#include "RZFiniteVolume.h"

#include <array>
#include <cmath>    // std::llround in surface-force reduction
#include <cstring>  // std::memcpy in surface MPI packing
#include <functional>
#include <fstream>
#include <iomanip>  // std::setprecision, std::fixed
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>

//===================================================================================
///-------------------------------- main class --------------------------------------
///
/// \brief ibm_solver_t is explicit geometry (triangulation based) immersed boundary method
/// class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
/// the geometry
///
template <typename wallmodel, typename param, typename cls_t>
class ibm_solver_t
{
public:
  // constant factor for image point
  static constexpr int  iorder_tparm = param::interp_order; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm = param::extrap_order; // number of image points used for ghost point extrapolation
  static_assert(param::interp_order == 1 || param::interp_order == 2,
      "IBM interp_order: 1 = bi/trilinear, 2 = WLS quadratic; no other value is implemented");
  static_assert(param::interp_order_surf == 1 || param::interp_order_surf == 2,
      "IBM interp_order_surf: 1 = bi/trilinear, 2 = WLS quadratic; no other value is implemented");
  static constexpr Real cim = param::alpha;
  static_assert(cim > Real(0.0), "IBM alpha must be positive");
  
  static constexpr int  iorder_tparm_surf = param::interp_order_surf; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm_surf = param::extrap_order_surf; // number of image points used for surface reconstruction
  static constexpr Real cim_surf  = param::alpha_surf;
  static_assert(cim_surf > Real(0.0), "IBM alpha_surf must be positive");
  static constexpr bool use_ns_noslip_pressure_compatibility =
      ibm_detail::uses_navier_stokes_noslip_pressure_compatibility<wallmodel>();
  static constexpr ibm_pressure_closure_t wall_pressure_closure = [] {
    if constexpr (requires { wallmodel::pressure_closure; }) {
      return wallmodel::pressure_closure;
    }
    return ibm_pressure_closure_t::zero_gradient;
  }();
  static constexpr bool supports_complete_viscous_surface_traction = [] {
    if constexpr (requires { wallmodel::stationary_no_slip; }) {
      return wallmodel::stationary_no_slip;
    }
    return false;
  }();
  static constexpr bool sharp_feature_pressure_limiter = [] {
    if constexpr (requires { param::sharp_feature_pressure_limiter; }) {
      return bool(param::sharp_feature_pressure_limiter);
    }
    return false;
  }();
  static constexpr bool sharp_feature_tangential_reconstruction = [] {
    if constexpr (requires { param::sharp_feature_tangential_reconstruction; }) {
      return bool(param::sharp_feature_tangential_reconstruction);
    }
    return false;
  }();
  static constexpr Real sharp_feature_angle_deg = [] {
    if constexpr (requires { param::sharp_feature_angle_deg; }) {
      return Real(param::sharp_feature_angle_deg);
    }
    return Real(30.0);
  }();
  static constexpr Real sharp_feature_radius_cells = [] {
    if constexpr (requires { param::sharp_feature_radius_cells; }) {
      return Real(param::sharp_feature_radius_cells);
    }
    return Real(3.0);
  }();
  static constexpr Real sharp_feature_tangential_radius_cells = [] {
    if constexpr (requires { param::sharp_feature_tangential_radius_cells; }) {
      return Real(param::sharp_feature_tangential_radius_cells);
    }
    return Real(3.0);
  }();
  // Support visibility is an initialization-time geometry check:
  //   0 = disabled (legacy, bit-identical interpolation stencils)
  //   1 = audit occluded fluid supports but keep their weights
  //   2 = audit and remove occluded supports before the WLS fit
  static constexpr int support_visibility_mode = [] {
    if constexpr (requires { param::support_visibility_mode; }) {
      return int(param::support_visibility_mode);
    }
    return 0;
  }();
  static constexpr Real support_visibility_condition_max = [] {
    if constexpr (requires { param::support_visibility_condition_max; }) {
      return Real(param::support_visibility_condition_max);
    }
    return Real(1.0e10);
  }();
  static constexpr ibm_shared_gp_reconstruction_t
      shared_gp_reconstruction = [] {
    if constexpr (requires { param::shared_gp_reconstruction; }) {
      return param::shared_gp_reconstruction;
    }
    return ibm_shared_gp_reconstruction_t::image_point;
  }();
  static constexpr bool use_bi_constrained_shared_gp =
      shared_gp_reconstruction ==
      ibm_shared_gp_reconstruction_t::boundary_intercept_constrained;
  static constexpr bool expanded_bic_support = [] {
    if constexpr (requires { param::expanded_bic_support; }) {
      return bool(param::expanded_bic_support);
    }
    return false;
  }();
  // Opt-in finite-volume correction for stationary-wall normal momentum. The
  // BI-constrained polynomial consumes support-cell averages and returns the
  // Cartesian ghost-cell average of rho*u_n.
  static constexpr bool cell_average_bic_normal_momentum = [] {
    if constexpr (requires { param::cell_average_bic_normal_momentum; }) {
      return bool(param::cell_average_bic_normal_momentum);
    }
    return false;
  }();
  // Opt-in auxiliary entropy extension.  No Euler wall condition prescribes
  // d(sigma)/dn, so sigma is extrapolated from a complete one-sided fluid jet
  // instead of imposing a homogeneous Neumann condition.
  static constexpr bool one_sided_entropy_jet = [] {
    if constexpr (requires { param::one_sided_entropy_jet; }) {
      return bool(param::one_sided_entropy_jet);
    }
    return false;
  }();
  // At an R-Z wall-axis junction, a regular scalar is even in signed radius.
  // The first-ring entropy jet therefore uses an axis-regular quadratic basis
  // instead of fitting an unphysical linear-r component from r>=0 supports.
  static constexpr bool rz_axis_regular_entropy_jet = [] {
    if constexpr (requires { param::rz_axis_regular_entropy_jet; }) {
      return bool(param::rz_axis_regular_entropy_jet);
    }
    return false;
  }();
  // R-Z production data semantics. Fluid supports are recovered as pointwise
  // centre states from annular conservative cell averages, while the stored
  // BI-CWLS target is the full annular ghost-cell average. This option is
  // deliberately compile-time so the Cartesian path remains unchanged.
  static constexpr bool rz_annular_bic_cell_average = [] {
    if constexpr (requires { param::rz_annular_bic_cell_average; }) {
      return bool(param::rz_annular_bic_cell_average);
    }
    return false;
  }();
  // Cartesian no-slip finite-volume semantics. Fluid support values are
  // recovered from conservative cell averages before BI-CWLS evaluation, and
  // the unique solid-side GP state is stored as a conservative cell average.
  // This remains a pure shared-GP/full-Cartesian method.
  static constexpr bool ns_cartesian_bic_cell_average = [] {
    if constexpr (requires { param::ns_cartesian_bic_cell_average; }) {
      return bool(param::ns_cartesian_bic_cell_average);
    }
    return false;
  }();
  static constexpr bool bic_point_target_functionals =
      rz_annular_bic_cell_average ||
      ns_cartesian_bic_cell_average;
  // Recover meridional point states for BI-CWLS from the annular conservative
  // averages using same-fluid, quadratic finite-volume stencils. The bulk
  // R-Z WENO field retains its independently verified centre-state semantics;
  // this opt-in changes only the fluid support data seen by shared-GP BIC.
  static constexpr bool rz_bic_support_centre_recovery = [] {
    if constexpr (requires { param::rz_bic_support_centre_recovery; }) {
      return bool(param::rz_bic_support_centre_recovery);
    }
    return false;
  }();
  // Hierarchical consistency scale for the auxiliary extensions. Compare the
  // quadratic state with an embedded linear state on exactly the same visible
  // fluid support. The local support range is regularized by
  // sqrt(|P2-P1|*q_scale): the resulting ratio vanishes at ordinary smooth
  // points and smooth extrema, but remains O(1) at an unresolved jump.
  // Non-smooth support variation is handled by an independent sensor below.
  static constexpr Real bic_auxiliary_hierarchy_tolerance = [] {
    if constexpr (requires { param::bic_auxiliary_hierarchy_tolerance; }) {
      return Real(param::bic_auxiliary_hierarchy_tolerance);
    }
    return Real(0.10);
  }();
  // Entropy-jet smoothness scale in log(p/rho^gamma).  Smooth support data
  // have Delta sigma=O(h), so the quartic limiter below is asymptotically
  // inactive.  An unresolved shock or entropy layer has Delta sigma=O(1) and
  // reduces only the auxiliary entropy extension to the homogeneous-Neumann
  // BIC reference; pressure and velocity closures remain independent.
  static constexpr Real bic_entropy_jet_variation_tolerance = [] {
    if constexpr (requires {
                    param::bic_entropy_jet_variation_tolerance;
                  }) {
      return Real(param::bic_entropy_jet_variation_tolerance);
    }
    return Real(0.10);
  }();
  // Characteristic-speed-normalized variation used to disable the
  // finite-volume normal-momentum correction across unresolved wave layers.
  // The reference speed contains the local sound and tangential speeds, so a
  // smooth slip wall gives S_u=O(h) even though u_n itself vanishes at the BI.
  static constexpr Real bic_normal_momentum_variation_tolerance = [] {
    if constexpr (requires {
                    param::bic_normal_momentum_variation_tolerance;
                  }) {
      return Real(param::bic_normal_momentum_variation_tolerance);
    }
    return Real(0.10);
  }();
  // Optional fixed, smooth-wall geometry contract. The discrete polygon/STL
  // remains the topology carrier used by the BVH (inside/outside, first hit,
  // and owning element), while all boundary closure operations use one
  // problem-supplied boundary intercept and orthonormal local frame. This is
  // intended for analytic or spline-defined 2-D walls; curvature is never
  // inferred from adjacent facets.
  static constexpr bool smooth_wall_geometry = [] {
    if constexpr (requires { param::smooth_wall_geometry; }) {
      return bool(param::smooth_wall_geometry);
    }
    return false;
  }();
  static constexpr bool stationary_slip_wall = [] {
    if constexpr (requires { wallmodel::stationary_slip; }) {
      return bool(wallmodel::stationary_slip);
    }
    return false;
  }();
  static constexpr bool stationary_no_slip_wall = [] {
    if constexpr (requires { wallmodel::stationary_no_slip; }) {
      return bool(wallmodel::stationary_no_slip);
    }
    return false;
  }();
  // A problem-specific fixed mixed wall (for example, a stationary slip wall
  // containing a prescribed injection patch) must not claim stationary_slip:
  // that flag selects the generic slip closure and would bypass its custom
  // boundary state.  It may instead opt in to the full-cell positivity
  // baseline explicitly.  The default remains fail-closed.
  static constexpr bool positivity_compatible_fixed_wall = [] {
    if constexpr (requires {
                    wallmodel::positivity_compatible_fixed_wall;
                  }) {
      return bool(wallmodel::positivity_compatible_fixed_wall);
    }
    return stationary_slip_wall;
  }();
  static constexpr bool isothermal_wall = [] {
    if constexpr (requires { wallmodel::isothermal_wall; }) {
      return bool(wallmodel::isothermal_wall);
    }
    return false;
  }();
  static constexpr bool adiabatic_wall = [] {
    if constexpr (requires { wallmodel::adiabatic_wall; }) {
      return bool(wallmodel::adiabatic_wall);
    }
    return false;
  }();
  static constexpr bool entropy_wall_extension = [] {
    if constexpr (requires { wallmodel::entropy_extension; }) {
      return bool(wallmodel::entropy_extension);
    }
    return false;
  }();
  static constexpr Real wall_temperature = [] {
    if constexpr (requires { wallmodel::Twall; }) {
      return Real(wallmodel::Twall);
    }
    return Real(0.0);
  }();
  static constexpr bool uses_sharp_feature_annotation =
      sharp_feature_pressure_limiter ||
      sharp_feature_tangential_reconstruction;
  static constexpr bool uses_shock_aware_pressure_closure = [] {
    if constexpr (requires { wallmodel::pressure_closure; }) {
      return wallmodel::pressure_closure ==
             ibm_pressure_closure_t::shock_aware_fluid_extrapolation;
    }
    return false;
  }();
  static constexpr bool uses_zero_gradient_pressure_closure = [] {
    if constexpr (requires { wallmodel::pressure_closure; }) {
      return wallmodel::pressure_closure ==
             ibm_pressure_closure_t::zero_gradient;
    }
    return true;
  }();
  static constexpr bool uses_prescribed_gradient_pressure_closure = [] {
    if constexpr (requires { wallmodel::pressure_closure; }) {
      return wallmodel::pressure_closure ==
             ibm_pressure_closure_t::prescribed_gradient;
    }
    return false;
  }();
  static constexpr bool uses_euler_slip_curvature_pressure_closure = [] {
    if constexpr (requires { wallmodel::pressure_closure; }) {
      return wallmodel::pressure_closure ==
             ibm_pressure_closure_t::euler_slip_analytic_curvature;
    }
    return false;
  }();
  static constexpr bool uses_entropy_compatible_pressure_closure =
      uses_prescribed_gradient_pressure_closure ||
      uses_euler_slip_curvature_pressure_closure;
  static_assert(!sharp_feature_pressure_limiter || AMREX_SPACEDIM == 2,
      "sharp_feature_pressure_limiter is currently implemented only for 2-D polygons");
  static_assert(!sharp_feature_tangential_reconstruction ||
                    AMREX_SPACEDIM == 2,
      "sharp_feature_tangential_reconstruction is currently implemented only for 2-D polygons");
  static_assert(!sharp_feature_pressure_limiter || uses_shock_aware_pressure_closure,
      "sharp_feature_pressure_limiter requires shock_aware_fluid_extrapolation");
  static_assert(!uses_sharp_feature_annotation ||
                    (sharp_feature_angle_deg > Real(0.0) &&
                     sharp_feature_angle_deg < Real(180.0)),
      "sharp_feature_angle_deg must lie in (0, 180)");
  static_assert(!sharp_feature_pressure_limiter ||
                    sharp_feature_radius_cells > Real(0.0),
      "sharp_feature_radius_cells must be positive");
  static_assert(!sharp_feature_tangential_reconstruction ||
                    sharp_feature_tangential_radius_cells > Real(0.0),
      "sharp_feature_tangential_radius_cells must be positive");
  static_assert(!sharp_feature_tangential_reconstruction ||
                    (iorder_tparm == 2 && iorder_tparm_surf == 2),
      "sharp-feature tangential reconstruction requires quadratic WLS for "
      "both GP and surface image points");

  // Convert the backend-specific Point representation from the common IBM
  // coordinate array used by smooth-wall callbacks.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static Point surface_point_from_array(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& xyz) noexcept
  {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
    return Point(xyz(0), xyz(1));
#else
    return Point(xyz(0), xyz(1), xyz(2));
#endif
#else
    Point point{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) point[d] = xyz(d);
    return point;
#endif
  }

  // Apply the fixed smooth-wall local-frame callback. The caller first
  // obtains an owning element from the faceted BVH, then evaluates this frame
  // at the closure boundary intercept. A false/non-finite callback result is
  // a geometry-contract violation, not a reason to mix faceted and smooth
  // closure data silently.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void apply_volume_surface_frame_override(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& xyz,
      LocalFrame& localframe) noexcept
  {
    if constexpr (smooth_wall_geometry) {
      LocalFrame candidate = localframe;
      const bool valid = param::smooth_wall_local_frame(xyz, candidate);
      Real norm2 = Real(0.0);
      Real tangent2 = Real(0.0);
      Real orthogonality = Real(0.0);
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        norm2 += candidate.normal[d] * candidate.normal[d];
        tangent2 += candidate.tangent1[d] * candidate.tangent1[d];
        orthogonality += candidate.normal[d] * candidate.tangent1[d];
      }
      const bool finite = amrex::Math::isfinite(
          norm2 + tangent2 + orthogonality);
      AMREX_ASSERT(valid && finite && norm2 > Real(0.0) &&
                   tangent2 > Real(0.0));
      if (valid && finite && norm2 > Real(0.0) && tangent2 > Real(0.0)) {
        const Real inv_norm = amrex::Math::rsqrt(norm2);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          candidate.normal[d] *= inv_norm;
        }
        Real projection = Real(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          projection += candidate.tangent1[d] * candidate.normal[d];
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          candidate.tangent1[d] -= projection * candidate.normal[d];
        }
        tangent2 = Real(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          tangent2 += candidate.tangent1[d] * candidate.tangent1[d];
        }
        AMREX_ASSERT(tangent2 > Real(0.0) &&
                     amrex::Math::isfinite(tangent2));
        const Real inv_tangent = amrex::Math::rsqrt(tangent2);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          candidate.tangent1[d] *= inv_tangent;
        }
#if (AMREX_SPACEDIM == 3)
        candidate.tangent2[0] =
            candidate.normal[1] * candidate.tangent1[2] -
            candidate.normal[2] * candidate.tangent1[1];
        candidate.tangent2[1] =
            candidate.normal[2] * candidate.tangent1[0] -
            candidate.normal[0] * candidate.tangent1[2];
        candidate.tangent2[2] =
            candidate.normal[0] * candidate.tangent1[1] -
            candidate.normal[1] * candidate.tangent1[0];
#endif
        localframe = candidate;
      }
    } else {
      amrex::ignore_unused(xyz, localframe);
    }
  }

  // Keep the raw BVH hit immutable for topology and owning-element audit, but
  // evaluate every wall closure at one problem-supplied smooth boundary
  // intercept when the contract is enabled.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void apply_volume_surface_point_override(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& raw_point,
      Array1D<Real, 0, AMREX_SPACEDIM - 1>& closure_point) noexcept
  {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      closure_point(d) = raw_point(d);
    }
    if constexpr (smooth_wall_geometry) {
      Array1D<Real, 0, AMREX_SPACEDIM - 1> candidate = closure_point;
      const bool valid =
          param::smooth_wall_boundary_intercept(raw_point, candidate);
      bool finite = true;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        finite = finite && amrex::Math::isfinite(candidate(d));
      }
      AMREX_ASSERT(valid && finite);
      if (valid && finite) closure_point = candidate;
    }
  }

  static_assert(!smooth_wall_geometry || AMREX_SPACEDIM == 2,
      "The production smooth-wall geometry contract is currently restricted "
      "to fixed two-dimensional walls");
#ifdef CNS_USE_FSI
  static_assert(!smooth_wall_geometry,
      "The production smooth-wall geometry contract does not yet support FSI motion");
#endif
  static_assert(
      !smooth_wall_geometry ||
          requires(const Array1D<Real, 0, AMREX_SPACEDIM - 1>& raw,
                   Array1D<Real, 0, AMREX_SPACEDIM - 1>& point) {
            param::smooth_wall_boundary_intercept(raw, point);
          },
      "smooth_wall_geometry requires "
      "param::smooth_wall_boundary_intercept(raw_point, boundary_intercept)");
  static_assert(
      !smooth_wall_geometry ||
          requires(const Array1D<Real, 0, AMREX_SPACEDIM - 1>& point,
                   LocalFrame& frame) {
            param::smooth_wall_local_frame(point, frame);
          },
      "smooth_wall_geometry requires "
      "param::smooth_wall_local_frame(boundary_intercept, frame)");
  static_assert(
      !smooth_wall_geometry ||
          requires(const Array1D<Real, 0, AMREX_SPACEDIM - 1>& point,
                   const Array1D<Real, 0, AMREX_SPACEDIM - 1>& normal,
                   const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tangent1,
                   const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tangent2) {
            param::wall_shape_operator(point, normal, tangent1, tangent2);
          },
      "smooth_wall_geometry requires param::wall_shape_operator at the "
      "same boundary intercept and local frame");
  static_assert(!smooth_wall_geometry || !uses_sharp_feature_annotation,
      "A smooth-wall geometry contract cannot be combined with sharp-feature annotations");

  static_assert(!use_ns_noslip_pressure_compatibility ||
                    (iorder_tparm == 2 && iorder_tparm_surf == 2),
      "Navier-Stokes no-slip pressure compatibility requires interp_order=2 "
      "and interp_order_surf=2 quadratic WLS");
  static_assert(!use_ns_noslip_pressure_compatibility ||
                    (eorder_tparm >= 2 && eorder_tparm_surf >= 2),
      "Navier-Stokes no-slip pressure compatibility requires extrap_order>=2 "
      "and extrap_order_surf>=2 for second-order wall pressure");
  static_assert(support_visibility_mode >= 0 && support_visibility_mode <= 2,
      "IBM support_visibility_mode must be 0 (off), 1 (audit), or 2 (filter)");
  static_assert(support_visibility_mode != 2 ||
                    (iorder_tparm == 2 && iorder_tparm_surf == 2),
      "IBM visibility filtering requires quadratic least-squares interpolation");
  static_assert(support_visibility_condition_max > Real(1.0),
      "IBM support_visibility_condition_max must exceed one");
  static_assert(!use_bi_constrained_shared_gp || iorder_tparm == 2,
      "BI-constrained shared-GP reconstruction currently requires the 3^D "
      "support storage selected by interp_order=2");
#ifdef AMREX_USE_CGAL
  static_assert(!use_bi_constrained_shared_gp,
      "BI-constrained shared-GP reconstruction currently requires the native "
      "BVH visibility backend; the CGAL support-functional path is not implemented");
#endif
  static_assert(!use_bi_constrained_shared_gp ||
                    support_visibility_mode == 2,
      "BI-constrained shared-GP reconstruction requires visibility-filtered "
      "fluid support (support_visibility_mode=2)");
  static_assert(!use_bi_constrained_shared_gp ||
                    (stationary_slip_wall || stationary_no_slip_wall),
      "BI-constrained shared-GP reconstruction currently supports stationary "
      "slip or no-slip walls only");
  static_assert(!use_bi_constrained_shared_gp ||
                    uses_zero_gradient_pressure_closure ||
                    uses_prescribed_gradient_pressure_closure ||
                    uses_euler_slip_curvature_pressure_closure ||
                    use_ns_noslip_pressure_compatibility,
      "BI-constrained shared-GP reconstruction supports zero-gradient or "
      "compatible normal pressure gradients");
  static_assert(!use_bi_constrained_shared_gp ||
                    !uses_prescribed_gradient_pressure_closure ||
                    requires(const Array1D<Real, 0, AMREX_SPACEDIM - 1>& xyz,
                             const Array1D<Real, 0, AMREX_SPACEDIM - 1>& normal) {
                      param::pressure_normal_derivative(xyz, normal);
                    },
      "BI-constrained prescribed-gradient pressure requires "
      "param::pressure_normal_derivative(xyz, normal)");
  static_assert(!use_bi_constrained_shared_gp ||
                    !uses_euler_slip_curvature_pressure_closure ||
                    stationary_slip_wall,
      "BI-constrained Euler curvature pressure compatibility requires a "
      "stationary slip wall");
  static_assert(!use_bi_constrained_shared_gp ||
                    !uses_euler_slip_curvature_pressure_closure ||
                    AMREX_SPACEDIM == 2,
      "BI-constrained Euler curvature pressure compatibility is currently "
      "restricted to two-dimensional smooth walls");
  static_assert(!use_bi_constrained_shared_gp ||
                    !uses_euler_slip_curvature_pressure_closure ||
                    expanded_bic_support,
      "BI-constrained Euler curvature pressure compatibility requires the "
      "expanded visible support candidate");
  static_assert(
      !use_bi_constrained_shared_gp ||
          !uses_euler_slip_curvature_pressure_closure ||
          requires(const Array1D<Real, 0, AMREX_SPACEDIM - 1>& xyz,
                   const Array1D<Real, 0, AMREX_SPACEDIM - 1>& normal,
                   const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tangent1,
                   const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tangent2) {
            param::wall_shape_operator(xyz, normal, tangent1, tangent2);
          },
      "BI-constrained Euler curvature pressure compatibility requires "
      "param::wall_shape_operator");
  static_assert(!use_bi_constrained_shared_gp || cls_t::NCONS == 5,
      "BI-constrained shared-GP reconstruction currently supports the "
      "single-species five-equation ideal-gas system only");
  static_assert(!entropy_wall_extension || use_bi_constrained_shared_gp,
      "Entropy-based Euler-slip extension requires BI-constrained shared-GP reconstruction");
  static_assert(!entropy_wall_extension || stationary_slip_wall,
      "Entropy-based thermodynamic extension requires a stationary slip wall");
  static_assert(!entropy_wall_extension ||
                    uses_entropy_compatible_pressure_closure,
      "Entropy-based Euler-slip extension requires a compatible prescribed "
      "or curvature-derived normal pressure gradient");
  static_assert(
      !entropy_wall_extension ||
          requires(const cls_t* closure) {
            closure->gamma;
            closure->Rspec;
          },
      "Entropy-based Euler-slip extension requires a calorically perfect "
      "ideal-gas closure exposing gamma and Rspec");
  static_assert(
      !cell_average_bic_normal_momentum ||
          (use_bi_constrained_shared_gp &&
           stationary_slip_wall && entropy_wall_extension &&
           uses_entropy_compatible_pressure_closure),
      "Cell-average-aware BI normal momentum requires the "
      "pressure-compatible Euler-slip BIC path");
  static_assert(!one_sided_entropy_jet ||
                    (use_bi_constrained_shared_gp && stationary_slip_wall &&
                     entropy_wall_extension &&
                     uses_entropy_compatible_pressure_closure),
      "The one-sided entropy jet currently requires the single-species "
      "pressure-compatible Euler-slip BIC path");
  static_assert(!rz_annular_bic_cell_average || AMREX_SPACEDIM == 2,
      "R-Z annular BI-CWLS is currently implemented only in two dimensions");
  static_assert(!rz_annular_bic_cell_average ||
                    (use_bi_constrained_shared_gp &&
                     cell_average_bic_normal_momentum &&
                     one_sided_entropy_jet && stationary_slip_wall),
      "R-Z annular BI-CWLS requires the frozen pressure-compatible, "
      "cell-average-aware shared-GP Euler-slip path");
  static_assert(!rz_bic_support_centre_recovery ||
                    rz_annular_bic_cell_average,
      "R-Z BI-CWLS support centre recovery requires annular shared-GP "
      "cell-average semantics");
  static_assert(!rz_axis_regular_entropy_jet ||
                    (rz_annular_bic_cell_average &&
                     rz_bic_support_centre_recovery &&
                     one_sided_entropy_jet),
      "R-Z axis-regular entropy extension requires the annular shared-GP "
      "one-sided entropy path");
  static_assert(bic_auxiliary_hierarchy_tolerance > Real(0.0),
      "IBM BIC auxiliary hierarchy tolerance must be positive");
  static_assert(bic_entropy_jet_variation_tolerance > Real(0.0),
      "IBM one-sided entropy-jet variation tolerance must be positive");
  static_assert(bic_normal_momentum_variation_tolerance > Real(0.0),
      "IBM cell-average normal-momentum variation tolerance must be positive");
  static_assert(!expanded_bic_support || use_bi_constrained_shared_gp,
      "Expanded BIC support requires BI-constrained shared-GP reconstruction");
  static_assert(!expanded_bic_support || AMREX_SPACEDIM == 2,
      "Expanded BIC support is currently implemented only in 2-D");
#if NUM_SPECIES > 1
  static_assert(!use_bi_constrained_shared_gp,
      "BI-constrained shared-GP reconstruction is not yet implemented for "
      "multiple species");
#endif
#ifdef AMREX_USE_CGAL
  static_assert(support_visibility_mode == 0,
      "IBM support visibility currently requires the CPU/GPU BVH backend");
#endif

  // number of ghost layers needed for IB method
  static constexpr int  ghost_layers = param::ghost_layers;
  static_assert(ghost_layers >= 1 && ghost_layers <= cls_t::NGHOST,
      "IBM ghost_layers must lie in [1, cls_t::NGHOST]");
  
  // true: interior of closed geometry is solid; false: interior is fluid
  static constexpr bool interior_is_solid = param::interior_is_solid; 

  // ideal number of interpolation points for each image point(ghost point extrapolation and surface reconstruction)
  static constexpr int  N_InterP      = ipow(iorder_tparm + 1, AMREX_SPACEDIM);
  static constexpr int  N_InterP_surf = ipow(iorder_tparm_surf + 1, AMREX_SPACEDIM);

  using SURFIMP = surfImp_t<eorder_tparm_surf, iorder_tparm_surf>;
  using GPSTORE = GPStore<eorder_tparm, iorder_tparm, cls_t::NCONS>;
  using GPSTOREVIEW =
      GPStoreView<eorder_tparm, iorder_tparm, cls_t::NCONS>;

  // MultiFabs pointer to Amr class instance
  Amr* amr_p;
  // Per-level marker MultiFabs. Production stores solid + reconstructed-GP
  // flags. Validation adds the actual-WENO-read flag and topological GP layer.
  static constexpr int marker_ncomp = 2;
  Vector<IBMultiFab<uint8_t>*> bmf_a;
  // Level-wide flattened ghost-point storage (CSR indexed by local FAB)
  Vector<GPSTORE> gpstore_a;
  // parameters for cell size and refinement ratio
  Vector<IntVect> rratio_a;                                 // vector of refinement ratio per level in each direction
  Vector<GpuArray<Real, AMREX_SPACEDIM>> dx_a;              // vector of cell sizes per level in each direction
  Vector<Real> diag_a;                                      // vector of cell diagonal length per level
  Vector<Real> di_a;                                        // image point distance per level for ghost point extrapolation
  Vector<Real> di_a_surf;                                   // image point distance per level for surface reconstruction
  // Image-point interpolation reaches farther than the numerical RHS stencil.
  // Keep these halos separate from cls_t::NGHOST so eorder/alpha changes do
  // not make GP or surface reconstruction depend on FAB/MPI partitioning.
  Vector<int> volume_interp_nghost_a;
  Vector<int> surface_interp_nghost_a;

  // Image-point halos and coarse-fine support coverage.
  #include "ibm_solver_amr_support.h"
  // geometry related data
  int ngeom = 0;                                            // number of geometries
  Vector<GeomType> geom_a;                                  // IB explicit geometry
#ifdef AMREX_USE_CGAL
  Vector<Tree> tree_a;                                      // CGAL AABB tree per geometry
  Vector<PrimitiveIndexMap> idxmap_a;                       // CGAL primitive-to-index map per geometry
#else
  Vector<BVH> bvh_a;                                        // BVH per geometry (replaces CGAL AABB tree)
#endif
  Vector<std::unique_ptr<inside_t>> inout_fa;                // in out testing function per geometry
  Vector<Bbox> bbox_a;                                      // bounding box per geometry (world-frame, for fast rejection)
  Vector<Bbox> bbox_body_a;                                 // bounding box per geometry (body-frame, constant after init)
  Vector<RigidTransform> transform_a;                       // rigid-body transform per geometry (body→world)
  Vector<std::string> geom_names;                           // geometry names stripped from input filenames

  Gpu::ManagedVector<LocalFrame> LocalFrame_a;              // local orthonormal frame matrix (flattened, body-frame)
  Gpu::ManagedVector<SurfElem> SurfElem_a;                  // surface element area and coordinates (flattened, body-frame)
  Gpu::ManagedVector<int> geom_offsets;                     // Start index for each geometry in flattened arrays
 
  // surface related data
  int ntotalfaces = 0;                                      // number of faces/edges across all geometries
  SURFIMP surfimp_soa;                                      // face/edge image point information (SoA structure)
  surfPhys_t surfphys_soa;                                  // face/edge physical data and identification (SoA structure)
  Vector<FaceCSR> faces_per_level;                          // faces integers per fab and per level [lev] (CSR format)
  /** 
   * \brief Destructor to release allocated memory
   */
  ~ibm_solver_t() noexcept
  {
    // Release per-level IBMultiFab pointers if any remain
    for (auto*& p : bmf_a) {
      if (p) { delete p; p = nullptr; }
    }

    // Release inside/outside testers (unique_ptr handles cleanup automatically)
    inout_fa.clear();
  }

  /**
   * \brief Explicitly release all GPU-managed memory before amrex::Finalize().
   *
   * Must be called while AMReX arenas are still alive. The global inline
   * IBM::ib outlives amrex::Finalize(), so its implicit destructor would
   * free arena memory after the arena is destroyed (static destruction
   * order problem).  Calling cleanup() first leaves the destructor with
   * nothing to free.
   */
  void cleanup() noexcept
  {
    // Delete raw-pointer members first (their internals use arena memory)
    for (auto*& p : bmf_a)    { if (p) { delete p; p = nullptr; } }
    inout_fa.clear();  // unique_ptr handles cleanup

    // Swap-with-empty idiom: guarantees capacity→0 and arena memory freed NOW,
    // so the post-Finalize destructor finds nothing to deallocate.
    { decltype(bmf_a)        tmp; tmp.swap(bmf_a);        }
    { decltype(inout_fa)     tmp; tmp.swap(inout_fa);     }
    { decltype(LocalFrame_a) tmp; tmp.swap(LocalFrame_a); }
    { decltype(SurfElem_a)   tmp; tmp.swap(SurfElem_a);   }
    { decltype(geom_offsets) tmp; tmp.swap(geom_offsets);  }

    // Compound types: destroying elements calls their ManagedVector destructors
    { decltype(geom_a) tmp; tmp.swap(geom_a); }
#ifdef AMREX_USE_CGAL
    { decltype(tree_a)   tmp; tmp.swap(tree_a);   }
    { decltype(idxmap_a) tmp; tmp.swap(idxmap_a); }
#else
    { decltype(bvh_a) tmp; tmp.swap(bvh_a); }
#endif
    { decltype(bbox_a)      tmp; tmp.swap(bbox_a);      }
    { decltype(bbox_body_a) tmp; tmp.swap(bbox_body_a); }
    { decltype(transform_a) tmp; tmp.swap(transform_a); }
    { decltype(geom_names)  tmp; tmp.swap(geom_names);  }
    { decltype(gpstore_a)  tmp; tmp.swap(gpstore_a);  }
    // These structs have their own clear() with shrink_to_fit()
    surfimp_soa.clear();
    surfphys_soa.clear();
    { decltype(faces_per_level) tmp; tmp.swap(faces_per_level); }

    amr_p = nullptr;
  }

  /**
   * Opt-in wall-clock diagnostics for geometry/IBM setup.  These are kept
   * separate from TinyProfiler so a long geometry phase is reported as soon
   * as it finishes (and before a short capacity job reaches Finalize).
   *
   * Runtime controls:
   *   ib.timing = 1
   *   ib.timing_detail = 1        (CPU-only GP sub-phase sampling)
   *   ib.timing_warn_seconds = 60
   */
  static bool performanceTimingEnabled()
  {
    static const bool enabled = [] {
      int value = 0;
      ParmParse pp("ib");
      pp.query("timing", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1, "ib.timing must be 0 or 1");
      return value != 0;
    }();
    return enabled;
  }

  /**
   * Opt-in finite-volume support semantics for the no-slip wall surface
   * observer. This changes only BI/IP surface recovery used for Cp, shear and
   * heat-flux output; it never changes a shared GP state or the full-cell RHS.
   */
  static bool nsSurfaceCellAverageRecoveryEnabled()
  {
    static const bool enabled = [] {
      int value = 0;
      ParmParse pp("ib");
      pp.query("ns_surface_cell_average_recovery", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1,
          "ib.ns_surface_cell_average_recovery must be 0 or 1");
      return value != 0;
    }();
    return enabled;
  }

  /**
   * Opt-in finite-volume support and target semantics for the evolved no-slip
   * shared-GP state. The option changes only the auxiliary solid-side state
   * consumed by the existing Cartesian flux/viscous operators.
   */
  static bool nsGPCellAverageRecoveryEnabled()
  {
    static const bool enabled = [] {
      int value = 0;
      ParmParse pp("ib");
      pp.query("ns_gp_cell_average_recovery", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1,
          "ib.ns_gp_cell_average_recovery must be 0 or 1");
      return value != 0;
    }();
    return enabled;
  }

  void validateAndPrintNSGPRecoveryConfiguration() const
  {
    const bool enabled = nsGPCellAverageRecoveryEnabled();
    amrex::Print()
        << "[IBM-NS-GP] cell_average_point_recovery=" << int(enabled)
        << " target=conservative_cell_average"
        << " extension=unique_shared_gp"
        << " rhs=full_cartesian\n";
    if (!enabled) return;

    if constexpr (!ns_cartesian_bic_cell_average) {
      amrex::Abort(
          "ib.ns_gp_cell_average_recovery requires a build with "
          "ProbIB::ns_cartesian_bic_cell_average=true");
    }
    if constexpr (AMREX_SPACEDIM != 2) {
      amrex::Abort(
          "ib.ns_gp_cell_average_recovery is currently qualified only for "
          "two-dimensional Cartesian no-slip flow");
    }
    if constexpr (!use_bi_constrained_shared_gp ||
                  !stationary_no_slip_wall) {
      amrex::Abort(
          "ib.ns_gp_cell_average_recovery requires stationary no-slip "
          "BI-CWLS shared-GP reconstruction");
    }
    if constexpr (cls_t::NCONS != 5) {
      amrex::Abort(
          "ib.ns_gp_cell_average_recovery is not qualified for multispecies "
          "conservative states");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !amr_p->Geom(0).IsRZ(),
        "ib.ns_gp_cell_average_recovery requires Cartesian geometry");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amr_p->maxLevel() == 0,
        "ib.ns_gp_cell_average_recovery is currently qualified only on a "
        "uniform level-0 grid");
  }

  void validateAndPrintNSSurfaceRecoveryConfiguration() const
  {
    const bool enabled = nsSurfaceCellAverageRecoveryEnabled();
    amrex::Print()
        << "[IBM-NS-Surface] cell_average_point_recovery="
        << int(enabled)
        << " evolved_state=unchanged rhs=full_cartesian\n";
    if (!enabled) return;

    if constexpr (AMREX_SPACEDIM != 2) {
      amrex::Abort(
          "ib.ns_surface_cell_average_recovery is currently qualified only "
          "for two-dimensional Cartesian no-slip verification");
    }
    if constexpr (!supports_complete_viscous_surface_traction) {
      amrex::Abort(
          "ib.ns_surface_cell_average_recovery requires a stationary "
          "no-slip wall model");
    }
    if constexpr (iorder_tparm_surf != 2 || eorder_tparm_surf < 2) {
      amrex::Abort(
          "ib.ns_surface_cell_average_recovery requires quadratic surface "
          "interpolation and at least two image points");
    }
    if constexpr (cls_t::NCONS != 5) {
      amrex::Abort(
          "ib.ns_surface_cell_average_recovery is not qualified for "
          "multispecies conservative states");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !amr_p->Geom(0).IsRZ(),
        "ib.ns_surface_cell_average_recovery currently requires Cartesian "
        "geometry");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amr_p->maxLevel() == 0,
        "ib.ns_surface_cell_average_recovery is currently qualified only on "
        "a uniform level-0 grid");
  }

  static bool detailedPerformanceTimingEnabled()
  {
    static const bool enabled = [] {
      int value = 0;
      ParmParse pp("ib");
      pp.query("timing_detail", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1, "ib.timing_detail must be 0 or 1");
      return value != 0;
    }();
    return performanceTimingEnabled() && enabled;
  }

  static Real performanceTimingWarningSeconds()
  {
    static const Real threshold = [] {
      Real value = Real(60.0);
      ParmParse pp("ib");
      pp.query("timing_warn_seconds", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value > Real(0.0), "ib.timing_warn_seconds must be positive");
      return value;
    }();
    return threshold;
  }

  static Real beginPerformanceTiming()
  {
    if (!performanceTimingEnabled()) return Real(0.0);
    Gpu::streamSynchronize();
    return amrex::second();
  }

  static void reportPerformanceDuration(const char* phase, int lev,
                                        Real local_seconds,
                                        Long local_work_items = Long(-1),
                                        bool replicated_work = false)
  {
    if (!performanceTimingEnabled()) return;
    Gpu::streamSynchronize();

    Real minimum_seconds = local_seconds;
    Real maximum_seconds = local_seconds;
    Real sum_seconds = local_seconds;
    ParallelDescriptor::ReduceRealMin(minimum_seconds);
    ParallelDescriptor::ReduceRealMax(maximum_seconds);
    ParallelDescriptor::ReduceRealSum(sum_seconds);
    const Real average_seconds =
        sum_seconds / Real(ParallelDescriptor::NProcs());

    Long total_work_items = local_work_items;
    Long maximum_local_work_items = local_work_items;
    if (local_work_items >= Long(0)) {
      ParallelDescriptor::ReduceLongSum(total_work_items);
      ParallelDescriptor::ReduceLongMax(maximum_local_work_items);
    }

    std::ostringstream line;
    line << std::setprecision(8)
         << "[IBM-Timing] phase=" << phase
         << " level=" << lev
         << " min_s=" << minimum_seconds
         << " avg_s=" << average_seconds
         << " max_s=" << maximum_seconds
         << " imbalance="
         << (average_seconds > Real(0.0)
                 ? maximum_seconds / average_seconds
                 : Real(0.0));
    if (local_work_items >= Long(0)) {
      line << " work_global=" << total_work_items
           << " work_local_max=" << maximum_local_work_items
           << " work_scope=" << (replicated_work ? "replicated" : "distributed");
      if (replicated_work && maximum_local_work_items > Long(0)) {
        line << " maxrank_us_per_item="
             << maximum_seconds * Real(1.0e6) /
                    Real(maximum_local_work_items);
      } else if (!replicated_work && total_work_items > Long(0)) {
        // Sum of per-rank phase times divided by distributed work is a useful
        // decomposition-independent CPU cost per item.  max_s/total_work
        // would instead improve artificially as ranks are added.
        line << " cpu_us_per_item="
             << sum_seconds * Real(1.0e6) /
                    Real(total_work_items);
      }
    }
    line << "\n";
    amrex::Print() << line.str();

    if (maximum_seconds > performanceTimingWarningSeconds()) {
      amrex::Print()
          << "[IBM-Timing-WARN] phase=" << phase
          << " level=" << lev
          << " max_s=" << maximum_seconds
          << " threshold_s=" << performanceTimingWarningSeconds()
          << " -- IBM geometry processing is consuming excessive wall time\n";
    }
  }

  static void finishPerformanceTiming(const char* phase, int lev,
                                      Real start_seconds,
                                      Long local_work_items = Long(-1),
                                      bool replicated_work = false)
  {
    if (!performanceTimingEnabled()) return;
    Gpu::streamSynchronize();
    reportPerformanceDuration(phase, lev,
                              amrex::second() - start_seconds,
                              local_work_items, replicated_work);
  }

  /**
   * \brief Initializes the Immersed Boundary (IB) method structures and geometry.
   *
   * This function sets up the AMR pointer, resizes internal data structures based on the maximum AMR level,
   * computes grid metrics (cell sizes, diagonals) for all levels, and loads the IB geometry from files.
   *
   * \param pointer_amr Pointer to the main Amr class instance.
   */
  void init(Amr* pointer_amr)
  {
    amr_p = pointer_amr;
    IBM::method_config::validate_and_print_runtime_contract();
    validateAndPrintInviscidConfiguration();
    validateAndPrintNSGPRecoveryConfiguration();
    validateAndPrintNSSurfaceRecoveryConfiguration();
    if constexpr (rz_annular_bic_cell_average) {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          amr_p->Geom(0).IsRZ(),
          "rz_annular_bic_cell_average requires R-Z geometry");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          amr_p->maxLevel() == 0,
          "RZ-2 annular shared-GP semantics are currently certified only on "
          "a uniform level-0 grid");
    }
    rratio_a = amr_p->refRatio();
    int lmax = amr_p->maxLevel();

    bmf_a.resize(lmax + 1);
    gpstore_a.resize(lmax + 1);
    faces_per_level.resize(lmax + 1);

    dx_a.resize(lmax + 1);
    dx_a[0] = amr_p->Geom(0).CellSizeArray();
    for (int i = 1; i <= lmax; i++) {
      for (int j = 0; j < AMREX_SPACEDIM; j++) {
        dx_a[i][j] = dx_a[i - 1][j] / rratio_a[i - 1][j];
      }
    }

    di_a.resize(lmax + 1);
    di_a_surf.resize(lmax + 1);
    diag_a.resize(lmax + 1);
    volume_interp_nghost_a.resize(lmax + 1);
    surface_interp_nghost_a.resize(lmax + 1);

    for (int i = 0; i <= lmax; i++) {
      diag_a[i] = std::sqrt(
        AMREX_D_TERM( std::pow(dx_a[i][0], 2),
                      + std::pow(dx_a[i][1], 2),
                      + std::pow(dx_a[i][2], 2)) );

      di_a[i] = cim * diag_a[i];
      di_a_surf[i] = cim_surf * diag_a[i];

      Real min_dx = dx_a[i][0];
      for (int d = 1; d < AMREX_SPACEDIM; ++d) {
        min_dx = amrex::min(min_dx, dx_a[i][d]);
      }
      const Real diagonal_cells = diag_a[i] / min_dx;
      // A centred interpolation block can extend one index beyond the cell
      // containing its image point.  The entropy jet may translate the full
      // block by two additional cells to obtain six visible 2-D supports.
      // Include that cold geometry-search reach in the stage primitive halo.
      const int interpolation_padding =
          1 + (sharp_feature_tangential_reconstruction ? 1 : 0) +
          (one_sided_entropy_jet ? 2 : 0);

      // The first volume IP may move as far as 2*di and each later IP adds
      // one di.  Include the GP-to-wall distance (<= one cell diagonal) and
      // the interpolation block's outer support.
      const Real volume_reach =
          (Real(1.0) + Real(eorder_tparm + 1) * cim) *
          diagonal_cells;
      volume_interp_nghost_a[i] = amrex::max(
          cls_t::NGHOST,
          static_cast<int>(std::ceil(volume_reach)) + interpolation_padding);

      // A surface point may be half a cell diagonal from the centre of its
      // owning cell. Its first IP search may move as far as 3*di; subsequent
      // points again add one di.
      const Real surface_reach =
          (Real(0.5) + Real(eorder_tparm_surf + 2) * cim_surf) *
          diagonal_cells;
      surface_interp_nghost_a[i] = amrex::max(
          cls_t::NGHOST,
          static_cast<int>(std::ceil(surface_reach)) + interpolation_padding);
    }

    // Read and build geometry independently on every MPI rank.  Report the
    // slowest rank because the first collective after initialization cannot
    // progress faster than that straggler.
    const Real t_geom = beginPerformanceTiming();
    read_geom();
    finishPerformanceTiming("read_geom_total", -1, t_geom,
                            Long(ntotalfaces), true);

    // One-time notice: at extrap_order>=2 a custom wall model that only
    // implements a legacy compute_surfIB signature never receives disIM /
    // n_valid, so its Neumann quantities (zero-grad P/T, slip u_t) silently
    // stay 1st-order while ghost extrapolation runs at 2nd. Make that loud.
    if constexpr (eorder_tparm >= 2) {
      using eo_tag  = std::integral_constant<int, eorder_tparm>;
      using Vec1D   = Array1D<Real, 0, AMREX_SPACEDIM - 1>;
      using DisArr  = Array1D<Real, 0, eorder_tparm - 1>;
      using Prims2D = Array2D<Real, 0, eorder_tparm + 1, 0, cls_t::NPRIM - 1>;
      constexpr bool wm_compat = ibm_detail::is_detected_v<
          ibm_detail::compute_surfIB_expr, wallmodel,
          eo_tag, const Vec1D&, const Vec1D&, const Vec1D&, const Vec1D&,
          const DisArr&, int, Prims2D&, const int, const cls_t*,
          const ibm_pressure_compatibility_t&>;
      constexpr bool wm_disim = ibm_detail::is_detected_v<
          ibm_detail::compute_surfIB_expr, wallmodel,
          eo_tag, const Vec1D&, const Vec1D&, const Vec1D&, const Vec1D&,
          const DisArr&, int, Prims2D&, const int, const cls_t*>;
      constexpr bool wm_2nd = wm_compat || wm_disim;
      if constexpr (!wm_2nd) {
        amrex::Print() << "  IBM note: extrap_order=" << eorder_tparm
                       << " but the wall model has only a legacy compute_surfIB signature;\n"
                       << "            wall Neumann BC stays 1st-order (see ibm_walltypes.h for the\n"
                       << "            disIM-aware template to adopt).\n";
      }
    }
  }

  /** \brief Validate and print the frozen shared-GP/full-Cartesian method. */
  static void validateAndPrintInviscidConfiguration()
  {
    ParmParse pp("cns");
    constexpr std::array<const char*, 11> retired_experimental_keys{{
        "llf_ibm_face_local_crossing",
        "llf_ibm_fluid_shell",
        "llf_ibm_shifted_shell",
        "llf_ibm_wall_matched_crossing",
        "llf_ibm_conservative_crossing",
        "afd_ibm_face_local_crossing",
        "afd_ibm_crossing_wls",
        "afd_ibm_fluid_shell",
        "afd_ibm_shifted_shell",
        "afd_ibm_conservative_crossing",
        "afd_ibm_patch_conservative_crossing"}};
    for (const char* key : retired_experimental_keys) {
      int enabled = 0;
      pp.query(key, enabled);
      if (enabled != 0) {
        amrex::Abort(
            std::string("cns.") + key +
            " belongs to a retired face-local/direct-wall IBM experiment; "
            "the production method is pure shared-GP/full-Cartesian");
      }
    }

    int rz_annular_bic = 0;
    int rz_support_centre_recovery = 0;
    int rz_axis_regular_entropy = 0;
    pp.query("rz_gp_annular_bic", rz_annular_bic);
    pp.query("rz_gp_support_centre_recovery",
             rz_support_centre_recovery);
    pp.query("rz_gp_axis_regular_entropy_jet",
             rz_axis_regular_entropy);
    const std::array<int, 3> rz_options{{
        rz_annular_bic, rz_support_centre_recovery,
        rz_axis_regular_entropy}};
    for (int value : rz_options) {
      if (value != 0 && value != 1) {
        amrex::Abort("R-Z shared-GP options must be 0 or 1");
      }
    }
    if (rz_annular_bic != int(rz_annular_bic_cell_average)) {
      amrex::Abort(
          "cns.rz_gp_annular_bic must agree with the compile-time "
          "ProbIB::rz_annular_bic_cell_average option");
    }
    if (rz_support_centre_recovery !=
        int(rz_bic_support_centre_recovery)) {
      amrex::Abort(
          "cns.rz_gp_support_centre_recovery must agree with the "
          "compile-time ProbIB::rz_bic_support_centre_recovery option");
    }
    if (rz_axis_regular_entropy !=
        int(rz_axis_regular_entropy_jet)) {
      amrex::Abort(
          "cns.rz_gp_axis_regular_entropy_jet must agree with the "
          "compile-time ProbIB::rz_axis_regular_entropy_jet option");
    }

    amrex::Print()
        << "[IBM-Config] method=pure_shared_gp_full_cartesian"
        << " shared_gp_reconstruction="
        << int(shared_gp_reconstruction)
        << " ghost_layers=" << ghost_layers
        << " interp/extrap=" << iorder_tparm << '/' << eorder_tparm
        << " alpha=" << cim
        << " visibility_mode=" << support_visibility_mode
        << " expanded_bic_support="
        << int(expanded_bic_support)
        << " cell_average_bic_normal_momentum="
        << int(cell_average_bic_normal_momentum)
        << " one_sided_entropy_jet=" << int(one_sided_entropy_jet)
        << " rz_annular_bic_cell_average="
        << int(rz_annular_bic_cell_average)
        << " rz_bic_support_centre_recovery="
        << int(rz_bic_support_centre_recovery)
        << " rz_axis_regular_entropy_jet="
        << int(rz_axis_regular_entropy_jet)
        << " auxiliary_hierarchy_tolerance="
        << bic_auxiliary_hierarchy_tolerance
        << " entropy_jet_variation_tolerance="
        << bic_entropy_jet_variation_tolerance
        << " normal_momentum_variation_tolerance="
        << bic_normal_momentum_variation_tolerance
        << " smooth_wall_geometry=" << int(smooth_wall_geometry)
        << " entropy_wall_extension=" << int(entropy_wall_extension)
        << " pressure_closure="
        << int(wall_pressure_closure)
        << '\n';
#ifdef CERISSE_GIT_VERSION
    amrex::Print() << "[IBM-Config] cerisse_git="
                   << CERISSE_GIT_VERSION;
#else
    amrex::Print() << "[IBM-Config] cerisse_git=not-embedded";
#endif
#ifdef AMREX_GIT_VERSION
    amrex::Print() << " amrex_git=" << AMREX_GIT_VERSION;
#endif
#ifdef IBM_DIAGNOSTIC_BUILD_ID
    amrex::Print() << " build_id=" << IBM_DIAGNOSTIC_BUILD_ID;
#endif
#ifdef IBM_DIAGNOSTIC_SOURCE_HASH
    amrex::Print() << " source_hash=" << IBM_DIAGNOSTIC_SOURCE_HASH;
#endif
#if defined(__CUDACC_VER_MAJOR__)
    amrex::Print() << " cuda_compiler=" << __CUDACC_VER_MAJOR__ << '.'
                   << __CUDACC_VER_MINOR__;
#else
    amrex::Print() << " cuda_compiler=none";
#endif
    amrex::Print() << '\n';
#ifdef AMREX_USE_CUDA
    int cuda_runtime_version = 0;
    int cuda_driver_version = 0;
    int cuda_device = -1;
    cudaDeviceProp cuda_properties{};
    const cudaError_t runtime_status =
        cudaRuntimeGetVersion(&cuda_runtime_version);
    const cudaError_t driver_status =
        cudaDriverGetVersion(&cuda_driver_version);
    const cudaError_t device_status = cudaGetDevice(&cuda_device);
    const cudaError_t properties_status =
        device_status == cudaSuccess
            ? cudaGetDeviceProperties(&cuda_properties, cuda_device)
            : device_status;
    amrex::Print()
        << "[IBM-Config] cuda_runtime="
        << cuda_runtime_version / 1000 << '.'
        << (cuda_runtime_version % 1000) / 10
        << " cuda_driver=" << cuda_driver_version / 1000 << '.'
        << (cuda_driver_version % 1000) / 10
        << " device="
        << (properties_status == cudaSuccess ? cuda_properties.name
                                             : "unavailable")
        << " compute_capability="
        << (properties_status == cudaSuccess ? cuda_properties.major : -1)
        << '.'
        << (properties_status == cudaSuccess ? cuda_properties.minor : -1)
        << " status(runtime,driver,device,properties)="
        << int(runtime_status) << ',' << int(driver_status) << ','
        << int(device_status) << ',' << int(properties_status) << '\n';
#endif
  }

  /**
   * \brief create IBMultiFab at a level and store pointers to it
   * \param bxa BoxArray for the level
   * \param dm DistributionMapping for the level
   * \param lev The AMR level
   */
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    // Prevent memory leak if MF already exists
    if (lev < 0 || lev >= static_cast<int>(bmf_a.size())) {
        amrex::Abort("ibm_solver_t::build_mf: lev is out of bounds");
    }

    // Prevent memory leak if MF already exists
    if (bmf_a[lev] != nullptr) {
        destroy_mf(lev);
    }

    // Use default MFInfo (configured to use The_Managed_Arena in IBMultiFab.h)
    bmf_a[lev] = new IBMultiFab<uint8_t>(
        bxa, dm, marker_ncomp, interpolationMarkerNghost(lev));
  }

  /**
   * \brief destroy IBMultiFab at a level
   * \param lev The AMR level of the IBMultiFab to be destroyed.
   */
  void destroy_mf(int lev)
  {
    if (lev >= 0 && lev < static_cast<int>(bmf_a.size())) {
      // safe delete: delete nullptr is valid in C++
      delete bmf_a[lev];
      bmf_a[lev] = nullptr;  // Prevent dangling pointer
    }
    if (lev >= 0 && lev < static_cast<int>(gpstore_a.size())) {
      gpstore_a[lev].clear();
    }
  }

  /**
   * \brief Rebuild all geometry-level data from current vertex positions.
   *
   * After externally modifying vertex positions in geom_a (e.g. for moving
   * geometry), call this to reconstruct BVH trees, inside/outside testers,
   * bounding boxes, and the surface cache (LocalFrame_a, SurfElem_a).
   * Must be called BEFORE rebuildIBM().
   */
  /**
   * \brief Lightweight rigid-body transform update (replaces full BVH rebuild).
   *
   * For rigid-body FSI, the geometry shape never changes — only its position
   * and orientation.  Instead of moving all vertices and rebuilding the BVH
   * tree, InsideTester, LocalFrame, and SurfElem every time step, we store
   * the geometry in its reference (body) frame and only update the transform.
   *
   * All query functions (computeMarkers, initialiseGPs, computeAllGPs, etc.)
   * apply the inverse transform to query points before using the body-frame
   * BVH, then forward-transform the results back to the world frame.
   * The results are mathematically identical to a full rebuild.
   *
   * \param geomIdx  Index of the geometry to update.
   * \param T        New rigid-body transform (body → world).
   */
  void updateRigidTransform(int geomIdx, const RigidTransform& T)
  {
      transform_a[geomIdx] = T;
      // Update world-frame bounding box from body-frame bbox + new transform
      bbox_a[geomIdx] = T.transform_bbox(bbox_body_a[geomIdx]);
  }

  /**
   * \brief Annotate 2-D edges with arc distance to the nearest sharp vertex.
   *
   * Refined polygon vertices on a straight edge have zero turning angle and
   * are ignored.  True corners are multi-source seeds on the cyclic edge
   * graph; Dijkstra propagation gives a topology-aware distance at every
   * endpoint.  Runtime GP/surface kernels then need only an O(1) edge-local
   * interpolation to construct the compact pressure fallback.
   */
  void annotateSharpFeatures2D()
  {
#if (AMREX_SPACEDIM == 2)
    if constexpr (!uses_sharp_feature_annotation) return;

    constexpr Real pi =
        Real(3.141592653589793238462643383279502884L);
    using QueueEntry = std::pair<Real, int>;

    for (int g = 0; g < ngeom; ++g) {
      const int n = static_cast<int>(geom_a[g].size());
      if (n < 3) continue;
      const int offset = geom_offsets[g];
      AMREX_ALWAYS_ASSERT(offset >= 0);
      AMREX_ALWAYS_ASSERT(offset + n <= static_cast<int>(SurfElem_a.size()));

      std::vector<Real> x(n), y(n);
      for (int i = 0; i < n; ++i) {
        const auto& p = geom_a[g].vertex(i);
#ifdef AMREX_USE_CGAL
        x[i] = p.x();
        y[i] = p.y();
#else
        x[i] = p[0];
        y[i] = p[1];
#endif
      }

      Real xmin = x[0];
      Real xmax = x[0];
      Real ymin = y[0];
      Real ymax = y[0];
      for (int i = 1; i < n; ++i) {
        xmin = std::min(xmin, x[i]);
        xmax = std::max(xmax, x[i]);
        ymin = std::min(ymin, y[i]);
        ymax = std::max(ymax, y[i]);
      }
      const Real geom_scale = std::max(xmax - xmin, ymax - ymin);
      const Real edge_tol = std::max(
          Real(64.0) * std::numeric_limits<Real>::epsilon() * geom_scale,
          std::numeric_limits<Real>::min());

      std::vector<Real> distance(
          n, std::numeric_limits<Real>::max());
      std::priority_queue<QueueEntry, std::vector<QueueEntry>,
                          std::greater<QueueEntry>> queue;
      int sharp_vertices = 0;

      for (int i = 0; i < n; ++i) {
        const int im = (i + n - 1) % n;
        const int ip = (i + 1) % n;
        const Real ax = x[i] - x[im];
        const Real ay = y[i] - y[im];
        const Real bx = x[ip] - x[i];
        const Real by = y[ip] - y[i];
        const Real la = std::sqrt(ax * ax + ay * ay);
        const Real lb = std::sqrt(bx * bx + by * by);
        if (la <= edge_tol || lb <= edge_tol) continue;
        const Real dot = (ax * bx + ay * by) / (la * lb);
        const Real cross = (ax * by - ay * bx) / (la * lb);
        const Real turn_deg = std::atan2(std::abs(cross),
                                         std::max(Real(-1.0),
                                                  std::min(Real(1.0), dot))) *
                              Real(180.0) / pi;
        if (turn_deg >= sharp_feature_angle_deg) {
          distance[i] = Real(0.0);
          queue.emplace(Real(0.0), i);
          ++sharp_vertices;
        }
      }

      while (!queue.empty()) {
        const auto [du, u] = queue.top();
        queue.pop();
        if (du > distance[u]) continue;

        const int up = (u + 1) % n;
        const int um = (u + n - 1) % n;
        const Real wp = SurfElem_a[offset + u].measure;
        const Real wm = SurfElem_a[offset + um].measure;
        if (du + wp < distance[up]) {
          distance[up] = du + wp;
          queue.emplace(distance[up], up);
        }
        if (du + wm < distance[um]) {
          distance[um] = du + wm;
          queue.emplace(distance[um], um);
        }
      }

      if (sharp_vertices > 0) {
        for (int i = 0; i < n; ++i) {
          SurfElem_a[offset + i].sharp_distance_source = distance[i];
          SurfElem_a[offset + i].sharp_distance_target =
              distance[(i + 1) % n];
        }
      }

      amrex::Print() << "[IBM] 2-D sharp-feature annotation: geometry=" << g
                     << " vertices=" << sharp_vertices
                     << " threshold_deg=" << sharp_feature_angle_deg
                     << " radius_cells=" << sharp_feature_radius_cells
                     << " tangential_wls="
                     << int(sharp_feature_tangential_reconstruction)
                     << " tangential_radius_cells="
                     << sharp_feature_tangential_radius_cells << "\n";
    }
#endif
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static Real sharpFeaturePressureFallback(
      const SurfElem& edge, Real coordinate_from_source, Real h)
  {
#if (AMREX_SPACEDIM == 2)
    if constexpr (sharp_feature_pressure_limiter) {
      if (edge.sharp_distance_source < Real(0.0) ||
          edge.sharp_distance_target < Real(0.0) ||
          !(h > Real(0.0))) {
        return Real(0.0);
      }
      const Real s = amrex::max(
          Real(0.0), amrex::min(edge.measure, coordinate_from_source));
      const Real distance = amrex::min(
          edge.sharp_distance_source + s,
          edge.sharp_distance_target + edge.measure - s);
      const Real radius = sharp_feature_radius_cells * h;
      const Real eta = Real(1.0) - amrex::max(
          Real(0.0), amrex::min(Real(1.0), distance / radius));
      return eta * eta * (Real(3.0) - Real(2.0) * eta);
    }
#endif
    amrex::ignore_unused(edge, coordinate_from_source, h);
    return Real(0.0);
  }

#ifndef AMREX_USE_CGAL
  struct WorldSegmentHit {
    int geometry = -1;
    int global_element = -1;
    Real fraction = std::numeric_limits<Real>::max();
    Point point{};

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool hit() const noexcept { return global_element >= 0; }
  };

  /// Find the first crossing of a world-frame open segment against every IBM
  /// body.  Global element id breaks equal-fraction ties deterministically, so
  /// overlapping FAB copies of one Cartesian face select the same owner.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static WorldSegmentHit firstSurfaceHitWorld(
      const Point& start, const Point& end,
      const GpuArray<BVH4QueryView, MAX_NGEOM>& queries,
      const GpuArray<RigidTransform, MAX_NGEOM>& transforms,
      const GpuArray<int, MAX_NGEOM + 1>& geometry_offsets,
      int geometry_count)
  {
    WorldSegmentHit result;
    for (int geometry = 0; geometry < geometry_count; ++geometry) {
      const Point start_body = transforms[geometry].to_body(start);
      const Point end_body = transforms[geometry].to_body(end);
      const SegmentHitResult candidate =
          queries[geometry].first_segment_hit(start_body, end_body);
      if (!candidate.hit()) continue;

      const int global_element =
          geometry_offsets[geometry] + candidate.prim_id;
      const Real tie_tolerance = Real(64.0) *
          std::numeric_limits<Real>::epsilon();
      const bool first = !result.hit();
      const bool earlier =
          candidate.fraction < result.fraction - tie_tolerance;
      const bool tied =
          amrex::Math::abs(candidate.fraction - result.fraction) <=
              tie_tolerance &&
          global_element < result.global_element;
      if (!(first || earlier || tied)) continue;

      result.geometry = geometry;
      result.global_element = global_element;
      result.fraction = candidate.fraction;
    }
    if (result.hit()) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        result.point[d] =
            start[d] + result.fraction * (end[d] - start[d]);
      }
    }
    return result;
  }

  /**
   * Build the line-of-sight mask for all interpolation supports of one GP or
   * surface element.  Both endpoints are in cells classified as fluid, but a
   * finite-segment surface hit proves only that the straight support segment is
   * obstructed.  It is not, by itself, a fluid connected-component test.
   */
  template <int eorder_t, int iorder_t,
            int NIP = ipow(iorder_t + 1, AMREX_SPACEDIM)>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void buildSupportVisibilityMask(
      Array2D<uint8_t, 0, eorder_t - 1, 0, NIP - 1>& support_visible,
      Array2D<int, 0, eorder_t - 1, 0, NIP - 1>&
          first_hit_element,
      Array1D<int, 0, eorder_t - 1>& candidate_count,
      Array1D<int, 0, eorder_t - 1>& visible_count,
      const Array2D<Real, 0, eorder_t - 1, 0,
                    AMREX_SPACEDIM - 1>& imp_xyz,
      const Array2D<int, 0, eorder_t - 1, 0,
                    AMREX_SPACEDIM - 1>& imp_ijk,
      const Array1D<int, 0, eorder_t - 1>& imp_ninterp,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
      const Array4<uint8_t const>& ibMarkers,
      const GpuArray<BVH4QueryView, MAX_NGEOM>& queries,
      const GpuArray<RigidTransform, MAX_NGEOM>& transforms,
      const GpuArray<int, MAX_NGEOM + 1>& geometry_offsets,
      int geometry_count, int& fluid_supports, int& occluded_supports,
      int& representative_first_hit)
  {
    fluid_supports = 0;
    occluded_supports = 0;
    representative_first_hit = -1;
    for (int image = 0; image < eorder_t; ++image) {
      candidate_count(image) = 0;
      visible_count(image) = 0;
      if (imp_ninterp(image) <= 0) {
        for (int support = 0; support < NIP; ++support) {
          support_visible(image, support) = uint8_t(0);
          first_hit_element(image, support) = -1;
        }
        continue;
      }
      Point image_point;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        image_point[d] = imp_xyz(image, d);
      }
      for (int support = 0; support < NIP; ++support) {
        support_visible(image, support) = uint8_t(0);
        first_hit_element(image, support) = -1;

        int remainder = support;
        int ijk[AMREX_SPACEDIM];
        Point support_point;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          const int offset = remainder % (iorder_t + 1);
          remainder /= (iorder_t + 1);
          ijk[d] = imp_ijk(image, d) + offset;
          support_point[d] =
              prob_lo[d] + (Real(ijk[d]) + Real(0.5)) * dxyz[d];
        }
#if (AMREX_SPACEDIM == 3)
        const int kk = ijk[2];
#else
        const int kk = 0;
#endif
        if (ibMarkers(ijk[0], ijk[1], kk, 0) != 0) continue;
        ++fluid_supports;
        ++candidate_count(image);

        const WorldSegmentHit hit = firstSurfaceHitWorld(
            image_point, support_point, queries, transforms,
            geometry_offsets, geometry_count);
        if (!hit.hit()) {
          support_visible(image, support) = uint8_t(1);
          ++visible_count(image);
          continue;
        }

        first_hit_element(image, support) = hit.global_element;
        ++occluded_supports;
        if (representative_first_hit < 0 ||
            hit.global_element < representative_first_hit) {
          representative_first_hit = hit.global_element;
        }
      }
    }
  }

  template <typename IPDATA,
            int NIP = ipow(3, AMREX_SPACEDIM)>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static int evaluateBIConstrainedSupportBlock(
      const Array1D<int, 0, AMREX_SPACEDIM - 1>& support_base,
      const Array2D<Real, 0, 0, 0, AMREX_SPACEDIM - 1>& support_anchor,
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& wall_point,
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& ghost_point,
      const LocalFrame& localframe, int lev,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>& dxyz, Real length_scale,
      const Box& bxg, const Array4<uint8_t const>& ibMarkers,
      const IPDATA& ipData, int data_idx,
      const GpuArray<BVH4QueryView, MAX_NGEOM>& queries,
      const GpuArray<RigidTransform, MAX_NGEOM>& transforms,
      const GpuArray<int, MAX_NGEOM + 1>& geometry_offsets,
      int geometry_count,
      Array2D<int, 0, NIP - 1, 0, AMREX_SPACEDIM - 1>& support_indices,
      Array1D<uint8_t, 0, NIP - 1>& support_mask,
      Array1D<Real, 0, NIP - 1>& dirichlet_weights,
      Real& dirichlet_boundary_weight,
      Array1D<Real, 0, NIP - 1>& neumann_weights,
      Real& neumann_gradient_weight,
      Array1D<Real, 0, NIP - 1>& fluid_trace_weights,
      int& fluid_trace_order, Real& fluid_trace_condition,
      Real& fluid_trace_weight_l1, int& support_count,
      Real& accepted_condition, Real& weight_l1)
  {
#if (AMREX_SPACEDIM == 3)
    const int support_k = support_base(2);
#else
    const int support_k = 0;
#endif
    const bool support_in_box =
        check_interpolation_stencil<IPDATA, 2>(
            support_base(0), support_base(1), support_k, bxg, lev, ipData,
            data_idx, CheckMode::Silent);

    for (int support = 0; support < NIP; ++support) {
      int remainder = support;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int offset = remainder % 3;
        remainder /= 3;
        support_indices(support, d) =
            support_in_box ? support_base(d) + offset : -99;
      }
      dirichlet_weights(support) = Real(0.0);
      neumann_weights(support) = Real(0.0);
      fluid_trace_weights(support) = Real(0.0);
      support_mask(support) = uint8_t(0);
    }
    dirichlet_boundary_weight = Real(0.0);
    neumann_gradient_weight = Real(0.0);
    fluid_trace_order = -1;
    fluid_trace_condition = std::numeric_limits<Real>::infinity();
    fluid_trace_weight_l1 = std::numeric_limits<Real>::infinity();
    support_count = 0;
    accepted_condition = std::numeric_limits<Real>::infinity();
    weight_l1 = std::numeric_limits<Real>::infinity();
    if (!support_in_box) return -1;

    Array2D<uint8_t, 0, 0, 0, NIP - 1> visible{};
    Array2D<int, 0, 0, 0, NIP - 1> first_hit{};
    Array1D<int, 0, 0> target_count{};
    Array1D<int, 0, 0> candidate_count{};
    Array1D<int, 0, 0> visible_count{};
    target_count(0) = valid_mirror<2>(
        support_base(0), support_base(1), support_k, ibMarkers);
    int checked = 0;
    int occluded = 0;
    int representative_hit = -1;
    Array2D<int, 0, 0, 0, AMREX_SPACEDIM - 1> support_base_array{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      support_base_array(0, d) = support_base(d);
    }
    buildSupportVisibilityMask<1, 2>(
        visible, first_hit, candidate_count, visible_count, support_anchor,
        support_base_array, target_count, prob_lo, dxyz, ibMarkers, queries,
        transforms, geometry_offsets, geometry_count, checked, occluded,
        representative_hit);
    amrex::ignore_unused(checked, occluded, representative_hit, first_hit,
                         candidate_count, visible_count);

    for (int support = 0; support < NIP; ++support) {
      support_mask(support) = visible(0, support);
    }
    Array1D<Real, 0, AMREX_SPACEDIM - 1> normal{};
    Array1D<Real, 0, AMREX_SPACEDIM - 1> tangent1{};
    Array1D<Real, 0, AMREX_SPACEDIM - 1> tangent2{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      normal(d) = localframe.normal[d];
      tangent1(d) = localframe.tangent1[d];
#if (AMREX_SPACEDIM == 3)
      tangent2(d) = localframe.tangent2[d];
#endif
    }
    const int constrained_order = ibm_bi_constrained_gp_functionals<NIP>(
        support_indices, support_mask, wall_point, ghost_point, normal,
        tangent1, tangent2, prob_lo, dxyz, length_scale,
        support_visibility_condition_max, dirichlet_weights,
        dirichlet_boundary_weight, neumann_weights, neumann_gradient_weight,
        support_count, accepted_condition, weight_l1);
    if constexpr (!uses_euler_slip_curvature_pressure_closure) {
      return constrained_order;
    } else {
      fluid_trace_order = ibm_unconstrained_wall_trace_functional<NIP>(
          support_indices, support_mask, wall_point, normal, tangent1,
          tangent2, prob_lo, dxyz, length_scale,
          support_visibility_condition_max, fluid_trace_weights,
          fluid_trace_condition, fluid_trace_weight_l1);
      // The current-stage wall trace supplies only the Euler compatibility
      // datum. It must not alter the support selected for the primary
      // BI-constrained ghost reconstruction; otherwise changing the pressure
      // closure also changes every reconstructed primitive variable. An
      // invalid trace is handled as a local pzero fallback at evaluation time.
      return constrained_order;
    }
  }

  template <typename PointArray>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static bool isFirstRZAxisRingTarget(
      const PointArray& ghost_point,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>& dxyz)
  {
    if constexpr (!rz_axis_regular_entropy_jet ||
                  AMREX_SPACEDIM != 2) {
      amrex::ignore_unused(ghost_point, prob_lo, dxyz);
      return false;
    } else {
      const Real first_ring_centre =
          prob_lo[0] + Real(0.5) * dxyz[0];
      const Real tolerance =
          Real(64.0) * std::numeric_limits<Real>::epsilon() *
          amrex::max(
              Real(1.0),
              amrex::max(amrex::Math::abs(ghost_point(0)),
                         amrex::Math::abs(first_ring_centre)));
      return amrex::Math::abs(
                 ghost_point(0) - first_ring_centre) <= tolerance;
    }
  }

  template <int NIP, bool AnnularTarget = rz_annular_bic_cell_average,
            typename IndexArray, typename MaskArray, typename PointArray,
            typename VectorArray, typename WeightArray>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static int buildOneSidedEntropyJetFunctional(
      const IndexArray& support_indices, const MaskArray& support_mask,
      const PointArray& wall_point, const PointArray& target_point,
      const VectorArray& wall_normal, const VectorArray& tangent1,
      const VectorArray& tangent2,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
      Real length_scale, Real maximum_condition,
      bool enforce_axis_regularity, WeightArray& weights,
      WeightArray& linear_weights, Real& accepted_condition,
      Real& weight_l1, Real& linear_condition, Real& linear_weight_l1)
  {
    if constexpr (rz_axis_regular_entropy_jet &&
                  AMREX_SPACEDIM == 2) {
      if (enforce_axis_regularity) {
        return ibm_rz_axis_regular_ghost_jet_functional<
            NIP, AnnularTarget>(
            support_indices, support_mask, wall_point, target_point,
            prob_lo, dxyz, length_scale, maximum_condition, weights,
            linear_weights, accepted_condition, weight_l1,
            linear_condition, linear_weight_l1);
      }
    }
    return ibm_one_sided_ghost_jet_functional<NIP, AnnularTarget>(
        support_indices, support_mask, wall_point, target_point,
        wall_normal, tangent1, tangent2, prob_lo, dxyz, length_scale,
        maximum_condition, weights, linear_weights, accepted_condition,
        weight_l1, linear_condition, linear_weight_l1);
  }

  /// Build one stored BI-constrained functional. The ordinary reflected 3^D
  /// block is bit-for-bit the only candidate when support expansion is off.
  /// The opt-in search examines adjacent blocks only when that block has already
  /// reduced below quadratic order; every candidate is independently checked
  /// for active-fluid membership, line of sight, rank, and condition number.
  template <typename IPDATA,
            int NIP = ipow(3, AMREX_SPACEDIM)>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static int buildBIConstrainedGPFunctionals(
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& wall_point,
      const Array1D<Real, 0, AMREX_SPACEDIM - 1>& ghost_point,
      const LocalFrame& localframe, int lev,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>& dxyz, Real length_scale,
      const Box& bxg, const Array4<uint8_t const>& ibMarkers,
      const IPDATA& ipData, int data_idx,
      const GpuArray<BVH4QueryView, MAX_NGEOM>& queries,
      const GpuArray<RigidTransform, MAX_NGEOM>& transforms,
      const GpuArray<int, MAX_NGEOM + 1>& geometry_offsets,
      int geometry_count,
      Array2D<int, 0, NIP - 1, 0, AMREX_SPACEDIM - 1>& support_indices,
      Array1D<uint8_t, 0, NIP - 1>& support_mask,
      Array1D<Real, 0, NIP - 1>& dirichlet_weights,
      Real& dirichlet_boundary_weight,
      Array1D<Real, 0, NIP - 1>& neumann_weights,
      Real& neumann_gradient_weight,
      Array1D<Real, 0, NIP - 1>& fluid_trace_weights,
      int& fluid_trace_order, Real& fluid_trace_condition,
      Real& fluid_trace_weight_l1, int& support_count,
      Real& accepted_condition, Real& weight_l1,
      Array1D<Real, 0, NIP - 1>& fv_dirichlet_weights,
      Real& fv_dirichlet_boundary_weight, int& fv_dirichlet_order,
      Real& fv_dirichlet_condition, Real& fv_dirichlet_weight_l1,
      Array1D<Real, 0, NIP - 1>& entropy_jet_weights,
      Array1D<Real, 0, NIP - 1>& entropy_jet_linear_weights,
      int& entropy_jet_order, Real& entropy_jet_condition,
      Real& entropy_jet_weight_l1, Real& entropy_jet_linear_condition,
      Real& entropy_jet_linear_weight_l1,
      RZBICPointFunctionals<NIP>& rz_point_functionals)
  {
    const bool axis_regular_entropy =
        isFirstRZAxisRingTarget(ghost_point, prob_lo, dxyz);
    Array2D<Real, 0, 0, 0, AMREX_SPACEDIM - 1> support_anchor{};
    Array1D<int, 0, AMREX_SPACEDIM - 1> original_base{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      support_anchor(0, d) =
          Real(2.0) * wall_point(d) - ghost_point(d);
      original_base(d) = ibm_interp_base_index<2>(
          support_anchor(0, d), prob_lo[d], dxyz[d]);
    }

    int best_order = evaluateBIConstrainedSupportBlock<IPDATA, NIP>(
        original_base, support_anchor, wall_point, ghost_point, localframe,
        lev, prob_lo, dxyz, length_scale, bxg, ibMarkers, ipData, data_idx,
        queries, transforms, geometry_offsets, geometry_count, support_indices,
        support_mask, dirichlet_weights, dirichlet_boundary_weight, neumann_weights,
        neumann_gradient_weight, fluid_trace_weights, fluid_trace_order,
        fluid_trace_condition, fluid_trace_weight_l1, support_count,
        accepted_condition, weight_l1);

    Array1D<Real, 0, AMREX_SPACEDIM - 1> normal{};
    Array1D<Real, 0, AMREX_SPACEDIM - 1> tangent1{};
    Array1D<Real, 0, AMREX_SPACEDIM - 1> tangent2{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      normal(d) = localframe.normal[d];
      tangent1(d) = localframe.tangent1[d];
#if (AMREX_SPACEDIM == 3)
      tangent2(d) = localframe.tangent2[d];
#endif
    }
    if constexpr (cell_average_bic_normal_momentum) {
      fv_dirichlet_order =
          ibm_bi_cell_average_dirichlet_functional<NIP>(
              support_indices, support_mask, wall_point, ghost_point, normal,
              tangent1, tangent2, prob_lo, dxyz, length_scale,
              support_visibility_condition_max, fv_dirichlet_weights,
              fv_dirichlet_boundary_weight, fv_dirichlet_condition,
              fv_dirichlet_weight_l1);
    }
    if constexpr (one_sided_entropy_jet) {
      entropy_jet_order = buildOneSidedEntropyJetFunctional<NIP>(
          support_indices, support_mask, wall_point, ghost_point, normal,
          tangent1, tangent2, prob_lo, dxyz, length_scale,
          support_visibility_condition_max, axis_regular_entropy,
          entropy_jet_weights,
          entropy_jet_linear_weights, entropy_jet_condition,
          entropy_jet_weight_l1, entropy_jet_linear_condition,
          entropy_jet_linear_weight_l1);
    }

    int best_required_order = best_order;
    if constexpr (cell_average_bic_normal_momentum) {
      best_required_order = amrex::min(
          best_required_order, fv_dirichlet_order);
    }
    if constexpr (one_sided_entropy_jet) {
      best_required_order = amrex::min(
          best_required_order, entropy_jet_order);
    }

    if constexpr (expanded_bic_support) {
      if (best_required_order < 2) {
        // A 2-D unconstrained quadratic jet has six coefficients, whereas
        // the wall-constrained quadratic BIC has only five.  A cut 3x3 block
        // can therefore be complete for the primary BIC but rank deficient
        // for the auxiliary entropy jet.  Permit one extra outward block
        // translation only when that jet is part of the production state.
        constexpr int shift_radius = one_sided_entropy_jet ? 2 : 1;
        constexpr int shift_width = 2 * shift_radius + 1;
        constexpr int candidate_blocks =
            ipow(shift_width, AMREX_SPACEDIM);
        const Real projection_tolerance = Real(64.0) *
            std::numeric_limits<Real>::epsilon() * length_scale;
        for (int block = 0; block < candidate_blocks; ++block) {
          int remainder = block;
          bool shifted = false;
          Real outward_projection = Real(0.0);
          Array1D<int, 0, AMREX_SPACEDIM - 1> candidate_base{};
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int shift = remainder % shift_width - shift_radius;
            remainder /= shift_width;
            shifted = shifted || shift != 0;
            outward_projection +=
                Real(shift) * dxyz[d] * localframe.normal[d];
            candidate_base(d) = original_base(d) + shift;
          }
          if (!shifted || !(outward_projection > projection_tolerance)) {
            continue;
          }

          Array2D<int, 0, NIP - 1, 0, AMREX_SPACEDIM - 1>
              candidate_indices{};
          Array1D<uint8_t, 0, NIP - 1> candidate_support_mask{};
          Array1D<Real, 0, NIP - 1> candidate_dirichlet{};
          Array1D<Real, 0, NIP - 1> candidate_neumann{};
          Array1D<Real, 0, NIP - 1> candidate_fluid_trace{};
          Real candidate_boundary_weight = Real(0.0);
          Real candidate_gradient_weight = Real(0.0);
          int candidate_fluid_trace_order = -1;
          Real candidate_fluid_trace_condition =
              std::numeric_limits<Real>::infinity();
          Real candidate_fluid_trace_l1 =
              std::numeric_limits<Real>::infinity();
          int candidate_support_count = 0;
          Real candidate_condition =
              std::numeric_limits<Real>::infinity();
          Real candidate_l1 = std::numeric_limits<Real>::infinity();
          const int candidate_order =
              evaluateBIConstrainedSupportBlock<IPDATA, NIP>(
                  candidate_base, support_anchor, wall_point, ghost_point,
                  localframe, lev, prob_lo, dxyz, length_scale, bxg, ibMarkers,
                  ipData, data_idx, queries, transforms, geometry_offsets,
                  geometry_count, candidate_indices, candidate_support_mask,
                  candidate_dirichlet,
                  candidate_boundary_weight, candidate_neumann,
                  candidate_gradient_weight, candidate_fluid_trace,
                  candidate_fluid_trace_order,
                  candidate_fluid_trace_condition,
                  candidate_fluid_trace_l1, candidate_support_count,
                  candidate_condition, candidate_l1);
          if (candidate_order != 2) continue;

          Array1D<Real, 0, NIP - 1> candidate_fv_dirichlet{};
          Array1D<Real, 0, NIP - 1> candidate_entropy_jet{};
          Array1D<Real, 0, NIP - 1> candidate_entropy_jet_linear{};
          Real candidate_fv_boundary_weight = Real(0.0);
          Real candidate_fv_condition =
              std::numeric_limits<Real>::infinity();
          Real candidate_fv_l1 = std::numeric_limits<Real>::infinity();
          int candidate_fv_order = -1;
          if constexpr (cell_average_bic_normal_momentum) {
            candidate_fv_order =
                ibm_bi_cell_average_dirichlet_functional<NIP>(
                    candidate_indices, candidate_support_mask, wall_point,
                    ghost_point, normal, tangent1, tangent2, prob_lo, dxyz,
                    length_scale, support_visibility_condition_max,
                    candidate_fv_dirichlet, candidate_fv_boundary_weight,
                    candidate_fv_condition, candidate_fv_l1);
          }
          Real candidate_entropy_condition =
              std::numeric_limits<Real>::infinity();
          Real candidate_entropy_l1 =
              std::numeric_limits<Real>::infinity();
          Real candidate_entropy_linear_condition =
              std::numeric_limits<Real>::infinity();
          Real candidate_entropy_linear_l1 =
              std::numeric_limits<Real>::infinity();
          int candidate_entropy_order = -1;
          if constexpr (one_sided_entropy_jet) {
            candidate_entropy_order =
                buildOneSidedEntropyJetFunctional<NIP>(
                    candidate_indices, candidate_support_mask, wall_point,
                    ghost_point, normal, tangent1, tangent2, prob_lo, dxyz,
                    length_scale, support_visibility_condition_max,
                    axis_regular_entropy,
                    candidate_entropy_jet, candidate_entropy_jet_linear,
                    candidate_entropy_condition, candidate_entropy_l1,
                    candidate_entropy_linear_condition,
                    candidate_entropy_linear_l1);
          }

          int candidate_required_order = candidate_order;
          if constexpr (cell_average_bic_normal_momentum) {
            candidate_required_order = amrex::min(
                candidate_required_order, candidate_fv_order);
          }
          if constexpr (one_sided_entropy_jet) {
            candidate_required_order = amrex::min(
                candidate_required_order, candidate_entropy_order);
          }
          const Real l1_tolerance = Real(64.0) *
              std::numeric_limits<Real>::epsilon() *
              amrex::max(Real(1.0),
                         amrex::max(amrex::Math::abs(candidate_l1),
                                    amrex::Math::abs(weight_l1)));
          const bool better =
              candidate_required_order > best_required_order ||
              (candidate_required_order == best_required_order &&
               (candidate_order > best_order ||
                (candidate_order == best_order &&
                 (candidate_l1 < weight_l1 - l1_tolerance ||
                  (amrex::Math::abs(candidate_l1 - weight_l1) <=
                       l1_tolerance &&
                   candidate_condition < accepted_condition)))));
          if (!better) continue;

          best_order = candidate_order;
          best_required_order = candidate_required_order;
          support_indices = candidate_indices;
          support_mask = candidate_support_mask;
          dirichlet_weights = candidate_dirichlet;
          neumann_weights = candidate_neumann;
          fluid_trace_weights = candidate_fluid_trace;
          fv_dirichlet_weights = candidate_fv_dirichlet;
          entropy_jet_weights = candidate_entropy_jet;
          entropy_jet_linear_weights = candidate_entropy_jet_linear;
          dirichlet_boundary_weight = candidate_boundary_weight;
          neumann_gradient_weight = candidate_gradient_weight;
          fv_dirichlet_boundary_weight = candidate_fv_boundary_weight;
          fluid_trace_order = candidate_fluid_trace_order;
          fluid_trace_condition = candidate_fluid_trace_condition;
          fluid_trace_weight_l1 = candidate_fluid_trace_l1;
          fv_dirichlet_order = candidate_fv_order;
          fv_dirichlet_condition = candidate_fv_condition;
          fv_dirichlet_weight_l1 = candidate_fv_l1;
          entropy_jet_order = candidate_entropy_order;
          entropy_jet_condition = candidate_entropy_condition;
          entropy_jet_weight_l1 = candidate_entropy_l1;
          entropy_jet_linear_condition =
              candidate_entropy_linear_condition;
          entropy_jet_linear_weight_l1 = candidate_entropy_linear_l1;
          support_count = candidate_support_count;
          accepted_condition = candidate_condition;
          weight_l1 = candidate_l1;
        }
      }
    }
    if constexpr (bic_point_target_functionals) {
      // Primitive averages cannot be passed through a nonlinear EOS to form
      // a conservative cell average. Cache point-evaluation functionals for
      // tensor-product 2x2 Gauss rules. R-Z uses three adjacent radial cells
      // for centre-state recovery; Cartesian no-slip uses the four targets in
      // the real GP cell.
      constexpr Real gauss_node =
          Real(0.57735026918962576450914878050195746);
      constexpr int target_count =
          rz_annular_bic_cell_average ? IBM_RZ_BIC_TARGETS
                                      : IBM_RZ_BIC_GAUSS_POINTS;
      for (int target = 0; target < target_count; ++target) {
        auto target_point = ghost_point;
        if constexpr (rz_annular_bic_cell_average) {
          const int radial_slot = target / IBM_RZ_BIC_GAUSS_POINTS;
          const int quadrature_point =
              target % IBM_RZ_BIC_GAUSS_POINTS;
          const int radial_cell_offset = radial_slot - 1;
          const Real radial_sign =
              (quadrature_point & 1) == 0 ? Real(-1.0) : Real(1.0);
          const Real axial_sign =
              (quadrature_point & 2) == 0 ? Real(-1.0) : Real(1.0);
          target_point(0) +=
              (Real(radial_cell_offset) +
               Real(0.5) * radial_sign * gauss_node) * dxyz[0];
          target_point(1) +=
              Real(0.5) * axial_sign * gauss_node * dxyz[1];
        } else {
          const Real x_sign =
              (target & 1) == 0 ? Real(-1.0) : Real(1.0);
          const Real y_sign =
              (target & 2) == 0 ? Real(-1.0) : Real(1.0);
          target_point(0) +=
              Real(0.5) * x_sign * gauss_node * dxyz[0];
          target_point(1) +=
              Real(0.5) * y_sign * gauss_node * dxyz[1];
        }

        Array1D<Real, 0, NIP - 1> point_dirichlet{};
        Array1D<Real, 0, NIP - 1> point_neumann{};
        Real point_boundary_weight = Real(0.0);
        Real point_gradient_weight = Real(0.0);
        int point_support_count = 0;
        Real point_condition = std::numeric_limits<Real>::infinity();
        Real point_l1 = std::numeric_limits<Real>::infinity();
        const int point_order =
            ibm_bi_constrained_gp_functionals<NIP, false>(
                support_indices, support_mask, wall_point, target_point,
                normal, tangent1, tangent2, prob_lo, dxyz, length_scale,
                support_visibility_condition_max, point_dirichlet,
                point_boundary_weight, point_neumann,
                point_gradient_weight, point_support_count, point_condition,
                point_l1);
        if (point_order != best_order) return -1;

        Array1D<Real, 0, NIP - 1> point_fv_dirichlet{};
        Real point_fv_boundary_weight = Real(0.0);
        Real point_fv_condition = std::numeric_limits<Real>::infinity();
        Real point_fv_l1 = std::numeric_limits<Real>::infinity();
        int point_fv_order = best_order;
        if constexpr (cell_average_bic_normal_momentum) {
          point_fv_order =
              ibm_bi_cell_average_dirichlet_functional<NIP, false>(
                  support_indices, support_mask, wall_point, target_point,
                  normal, tangent1, tangent2, prob_lo, dxyz, length_scale,
                  support_visibility_condition_max, point_fv_dirichlet,
                  point_fv_boundary_weight, point_fv_condition,
                  point_fv_l1);
          if (point_fv_order != fv_dirichlet_order) return -1;
        }

        Array1D<Real, 0, NIP - 1> point_entropy{};
        Array1D<Real, 0, NIP - 1> point_entropy_linear{};
        Real point_entropy_condition =
            std::numeric_limits<Real>::infinity();
        Real point_entropy_l1 = std::numeric_limits<Real>::infinity();
        Real point_entropy_linear_condition =
            std::numeric_limits<Real>::infinity();
        Real point_entropy_linear_l1 =
            std::numeric_limits<Real>::infinity();
        int point_entropy_order = best_order;
        if constexpr (one_sided_entropy_jet) {
          point_entropy_order =
              buildOneSidedEntropyJetFunctional<NIP, false>(
                  support_indices, support_mask, wall_point, target_point,
                  normal, tangent1, tangent2, prob_lo, dxyz, length_scale,
                  support_visibility_condition_max, axis_regular_entropy,
                  point_entropy,
                  point_entropy_linear, point_entropy_condition,
                  point_entropy_l1, point_entropy_linear_condition,
                  point_entropy_linear_l1);
          if (point_entropy_order != entropy_jet_order) return -1;
        }

        for (int support = 0; support < NIP; ++support) {
          rz_point_functionals.dirichlet_weights(target, support) =
              point_dirichlet(support);
          rz_point_functionals.neumann_weights(target, support) =
              point_neumann(support);
          rz_point_functionals.fv_dirichlet_weights(target, support) =
              point_fv_dirichlet(support);
          rz_point_functionals.entropy_jet_weights(target, support) =
              point_entropy(support);
          rz_point_functionals.entropy_jet_linear_weights(target, support) =
              point_entropy_linear(support);
        }
        rz_point_functionals.dirichlet_boundary_weight(target) =
            point_boundary_weight;
        rz_point_functionals.neumann_gradient_weight(target) =
            point_gradient_weight;
        rz_point_functionals.fv_dirichlet_boundary_weight(target) =
            point_fv_boundary_weight;
      }
    } else {
      amrex::ignore_unused(rz_point_functionals);
    }
    return best_order;
  }
#endif

  /**
   * \brief Bias a centred 3x3 WLS block along the owning polygon edge.
   *
   * The image point and its wall-normal distance are unchanged.  Only the
   * Cartesian candidate block used to fit the image-point value is translated
   * by one cell in the direction of the owning edge that points away from the
   * nearest annotated corner.  The translated block is accepted only when it
   * remains inside the owning FAB and supports a full-rank quadratic fit.
   * Consequently this operation changes the selected trace near a singular
   * polygon vertex without turning the normal reconstruction into an oblique
   * image-point construction.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static Real sourceDirectedTangent(const LocalFrame& localframe, int d)
  {
    // Geometry edges are stored source -> target.  convert_inout() reverses
    // tangent1 together with the normal for interior-is-fluid geometries, but
    // the sharp-feature arc distances retain the original edge numbering.
    return (interior_is_solid ? Real(1.0) : Real(-1.0)) *
           localframe.tangent1[d];
  }

  template <int eorder_t, int iorder_t, typename IPDATA>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static int applySharpFeatureTangentialStencil(
      const SurfElem& edge, Real coordinate_from_source, Real h,
      const LocalFrame& localframe, int lev,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>& dxyz, const Box& bxg,
      const Array4<uint8_t const>& ibMarkers, const IPDATA& ipData,
      int data_idx,
      const Array2D<Real, 0, eorder_t - 1, 0,
                    AMREX_SPACEDIM - 1>& imp_xyz,
      Array2D<int, 0, eorder_t - 1, 0,
              AMREX_SPACEDIM - 1>& imp_ijk,
      Array1D<int, 0, eorder_t - 1>& imp_ninterp)
  {
#if (AMREX_SPACEDIM == 2)
    if constexpr (sharp_feature_tangential_reconstruction) {
      static_assert(iorder_t == 2,
                    "one-sided sharp-feature stencil requires 3x3 WLS");
      if (edge.sharp_distance_source < Real(0.0) ||
          edge.sharp_distance_target < Real(0.0) || !(h > Real(0.0))) {
        return 0;
      }

      const Real s = amrex::max(
          Real(0.0), amrex::min(edge.measure, coordinate_from_source));
      const Real distance_from_source = edge.sharp_distance_source + s;
      const Real distance_from_target =
          edge.sharp_distance_target + edge.measure - s;
      const Real feature_distance =
          amrex::min(distance_from_source, distance_from_target);
      if (!(feature_distance < sharp_feature_tangential_radius_cells * h)) {
        return 0;
      }
      const Real away_sign =
          distance_from_source <= distance_from_target ? Real(1.0)
                                                       : Real(-1.0);

      // Select the neighbouring Cartesian block whose physical translation
      // is most closely aligned with the owner-edge tangent.  Searching all
      // eight neighbours handles diagonal edges without introducing a
      // preferred coordinate direction.
      int shift_i = 0;
      int shift_j = 0;
      Real best_alignment = Real(-2.0);
      for (int sj = -1; sj <= 1; ++sj) {
        for (int si = -1; si <= 1; ++si) {
          if (si == 0 && sj == 0) continue;
          const Real vx = Real(si) * dxyz[0];
          const Real vy = Real(sj) * dxyz[1];
          const Real length = std::sqrt(vx * vx + vy * vy);
          const Real alignment =
              away_sign *
              (vx * sourceDirectedTangent(localframe, 0) +
               vy * sourceDirectedTangent(localframe, 1)) /
              length;
          if (alignment > best_alignment) {
            best_alignment = alignment;
            shift_i = si;
            shift_j = sj;
          }
        }
      }
      if (!(best_alignment > Real(0.0))) return 0;

      constexpr int interp_threshold =
          is_gpData_t<IPDATA>::value ? INTERP_THRESHOLD_GP
                                     : INTERP_THRESHOLD_SURF;
      int shifted = 0;
      for (int iim = 0; iim < eorder_t; ++iim) {
        Array1D<int, 0, AMREX_SPACEDIM - 1> candidate_ijk;
        Array1D<Real, 0, AMREX_SPACEDIM - 1> candidate_xyz;
        candidate_ijk(0) = imp_ijk(iim, 0) + shift_i;
        candidate_ijk(1) = imp_ijk(iim, 1) + shift_j;
        candidate_xyz(0) = imp_xyz(iim, 0);
        candidate_xyz(1) = imp_xyz(iim, 1);

        const bool in_box =
            check_interpolation_stencil<IPDATA, iorder_t>(
                candidate_ijk(0), candidate_ijk(1), 0, bxg, lev, ipData,
                data_idx, CheckMode::Silent);
        if (!in_box) continue;
        const int n_fluid = valid_mirror<iorder_t>(
            candidate_ijk(0), candidate_ijk(1), 0, ibMarkers);
        if (n_fluid < interp_threshold) continue;
        const int fit_order =
            ibm_interp_candidate_fit_order<iorder_t>(
                candidate_ijk, candidate_xyz, prob_lo, dxyz, ibMarkers);
        if (fit_order != 2) continue;

        imp_ijk(iim, 0) = candidate_ijk(0);
        imp_ijk(iim, 1) = candidate_ijk(1);
        imp_ninterp(iim) = n_fluid;
        ++shifted;
      }
      return shifted;
    }
#endif
    amrex::ignore_unused(edge, coordinate_from_source, h, localframe, lev,
                         prob_lo, dxyz, bxg, ibMarkers, ipData, data_idx,
                         imp_xyz, imp_ijk, imp_ninterp);
    return 0;
  }

  /**
   * \brief Sanity-check that surface normals point from solid to fluid.
   *
   * Samples N faces per geometry, computes a probe point at
   * centroid + eps*normal (body frame), and tests inside/outside.  If the
   * probe lands inside the body, the normal is pointing the wrong way,
   * which usually means the user's `interior_is_solid` parameter does not
   * match the mesh's actual orientation (e.g., STL exported with reversed
   * face winding, or polygon file in the wrong CCW/CW convention).
   *
   * Warns on any failure; aborts if the majority of samples fail (clear
   * miscofiguration vs. occasional boundary-grazing false positive).
   *
   * Called once at the end of rebuildGeometryData() — after convert_inout
   * has had its chance to flip normals based on `interior_is_solid`.
   */
  void verify_normal_orientation()
  {
    BL_PROFILE("IBM::verify_normal_orientation");
    constexpr int N_SAMPLES_PER_GEOM = 20;
    constexpr Real EPS_FACTOR        = Real(0.01);  // 1% of body diagonal

    for (int ii = 0; ii < ngeom; ++ii) {
      const int n_faces = geom_offsets[ii + 1] - geom_offsets[ii];
      if (n_faces == 0) continue;

      // Length scale: body-frame bounding box diagonal.
      Real diag2 = Real(0.0);
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real ext = bbox_max_d(bbox_body_a[ii], d)
                       - bbox_min_d(bbox_body_a[ii], d);
        diag2 += ext * ext;
      }
      const Real eps = EPS_FACTOR * std::sqrt(diag2);

      const int n_samples = std::min(N_SAMPLES_PER_GEOM, n_faces);
      const int step      = std::max(1, n_faces / n_samples);

      int n_failures = 0;
      int n_checked  = 0;
      for (int local_f = 0; local_f < n_faces; local_f += step) {
        const int f_idx = geom_offsets[ii] + local_f;
        const auto& se  = SurfElem_a[f_idx];
        const auto& lf  = LocalFrame_a[f_idx];

        Point probe;
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        probe = Point(se.centroid[0] + eps * lf.normal[0],
                      se.centroid[1] + eps * lf.normal[1]);
#else
        probe = Point(se.centroid[0] + eps * lf.normal[0],
                      se.centroid[1] + eps * lf.normal[1],
                      se.centroid[2] + eps * lf.normal[2]);
#endif
#else
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          probe[d] = se.centroid[d] + eps * lf.normal[d];
        }
#endif

        BoundedSide side = (*inout_fa[ii])(probe);
        if (side == BoundedSide::Inside) ++n_failures;
        ++n_checked;
      }

      if (n_failures > 0) {
        amrex::Print()
            << "[IBM] WARNING: normal-orientation sanity check failed for geometry " << ii
            << ": " << n_failures << "/" << n_checked
            << " sample probes (centroid + eps*normal) landed inside the body.\n"
            << "       eps = " << eps << " (body-frame, ~"
            << int(EPS_FACTOR * Real(100)) << "% of diag).\n"
            << "       Likely cause: the 'interior_is_solid' parameter does not match\n"
            << "       the mesh's actual outward-normal orientation. STL files with\n"
            << "       inverted face windings, or polygon files in the wrong CCW/CW\n"
            << "       convention, will trigger this. Visualise the mesh in ParaView\n"
            << "       and verify face normals point away from the solid interior.\n";

        if (n_failures > n_checked / 2) {
          amrex::Abort("[IBM] verify_normal_orientation: majority of samples failed; "
                       "check `interior_is_solid` parameter and mesh orientation.");
        }
      }
    }
  }

  /**
   * \brief Full geometry rebuild (for deformable bodies or re-initialization).
   *
   * This is the expensive path that rebuilds BVH, InsideTester, LocalFrame,
   * and SurfElem from scratch.  For rigid-body FSI, use updateRigidTransform()
   * instead — it is O(1) per geometry.
   */
  void rebuildGeometryData()
  {
    LocalFrame_a.clear();
    SurfElem_a.clear();

    for (int i = 0; i < ngeom; i++) {
#ifdef AMREX_USE_CGAL
      tree_a[i].clear();
#if (AMREX_SPACEDIM == 2)
      tree_a[i].insert(geom_a[i].edges_begin(), geom_a[i].edges_end());
#elif (AMREX_SPACEDIM == 3)
      tree_a[i].insert(faces(geom_a[i]).first, faces(geom_a[i]).second, geom_a[i]);
#endif
      tree_a[i].build();

      inout_fa[i] = std::make_unique<inside_t>(geom_a[i]);
#else
      bvh_a[i].build(geom_a[i]);

      inout_fa[i] = std::make_unique<inside_t>(geom_a[i], bvh_a[i]);
#endif

#if defined(AMREX_USE_CGAL) && (AMREX_SPACEDIM == 3)
      bbox_body_a[i] = PMP::bbox(geom_a[i]);
#else
      bbox_body_a[i] = geom_a[i].bbox();
#endif
      bbox_a[i] = transform_a[i].transform_bbox(bbox_body_a[i]);

      geom_offsets[i] = static_cast<int>(LocalFrame_a.size());
#ifdef AMREX_USE_CGAL
      build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                           idxmap_a[i], geom_offsets[i], i);
#elif defined(AMREX_USE_GPU)
      build_geometry_cache_gpu(geom_a[i], SurfElem_a, LocalFrame_a, i);
#else
      build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                           geom_offsets[i], i);
#endif
    }
    geom_offsets[ngeom] = static_cast<int>(LocalFrame_a.size());

    annotateSharpFeatures2D();

    if (!interior_is_solid) {
      convert_inout(LocalFrame_a);
    }

    // Defensive check: detect mismatched `interior_is_solid` setting or
    // reversed mesh winding before downstream IBM machinery reads
    // LocalFrame_a / inout_fa.
    verify_normal_orientation();
  }

  /**
   * \brief Compute solid/fluid markers and identify ghost points on a level.
   *
   * Per cell, writes two uint8_t components into the IBMultiFab:
   *   comp 0 — solid marker:  0 = fluid, ii+1 = solid in geometry ii
   *   comp 1 — ghost marker:  0 = not a ghost point, otherwise carries the
   *                           geometry id of the owning solid cell
   *
   * Pipeline (per FAB):
   *   FAST-SKIP : reject FABs whose working region misses every body bbox.
   *   Step 1    : solid markers via fast bbox-reject + BVH inside test.
   *   Step 2    : ghost markers + ghost-point count (fused reduction).
   *
   * Per-FAB ghost-point counts feed the level-wide CSR GPStore allocation
   * at the end of the function.
   *
   * \param lev AMR level to process.
   */
  void computeMarkers(int lev)
  {
    BL_PROFILE("IBM::computeMarkers");
    const Real timing_total_start = beginPerformanceTiming();
    const bool timing_detail = detailedPerformanceTimingEnabled();
    Real timing_inout_seconds = Real(0.0);
    Real timing_ghost_seconds = Real(0.0);
    Long timing_inout_cells = Long(0);
    Long timing_ghost_cells = Long(0);
    Long timing_marker_query_cells = Long(0);

    // Step 2 reads markers in grow(bx, GP_BOX_EXTRA).  The marker allocation
    // now covers the farther image-point halo; retain this lower-bound check
    // because every RHS configuration still guarantees cls_t::NGHOST.
    static_assert(GP_BOX_EXTRA <= cls_t::NGHOST,
        "GP_BOX_EXTRA must not exceed cls_t::NGHOST: Step 2 would read "
        "marker cells that Step 1 never initialised.");

    // ---- Level-wide handles ----------------------------------------------
    auto& mfab        = *bmf_a[lev];
    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto dx_lev  = dx_a[lev];
    const int marker_nghost = interpolationMarkerNghost(lev);
    const Box physical_domain = amr_p->Geom(lev).Domain();
    GpuArray<int, AMREX_SPACEDIM> domain_lo{};
    GpuArray<int, AMREX_SPACEDIM> domain_hi{};
    GpuArray<int, AMREX_SPACEDIM> periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      domain_lo[d] = physical_domain.smallEnd(d);
      domain_hi[d] = physical_domain.bigEnd(d);
      periodic[d] = amr_p->Geom(lev).isPeriodic(d) ? 1 : 0;
    }

    // ---- Per-FAB GP counts (drive level-wide CSR allocation below) -------
    const int nfabs_local = mfab.local_size();
    Vector<int> gp_counts(nfabs_local, 0);
    int ifab_local = 0;

#ifdef AMREX_USE_GPU
    // Batched per-FAB GP counters: Step 2 accumulates into one device array
    // and a single D2H copy after the MFIter loop replaces the previous
    // blocking ReduceData::value() (4 B D2H + stream sync) per FAB.
    // Pre-zeroed from gp_counts (all zeros) so fast-skipped FABs, whose
    // slot no kernel ever touches, read back 0.
    Gpu::DeviceVector<int> d_gp_counts(nfabs_local);
    Gpu::copyAsync(Gpu::hostToDevice, gp_counts.begin(), gp_counts.end(),
                   d_gp_counts.begin());
    int* const p_gp_counts = d_gp_counts.data();
#endif

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      const Box&  bx        = mfi.tilebox();
      const auto& ibMarkers = mfab.array(mfi);

      // ===================================================================
      // FAST-SKIP — short-circuit FABs that are entirely freestream.
      //
      // If the FAB's working region (valid box + the ghost cells consumed
      // downstream) is disjoint from every body's world-frame AABB, every
      // cell must be fluid; setVal(0) is mathematically equivalent to the
      // full kernel but skips the GPU launch and per-cell BVH query.
      //
      // Working region grows valid box by max(NGHOST, GP_BOX_EXTRA) + 2:
      //   NGHOST       — stencil width read by compute_rhs
      //   GP_BOX_EXTRA — extra growth for ghost-point detection (Step 2)
      //   +2           — safety margin for moving-geometry bbox inflation
      //
      // bbox_a[ii] is in the world frame: built once for static geometry,
      // refreshed by the caller every regrid for FSI before this runs.
      // ===================================================================
      {
        const int check_grow = amrex::max(marker_nghost, GP_BOX_EXTRA) + 2;
        const Box bxcheck = amrex::grow(bx, check_grow);

        bool any_touches = false;
        for (int ii = 0; ii < ngeom; ++ii) {
          bool intersects = true;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const Real box_lo_d = prob_lo[d] +  bxcheck.smallEnd(d)      * dx_lev[d];
            const Real box_hi_d = prob_lo[d] + (bxcheck.bigEnd(d) + 1)   * dx_lev[d];
            if (box_hi_d < bbox_min_d(bbox_a[ii], d) ||
                box_lo_d > bbox_max_d(bbox_a[ii], d)) {
              intersects = false;
              break;
            }
          }
          if (intersects) { any_touches = true; break; }
        }

        if (!any_touches) {
          // Pure freestream: clear both marker components over the FAB's
          // full allocated region (valid + ghost), then move on.
          mfab.get(mfi).template setVal<amrex::RunOn::Device>(uint8_t(0));
          gp_counts[ifab_local] = 0;
          continue;
        }
        if (performanceTimingEnabled()) {
          timing_marker_query_cells +=
              amrex::grow(bx, marker_nghost).numPts();
        }
      }

#ifdef AMREX_USE_GPU
      // ===================================================================
      // GPU path
      // ===================================================================
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
          "ngeom exceeds MAX_NGEOM; raise MAX_NGEOM in ibm_containers.h");
      const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);

      // Pre-pack trivially-copyable POD views for device capture.
      GpuArray<InsideTesterView, MAX_NGEOM> inout_views;
      GpuArray<AABB,             MAX_NGEOM> bboxes;
      GpuArray<RigidTransform,   MAX_NGEOM> transforms;
      for (int ii = 0; ii < ngeom_local; ++ii) {
        inout_views[ii] = inout_fa[ii]->view();
        bboxes[ii]      = bbox_a[ii];          // world-frame bbox
        transforms[ii]  = transform_a[ii];     // body→world transform
      }

      // --- Step 1: solid markers over the complete IBM interpolation halo -
      // Query points are in world frame: fast-reject vs world bbox, then
      // inverse-transform to body frame for the BVH inside/outside test.
      const Box& bxg = amrex::grow(bx, marker_nghost);
      amrex::ParallelFor(bxg,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        ibMarkers(i, j, k, 0) = uint8_t(0);
        ibMarkers(i, j, k, 1) = uint8_t(0);

        Point gridpoint = make_grid_point(prob_lo, dx_lev, i, j, k);
        for (int ii = 0; ii < ngeom_local; ++ii) {
          if (!bbox_contains(bboxes[ii], gridpoint)) continue;
          Point gp_body = transforms[ii].to_body(gridpoint);
          BoundedSide result = inout_views[ii](gp_body);
          if (result == BoundedSide::Inside || result == BoundedSide::OnBoundary) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      });

      // --- Step 2: ghost markers (comp 1) + count ------------------------
      // A solid cell is a ghost point iff at least one neighbour within
      // ±ghost_layers along an axis is fluid. Comp 1 stores the geometry id
      // (== comp 0) on ghost points; non-ghost cells keep the comp-1 zero
      // written by Step 1 (relies on GP_BOX_EXTRA <= NGHOST static_assert).
      // Count via per-FAB device atomic slot instead of a fused ReduceOps:
      // ReduceData::value() forced a blocking D2H + stream sync per FAB;
      // the slots are read back once for all FABs after the loop.
      {
        constexpr int gl   = ghost_layers;
        const Box&    bxgp = amrex::grow(bx, GP_BOX_EXTRA);
        int* const p_count_fab = p_gp_counts + ifab_local;

        amrex::ParallelFor(bxgp,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (!ibMarkers(i, j, k, 0)) return;  // fluid: not a candidate

          bool ghost = false;
          int layer = 0;
          for (int d = 1; d <= gl; ++d) {
            const bool i_here = periodic[0] ||
                (i >= domain_lo[0] && i <= domain_hi[0]);
            const bool j_here = periodic[1] ||
                (j >= domain_lo[1] && j <= domain_hi[1]);
            const bool i_minus = periodic[0] ||
                (i-d >= domain_lo[0] && i-d <= domain_hi[0]);
            const bool i_plus = periodic[0] ||
                (i+d >= domain_lo[0] && i+d <= domain_hi[0]);
            const bool j_minus = periodic[1] ||
                (j-d >= domain_lo[1] && j-d <= domain_hi[1]);
            const bool j_plus = periodic[1] ||
                (j+d >= domain_lo[1] && j+d <= domain_hi[1]);
#if (AMREX_SPACEDIM == 3)
            const bool k_here = periodic[2] ||
                (k >= domain_lo[2] && k <= domain_hi[2]);
            const bool k_minus = periodic[2] ||
                (k-d >= domain_lo[2] && k-d <= domain_hi[2]);
            const bool k_plus = periodic[2] ||
                (k+d >= domain_lo[2] && k+d <= domain_hi[2]);
            ghost = ghost ||
                (i_minus && j_here && k_here &&
                 !ibMarkers(i-d, j, k, 0)) ||
                (i_plus && j_here && k_here &&
                 !ibMarkers(i+d, j, k, 0)) ||
                (i_here && j_minus && k_here &&
                 !ibMarkers(i, j-d, k, 0)) ||
                (i_here && j_plus && k_here &&
                 !ibMarkers(i, j+d, k, 0)) ||
                (i_here && j_here && k_minus &&
                 !ibMarkers(i, j, k-d, 0)) ||
                (i_here && j_here && k_plus &&
                 !ibMarkers(i, j, k+d, 0));
#else
            ghost = ghost ||
                (i_minus && j_here && !ibMarkers(i-d, j, k, 0)) ||
                (i_plus && j_here && !ibMarkers(i+d, j, k, 0)) ||
                (i_here && j_minus && !ibMarkers(i, j-d, k, 0)) ||
                (i_here && j_plus && !ibMarkers(i, j+d, k, 0));
#endif
            if (ghost) {
              layer = d;
              break;
            }
          }

          if (ghost) {
            ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
            Gpu::Atomic::Add(p_count_fab, 1);
          }
        });
      }

#else
      // ===================================================================
      // CPU path — semantically identical to GPU. CPU additionally emits
      // an OnBoundary diagnostic via IB_WarnOnBoundary (cannot print on GPU).
      // ===================================================================

      // --- Step 1: solid markers over the complete IBM interpolation halo -
      const Box marker_box = amrex::grow(bx, marker_nghost);
      const Real inout_start =
          timing_detail ? amrex::second() : Real(0.0);
      amrex::LoopOnCpu(marker_box,
        [&](int i, int j, int k)
      {
        ibMarkers(i, j, k, 0) = uint8_t(0);
        ibMarkers(i, j, k, 1) = uint8_t(0);

        Point gridpoint = make_grid_point(prob_lo, dx_lev, i, j, k);
        for (int ii = 0; ii < ngeom; ++ii) {
          if (!bbox_contains(bbox_a[ii], gridpoint)) continue;
          Point gp_body = transform_a[ii].to_body(gridpoint);
          inside_t& inside = *inout_fa[ii];
          BoundedSide result = inside(gp_body);
          IB_WarnOnBoundary(ii, lev, i, j, k, result, gridpoint);
          if (result == BoundedSide::Inside || result == BoundedSide::OnBoundary) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      });
      if (timing_detail) {
        timing_inout_seconds += amrex::second() - inout_start;
        timing_inout_cells += marker_box.numPts();
      }

      // --- Step 2: ghost markers (comp 1) + count ------------------------
      int ngps_fab = 0;
      const Box ghost_box = amrex::grow(bx, GP_BOX_EXTRA);
      const Real ghost_start =
          timing_detail ? amrex::second() : Real(0.0);
      amrex::LoopOnCpu(ghost_box,
        [&](int i, int j, int k)
      {
        if (!ibMarkers(i, j, k, 0)) return;  // fluid: not a candidate

        bool ghost = false;
        int layer = 0;
        for (int d = 1; d <= ghost_layers; ++d) {
          const bool i_here = periodic[0] ||
              (i >= domain_lo[0] && i <= domain_hi[0]);
          const bool j_here = periodic[1] ||
              (j >= domain_lo[1] && j <= domain_hi[1]);
          const bool i_minus = periodic[0] ||
              (i-d >= domain_lo[0] && i-d <= domain_hi[0]);
          const bool i_plus = periodic[0] ||
              (i+d >= domain_lo[0] && i+d <= domain_hi[0]);
          const bool j_minus = periodic[1] ||
              (j-d >= domain_lo[1] && j-d <= domain_hi[1]);
          const bool j_plus = periodic[1] ||
              (j+d >= domain_lo[1] && j+d <= domain_hi[1]);
#if (AMREX_SPACEDIM == 3)
          const bool k_here = periodic[2] ||
              (k >= domain_lo[2] && k <= domain_hi[2]);
          const bool k_minus = periodic[2] ||
              (k-d >= domain_lo[2] && k-d <= domain_hi[2]);
          const bool k_plus = periodic[2] ||
              (k+d >= domain_lo[2] && k+d <= domain_hi[2]);
          ghost = ghost ||
              (i_minus && j_here && k_here &&
               !ibMarkers(i-d, j, k, 0)) ||
              (i_plus && j_here && k_here &&
               !ibMarkers(i+d, j, k, 0)) ||
              (i_here && j_minus && k_here &&
               !ibMarkers(i, j-d, k, 0)) ||
              (i_here && j_plus && k_here &&
               !ibMarkers(i, j+d, k, 0)) ||
              (i_here && j_here && k_minus &&
               !ibMarkers(i, j, k-d, 0)) ||
              (i_here && j_here && k_plus &&
               !ibMarkers(i, j, k+d, 0));
#else
          ghost = ghost ||
              (i_minus && j_here && !ibMarkers(i-d, j, k, 0)) ||
              (i_plus && j_here && !ibMarkers(i+d, j, k, 0)) ||
              (i_here && j_minus && !ibMarkers(i, j-d, k, 0)) ||
              (i_here && j_plus && !ibMarkers(i, j+d, k, 0));
#endif
          if (ghost) {
            layer = d;
            break;
          }
        }

        if (ghost) {
          ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
          ++ngps_fab;
        }
        // Non-ghost solid cells keep the comp-1 zero from Step 1 init
        // (relies on GP_BOX_EXTRA <= NGHOST static_assert).
      });
      if (timing_detail) {
        timing_ghost_seconds += amrex::second() - ghost_start;
        timing_ghost_cells += ghost_box.numPts();
      }
      gp_counts[ifab_local] = ngps_fab;

#endif // AMREX_USE_GPU

    } // end MFIter

#ifdef AMREX_USE_GPU
    // Single batched readback of all per-FAB GP counts (one sync for the
    // whole level; replaces one blocking readback per FAB).
    Gpu::copy(Gpu::deviceToHost, d_gp_counts.begin(), d_gp_counts.end(),
              gp_counts.begin());
#endif

    // ---- Build the level-wide flattened GPStore (CSR layout) -------------
    gpstore_a[lev].allocate(gp_counts);
    if constexpr (use_bi_constrained_shared_gp) {
      gpstore_a[lev].allocate_bi_constrained();
      if constexpr (bic_point_target_functionals) {
        gpstore_a[lev].allocate_rz_bic_point_functionals();
      }
    }

    if (timing_detail) {
      reportPerformanceDuration("markers_inout_test", lev,
                                timing_inout_seconds,
                                timing_inout_cells, false);
      reportPerformanceDuration("markers_ghost_detect", lev,
                                timing_ghost_seconds,
                                timing_ghost_cells, false);
    }
    finishPerformanceTiming("markers_total", lev, timing_total_start,
                            timing_marker_query_cells, false);
  }

  /**
   * \brief Initialises geometric and interpolation data for Ghost Points (GPs).
   *
   * This function iterates over all ghost points identified in the `computeMarkers` step.
   * For each ghost point, it performs the following operations:
   * 1. Identifies the specific geometry (body) the ghost point belongs to.
   * 2. Finds the closest point on the surface (IB point) using BVH trees.
   * 3. Computes the normal distance from the ghost point to the surface.
   * 4. Constructs a local orthonormal frame (normal, tangent1, tangent2) at the IB point.
   * 5. Projects "Image Points" (IMs) into the fluid domain along the surface normal.
   * 6. Computes trilinear interpolation weights and indices for these Image Points.
   *
   * All computed data is stored in the level-wide `GPStore` (CSR layout) for use in boundary condition reconstruction.
   *
   * \param lev The current AMR level index.
   */
  void initialiseGPs(int lev) {
    BL_PROFILE("IBM::initialiseGPs");
    const Real timing_total_start = beginPerformanceTiming();
    auto& mfab = *bmf_a[lev];
    auto& gpstore = gpstore_a[lev];
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

#ifdef AMREX_USE_GPU
    // ================================================================
    // GPU path: ParallelFor with atomic counter for GP index
    // ================================================================

    // ------------------------------------------------------------------
    // Pack per-geometry data into trivially-copyable POD GpuArrays so the
    // device lambda can capture-by-value.  Stack-allocated, fixed length
    // MAX_NGEOM (16) — multi-body limit; raise in ibm_containers.h if hit.
    // ------------------------------------------------------------------
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
        "ngeom exceeds MAX_NGEOM; increase MAX_NGEOM in ibm_containers.h");
    const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);

    // bqv[ii] : BVH4 query view for geometry ii (POD: raw ptrs to nodes/
    //           verts/faces). Used on device for closest-point queries in
    //           body frame. Owning BVH lives in bvh_a[ii].
    GpuArray<BVH4QueryView, MAX_NGEOM> bqv;

    // geom_off[ii] : level-wide flat-id offset for geometry ii's elements.
    //                Body-local prim_id + geom_off[ii] = global element id
    //                (used to index into LocalFrame_a / SurfElem_a).
    //                geom_off[ngeom] is the sentinel = total #elements.
    GpuArray<int, MAX_NGEOM + 1> geom_off;

    // transforms[ii] : rigid-body transform (body→world) for geometry ii.
    //                  Identity for static bodies; refreshed by FSI before
    //                  this kernel via updateRigidTransform().
    GpuArray<RigidTransform, MAX_NGEOM> transforms;

    for (int ii = 0; ii < ngeom_local; ii++) {
        bqv[ii]        = bvh_a[ii].query_view(geom_a[ii]);
        geom_off[ii]   = geom_offsets[ii];
        transforms[ii] = transform_a[ii];
    }
    geom_off[ngeom_local] = geom_offsets[ngeom_local];

    // ------------------------------------------------------------------
    // Level-scalar constants captured by value into the device lambda.
    // ------------------------------------------------------------------

    // Cell sizes on this level (dx, dy[, dz]).
    const auto dx_lev       = dx_a[lev];

    // Image-point spacing along the surface normal (= alpha * cell_diagonal).
    // search_*_image_point uses multiples of this to place IPs.
    const Real di_lev       = di_a[lev];

    // Cell diagonal length on this level — used as upper-bound sanity
    // check for ghost-point ↔ closest-surface-point distance.
    const Real diag_lev     = diag_a[lev];

    // Raw pointer to level-wide LocalFrame array (per-element body-frame
    // normal + tangents). Indexed by global f_idx = prim_id + geom_off[ii].
    const auto* lf_ptr      = LocalFrame_a.data();
    const auto* se_ptr      = SurfElem_a.data();

    // GPU-capturable POD view of the CSR GPStore. Holds both const read
    // pointers and writable _w pointers used during this init pass only.
    auto gpview = gpstore.view();

    // All host-visible counters batched into ONE device array read back with
    // a single copy+sync after the MFIter loop (previously: one blocking
    // dataValue() readback + stream sync PER FAB for the count verification,
    // plus a reset kernel per FAB):
    //   slots [0, nfabs)  — per-FAB atomic GP counters (assign GP slots);
    //                       pre-zeroed here, so no per-FAB reset is needed.
    //   slot  nfabs       — first-IP stencil-out-of-box failure count
    //   slot  nfabs+1     — first-IP below-threshold placement failure count
    //   slot  nfabs+2     — first-IP interpolation fit failure count
    //   slot  nfabs+3     — first-IP constant-fit fallback count
    //   slot  nfabs+4     — image-point WLS blocks shifted tangentially
    //   slot  nfabs+5     — fluid supports checked for line of sight
    //   slot  nfabs+6     — supports occluded by an IBM surface
    //   slot  nfabs+7     — first-IP fits demoted by visibility filtering
    //   slot  nfabs+8     — minimum representative first-hit element id
    //   slot  nfabs+9     — BI-constrained GP fit failures
    //   slot  nfabs+10    — BI-constrained GP linear fallbacks
    // (failure counters give GPU parity with the CPU path's abort/print
    // semantics: device printf is lossy and asserts are stripped in release,
    // so failures are aggregated and reported host-side).
    const int nfabs_local = mfab.local_size();

    // Host mirror of the managed CSR offsets. fab_offsets is managed memory,
    // and on devices with concurrentManagedAccess==0 (WSL) ANY host touch of
    // managed memory while device work is in flight faults. The old per-FAB
    // blocking readbacks serialized the device before every host read by
    // accident; with batched counters the launch loops below run ahead
    // asynchronously, so every host-side offset read must go through this
    // plain host copy instead (device kernels keep the managed copy via
    // gpview). The one sync here replaces this function's former 2-per-FAB.
    Gpu::streamSynchronize();
    Vector<int> h_fab_offsets(nfabs_local + 1);
    std::copy(gpstore.fab_offsets.begin(), gpstore.fab_offsets.end(),
              h_fab_offsets.begin());

    Vector<int> h_counts(nfabs_local + 11, 0);
    h_counts[nfabs_local + 8] = std::numeric_limits<int>::max();
    Gpu::DeviceVector<int> d_counts(nfabs_local + 11);
    Gpu::copyAsync(Gpu::hostToDevice, h_counts.begin(), h_counts.end(),
                   d_counts.begin());
    int* p_gp_count    = d_counts.data();                    // per-FAB slots
    int* p_err_stencil = d_counts.data() + nfabs_local;      // first-IP stencil out of box
    int* p_err_place   = d_counts.data() + nfabs_local + 1;  // first-IP below threshold
    int* p_err_fit     = d_counts.data() + nfabs_local + 2;  // first-IP LS/weight fit failure
    int* p_fit_constant = d_counts.data() + nfabs_local + 3; // explicit rank-deficient fallback
    int* p_tangent_shift = d_counts.data() + nfabs_local + 4;
    int* p_visibility_supports = d_counts.data() + nfabs_local + 5;
    int* p_visibility_occluded = d_counts.data() + nfabs_local + 6;
    int* p_visibility_demoted = d_counts.data() + nfabs_local + 7;
    int* p_visibility_first_hit = d_counts.data() + nfabs_local + 8;
    int* p_constrained_fit_fail = d_counts.data() + nfabs_local + 9;
    int* p_constrained_linear = d_counts.data() + nfabs_local + 10;

    int ifab_local = 0;
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      const int gp_offset = h_fab_offsets[ifab_local];
      const int ngps_fab  = h_fab_offsets[ifab_local + 1] - gp_offset;

      if (ngps_fab == 0) continue;

      // Image-point support is wider than the RHS stencil when eorder>1.
      // The bounds match the enlarged primitive/marker exchange halo, while
      // interpolationReadBox retains the physical FillPatch limit.
      const Box bxg = interpolationReadBox(
          mfi, lev, volumeInterpolationNghost(lev));
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);

      // This FAB's dedicated counter slot (pre-zeroed at allocation above,
      // so the old per-FAB single-thread reset kernel is gone).
      int* const p_fab_count = p_gp_count + ifab_local;

      const Box bxgp = amrex::grow(bx, GP_BOX_EXTRA);
      amrex::ParallelFor(bxgp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        // NVCC extended lambdas cannot first-capture a variable from inside an
        // if-constexpr branch.  Anchor the optional visibility inputs here;
        // mode 0 still compiles them away with the rest of the audit path.
        amrex::ignore_unused(
            bqv, transforms, geom_off, ngeom_local, diag_lev,
            p_visibility_supports, p_visibility_occluded,
            p_visibility_demoted, p_visibility_first_hit,
            p_constrained_fit_fail, p_constrained_linear);
        if (!ibMarkers(i, j, k, 1)) return;

        // Atomically claim a slot in the flat GP array.
        // Note: local_idx assignment order is non-deterministic on GPU —
        // the GP at the same physical cell (i,j,k) may receive different
        // gidx across runs. This does not affect physics: downstream reads
        // gpview.gp_ijk[ii] to recover (i,j,k), and final ghost-cell
        // primitives are written back to those cells regardless of how
        // gidx is permuted.
        int local_idx = Gpu::Atomic::Add(p_fab_count, 1);
        const int gidx = gp_offset + local_idx;

        Point gp = make_grid_point(prob_lo, dx_lev, i, j, k);
        gpview.gp_ijk_w[gidx] = make_vec<int>(i, j, k);  // dim-aware: 2D ignores k

        int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;
        AMREX_ASSERT(geomIdx >= 0 && geomIdx < ngeom_local);
        gpview.geomIdx_w[gidx] = geomIdx;

        // Transform query point to body frame, then BVH closest-point query
        const auto& T = transforms[geomIdx];
        Point gp_body = T.to_body(gp);
        ClosestPointResult closest_elem = bqv[geomIdx].closest_point_query(gp_body);
        int f_idx = closest_elem.prim_id + geom_off[geomIdx];
        gpview.elemIdx_w[gidx] = f_idx;

        // Transform the faceted closest point back to world coordinates. The
        // owning element remains the BVH result, while an enabled smooth-wall
        // contract supplies the unique closure boundary intercept.
        const Point raw_cp = T.to_world(closest_elem.point);
        const auto raw_cp_vec = make_vec<Real>(raw_cp);
        Array1D<Real, 0, AMREX_SPACEDIM - 1> cp_vec;
        apply_volume_surface_point_override(raw_cp_vec, cp_vec);
        const Point cp = surface_point_from_array(cp_vec);

        // Distance (invariant under rigid transform, but computed in world frame)
        Real disGP_val = std::sqrt(point_distance_sq(gp, cp));

        // Sanity: a layer-k ghost point is at most k Cartesian cells from a
        // fluid cell.  WENO5/TENO5 may use three reconstructed GP layers to
        // retain a complete crossing stencil, so scale this bound with the
        // configured GP band rather than assuming ghost_layers=1.
        // Failure here means one of:
        //   (a) BVH returned a wrong primitive (geometry data corruption)
        //   (b) marker comp 1 was set on a cell not adjacent to fluid
        //   (c) FAB layout / DistributionMapping inconsistency between
        //       computeMarkers and initialiseGPs (shouldn't happen, but
        //       would manifest here first if a regrid slipped between).
        AMREX_ASSERT(disGP_val < Real(ghost_layers) * diag_lev);
        gpview.disGP_w[gidx] = disGP_val;
        gpview.ib_xyz_w[gidx] = cp_vec;

        // Rotate body-frame LocalFrame to world frame
        const LocalFrame& lf_body = lf_ptr[f_idx];
        LocalFrame localframe;
        T.rotate_to_world(lf_body.normal,   localframe.normal);
        T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
        T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif
        apply_volume_surface_frame_override(cp_vec, localframe);

        // Image points: stack-local SoA, written into GPStore once at the end.
        // Value-initialisation keeps diagnostics deterministic before the
        // aggregated host-side failure check aborts an invalid first IP.
        Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz{};
        Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk{};
        Array1D<Real, 0, eorder_tparm - 1> disIM{};
        Array1D< int, 0, eorder_tparm - 1> imp_ninterp{};
        Array1D< int, 0, eorder_tparm - 1> imp_fit_order{};

        // Chained walk along the outward normal — jj=0 from IB point,
        // jj>=1 from previous IP. Shared with CPU path.
        const int place_status = place_image_points<eorder_tparm, iorder_tparm>(
            cp, localframe,
            lev, prob_lo, dx_lev, di_lev,
            bxg, ibMarkers,
            gpview, gidx,
            imp_xyz, imp_ijk, disIM, imp_ninterp);
        if (place_status & 1) { Gpu::Atomic::Add(p_err_stencil, 1); }
        if (place_status & 2) { Gpu::Atomic::Add(p_err_place, 1); }

        if (sharp_feature_tangential_reconstruction) {
          const SurfElem& edge = se_ptr[f_idx];
          Real coordinate = Real(0.5) * edge.measure;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real centroid_world = T.d[d];
            for (int q = 0; q < AMREX_SPACEDIM; ++q) {
              centroid_world += T.R[d][q] * edge.centroid[q];
            }
            coordinate +=
                (gpview.ib_xyz_w[gidx](d) - centroid_world) *
                sourceDirectedTangent(localframe, d);
          }
          const Real h = amrex::min(dx_lev[0], dx_lev[1]);
          const int shifted =
              applySharpFeatureTangentialStencil<eorder_tparm,
                                                  iorder_tparm>(
                  edge, coordinate, h, localframe, lev, prob_lo, dx_lev,
                  bxg, ibMarkers, gpview, gidx, imp_xyz, imp_ijk,
                  imp_ninterp);
          if (shifted > 0) Gpu::Atomic::Add(p_tangent_shift, shifted);
        }
        gpview.imp_xyz_w[gidx]       = imp_xyz;
        gpview.imp_ijk_w[gidx]       = imp_ijk;
        gpview.disIM_w[gidx]         = disIM;

        Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1> imp_ipweights;
        Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;
        Array2D<uint8_t, 0, eorder_tparm - 1, 0, N_InterP - 1>
            visibility_mask{};
        Array2D<uint8_t, 0, eorder_tparm - 1, 0, N_InterP - 1>
            weight_mask;
        Array2D<int, 0, eorder_tparm - 1, 0, N_InterP - 1>
            first_hit_element{};
        Array1D<int, 0, eorder_tparm - 1> candidate_count{};
        Array1D<int, 0, eorder_tparm - 1> visible_count{};
        for (int image = 0; image < eorder_tparm; ++image) {
          for (int support = 0; support < N_InterP; ++support) {
            weight_mask(image, support) = uint8_t(1);
          }
        }

        if constexpr (support_visibility_mode > 0) {
          int visible_supports = 0;
          int occluded_supports = 0;
          int first_hit = -1;
          buildSupportVisibilityMask<eorder_tparm, iorder_tparm>(
              visibility_mask, first_hit_element, candidate_count,
              visible_count, imp_xyz, imp_ijk, imp_ninterp, prob_lo, dx_lev,
              ibMarkers, bqv, transforms, geom_off, ngeom_local,
              visible_supports, occluded_supports, first_hit);
          Gpu::Atomic::Add(p_visibility_supports, visible_supports);
          Gpu::Atomic::Add(p_visibility_occluded, occluded_supports);
          if (first_hit >= 0) {
            Gpu::Atomic::Min(p_visibility_first_hit, first_hit);
          }
          for (int image = 0; image < eorder_tparm; ++image) {

            for (int support = 0; support < N_InterP; ++support) {
              if constexpr (support_visibility_mode == 2) {
                weight_mask(image, support) =
                    visibility_mask(image, support);
              }
            }
            if constexpr (support_visibility_mode == 2) {
              imp_ninterp(image) = visible_count(image);
            }
          }
        }

        computeIPweights<eorder_tparm, iorder_tparm, GPSTOREVIEW>(
            imp_ipweights, imp_ip_ijk,
            imp_xyz, imp_ijk, imp_ninterp, imp_fit_order,
            weight_mask, prob_lo, dx_lev, ibMarkers,
            support_visibility_mode == 2
                ? support_visibility_condition_max
                : std::numeric_limits<Real>::infinity());
        if (place_status == 0 && imp_ninterp(0) == 0) {
          Gpu::Atomic::Add(p_err_fit, 1);
        }
        if (place_status == 0 && imp_fit_order(0) == 0) {
          Gpu::Atomic::Add(p_fit_constant, 1);
        }

        // Store imp_ninterp AFTER computeIPweights: the iorder=2 WLS path
        // demotes an unfittable image point by setting imp_ninterp(iim)=0,
        // and n_valid/eff_order downstream must see that demotion.
        gpview.imp_ninterp_w[gidx]   = imp_ninterp;
        gpview.imp_ipweights_w[gidx] = imp_ipweights;
        gpview.imp_ip_ijk_w[gidx]    = imp_ip_ijk;

        if constexpr (use_bi_constrained_shared_gp) {
          const auto gp_vec = make_vec<Real>(gp);
          Array2D<int, 0, N_InterP - 1, 0, AMREX_SPACEDIM - 1>
              support_indices{};
          Array1D<uint8_t, 0, N_InterP - 1> support_mask{};
          Array1D<Real, 0, N_InterP - 1> dirichlet_weights{};
          Array1D<Real, 0, N_InterP - 1> neumann_weights{};
          Array1D<Real, 0, N_InterP - 1> fluid_trace_weights{};
          Array1D<Real, 0, N_InterP - 1> fv_dirichlet_weights{};
          Array1D<Real, 0, N_InterP - 1> entropy_jet_weights{};
          Array1D<Real, 0, N_InterP - 1> entropy_jet_linear_weights{};
          RZBICPointFunctionals<N_InterP> rz_point_functionals{};
          Real dirichlet_boundary_weight = Real(0.0);
          Real neumann_gradient_weight = Real(0.0);
          Real fv_dirichlet_boundary_weight = Real(0.0);
          int fluid_trace_order = -1;
          Real fluid_trace_condition =
              std::numeric_limits<Real>::infinity();
          Real fluid_trace_weight_l1 =
              std::numeric_limits<Real>::infinity();
          int constrained_support_count = 0;
          Real constrained_condition =
              std::numeric_limits<Real>::infinity();
          Real constrained_weight_l1 =
              std::numeric_limits<Real>::infinity();
          Real fv_dirichlet_condition =
              std::numeric_limits<Real>::infinity();
          Real fv_dirichlet_weight_l1 =
              std::numeric_limits<Real>::infinity();
          Real entropy_jet_condition =
              std::numeric_limits<Real>::infinity();
          Real entropy_jet_weight_l1 =
              std::numeric_limits<Real>::infinity();
          Real entropy_jet_linear_condition =
              std::numeric_limits<Real>::infinity();
          Real entropy_jet_linear_weight_l1 =
              std::numeric_limits<Real>::infinity();
          int fv_dirichlet_order = -1;
          int entropy_jet_order = -1;
          const int constrained_order =
              buildBIConstrainedGPFunctionals<GPSTOREVIEW, N_InterP>(
                  cp_vec, gp_vec, localframe, lev, prob_lo, dx_lev, diag_lev,
                  bxg, ibMarkers, gpview, gidx, bqv, transforms, geom_off,
                  ngeom_local, support_indices, support_mask,
                  dirichlet_weights,
                  dirichlet_boundary_weight, neumann_weights,
                  neumann_gradient_weight, fluid_trace_weights,
                  fluid_trace_order, fluid_trace_condition,
                  fluid_trace_weight_l1, constrained_support_count,
                  constrained_condition, constrained_weight_l1,
                  fv_dirichlet_weights, fv_dirichlet_boundary_weight,
                  fv_dirichlet_order, fv_dirichlet_condition,
                  fv_dirichlet_weight_l1, entropy_jet_weights,
                  entropy_jet_linear_weights, entropy_jet_order,
                  entropy_jet_condition, entropy_jet_weight_l1,
                  entropy_jet_linear_condition,
                  entropy_jet_linear_weight_l1, rz_point_functionals);

          gpview.constrained_support_ijk_w[gidx] = support_indices;
          gpview.constrained_dirichlet_weights_w[gidx] =
              dirichlet_weights;
          gpview.constrained_neumann_weights_w[gidx] = neumann_weights;
          gpview.constrained_fv_dirichlet_weights_w[gidx] =
              fv_dirichlet_weights;
          gpview.constrained_entropy_jet_weights_w[gidx] =
              entropy_jet_weights;
          gpview.constrained_entropy_jet_linear_weights_w[gidx] =
              entropy_jet_linear_weights;
          gpview.constrained_fluid_trace_weights_w[gidx] =
              fluid_trace_weights;
          gpview.constrained_dirichlet_boundary_weight_w[gidx] =
              dirichlet_boundary_weight;
          gpview.constrained_neumann_gradient_weight_w[gidx] =
              neumann_gradient_weight;
          gpview.constrained_fv_dirichlet_boundary_weight_w[gidx] =
              fv_dirichlet_boundary_weight;
          gpview.constrained_order_w[gidx] = constrained_order;
          gpview.constrained_fv_dirichlet_order_w[gidx] =
              fv_dirichlet_order;
          gpview.constrained_entropy_jet_order_w[gidx] = entropy_jet_order;
          gpview.constrained_fluid_trace_order_w[gidx] = fluid_trace_order;
          gpview.constrained_support_count_w[gidx] =
              constrained_support_count;
          gpview.constrained_condition_w[gidx] = constrained_condition;
          gpview.constrained_weight_l1_w[gidx] = constrained_weight_l1;
          gpview.constrained_fv_dirichlet_condition_w[gidx] =
              fv_dirichlet_condition;
          gpview.constrained_fv_dirichlet_weight_l1_w[gidx] =
              fv_dirichlet_weight_l1;
          gpview.constrained_entropy_jet_condition_w[gidx] =
              entropy_jet_condition;
          gpview.constrained_entropy_jet_weight_l1_w[gidx] =
              entropy_jet_weight_l1;
          gpview.constrained_entropy_jet_linear_condition_w[gidx] =
              entropy_jet_linear_condition;
          gpview.constrained_entropy_jet_linear_weight_l1_w[gidx] =
              entropy_jet_linear_weight_l1;
          gpview.constrained_fluid_trace_condition_w[gidx] =
              fluid_trace_condition;
          gpview.constrained_fluid_trace_weight_l1_w[gidx] =
              fluid_trace_weight_l1;
          if constexpr (bic_point_target_functionals) {
            gpview.rz_bic_point_functionals_w[gidx] =
                rz_point_functionals;
          }
          if (constrained_order < 1) {
            Gpu::Atomic::Add(p_constrained_fit_fail, 1);
          } else if (constrained_order == 1) {
            Gpu::Atomic::Add(p_constrained_linear, 1);
          }
        }

        ibMarkers(i, j, k, 1) = static_cast<uint8_t>(imp_ninterp(0));
      }); // end ParallelFor

    } // end MFIter

    // Single batched readback of per-FAB counts + error counters (one sync
    // for the whole level). Verifying after all FABs launched instead of
    // per-FAB does not add risk: any CSR overflow from a miscounted FAB
    // already happened inside that FAB's own kernel (the atomic claims
    // slots during the launch, before any readback could run), and in both
    // versions we abort before the GP data is consumed.
    Gpu::copy(Gpu::deviceToHost, d_counts.begin(), d_counts.end(),
              h_counts.begin());
    for (int f = 0; f < nfabs_local; ++f) {
      const int gp_offset = h_fab_offsets[f];
      const int ngps_fab  = h_fab_offsets[f + 1] - gp_offset;
      // Fast-skipped FABs (ngps_fab==0) never touch their pre-zeroed slot.
      if (h_counts[f] != ngps_fab) {
        amrex::Print() << "Error in initialiseGPs (GPU): GP count mismatch\n"
                       << "  level=" << lev
                       << "  FAB=" << f
                       << "  expected=" << ngps_fab
                       << "  got=" << h_counts[f]
                       << "  gp_offset=" << gp_offset << "\n";
        amrex::Abort("initialiseGPs: ghost point count mismatch");
      }
    }

    // Aggregated failure report (parity with the CPU path, which aborts on a
    // first-IP stencil escape and prints per-GP placement failures).
    // Counters were read back in the batched copy above — no extra sync.
    {
      const int n_stencil = h_counts[nfabs_local];
      const int n_place   = h_counts[nfabs_local + 1];
      const int n_fit     = h_counts[nfabs_local + 2];
      const int n_constant = h_counts[nfabs_local + 3];
      const int n_tangent_shift = h_counts[nfabs_local + 4];
      if (n_stencil > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_stencil
                       << " ghost point(s) on level " << lev
                       << " have their first image-point stencil outside the grown box.\n";
        amrex::Abort("initialiseGPs: interpolation stencil out of box (GPU)");
      }
      if (n_place > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_place
                       << " ghost point(s) on level " << lev
                       << " failed first image-point placement (below threshold).\n";
        amrex::Abort("initialiseGPs: invalid first image point (GPU)");
      }
      if (n_fit > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_fit
                       << " ghost point(s) on level " << lev
                       << " have no resolvable interpolation fit at the first image point.\n";
        amrex::Abort("initialiseGPs: invalid first image-point fit (GPU)");
      }
      if (n_constant > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_constant
                       << " ghost point(s) on level " << lev
                       << " use the explicit constant-fit interpolation fallback.\n";
      }
      if constexpr (sharp_feature_tangential_reconstruction) {
        amrex::Print() << "initialiseGPs (GPU): " << n_tangent_shift
                       << " image-point WLS block(s) on level " << lev
                       << " use same-edge-biased tangential reconstruction.\n";
      }
      if constexpr (support_visibility_mode > 0) {
        const int checked = h_counts[nfabs_local + 5];
        const int occluded = h_counts[nfabs_local + 6];
        const int demoted = h_counts[nfabs_local + 7];
        const int first_hit = h_counts[nfabs_local + 8];
        const Real fraction = checked > 0
            ? Real(100.0) * Real(occluded) / Real(checked)
            : Real(0.0);
        amrex::Print()
            << "[IBM-Visibility] GP level " << lev
            << " mode=" << support_visibility_mode
            << " checked=" << checked
            << " occluded=" << occluded
            << " (" << fraction << "%)";
        amrex::ignore_unused(demoted);
        if (first_hit != std::numeric_limits<int>::max()) {
          amrex::Print() << " representative-element=" << first_hit;
        }
        amrex::Print() << "\n";
      }
      if constexpr (use_bi_constrained_shared_gp) {
        const int fit_failures = h_counts[nfabs_local + 9];
        const int linear_fallbacks = h_counts[nfabs_local + 10];
        if (fit_failures > 0) {
          amrex::Print()
              << "initialiseGPs (GPU): " << fit_failures
              << " BI-constrained shared-GP fit(s) failed on level " << lev
              << ".\n";
          amrex::Abort(
              "initialiseGPs: invalid BI-constrained shared-GP fit (GPU)");
        }
        amrex::Print()
            << "[IBM-BI-Constrained-GP] level=" << lev
            << " total=" << gpstore.total_ngps
            << " quadratic=" << gpstore.total_ngps - linear_fallbacks
            << " linear=" << linear_fallbacks << "\n";
      }
    }

#else
    // ================================================================
    // CPU path: original LoopOnCpu implementation
    // ================================================================
#ifndef AMREX_USE_CGAL
    GpuArray<BVH4QueryView, MAX_NGEOM> visibility_queries{};
    GpuArray<RigidTransform, MAX_NGEOM> visibility_transforms{};
    GpuArray<int, MAX_NGEOM + 1> visibility_geometry_offsets{};
    if constexpr (support_visibility_mode > 0) {
      // geom_offsets is managed storage populated during geometry setup. On
      // WSL2, host access before the asynchronous initialization completes can
      // fault even though ordinary unified-memory systems tolerate it.
      Gpu::streamSynchronize();
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
          "ngeom exceeds MAX_NGEOM; increase MAX_NGEOM in ibm_containers.h");
      for (int geometry = 0; geometry < ngeom; ++geometry) {
        visibility_queries[geometry] =
            bvh_a[geometry].query_view(geom_a[geometry]);
        visibility_transforms[geometry] = transform_a[geometry];
        visibility_geometry_offsets[geometry] = geom_offsets[geometry];
      }
      visibility_geometry_offsets[ngeom] = geom_offsets[ngeom];
    }
#endif
    int ifab_local = 0;
    int n_stencil_fail = 0;
    int n_place_fail = 0;
    int n_fit_fail = 0;
    int n_constant_fit = 0;
    int n_tangent_shift = 0;
    int n_visibility_supports = 0;
    int n_visibility_occluded = 0;
    int n_visibility_demoted = 0;
    int visibility_first_hit = std::numeric_limits<int>::max();
    int n_constrained_fit_fail = 0;
    int n_constrained_linear = 0;
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      // Get the pre-allocated range in the flattened GPStore for this FAB.
      // (CPU build: ManagedVector is plain host memory — direct read is safe;
      //  the h_fab_offsets host mirror exists only in the GPU branch above.)
      const int gp_offset = gpstore.fab_offsets[ifab_local];
      const int ngps_fab  = gpstore.fab_offsets[ifab_local + 1] - gp_offset;

      if (ngps_fab == 0) continue;

      const Box bxg = interpolationReadBox(
          mfi, lev, volumeInterpolationNghost(lev));
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);  // uint8_t array

      // CPU loop: Write directly into gpstore at pre-allocated offsets (no push_back).
      int gp_count = 0; 
      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA), [&](int i, int j, int k) {
        // for each ghost point
        if (ibMarkers(i, j, k, 1)) {
          const int gidx = gp_offset + gp_count;  // global flat index in GPStore
          gp_count++;

          Point gp = make_grid_point(prob_lo, dx_a[lev], i, j, k);
          gpstore.gp_ijk[gidx] = make_vec<int>(i, j, k);
          
          int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;
          AMREX_ASSERT_WITH_MESSAGE(geomIdx >= 0 && geomIdx < ngeom, 
                  "Invalid geometry index in initialiseGPs");
          gpstore.geomIdx[gidx] = geomIdx;

          // Rigid-body transform: geometry trees/frames live in BODY frame;
          // query in body frame, return results to world frame (mirrors the
          // GPU path — previously missing here, which made CPU builds use
          // wrong closest points/normals once an FSI body moved).
          const auto& T = transform_a[geomIdx];
          Point gp_body = T.to_body(gp);

#ifdef AMREX_USE_CGAL
          ClosestPointResult closest_elem =
              cgal_closest_point_query(tree_a[geomIdx], idxmap_a[geomIdx],
                                       gp_body, this->geom_offsets[geomIdx]);
#else
          ClosestPointResult closest_elem =
              bvh_a[geomIdx].closest_point_query(gp_body, geom_a[geomIdx]);
#endif

          int f_idx = closest_elem.prim_id + this->geom_offsets[geomIdx];
          gpstore.elemIdx[gidx] = f_idx;

          const Point raw_cp = T.to_world(closest_elem.point);
          const auto raw_cp_vec = make_vec<Real>(raw_cp);
          Array1D<Real, 0, AMREX_SPACEDIM - 1> cp_vec;
          apply_volume_surface_point_override(raw_cp_vec, cp_vec);
          const Point cp = surface_point_from_array(cp_vec);

          Real disGP = std::sqrt(point_distance_sq(gp, cp));

          AMREX_ASSERT_WITH_MESSAGE(
              disGP < Real(ghost_layers) * diag_a[lev],
              "Ghost point and boundary-intercept distance exceeds the "
              "configured ghost-layer diagonal bound");
          gpstore.disGP[gidx] = disGP;

          gpstore.ib_xyz[gidx] = cp_vec;

          const LocalFrame& lf_body = LocalFrame_a[f_idx];
          LocalFrame localframe;
          T.rotate_to_world(lf_body.normal,   localframe.normal);
          T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
          T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif
          apply_volume_surface_frame_override(cp_vec, localframe);

          // Value-initialised ({}) so a failed search remains deterministic
          // until the aggregated host-side failure check below aborts.
          Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz{};
          Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk{};
          Array1D<Real, 0, eorder_tparm - 1> disIM{};
          Array1D< int, 0, eorder_tparm - 1> imp_ninterp{};
          Array1D< int, 0, eorder_tparm - 1> imp_fit_order{};

          // Chained walk along the outward normal — shared with GPU path.
          const int place_status = place_image_points<eorder_tparm, iorder_tparm>(
              cp, localframe,
              lev, prob_lo, dx_a[lev], di_a[lev],
              bxg, ibMarkers,
              gpstore, gidx,
              imp_xyz, imp_ijk, disIM, imp_ninterp);
          if (place_status & 1) { ++n_stencil_fail; }
          if (place_status & 2) { ++n_place_fail; }

          if constexpr (sharp_feature_tangential_reconstruction) {
            const SurfElem& edge = SurfElem_a[f_idx];
            Real coordinate = Real(0.5) * edge.measure;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              Real centroid_world = T.d[d];
              for (int q = 0; q < AMREX_SPACEDIM; ++q) {
                centroid_world += T.R[d][q] * edge.centroid[q];
              }
              coordinate +=
                  (gpstore.ib_xyz[gidx](d) - centroid_world) *
                  sourceDirectedTangent(localframe, d);
            }
            const Real h = amrex::min(dx_a[lev][0], dx_a[lev][1]);
            n_tangent_shift +=
                applySharpFeatureTangentialStencil<eorder_tparm,
                                                    iorder_tparm>(
                    edge, coordinate, h, localframe, lev, prob_lo,
                    dx_a[lev], bxg, ibMarkers, gpstore, gidx, imp_xyz,
                    imp_ijk, imp_ninterp);
          }
          gpstore.imp_xyz[gidx] = imp_xyz;
          gpstore.imp_ijk[gidx] = imp_ijk;
          gpstore.disIM[gidx] = disIM;

          Array2D<Real, 0, eorder_tparm - 1 , 0, N_InterP -1 > imp_ipweights;
          Array3D< int, 0, eorder_tparm - 1 , 0, N_InterP -1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;
          Array2D<uint8_t, 0, eorder_tparm - 1, 0, N_InterP - 1>
              visibility_mask{};
          Array2D<uint8_t, 0, eorder_tparm - 1, 0, N_InterP - 1>
              weight_mask;
          Array2D<int, 0, eorder_tparm - 1, 0, N_InterP - 1>
              first_hit_element{};
          Array1D<int, 0, eorder_tparm - 1> candidate_count{};
          Array1D<int, 0, eorder_tparm - 1> visible_count{};
          for (int image = 0; image < eorder_tparm; ++image) {
            for (int support = 0; support < N_InterP; ++support) {
              weight_mask(image, support) = uint8_t(1);
            }
          }

#ifndef AMREX_USE_CGAL
          if constexpr (support_visibility_mode > 0) {
            int checked = 0;
            int occluded = 0;
            int first_hit = -1;
            buildSupportVisibilityMask<eorder_tparm, iorder_tparm>(
                visibility_mask, first_hit_element, candidate_count,
                visible_count, imp_xyz, imp_ijk, imp_ninterp, prob_lo,
                dx_a[lev], ibMarkers, visibility_queries,
                visibility_transforms, visibility_geometry_offsets, ngeom,
                checked, occluded, first_hit);
            n_visibility_supports += checked;
            n_visibility_occluded += occluded;
            if (first_hit >= 0) {
              visibility_first_hit =
                  amrex::min(visibility_first_hit, first_hit);
            }
            for (int image = 0; image < eorder_tparm; ++image) {

              for (int support = 0; support < N_InterP; ++support) {
                if constexpr (support_visibility_mode == 2) {
                  weight_mask(image, support) =
                      visibility_mask(image, support);
                }
              }
              if constexpr (support_visibility_mode == 2) {
                imp_ninterp(image) = visible_count(image);
              }
            }
          }
#endif

          computeIPweights<eorder_tparm, iorder_tparm, GPSTORE>(
              imp_ipweights,
              imp_ip_ijk,
              imp_xyz,
              imp_ijk,
              imp_ninterp,
              imp_fit_order,
              weight_mask, prob_lo, dx_a[lev], ibMarkers,
              support_visibility_mode == 2
                  ? support_visibility_condition_max
                  : std::numeric_limits<Real>::infinity());
          if (place_status == 0 && imp_ninterp(0) == 0) {
            ++n_fit_fail;
          }
          if (place_status == 0 && imp_fit_order(0) == 0) {
            ++n_constant_fit;
          }

          // Store imp_ninterp AFTER computeIPweights (WLS demotion visibility).
          gpstore.imp_ninterp[gidx] = imp_ninterp;
          gpstore.imp_ipweights[gidx] = imp_ipweights;
          gpstore.imp_ip_ijk[gidx] = imp_ip_ijk;

#ifndef AMREX_USE_CGAL
          if constexpr (use_bi_constrained_shared_gp) {
            const auto gp_vec = make_vec<Real>(gp);
            Array2D<int, 0, N_InterP - 1, 0, AMREX_SPACEDIM - 1>
                support_indices{};
            Array1D<uint8_t, 0, N_InterP - 1> support_mask{};
            Array1D<Real, 0, N_InterP - 1> dirichlet_weights{};
            Array1D<Real, 0, N_InterP - 1> neumann_weights{};
            Array1D<Real, 0, N_InterP - 1> fluid_trace_weights{};
            Array1D<Real, 0, N_InterP - 1> fv_dirichlet_weights{};
            Array1D<Real, 0, N_InterP - 1> entropy_jet_weights{};
            Array1D<Real, 0, N_InterP - 1> entropy_jet_linear_weights{};
            RZBICPointFunctionals<N_InterP> rz_point_functionals{};
            Real dirichlet_boundary_weight = Real(0.0);
            Real neumann_gradient_weight = Real(0.0);
            Real fv_dirichlet_boundary_weight = Real(0.0);
            int fluid_trace_order = -1;
            Real fluid_trace_condition =
                std::numeric_limits<Real>::infinity();
            Real fluid_trace_weight_l1 =
                std::numeric_limits<Real>::infinity();
            int constrained_support_count = 0;
            Real constrained_condition =
                std::numeric_limits<Real>::infinity();
            Real constrained_weight_l1 =
                std::numeric_limits<Real>::infinity();
            Real fv_dirichlet_condition =
                std::numeric_limits<Real>::infinity();
            Real fv_dirichlet_weight_l1 =
                std::numeric_limits<Real>::infinity();
            Real entropy_jet_condition =
                std::numeric_limits<Real>::infinity();
            Real entropy_jet_weight_l1 =
                std::numeric_limits<Real>::infinity();
            Real entropy_jet_linear_condition =
                std::numeric_limits<Real>::infinity();
            Real entropy_jet_linear_weight_l1 =
                std::numeric_limits<Real>::infinity();
            int fv_dirichlet_order = -1;
            int entropy_jet_order = -1;
            const int constrained_order =
                buildBIConstrainedGPFunctionals<GPSTORE, N_InterP>(
                    cp_vec, gp_vec, localframe, lev, prob_lo, dx_a[lev],
                    diag_a[lev], bxg, ibMarkers, gpstore, gidx,
                    visibility_queries, visibility_transforms,
                    visibility_geometry_offsets, ngeom, support_indices,
                    support_mask, dirichlet_weights,
                    dirichlet_boundary_weight,
                    neumann_weights, neumann_gradient_weight,
                    fluid_trace_weights, fluid_trace_order,
                    fluid_trace_condition, fluid_trace_weight_l1,
                    constrained_support_count, constrained_condition,
                    constrained_weight_l1, fv_dirichlet_weights,
                    fv_dirichlet_boundary_weight, fv_dirichlet_order,
                    fv_dirichlet_condition, fv_dirichlet_weight_l1,
                    entropy_jet_weights, entropy_jet_linear_weights,
                    entropy_jet_order, entropy_jet_condition,
                    entropy_jet_weight_l1, entropy_jet_linear_condition,
                    entropy_jet_linear_weight_l1, rz_point_functionals);

            gpstore.constrained_support_ijk[gidx] = support_indices;
            gpstore.constrained_dirichlet_weights[gidx] =
                dirichlet_weights;
            gpstore.constrained_neumann_weights[gidx] = neumann_weights;
            gpstore.constrained_fv_dirichlet_weights[gidx] =
                fv_dirichlet_weights;
            gpstore.constrained_entropy_jet_weights[gidx] =
                entropy_jet_weights;
            gpstore.constrained_entropy_jet_linear_weights[gidx] =
                entropy_jet_linear_weights;
            gpstore.constrained_fluid_trace_weights[gidx] =
                fluid_trace_weights;
            gpstore.constrained_dirichlet_boundary_weight[gidx] =
                dirichlet_boundary_weight;
            gpstore.constrained_neumann_gradient_weight[gidx] =
                neumann_gradient_weight;
            gpstore.constrained_fv_dirichlet_boundary_weight[gidx] =
                fv_dirichlet_boundary_weight;
            gpstore.constrained_order[gidx] = constrained_order;
            gpstore.constrained_fv_dirichlet_order[gidx] =
                fv_dirichlet_order;
            gpstore.constrained_entropy_jet_order[gidx] = entropy_jet_order;
            gpstore.constrained_fluid_trace_order[gidx] = fluid_trace_order;
            gpstore.constrained_support_count[gidx] =
                constrained_support_count;
            gpstore.constrained_condition[gidx] = constrained_condition;
            gpstore.constrained_weight_l1[gidx] = constrained_weight_l1;
            gpstore.constrained_fv_dirichlet_condition[gidx] =
                fv_dirichlet_condition;
            gpstore.constrained_fv_dirichlet_weight_l1[gidx] =
                fv_dirichlet_weight_l1;
            gpstore.constrained_entropy_jet_condition[gidx] =
                entropy_jet_condition;
            gpstore.constrained_entropy_jet_weight_l1[gidx] =
                entropy_jet_weight_l1;
            gpstore.constrained_entropy_jet_linear_condition[gidx] =
                entropy_jet_linear_condition;
            gpstore.constrained_entropy_jet_linear_weight_l1[gidx] =
                entropy_jet_linear_weight_l1;
            gpstore.constrained_fluid_trace_condition[gidx] =
                fluid_trace_condition;
            gpstore.constrained_fluid_trace_weight_l1[gidx] =
                fluid_trace_weight_l1;
            if constexpr (bic_point_target_functionals) {
              gpstore.rz_bic_point_functionals[gidx] =
                  rz_point_functionals;
            }
            if (constrained_order < 1) {
              ++n_constrained_fit_fail;
            } else if (constrained_order == 1) {
              ++n_constrained_linear;
            }
          }
#endif

          ibMarkers(i, j, k, 1) = static_cast<uint8_t>(imp_ninterp(0));

          } //end if (ibMarkers(i,j,k,1))
      });//end loop on bx

      if(gp_count != ngps_fab) {
        amrex::Abort("Error in initialiseGPs: mismatch in ghost point count");
      }
    } //end MFIter

    if (n_stencil_fail > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_stencil_fail
                     << " ghost point(s) on level " << lev
                     << " have their first image-point stencil outside the grown box.\n";
      amrex::Abort("initialiseGPs: interpolation stencil out of box (CPU)");
    }
    if (n_place_fail > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_place_fail
                     << " ghost point(s) on level " << lev
                     << " failed first image-point placement (below threshold).\n";
      amrex::Abort("initialiseGPs: invalid first image point (CPU)");
    }
    if (n_fit_fail > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_fit_fail
                     << " ghost point(s) on level " << lev
                     << " have no resolvable interpolation fit at the first image point.\n";
      amrex::Abort("initialiseGPs: invalid first image-point fit (CPU)");
    }
    if (n_constant_fit > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_constant_fit
                     << " ghost point(s) on level " << lev
                     << " use the explicit constant-fit interpolation fallback.\n";
    }
    if constexpr (use_bi_constrained_shared_gp) {
      int constrained_total = gpstore.total_ngps;
      ParallelDescriptor::ReduceIntSum(constrained_total);
      ParallelDescriptor::ReduceIntSum(n_constrained_fit_fail);
      ParallelDescriptor::ReduceIntSum(n_constrained_linear);
      if (n_constrained_fit_fail > 0) {
        amrex::Print()
            << "initialiseGPs (CPU): " << n_constrained_fit_fail
            << " BI-constrained shared-GP fit(s) failed on level " << lev
            << ".\n";
        amrex::Abort(
            "initialiseGPs: invalid BI-constrained shared-GP fit (CPU)");
      }
      amrex::Print()
          << "[IBM-BI-Constrained-GP] level=" << lev
          << " total=" << constrained_total
          << " quadratic=" << constrained_total - n_constrained_linear
          << " linear=" << n_constrained_linear << "\n";
    }
    if constexpr (sharp_feature_tangential_reconstruction) {
      amrex::Print() << "initialiseGPs (CPU): " << n_tangent_shift
                     << " image-point WLS block(s) on level " << lev
                     << " use same-edge-biased tangential reconstruction.\n";
    }
    if constexpr (support_visibility_mode > 0) {
      ParallelDescriptor::ReduceIntSum(n_visibility_supports);
      ParallelDescriptor::ReduceIntSum(n_visibility_occluded);
      ParallelDescriptor::ReduceIntSum(n_visibility_demoted);
      ParallelDescriptor::ReduceIntMin(visibility_first_hit);
      const Real fraction = n_visibility_supports > 0
          ? Real(100.0) * Real(n_visibility_occluded) /
                Real(n_visibility_supports)
          : Real(0.0);
      amrex::Print()
          << "[IBM-Visibility] GP level " << lev
          << " mode=" << support_visibility_mode
          << " checked=" << n_visibility_supports
          << " occluded=" << n_visibility_occluded
          << " (" << fraction << "%)";
      if (visibility_first_hit != std::numeric_limits<int>::max()) {
        amrex::Print() << " representative-element=" << visibility_first_hit;
      }
      amrex::Print() << "\n";
    }

#endif // AMREX_USE_GPU

    if (performanceTimingEnabled()) {
      Gpu::streamSynchronize();
      reportPerformanceDuration(
          "gp_geometry_imagepoint_fused", lev,
          amrex::second() - timing_total_start,
          Long(gpstore.total_ngps), false);
    }

    // initialiseGPs rewrites marker component 1 from the temporary geometry id
    // to the final first-IP interpolation count.  Only valid cells own that
    // value; same-level ghost copies must be refreshed before any flux stencil
    // inspects them.  Exchange both components so the solid and GP views have
    // one authoritative, box-layout-independent representation.
    mfab.FillBoundary(0, 2, amr_p->Geom(lev).periodicity());

    // A geometrically valid marker halo is not proof that the corresponding
    // fine-level state exists.  Audit the exact nonzero interpolation supports
    // against the global valid BoxArray before any RHS evaluation can consume
    // an unfilled coarse-fine gap.
    auditGPAMRSupportCoverage(lev);

    // NOTE: gpstore.shrink() disabled — shrink_to_fit() on PODVector
    // corrupts data on some platforms (observed with AMReX ManagedVector on CPU).

    finishPerformanceTiming("initialise_gps_total", lev,
                            timing_total_start,
                            Long(gpstore.total_ngps), false);
  }


  bool rzAnnularCellAverageEnabled(int lev) const
  {
    if constexpr (rz_annular_bic_cell_average) {
      AMREX_ALWAYS_ASSERT(lev >= 0 && lev < static_cast<int>(bmf_a.size()));
      return amr_p->Geom(lev).IsRZ();
    }
    amrex::ignore_unused(lev);
    return false;
  }

  /**
   * Refresh the cell-centred primitive ghosts across the physical R-Z axis.
   *
   * The production state at every solid-side cell is still the one unique
   * shared GP state.  This routine only applies the coordinate parity of that
   * already-published primitive field to the negative-radius physical ghost
   * cells used by the full-background-grid flux stencil.  It neither creates
   * a second GP state nor changes a Cartesian face area or cell volume.
   */
  void fillRZAxisPrimitiveParity(MultiFab& prims_mf, int lev) const
  {
    const auto& geom = amr_p->Geom(lev);
    if (!geom.IsRZ() || geom.ProbLo(0) != Real(0.0)) return;

    const int axis_cell = geom.Domain().smallEnd(0);
    for (MFIter mfi(prims_mf, false); mfi.isValid(); ++mfi) {
      Box axis_ghosts = amrex::grow(mfi.validbox(), prims_mf.nGrow());
      axis_ghosts &= prims_mf[mfi].box();
      if (axis_ghosts.smallEnd(0) >= axis_cell) continue;
      axis_ghosts.setBig(
          0, amrex::min(axis_ghosts.bigEnd(0), axis_cell - 1));

      const int deepest_mirror =
          2 * axis_cell - 1 - axis_ghosts.smallEnd(0);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          deepest_mirror <= prims_mf[mfi].box().bigEnd(0),
          "R-Z axis parity requires every negative-radius primitive ghost "
          "to have a readable positive-radius mirror");

      const auto prims = prims_mf.array(mfi);
      ParallelFor(
          axis_ghosts, cls_t::NPRIM,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            const int mirror_i = 2 * axis_cell - 1 - i;
            Real value = prims(mirror_i, j, k, n);
            // The present two-dimensional R-Z closure treats only the radial
            // velocity as the symmetry-normal component.  This deliberately
            // matches CNS_setup's physical symmetry BC: QW is the dormant
            // Cartesian third component here, not an enabled swirl model.
            if (n == cls_t::QU) value = -value;
            prims(i, j, k, n) = value;
          });
    }
  }

  /// Recover pointwise radial centre states from annular conservative cell
  /// averages for both active-fluid and shared ghost cells.  During the
  /// active-fluid pass, a quadratic recovery is accepted only when all three
  /// annular inputs belong to the same active-fluid extension.  A wall-adjacent
  /// cell without such a stencil retains its local Q(Ubar) value until the
  /// unique shared GP has been reconstructed; no deep-solid state is sampled.
  void recoverRZCentrePrimitives(MultiFab& annular_state,
                                 MultiFab& prims_mf,
                                 const cls_t* cls,
                                 int lev,
                                 bool active_fluid_only) const
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        rzAnnularCellAverageEnabled(lev),
        "recoverRZCentrePrimitives requires the R-Z annular shared-GP mode");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        annular_state.boxArray() == prims_mf.boxArray() &&
            annular_state.DistributionMap() == prims_mf.DistributionMap(),
        "R-Z centre-state recovery requires matching state/primitive layouts");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        active_fluid_only,
        "R-Z production centre recovery is defined only on active fluid; "
        "the unique shared GP is constructed in the subsequent BI-CWLS "
        "phase");

    annular_state.FillBoundary(amr_p->Geom(lev).periodicity());
    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto dxyz = amr_p->Geom(lev).CellSizeArray();
    const Real radial_origin_over_dr = prob_lo[0] / dxyz[0];
    const int radial_domain_lo =
        amr_p->Geom(lev).Domain().smallEnd(0);
    const int grow = amrex::min(annular_state.nGrow(), prims_mf.nGrow());
    auto& markers = *bmf_a[lev];
    markers.FillBoundary(amr_p->Geom(lev).periodicity());
    for (MFIter mfi(prims_mf, false); mfi.isValid(); ++mfi) {
      Box bx = amrex::grow(mfi.validbox(), grow);
      bx &= annular_state[mfi].box();
      bx &= prims_mf[mfi].box();
      bx &= markers[mfi].box();
      if (active_fluid_only) {
        // Physical negative-radius ghosts are populated by parity only after
        // every positive-radius fluid/GP state has reached its stage-local
        // final value.  Never deconvolve those ghosts independently.
        bx.setSmall(0, amrex::max(bx.smallEnd(0), radial_domain_lo));
      }
      const auto state = annular_state.const_array(mfi);
      const auto prims = prims_mf.array(mfi);
      const auto marker = markers.const_array(mfi);
      const int radial_data_lo = amrex::max(
          active_fluid_only ? radial_domain_lo
                            : annular_state[mfi].box().smallEnd(0),
          amrex::max(annular_state[mfi].box().smallEnd(0),
                     markers[mfi].box().smallEnd(0)));
      const int radial_data_hi = amrex::min(
          annular_state[mfi].box().bigEnd(0),
          markers[mfi].box().bigEnd(0));
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          radial_data_hi - radial_data_lo >= 2,
          "R-Z centre-state recovery requires three radial annular cells");
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (active_fluid_only && marker(i, j, k, 0) != 0) return;
        Real point_state[cls_t::NCONS];
        Real point_primitives[cls_t::NPRIM];
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          point_primitives[n] = prims(i, j, k, n);
        }
        if (active_fluid_only) {
          // Prefer centred, then forward, then backward quadratic recovery.
          // Every selected state is an active background-grid fluid cell;
          // neither a distinct face state nor cut geometry is introduced.
          const int candidate_lo[3] = {i - 1, i, i - 2};
          int stencil_lo = i;
          bool found_same_fluid_stencil = false;
          for (int candidate = 0;
               candidate < 3 && !found_same_fluid_stencil;
               ++candidate) {
            const int lo = candidate_lo[candidate];
            if (lo < radial_data_lo || lo + 2 > radial_data_hi) continue;
            bool all_active_fluid = true;
            for (int q = 0; q < 3; ++q) {
              all_active_fluid =
                  all_active_fluid && marker(lo + q, j, k, 0) == 0;
            }
            if (all_active_fluid) {
              stencil_lo = lo;
              found_same_fluid_stencil = true;
            }
          }
          if (!found_same_fluid_stencil) return;
          cerisse::rz_fv::annular_average_to_radial_center_from_stencil<
              cls_t::NCONS>(i, j, k, stencil_lo, radial_origin_over_dr,
                            state, point_state);
        } else {
          cerisse::rz_fv::annular_average_to_radial_center<cls_t::NCONS>(
              i, j, k, radial_data_lo, radial_data_hi,
              radial_origin_over_dr, state, point_state);
        }
        cls->cons2prims_point(point_state, point_primitives);
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(i, j, k, n) = point_primitives[n];
        }
      });
    }
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
  }

  /// Recover Cartesian point-valued BI-CWLS support states from conservative
  /// full-cell averages without crossing the immersed boundary. Two separable
  /// quadratic finite-volume passes remove the x and y cell averages. A cell
  /// without a complete same-fluid stencil retains the existing Q(Ubar)
  /// full-state fallback in support_prims.
  std::array<Long, 4> recoverCartesianBICSupportCentrePrimitives(
      MultiFab& conservative_state, MultiFab& support_prims,
      const cls_t* cls, int lev, bool collect_stats = false) const
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nsGPCellAverageRecoveryEnabled(),
        "Cartesian BI-CWLS support recovery requires "
        "ib.ns_gp_cell_average_recovery=1");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !amr_p->Geom(lev).IsRZ() && AMREX_SPACEDIM == 2,
        "Cartesian BI-CWLS support recovery requires two-dimensional "
        "Cartesian geometry");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        conservative_state.boxArray() == support_prims.boxArray() &&
            conservative_state.DistributionMap() ==
                support_prims.DistributionMap(),
        "Cartesian BI-CWLS support recovery requires matching layouts");

    auto& geom = amr_p->Geom(lev);
    conservative_state.FillBoundary(geom.periodicity());
    auto& markers = *bmf_a[lev];
    markers.FillBoundary(geom.periodicity());

    constexpr int x_success_comp = cls_t::NCONS;
    constexpr int y_stencil_success_comp = cls_t::NCONS + 1;
    constexpr int point_success_comp = cls_t::NCONS + 2;
    MultiFab recovered_centres(
        support_prims.boxArray(), support_prims.DistributionMap(),
        cls_t::NCONS + 3, 2, MFInfo().SetArena(The_Async_Arena()));
    recovered_centres.setVal(Real(0.0));

    const Box domain = geom.Domain();
    const int domain_x_lo = domain.smallEnd(0);
    const int domain_x_hi = domain.bigEnd(0);
    const int domain_y_lo = domain.smallEnd(1);
    const int domain_y_hi = domain.bigEnd(1);

    // Remove the x average at fixed y.
    for (MFIter mfi(recovered_centres, false); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();
      box &= conservative_state[mfi].box();
      box &= markers[mfi].box();
      const auto state = conservative_state.const_array(mfi);
      const auto recovered = recovered_centres.array(mfi);
      const auto marker = markers.const_array(mfi);
      const Box state_box = conservative_state[mfi].box();
      const Box marker_box = markers[mfi].box();
      const int data_lo = amrex::max(
          domain_x_lo,
          amrex::max(state_box.smallEnd(0), marker_box.smallEnd(0)));
      const int data_hi = amrex::min(
          domain_x_hi,
          amrex::min(state_box.bigEnd(0), marker_box.bigEnd(0)));
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (marker(i, j, k, 0) != 0) return;

        const int candidate_lo[3] = {i - 1, i, i - 2};
        int stencil_lo = i;
        bool found = false;
        for (int candidate = 0; candidate < 3 && !found; ++candidate) {
          const int lo = candidate_lo[candidate];
          if (lo < data_lo || lo + 2 > data_hi) continue;
          bool all_fluid = true;
          for (int q = 0; q < 3; ++q) {
            all_fluid =
                all_fluid && marker(lo + q, j, k, 0) == 0;
          }
          if (all_fluid) {
            stencil_lo = lo;
            found = true;
          }
        }
        if (!found) return;

        Real weight[3];
        if (stencil_lo == i - 1) {
          weight[0] = Real(-1.0 / 24.0);
          weight[1] = Real(13.0 / 12.0);
          weight[2] = Real(-1.0 / 24.0);
        } else if (stencil_lo == i) {
          weight[0] = Real(23.0 / 24.0);
          weight[1] = Real(1.0 / 12.0);
          weight[2] = Real(-1.0 / 24.0);
        } else {
          weight[0] = Real(-1.0 / 24.0);
          weight[1] = Real(1.0 / 12.0);
          weight[2] = Real(23.0 / 24.0);
        }
        for (int n = 0; n < cls_t::NCONS; ++n) {
          recovered(i, j, k, n) =
              weight[0] * state(stencil_lo, j, k, n) +
              weight[1] * state(stencil_lo + 1, j, k, n) +
              weight[2] * state(stencil_lo + 2, j, k, n);
        }
        recovered(i, j, k, x_success_comp) = Real(1.0);
      });
    }
    recovered_centres.FillBoundary(geom.periodicity());

    // Remove the y average from the x-centre states, then apply the EOS only
    // after a complete conservative point state has been recovered.
    for (MFIter mfi(recovered_centres, false); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();
      box &= support_prims[mfi].box();
      box &= markers[mfi].box();
      const auto recovered = recovered_centres.array(mfi);
      const auto support = support_prims.array(mfi);
      const auto marker = markers.const_array(mfi);
      const Box recovered_box = recovered_centres[mfi].box();
      const Box marker_box = markers[mfi].box();
      const int data_lo = amrex::max(
          domain_y_lo,
          amrex::max(recovered_box.smallEnd(1), marker_box.smallEnd(1)));
      const int data_hi = amrex::min(
          domain_y_hi,
          amrex::min(recovered_box.bigEnd(1), marker_box.bigEnd(1)));
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (marker(i, j, k, 0) != 0 ||
            recovered(i, j, k, x_success_comp) < Real(0.5)) {
          return;
        }

        const int candidate_lo[3] = {j - 1, j, j - 2};
        int stencil_lo = j;
        bool found = false;
        for (int candidate = 0; candidate < 3 && !found; ++candidate) {
          const int lo = candidate_lo[candidate];
          if (lo < data_lo || lo + 2 > data_hi) continue;
          bool all_fluid = true;
          for (int q = 0; q < 3; ++q) {
            all_fluid =
                all_fluid && marker(i, lo + q, k, 0) == 0 &&
                recovered(i, lo + q, k, x_success_comp) > Real(0.5);
          }
          if (all_fluid) {
            stencil_lo = lo;
            found = true;
          }
        }
        if (!found) return;
        recovered(i, j, k, y_stencil_success_comp) = Real(1.0);

        Real weight[3];
        if (stencil_lo == j - 1) {
          weight[0] = Real(-1.0 / 24.0);
          weight[1] = Real(13.0 / 12.0);
          weight[2] = Real(-1.0 / 24.0);
        } else if (stencil_lo == j) {
          weight[0] = Real(23.0 / 24.0);
          weight[1] = Real(1.0 / 12.0);
          weight[2] = Real(-1.0 / 24.0);
        } else {
          weight[0] = Real(-1.0 / 24.0);
          weight[1] = Real(1.0 / 12.0);
          weight[2] = Real(23.0 / 24.0);
        }

        Real point_state[cls_t::NCONS];
        bool admissible = true;
        for (int n = 0; n < cls_t::NCONS; ++n) {
          point_state[n] =
              weight[0] * recovered(i, stencil_lo, k, n) +
              weight[1] * recovered(i, stencil_lo + 1, k, n) +
              weight[2] * recovered(i, stencil_lo + 2, k, n);
          admissible =
              admissible && amrex::Math::isfinite(point_state[n]);
        }
        const Real density = point_state[cls_t::URHO];
        if (!(density > Real(0.0))) admissible = false;
        if (admissible) {
          const Real kinetic =
              Real(0.5) *
              (point_state[cls_t::UMX] * point_state[cls_t::UMX] +
               point_state[cls_t::UMY] * point_state[cls_t::UMY] +
               point_state[cls_t::UMZ] * point_state[cls_t::UMZ]) /
              density;
          const Real internal_energy_density =
              point_state[cls_t::UET] - kinetic;
          admissible =
              amrex::Math::isfinite(internal_energy_density) &&
              internal_energy_density > Real(0.0);
        }
        if (!admissible) return;

        Real point_primitives[cls_t::NPRIM] = {Real(0.0)};
        cls->cons2prims_point(point_state, point_primitives);
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          if (!amrex::Math::isfinite(point_primitives[n])) return;
        }
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          support(i, j, k, n) = point_primitives[n];
        }
        recovered(i, j, k, point_success_comp) = Real(1.0);
      });
    }
    support_prims.FillBoundary(geom.periodicity());

    std::array<Long, 4> counts{{Long(0), Long(0), Long(0), Long(0)}};
    amrex::ignore_unused(collect_stats);
    return counts;
  }

  /// Evaluate the existing BI-constrained no-slip extension at four Cartesian
  /// Gauss points and publish one thermodynamically consistent conservative
  /// average at each real solid-side ghost-cell centre.
  void reconstructCartesianConservativeGhostAverageFromPointSupport(
      MultiFab& point_support_prims, const cls_t* cls, int lev)
  {
    if constexpr (!bic_point_target_functionals) {
      amrex::Abort(
          "Cartesian conservative GP averaging requires cached BI-CWLS point "
          "functionals");
    } else {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          !amr_p->Geom(lev).IsRZ() && AMREX_SPACEDIM == 2,
          "Cartesian conservative GP averaging requires two-dimensional "
          "Cartesian geometry");

      auto& gpstore = gpstore_a[lev];
      if (gpstore.total_ngps == 0) return;

      MultiFab conservative_average(
          point_support_prims.boxArray(),
          point_support_prims.DistributionMap(), cls_t::NCONS,
          point_support_prims.nGrow(),
          MFInfo().SetArena(The_Async_Arena()));
      conservative_average.setVal(Real(0.0));

      const int nfabs_local = gpstore.nfabs;
      Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
      Gpu::DeviceVector<Array4<Real>> d_average(nfabs_local);
      {
        Vector<Array4<Real>> h_prims(nfabs_local);
        Vector<Array4<Real>> h_average(nfabs_local);
        int ifab = 0;
        for (MFIter mfi(*bmf_a[lev], false); mfi.isValid(); ++mfi, ++ifab) {
          h_prims[ifab] = point_support_prims.array(mfi);
          h_average[ifab] = conservative_average.array(mfi);
        }
        AMREX_ALWAYS_ASSERT(ifab == nfabs_local);
        Gpu::copyAsync(
            Gpu::hostToDevice, h_prims.begin(), h_prims.end(),
            d_prims.begin());
        Gpu::copyAsync(
            Gpu::hostToDevice, h_average.begin(), h_average.end(),
            d_average.begin());
      }

      auto gpview = gpstore.view();
      auto* prims_arr = d_prims.data();
      auto* average_arr = d_average.data();
      const int total_ngps = gpview.total_ngps;
      constexpr Real quadrature_weight = Real(0.25);
      for (int target = 0; target < IBM_RZ_BIC_GAUSS_POINTS; ++target) {
        computeAllGPs(point_support_prims, cls, lev, target);
        ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE(int ii) noexcept {
          const int ifab = gpview.gp_fab[ii];
          const auto prims = prims_arr[ifab];
          const auto average = average_arr[ifab];
          const int i = gpview.gp_ijk[ii](0);
          const int j = gpview.gp_ijk[ii](1);
          const int k = 0;
          const IntVect cell(AMREX_D_DECL(i, j, k));
          Real conservative[cls_t::NCONS];
          cls->prims2cons(cell, prims, conservative);
          for (int n = 0; n < cls_t::NCONS; ++n) {
            average(i, j, k, n) +=
                quadrature_weight * conservative[n];
          }
        });
      }

      Gpu::DeviceScalar<int> invalid_average(0);
      int* invalid_average_ptr = invalid_average.dataPtr();
      ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE(int ii) noexcept {
        const int ifab = gpview.gp_fab[ii];
        const auto prims = prims_arr[ifab];
        const auto average = average_arr[ifab];
        const int i = gpview.gp_ijk[ii](0);
        const int j = gpview.gp_ijk[ii](1);
        const int k = 0;
        Real conservative[cls_t::NCONS];
        bool admissible = true;
        for (int n = 0; n < cls_t::NCONS; ++n) {
          conservative[n] = average(i, j, k, n);
          admissible =
              admissible && amrex::Math::isfinite(conservative[n]);
        }
        const Real density = conservative[cls_t::URHO];
        if (!(density > Real(0.0))) admissible = false;
        if (admissible) {
          const Real kinetic =
              Real(0.5) *
              (conservative[cls_t::UMX] * conservative[cls_t::UMX] +
               conservative[cls_t::UMY] * conservative[cls_t::UMY] +
               conservative[cls_t::UMZ] * conservative[cls_t::UMZ]) /
              density;
          const Real internal_energy_density =
              conservative[cls_t::UET] - kinetic;
          admissible =
              amrex::Math::isfinite(internal_energy_density) &&
              internal_energy_density > Real(0.0);
        }
        if (!admissible) {
          Gpu::Atomic::Add(invalid_average_ptr, 1);
          return;
        }
        Real primitive[cls_t::NPRIM] = {Real(0.0)};
        cls->cons2prims_point(conservative, primitive);
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(i, j, k, n) = primitive[n];
        }
      });
      int invalid_average_count = invalid_average.dataValue();
      ParallelDescriptor::ReduceIntSum(invalid_average_count);
      if (invalid_average_count > 0) {
        amrex::Abort(
            "Cartesian conservative GP quadrature produced " +
            std::to_string(invalid_average_count) +
            " inadmissible shared ghost-cell average(s)");
      }
      point_support_prims.FillBoundary(amr_p->Geom(lev).periodicity());
    }
  }

  /// Complete Cartesian no-slip finite-volume shared-GP data path.
  void computeAllGPsCartesianConservativeAverage(
      MultiFab& prims_mf, MultiFab& conservative_state,
      const cls_t* cls, int lev)
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nsGPCellAverageRecoveryEnabled(),
        "Cartesian conservative shared-GP path requires "
        "ib.ns_gp_cell_average_recovery=1");

    MultiFab support_prims(
        prims_mf.boxArray(), prims_mf.DistributionMap(), prims_mf.nComp(),
        prims_mf.nGrow(), MFInfo().SetArena(The_Async_Arena()));
    MultiFab::Copy(
        support_prims, prims_mf, 0, 0, prims_mf.nComp(),
        prims_mf.nGrow());
    const auto recovery_counts =
        recoverCartesianBICSupportCentrePrimitives(
            conservative_state, support_prims, cls, lev,
            false
        );
    amrex::ignore_unused(recovery_counts);

    reconstructCartesianConservativeGhostAverageFromPointSupport(
        support_prims, cls, lev);

    auto& markers = *bmf_a[lev];
    for (MFIter mfi(prims_mf, false); mfi.isValid(); ++mfi) {
      Box box = amrex::grow(mfi.validbox(), prims_mf.nGrow());
      box &= prims_mf[mfi].box();
      box &= support_prims[mfi].box();
      box &= markers[mfi].box();
      const auto prims = prims_mf.array(mfi);
      const auto support = support_prims.const_array(mfi);
      const auto marker = markers.const_array(mfi);
      ParallelFor(
          box, cls_t::NPRIM,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (marker(i, j, k, 0) != 0 &&
                marker(i, j, k, 1) != 0) {
              prims(i, j, k, n) = support(i, j, k, n);
            }
          });
    }
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
  }

  /// Recover point-valued BI-CWLS support states from annular conservative
  /// cell averages without crossing the immersed boundary. The two separable
  /// quadratic finite-volume passes use centred, forward, then backward
  /// three-cell stencils. A cell with no complete same-fluid stencil retains
  /// the radial-centre production value already present in support_prims.
  std::array<Long, 3> recoverRZBICSupportCentrePrimitives(
      MultiFab& annular_state, MultiFab& support_prims,
      const cls_t* cls, int lev,
      bool collect_stats = false) const
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        rzAnnularCellAverageEnabled(lev),
        "R-Z BI-CWLS support recovery requires annular shared-GP semantics");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        annular_state.boxArray() == support_prims.boxArray() &&
            annular_state.DistributionMap() ==
                support_prims.DistributionMap(),
        "R-Z BI-CWLS support recovery requires matching layouts");

    auto& geom = amr_p->Geom(lev);
    annular_state.FillBoundary(geom.periodicity());
    auto& markers = *bmf_a[lev];
    markers.FillBoundary(geom.periodicity());

    constexpr int radial_success_comp = cls_t::NCONS;
    constexpr int axial_success_comp = cls_t::NCONS + 1;
    MultiFab recovered_centres(
        support_prims.boxArray(), support_prims.DistributionMap(),
        cls_t::NCONS + 2, 2, MFInfo().SetArena(The_Async_Arena()));
    recovered_centres.setVal(Real(0.0));

    const auto prob_lo = geom.ProbLoArray();
    const auto dxyz = geom.CellSizeArray();
    const Real radial_origin_over_dr = prob_lo[0] / dxyz[0];
    const Box domain = geom.Domain();
    const int domain_radial_lo = domain.smallEnd(0);
    const int domain_radial_hi = domain.bigEnd(0);
    const int domain_axial_lo = domain.smallEnd(1);
    const int domain_axial_hi = domain.bigEnd(1);

    // First remove the r-weighted annular average at fixed axial index.
    for (MFIter mfi(recovered_centres, false); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();
      box &= annular_state[mfi].box();
      box &= markers[mfi].box();
      const auto state = annular_state.const_array(mfi);
      const auto recovered = recovered_centres.array(mfi);
      const auto marker = markers.const_array(mfi);
      const Box state_box = annular_state[mfi].box();
      const Box marker_box = markers[mfi].box();
      const int radial_data_lo =
          amrex::max(domain_radial_lo,
                     amrex::max(state_box.smallEnd(0),
                                marker_box.smallEnd(0)));
      const int radial_data_hi =
          amrex::min(domain_radial_hi,
                     amrex::min(state_box.bigEnd(0),
                                marker_box.bigEnd(0)));
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (marker(i, j, k, 0) != 0) return;

        const int candidate_lo[3] = {i - 1, i, i - 2};
        int stencil_lo = i;
        bool found = false;
        for (int candidate = 0; candidate < 3 && !found; ++candidate) {
          const int lo = candidate_lo[candidate];
          if (lo < radial_data_lo || lo + 2 > radial_data_hi) continue;
          bool all_fluid = true;
          for (int q = 0; q < 3; ++q) {
            all_fluid = all_fluid && marker(lo + q, j, k, 0) == 0;
          }
          if (all_fluid) {
            stencil_lo = lo;
            found = true;
          }
        }
        if (!found) return;

        Real radial_centre[cls_t::NCONS];
        cerisse::rz_fv::annular_average_to_radial_center_from_stencil<
            cls_t::NCONS>(i, j, k, stencil_lo, radial_origin_over_dr,
                          state, radial_centre);
        for (int n = 0; n < cls_t::NCONS; ++n) {
          recovered(i, j, k, n) = radial_centre[n];
        }
        recovered(i, j, k, radial_success_comp) = Real(1.0);
      });
    }
    recovered_centres.FillBoundary(geom.periodicity());

    // Then remove the uniform axial average from the radial-centre values.
    for (MFIter mfi(recovered_centres, false); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();
      box &= support_prims[mfi].box();
      box &= markers[mfi].box();
      const auto recovered = recovered_centres.array(mfi);
      const auto support = support_prims.array(mfi);
      const auto marker = markers.const_array(mfi);
      const Box recovered_box = recovered_centres[mfi].box();
      const Box marker_box = markers[mfi].box();
      const int axial_data_lo =
          amrex::max(domain_axial_lo,
                     amrex::max(recovered_box.smallEnd(1),
                                marker_box.smallEnd(1)));
      const int axial_data_hi =
          amrex::min(domain_axial_hi,
                     amrex::min(recovered_box.bigEnd(1),
                                marker_box.bigEnd(1)));
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (marker(i, j, k, 0) != 0 ||
            recovered(i, j, k, radial_success_comp) < Real(0.5)) {
          return;
        }

        const int candidate_lo[3] = {j - 1, j, j - 2};
        int stencil_lo = j;
        bool found = false;
        for (int candidate = 0; candidate < 3 && !found; ++candidate) {
          const int lo = candidate_lo[candidate];
          if (lo < axial_data_lo || lo + 2 > axial_data_hi) continue;
          bool all_fluid = true;
          for (int q = 0; q < 3; ++q) {
            all_fluid =
                all_fluid && marker(i, lo + q, k, 0) == 0 &&
                recovered(i, lo + q, k, radial_success_comp) > Real(0.5);
          }
          if (all_fluid) {
            stencil_lo = lo;
            found = true;
          }
        }
        if (!found) return;

        Real weight[3];
        if (stencil_lo == j - 1) {
          weight[0] = Real(-1.0 / 24.0);
          weight[1] = Real(13.0 / 12.0);
          weight[2] = Real(-1.0 / 24.0);
        } else if (stencil_lo == j) {
          weight[0] = Real(23.0 / 24.0);
          weight[1] = Real(1.0 / 12.0);
          weight[2] = Real(-1.0 / 24.0);
        } else {
          weight[0] = Real(-1.0 / 24.0);
          weight[1] = Real(1.0 / 12.0);
          weight[2] = Real(23.0 / 24.0);
        }

        Real point_state[cls_t::NCONS];
        Real point_primitives[cls_t::NPRIM] = {Real(0.0)};
        for (int n = 0; n < cls_t::NCONS; ++n) {
          point_state[n] =
              weight[0] * recovered(i, stencil_lo, k, n) +
              weight[1] * recovered(i, stencil_lo + 1, k, n) +
              weight[2] * recovered(i, stencil_lo + 2, k, n);
        }
        cls->cons2prims_point(point_state, point_primitives);
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          support(i, j, k, n) = point_primitives[n];
        }
        recovered(i, j, k, axial_success_comp) = Real(1.0);
      });
    }
    support_prims.FillBoundary(geom.periodicity());

    std::array<Long, 3> counts{{Long(0), Long(0), Long(0)}};
    amrex::ignore_unused(collect_stats);
    return counts;
  }

  /// Accumulate one thermodynamically consistent Gauss-point state into one
  /// of the three annular target cells attached to each real shared GP.
  /// These target averages are auxiliary evaluations of the same continuous
  /// BI-CWLS extension; they are not additional ghost cells or face states.
  void accumulateRZAnnularGhostState(
      Array4<Real>* prims_arr, const cls_t* cls, int lev,
      int point_target) const
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        rzAnnularCellAverageEnabled(lev),
        "accumulateRZAnnularGhostState requires the R-Z annular mode");
    AMREX_ALWAYS_ASSERT(
        point_target >= 0 && point_target < IBM_RZ_BIC_ANNULAR_TARGETS);
    auto& gpstore = gpstore_a[lev];
    auto gpview = gpstore.view();
    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto dxyz = amr_p->Geom(lev).CellSizeArray();
    constexpr Real gauss_node =
        Real(0.57735026918962576450914878050195746);
    const int radial_slot = point_target / IBM_RZ_BIC_GAUSS_POINTS;
    const int quadrature_point =
        point_target % IBM_RZ_BIC_GAUSS_POINTS;
    const int radial_cell_offset = radial_slot - 1;
    const int total_ngps = gpview.total_ngps;
    ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE(int ii) noexcept {
        const int ifab = gpview.gp_fab[ii];
        const auto prims = prims_arr[ifab];
        const int i = gpview.gp_ijk[ii](0);
        const int j = gpview.gp_ijk[ii](1);
#if (AMREX_SPACEDIM == 3)
        const int k = gpview.gp_ijk[ii](2);
#else
        const int k = 0;
#endif
        const IntVect iv(AMREX_D_DECL(i, j, k));
        Real conservative[cls_t::NCONS];
        cls->prims2cons(iv, prims, conservative);
        const Real radial_sign =
            (quadrature_point & 1) == 0 ? Real(-1.0) : Real(1.0);
        const Real radial_centre =
            prob_lo[0] +
            (Real(i + radial_cell_offset) + Real(0.5)) * dxyz[0];
        const Real radial_quadrature = radial_centre +
            Real(0.5) * radial_sign * gauss_node * dxyz[0];
        const Real weight =
            radial_quadrature / (Real(4.0) * radial_centre);
        for (int n = 0; n < cls_t::NCONS; ++n) {
          if (quadrature_point == 0) {
            gpview.rz_bic_annular_stencil_w[ii](radial_slot, n) =
                weight * conservative[n];
          } else {
            gpview.rz_bic_annular_stencil_w[ii](radial_slot, n) +=
                weight * conservative[n];
          }
        }
    });
  }

  /// Recover the unique GP-centre state from the three annular averages and
  /// publish the central annular average to the conservative ghost cell.
  void finalizeRZAnnularGhostStates(
      Array4<Real>* prims_arr, Array4<Real>* annular_arr,
      const cls_t* cls, int lev)
  {
    auto& gpstore = gpstore_a[lev];
    auto gpview = gpstore.view();
    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto dxyz = amr_p->Geom(lev).CellSizeArray();
    const Real radial_origin_over_dr = prob_lo[0] / dxyz[0];
    const int total_ngps = gpview.total_ngps;
    ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE(int ii) noexcept {
      const int ifab = gpview.gp_fab[ii];
      const auto prims = prims_arr[ifab];
      const auto state = annular_arr[ifab];
      const int i = gpview.gp_ijk[ii](0);
      const int j = gpview.gp_ijk[ii](1);
#if (AMREX_SPACEDIM == 3)
      const int k = gpview.gp_ijk[ii](2);
#else
      const int k = 0;
#endif
      const auto& annular_stencil =
          gpview.rz_bic_annular_stencil[ii];
      Real point_state[cls_t::NCONS];
      Real point_primitives[cls_t::NPRIM] = {Real(0.0)};
      cerisse::rz_fv::annular_stencil_to_radial_center<cls_t::NCONS>(
          i, radial_origin_over_dr, annular_stencil, point_state);
      cls->cons2prims_point(point_state, point_primitives);
      for (int n = 0; n < cls_t::NPRIM; ++n) {
        prims(i, j, k, n) = point_primitives[n];
      }
      for (int n = 0; n < cls_t::NCONS; ++n) {
        state(i, j, k, n) = annular_stencil(1, n);
      }
    });
  }

  /// Reconstruct the R-Z shared GP from the primitive support field currently
  /// stored in prims_mf using the production target and centre-recovery
  /// operators.
  void reconstructRZAnnularGhostStates(MultiFab& prims_mf,
                                       MultiFab& annular_state,
                                       const cls_t* cls,
                                       int lev)
  {
    auto& gpstore = gpstore_a[lev];
    const int nfabs_local = gpstore.nfabs;
    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    Gpu::DeviceVector<Array4<Real>> d_annular(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      Vector<Array4<Real>> h_annular(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(*bmf_a[lev], false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
        h_annular[ifab] = annular_state.array(mfi);
      }
      AMREX_ALWAYS_ASSERT(ifab == nfabs_local);
      Gpu::copyAsync(
          Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      Gpu::copyAsync(Gpu::hostToDevice, h_annular.begin(), h_annular.end(),
                     d_annular.begin());
    }

    for (int target = 0; target < IBM_RZ_BIC_ANNULAR_TARGETS; ++target) {
      computeAllGPs(prims_mf, cls, lev, target);
      accumulateRZAnnularGhostState(d_prims.data(), cls, lev, target);
    }
    finalizeRZAnnularGhostStates(
        d_prims.data(), d_annular.data(), cls, lev);
    annular_state.FillBoundary(amr_p->Geom(lev).periodicity());
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
  }

  /// Complete RZ-2 data path. Fluid supports and the unique shared GP state
  /// both use annular conservative averages followed by the same radial
  /// centre-state recovery.
  void computeAllGPsRZAnnular(MultiFab& prims_mf,
                              MultiFab& annular_state,
                              const cls_t* cls,
                              int lev)
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        rzAnnularCellAverageEnabled(lev),
        "computeAllGPsRZAnnular called outside the R-Z annular mode");
    recoverRZCentrePrimitives(
        annular_state, prims_mf, cls, lev, /*active_fluid_only=*/true);
    // Phase A has changed positive-radius active-fluid centre values.  Refresh
    // their physical-axis mirrors before any BI support pass can read them.
    fillRZAxisPrimitiveParity(prims_mf, lev);
    if constexpr (!rz_bic_support_centre_recovery) {
      reconstructRZAnnularGhostStates(prims_mf, annular_state, cls, lev);
      fillRZAxisPrimitiveParity(prims_mf, lev);
      return;
    }

    // BI-CWLS consumes meridional point states, whereas the R-Z WENO field
    // retains the RZ-1 centre-state semantics above. Reconstruct the unique
    // shared GP from a separate support field, then publish only that GP back
    // to the production primitive field.
    MultiFab support_prims(
        prims_mf.boxArray(), prims_mf.DistributionMap(), prims_mf.nComp(),
        prims_mf.nGrow(), MFInfo().SetArena(The_Async_Arena()));
    MultiFab::Copy(support_prims, prims_mf, 0, 0, prims_mf.nComp(),
                   prims_mf.nGrow());
    recoverRZBICSupportCentrePrimitives(
        annular_state, support_prims, cls, lev);
    reconstructRZAnnularGhostStates(
        support_prims, annular_state, cls, lev);

    auto& markers = *bmf_a[lev];
    for (MFIter mfi(prims_mf, false); mfi.isValid(); ++mfi) {
      Box box = amrex::grow(mfi.validbox(), prims_mf.nGrow());
      box &= prims_mf[mfi].box();
      box &= support_prims[mfi].box();
      box &= markers[mfi].box();
      const auto prims = prims_mf.array(mfi);
      const auto support = support_prims.const_array(mfi);
      const auto marker = markers.const_array(mfi);
      ParallelFor(
          box, cls_t::NPRIM,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (marker(i, j, k, 0) != 0 &&
                marker(i, j, k, 1) != 0) {
              prims(i, j, k, n) = support(i, j, k, n);
            }
          });
    }
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
    // FillBoundary does not impose a physical-axis boundary condition.  Do
    // this last, after the unique shared GP has been copied into production,
    // so a mirrored solid-side GP cannot retain a pre-reconstruction value.
    fillRZAxisPrimitiveParity(prims_mf, lev);
  }

  // ========================================================================
  // Sparse first-hit crossing-face geometry and reconstruction metadata.
  // Strict crossing-face invariants remain part of production builds.
  void computeAllGPs(MultiFab& prims_mf,
                     const cls_t* cls,
                     int lev,
                     int rz_point_target = -1)
  {
    BL_PROFILE("IBM::computeAllGPs");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        prims_mf.nGrow() >= volumeInterpolationNghost(lev),
        "computeAllGPs requires the IBM volume interpolation halo");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        rz_point_target >= -1 && rz_point_target < IBM_RZ_BIC_TARGETS,
        "invalid BI-CWLS point target");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        rz_point_target < 0 ||
            rzAnnularCellAverageEnabled(lev) ||
            (ns_cartesian_bic_cell_average &&
             !amr_p->Geom(lev).IsRZ() &&
             rz_point_target < IBM_RZ_BIC_GAUSS_POINTS),
        "BI-CWLS point targets require R-Z annular mode or an enabled "
        "Cartesian conservative ghost-average path");

    // Publish valid fluid primitives into the wider image-point halo before
    // any GP reads them.  This same-level exchange is what makes the cached
    // global support indices independent of FAB and MPI ownership.
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());

    auto& gpstore = gpstore_a[lev];
    if (gpstore.total_ngps > 0) {

    auto gpview = gpstore.view();
    auto const* lf_ptr = LocalFrame_a.data();
    auto const* se_ptr = SurfElem_a.data();
    const int nfabs_local = gpview.nfabs;

    // Capture per-geometry transforms for rotating body-frame LocalFrame to world
    AMREX_ALWAYS_ASSERT(ngeom <= MAX_NGEOM);
    const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    for (int ii = 0; ii < ngeom_local; ii++) {
        transforms[ii] = transform_a[ii];
    }

    // Build device array of Array4 pointers (one per local FAB).
    // No snapshot needed: GPs write only to ghost cells (ibMarkers==1)
    // while image-point interpolation reads only from fluid cells (ibMarkers==0).
    auto& mfab = *bmf_a[lev];

    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
      }
      Gpu::copyAsync(Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      // No streamSynchronize: copyAsync and ParallelFor share the same stream.
    }

    auto* prims_arr = d_prims.data();

    // Single kernel launch over all ghost points on this level
    auto* copy = this;
    const int total_ngps = gpview.total_ngps;
    const auto prob_lo_lev = amr_p->Geom(lev).ProbLoArray();
    const auto dx_lev = dx_a[lev];
    const bool cartesian_geometry = !amr_p->Geom(lev).IsRZ();

    ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE (int ii) noexcept
    {
      // nvcc does not permit an extended device lambda to first-capture a
      // variable inside an if-constexpr branch.  These geometry values are
      // consumed only by the optional NS compatibility policy, so make their
      // capture unconditional without generating work.
      amrex::ignore_unused(prob_lo_lev, dx_lev, cls, cartesian_geometry,
                           se_ptr, rz_point_target);
      // Determine which FAB this GP belongs to
      int ifab = gpview.gp_fab[ii];
      auto prims = prims_arr[ifab];  // read image points & write ghost cells

      // 1) Reconstruct local orthonormal frame — rotate from body to world frame
      int elem_idx = gpview.elemIdx[ii];
      const auto& frame = lf_ptr[elem_idx];

      // Find which geometry this GP belongs to (for transform lookup)
      int geomIdx = gpview.geomIdx[ii];
      const auto& xform = transforms[geomIdx];

      int type_solid_bc = 0;

      // Rotate body-frame vectors to world frame, then apply an optional
      // analytic verification frame at the cached wall point.
      LocalFrame localframe;
      xform.rotate_to_world(frame.normal,   localframe.normal);
      xform.rotate_to_world(frame.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
      xform.rotate_to_world(frame.tangent2, localframe.tangent2);
#endif
      apply_volume_surface_frame_override(gpview.ib_xyz[ii], localframe);

      Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          nvec(d)  = localframe.normal[d];
          t1vec(d) = localframe.tangent1[d];
#if (AMREX_SPACEDIM == 3)
          t2vec(d) = localframe.tangent2[d];
#else
          t2vec(d) = Real(0.0);
#endif
      }

      // 2) Storage for primitive variables along the normal
      Array2D<Real, 0, eorder_tparm + 1, 0, cls_t::NPRIM - 1> primsNormal;
      for (int p = 0; p <= eorder_tparm + 1; ++p) {
          for (int n = 0; n < cls_t::NPRIM; ++n) {
              primsNormal(p,n) = Real(0.0);
          }
      }

      // 3) Interpolate primitive variables at all image points
      copy->template interpolateIMs<eorder_tparm, iorder_tparm>(
          gpview.imp_ip_ijk[ii], gpview.imp_ipweights[ii], prims, primsNormal);

      // 4) Transform velocities at image points to local frame
      for (int iip = 2; iip < 2 + eorder_tparm; ++iip) {
          copy->template global2local<eorder_tparm>(iip, primsNormal, nvec, t1vec, t2vec);
      }

      // 5) Apply wall model at IB surface.
      //    n_valid_bc = number of leading valid image points. The wall model
      //    uses the 2nd-order one-sided form only when >=2 are available and
      //    extrap_order>=2; otherwise it degrades to the 1st-order single-IP
      //    form (identical to the previous behaviour at extrap_order=1).
      int n_valid_bc = 0;
      for (int kk = 0; kk < eorder_tparm; ++kk) {
          if (gpview.imp_ninterp[ii](kk) < INTERP_THRESHOLD_GP) break;
          n_valid_bc = kk + 1;
      }
      ibm_pressure_compatibility_t pressure_compat{};
      if constexpr (use_ns_noslip_pressure_compatibility) {
        if (n_valid_bc >= 2) {
          pressure_compat =
              copy->template ibm_viscous_pressure_compatibility<
                  iorder_tparm, N_InterP>(
                      gpview.imp_ip_ijk[ii], gpview.imp_xyz[ii],
                      gpview.imp_ninterp[ii](0), prims, nvec,
                      prob_lo_lev, dx_lev, cls, cartesian_geometry);
        }
        if constexpr (use_bi_constrained_shared_gp &&
                      expanded_bic_support) {
          if (pressure_compat.valid == 0 &&
              gpview.constrained_order[ii] >= 2) {
            Array1D<uint8_t, 0, N_InterP - 1> compatibility_mask{};
            Array1D<Real, 0, AMREX_SPACEDIM - 1>
                compatibility_point{};
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              compatibility_point(d) = gpview.imp_xyz[ii](0, d);
            }
            for (int support = 0; support < N_InterP; ++support) {
              const bool used =
                  gpview.constrained_dirichlet_weights[ii](support) !=
                      Real(0.0) ||
                  gpview.constrained_neumann_weights[ii](support) !=
                      Real(0.0);
              compatibility_mask(support) = used ? uint8_t(1) : uint8_t(0);
            }
            pressure_compat =
                copy->template
                    ibm_viscous_pressure_compatibility_from_support<
                        N_InterP>(
                        gpview.constrained_support_ijk[ii],
                        compatibility_mask, compatibility_point, prims, nvec,
                        prob_lo_lev, dx_lev, cls, cartesian_geometry,
                        support_visibility_condition_max);
          }
        }
      }
      if constexpr (sharp_feature_pressure_limiter) {
        const auto& edge = se_ptr[elem_idx];
        Real coordinate = Real(0.5) * edge.measure;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          Real centroid_world = xform.d[d];
          for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            centroid_world += xform.R[d][j] * edge.centroid[j];
          }
          coordinate +=
              (gpview.ib_xyz[ii](d) - centroid_world) *
              localframe.tangent1[d];
        }
        const Real h = amrex::min(dx_lev[0], dx_lev[1]);
        pressure_compat.feature_pressure_fallback =
            sharpFeaturePressureFallback(edge, coordinate, h);
      }
      ibm_detail::dispatch_compute_surfIB<eorder_tparm, wallmodel>(
          gpview.ib_xyz[ii], nvec, t1vec, t2vec, gpview.disIM[ii], n_valid_bc,
          primsNormal, type_solid_bc, cls, pressure_compat);

      // 6) Extrapolate from surface/image points back to ghost point
      copy->template extrapolate<eorder_tparm>(
          primsNormal, gpview.imp_ninterp[ii], gpview.disGP[ii], gpview.disIM[ii]);

      // 7) Transform ghost-point velocity back to global coordinates
      int idx = 0;
      copy->template local2global<eorder_tparm>(idx, primsNormal, nvec, t1vec, t2vec);

      // 8) Extract primitive variables at ghost point (slot 0)
      Real P = primsNormal(0, cls_t::QPRES);
      Real T = primsNormal(0, cls_t::QT);

      Real Y[NUM_SPECIES] = { Real(0.0) };
#if NUM_SPECIES > 1
      for (int n = 0; n < NUM_SPECIES; ++n) {
          Y[n] = primsNormal(0, cls_t::QFS + n);
      }
#endif

      Real ux = primsNormal(0, cls_t::QU);
      Real uy = primsNormal(0, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
      Real uz = primsNormal(0, cls_t::QW);
#else
      Real uz = Real(0.0);
#endif


      if constexpr (use_bi_constrained_shared_gp) {
        // Evaluate one BI-constrained extension at this real ghost-cell centre.
        // Pressure uses the selected normal-momentum compatibility condition.
        // Euler-slip thermodynamics use the entropy proxy and EOS; the older
        // adiabatic wall retains homogeneous temperature Neumann data.
        // No-penetration is imposed as a Dirichlet condition on local normal
        // velocity. Slip tangential velocity uses homogeneous Neumann data,
        // while a no-slip velocity uses the same zero Dirichlet functional.
        P = Real(0.0);
        T = Real(0.0);
        Real velocity_normal = Real(0.0);
        Real normal_momentum = Real(0.0);
        Real velocity_tangent1 = Real(0.0);
        Real velocity_tangent2 = Real(0.0);
        Real entropy_proxy = Real(0.0);
        Real entropy_proxy_linear = Real(0.0);
        Real entropy_proxy_reference = Real(0.0);
        Real entropy_support_min = std::numeric_limits<Real>::infinity();
        Real entropy_support_max = -std::numeric_limits<Real>::infinity();
        Real normal_support_min = std::numeric_limits<Real>::infinity();
        Real normal_support_max = -std::numeric_limits<Real>::infinity();
        bool entropy_support_valid =
            !one_sided_entropy_jet ||
            gpview.constrained_entropy_jet_order[ii] >= 1;
        const bool normal_momentum_support_valid =
            !cell_average_bic_normal_momentum ||
            gpview.constrained_fv_dirichlet_order[ii] >= 1;
        Real fluid_trace_density = Real(0.0);
        Real fluid_trace_tangent1 = Real(0.0);
        Real fluid_trace_tangent2 = Real(0.0);
        Real fluid_trace_pressure = Real(0.0);
        Real nearest_support_pressure = Real(0.0);
        Real nearest_support_temperature = Real(0.0);
        Real nearest_support_normal = Real(0.0);
        Real nearest_support_tangent1 = Real(0.0);
        Real nearest_support_tangent2 = Real(0.0);
        Real nearest_support_distance2 =
            std::numeric_limits<Real>::infinity();
#if CLIP_TEMPERATURE_MIN
        constexpr Real minimum_admissible_temperature =
            CNSConstants::min_temp();
#else
        constexpr Real minimum_admissible_temperature = Real(0.0);
#endif

        for (int support = 0; support < N_InterP; ++support) {
          const int si = gpview.constrained_support_ijk[ii](support, 0);
          const int sj = gpview.constrained_support_ijk[ii](support, 1);
#if (AMREX_SPACEDIM == 3)
          const int sk = gpview.constrained_support_ijk[ii](support, 2);
#else
          const int sk = 0;
#endif
          const Real dirichlet_weight =
              rz_point_target >= 0
                  ? gpview.rz_bic_point_functionals[ii]
                        .dirichlet_weights(rz_point_target, support)
                  : gpview.constrained_dirichlet_weights[ii](support);
          const Real neumann_weight =
              rz_point_target >= 0
                  ? gpview.rz_bic_point_functionals[ii]
                        .neumann_weights(rz_point_target, support)
                  : gpview.constrained_neumann_weights[ii](support);
          Real fv_dirichlet_weight = Real(0.0);
          if constexpr (cell_average_bic_normal_momentum) {
            fv_dirichlet_weight =
                rz_point_target >= 0
                    ? gpview.rz_bic_point_functionals[ii]
                          .fv_dirichlet_weights(rz_point_target, support)
                    : gpview.constrained_fv_dirichlet_weights[ii](support);
          }
          Real entropy_weight = neumann_weight;
          Real entropy_linear_weight = neumann_weight;
          if constexpr (one_sided_entropy_jet) {
            entropy_weight = rz_point_target >= 0
                ? gpview.rz_bic_point_functionals[ii]
                      .entropy_jet_weights(rz_point_target, support)
                : gpview.constrained_entropy_jet_weights[ii](support);
            entropy_linear_weight = rz_point_target >= 0
                ? gpview.rz_bic_point_functionals[ii]
                      .entropy_jet_linear_weights(rz_point_target, support)
                : gpview.constrained_entropy_jet_linear_weights[ii](support);
          }
          Real fluid_trace_weight = Real(0.0);
          if constexpr (uses_euler_slip_curvature_pressure_closure) {
            fluid_trace_weight =
                gpview.constrained_fluid_trace_weights[ii](support);
          }
          const bool support_is_used =
              dirichlet_weight != Real(0.0) ||
              neumann_weight != Real(0.0) ||
              fv_dirichlet_weight != Real(0.0) ||
              entropy_weight != Real(0.0) ||
              entropy_linear_weight != Real(0.0) ||
              fluid_trace_weight != Real(0.0);
          if (!support_is_used) {
            continue;
          }

          const Real support_pressure =
              prims(si, sj, sk, cls_t::QPRES);
          const Real support_density =
              prims(si, sj, sk, cls_t::QRHO);
          P += neumann_weight * support_pressure;
          if constexpr (entropy_wall_extension) {
            if (entropy_weight != Real(0.0) ||
                entropy_linear_weight != Real(0.0) ||
                neumann_weight != Real(0.0)) {
              const bool entropy_state_valid =
                  support_pressure > CNSConstants::min_press() &&
                  support_density > Real(0.0) &&
                  amrex::Math::isfinite(
                      support_pressure + support_density);
              entropy_support_valid =
                  entropy_support_valid && entropy_state_valid;
              if (entropy_state_valid) {
                const Real support_entropy =
                    std::log(support_pressure) -
                    cls->gamma * std::log(support_density);
                entropy_proxy += entropy_weight * support_entropy;
                entropy_proxy_linear +=
                    entropy_linear_weight * support_entropy;
                entropy_proxy_reference +=
                    neumann_weight * support_entropy;
                entropy_support_min =
                    amrex::min(entropy_support_min, support_entropy);
                entropy_support_max =
                    amrex::max(entropy_support_max, support_entropy);
              }
            }
          } else if constexpr (adiabatic_wall) {
            T += neumann_weight * prims(si, sj, sk, cls_t::QT);
          } else {
            T += dirichlet_weight * prims(si, sj, sk, cls_t::QT);
          }

          const Real support_ux = prims(si, sj, sk, cls_t::QU);
          const Real support_uy = prims(si, sj, sk, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
          const Real support_uz = prims(si, sj, sk, cls_t::QW);
#endif
          const Real support_normal =
              support_ux * nvec(0) + support_uy * nvec(1)
#if (AMREX_SPACEDIM == 3)
              + support_uz * nvec(2)
#endif
              ;
          const Real support_tangent1 =
              support_ux * t1vec(0) + support_uy * t1vec(1)
#if (AMREX_SPACEDIM == 3)
              + support_uz * t1vec(2)
#endif
              ;
          const Real support_tangent2 =
#if (AMREX_SPACEDIM == 3)
              support_ux * t2vec(0) + support_uy * t2vec(1) +
              support_uz * t2vec(2);
#else
              Real(0.0);
#endif

          if (amrex::Math::isfinite(support_normal)) {
            normal_support_min =
                amrex::min(normal_support_min, support_normal);
            normal_support_max =
                amrex::max(normal_support_max, support_normal);
          }

          if constexpr (entropy_wall_extension) {
            const Real support_temperature =
                prims(si, sj, sk, cls_t::QT);
            const Real support_x =
                prob_lo_lev[0] + (Real(si) + Real(0.5)) * dx_lev[0];
            const Real support_y =
                prob_lo_lev[1] + (Real(sj) + Real(0.5)) * dx_lev[1];
            Real distance2 =
                (support_x - gpview.ib_xyz[ii](0)) *
                    (support_x - gpview.ib_xyz[ii](0)) +
                (support_y - gpview.ib_xyz[ii](1)) *
                    (support_y - gpview.ib_xyz[ii](1));
#if (AMREX_SPACEDIM == 3)
            const Real support_z =
                prob_lo_lev[2] + (Real(sk) + Real(0.5)) * dx_lev[2];
            distance2 +=
                (support_z - gpview.ib_xyz[ii](2)) *
                (support_z - gpview.ib_xyz[ii](2));
#endif
            const bool support_state_admissible =
                support_pressure > CNSConstants::min_press() &&
                support_temperature > minimum_admissible_temperature &&
                amrex::Math::isfinite(
                    support_pressure + support_temperature + support_normal +
                    support_tangent1 + support_tangent2);
            if (support_state_admissible &&
                distance2 < nearest_support_distance2) {
              nearest_support_distance2 = distance2;
              nearest_support_pressure = support_pressure;
              nearest_support_temperature = support_temperature;
              nearest_support_normal = support_normal;
              nearest_support_tangent1 = support_tangent1;
              nearest_support_tangent2 = support_tangent2;
            }
          }

          velocity_normal += dirichlet_weight * support_normal;
          if constexpr (cell_average_bic_normal_momentum) {
            // Cartesian mode consumes Q(Ubar); R-Z mode consumes the recovered
            // centre-state momentum and applies an annular target functional.
            // The stationary-wall boundary contribution is zero in either
            // representation.
            normal_momentum += fv_dirichlet_weight * support_density *
                               support_normal;
          }
          if constexpr (stationary_slip_wall) {
            velocity_tangent1 += neumann_weight * support_tangent1;
            velocity_tangent2 += neumann_weight * support_tangent2;
          } else {
            velocity_tangent1 += dirichlet_weight * support_tangent1;
            velocity_tangent2 += dirichlet_weight * support_tangent2;
          }
          if constexpr (uses_euler_slip_curvature_pressure_closure) {
            fluid_trace_density +=
                fluid_trace_weight * prims(si, sj, sk, cls_t::QRHO);
            fluid_trace_pressure +=
                fluid_trace_weight * prims(si, sj, sk, cls_t::QPRES);
            fluid_trace_tangent1 +=
                fluid_trace_weight * support_tangent1;
            fluid_trace_tangent2 +=
                fluid_trace_weight * support_tangent2;
          }
        }

        if constexpr (isothermal_wall) {
          const Real boundary_weight = rz_point_target >= 0
              ? gpview.rz_bic_point_functionals[ii]
                    .dirichlet_boundary_weight(rz_point_target)
              : gpview.constrained_dirichlet_boundary_weight[ii];
          T += boundary_weight * wall_temperature;
        }
        if constexpr (uses_prescribed_gradient_pressure_closure) {
          const Real dpdn = param::pressure_normal_derivative(
              gpview.ib_xyz[ii], nvec);
          const Real gradient_weight = rz_point_target >= 0
              ? gpview.rz_bic_point_functionals[ii]
                    .neumann_gradient_weight(rz_point_target)
              : gpview.constrained_neumann_gradient_weight[ii];
          P += gradient_weight * dpdn;
        }
        if constexpr (use_ns_noslip_pressure_compatibility) {
          // For a stationary no-slip wall, the normal momentum equation gives
          // dp/dn = n dot div(tau).  The compatibility value is reconstructed
          // from the current-stage fluid trace above; the same BI-CWLS
          // Neumann functional then evaluates pressure at this unique shared
          // ghost-cell centre.  If the quadratic compatibility reconstruction
          // is unavailable, retaining P is the explicit homogeneous-Neumann
          // fallback used by the legacy wall closure.
          if (pressure_compat.valid != 0 &&
              amrex::Math::isfinite(pressure_compat.dpdn)) {
            const Real gradient_weight = rz_point_target >= 0
                ? gpview.rz_bic_point_functionals[ii]
                      .neumann_gradient_weight(rz_point_target)
                : gpview.constrained_neumann_gradient_weight[ii];
            P += gradient_weight * pressure_compat.dpdn;
          }
        }
        if constexpr (entropy_wall_extension) {
          // The pressure-compatible candidate is auxiliary: if it is not
          // thermodynamically admissible, reject the complete ghost state and
          // use an admissible lower-order extension. This local order
          // reduction is inactive for a smooth resolved wall and avoids
          // passing a mixed pressure/entropy state to the EOS.
          int pressure_fallback = 5;
          Real pressure_candidate = P;
          if constexpr (uses_euler_slip_curvature_pressure_closure) {
            const auto shape = param::wall_shape_operator(
                gpview.ib_xyz[ii], nvec, t1vec, t2vec);
            const bool trace_valid =
                gpview.constrained_fluid_trace_order[ii] >= 1 &&
                fluid_trace_density > Real(0.0) &&
                amrex::Math::isfinite(
                    fluid_trace_density + fluid_trace_tangent1 +
                    fluid_trace_tangent2);
            const bool shape_valid =
                shape.valid != 0 && amrex::Math::isfinite(shape.k11) &&
                amrex::Math::isfinite(shape.k12) &&
                amrex::Math::isfinite(shape.k22);
            if (trace_valid && shape_valid) {
              const Real curvature_acceleration =
                  shape.k11 * fluid_trace_tangent1 * fluid_trace_tangent1 +
                  Real(2.0) * shape.k12 * fluid_trace_tangent1 *
                      fluid_trace_tangent2 +
                  shape.k22 * fluid_trace_tangent2 * fluid_trace_tangent2;
              const Real dpdn = fluid_trace_density * curvature_acceleration;
              if (amrex::Math::isfinite(dpdn)) {
                pressure_candidate = P +
                    (rz_point_target >= 0
                         ? gpview.rz_bic_point_functionals[ii]
                               .neumann_gradient_weight(rz_point_target)
                         : gpview.constrained_neumann_gradient_weight[ii]) *
                        dpdn;
                pressure_fallback = 0;
              }
            }
          } else if constexpr (uses_prescribed_gradient_pressure_closure) {
            // P already contains the prescribed Neumann contribution added
            // above. Re-evaluate only to validate the problem callback and to
            // expose the same diagnostic quantity as the curvature path.
            const Real dpdn = param::pressure_normal_derivative(
                gpview.ib_xyz[ii], nvec);
            if (amrex::Math::isfinite(dpdn) &&
                amrex::Math::isfinite(pressure_candidate)) {
              pressure_fallback = 0;
            }
          }
          Real temperature_candidate = T;
          Real density_candidate = Real(-1.0);
          const Real velocity_normal_reference = velocity_normal;
          Real velocity_normal_candidate = velocity_normal;
          Real velocity_normal_cell_average = velocity_normal;
          Real entropy_hierarchy_ratio = Real(0.0);
          Real entropy_variation_ratio = Real(0.0);
          Real entropy_variation_fraction = Real(1.0);
          Real normal_hierarchy_ratio = Real(0.0);
          Real normal_variation_ratio = Real(0.0);
          Real normal_variation_fraction = Real(1.0);
          const Real entropy_proxy_high = entropy_proxy;
          Real entropy_high_order_fraction = Real(1.0);
          Real normal_high_order_fraction = Real(1.0);
          int auxiliary_fallback = 0;
          bool entropy_hierarchy_consistent = true;
          if constexpr (one_sided_entropy_jet) {
            const bool entropy_range_valid =
                amrex::Math::isfinite(
                    entropy_support_min + entropy_support_max +
                    entropy_proxy_high + entropy_proxy_linear) &&
                entropy_support_max >= entropy_support_min;
            if (entropy_range_valid) {
              const Real entropy_magnitude = amrex::max(
                  Real(1.0),
                  amrex::max(
                      amrex::max(amrex::Math::abs(entropy_support_min),
                                 amrex::Math::abs(entropy_support_max)),
                      amrex::max(amrex::Math::abs(entropy_proxy_high),
                                 amrex::Math::abs(entropy_proxy_linear))));
              const Real entropy_hierarchy_difference =
                  amrex::Math::abs(
                      entropy_proxy_high - entropy_proxy_linear);
              const Real entropy_support_range =
                  entropy_support_max - entropy_support_min;
              // A raw |P2-P1|/range sensor is singular at a smooth extremum:
              // both numerator and range are O(h^2), so it remains active
              // under refinement.  The square-root regularisation preserves
              // the local-range shock sensor while making the ratio O(h) at
              // both ordinary smooth points and smooth extrema.
              const Real entropy_regularized_range =
                  entropy_support_range +
                  std::sqrt(
                      entropy_hierarchy_difference * entropy_magnitude);
              const Real entropy_roundoff_scale =
                  Real(64.0) *
                  std::sqrt(std::numeric_limits<Real>::epsilon()) *
                  entropy_magnitude;
              entropy_hierarchy_ratio = entropy_hierarchy_difference /
                  amrex::max(
                      entropy_regularized_range, entropy_roundoff_scale);
              const Real scaled_ratio = amrex::min(
                  entropy_hierarchy_ratio /
                      bic_auxiliary_hierarchy_tolerance,
                  Real(1.0e6));
              const Real scaled_ratio_sq = scaled_ratio * scaled_ratio;
              entropy_high_order_fraction = Real(1.0) /
                  (Real(1.0) + scaled_ratio_sq * scaled_ratio_sq);
              const Real entropy_hierarchy_value = entropy_proxy_linear +
                  entropy_high_order_fraction *
                      (entropy_proxy_high - entropy_proxy_linear);
              // P1 and P2 may agree while both extrapolate an unresolved
              // entropy jump.  Delta sigma distinguishes that case: it
              // vanishes linearly under smooth-grid refinement but remains
              // finite across a shock.  Reduce only sigma_g toward the
              // homogeneous-Neumann BIC reference; never blend the complete
              // primitive state with this scalar sensor.
              entropy_variation_ratio =
                  (entropy_support_max - entropy_support_min) /
                  bic_entropy_jet_variation_tolerance;
              const Real capped_variation_ratio =
                  amrex::min(entropy_variation_ratio, Real(1.0e6));
              const Real variation_ratio_sq =
                  capped_variation_ratio * capped_variation_ratio;
              entropy_variation_fraction = Real(1.0) /
                  (Real(1.0) + variation_ratio_sq * variation_ratio_sq);
              entropy_proxy = entropy_proxy_reference +
                  entropy_variation_fraction *
                      (entropy_hierarchy_value - entropy_proxy_reference);
              entropy_high_order_fraction *= entropy_variation_fraction;
              if (entropy_high_order_fraction < Real(0.5)) {
                auxiliary_fallback = 1;
              }
              entropy_hierarchy_consistent =
                  amrex::Math::isfinite(entropy_proxy);
            } else {
              entropy_hierarchy_ratio =
                  std::numeric_limits<Real>::infinity();
              entropy_hierarchy_consistent = false;
            }
          }
          if constexpr (entropy_wall_extension) {
            if (entropy_support_valid && cls->gamma > Real(1.0) &&
                cls->Rspec > Real(0.0) &&
                pressure_candidate > CNSConstants::min_press()) {
              density_candidate = std::exp(
                  (std::log(pressure_candidate) - entropy_proxy) /
                  cls->gamma);
              temperature_candidate =
                  pressure_candidate / (density_candidate * cls->Rspec);
            } else {
              temperature_candidate = Real(-1.0);
            }
          }
          bool normal_candidate_valid = true;
          bool normal_hierarchy_consistent = true;
          if constexpr (cell_average_bic_normal_momentum) {
            normal_candidate_valid = normal_momentum_support_valid &&
                density_candidate > Real(0.0) &&
                amrex::Math::isfinite(density_candidate + normal_momentum);
            if (normal_candidate_valid) {
              velocity_normal_candidate = normal_momentum / density_candidate;
              velocity_normal_cell_average = velocity_normal_candidate;
            }
            const bool normal_range_valid =
                normal_candidate_valid &&
                amrex::Math::isfinite(
                    normal_support_min + normal_support_max +
                    velocity_normal + velocity_normal_candidate) &&
                normal_support_max >= normal_support_min;
            if (normal_range_valid) {
              Real characteristic_speed = amrex::max(
                  Real(1.0),
                  amrex::max(
                      amrex::Math::abs(nearest_support_normal),
                      amrex::max(
                          amrex::Math::abs(nearest_support_tangent1),
                          amrex::Math::abs(nearest_support_tangent2))));
              if (nearest_support_temperature > Real(0.0) &&
                  cls->gamma > Real(1.0) && cls->Rspec > Real(0.0)) {
                const Real sound_speed = std::sqrt(
                    cls->gamma * cls->Rspec *
                    nearest_support_temperature);
                if (amrex::Math::isfinite(sound_speed)) {
                  characteristic_speed =
                      amrex::max(characteristic_speed, sound_speed);
                }
              }
              const Real normal_hierarchy_difference =
                  amrex::Math::abs(
                      velocity_normal_candidate - velocity_normal);
              const Real normal_support_range =
                  normal_support_max - normal_support_min;
              const Real normal_regularized_range =
                  normal_support_range +
                  std::sqrt(
                      normal_hierarchy_difference * characteristic_speed);
              const Real normal_roundoff_scale =
                  Real(64.0) *
                  std::sqrt(std::numeric_limits<Real>::epsilon()) *
                  characteristic_speed;
              normal_hierarchy_ratio = normal_hierarchy_difference /
                  amrex::max(
                      normal_regularized_range, normal_roundoff_scale);
              const Real scaled_ratio = amrex::min(
                  normal_hierarchy_ratio /
                      bic_auxiliary_hierarchy_tolerance,
                  Real(1.0e6));
              const Real scaled_ratio_sq = scaled_ratio * scaled_ratio;
              normal_high_order_fraction = Real(1.0) /
                  (Real(1.0) + scaled_ratio_sq * scaled_ratio_sq);
              velocity_normal_candidate = velocity_normal +
                  normal_high_order_fraction *
                      (velocity_normal_candidate - velocity_normal);
              normal_variation_ratio =
                  (normal_support_max - normal_support_min) /
                  (bic_normal_momentum_variation_tolerance *
                   characteristic_speed);
              const Real capped_variation_ratio =
                  amrex::min(normal_variation_ratio, Real(1.0e6));
              const Real variation_ratio_sq =
                  capped_variation_ratio * capped_variation_ratio;
              normal_variation_fraction = Real(1.0) /
                  (Real(1.0) + variation_ratio_sq * variation_ratio_sq);
              velocity_normal_candidate = velocity_normal +
                  normal_variation_fraction *
                      (velocity_normal_candidate - velocity_normal);
              normal_high_order_fraction *= normal_variation_fraction;
              if (normal_high_order_fraction < Real(0.5)) {
                auxiliary_fallback = 1;
              }
              normal_hierarchy_consistent =
                  amrex::Math::isfinite(velocity_normal_candidate);
            } else {
              normal_hierarchy_ratio =
                  std::numeric_limits<Real>::infinity();
              normal_hierarchy_consistent = false;
            }
          }
          const bool thermodynamic_candidate_valid =
              !entropy_wall_extension ||
              (density_candidate > Real(0.0) &&
               amrex::Math::isfinite(density_candidate));
          bool candidate_state_admissible =
              pressure_candidate > CNSConstants::min_press() &&
              temperature_candidate > minimum_admissible_temperature &&
              thermodynamic_candidate_valid &&
              normal_candidate_valid &&
              entropy_hierarchy_consistent &&
              normal_hierarchy_consistent &&
              amrex::Math::isfinite(
                  pressure_candidate + temperature_candidate +
                  velocity_normal_candidate + velocity_tangent1 +
                  velocity_tangent2);
          // If the high-order state is inadmissible even though the hierarchy
          // remains smooth, try one coherent lower-order BIC state before the
          // nearest-support fallback. This is order reduction, not clipping.
          if constexpr (cell_average_bic_normal_momentum ||
                        one_sided_entropy_jet) {
            if (!(candidate_state_admissible && pressure_fallback == 0) &&
                pressure_fallback == 0 &&
                pressure_candidate > CNSConstants::min_press()) {
              Real reduced_density_candidate = Real(-1.0);
              Real reduced_temperature_candidate = T;
              if constexpr (entropy_wall_extension) {
                if (entropy_support_valid && cls->gamma > Real(1.0) &&
                    cls->Rspec > Real(0.0) &&
                    amrex::Math::isfinite(entropy_proxy_reference)) {
                  reduced_density_candidate = std::exp(
                      (std::log(pressure_candidate) -
                       entropy_proxy_reference) /
                      cls->gamma);
                  reduced_temperature_candidate =
                      pressure_candidate /
                      (reduced_density_candidate * cls->Rspec);
                } else {
                  reduced_temperature_candidate = Real(-1.0);
                }
              }
              const bool reduced_thermodynamics_admissible =
                  !entropy_wall_extension ||
                  (reduced_density_candidate > Real(0.0) &&
                   amrex::Math::isfinite(reduced_density_candidate));
              const bool reduced_state_admissible =
                  reduced_thermodynamics_admissible &&
                  reduced_temperature_candidate >
                      minimum_admissible_temperature &&
                  amrex::Math::isfinite(
                      pressure_candidate + reduced_temperature_candidate +
                      velocity_normal + velocity_tangent1 +
                      velocity_tangent2);
              if (reduced_state_admissible) {
                density_candidate = reduced_density_candidate;
                temperature_candidate = reduced_temperature_candidate;
                velocity_normal_candidate = velocity_normal;
                candidate_state_admissible = true;
                auxiliary_fallback = 1;
              }
            }
          }
          const bool nearest_state_admissible =
              nearest_support_pressure > CNSConstants::min_press() &&
              nearest_support_temperature > minimum_admissible_temperature &&
              amrex::Math::isfinite(
                  nearest_support_pressure + nearest_support_temperature +
                  nearest_support_normal + nearest_support_tangent1 +
                  nearest_support_tangent2);
          if (candidate_state_admissible && pressure_fallback == 0) {
            P = pressure_candidate;
            T = temperature_candidate;
            velocity_normal = velocity_normal_candidate;
          } else if (nearest_state_admissible) {
            // A rejected pressure candidate identifies a non-smooth local
            // extension. Reduce the complete primitive state together rather
            // than combining a low-order pressure with high-order temperature
            // and velocity. The normal component is reflected across the
            // stationary slip wall; thermodynamic and tangential components
            // are copied from the nearest selected visible fluid support.
            P = nearest_support_pressure;
            T = nearest_support_temperature;
            velocity_normal = -nearest_support_normal;
            velocity_tangent1 = nearest_support_tangent1;
            velocity_tangent2 = nearest_support_tangent2;
            pressure_fallback = 1;
          } else if (P > CNSConstants::min_press()) {
            Real reduced_temperature = T;
            Real reduced_density = Real(-1.0);
            Real reduced_velocity_normal = velocity_normal;
            if constexpr (entropy_wall_extension) {
              if (entropy_support_valid && cls->gamma > Real(1.0) &&
                  cls->Rspec > Real(0.0)) {
                reduced_density = std::exp(
                    (std::log(P) - entropy_proxy) / cls->gamma);
                reduced_temperature = P / (reduced_density * cls->Rspec);
              } else {
                reduced_temperature = Real(-1.0);
              }
            }
            bool reduced_normal_valid = true;
            if constexpr (cell_average_bic_normal_momentum) {
              reduced_normal_valid = normal_momentum_support_valid &&
                  reduced_density > Real(0.0) &&
                  amrex::Math::isfinite(reduced_density + normal_momentum);
              if (reduced_normal_valid) {
                reduced_velocity_normal = normal_momentum / reduced_density;
              }
            }
            if (reduced_temperature > minimum_admissible_temperature &&
                reduced_normal_valid &&
                amrex::Math::isfinite(
                    P + reduced_temperature + reduced_velocity_normal +
                    velocity_tangent1 + velocity_tangent2)) {
              T = reduced_temperature;
              velocity_normal = reduced_velocity_normal;
              pressure_fallback = 2;
            } else if (nearest_state_admissible) {
              P = nearest_support_pressure;
              T = nearest_support_temperature;
              velocity_normal = -nearest_support_normal;
              velocity_tangent1 = nearest_support_tangent1;
              velocity_tangent2 = nearest_support_tangent2;
              pressure_fallback = 1;
            } else {
              P = pressure_candidate;
              T = temperature_candidate;
              pressure_fallback = 4;
            }
          } else if (
              gpview.constrained_fluid_trace_order[ii] >= 1 &&
              fluid_trace_pressure > CNSConstants::min_press() &&
              amrex::Math::isfinite(fluid_trace_pressure)) {
            P = fluid_trace_pressure;
            pressure_fallback = 3;
          } else {
            // The active-fluid solution itself is already inadmissible or no
            // valid support remains. Leave the candidate for the established
            // EOS floor and expose this exceptional path in diagnostics.
            P = pressure_candidate;
            pressure_fallback = 4;
          }
        }
        ux = velocity_normal * nvec(0) +
             velocity_tangent1 * t1vec(0) +
             velocity_tangent2 * t2vec(0);
        uy = velocity_normal * nvec(1) +
             velocity_tangent1 * t1vec(1) +
             velocity_tangent2 * t2vec(1);
#if (AMREX_SPACEDIM == 3)
        uz = velocity_normal * nvec(2) +
             velocity_tangent1 * t1vec(2) +
             velocity_tangent2 * t2vec(2);
#else
        uz = Real(0.0);
#endif
      }


      // 9) Enforce thermodynamic consistency
      Real Q[cls_t::NPRIM];
      cls->ensurePTYfillq(P, T, Y, ux, uy, uz, Q);


      // 10) Write ghost-cell primitive variables back into prims
      int i = gpview.gp_ijk[ii](0);
      int j = gpview.gp_ijk[ii](1);
#if (AMREX_SPACEDIM == 3)
      int k = gpview.gp_ijk[ii](2);
#else
      int k = 0;
#endif

      for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(i,j,k,n) = Q[n];
      }
    }); // end ParallelFor over all ghost points
    }

    // GP reconstruction writes the valid owner cell only.  A neighboring FAB
    // can read that cell through its private ghost copy in WENO/TENO or viscous
    // stencils, so publish the reconstructed values before returning.  This is
    // intentionally unconditional: ranks with zero local GPs must still
    // participate in the same-level exchange initiated by ranks that own GPs.
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
    // FillBoundary publishes the unique shared GP across FAB/MPI ownership,
    // but it does not impose a physical boundary condition.  In the generic
    // point-GP path, refresh negative-radius coordinate ghosts only after the
    // positive-radius GP has reached its final value.  Annular quadrature
    // target passes are closed by computeAllGPsRZAnnular after all targets are
    // accumulated, so do not alter their intermediate state here.
    if (rz_point_target < 0) {
      fillRZAxisPrimitiveParity(prims_mf, lev);
    }
  }


  /**
   * \brief Repair conservative state in cells exposed by moving geometry.
   *
   * A level-wide Jacobi front propagates data from cells that are fluid in
   * both topologies into cells that changed solid -> fluid.  Tags and state
   * ghosts are exchanged after every iteration, so the result is independent
   * of GPU scheduling, FAB decomposition and MPI rank boundaries.  Only a
   * one-cell halo is accessed; there are no unchecked expanding-ring reads.
   *
   * If an exposed component is disconnected from every old-fluid donor, the
   * state is physically undefined.  Abort instead of inventing a freestream
   * value or retaining stale solid data.
   */
  void fixExposedCells(const FabArray<BaseFab<uint8_t>>& old_markers,
                       MultiFab& state_mf,
                       int lev)
  {
    BL_PROFILE("IBM::fixExposedCells");

    auto& new_markers = *bmf_a[lev];
    AMREX_ALWAYS_ASSERT(old_markers.boxArray() == new_markers.boxArray());
    AMREX_ALWAYS_ASSERT(state_mf.boxArray() == new_markers.boxArray());
    AMREX_ALWAYS_ASSERT(state_mf.nGrowVect().allGE(IntVect(1)));

    constexpr int invalid = -1;
    constexpr int donor   = 0;
    constexpr int pending = 1;
    constexpr int fixed   = 2;

    // Two tag components are the old/new Jacobi buffers.  Ghosts start invalid
    // and become valid only through FillBoundary from an owned same-level cell.
    iMultiFab tags(new_markers.boxArray(), new_markers.DistributionMap(),
                   2, 1, MFInfo().SetArena(The_Async_Arena()));
    tags.setVal(invalid);

    Gpu::DeviceScalar<int> d_pending(0);
    int* const p_pending = d_pending.dataPtr();

    for (MFIter mfi(new_markers, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& old_mk = old_markers.const_array(mfi);
      auto const& new_mk = new_markers.const_array(mfi);
      auto const& tag = tags.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int t = invalid;
        if (new_mk(i,j,k,0) == 0) {
          if (old_mk(i,j,k,0) == 0) {
            t = donor;
          } else {
            t = pending;
            Gpu::Atomic::Add(p_pending, 1);
          }
        }
        tag(i,j,k,0) = t;
        tag(i,j,k,1) = t;
      });
    }

    int n_pending = d_pending.dataValue();
    ParallelDescriptor::ReduceIntSum(n_pending);
    if (n_pending == 0) return;

    const auto& periodicity = amr_p->Geom(lev).periodicity();
    state_mf.FillBoundary(periodicity);
    tags.FillBoundary(periodicity);

    int max_iters = 0;
    const Box& domain = amr_p->Geom(lev).Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      max_iters = amrex::max(max_iters, domain.length(d));
    }

    int previous_pending = n_pending;
    for (int iter = 0; iter < max_iters; ++iter) {
      const int told = iter & 1;
      const int tnew = 1 - told;
      const int zero = 0;
      Gpu::htod_memcpy(p_pending, &zero, sizeof(int));

      for (MFIter mfi(new_markers, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);
        auto const& state = state_mf.array(mfi);
        const int nc = cls_t::NCONS;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          const int old_tag = tag(i,j,k,told);
          if (old_tag != pending) {
            tag(i,j,k,tnew) = old_tag;
            return;
          }

          Real sum[cls_t::NCONS] = {};
          int count = 0;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            for (int s = -1; s <= 1; s += 2) {
              int ii = i;
              int jj = j;
              int kk = k;
              if (d == 0) ii += s;
              if (d == 1) jj += s;
#if (AMREX_SPACEDIM == 3)
              if (d == 2) kk += s;
#endif
              const int neighbour_tag = tag(ii,jj,kk,told);
              if (neighbour_tag != donor && neighbour_tag != fixed) continue;

              bool valid_state = state(ii,jj,kk,cls_t::URHO) > Real(0.0)
                              && state(ii,jj,kk,cls_t::UET)  > Real(0.0);
              for (int n = 0; n < nc; ++n) {
                valid_state = valid_state && std::isfinite(state(ii,jj,kk,n));
              }
              if (!valid_state) continue;

              for (int n = 0; n < nc; ++n) {
                sum[n] += state(ii,jj,kk,n);
              }
              ++count;
            }
          }

          if (count > 0) {
            const Real inv = Real(1.0) / Real(count);
            for (int n = 0; n < nc; ++n) {
              state(i,j,k,n) = sum[n] * inv;
            }
            tag(i,j,k,tnew) = fixed;
          } else {
            tag(i,j,k,tnew) = pending;
            Gpu::Atomic::Add(p_pending, 1);
          }
        });
      }

      // Publish newly fixed states and tags before the next Jacobi iteration.
      state_mf.FillBoundary(periodicity);
      tags.FillBoundary(tnew, 1, periodicity);

      n_pending = d_pending.dataValue();
      ParallelDescriptor::ReduceIntSum(n_pending);
      if (n_pending == 0) return;

      if (n_pending >= previous_pending) {
        amrex::Abort(
            "fixExposedCells: exposed fluid region is disconnected from all "
            "old-fluid donors on level " + std::to_string(lev) +
            " (" + std::to_string(n_pending) + " cell(s) remain)");
      }
      previous_pending = n_pending;
    }

    amrex::Abort("fixExposedCells: Jacobi propagation exceeded the level "
                 "domain extent on level " + std::to_string(lev));
  }

  /**
    * \brief Compute surface indices and interpolation data for all faces/edges at given level
    *
    * Algorithm:
    *  1. Build spatial lookup (global_fab_idx -> local_fab_idx)
    *  2. For each face: compute mirror point, find owning FAB, compute interpolation weights
    *  3. Build CSR structure for GPU-friendly access
    *
    * \param lev AMR level
    */
  void computeSurfIndices(int lev)
  {
    BL_PROFILE("IBM::computeSurfIndices");
    const Real timing_total_start = beginPerformanceTiming();
    const bool timing_detail = detailedPerformanceTimingEnabled();
    Real timing_locate_seconds = Real(0.0);
    Real timing_surface_image_seconds = Real(0.0);

    int myrank = amrex::ParallelDescriptor::MyProc();
    amrex::Print()  << "Compute Surface Index at LEVEL " << lev  << std::endl;

    auto& mfab = *bmf_a[lev];

    const BoxArray& ba_global     = mfab.boxArray();
    const DistributionMapping& dm = mfab.DistributionMap();

    const int nfab_local  = mfab.local_size();
    const int nfab_global = ba_global.size();

    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto& domain = amr_p->Geom(lev).Domain();

#ifndef AMREX_USE_CGAL
    GpuArray<BVH4QueryView, MAX_NGEOM> visibility_queries{};
    GpuArray<RigidTransform, MAX_NGEOM> visibility_transforms{};
    GpuArray<int, MAX_NGEOM + 1> visibility_geometry_offsets{};
    if constexpr (support_visibility_mode > 0) {
      // Geometry setup populates managed BVH metadata asynchronously.  Finish
      // that initialization before copying offsets into host-side launch data;
      // otherwise WSL2 can fault on the managed-memory read below.
      Gpu::streamSynchronize();
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
          "ngeom exceeds MAX_NGEOM; increase MAX_NGEOM in ibm_containers.h");
      for (int geometry = 0; geometry < ngeom; ++geometry) {
        visibility_queries[geometry] =
            bvh_a[geometry].query_view(geom_a[geometry]);
        visibility_transforms[geometry] = transform_a[geometry];
        visibility_geometry_offsets[geometry] = geom_offsets[geometry];
      }
      visibility_geometry_offsets[ngeom] = geom_offsets[ngeom];
    }
#endif

    // ========================================================================
    // Phase 0: Initialize surfimp_soa / surfphys_soa
    // ========================================================================
    if (lev == amr_p->finestLevel()) {
      const bool need_cell_average_surface_data =
          nsSurfaceCellAverageRecoveryEnabled();
      if ((surfimp_soa.elemIdx.size() != ntotalfaces) ||
          (surfphys_soa.elemIdx.size() != ntotalfaces) ||
          (need_cell_average_surface_data &&
           surfimp_soa.imp_ipweights_cell_average.size() != ntotalfaces)) {
        surfimp_soa.resize(
            ntotalfaces, need_cell_average_surface_data);
        surfphys_soa.resize(ntotalfaces);
      }
      surfphys_soa.reset();
    }

    // ========================================================================
    // Phase 1: Build per-FAB lookup structures
    // ========================================================================
    Vector<Array4<uint8_t const>> fab_markers(nfab_local);
    Vector<Box> fab_bxg(nfab_local);
    Vector<int> global_to_local_fab(nfab_global, -1);

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
        int lidx = mfi.LocalIndex();
        int gidx = mfi.index();
        fab_markers[lidx] = mfab.const_array(mfi);
        fab_bxg[lidx] = interpolationReadBox(
            mfi, lev, surfaceInterpolationNghost(lev));
        global_to_local_fab[gidx] = lidx;
    }

    // ========================================================================
    // Phase 2: Fast locate — assign each face to a FAB (lightweight)
    //
    // This pass only does centroid→cell→FAB lookup and writes metadata.
    // The heavy image-point computation is deferred to Phase 3.
    // ========================================================================
    // Per-face FAB assignment: >=0 means locally owned, -1 means not local
    std::vector<int> surface_owner_fab(ntotalfaces, -1);

    int faces_found = 0;
    int faces_out_domain = 0;
    int faces_out_level  = 0;
    int faces_out_rank   = 0;
    int faces_in_finer   = 0;

    // Precompute inverse dx for faster centroid→cell conversion
    amrex::GpuArray<Real, AMREX_SPACEDIM> inv_dx;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        inv_dx[d] = Real(1.0) / dx_a[lev][d];
    }

    const Real timing_locate_start =
        timing_detail ? amrex::second() : Real(0.0);
    {
      std::vector<std::pair<int, Box>> isects;

      for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {

        // Skip faces already assigned to a finer level
        if (surfphys_soa.elemfound[f_idx] && surfphys_soa.lev[f_idx] > lev) {
            faces_in_finer++;
            continue;
        }

        // Convert body-frame centroid to world frame, then to cell index
        const SurfElem& surfelem = SurfElem_a[f_idx];
        // RZ fix: skip non-physical surface faces lying exactly on the symmetry
        // axis.  Where a closed axisymmetric contour runs along r=0 through the
        // solid interior (e.g. a long centre body / plenum), those faces have
        // zero circumference (2*pi*r -> 0) and their image point falls inside
        // the solid, so the stencil search always fails and the failed-search
        // host printf floods stdout -- serializing one rank and stalling
        // refinement of long on-axis geometries.  They carry no valid surface
        // data, so skip them.  (RZ-gated: no effect on Cartesian.)
        //
        // Exactly on-axis contour-closure segments have zero 2*pi*r surface
        // measure and carry no physical pressure-force contribution. Physical
        // curved elements with a small but non-zero centroid radius must remain;
        // a grid-dependent r<h/2 test would delete an O(h) polar cap as the
        // fluid grid is refined against fixed geometry.
        if (amr_p->Geom(lev).IsRZ()
            && std::abs(surfelem.centroid[0]) <= IBM_EPS::GEOM) {
            faces_out_domain++;
            continue;
        }
        int gIdx = getGeomIdx(f_idx);
        const auto& T = transform_a[gIdx];
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(surfelem.centroid[0], surfelem.centroid[1]);
#else
        Point c_body(surfelem.centroid[0], surfelem.centroid[1], surfelem.centroid[2]);
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = surfelem.centroid[dd];
#endif
        const Point raw_c_world = T.to_world(c_body);
        const auto raw_c_vec = make_vec<Real>(raw_c_world);
        Array1D<Real, 0, AMREX_SPACEDIM - 1> c_vec;
        apply_volume_surface_point_override(raw_c_vec, c_vec);
        const Point c_world = surface_point_from_array(c_vec);

        IntVect iv;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv[d] = static_cast<int>(std::floor(
                    (c_world[d] - prob_lo[d]) * inv_dx[d]));
        }

        if (!domain.contains(iv)) { faces_out_domain++; continue; }

        ba_global.intersections(Box(iv, iv), isects, true, IntVect::TheZeroVector());
        if (isects.empty()) { faces_out_level++; continue; }

        int gidx       = isects[0].first;
        int owner_rank = dm[gidx];
        int local_fab  = global_to_local_fab[gidx];

        surfphys_soa.ifab[f_idx] = local_fab;
        surfphys_soa.rank[f_idx] = owner_rank;
        surfphys_soa.lev[f_idx]  = lev;

        if (owner_rank != myrank) {
            surfphys_soa.elemfound[f_idx] = -1;
            faces_out_rank++;
            continue;
        }

        surfphys_soa.elemfound[f_idx] = 1;
        surface_owner_fab[f_idx] = local_fab;
        faces_found++;
      }
    } // isects freed here
    if (timing_detail) {
      timing_locate_seconds = amrex::second() - timing_locate_start;
    }

    // ========================================================================
    // Phase 3: Group locally-owned faces by FAB for cache locality
    //
    // Faces in the same FAB access the same ibMarkers array.  Processing
    // them together keeps marker data in cache and avoids thrashing when
    // calling search_optimal_image_point / computeIPweights.
    // ========================================================================
    std::vector<std::vector<int>> fab_faces(nfab_local);
    for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {
        if (surface_owner_fab[f_idx] >= 0) {
            fab_faces[surface_owner_fab[f_idx]].push_back(f_idx);
        }
    }

    // ========================================================================
    // Phase 4: Compute image points and interpolation weights per FAB
    //
    // Each face writes exclusively to its own f_idx slot in surfimp_soa /
    // surfphys_soa — no cross-face data dependencies.  The inner helpers
    // (search_optimal_image_point, computeIPweights) are stateless and
    // thread-safe, so the outer FAB loop is safe to parallelize with OpenMP.
    // ========================================================================
    int surf_stencil_fail = 0;
    int surf_place_fail = 0;
    int surf_fit_fail = 0;
    int surf_cell_average_fit_fail = 0;
    int surf_cell_average_fit_reduced = 0;
    int surf_tangent_shift = 0;
    int surf_visibility_supports = 0;
    int surf_visibility_occluded = 0;
    int surf_visibility_demoted = 0;
    int surf_visibility_first_hit = std::numeric_limits<int>::max();
    const Real timing_surface_image_start =
        timing_detail ? amrex::second() : Real(0.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+:surf_stencil_fail,surf_place_fail,surf_fit_fail,surf_cell_average_fit_fail,surf_cell_average_fit_reduced,surf_tangent_shift,surf_visibility_supports,surf_visibility_occluded,surf_visibility_demoted) reduction(min:surf_visibility_first_hit)
#endif
    for (int lfab = 0; lfab < nfab_local; ++lfab) {
      const auto& flist = fab_faces[lfab];
      if (flist.empty()) continue;

      // Pin the marker array and box for this FAB — all faces in flist
      // will read from this same data, maximizing cache reuse.
      auto const ibMarkers = fab_markers[lfab];
      auto const bxg       = fab_bxg[lfab];

      for (int f_idx : flist) {
        surfimp_soa.elemIdx[f_idx]  = f_idx;
        surfphys_soa.elemIdx[f_idx] = f_idx;

        // Rotate body-frame LocalFrame to world frame
        const LocalFrame& lf_body = LocalFrame_a[f_idx];
        const SurfElem&   surfelem = SurfElem_a[f_idx];
        int gIdx = getGeomIdx(f_idx);
        const auto& T = transform_a[gIdx];

        LocalFrame localframe;
        T.rotate_to_world(lf_body.normal,   localframe.normal);
        T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
        T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif

        // Transform body-frame centroid to world frame
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(surfelem.centroid[0], surfelem.centroid[1]);
#else
        Point c_body(surfelem.centroid[0], surfelem.centroid[1], surfelem.centroid[2]);
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = surfelem.centroid[dd];
#endif
        const Point raw_surf_centroid = T.to_world(c_body);
        const auto raw_surf_vec = make_vec<Real>(raw_surf_centroid);
        Array1D<Real, 0, AMREX_SPACEDIM - 1> surf_vec;
        apply_volume_surface_point_override(raw_surf_vec, surf_vec);
        const Point surf_centroid = surface_point_from_array(surf_vec);
        apply_volume_surface_frame_override(surf_vec, localframe);

        // Temporary storage for this face's image point data.  Keep it
        // deterministic so an invalid first IP can be diagnosed and rejected.
        Array2D<Real, 0, eorder_tparm_surf - 1, 0, IDIM> imp_xyz{};
        Array2D< int, 0, eorder_tparm_surf - 1, 0, IDIM> imp_ijk{};
        Array1D<Real, 0, eorder_tparm_surf - 1> disIM{};
        Array1D< int, 0, eorder_tparm_surf - 1> imp_ninterp{};
        Array1D< int, 0, eorder_tparm_surf - 1> imp_fit_order{};
        int first_place_status = 0;

        for (int jj = 0; jj < eorder_tparm_surf; jj++) {
          Point cp_start;
          if (jj == 0) {
            cp_start = surf_centroid;
          } else {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
            cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1));
#else
            cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1), imp_xyz(jj - 1, 2));
#endif
#else
            for (int d = 0; d < AMREX_SPACEDIM; ++d) cp_start[d] = imp_xyz(jj - 1, d);
#endif
          }

          if (jj == 0) {
            first_place_status =
                search_optimal_image_point<eorder_tparm_surf, iorder_tparm_surf>(
                                          cp_start, localframe,
                                          lev, prob_lo, dx_a[lev], di_a_surf[lev],
                                          bxg, ibMarkers,
                                          surfimp_soa, f_idx,
                                          imp_xyz, imp_ijk, disIM, imp_ninterp);
            if (first_place_status & 1) { ++surf_stencil_fail; }
            if (first_place_status & 2) { ++surf_place_fail; }
            if (first_place_status != 0) { surfphys_soa.elemfound[f_idx] = 0; }
          }
          else {
            search_image_point<eorder_tparm_surf, iorder_tparm_surf>(
                                      jj, cp_start, localframe,
                                      lev, prob_lo, dx_a[lev], di_a_surf[lev],
                                      bxg, ibMarkers,
                                      surfimp_soa, f_idx,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);
          }
        } // end loop on image points

        if constexpr (sharp_feature_tangential_reconstruction) {
          const Real h = amrex::min(dx_a[lev][0], dx_a[lev][1]);
          surf_tangent_shift +=
              applySharpFeatureTangentialStencil<eorder_tparm_surf,
                                                  iorder_tparm_surf>(
                  surfelem, Real(0.5) * surfelem.measure, h, localframe, lev,
                  prob_lo, dx_a[lev], bxg, ibMarkers, surfimp_soa, f_idx,
                  imp_xyz, imp_ijk, imp_ninterp);
        }
        surfimp_soa.imp_xyz[f_idx]     = imp_xyz;
        surfimp_soa.imp_ijk[f_idx]     = imp_ijk;
        surfimp_soa.disIM[f_idx]       = disIM;

        Array3D< int, 0, eorder_tparm_surf - 1, 0, N_InterP_surf - 1, 0, IDIM> imp_ip_ijk;
        Array2D<Real, 0, eorder_tparm_surf - 1, 0, N_InterP_surf - 1> imp_ipweights;
        Array2D<uint8_t, 0, eorder_tparm_surf - 1,
                0, N_InterP_surf - 1> visibility_mask{};
        Array2D<uint8_t, 0, eorder_tparm_surf - 1,
                0, N_InterP_surf - 1> weight_mask;
        Array2D<int, 0, eorder_tparm_surf - 1,
                0, N_InterP_surf - 1> first_hit_element{};
        Array1D<int, 0, eorder_tparm_surf - 1> candidate_count{};
        Array1D<int, 0, eorder_tparm_surf - 1> visible_count{};
        for (int image = 0; image < eorder_tparm_surf; ++image) {
          for (int support = 0; support < N_InterP_surf; ++support) {
            weight_mask(image, support) = uint8_t(1);
          }
        }

#ifndef AMREX_USE_CGAL
        if constexpr (support_visibility_mode > 0) {
          int checked = 0;
          int occluded = 0;
          int first_hit = -1;
          buildSupportVisibilityMask<eorder_tparm_surf,
                                     iorder_tparm_surf>(
              visibility_mask, first_hit_element, candidate_count,
              visible_count, imp_xyz, imp_ijk, imp_ninterp, prob_lo,
              dx_a[lev], ibMarkers, visibility_queries,
              visibility_transforms, visibility_geometry_offsets, ngeom,
              checked, occluded, first_hit);
          surf_visibility_supports += checked;
          surf_visibility_occluded += occluded;
          if (first_hit >= 0) {
            surf_visibility_first_hit =
                amrex::min(surf_visibility_first_hit, first_hit);
          }
          for (int image = 0; image < eorder_tparm_surf; ++image) {

            for (int support = 0; support < N_InterP_surf; ++support) {
              if constexpr (support_visibility_mode == 2) {
                weight_mask(image, support) =
                    visibility_mask(image, support);
              }
            }
            if constexpr (support_visibility_mode == 2) {
              imp_ninterp(image) = visible_count(image);
            }
          }
        }
#endif

        const auto geometric_ninterp = imp_ninterp;
        computeIPweights<eorder_tparm_surf, iorder_tparm_surf, SURFIMP>(
            imp_ipweights, imp_ip_ijk,
            imp_xyz, imp_ijk, imp_ninterp, imp_fit_order,
            weight_mask, prob_lo, dx_a[lev], ibMarkers,
            support_visibility_mode == 2
                ? support_visibility_condition_max
                : std::numeric_limits<Real>::infinity());
        if (first_place_status == 0 && imp_ninterp(0) == 0) {
          ++surf_fit_fail;
          surfphys_soa.elemfound[f_idx] = 0;
        }

        Array1D<int, 0, eorder_tparm_surf - 1>
            imp_ninterp_cell_average{};
        Array1D<int, 0, eorder_tparm_surf - 1>
            imp_fit_order_cell_average{};
        Array3D<int, 0, eorder_tparm_surf - 1,
                0, N_InterP_surf - 1, 0, IDIM>
            imp_ip_ijk_cell_average{};
        Array2D<Real, 0, eorder_tparm_surf - 1,
                0, N_InterP_surf - 1>
            imp_ipweights_cell_average{};
        if (nsSurfaceCellAverageRecoveryEnabled()) {
          imp_ninterp_cell_average = geometric_ninterp;
          computeIPweights<
              eorder_tparm_surf, iorder_tparm_surf, SURFIMP>(
              imp_ipweights_cell_average, imp_ip_ijk_cell_average,
              imp_xyz, imp_ijk, imp_ninterp_cell_average,
              imp_fit_order_cell_average, weight_mask, prob_lo, dx_a[lev],
              ibMarkers,
              support_visibility_mode == 2
                  ? support_visibility_condition_max
                  : std::numeric_limits<Real>::infinity(),
              iorder_tparm_surf, true);
          if (first_place_status == 0 &&
              imp_ninterp_cell_average(0) == 0) {
            ++surf_cell_average_fit_fail;
          }
          for (int image = 0; image < eorder_tparm_surf; ++image) {
            if (imp_ninterp_cell_average(image) >=
                    INTERP_THRESHOLD_SURF &&
                imp_fit_order_cell_average(image) < 2) {
              ++surf_cell_average_fit_reduced;
            }
          }
          surfimp_soa.imp_ninterp_cell_average[f_idx] =
              imp_ninterp_cell_average;
          surfimp_soa.imp_fit_order_cell_average[f_idx] =
              imp_fit_order_cell_average;
          surfimp_soa.imp_ip_ijk_cell_average[f_idx] =
              imp_ip_ijk_cell_average;
          surfimp_soa.imp_ipweights_cell_average[f_idx] =
              imp_ipweights_cell_average;
        }

        // Store imp_ninterp / ip_quality AFTER computeIPweights (WLS demotion).
        surfimp_soa.imp_ninterp[f_idx] = imp_ninterp;
        surfphys_soa.ip_quality[f_idx] =
            nsSurfaceCellAverageRecoveryEnabled()
                ? imp_ninterp_cell_average(0)
                : imp_ninterp(0);
        surfphys_soa.ip_fit_order[f_idx] =
            nsSurfaceCellAverageRecoveryEnabled()
                ? imp_fit_order_cell_average(0)
                : imp_fit_order(0);
        surfimp_soa.imp_ip_ijk[f_idx]    = imp_ip_ijk;
        surfimp_soa.imp_ipweights[f_idx] = imp_ipweights;
      } // end loop over faces in this FAB
    } // end loop over FABs
    if (timing_detail) {
      timing_surface_image_seconds =
          amrex::second() - timing_surface_image_start;
    }

    ParallelDescriptor::ReduceIntSum(surf_stencil_fail);
    ParallelDescriptor::ReduceIntSum(surf_place_fail);
    ParallelDescriptor::ReduceIntSum(surf_fit_fail);
    ParallelDescriptor::ReduceIntSum(surf_cell_average_fit_fail);
    ParallelDescriptor::ReduceIntSum(surf_cell_average_fit_reduced);
    ParallelDescriptor::ReduceIntSum(surf_tangent_shift);
    if constexpr (support_visibility_mode > 0) {
      ParallelDescriptor::ReduceIntSum(surf_visibility_supports);
      ParallelDescriptor::ReduceIntSum(surf_visibility_occluded);
      ParallelDescriptor::ReduceIntSum(surf_visibility_demoted);
      ParallelDescriptor::ReduceIntMin(surf_visibility_first_hit);
    }
    if (surf_stencil_fail > 0) {
      amrex::Print() << "computeSurfIndices: " << surf_stencil_fail
                     << " surface first-IP stencil(s) on level " << lev
                     << " extend outside their owning FAB's grown box.\n";
      amrex::Abort("computeSurfIndices: surface interpolation stencil out of box");
    }
    if (surf_place_fail > 0) {
      amrex::Print() << "computeSurfIndices: " << surf_place_fail
                     << " surface element(s) on level " << lev
                     << " failed first image-point placement.\n";
      amrex::Abort("computeSurfIndices: invalid surface first image point");
    }
    if (surf_fit_fail > 0) {
      amrex::Print() << "computeSurfIndices: " << surf_fit_fail
                     << " surface element(s) on level " << lev
                     << " have no resolvable interpolation fit at the first image point.\n";
      amrex::Abort("computeSurfIndices: invalid surface first image-point fit");
    }
    if (surf_cell_average_fit_fail > 0) {
      amrex::Print()
          << "computeSurfIndices: " << surf_cell_average_fit_fail
          << " finite-volume surface observer(s) on level " << lev
          << " have no resolvable first image-point fit.\n";
      amrex::Abort(
          "computeSurfIndices: invalid cell-average surface recovery fit");
    }
    if (nsSurfaceCellAverageRecoveryEnabled()) {
      amrex::Print()
          << "[IBM-NS-Surface] level=" << lev
          << " cell_average_fit_reductions="
          << surf_cell_average_fit_reduced << '\n';
    }
    if constexpr (sharp_feature_tangential_reconstruction) {
      amrex::Print() << "computeSurfIndices: " << surf_tangent_shift
                     << " surface image-point WLS block(s) on level " << lev
                     << " use same-edge-biased tangential reconstruction.\n";
    }
    if constexpr (support_visibility_mode > 0) {
      const Real fraction = surf_visibility_supports > 0
          ? Real(100.0) * Real(surf_visibility_occluded) /
                Real(surf_visibility_supports)
          : Real(0.0);
      amrex::Print()
          << "[IBM-Visibility] surface level " << lev
          << " mode=" << support_visibility_mode
          << " checked=" << surf_visibility_supports
          << " occluded=" << surf_visibility_occluded
          << " (" << fraction << "%)";
      if (surf_visibility_first_hit != std::numeric_limits<int>::max()) {
        amrex::Print()
            << " representative-element=" << surf_visibility_first_hit;
      }
      amrex::Print() << "\n";
    }

    // Build CSR at the coarsest level (surface is built from finest to coarsest)
    if (lev == 0) buildCSR();

    if constexpr (support_visibility_mode > 0) {
    }
    auditSurfaceAMRSupportCoverage(lev);
    if (timing_detail) {
      reportPerformanceDuration("surface_face_locate", lev,
                                timing_locate_seconds,
                                Long(ntotalfaces), true);
      reportPerformanceDuration("surface_imagepoint_weights", lev,
                                timing_surface_image_seconds,
                                Long(faces_found), false);
    }
    finishPerformanceTiming("surface_indices_total", lev,
                            timing_total_start,
                            Long(faces_found), false);
  }

  /**
   * \brief Computes surface properties (pressure, temperature, gradients) for each surface face.
   *
   * This function iterates over all surface element owned by the current process. For each face:
   * 1. Interpolates primitive variables at image points using the pre-computed weights.
   * 2. Applies the wall model to determine surface state (e.g., no-slip, adiabatic/isothermal).
   * 3. Computes gradients (e.g., dT/dn) at the surface.
   * 4. Stores the results back into the Surface Data SoA structure.
   *
   * \param stateprops MultiFab containing the fluid state properties.
   * \param cls        Pointer to the physics/closure class.
   * \param lev        Current AMR level.
  */
  void computeSURFs(
      MultiFab& prims_mf, MultiFab* conservative_cell_averages,
      const cls_t* cls, int lev
      )
  {
    BL_PROFILE("IBM::computeSURFs");


    bool use_cell_average_surface_recovery =
        nsSurfaceCellAverageRecoveryEnabled();

    // Skip if no CSR data for this level
    if (lev >= static_cast<int>(faces_per_level.size())) return;
    auto& csr = faces_per_level[lev];
    const int nfaces_local = static_cast<int>(csr.face_indices.size());

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        prims_mf.nGrow() >= surfaceInterpolationNghost(lev),
        "computeSURFs requires the IBM surface interpolation halo");
    // FillBoundary is collective across the MultiFab communicator. A rank with
    // no locally owned STL element must still participate while another rank
    // reconstructs a surface element from its ghost data.
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
    if (use_cell_average_surface_recovery) {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          conservative_cell_averages != nullptr,
          "cell-average surface recovery requires conservative support data");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          conservative_cell_averages->nGrow() >=
              surfaceInterpolationNghost(lev),
          "cell-average surface recovery requires the IBM surface halo");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          conservative_cell_averages->nComp() >= cls_t::NCONS,
          "cell-average surface recovery has insufficient state components");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          conservative_cell_averages->boxArray() == prims_mf.boxArray(),
          "surface primitive and conservative support layouts must match");
      conservative_cell_averages->FillBoundary(
          amr_p->Geom(lev).periodicity());
    }
    if (nfaces_local == 0) {
      // The invalid-state reduction below is communicator-wide.  A rank with
      // no locally owned surface faces must contribute zero instead of
      // returning early; otherwise its next collective can be matched against
      // another rank's ReduceIntSum, which is an MPI collective-order
      // violation (and commonly reports MPI_ERR_TRUNCATE).
      int invalid_surface_state_count = 0;
      ParallelDescriptor::ReduceIntSum(invalid_surface_state_count);
      if (invalid_surface_state_count > 0) {
        amrex::Abort(
            "cell-average NS surface recovery produced " +
            std::to_string(invalid_surface_state_count) +
            " non-admissible image-point conservative state(s); surface "
            "output is fail-closed and the evolved full-cell state was not "
            "modified");
      }
      return;
    }

    auto& mfab = *bmf_a[lev];
    const int nfabs_local = mfab.local_size();

    // Build device array of Array4<Real> pointers (one per local FAB)
    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
      }
      Gpu::copyAsync(Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      // Sync the async H2D copy before the host-side reads of managed members
      // (geom_offsets) below: on WSL2, host access to managed memory while a
      // copy is in flight on the stream faults. (computeAllGPs reads these
      // BEFORE its copyAsync, so it is unaffected.)
      Gpu::streamSynchronize();
    }

    auto* prims_arr = d_prims.data();
    Gpu::DeviceVector<Array4<Real>> d_conservative;
    if (use_cell_average_surface_recovery) {
      d_conservative.resize(nfabs_local);
      Vector<Array4<Real>> h_conservative(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_conservative[ifab] =
            conservative_cell_averages->array(mfi);
      }
      Gpu::copy(
          Gpu::hostToDevice, h_conservative.begin(), h_conservative.end(),
          d_conservative.begin());
    }
    auto* conservative_arr = d_conservative.data();

    // Device-accessible SoA pointers
    const int* d_face_indices  = csr.face_indices.data();

    auto* sp_pressure    = surfphys_soa.pressure.data();
    auto* sp_temperature = surfphys_soa.temperature.data();
    auto* sp_dTdn        = surfphys_soa.dTdn.data();
    auto* sp_tau_normal  = surfphys_soa.tau_normal.data();
    auto* sp_tau1        = surfphys_soa.tau1.data();
    auto* sp_tau2        = surfphys_soa.tau2.data();
    auto* sp_pshock      = surfphys_soa.pressure_shock_sensor.data();
    auto* sp_pfallback   = surfphys_soa.pressure_fallback.data();
    const int* sp_ifab   = surfphys_soa.ifab.data();

    const auto* si_imp_ip_ijk    = surfimp_soa.imp_ip_ijk.data();
    const auto* si_imp_ipweights = surfimp_soa.imp_ipweights.data();
    const auto* si_imp_ninterp_cell_average =
        surfimp_soa.imp_ninterp_cell_average.data();
    const auto* si_imp_ip_ijk_cell_average =
        surfimp_soa.imp_ip_ijk_cell_average.data();
    const auto* si_imp_ipweights_cell_average =
        surfimp_soa.imp_ipweights_cell_average.data();
    const auto* si_imp_xyz       = surfimp_soa.imp_xyz.data();
    const auto* si_imp_ijk       = surfimp_soa.imp_ijk.data();
    const auto* si_disIM         = surfimp_soa.disIM.data();
    const auto* si_imp_ninterp   = surfimp_soa.imp_ninterp.data();

    auto const* lf_ptr = LocalFrame_a.data();
    auto const* se_ptr = SurfElem_a.data();

    // Capture transforms and geom_offsets for body→world rotation
    AMREX_ALWAYS_ASSERT(ngeom <= MAX_NGEOM);
    const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    GpuArray<int, MAX_NGEOM + 1> geom_off;
    for (int ii = 0; ii < ngeom_local; ii++) {
        transforms[ii] = transform_a[ii];
        geom_off[ii]   = geom_offsets[ii];
    }
    geom_off[ngeom_local] = geom_offsets[ngeom_local];

    auto* copy = this;
    const auto prob_lo_lev = amr_p->Geom(lev).ProbLoArray();
    const auto dx_lev = dx_a[lev];
    const bool cartesian_geometry = !amr_p->Geom(lev).IsRZ();
    Gpu::DeviceScalar<int> invalid_surface_state(0);
    int* invalid_surface_state_ptr = invalid_surface_state.dataPtr();
    ParallelFor(nfaces_local, [=] AMREX_GPU_DEVICE (int ii) noexcept
    {
      // See the GP path above: force the optional compatibility inputs into
      // the extended-lambda capture set before entering if constexpr.
      amrex::ignore_unused(
          si_imp_xyz, prob_lo_lev, dx_lev, cls, cartesian_geometry);
      // Global face index from CSR
      int f_idx = d_face_indices[ii];
      int ifab = sp_ifab[f_idx];
      auto prims = prims_arr[ifab];
      auto conservative = use_cell_average_surface_recovery
          ? conservative_arr[ifab]
          : prims_arr[ifab];

      // Find geometry index for this face (to get the right transform)
      int gIdx = 0;
      for (int g = 0; g < ngeom_local; ++g) {
          if (f_idx >= geom_off[g] && f_idx < geom_off[g + 1]) { gIdx = g; break; }
      }
      const auto& T = transforms[gIdx];

      // 1) Reconstruct local orthonormal frame — rotate body→world
      const auto& frame = lf_ptr[f_idx];

      amrex::Real n_w[AMREX_SPACEDIM], t1_w[AMREX_SPACEDIM];
#if (AMREX_SPACEDIM == 3)
      amrex::Real t2_w[AMREX_SPACEDIM];
#endif
      T.rotate_to_world(frame.normal,   n_w);
      T.rotate_to_world(frame.tangent1, t1_w);
#if (AMREX_SPACEDIM == 3)
      T.rotate_to_world(frame.tangent2, t2_w);
#endif

      Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          nvec(d)  = n_w[d];
          t1vec(d) = t1_w[d];
#if (AMREX_SPACEDIM == 3)
          t2vec(d) = t2_w[d];
#else
          t2vec(d) = Real(0.0);
#endif
      }

      // Surface centroid coordinates — transform body→world
      Array1D<Real, 0, AMREX_SPACEDIM - 1> xyz;
      {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(se_ptr[f_idx].centroid[0], se_ptr[f_idx].centroid[1]);
#else
        Point c_body(se_ptr[f_idx].centroid[0], se_ptr[f_idx].centroid[1], se_ptr[f_idx].centroid[2]);
#endif
        Point c_world = T.to_world(c_body);
        xyz(0) = c_world.x(); xyz(1) = c_world.y();
#if (AMREX_SPACEDIM == 3)
        xyz(2) = c_world.z();
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = se_ptr[f_idx].centroid[dd];
        Point c_world = T.to_world(c_body);
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) xyz(dd) = c_world[dd];
#endif
      }

      Array1D<Real, 0, AMREX_SPACEDIM - 1> raw_xyz = xyz;
      apply_volume_surface_point_override(raw_xyz, xyz);
      LocalFrame surface_frame;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        surface_frame.normal[d] = n_w[d];
        surface_frame.tangent1[d] = t1_w[d];
#if (AMREX_SPACEDIM == 3)
        surface_frame.tangent2[d] = t2_w[d];
#endif
      }
      apply_volume_surface_frame_override(xyz, surface_frame);
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        nvec(d) = surface_frame.normal[d];
        t1vec(d) = surface_frame.tangent1[d];
#if (AMREX_SPACEDIM == 3)
        t2vec(d) = surface_frame.tangent2[d];
#else
        t2vec(d) = Real(0.0);
#endif
      }

      // 2) Zero-initialize primsNormal
      Array2D<Real, 0, eorder_tparm_surf + 1, 0, cls_t::NPRIM - 1> primsNormal;
      for (int p = 0; p <= eorder_tparm_surf + 1; ++p) {
          for (int n = 0; n < cls_t::NPRIM; ++n) {
              primsNormal(p,n) = Real(0.0);
          }
      }

      // Number of leading valid image points.  The surface wall model and
      // one-sided gradient use the quadratic form only when both are valid.
      int n_valid_surf = 0;
      for (int kk = 0; kk < eorder_tparm_surf; ++kk) {
          const int support_count = use_cell_average_surface_recovery
              ? si_imp_ninterp_cell_average[f_idx](kk)
              : si_imp_ninterp[f_idx](kk);
          if (support_count < INTERP_THRESHOLD_SURF) break;
          n_valid_surf = kk + 1;
      }

      {
        if (use_cell_average_surface_recovery) {
          const bool admissible =
              copy->template
                  interpolateConservativeCellAveragesToPointPrimitives<
                      eorder_tparm_surf, iorder_tparm_surf>(
                      si_imp_ip_ijk_cell_average[f_idx],
                      si_imp_ipweights_cell_average[f_idx],
                      si_imp_ninterp_cell_average[f_idx],
                      conservative, cls, primsNormal);
          if (!admissible) {
            Gpu::Atomic::Add(invalid_surface_state_ptr, 1);
            return;
          }
        } else {
          // Surface B/C and the legacy production observer use cached point
          // interpolation. Mode B changes only the supplied support values.
          copy->template interpolateIMs<
              eorder_tparm_surf, iorder_tparm_surf>(
                  si_imp_ip_ijk[f_idx], si_imp_ipweights[f_idx],
                  prims, primsNormal);
        }
        for (int iip = 2; iip < 2 + eorder_tparm_surf; ++iip) {
          copy->template global2local<eorder_tparm_surf>(
              iip, primsNormal, nvec, t1vec, t2vec);
        }
      }

      // Apply the numerical wall model at the BI for the production path and
      // Surface B/C.  Surface A already contains the exact BI state.
      int type_solid_bc = 0;
      ibm_pressure_compatibility_t pressure_compat{};
      if constexpr (use_ns_noslip_pressure_compatibility) {
        if (
            n_valid_surf >= 2) {
          pressure_compat =
              copy->template ibm_viscous_pressure_compatibility<
                  iorder_tparm_surf, N_InterP_surf>(
                      si_imp_ip_ijk[f_idx], si_imp_xyz[f_idx],
                      si_imp_ninterp[f_idx](0), prims, nvec,
                      prob_lo_lev, dx_lev, cls, cartesian_geometry);
        }
      }
      if constexpr (uses_prescribed_gradient_pressure_closure) {
        pressure_compat.dpdn =
            param::pressure_normal_derivative(xyz, nvec);
        pressure_compat.valid =
            amrex::Math::isfinite(pressure_compat.dpdn) ? 1 : 0;
      }
      if constexpr (sharp_feature_pressure_limiter) {
        const auto& edge = se_ptr[f_idx];
        const Real h = amrex::min(dx_lev[0], dx_lev[1]);
        pressure_compat.feature_pressure_fallback =
            sharpFeaturePressureFallback(
                edge, Real(0.5) * edge.measure, h);
      }
      {
        ibm_detail::dispatch_compute_surfIB<eorder_tparm_surf, wallmodel>(
            xyz, nvec, t1vec, t2vec, si_disIM[f_idx], n_valid_surf,
            primsNormal, type_solid_bc, cls, pressure_compat);
      }

      // 6) Extract surface quantities
      Real P_surf = primsNormal(1, cls_t::QPRES);
      Real T_surf = primsNormal(1, cls_t::QT);
      Real pressure_shock_sensor = Real(0.0);
      Real pressure_fallback = Real(0.0);
      {
        ibm_detail::dispatch_pressure_closure_diagnostics<
            eorder_tparm_surf, wallmodel, cls_t>(
                primsNormal, si_disIM[f_idx], n_valid_surf,
                pressure_compat,
                pressure_shock_sensor, pressure_fallback);
      }


      // Wall-normal gradients via 2nd-order one-sided stencil (surface + IP1 +
      // IP2, general spacing) when >=2 image points are valid; else 1st-order
      // one-sided (= the previous (IP1 - surface)/dis formula).
      //   dT/dn = grad(T)·n   (heat flux)
      //   tau   = mu * d(u_tangential)/dn   (local frame: QU=normal,
      //           QV=tangent1, QW=tangent2)
      Real dTdn = ibm_wall_normal_deriv<eorder_tparm_surf>(
          primsNormal, cls_t::QT, si_disIM[f_idx], n_valid_surf);

      Real mu = cls->visc(T_surf);
      Real tau_normal = Real(0.0);
      if constexpr (supports_complete_viscous_surface_traction) {
        const Real xi = cls->xi(T_surf);
        const Real dun_dn = ibm_wall_normal_deriv<eorder_tparm_surf>(
            primsNormal, cls_t::QU, si_disIM[f_idx], n_valid_surf);
        tau_normal = (Real(4.0 / 3.0) * mu + xi) * dun_dn;
      }
      Real tau1 = mu * ibm_wall_normal_deriv<eorder_tparm_surf>(
          primsNormal, cls_t::QV, si_disIM[f_idx], n_valid_surf);
#if (AMREX_SPACEDIM == 3)
      Real tau2 = mu * ibm_wall_normal_deriv<eorder_tparm_surf>(
          primsNormal, cls_t::QW, si_disIM[f_idx], n_valid_surf);
#else
      Real tau2 = Real(0.0);
#endif

      // 7) Store results
      sp_pressure[f_idx]    = P_surf;
      sp_temperature[f_idx] = T_surf;
      sp_dTdn[f_idx]        = dTdn;
      sp_tau_normal[f_idx]   = tau_normal;
      sp_tau1[f_idx]        = tau1;
      sp_tau2[f_idx]        = tau2;
      sp_pshock[f_idx]      = pressure_shock_sensor;
      sp_pfallback[f_idx]   = pressure_fallback;
    }); // end ParallelFor

    // Ensure GPU writes are visible to CPU before gatherSurfData / plotSURF
    Gpu::streamSynchronize();
    int invalid_surface_state_count =
        invalid_surface_state.dataValue();
    ParallelDescriptor::ReduceIntSum(invalid_surface_state_count);
    if (invalid_surface_state_count > 0) {
      amrex::Abort(
          "cell-average NS surface recovery produced " +
          std::to_string(invalid_surface_state_count) +
          " non-admissible image-point conservative state(s); surface output "
          "is fail-closed and the evolved full-cell state was not modified");
    }

  }

  // Conservative crossing loads and physical/load-consistent surface forces.
  #include "ibm_solver_surface_load.h"


//============================================================================
///--------------------------- private functions -----------------------------
private:

  // Interpolation, extrapolation, coordinate-transform helpers
  #include "ibm_solver_interp.h"

  // Geometry I/O, VTK output, MPI gather, CSR builder
  #include "ibm_solver_io.h"

}; // end class ibm_solver_t

#endif // IBM_SOLVER_H_
