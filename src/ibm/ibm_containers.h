#ifndef IBM_CONTAINERS_H_
#define IBM_CONTAINERS_H_

// ============================================================================
// ibm_containers.h — Container types, constants, and helper types for the IBM solver
//
// Contains:
//   0. IBFab / IBMultiFab  : IBM-specific AMR containers (marker storage only)
//   1. ibm_detail namespace: SFINAE helpers for wall-model dispatch
//   2. ipow()              : Constexpr integer power
//   3. Constants           : Dimension indices, thresholds, image-point factors
//   4. surfImp_t           : SoA storage for surface image-point data
//   6. surfPhys_t          : SoA storage for reconstructed surface fields
//   7. is_gp_store_view    : Type trait selecting GP vs surface interp paths
//   8. GPStoreView/GPStore : Level-wide flattened ghost-point storage (CSR)
//   9. FaceCSR             : CSR structure for per-FAB face iteration
//  10. CheckMode           : Interpolation stencil check policy
//
// Per-FAB ghost-point data lived in a now-removed ``gpData_t`` member of
// ``IBFab``.  All ghost-point geometry/interpolation data is now flattened
// into the level-wide CSR ``GPStore``; ``IBFab`` only owns the marker FAB.
// ============================================================================

#include <algorithm>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

#include <AMReX_FabArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>

#ifdef AMREX_USE_CGAL
#include "ibm_backend_cgal.h"
#else
#include "ibm_backend_bvh.h"
#endif

// ============================================================================
// 0. IBFab / IBMultiFab — IBM-specific AMR containers
// ============================================================================

/// \brief IBFab holds the marker data for a single AMR box.
/// \tparam marker_t Type of the marker data (typically uint8_t)
///
/// Ghost-point geometry/interpolation data is NOT stored here — it lives in
/// the level-wide CSR ``GPStore``.  IBFab is a thin wrapper over
/// ``amrex::BaseFab<marker_t>`` that disables implicit deep-copy.
template<typename marker_t>
class IBFab : public amrex::BaseFab<marker_t> {
public:
  // FabArray communication creates temporary FABs before defining their
  // storage.  Match BaseFab's default-constructible interface so standard
  // FillBoundary/ParallelCopy paths work for marker data.
  IBFab() = default;

  explicit IBFab(const amrex::Box& b, int ncomp,
                 bool alloc = true, bool shared = false, amrex::Arena* ar = nullptr)
      : amrex::BaseFab<marker_t>(b, ncomp, alloc, shared, ar) {}

  // MakeType ctor (alias/deep-copy): forwarded straight to BaseFab.
  explicit IBFab(const IBFab<marker_t>& rhs, amrex::MakeType make_type, int scomp, int ncomp)
      : amrex::BaseFab<marker_t>(rhs, make_type, scomp, ncomp) {}

  ~IBFab() = default;

  IBFab(const IBFab&) = delete;
  IBFab& operator=(const IBFab&) = delete;

  IBFab(IBFab&&) noexcept = default;
  IBFab& operator=(IBFab&&) noexcept = default;
};

/// \brief IBMultiFab holds an array of IBFab on a level
/// \tparam marker_t Type of the marker data
template<typename marker_t>
class IBMultiFab : public amrex::FabArray<IBFab<marker_t>> {
public:
  using Fab  = IBFab<marker_t>;
  using Base = amrex::FabArray<Fab>;

  /// \brief Default MFInfo that routes allocation to managed (CPU+GPU) memory.
  static amrex::MFInfo DataMFInfo() {
      amrex::MFInfo info;
      info.SetArena(amrex::The_Managed_Arena());
      return info;
  }

  explicit IBMultiFab(
      const amrex::BoxArray& bxs, const amrex::DistributionMapping& dm, int nvar, int ngrow,
      const amrex::MFInfo& info = DataMFInfo(),
      const amrex::FabFactory<Fab>& factory = amrex::DefaultFabFactory<Fab>())
      : Base(bxs, dm, nvar, ngrow, info, factory) {}

  ~IBMultiFab() = default;

  IBMultiFab(IBMultiFab&&) noexcept = default;
  IBMultiFab& operator=(IBMultiFab&&) noexcept = default;

  IBMultiFab(const IBMultiFab&) = delete;
  IBMultiFab& operator=(const IBMultiFab&) = delete;

  /// \brief Copy marker components from this IBMultiFab into a Real MultiFab (e.g. for plotfile).
  void copytoRealMF(amrex::MultiFab& mf, int ibcomp, int mfcomp) {
    const int ncomp_copy = std::min(this->nComp() - ibcomp, mf.nComp() - mfcomp);
    if (ncomp_copy <= 0) return;

    for (amrex::MFIter mfi(*this, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& src = this->get(mfi).array();
      auto const& dst = mf.array(mfi);

      amrex::ParallelFor(bx, ncomp_copy,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          dst(i, j, k, mfcomp + n) = static_cast<amrex::Real>(src(i, j, k, ibcomp + n));
        });
    }
  }
};

// ============================================================================
// 1. SFINAE helpers for wall-model dispatch
// ============================================================================

struct ibm_pressure_compatibility_t;

namespace ibm_detail {
template <class...>
using void_t = void;

// Robust detection idiom: the void_t<Op<Args...>> probe lives ONLY in the
// partial-specialisation's template-argument list (a deduced/SFINAE context),
// so a fully-failing probe (Op<Args...> ill-formed — e.g. a wall model that
// lacks the signature being probed) resolves to false_type instead of a hard
// error.  The _v alias passes a plain `void` for the AlwaysVoid slot.
template <class AlwaysVoid, template <class...> class Op, class... Args>
struct is_detected_impl : std::false_type {};

template <template <class...> class Op, class... Args>
struct is_detected_impl<void_t<Op<Args...>>, Op, Args...> : std::true_type {};

template <template <class...> class Op, class... Args>
constexpr bool is_detected_v = is_detected_impl<void, Op, Args...>::value;

template <class WM, class... Args>
using compute_surfIB_expr = decltype(WM::compute_surfIB(std::declval<Args>()...));

/// nvcc workaround: dispatch compute_surfIB outside constexpr-if in __device__ lambda.
/// Priority: pressure-compatibility + disIM > disIM-aware (2nd-order Neumann)
/// > full (xyz,n,t1,t2,q,...) > reduced (xyz,n,q,...).  Custom wall models
/// that define only an older signature keep working unchanged.
///
/// The disIM-aware overload is order-templated; its extrapolation order EO
/// appears only as EO+1 / EO-1 in the parameter array bounds (non-deduced
/// contexts), so it is made deducible by a leading std::integral_constant<int,EO>
/// tag.  This keeps detection on the standard argument-deduction path (no
/// `::template` on a possibly-non-template member, which is not portably
/// SFINAE-friendly).
template <int eorder, class wallmodel_t, class cls_type,
          class PointArr, class Vec1D, class DisArr, class Prims2D>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dispatch_compute_surfIB(
    const PointArr& xyz,
    const Vec1D& nvec,
    const Vec1D& t1vec,
    const Vec1D& t2vec,
    const DisArr& disIM,
    int n_valid,
    Prims2D& primsNormal,
    int type_solid_bc,
    const cls_type* cls,
    const ibm_pressure_compatibility_t& pressure_compat)
{
    using eo_tag = std::integral_constant<int, eorder>;
    constexpr bool has_compat = is_detected_v<
        compute_surfIB_expr,
        wallmodel_t,
        eo_tag,
        const Vec1D&, const Vec1D&, const Vec1D&, const Vec1D&,
        const DisArr&, int, Prims2D&, const int, const cls_type*,
        const ibm_pressure_compatibility_t&>;
    constexpr bool has_disim = is_detected_v<
        compute_surfIB_expr,
        wallmodel_t,
        eo_tag,
        const Vec1D&, const Vec1D&, const Vec1D&, const Vec1D&,
        const DisArr&, int, Prims2D&, const int, const cls_type*>;
    constexpr bool has_full = is_detected_v<
        compute_surfIB_expr,
        wallmodel_t,
        const Vec1D&,
        const Vec1D&,
        const Vec1D&,
        const Vec1D&,
        Prims2D&,
        const int,
        const cls_type*>;
    if constexpr (has_compat) {
        wallmodel_t::compute_surfIB(
            eo_tag{}, xyz, nvec, t1vec, t2vec, disIM, n_valid, primsNormal,
            type_solid_bc, cls, pressure_compat);
    } else if constexpr (has_disim) {
        amrex::ignore_unused(pressure_compat);
        wallmodel_t::compute_surfIB(
            eo_tag{}, xyz, nvec, t1vec, t2vec, disIM, n_valid, primsNormal, type_solid_bc, cls);
    } else if constexpr (has_full) {
        amrex::ignore_unused(disIM, n_valid, pressure_compat);
        wallmodel_t::compute_surfIB(xyz, nvec, t1vec, t2vec, primsNormal, type_solid_bc, cls);
    } else {
        amrex::ignore_unused(t1vec, t2vec, disIM, n_valid, pressure_compat);
        wallmodel_t::compute_surfIB(xyz, nvec, primsNormal, type_solid_bc, cls);
    }
}

} // namespace ibm_detail

// ============================================================================
// 2. Constexpr integer power
// ============================================================================

AMREX_GPU_HOST_DEVICE constexpr int ipow(int base, int exp) {
    return (exp == 0) ? 1 : base * ipow(base, exp - 1);
}

// ============================================================================
// 2b. Second-order one-sided wall-normal helpers (image-point IBM)
//
// Geometry along the outward wall normal: the surface (boundary intercept) sits
// at normal-distance x = 0, the image points IP1, IP2 at x1 = disIM(0),
// x2 = disIM(1).  Slot layout in the primsNormal array:
//   q(1,.) = surface (wall)   q(2,.) = IP1   q(3,.) = IP2
//
//   ibm_normal_grad_surface : surface value for a prescribed dphi/dn|_0
//   ibm_zero_grad_surface   : prescribed-gradient helper with dphi/dn|_0 = 0
//   ibm_fluid_extrap_surface: unconstrained one-sided extrapolation to the wall
//   ibm_wall_normal_deriv   : one-sided dphi/dn|_0 (heat flux, shear)
//
// The derivative-based helpers use the general (non-uniform) spacing x1, x2.
// Every helper reduces to the previous single-image-point form when EO < 2,
// when fewer than two image points are valid (n_valid < 2), or when the image
// points are nearly coincident.  q(3,.)/disIM(1) are read only inside
// `if constexpr (EO >= 2)`, so they are never instantiated for EO = 1.
//
// Coefficients are the derivatives of the Lagrange basis through {0, x1, x2}:
//   cs = -(x1+x2)/(x1 x2),  c1 = x2/(x1 (x2-x1)),  c2 = -x1/(x2 (x2-x1))
// Zero-gradient surface value: phi_s = -(c1 phi1 + c2 phi2)/cs,
// with -1/cs = x1 x2/(x1+x2).  For x2 = 2 x1 this collapses to (4 phi1 - phi2)/3.
// ============================================================================

enum class ibm_pressure_closure_t : int {
    zero_gradient = 0,
    fluid_extrapolation = 1,
    prescribed_gradient = 2,
    shock_aware_fluid_extrapolation = 3,
    euler_slip_analytic_curvature = 4,
    navier_stokes_noslip_viscous_compatibility = 5
};

/// Volume ghost-state reconstruction used by the shared-GP path.  The legacy
/// method first interpolates one or more image points and then extrapolates
/// through the boundary intercept.  The BI-constrained method fits the visible
/// fluid support subject to the wall condition and evaluates that same
/// polynomial directly at the real ghost-cell centre.
enum class ibm_shared_gp_reconstruction_t : int {
    image_point = 0,
    boundary_intercept_constrained = 1
};

/// Optional normal pressure derivative supplied by a solver-side momentum
/// compatibility reconstruction.  `valid=0` is a deliberate request for the
/// wall model to use its documented fallback; dpdn is then ignored.
struct ibm_pressure_compatibility_t {
    Real dpdn = Real(0.0);
    int valid = 0;
    // Minimum return-to-IP1 fraction requested by the geometry path.  The
    // 2-D sharp-feature treatment uses this to regularise the non-smooth
    // corner neighbourhood while leaving the pressure closure unchanged
    // away from the O(h) feature band.
    Real feature_pressure_fallback = Real(0.0);
};

/// Shape operator in the local surface basis used by primsNormal:
///   B_ab = t_a . (grad n) t_b,  a,b in {1,2}.
/// The IBM normal points from solid into fluid, so a convex sphere/circle has
/// positive principal curvature (+1/R).  In 2-D only k11 is used.  Analytic
/// geometry callbacks return valid=0 outside their declared smooth region;
/// the wall closure then falls back to fluid-side pressure extrapolation.
struct ibm_shape_operator_t {
    Real k11;
    Real k12;
    Real k22;
    int valid;
};

struct ibm_pressure_reconstruction_t {
    Real value;
    Real unlimited_value;
    Real shock_sensor;
    Real fallback_fraction;
};

struct ibm_surface_force_audit_t {
    GpuArray<Real, 3> force{{Real(0.0), Real(0.0), Real(0.0)}};
    GpuArray<Real, 3> pressure_force{{Real(0.0), Real(0.0), Real(0.0)}};
    GpuArray<Real, 3> viscous_force{{Real(0.0), Real(0.0), Real(0.0)}};
    GpuArray<Real, 3> normal_closure{{Real(0.0), Real(0.0), Real(0.0)}};
    Real measure = Real(0.0);
    Long faces = 0;
    Long nonfinite_faces = 0;
};

/// Reconstruct the surface value from image-point values and a prescribed
/// outward-normal derivative.  The derivative and disIM must use the same
/// normal orientation (surface -> fluid).  For EO>=2 this inverts the same
/// one-sided Lagrange derivative used by ibm_wall_normal_deriv; EO=1 uses the
/// linear relation phi_s = phi_1 - x_1 dphi/dn.
template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_normal_grad_surface(
    const QArr& q, int n, const DisArr& disIM, int n_valid,
    Real normal_derivative)
{
    const Real u1 = q(2, n);   // IP1
    constexpr Real eps = Real(1.0e-12);
    if constexpr (EO >= 3) {
        // Higher-order prescribed-gradient form on {0, x1, x2, x3}, using
        // Lagrange-basis derivatives at s=0:
        //   cs = -(1/x1 + 1/x2 + 1/x3)
        //   c1 = x2 x3/(x1 (x1-x2)(x1-x3)),  c2/c3 by cyclic index swap.
        // phi_s = (dphi/dn - c1 u1 - c2 u2 - c3 u3)/cs.
        if (n_valid >= 3) {
            const Real x1 = disIM(0);
            const Real x2 = disIM(1);
            const Real x3 = disIM(2);
            const bool distinct =
                x1 > Real(0.0) &&
                amrex::Math::abs(x2 - x1) > eps * amrex::max(x1, x2) &&
                amrex::Math::abs(x3 - x2) > eps * amrex::max(x2, x3) &&
                amrex::Math::abs(x3 - x1) > eps * amrex::max(x1, x3);
            if (distinct) {
                const Real u2 = q(3, n);
                const Real u3 = q(4, n);
                const Real cs = -(Real(1.0)/x1 + Real(1.0)/x2 + Real(1.0)/x3);
                const Real c1 = x2 * x3 / (x1 * (x1 - x2) * (x1 - x3));
                const Real c2 = x1 * x3 / (x2 * (x2 - x1) * (x2 - x3));
                const Real c3 = x1 * x2 / (x3 * (x3 - x1) * (x3 - x2));
                return (normal_derivative -
                        (c1 * u1 + c2 * u2 + c3 * u3)) / cs;
            }
        }
    }
    if constexpr (EO >= 2) {
        if (n_valid >= 2) {
            const Real x1 = disIM(0);
            const Real x2 = disIM(1);
            if (x1 > Real(0.0) && amrex::Math::abs(x2 - x1) > eps * amrex::max(x1, x2)) {
                const Real u2 = q(3, n);   // IP2
                const Real c1 =  x2 / (x1 * (x2 - x1));
                const Real c2 = -x1 / (x2 * (x2 - x1));
                const Real cs = -(x1 + x2) / (x1 * x2);
                return (normal_derivative - (c1 * u1 + c2 * u2)) / cs;
            }
        }
    } else {
        amrex::ignore_unused(disIM, n_valid);
    }
    const Real x1 = disIM(0);
    return (n_valid >= 1 && x1 > Real(0.0))
               ? u1 - x1 * normal_derivative
               : u1;
}

template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_zero_grad_surface(const QArr& q, int n, const DisArr& disIM, int n_valid)
{
    return ibm_normal_grad_surface<EO>(
        q, n, disIM, n_valid, Real(0.0));
}

/// Reconstruct a wall value without imposing a Neumann condition.  Pressure at
/// an inviscid wall is determined by the interior solution, not by a universal
/// dp/dn=0 condition.  A linear extrapolation through IP1/IP2 gives an O(h^2)
/// wall value; IP1/IP2/IP3 give an O(h^3) value.  The EO=1 fallback is the
/// nearest image-point value and therefore retains the legacy behaviour.
template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_fluid_extrap_surface(
    const QArr& q, int n, const DisArr& disIM, int n_valid)
{
    const Real u1 = q(2, n);
    constexpr Real eps = Real(1.0e-12);
    if constexpr (EO >= 3) {
        if (n_valid >= 3) {
            const Real x1 = disIM(0);
            const Real x2 = disIM(1);
            const Real x3 = disIM(2);
            const bool distinct =
                x1 > Real(0.0) &&
                amrex::Math::abs(x2 - x1) > eps * amrex::max(x1, x2) &&
                amrex::Math::abs(x3 - x2) > eps * amrex::max(x2, x3) &&
                amrex::Math::abs(x3 - x1) > eps * amrex::max(x1, x3);
            if (distinct) {
                const Real u2 = q(3, n);
                const Real u3 = q(4, n);
                const Real l1 = x2 * x3 / ((x1 - x2) * (x1 - x3));
                const Real l2 = x1 * x3 / ((x2 - x1) * (x2 - x3));
                const Real l3 = x1 * x2 / ((x3 - x1) * (x3 - x2));
                return l1 * u1 + l2 * u2 + l3 * u3;
            }
        }
    }
    if constexpr (EO >= 2) {
        if (n_valid >= 2) {
            const Real x1 = disIM(0);
            const Real x2 = disIM(1);
            if (x1 > Real(0.0) &&
                amrex::Math::abs(x2 - x1) > eps * amrex::max(x1, x2)) {
                return (x2 * u1 - x1 * q(3, n)) / (x2 - x1);
            }
        }
    } else {
        amrex::ignore_unused(disIM, n_valid);
    }
    return u1;
}

/// Irreversible-compression sensor on the first two image points.  IP1 is
/// nearest the wall and IP2 is farther into the fluid.  The reversible part
/// of a compression is estimated by a trapezoidal integral of c^2 d(rho),
/// leaving only pressure rise in excess of the local acoustic/isentropic
/// relation.  Expansions and smooth isentropic curvature gradients therefore
/// do not trigger the shock fallback.  The normalized sensor is clamped to two.
/// A value of zero also denotes an unavailable two-point stencil; callers use
/// n_valid separately to distinguish that case.
template <int EO, typename QArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_pressure_shock_sensor(
    const QArr& q, int np, int nrho, int nc, int n_valid)
{
    if constexpr (EO >= 2) {
        if (n_valid >= 2) {
            const Real p1 = q(2, np);
            const Real p2 = q(3, np);
            const Real rho1 = q(2, nrho);
            const Real rho2 = q(3, nrho);
            const Real c1 = q(2, nc);
            const Real c2 = q(3, nc);
            if (!(p1 > Real(0.0)) || !(p2 > Real(0.0)) ||
                !(rho1 > Real(0.0)) || !(rho2 > Real(0.0)) ||
                !(c1 > Real(0.0)) || !(c2 > Real(0.0)) ||
                !amrex::Math::isfinite(p1) || !amrex::Math::isfinite(p2) ||
                !amrex::Math::isfinite(rho1) || !amrex::Math::isfinite(rho2) ||
                !amrex::Math::isfinite(c1) || !amrex::Math::isfinite(c2)) {
                return Real(2.0);
            }
            const Real dp = p1 - p2;
            const Real drho = rho1 - rho2;
            if (dp <= Real(0.0) || drho <= Real(0.0)) return Real(0.0);
            const Real c2_bar = Real(0.5) * (c1 * c1 + c2 * c2);
            const Real irreversible = dp - c2_bar * drho;
            const Real sensor = Real(2.0) *
                amrex::max(irreversible, Real(0.0)) / (p1 + p2);
            return amrex::min(sensor, Real(2.0));
        }
    } else {
        amrex::ignore_unused(q, np, nrho, nc, n_valid);
    }
    return Real(0.0);
}

/// Fluid-side pressure extrapolation with a continuous shock-triggered return
/// to the nearest image-point value.  The unlimited reconstruction is retained
/// exactly for sensor <= sensor_low.  Between sensor_low and sensor_high a
/// cubic smoothstep blends to IP1; at and above sensor_high the result is IP1.
/// Invalid/non-positive extrapolates also return IP1 independently of the
/// sensor.  The two thresholds use the shock sensor above and must
/// satisfy 0 <= sensor_low < sensor_high <= 2.
template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ibm_pressure_reconstruction_t ibm_shock_aware_pressure_surface(
    const QArr& q, int np, int nrho, int nc,
    const DisArr& disIM, int n_valid,
    Real sensor_low, Real sensor_high,
    Real feature_fallback = Real(0.0))
{
    const Real p1 = q(2, np);
    const Real unlimited = ibm_fluid_extrap_surface<EO>(
        q, np, disIM, n_valid);

    if constexpr (EO < 2) {
        amrex::ignore_unused(sensor_low, sensor_high);
        return {p1, unlimited, Real(0.0), Real(1.0)};
    }

    if (n_valid < 2) {
        return {p1, unlimited, Real(0.0), Real(1.0)};
    }

    const Real sensor = ibm_pressure_shock_sensor<EO>(
        q, np, nrho, nc, n_valid);
    if (!(p1 > Real(0.0)) || !amrex::Math::isfinite(p1)) {
        Real safe_value = p1;
        if constexpr (EO >= 2) {
            const Real p2 = q(3, np);
            if (p2 > Real(0.0) && amrex::Math::isfinite(p2)) safe_value = p2;
        }
        return {safe_value, unlimited, sensor, Real(1.0)};
    }
    if (!(unlimited > Real(0.0)) || !amrex::Math::isfinite(unlimited)) {
        return {p1, unlimited, sensor, Real(1.0)};
    }

    constexpr Real eps = Real(1.0e-12);
    const Real width = amrex::max(sensor_high - sensor_low, eps);
    const Real xi = amrex::max(
        Real(0.0), amrex::min(Real(1.0), (sensor - sensor_low) / width));
    Real fallback = xi * xi * (Real(3.0) - Real(2.0) * xi);
    // A sharp corner is a genuine singular point of the inviscid wall
    // solution: a normal polynomial through IP1/IP2 is not a consistent
    // smooth extrapolation there.  The geometry path supplies a compact O(h)
    // fallback floor.  Taking the maximum preserves any stronger shock
    // fallback and makes the transition exactly recover the unlimited
    // closure outside the feature band.
    const Real bounded_feature_fallback = amrex::max(
        Real(0.0), amrex::min(Real(1.0), feature_fallback));
    fallback = amrex::max(fallback, bounded_feature_fallback);
    Real value = (Real(1.0) - fallback) * unlimited + fallback * p1;
    if (!(value > Real(0.0)) || !amrex::Math::isfinite(value)) {
        value = p1;
        fallback = Real(1.0);
    }
    return {value, unlimited, sensor, fallback};
}

namespace ibm_detail {

template <typename WallModel>
using pressure_closure_expr = decltype(WallModel::pressure_closure);

template <typename WallModel>
AMREX_GPU_HOST_DEVICE constexpr bool
uses_navier_stokes_noslip_pressure_compatibility()
{
    if constexpr (is_detected_v<pressure_closure_expr, WallModel>) {
        return WallModel::pressure_closure ==
               ibm_pressure_closure_t::navier_stokes_noslip_viscous_compatibility;
    }
    return false;
}

/// Surface-only pressure diagnostics with a no-op fallback for custom wall
/// models.  Built-in walls expose the compile-time policy and thresholds, but
/// ibm_solver.h does not need to include or know their concrete definitions.
template <int EO, typename WallModel, typename Cls, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dispatch_pressure_closure_diagnostics(
    const QArr& q, const DisArr& disIM, int n_valid,
    const ibm_pressure_compatibility_t& pressure_compat,
    Real& shock_sensor, Real& fallback_fraction)
{
    shock_sensor = Real(-1.0);
    fallback_fraction = Real(-1.0);
    if constexpr (is_detected_v<pressure_closure_expr, WallModel>) {
        shock_sensor = ibm_pressure_shock_sensor<EO>(
            q, Cls::QPRES, Cls::QRHO, Cls::QC, n_valid);
        if constexpr (
            WallModel::pressure_closure ==
            ibm_pressure_closure_t::shock_aware_fluid_extrapolation) {
            const auto result = ibm_shock_aware_pressure_surface<EO>(
                q, Cls::QPRES, Cls::QRHO, Cls::QC, disIM, n_valid,
                WallModel::pressure_sensor_low,
                WallModel::pressure_sensor_high,
                pressure_compat.feature_pressure_fallback);
            shock_sensor = result.shock_sensor;
            fallback_fraction = result.fallback_fraction;
        } else if constexpr (
            WallModel::pressure_closure ==
            ibm_pressure_closure_t::prescribed_gradient) {
            // The positive prescribed-gradient helper falls back to IP1 when
            // its candidate is non-positive or non-finite.  Mirror that exact
            // branch in the output-only diagnostic instead of leaving the
            // generic -1 sentinel, which would hide a real order reduction.
            shock_sensor = Real(-1.0);
            const Real candidate = ibm_normal_grad_surface<EO>(
                q, Cls::QPRES, disIM, n_valid, pressure_compat.dpdn);
            const bool used = pressure_compat.valid != 0 &&
                              amrex::Math::isfinite(pressure_compat.dpdn) &&
                              candidate > Real(0.0) &&
                              amrex::Math::isfinite(candidate);
            fallback_fraction = used ? Real(0.0) : Real(1.0);
        } else if constexpr (
            WallModel::pressure_closure ==
            ibm_pressure_closure_t::navier_stokes_noslip_viscous_compatibility) {
            shock_sensor = Real(-1.0);
            const Real candidate = ibm_normal_grad_surface<EO>(
                q, Cls::QPRES, disIM, n_valid, pressure_compat.dpdn);
            const bool used = pressure_compat.valid != 0 &&
                              amrex::Math::isfinite(pressure_compat.dpdn) &&
                              candidate > Real(0.0) &&
                              amrex::Math::isfinite(candidate);
            fallback_fraction = used ? Real(0.0) : Real(1.0);
        }
    } else {
        amrex::ignore_unused(q, disIM, n_valid, pressure_compat);
    }
}

} // namespace ibm_detail

/// Positivity-guarded variant for positive-definite primitives (P, T, Y).
/// The 2nd-order two-point form (x2^2 u1 - x1^2 u2)/(x2^2 - x1^2) goes
/// non-positive when u2/u1 > (x2/x1)^2 — a strong gradient between IP1 and IP2
/// (e.g. a shock crossing the stencil). Fall back to the 1st-order value u1,
/// which is a convex combination of fluid-cell values and therefore positive.
/// Mirrors the documented ghost positivity floor in extrapolate(); exact no-op
/// at EO=1 (base helper already returns u1). NOT for velocities (legitimately
/// signed).
template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_zero_grad_surface_pos(const QArr& q, int n, const DisArr& disIM, int n_valid)
{
    const Real v = ibm_zero_grad_surface<EO>(q, n, disIM, n_valid);
    return (v > Real(0.0) && amrex::Math::isfinite(v)) ? v : q(2, n);
}

template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_normal_grad_surface_pos(
    const QArr& q, int n, const DisArr& disIM, int n_valid,
    Real normal_derivative)
{
    const Real v = ibm_normal_grad_surface<EO>(
        q, n, disIM, n_valid, normal_derivative);
    return (v > Real(0.0) && amrex::Math::isfinite(v)) ? v : q(2, n);
}

template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_fluid_extrap_surface_pos(
    const QArr& q, int n, const DisArr& disIM, int n_valid)
{
    const Real v = ibm_fluid_extrap_surface<EO>(
        q, n, disIM, n_valid);
    return (v > Real(0.0) && amrex::Math::isfinite(v)) ? v : q(2, n);
}

template <int EO, typename QArr, typename DisArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_wall_normal_deriv(const QArr& q, int n, const DisArr& disIM, int n_valid)
{
    const Real phi_s = q(1, n);   // surface
    const Real phi1  = q(2, n);   // IP1
    const Real x1    = disIM(0);
    constexpr Real eps = Real(1.0e-12);
    if constexpr (EO >= 3) {
        // 3rd-order one-sided derivative on {0, x1, x2, x3} (exact for cubics;
        // coefficients are the Lagrange-basis derivatives at s = 0).
        if (n_valid >= 3) {
            const Real x2 = disIM(1);
            const Real x3 = disIM(2);
            const bool distinct =
                x1 > Real(0.0) &&
                amrex::Math::abs(x2 - x1) > eps * amrex::max(x1, x2) &&
                amrex::Math::abs(x3 - x2) > eps * amrex::max(x2, x3) &&
                amrex::Math::abs(x3 - x1) > eps * amrex::max(x1, x3);
            if (distinct) {
                const Real phi2 = q(3, n);
                const Real phi3 = q(4, n);
                const Real cs = -(Real(1.0)/x1 + Real(1.0)/x2 + Real(1.0)/x3);
                const Real c1 = x2 * x3 / (x1 * (x1 - x2) * (x1 - x3));
                const Real c2 = x1 * x3 / (x2 * (x2 - x1) * (x2 - x3));
                const Real c3 = x1 * x2 / (x3 * (x3 - x1) * (x3 - x2));
                return cs * phi_s + c1 * phi1 + c2 * phi2 + c3 * phi3;
            }
        }
    }
    if constexpr (EO >= 2) {
        if (n_valid >= 2) {
            const Real x2 = disIM(1);
            if (x1 > Real(0.0) && amrex::Math::abs(x2 - x1) > eps * amrex::max(x1, x2)) {
                const Real phi2 = q(3, n);   // IP2
                const Real cs = -(x1 + x2) / (x1 * x2);
                const Real c1 =  x2 / (x1 * (x2 - x1));
                const Real c2 = -x1 / (x2 * (x2 - x1));
                return cs * phi_s + c1 * phi1 + c2 * phi2;
            }
        }
    } else {
        amrex::ignore_unused(n_valid);
    }
    const Real inv = (x1 > Real(0.0)) ? Real(1.0) / x1 : Real(0.0);
    return (phi1 - phi_s) * inv;   // 1st-order one-sided
}

// ============================================================================
// 3. Constants
// ============================================================================

// Index dimension
static constexpr int IDIM = AMREX_SPACEDIM - 1;     

// Box extra width for ghost point search, no more than cls_t::NGHOST - 1
static constexpr int GP_BOX_EXTRA = 0;

// Maximum number of IB geometries (for GPU GpuArray captures)
static constexpr int MAX_NGEOM = 16;

// Minimum number of valid fluid points required in the interpolation stencil
//
// For the first image point, if the number of available interpolation points is
// less than INTERP_THRESHOLD, an error will be raised and the corresponding IB
// information will be reported.
//
// For all subsequent image points, if the number of interpolation points is
// insufficient, the image point will be discarded and the order of the
// extrapolation will be reduced accordingly.
// In 2D (4 points total), we require 2.
// In 3D (8 points total), we require 3.
#if (AMREX_SPACEDIM == 2)
static constexpr int INTERP_THRESHOLD_GP   = 2;
static constexpr int INTERP_THRESHOLD_SURF = 1; // for surface data reconstruction, less strict
#else
static constexpr int INTERP_THRESHOLD_GP   = 3;
static constexpr int INTERP_THRESHOLD_SURF = 2; // for surface data reconstruction, less strict
#endif

// Number of attempts for the first image point placement
static constexpr int N_ATTEMPTS_GP   = 3;  
static constexpr int N_ATTEMPTS_SURF = 5;   

// ============================================================================
// 4. surfImp_t — Surface image-point data (SoA)
// ============================================================================

template <int eorder_tparm_surf, int iorder_tparm_surf>
struct surfImp_t {
  // ideal number of interpolation points for each image point
  static constexpr int  N_InterP = ipow(iorder_tparm_surf + 1, AMREX_SPACEDIM);

  // Surface identification
  Gpu::ManagedVector<int> elemIdx;      // Global face index across all geometries

  // Image point data (per face)
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1, 0, IDIM>> imp_xyz;                       // Physical-space coordinates of image points placed along the outward normal
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm_surf - 1, 0, IDIM>> imp_ijk;                       // Index of the "bottom-left" grid cell associated with each image point
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm_surf - 1>> disIM;                                  // Normal distances from the surface (IB point) to each image point
  
  // Interpolation data for image points (per face)
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm_surf - 1>> imp_ninterp;                            // Actual number of interpolation points used for each image point
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm_surf - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;   // Indices of the 8-point interpolation stencil for each image point
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1, 0, N_InterP - 1>> imp_ipweights;         // Trilinear interpolation weights for the 8-point stencil of each image point
  // Optional finite-volume surface-observer data. These weights map full
  // Cartesian conservative cell averages to point values at each image point.
  // They are never used by the GP extension or the evolved full-cell RHS.
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm_surf - 1>>
      imp_ninterp_cell_average;
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm_surf - 1>>
      imp_fit_order_cell_average;
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm_surf - 1,
                              0, N_InterP - 1, 0, IDIM>>
      imp_ip_ijk_cell_average;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1,
                              0, N_InterP - 1>>
      imp_ipweights_cell_average;

  void resize(int n, bool allocate_cell_average = false) {
      elemIdx.resize(n);
      imp_xyz.resize(n);
      imp_ijk.resize(n);
      disIM.resize(n);
      imp_ninterp.resize(n);
      imp_ip_ijk.resize(n);
      imp_ipweights.resize(n);
      if (allocate_cell_average) {
          imp_ninterp_cell_average.resize(n);
          imp_fit_order_cell_average.resize(n);
          imp_ip_ijk_cell_average.resize(n);
          imp_ipweights_cell_average.resize(n);
      } else {
          imp_ninterp_cell_average.clear();
          imp_fit_order_cell_average.clear();
          imp_ip_ijk_cell_average.clear();
          imp_ipweights_cell_average.clear();
      }
  }

  void clear() {
      elemIdx.clear();       elemIdx.shrink_to_fit();
      imp_xyz.clear();       imp_xyz.shrink_to_fit();
      imp_ijk.clear();       imp_ijk.shrink_to_fit();
      disIM.clear();         disIM.shrink_to_fit();
      imp_ninterp.clear();   imp_ninterp.shrink_to_fit();
      imp_ip_ijk.clear();    imp_ip_ijk.shrink_to_fit();
      imp_ipweights.clear(); imp_ipweights.shrink_to_fit();
      imp_ninterp_cell_average.clear();
      imp_ninterp_cell_average.shrink_to_fit();
      imp_fit_order_cell_average.clear();
      imp_fit_order_cell_average.shrink_to_fit();
      imp_ip_ijk_cell_average.clear();
      imp_ip_ijk_cell_average.shrink_to_fit();
      imp_ipweights_cell_average.clear();
      imp_ipweights_cell_average.shrink_to_fit();
  }

  void shrink() {
      // Only shrink if capacity is significantly larger than size (e.g. > 4x)
      // to avoid frequent reallocations (memory jitter).
      if (elemIdx.capacity() > static_cast<std::size_t>(4 * elemIdx.size())) {
          elemIdx.shrink_to_fit();
          imp_xyz.shrink_to_fit();
          imp_ijk.shrink_to_fit();
          disIM.shrink_to_fit();
          imp_ninterp.shrink_to_fit();
          imp_ip_ijk.shrink_to_fit();
          imp_ipweights.shrink_to_fit();
          imp_ninterp_cell_average.shrink_to_fit();
          imp_fit_order_cell_average.shrink_to_fit();
          imp_ip_ijk_cell_average.shrink_to_fit();
          imp_ipweights_cell_average.shrink_to_fit();
      }
  }
};

// ============================================================================
// 7. surfPhys_t — Surface physical data (SoA)
// ============================================================================

struct surfPhys_t {

  // CPU only attributes
  surfPhys_t() : filled_elems(0) {}
  int filled_elems;

  // Indexing info
  Gpu::ManagedVector<int> elemIdx;      // Global face index across all geometries
  Gpu::ManagedVector<int> ifab;         // Local FAB index on this MPI rank
  Gpu::ManagedVector<int> lev;          // AMR level
  Gpu::ManagedVector<int> rank;         // Owning MPI rank
  Gpu::ManagedVector<int> elemfound;    // Whether this face has been located (int for GPU compatibility)
  Gpu::ManagedVector<int> ip_quality;   // number of fluid interpolation points used for the first image point
  Gpu::ManagedVector<int> ip_fit_order; // first-IP polynomial fit: 2/1/0, or -1 invalid
  //elemfound: 0 = outside this level, 1 = owned/valid, -1 = reserved by other rank

  // Surface fields (per face)
  Gpu::ManagedVector<Real> pressure;    // reconstructed local surface pressure
  Gpu::ManagedVector<Real> tau_normal;  // reconstructed normal viscous traction
  Gpu::ManagedVector<Real> tau1;        // reconstructed local surface shear stress 1
  Gpu::ManagedVector<Real> tau2;        // reconstructed local surface shear stress 2
  Gpu::ManagedVector<Real> temperature; // reconstructed temperature
  Gpu::ManagedVector<Real> dTdn;        // reconstructed grad(T)·n
  Gpu::ManagedVector<Real> pressure_shock_sensor; // irreversible IP1/IP2 compression
  Gpu::ManagedVector<Real> pressure_fallback;    // shock-aware return-to-IP1 fraction
  
  // Helper function to resize all vectors
  // Note: resize() initializes new elements to 0 / default constructor.
  // If you need specific default values (e.g. -1 for indices), set them manually after resize.
  void resize(int n) {
    int old_n = elemIdx.size();
    
    elemIdx.resize(n);
    ifab.resize(n);
    lev.resize(n);
    rank.resize(n);
    elemfound.resize(n);
    
    pressure.resize(n);
    tau_normal.resize(n);
    tau1.resize(n);
    tau2.resize(n);
    temperature.resize(n);
    dTdn.resize(n);
    pressure_shock_sensor.resize(n);
    pressure_fallback.resize(n);
    ip_quality.resize(n);
    ip_fit_order.resize(n);

    // ManagedVector resize/initialization may still be in flight on WSL2 GPUs.
    // The metadata defaults below are written on the host, so establish a
    // migration-safe synchronization point before touching managed storage.
    Gpu::streamSynchronize();

    // Initialize new elements with specific defaults if n > old_n
    if (n > old_n) {
        for (int i = old_n; i < n; ++i) {
            ifab[i] = -1;
            lev[i]  = -1;
            rank[i] = -1;
            elemfound[i] = 0; // false
            ip_quality[i] = -1;
            ip_fit_order[i] = -1;
            pressure_shock_sensor[i] = Real(-1.0);
            pressure_fallback[i] = Real(-1.0);
        }
    }
  }

  // Clear and free memory
  void clear() {
      filled_elems = 0;

      elemIdx.clear();     elemIdx.shrink_to_fit();
      ifab.clear();        ifab.shrink_to_fit();
      lev.clear();         lev.shrink_to_fit();
      rank.clear();        rank.shrink_to_fit();
      elemfound.clear();   elemfound.shrink_to_fit();

      pressure.clear();    pressure.shrink_to_fit();
      tau_normal.clear();  tau_normal.shrink_to_fit();
      tau1.clear();        tau1.shrink_to_fit();
      tau2.clear();        tau2.shrink_to_fit();
      temperature.clear(); temperature.shrink_to_fit();
      dTdn.clear();        dTdn.shrink_to_fit();
      pressure_shock_sensor.clear(); pressure_shock_sensor.shrink_to_fit();
      pressure_fallback.clear();    pressure_fallback.shrink_to_fit();
      ip_quality.clear();  ip_quality.shrink_to_fit();
      ip_fit_order.clear(); ip_fit_order.shrink_to_fit();
  }

  // Explicitly release memory
  void shrink() {
    if (elemIdx.capacity() > static_cast<std::size_t>(4 * elemIdx.size())) {
        elemIdx.shrink_to_fit();
        ifab.shrink_to_fit();
        lev.shrink_to_fit();
        rank.shrink_to_fit();
        elemfound.shrink_to_fit();
        
        pressure.shrink_to_fit();
        tau_normal.shrink_to_fit();
        tau1.shrink_to_fit();
        tau2.shrink_to_fit();
        temperature.shrink_to_fit();
        dTdn.shrink_to_fit();
        pressure_shock_sensor.shrink_to_fit();
        pressure_fallback.shrink_to_fit();
        ip_quality.shrink_to_fit();
        ip_fit_order.shrink_to_fit();
    }
  }

  // Reset metadata for regrid
  void reset() {
    // Surface arrays are managed storage and can have outstanding device work
    // from allocation or the preceding level.  Host writes require an explicit
    // synchronization on devices without concurrent managed access.
    Gpu::streamSynchronize();
    int n = elemIdx.size();
    for (int i = 0; i < n; ++i) {
        ifab[i] = -1;
        lev[i]  = -1;
        rank[i] = -1;
        elemfound[i] = 0; // false
        ip_quality[i] = -1;
        ip_fit_order[i] = -1;
        pressure_shock_sensor[i] = Real(-1.0);
        pressure_fallback[i] = Real(-1.0);
    }
  }

  AMREX_FORCE_INLINE
  bool owned(int f_idx, int lev, int ifab) const noexcept
  {
        return elemfound[f_idx] == 1 &&
                this->lev[f_idx]  == lev &&
                this->ifab[f_idx] == ifab;
  }

};

// ============================================================================
// 6. Type traits and helpers
// ============================================================================

// Detects the ghost-point storage view (has ``gp_ijk``) vs the surface storage
// (no ``gp_ijk``).  Used by interp templates in ibm_solver_interp.h to pick
// GP vs surface code paths at compile time.
//
// The historical name ``is_gpData_t`` predates the CSR refactor that removed
// the per-FAB gpData_t struct; it now matches GPStoreView (the live storage).
template <typename T, typename = void>
struct is_gpData_t : std::false_type {};

template <typename T>
struct is_gpData_t<T, std::void_t<decltype(std::declval<T>().gp_ijk)>> : std::true_type {};

// ============================================================================
// 7. GPStoreView / GPStore — Level-wide flattened ghost-point storage (CSR)
// ============================================================================

// Four tensor-product Gauss points integrate a smooth conservative state over
// one 2-D annular cell to the accuracy required by the quadratic BI-CWLS
// extension. Three adjacent radial target cells provide the same centre-state
// recovery stencil used by the fluid finite-volume data.
static constexpr int IBM_RZ_BIC_GAUSS_POINTS = 4;
static constexpr int IBM_RZ_BIC_RADIAL_STENCIL_CELLS = 3;
static constexpr int IBM_RZ_BIC_ANNULAR_TARGETS =
    IBM_RZ_BIC_RADIAL_STENCIL_CELLS * IBM_RZ_BIC_GAUSS_POINTS;
static constexpr int IBM_RZ_BIC_TARGETS = IBM_RZ_BIC_ANNULAR_TARGETS;

template <int N_InterP>
struct RZBICPointFunctionals {
  Array2D<Real, 0, IBM_RZ_BIC_TARGETS - 1, 0, N_InterP - 1>
      dirichlet_weights{};
  Array2D<Real, 0, IBM_RZ_BIC_TARGETS - 1, 0, N_InterP - 1>
      neumann_weights{};
  Array2D<Real, 0, IBM_RZ_BIC_TARGETS - 1, 0, N_InterP - 1>
      fv_dirichlet_weights{};
  Array2D<Real, 0, IBM_RZ_BIC_TARGETS - 1, 0, N_InterP - 1>
      entropy_jet_weights{};
  Array2D<Real, 0, IBM_RZ_BIC_TARGETS - 1, 0, N_InterP - 1>
      entropy_jet_linear_weights{};
  Array1D<Real, 0, IBM_RZ_BIC_TARGETS - 1> dirichlet_boundary_weight{};
  Array1D<Real, 0, IBM_RZ_BIC_TARGETS - 1> neumann_gradient_weight{};
  Array1D<Real, 0, IBM_RZ_BIC_TARGETS - 1>
      fv_dirichlet_boundary_weight{};
};

/// \brief GPU-capturable POD view into GPStore.  Holds raw pointers only.
///        This struct is trivially copyable and can be captured by GPU lambdas.
template <int eorder_tparm, int iorder_tparm, int ncons_tparm>
struct GPStoreView {
  static constexpr int N_InterP = ipow(iorder_tparm + 1, AMREX_SPACEDIM);

  int  total_ngps;                                   // total ghost points on this level
  int  nfabs;                                        // number of local FABs
  const int*  fab_offsets;                            // CSR offsets: fab_offsets[ifab] .. fab_offsets[ifab+1]

  // Per-GP arrays (indexed 0..total_ngps-1)
  const int*                    gp_fab;                   // local FAB index for each GP
  const Array1D< int, 0, IDIM>* gp_ijk;
  const Array1D<Real, 0, IDIM>* ib_xyz;
  const Real*                   disGP;
  const int*                    geomIdx;
  const int*                    elemIdx;
  const Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>* imp_xyz;
  const Array2D< int, 0, eorder_tparm - 1, 0, IDIM>* imp_ijk;
  const Array1D<Real, 0, eorder_tparm - 1>*           disIM;
  const Array1D< int, 0, eorder_tparm - 1>*           imp_ninterp;
  const Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>* imp_ip_ijk;
  const Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>*          imp_ipweights;
  // Optional BI-constrained shared-GP metadata.  The support block has the
  // same 3^D footprint as iorder=2, but its anchor is the reflected real GP;
  // no image-point value is formed.  The two functionals impose either a
  // Dirichlet value or a wall-normal derivative at the BI.
  const Array2D<int, 0, N_InterP - 1, 0, IDIM>* constrained_support_ijk;
  const Array1D<Real, 0, N_InterP - 1>* constrained_dirichlet_weights;
  const Array1D<Real, 0, N_InterP - 1>* constrained_neumann_weights;
  // Finite-volume Dirichlet functional: support and target rows are exact
  // Cartesian cell averages of the local polynomial basis.
  const Array1D<Real, 0, N_InterP - 1>*
      constrained_fv_dirichlet_weights;
  // Unconstrained one-sided fluid jet evaluated at the real GP centre.
  const Array1D<Real, 0, N_InterP - 1>* constrained_entropy_jet_weights;
  // Embedded linear jet on the identical visible support.  Runtime compares
  // this state with the quadratic jet to detect loss of smoothness without
  // changing the geometry, support set, or wall constraints.
  const Array1D<Real, 0, N_InterP - 1>*
      constrained_entropy_jet_linear_weights;
  // Unconstrained one-sided value at the boundary intercept.  This is used
  // only for compatibility data that depend on the current fluid-side trace;
  // it never replaces the constrained ghost value itself.
  const Array1D<Real, 0, N_InterP - 1>* constrained_fluid_trace_weights;
  const Real* constrained_dirichlet_boundary_weight;
  const Real* constrained_neumann_gradient_weight;
  const Real* constrained_fv_dirichlet_boundary_weight;
  const int* constrained_order;
  const int* constrained_fv_dirichlet_order;
  const int* constrained_entropy_jet_order;
  const int* constrained_fluid_trace_order;
  const int* constrained_support_count;
  const Real* constrained_condition;
  const Real* constrained_weight_l1;
  const Real* constrained_fv_dirichlet_condition;
  const Real* constrained_fv_dirichlet_weight_l1;
  const Real* constrained_entropy_jet_condition;
  const Real* constrained_entropy_jet_weight_l1;
  const Real* constrained_entropy_jet_linear_condition;
  const Real* constrained_entropy_jet_linear_weight_l1;
  const Real* constrained_fluid_trace_condition;
  const Real* constrained_fluid_trace_weight_l1;
  // R-Z only: point-evaluation functionals used to integrate a
  // thermodynamically consistent annular conservative ghost-cell average.
  const RZBICPointFunctionals<N_InterP>* rz_bic_point_functionals;
  const Array2D<Real, 0, IBM_RZ_BIC_RADIAL_STENCIL_CELLS - 1,
                0, ncons_tparm - 1>* rz_bic_annular_stencil;

  // Non-const data pointers for initialiseGPs (write pass)
  Array1D< int, 0, IDIM>* gp_ijk_w;
  Array1D<Real, 0, IDIM>* ib_xyz_w;
  Real*                   disGP_w;
  int*                    geomIdx_w;
  int*                    elemIdx_w;
  Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>* imp_xyz_w;
  Array2D< int, 0, eorder_tparm - 1, 0, IDIM>* imp_ijk_w;
  Array1D<Real, 0, eorder_tparm - 1>*           disIM_w;
  Array1D< int, 0, eorder_tparm - 1>*           imp_ninterp_w;
  Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>* imp_ip_ijk_w;
  Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>*          imp_ipweights_w;
  Array2D<int, 0, N_InterP - 1, 0, IDIM>* constrained_support_ijk_w;
  Array1D<Real, 0, N_InterP - 1>* constrained_dirichlet_weights_w;
  Array1D<Real, 0, N_InterP - 1>* constrained_neumann_weights_w;
  Array1D<Real, 0, N_InterP - 1>* constrained_fv_dirichlet_weights_w;
  Array1D<Real, 0, N_InterP - 1>* constrained_entropy_jet_weights_w;
  Array1D<Real, 0, N_InterP - 1>*
      constrained_entropy_jet_linear_weights_w;
  Array1D<Real, 0, N_InterP - 1>* constrained_fluid_trace_weights_w;
  Real* constrained_dirichlet_boundary_weight_w;
  Real* constrained_neumann_gradient_weight_w;
  Real* constrained_fv_dirichlet_boundary_weight_w;
  int* constrained_order_w;
  int* constrained_fv_dirichlet_order_w;
  int* constrained_entropy_jet_order_w;
  int* constrained_fluid_trace_order_w;
  int* constrained_support_count_w;
  Real* constrained_condition_w;
  Real* constrained_weight_l1_w;
  Real* constrained_fv_dirichlet_condition_w;
  Real* constrained_fv_dirichlet_weight_l1_w;
  Real* constrained_entropy_jet_condition_w;
  Real* constrained_entropy_jet_weight_l1_w;
  Real* constrained_entropy_jet_linear_condition_w;
  Real* constrained_entropy_jet_linear_weight_l1_w;
  Real* constrained_fluid_trace_condition_w;
  Real* constrained_fluid_trace_weight_l1_w;
  RZBICPointFunctionals<N_InterP>* rz_bic_point_functionals_w;
  Array2D<Real, 0, IBM_RZ_BIC_RADIAL_STENCIL_CELLS - 1,
          0, ncons_tparm - 1>* rz_bic_annular_stencil_w;

  /// Get GP range for a local FAB index
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int gp_begin(int ifab) const { return fab_offsets[ifab]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int gp_end(int ifab) const { return fab_offsets[ifab + 1]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int ngps(int ifab) const { return fab_offsets[ifab + 1] - fab_offsets[ifab]; }
};

/// \brief Level-wide flattened ghost-point storage.  Owns memory via ManagedVector.
///        CSR layout: fab_offsets[ifab] gives the first GP index for local fab ifab.
template <int eorder_tparm, int iorder_tparm, int ncons_tparm>
struct GPStore {
  static constexpr int N_InterP = ipow(iorder_tparm + 1, AMREX_SPACEDIM);

  int total_ngps = 0;   // total ghost points on this level (sum over all FABs)
  int nfabs      = 0;   // number of local FABs

  // CSR index: size = nfabs + 1.  fab_offsets[i] = start index of FAB i's GPs.
  Gpu::ManagedVector<int> fab_offsets;

  // Per-GP arrays (flattened, size = total_ngps each)
  Gpu::ManagedVector<int>                    gp_fab;   // gp_fab[ii] = local FAB index for GP ii
  Gpu::ManagedVector<Array1D< int, 0, IDIM>> gp_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, IDIM>> ib_xyz;
  Gpu::ManagedVector<Real>                   disGP;
  Gpu::ManagedVector<int>                    geomIdx;
  Gpu::ManagedVector<int>                    elemIdx;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>> imp_xyz;
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>> imp_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm - 1>>           disIM;
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm - 1>>           imp_ninterp;
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>>          imp_ipweights;
  Gpu::ManagedVector<Array2D<int, 0, N_InterP - 1, 0, IDIM>>
      constrained_support_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, N_InterP - 1>>
      constrained_dirichlet_weights;
  Gpu::ManagedVector<Array1D<Real, 0, N_InterP - 1>>
      constrained_neumann_weights;
  Gpu::ManagedVector<Array1D<Real, 0, N_InterP - 1>>
      constrained_fv_dirichlet_weights;
  Gpu::ManagedVector<Array1D<Real, 0, N_InterP - 1>>
      constrained_entropy_jet_weights;
  Gpu::ManagedVector<Array1D<Real, 0, N_InterP - 1>>
      constrained_entropy_jet_linear_weights;
  Gpu::ManagedVector<Array1D<Real, 0, N_InterP - 1>>
      constrained_fluid_trace_weights;
  Gpu::ManagedVector<Real> constrained_dirichlet_boundary_weight;
  Gpu::ManagedVector<Real> constrained_neumann_gradient_weight;
  Gpu::ManagedVector<Real> constrained_fv_dirichlet_boundary_weight;
  Gpu::ManagedVector<int> constrained_order;
  Gpu::ManagedVector<int> constrained_fv_dirichlet_order;
  Gpu::ManagedVector<int> constrained_entropy_jet_order;
  Gpu::ManagedVector<int> constrained_fluid_trace_order;
  Gpu::ManagedVector<int> constrained_support_count;
  Gpu::ManagedVector<Real> constrained_condition;
  Gpu::ManagedVector<Real> constrained_weight_l1;
  Gpu::ManagedVector<Real> constrained_fv_dirichlet_condition;
  Gpu::ManagedVector<Real> constrained_fv_dirichlet_weight_l1;
  Gpu::ManagedVector<Real> constrained_entropy_jet_condition;
  Gpu::ManagedVector<Real> constrained_entropy_jet_weight_l1;
  Gpu::ManagedVector<Real> constrained_entropy_jet_linear_condition;
  Gpu::ManagedVector<Real> constrained_entropy_jet_linear_weight_l1;
  Gpu::ManagedVector<Real> constrained_fluid_trace_condition;
  Gpu::ManagedVector<Real> constrained_fluid_trace_weight_l1;
  Gpu::ManagedVector<RZBICPointFunctionals<N_InterP>>
      rz_bic_point_functionals;
  Gpu::ManagedVector<
      Array2D<Real, 0, IBM_RZ_BIC_RADIAL_STENCIL_CELLS - 1,
              0, ncons_tparm - 1>>
      rz_bic_annular_stencil;

  /// Allocate flat arrays from per-fab counts.
  /// \param counts  Vector of GP counts per local FAB (size = nfabs_in).
  void allocate(const amrex::Vector<int>& counts) {
    nfabs = static_cast<int>(counts.size());

    // Build CSR offsets via exclusive prefix sum
    fab_offsets.resize(nfabs + 1);
    fab_offsets[0] = 0;
    for (int i = 0; i < nfabs; ++i) {
      fab_offsets[i + 1] = fab_offsets[i] + counts[i];
    }
    total_ngps = fab_offsets[nfabs];

    // Resize all per-GP arrays
    gp_fab.resize(total_ngps);
    gp_ijk.resize(total_ngps);

    // Fill gp_fab: expand CSR offsets to per-GP FAB index
    for (int f = 0; f < nfabs; ++f) {
      for (int g = fab_offsets[f]; g < fab_offsets[f + 1]; ++g) {
        gp_fab[g] = f;
      }
    }
    ib_xyz.resize(total_ngps);
    disGP.resize(total_ngps);
    geomIdx.resize(total_ngps);
    elemIdx.resize(total_ngps);
    imp_xyz.resize(total_ngps);
    imp_ijk.resize(total_ngps);
    disIM.resize(total_ngps);
    imp_ninterp.resize(total_ngps);
    imp_ip_ijk.resize(total_ngps);
    imp_ipweights.resize(total_ngps);
  }

  /// Allocate metadata used only by the opt-in BI-constrained shared-GP path.
  /// Keeping this separate avoids a permanent memory cost for legacy cases.
  void allocate_bi_constrained() {
    constrained_support_ijk.resize(total_ngps);
    constrained_dirichlet_weights.resize(total_ngps);
    constrained_neumann_weights.resize(total_ngps);
    constrained_fv_dirichlet_weights.resize(total_ngps);
    constrained_entropy_jet_weights.resize(total_ngps);
    constrained_entropy_jet_linear_weights.resize(total_ngps);
    constrained_fluid_trace_weights.resize(total_ngps);
    constrained_dirichlet_boundary_weight.resize(total_ngps);
    constrained_neumann_gradient_weight.resize(total_ngps);
    constrained_fv_dirichlet_boundary_weight.resize(total_ngps);
    constrained_order.resize(total_ngps);
    constrained_fv_dirichlet_order.resize(total_ngps);
    constrained_entropy_jet_order.resize(total_ngps);
    constrained_fluid_trace_order.resize(total_ngps);
    constrained_support_count.resize(total_ngps);
    constrained_condition.resize(total_ngps);
    constrained_weight_l1.resize(total_ngps);
    constrained_fv_dirichlet_condition.resize(total_ngps);
    constrained_fv_dirichlet_weight_l1.resize(total_ngps);
    constrained_entropy_jet_condition.resize(total_ngps);
    constrained_entropy_jet_weight_l1.resize(total_ngps);
    constrained_entropy_jet_linear_condition.resize(total_ngps);
    constrained_entropy_jet_linear_weight_l1.resize(total_ngps);
    constrained_fluid_trace_condition.resize(total_ngps);
    constrained_fluid_trace_weight_l1.resize(total_ngps);
  }

  /// Allocate the extra target functionals only for the opt-in R-Z annular
  /// shared-GP path. Cartesian BI-CWLS cases retain their previous footprint.
  void allocate_rz_bic_point_functionals() {
    rz_bic_point_functionals.resize(total_ngps);
    rz_bic_annular_stencil.resize(total_ngps);
  }

  /// Return a GPU-capturable read/write view of this store.
  GPStoreView<eorder_tparm, iorder_tparm, ncons_tparm> view() const {
    GPStoreView<eorder_tparm, iorder_tparm, ncons_tparm> v{};
    v.total_ngps   = total_ngps;
    v.nfabs        = nfabs;
    v.fab_offsets   = fab_offsets.data();

    v.gp_fab        = gp_fab.data();
    v.gp_ijk        = gp_ijk.data();
    v.ib_xyz        = ib_xyz.data();
    v.disGP         = disGP.data();
    v.geomIdx       = geomIdx.data();
    v.elemIdx       = elemIdx.data();
    v.imp_xyz       = imp_xyz.data();
    v.imp_ijk       = imp_ijk.data();
    v.disIM         = disIM.data();
    v.imp_ninterp   = imp_ninterp.data();
    v.imp_ip_ijk    = imp_ip_ijk.data();
    v.imp_ipweights = imp_ipweights.data();
    v.constrained_support_ijk = constrained_support_ijk.data();
    v.constrained_dirichlet_weights = constrained_dirichlet_weights.data();
    v.constrained_neumann_weights = constrained_neumann_weights.data();
    v.constrained_fv_dirichlet_weights =
        constrained_fv_dirichlet_weights.data();
    v.constrained_entropy_jet_weights =
        constrained_entropy_jet_weights.data();
    v.constrained_entropy_jet_linear_weights =
        constrained_entropy_jet_linear_weights.data();
    v.constrained_fluid_trace_weights =
        constrained_fluid_trace_weights.data();
    v.constrained_dirichlet_boundary_weight =
        constrained_dirichlet_boundary_weight.data();
    v.constrained_neumann_gradient_weight =
        constrained_neumann_gradient_weight.data();
    v.constrained_fv_dirichlet_boundary_weight =
        constrained_fv_dirichlet_boundary_weight.data();
    v.constrained_order = constrained_order.data();
    v.constrained_fv_dirichlet_order =
        constrained_fv_dirichlet_order.data();
    v.constrained_entropy_jet_order = constrained_entropy_jet_order.data();
    v.constrained_fluid_trace_order =
        constrained_fluid_trace_order.data();
    v.constrained_support_count = constrained_support_count.data();
    v.constrained_condition = constrained_condition.data();
    v.constrained_weight_l1 = constrained_weight_l1.data();
    v.constrained_fv_dirichlet_condition =
        constrained_fv_dirichlet_condition.data();
    v.constrained_fv_dirichlet_weight_l1 =
        constrained_fv_dirichlet_weight_l1.data();
    v.constrained_entropy_jet_condition =
        constrained_entropy_jet_condition.data();
    v.constrained_entropy_jet_weight_l1 =
        constrained_entropy_jet_weight_l1.data();
    v.constrained_entropy_jet_linear_condition =
        constrained_entropy_jet_linear_condition.data();
    v.constrained_entropy_jet_linear_weight_l1 =
        constrained_entropy_jet_linear_weight_l1.data();
    v.constrained_fluid_trace_condition =
        constrained_fluid_trace_condition.data();
    v.constrained_fluid_trace_weight_l1 =
        constrained_fluid_trace_weight_l1.data();
    v.rz_bic_point_functionals = rz_bic_point_functionals.data();
    v.rz_bic_annular_stencil = rz_bic_annular_stencil.data();

    // const_cast for writable pointers (initialiseGPs write pass)
    v.gp_ijk_w        = const_cast<Array1D< int, 0, IDIM>*>(gp_ijk.data());
    v.ib_xyz_w        = const_cast<Array1D<Real, 0, IDIM>*>(ib_xyz.data());
    v.disGP_w         = const_cast<Real*>(disGP.data());
    v.geomIdx_w       = const_cast<int*>(geomIdx.data());
    v.elemIdx_w       = const_cast<int*>(elemIdx.data());
    v.imp_xyz_w       = const_cast<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>*>(imp_xyz.data());
    v.imp_ijk_w       = const_cast<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>*>(imp_ijk.data());
    v.disIM_w         = const_cast<Array1D<Real, 0, eorder_tparm - 1>*>(disIM.data());
    v.imp_ninterp_w   = const_cast<Array1D< int, 0, eorder_tparm - 1>*>(imp_ninterp.data());
    v.imp_ip_ijk_w    = const_cast<Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>*>(imp_ip_ijk.data());
    v.imp_ipweights_w = const_cast<Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>*>(imp_ipweights.data());
    v.constrained_support_ijk_w = const_cast<
        Array2D<int, 0, N_InterP - 1, 0, IDIM>*>(
            constrained_support_ijk.data());
    v.constrained_dirichlet_weights_w = const_cast<
        Array1D<Real, 0, N_InterP - 1>*>(
            constrained_dirichlet_weights.data());
    v.constrained_neumann_weights_w = const_cast<
        Array1D<Real, 0, N_InterP - 1>*>(
            constrained_neumann_weights.data());
    v.constrained_fv_dirichlet_weights_w = const_cast<
        Array1D<Real, 0, N_InterP - 1>*>(
            constrained_fv_dirichlet_weights.data());
    v.constrained_entropy_jet_weights_w = const_cast<
        Array1D<Real, 0, N_InterP - 1>*>(
            constrained_entropy_jet_weights.data());
    v.constrained_entropy_jet_linear_weights_w = const_cast<
        Array1D<Real, 0, N_InterP - 1>*>(
            constrained_entropy_jet_linear_weights.data());
    v.constrained_fluid_trace_weights_w = const_cast<
        Array1D<Real, 0, N_InterP - 1>*>(
            constrained_fluid_trace_weights.data());
    v.constrained_dirichlet_boundary_weight_w =
        const_cast<Real*>(constrained_dirichlet_boundary_weight.data());
    v.constrained_neumann_gradient_weight_w =
        const_cast<Real*>(constrained_neumann_gradient_weight.data());
    v.constrained_fv_dirichlet_boundary_weight_w = const_cast<Real*>(
        constrained_fv_dirichlet_boundary_weight.data());
    v.constrained_order_w = const_cast<int*>(constrained_order.data());
    v.constrained_fv_dirichlet_order_w =
        const_cast<int*>(constrained_fv_dirichlet_order.data());
    v.constrained_entropy_jet_order_w =
        const_cast<int*>(constrained_entropy_jet_order.data());
    v.constrained_fluid_trace_order_w =
        const_cast<int*>(constrained_fluid_trace_order.data());
    v.constrained_support_count_w =
        const_cast<int*>(constrained_support_count.data());
    v.constrained_condition_w = const_cast<Real*>(constrained_condition.data());
    v.constrained_weight_l1_w = const_cast<Real*>(constrained_weight_l1.data());
    v.constrained_fv_dirichlet_condition_w =
        const_cast<Real*>(constrained_fv_dirichlet_condition.data());
    v.constrained_fv_dirichlet_weight_l1_w =
        const_cast<Real*>(constrained_fv_dirichlet_weight_l1.data());
    v.constrained_entropy_jet_condition_w =
        const_cast<Real*>(constrained_entropy_jet_condition.data());
    v.constrained_entropy_jet_weight_l1_w =
        const_cast<Real*>(constrained_entropy_jet_weight_l1.data());
    v.constrained_entropy_jet_linear_condition_w =
        const_cast<Real*>(constrained_entropy_jet_linear_condition.data());
    v.constrained_entropy_jet_linear_weight_l1_w =
        const_cast<Real*>(constrained_entropy_jet_linear_weight_l1.data());
    v.constrained_fluid_trace_condition_w =
        const_cast<Real*>(constrained_fluid_trace_condition.data());
    v.constrained_fluid_trace_weight_l1_w =
        const_cast<Real*>(constrained_fluid_trace_weight_l1.data());
    v.rz_bic_point_functionals_w = const_cast<
        RZBICPointFunctionals<N_InterP>*>(rz_bic_point_functionals.data());
    v.rz_bic_annular_stencil_w = const_cast<
        Array2D<Real, 0, IBM_RZ_BIC_RADIAL_STENCIL_CELLS - 1,
                0, ncons_tparm - 1>*>(rz_bic_annular_stencil.data());

    return v;
  }

  void clear() {
    total_ngps = 0;
    nfabs      = 0;
    fab_offsets.clear();
    gp_fab.clear();
    gp_ijk.clear();
    ib_xyz.clear();
    disGP.clear();
    geomIdx.clear();
    elemIdx.clear();
    imp_xyz.clear();
    imp_ijk.clear();
    disIM.clear();
    imp_ninterp.clear();
    imp_ip_ijk.clear();
    imp_ipweights.clear();
    constrained_support_ijk.clear();
    constrained_dirichlet_weights.clear();
    constrained_neumann_weights.clear();
    constrained_fv_dirichlet_weights.clear();
    constrained_entropy_jet_weights.clear();
    constrained_entropy_jet_linear_weights.clear();
    constrained_fluid_trace_weights.clear();
    constrained_dirichlet_boundary_weight.clear();
    constrained_neumann_gradient_weight.clear();
    constrained_fv_dirichlet_boundary_weight.clear();
    constrained_order.clear();
    constrained_fv_dirichlet_order.clear();
    constrained_entropy_jet_order.clear();
    constrained_fluid_trace_order.clear();
    constrained_support_count.clear();
    constrained_condition.clear();
    constrained_weight_l1.clear();
    constrained_fv_dirichlet_condition.clear();
    constrained_fv_dirichlet_weight_l1.clear();
    constrained_entropy_jet_condition.clear();
    constrained_entropy_jet_weight_l1.clear();
    constrained_entropy_jet_linear_condition.clear();
    constrained_entropy_jet_linear_weight_l1.clear();
    constrained_fluid_trace_condition.clear();
    constrained_fluid_trace_weight_l1.clear();
    rz_bic_point_functionals.clear();
    rz_bic_annular_stencil.clear();
  }

  void shrink() {
    auto cap = gp_ijk.capacity();
    if (cap > static_cast<std::size_t>(4 * total_ngps) && cap > 1000u) {
      fab_offsets.shrink_to_fit();
      gp_fab.shrink_to_fit();
      gp_ijk.shrink_to_fit();
      ib_xyz.shrink_to_fit();
      disGP.shrink_to_fit();
      geomIdx.shrink_to_fit();
      elemIdx.shrink_to_fit();
      imp_xyz.shrink_to_fit();
      imp_ijk.shrink_to_fit();
      disIM.shrink_to_fit();
      imp_ninterp.shrink_to_fit();
      imp_ip_ijk.shrink_to_fit();
      imp_ipweights.shrink_to_fit();
      constrained_support_ijk.shrink_to_fit();
      constrained_dirichlet_weights.shrink_to_fit();
      constrained_neumann_weights.shrink_to_fit();
      constrained_fv_dirichlet_weights.shrink_to_fit();
      constrained_entropy_jet_weights.shrink_to_fit();
      constrained_entropy_jet_linear_weights.shrink_to_fit();
      constrained_fluid_trace_weights.shrink_to_fit();
      constrained_dirichlet_boundary_weight.shrink_to_fit();
      constrained_neumann_gradient_weight.shrink_to_fit();
      constrained_fv_dirichlet_boundary_weight.shrink_to_fit();
      constrained_order.shrink_to_fit();
      constrained_fv_dirichlet_order.shrink_to_fit();
      constrained_entropy_jet_order.shrink_to_fit();
      constrained_fluid_trace_order.shrink_to_fit();
      constrained_support_count.shrink_to_fit();
      constrained_condition.shrink_to_fit();
      constrained_weight_l1.shrink_to_fit();
      constrained_fv_dirichlet_condition.shrink_to_fit();
      constrained_fv_dirichlet_weight_l1.shrink_to_fit();
      constrained_entropy_jet_condition.shrink_to_fit();
      constrained_entropy_jet_weight_l1.shrink_to_fit();
      constrained_entropy_jet_linear_condition.shrink_to_fit();
      constrained_entropy_jet_linear_weight_l1.shrink_to_fit();
      constrained_fluid_trace_condition.shrink_to_fit();
      constrained_fluid_trace_weight_l1.shrink_to_fit();
      rz_bic_point_functionals.shrink_to_fit();
      rz_bic_annular_stencil.shrink_to_fit();
    }
  }
};

// ============================================================================
// 9. FaceCSR — CSR structure for per-FAB face iteration
// ============================================================================

struct FaceCSR {
  Gpu::ManagedVector<int> fab_offsets;   // Offsets for each FAB (size = nfab + 1)
  Gpu::ManagedVector<int> face_indices;  // Contiguous array of face indices

  void clear() {
        fab_offsets.clear();  
        face_indices.clear();
  }

  void shrink() {
        fab_offsets.shrink_to_fit();
        face_indices.shrink_to_fit();
  }
};

// ============================================================================
// 10. CheckMode — Interpolation stencil check policy
// ============================================================================

enum class CheckMode {
    Silent,      // Do not output anything, just return status
    Warn,        // Output a warning message, return status
    Abort        // Abort execution immediately on failure
};

#endif // IBM_CONTAINERS_H_
