#include <AMReX_FluxRegister.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>
#include <CNS.h>
#include <prob.h>
#include <RZFiniteVolume.h>
#include <limits>  // for std::numeric_limits (RZ divergence floor)

// IBM and EB marker paths use different array types; combined mode is not yet supported.
#if defined(AMREX_USE_GPIBM) && defined(CNS_USE_EB)
#error "compute_rhs.cpp expects exactly one of AMREX_USE_GPIBM or CNS_USE_EB to be enabled"
#endif

#ifdef AMREX_USE_GPIBM
#include <ibm_solver.h>
#endif

using namespace amrex;

namespace {

using FaceFluxArray = std::array<FArrayBox*, AMREX_SPACEDIM>;
using FaceFallbackMaskArray =
    std::array<Array4<int>, AMREX_SPACEDIM>;

template <typename T, typename = void>
struct local_llf_fallback_mask_capable : std::false_type {};
template <typename T>
struct local_llf_fallback_mask_capable<
    T, std::void_t<decltype(T::local_llf_fallback_mask_capable)>>
    : std::bool_constant<T::local_llf_fallback_mask_capable> {};

// Every Euler operator uses the standard cylindrical +p/r radial-momentum
// source unless it explicitly opts out.  no_euler_t opts out so a
// diffusion-only R-Z RHS is not contaminated by a pressure source without a
// matching pressure face flux.
template <typename T, typename = void>
struct rz_euler_geometric_source_active : std::true_type {};
template <typename T>
struct rz_euler_geometric_source_active<
    T, std::void_t<decltype(T::rz_euler_geometric_source_active)>>
    : std::bool_constant<T::rz_euler_geometric_source_active> {};

template <typename T, typename = void>
struct rz_annular_pressure_consistency_capable : std::false_type {};
template <typename T>
struct rz_annular_pressure_consistency_capable<
    T, std::void_t<decltype(T::rz_annular_pressure_consistency_capable)>>
    : std::bool_constant<T::rz_annular_pressure_consistency_capable> {};

template <typename T, typename = void>
struct rz_radial_axis_face_flux_is_metric : std::false_type {};
template <typename T>
struct rz_radial_axis_face_flux_is_metric<
    T, std::void_t<decltype(T::rz_radial_axis_face_flux_is_metric)>>
    : std::bool_constant<T::rz_radial_axis_face_flux_is_metric> {};

// Most paired-pressure operators impose zero radial-momentum advection at
// r=0.  An operator may instead provide an auxiliary metric flux in the axis
// face slot, for example to close a parity-aware low-order rescue at the first
// interior radial face.  This is deliberately a separate capability from the
// generic axis-metric storage trait so existing WENO/AFD paths are unchanged.
template <typename T, typename = void>
struct rz_paired_axis_advective_flux_is_metric : std::false_type {};
template <typename T>
struct rz_paired_axis_advective_flux_is_metric<
    T,
    std::void_t<decltype(T::rz_paired_axis_advective_flux_is_metric)>>
    : std::bool_constant<T::rz_paired_axis_advective_flux_is_metric> {};

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

template <typename RhsT, typename PrimitiveView, typename StateView,
          typename ClosureT>
void compute_all_fluid_euler_fluxes(
    RhsT& rhs_operator, const Geometry& geom, const MFIter& mfi,
    const PrimitiveView& prims, const FaceFluxArray& face_fluxes,
    const StateView& state, const ClosureT* closure,
    const Array4<const Real>& pressure_reconstruction_states,
    const int pressure_reconstruction_component,
    const Array4<Real>& radial_pressure_face_flux,
    const bool capture_llf_fallback_mask,
    const FaceFallbackMaskArray& llf_fallback_mask)
{
  if (capture_llf_fallback_mask) {
    if constexpr (local_llf_fallback_mask_capable<RhsT>::value) {
      if (geom.IsRZ()) {
        rhs_operator
            .eflux_with_rz_paired_pressure_and_local_llf_fallback_mask(
                geom, mfi, prims, face_fluxes, state, closure,
                pressure_reconstruction_states,
                pressure_reconstruction_component,
                radial_pressure_face_flux, llf_fallback_mask);
      } else {
        rhs_operator.eflux_with_local_llf_fallback_mask(
            geom, mfi, prims, face_fluxes, state, closure,
            llf_fallback_mask);
      }
      return;
    } else {
      amrex::Abort(
          "the selected Euler operator cannot record an LLF fallback mask");
    }
  }

  if constexpr (rz_paired_pressure_flux_capable<RhsT>::value) {
    if (geom.IsRZ()) {
      rhs_operator.eflux_with_rz_paired_pressure(
          geom, mfi, prims, face_fluxes, state, closure,
          pressure_reconstruction_states,
          pressure_reconstruction_component,
          radial_pressure_face_flux);
      return;
    }
  }

  rhs_operator.eflux(
      geom, mfi, prims, face_fluxes, state, closure);
}

#ifdef AMREX_USE_GPIBM
template <typename RhsT, typename PrimitiveView, typename StateView,
          typename ClosureT, typename MarkerView>
void compute_shared_gp_euler_fluxes(
    RhsT& rhs_operator, const Geometry& geom, const MFIter& mfi,
    const PrimitiveView& prims, const FaceFluxArray& face_fluxes,
    const StateView& state, const ClosureT* closure,
    const MarkerView& markers,
    const Array4<const Real>& pressure_reconstruction_states,
    const int pressure_reconstruction_component,
    const Array4<Real>& radial_pressure_face_flux,
    const bool capture_llf_fallback_mask,
    const FaceFallbackMaskArray& llf_fallback_mask)
{
  if (capture_llf_fallback_mask) {
    if constexpr (local_llf_fallback_mask_capable<RhsT>::value) {
      if (geom.IsRZ()) {
        rhs_operator
            .eflux_ibm_with_rz_paired_pressure_and_local_llf_fallback_mask(
                geom, mfi, prims, face_fluxes, state, closure, markers,
                pressure_reconstruction_states,
                pressure_reconstruction_component,
                radial_pressure_face_flux, llf_fallback_mask);
      } else {
        rhs_operator.eflux_ibm_with_local_llf_fallback_mask(
            geom, mfi, prims, face_fluxes, state, closure, markers,
            llf_fallback_mask);
      }
      return;
    } else {
      amrex::Abort(
          "the selected GP Euler operator cannot record an LLF fallback mask");
    }
  }

  if constexpr (rz_paired_pressure_flux_capable<RhsT>::value) {
    if (geom.IsRZ()) {
      rhs_operator.eflux_ibm_with_rz_paired_pressure(
          geom, mfi, prims, face_fluxes, state, closure, markers,
          pressure_reconstruction_states,
          pressure_reconstruction_component,
          radial_pressure_face_flux);
      return;
    }
  }
  rhs_operator.eflux_ibm(
      geom, mfi, prims, face_fluxes, state, closure, markers);
}
#endif

template <typename RhsT, typename PrimArrayT>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
dispatch_rz_annular_pressure_from_cell_average(
    const IntVect& cell, const PrimArrayT& prims, const Real radial,
    const Real dr, const int pressure_component) noexcept
{
  if constexpr (rz_annular_pressure_consistency_capable<RhsT>::value) {
    return RhsT::rz_annular_pressure_from_cell_average(
        cell, prims, radial, dr);
  } else {
    return prims(cell, pressure_component);
  }
}

template <typename PrimArrayT, typename RhsArrayT, typename ClosuresT,
          typename ProbParmT>
auto try_rhs_nscbc(int, const Geometry& geom, const MFIter& mfi,
                   PrimArrayT const& prims, RhsArrayT const& rhs,
                   const ClosuresT* closures, const ProbParmT* pparm,
                   const Real dt, const Real time,
                   const bool forbid_active_rhs_override)
    -> decltype(rhs_nscbc(geom, mfi, prims, rhs, closures, pparm, dt, time),
                void())
{
  if (forbid_active_rhs_override) {
    amrex::Abort(
        "the pure-GP shared-face positivity limiter cannot be combined with "
        "a problem rhs_nscbc hook that may overwrite active-cell RHS; use "
        "the certified stage-local ghost-cell NSCBC path instead");
  }
  rhs_nscbc(geom, mfi, prims, rhs, closures, pparm, dt, time);
}

template <typename PrimArrayT, typename RhsArrayT, typename ClosuresT,
          typename ProbParmT>
void try_rhs_nscbc(long, const Geometry&, const MFIter&,
                   PrimArrayT const&, RhsArrayT const&,
                   const ClosuresT*, const ProbParmT*,
                   const Real, const Real, const bool)
{
}

// The numerical schemes fill face_fluxes. These helpers only assemble their
// conservative flux difference into rhs.
void assemble_cartesian_flux_divergence(
    const Geometry& geom, const Box& cell_box,
    const FaceFluxArray& face_fluxes, const Array4<Real>& rhs,
    const int ncons)
{
  const auto dx = geom.CellSizeArray();
  auto const& fx = face_fluxes[0]->array();
#if (AMREX_SPACEDIM >= 2)
  auto const& fy = face_fluxes[1]->array();
#endif
#if (AMREX_SPACEDIM == 3)
  auto const& fz = face_fluxes[2]->array();
#endif

#if (AMREX_SPACEDIM == 1)
  const Real invdx = Real(1.0) / dx[0];
  ParallelFor(
      cell_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int n = 0; n < ncons; ++n) {
          rhs(i, j, k, n) +=
              (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
        }
      });
#elif (AMREX_SPACEDIM == 2)
  const Real invdx = Real(1.0) / dx[0];
  const Real invdy = Real(1.0) / dx[1];
  ParallelFor(
      cell_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int n = 0; n < ncons; ++n) {
          rhs(i, j, k, n) +=
              (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
          rhs(i, j, k, n) +=
              (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
        }
      });
#else
  const Real invdx = Real(1.0) / dx[0];
  const Real invdy = Real(1.0) / dx[1];
  const Real invdz = Real(1.0) / dx[2];
  ParallelFor(
      cell_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int n = 0; n < ncons; ++n) {
          rhs(i, j, k, n) +=
              (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
          rhs(i, j, k, n) +=
              (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
          rhs(i, j, k, n) +=
              (fz(i, j, k, n) - fz(i, j, k + 1, n)) * invdz;
        }
      });
#endif
}

#if (AMREX_SPACEDIM == 2)
template <typename RhsT>
void assemble_rz_flux_divergence(
    const Geometry& geom, const Box& cell_box,
    const Array4<Real>& prims, const FaceFluxArray& face_fluxes,
    const Array4<const Real>& rz_center_pressure,
    const bool use_cell_average_deconvolution,
    const Array4<Real>& rz_pressure_face_flux,
    const bool use_rz_paired_pressure,
    const Array4<Real>& rhs, const int ncons)
{
  const auto dx = geom.CellSizeArray();
  const auto prob_lo = geom.ProbLoArray();
  auto const& radial_flux = face_fluxes[0]->array();
  auto const& axial_flux = face_fluxes[1]->array();

  const Real dr = dx[0];
  const Real dz = dx[1];
  const Real radial_origin = prob_lo[0];
  const int pressure_component = PROB::ProbClosures::QPRES;
  const int radial_momentum_component = PROB::ProbClosures::UMX;
  const Real inverse_dz = Real(1.0) / dz;
  const Real inverse_dr = Real(1.0) / dr;

  static const int s_annular_pressure_consistency = [] {
    int value = 0;
    ParmParse pp("cns");
    pp.query("rz_euler_annular_pressure_consistency", value);
    return value;
  }();
  static const int s_euler_point_flux = [] {
    int value = 0;
    ParmParse pp("cns");
    pp.query("rz_euler_point_flux", value);
    return value;
  }();

  const bool requested_annular_pressure_consistency =
      s_annular_pressure_consistency != 0;
  const bool use_annular_pressure_consistency =
      requested_annular_pressure_consistency && !use_rz_paired_pressure;
  constexpr bool euler_geometric_source_active =
      rz_euler_geometric_source_active<RhsT>::value;
  constexpr bool axis_face_flux_is_metric =
      rz_radial_axis_face_flux_is_metric<RhsT>::value;
  constexpr bool paired_axis_advective_flux_is_metric =
      rz_paired_axis_advective_flux_is_metric<RhsT>::value;

  if (use_annular_pressure_consistency && s_euler_point_flux == 0) {
    amrex::Abort(
        "cns.rz_euler_annular_pressure_consistency=1 requires explicit "
        "cns.rz_euler_point_flux=1");
  }
  if (use_annular_pressure_consistency &&
      !rz_annular_pressure_consistency_capable<RhsT>::value) {
    amrex::Abort(
        "cns.rz_euler_annular_pressure_consistency requires a compatible "
        "R-Z WENO/TENO Euler flux");
  }
  if (use_cell_average_deconvolution && !use_rz_paired_pressure &&
      !use_annular_pressure_consistency) {
    amrex::Abort(
        "cns.rz_euler_cell_average_deconvolution requires "
        "cns.rz_euler_annular_pressure_consistency=1");
  }

  ParallelFor(
      cell_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real r_lo = radial_origin + Real(i) * dr;
        const Real r_hi = r_lo + dr;
        const Real r_center = r_lo + Real(0.5) * dr;
        const Real radial_volume_factor = (r_hi + r_lo) * dr;
        const Real radial_volume_floor =
            amrex::max(std::numeric_limits<Real>::min(),
                       Real(1.0e-14) * dr * dr);
        const Real inverse_radial_volume =
            Real(1.0) /
            amrex::max(radial_volume_factor, radial_volume_floor);
        const Real inverse_radius =
            Real(1.0) /
            amrex::max(r_center, Real(1.0e-14) * dr);

        for (int n = 0; n < ncons; ++n) {
          if (n == radial_momentum_component &&
              use_rz_paired_pressure) {
            const Real pressure_lo = rz_pressure_face_flux(i, j, k, 0);
            const Real pressure_hi =
                rz_pressure_face_flux(i + 1, j, k, 0);
            const Real metric_advective_lo =
                r_lo > Real(0.0)
                    ? r_lo *
                          (radial_flux(i, j, k, n) - pressure_lo)
                    : (paired_axis_advective_flux_is_metric
                           ? radial_flux(i, j, k, n)
                           : Real(0.0));
            const Real metric_advective_hi =
                r_hi *
                (radial_flux(i + 1, j, k, n) - pressure_hi);
            rhs(i, j, k, n) +=
                Real(2.0) *
                    (metric_advective_lo - metric_advective_hi) *
                    inverse_radial_volume +
                (pressure_lo - pressure_hi) * inverse_dr;
          } else {
            const Real radial_low_metric_flux =
                axis_face_flux_is_metric && !(r_lo > Real(0.0))
                    ? radial_flux(i, j, k, n)
                    : r_lo * radial_flux(i, j, k, n);
            rhs(i, j, k, n) +=
                Real(2.0) *
                (radial_low_metric_flux -
                 r_hi * radial_flux(i + 1, j, k, n)) *
                inverse_radial_volume;
          }

          if (n == radial_momentum_component &&
              euler_geometric_source_active &&
              !use_rz_paired_pressure) {
            if (use_annular_pressure_consistency) {
              Real pressure_center;
              Real pressure_hi;
              Real pressure_lo;
              if (use_cell_average_deconvolution) {
                pressure_center = rz_center_pressure(i, j, k);
                pressure_hi = rz_center_pressure(i + 1, j, k);
                pressure_lo = rz_center_pressure(i - 1, j, k);
              } else {
                const IntVect cell(AMREX_D_DECL(i, j, k));
                const IntVect radial_cell =
                    IntVect::TheDimensionVector(0);
                pressure_center =
                    dispatch_rz_annular_pressure_from_cell_average<RhsT>(
                        cell, prims, r_center, dr, pressure_component);
                pressure_hi =
                    dispatch_rz_annular_pressure_from_cell_average<RhsT>(
                        cell + radial_cell, prims, r_center + dr, dr,
                        pressure_component);
                pressure_lo =
                    dispatch_rz_annular_pressure_from_cell_average<RhsT>(
                        cell - radial_cell, prims, r_center - dr, dr,
                        pressure_component);
              }
              const Real pressure_source_correction =
                  (pressure_hi - Real(2.0) * pressure_center + pressure_lo) *
                  (inverse_radius / Real(24.0));
              rhs(i, j, k, n) +=
                  pressure_center * inverse_radius +
                  pressure_source_correction;
            } else {
              rhs(i, j, k, n) +=
                  prims(i, j, k, pressure_component) * inverse_radius;
            }
          }

          rhs(i, j, k, n) +=
              (axial_flux(i, j, k, n) -
               axial_flux(i, j + 1, k, n)) *
              inverse_dz;
        }
      });
}
#endif

template <typename RhsT>
void assemble_face_flux_divergence(
    const Geometry& geom, const Box& cell_box,
    const Array4<Real>& prims, const FaceFluxArray& face_fluxes,
    const Array4<const Real>& rz_center_pressure,
    const bool use_cell_average_deconvolution,
    const Array4<Real>& rz_pressure_face_flux,
    const bool use_rz_paired_pressure,
    const Array4<Real>& rhs, const int ncons)
{
#if (AMREX_SPACEDIM == 2)
  if (geom.IsRZ()) {
    assemble_rz_flux_divergence<RhsT>(
        geom, cell_box, prims, face_fluxes, rz_center_pressure,
        use_cell_average_deconvolution, rz_pressure_face_flux,
        use_rz_paired_pressure, rhs, ncons);
    return;
  }
#else
  amrex::ignore_unused(
      prims, rz_center_pressure, use_cell_average_deconvolution,
      rz_pressure_face_flux, use_rz_paired_pressure);
#endif

  assemble_cartesian_flux_divergence(
      geom, cell_box, face_fluxes, rhs, ncons);
}

}  // namespace

// Assemble one Runge--Kutta stage RHS. Device kernels are launched once per FAB
// from a single untiled MFIter loop.

void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxRegister* fr_as_crse,
                      FluxRegister* fr_as_fine, Real stage_time,
                      Real reflux_dt,
                      std::array<MultiFab*, AMREX_SPACEDIM>
                          captured_face_flux,
                      MultiFab* captured_rz_pressure_face_flux,
                      std::array<iMultiFab*, AMREX_SPACEDIM>
                          captured_llf_fallback_mask) {
  BL_PROFILE("CNS::compute_rhs()");

  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  const PROB::ProbParm* pparm_d = CNS::d_prob_parm;

  // Reflux must integrate the actual numerical face flux with the final RK
  // quadrature weight, which is generally not the same as the local substep
  // scale passed as dt (SSPRK33 is the simplest counterexample).  Keep this
  // allocation entirely off the default do_reflux=0 path.
  const bool register_fluxes =
      reflux_dt != Real(0.0) && (fr_as_crse != nullptr || fr_as_fine != nullptr);
  bool capture_stage_flux = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    capture_stage_flux = capture_stage_flux || captured_face_flux[dir] != nullptr;
  }
  if (capture_stage_flux) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (captured_face_flux[dir] == nullptr) {
        amrex::Abort(
            "compute_rhs face-flux capture requires every direction");
      }
      const BoxArray expected_boxes = amrex::convert(
          statemf.boxArray(), IntVect::TheDimensionVector(dir));
      if (captured_face_flux[dir]->boxArray() != expected_boxes ||
          captured_face_flux[dir]->DistributionMap() !=
              statemf.DistributionMap() ||
          captured_face_flux[dir]->nComp() != cls_h.NCONS) {
        amrex::Abort("compute_rhs face-flux capture layout mismatch");
      }
    }
  }
  if (captured_rz_pressure_face_flux != nullptr) {
#if (AMREX_SPACEDIM == 2)
    const BoxArray expected_boxes = amrex::convert(
        statemf.boxArray(), IntVect::TheDimensionVector(0));
    if (!geom.IsRZ() ||
        captured_rz_pressure_face_flux->boxArray() != expected_boxes ||
        captured_rz_pressure_face_flux->DistributionMap() !=
            statemf.DistributionMap() ||
        captured_rz_pressure_face_flux->nComp() != 1) {
      amrex::Abort(
          "compute_rhs R-Z pressure-face capture layout mismatch");
    }
#else
    amrex::Abort(
        "compute_rhs R-Z pressure-face capture requires a 2-D build");
#endif
  }
  const bool capture_rz_pressure =
      captured_rz_pressure_face_flux != nullptr;
  bool capture_llf_fallback = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    capture_llf_fallback = capture_llf_fallback ||
                           captured_llf_fallback_mask[dir] != nullptr;
  }
  if (capture_llf_fallback) {
    if constexpr (!local_llf_fallback_mask_capable<PROB::ProbRHS>::value) {
      amrex::Abort(
          "cns.llf_fallback_mask_diagnostics requires a WENO/TENO LLF "
          "Euler operator");
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (captured_llf_fallback_mask[dir] == nullptr) {
        amrex::Abort(
            "LLF fallback-mask capture requires every coordinate direction");
      }
      const BoxArray expected_boxes = amrex::convert(
          statemf.boxArray(), IntVect::TheDimensionVector(dir));
      if (captured_llf_fallback_mask[dir]->boxArray() != expected_boxes ||
          captured_llf_fallback_mask[dir]->DistributionMap() !=
              statemf.DistributionMap() ||
          captured_llf_fallback_mask[dir]->nComp() != 1) {
        amrex::Abort("compute_rhs LLF fallback-mask layout mismatch");
      }
    }
  }
  const bool retain_level_face_fluxes = register_fluxes || capture_stage_flux;
  std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> level_face_fluxes;
  std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> level_face_areas;
  if (retain_level_face_fluxes) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      BoxArray face_boxes = amrex::convert(
          statemf.boxArray(), IntVect::TheDimensionVector(dir));
      level_face_fluxes[dir] = std::make_unique<MultiFab>(
          face_boxes, statemf.DistributionMap(), cls_h.NCONS, 0,
          MFInfo().SetArena(The_Async_Arena()));
      if (register_fluxes) {
        level_face_areas[dir] = std::make_unique<MultiFab>(
            face_boxes, statemf.DistributionMap(), 1, 0,
            MFInfo().SetArena(The_Async_Arena()));
        geom.GetFaceArea(*level_face_areas[dir], dir);
      }
    }
  }
  std::unique_ptr<MultiFab> level_rz_pressure_face_flux;
  if (capture_rz_pressure) {
    const BoxArray radial_face_boxes = amrex::convert(
        statemf.boxArray(), IntVect::TheDimensionVector(0));
    level_rz_pressure_face_flux = std::make_unique<MultiFab>(
        radial_face_boxes, statemf.DistributionMap(), 1, 0,
        MFInfo().SetArena(The_Async_Arena()));
  }

  // RK stage time is explicit.  StateData::curTime() is not reliable for the
  // first stage after swapTimeLevels(): it already points at t^{n+1} while
  // FillPatch reads U^n at t^n.
  const Real cur_time = stage_time;

  static const int s_rz_euler_cell_average_deconvolution = [] {
    int value = 0;
    ParmParse pp("cns");
    pp.query("rz_euler_cell_average_deconvolution", value);
    return value;
  }();
  const bool rz_euler_cell_average_deconvolution =
      s_rz_euler_cell_average_deconvolution != 0;
#if (AMREX_SPACEDIM == 2) && !defined(CNS_USE_EB)
  const bool use_rz_paired_pressure =
      geom.IsRZ() &&
      rz_paired_pressure_flux_capable<PROB::ProbRHS>::value;
#else
  const bool use_rz_paired_pressure = false;
#endif
  if (geom.IsRZ() &&
      rz_paired_pressure_flux_required<PROB::ProbRHS>::value &&
      !use_rz_paired_pressure) {
    amrex::Abort(
        "The selected R-Z inviscid scheme requires the paired metric-"
        "advection/ordinary-pressure operator.  The selected diffusive or "
        "geometry path is not compatible; refusing a silent direct-p/r "
        "fallback.");
  }
  if (capture_rz_pressure && !use_rz_paired_pressure) {
    amrex::Abort(
        "R-Z pressure-face capture requires the paired pressure operator");
  }
  static const int s_rz_gp_annular_bic = [] {
    int value = 0;
    ParmParse pp("cns");
    pp.query("rz_gp_annular_bic", value);
    return value;
  }();
  const bool rz_gp_annular_bic = s_rz_gp_annular_bic != 0;
  if (rz_euler_cell_average_deconvolution && !geom.IsRZ()) {
    amrex::Abort(
        "cns.rz_euler_cell_average_deconvolution requires R-Z geometry");
  }
#ifdef CNS_USE_EB
  if (rz_euler_cell_average_deconvolution || rz_gp_annular_bic) {
    amrex::Abort(
        "R-Z annular shared-GP semantics are not an EB flow path");
  }
#endif
#ifdef AMREX_USE_GPIBM
  if (rz_euler_cell_average_deconvolution != rz_gp_annular_bic) {
    amrex::Abort(
        "R-Z GP-IBM requires cns.rz_euler_cell_average_deconvolution and "
        "cns.rz_gp_annular_bic to be enabled together");
  }
  if (use_rz_paired_pressure && CNS::ibm_positivity_flux_limiter &&
      !capture_rz_pressure) {
    amrex::Abort(
        "R-Z paired-pressure GP-IBM positivity limiting requires stage-local "
        "capture of the radial pressure companion");
  }
#else
  if (rz_gp_annular_bic) {
    amrex::Abort("cns.rz_gp_annular_bic requires AMREX_USE_GPIBM");
  }
#endif

#if (AMREX_SPACEDIM < 3)
  // UMZ (theta/z-momentum) is identically zero in 2D but the component exists
  // and cons2prims reads it unconditionally. CPU malloc zero-pages hid this;
  // CUDA arena memory is recycled and NOT zeroed. Sanitize on EVERY build
  // (was IBM-only; GPU-gate fix 2026-07-15).
  statemf.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, statemf.nGrow());
#endif

#ifdef AMREX_USE_GPIBM
  MultiFab& conservative_flux_state = statemf;

  // Convert conserved variables to primitives level-wide, then apply IBM
  // ghost-point corrections in a single pass before the MFIter loop.
  // This avoids redundant per-fab conversions and ensures all ghost-point
  // data are consistent when each fab's flux kernel executes.
  // WENO needs NGHOST cells, while IBM image points can reach farther when
  // alpha/eorder > 1.  The wider primitive-only scratch halo is exchanged in
  // computeAllGPs; the conservative state allocation remains unchanged.
  const int ibm_primitive_nghost =
      IBM::ib.volumeInterpolationNghost(level);
  MultiFab prims_mf(statemf.boxArray(), statemf.DistributionMap(),
                    cls_h.NPRIM, ibm_primitive_nghost,
                    MFInfo().SetArena(The_Async_Arena()));

  {
    BL_PROFILE_VAR("CNS::compute_rhs::cons2prims", prof_cons2prims);
    for (MFIter mfi(conservative_flux_state, false); mfi.isValid(); ++mfi) {
      cls_h.cons2prims(mfi, conservative_flux_state.array(mfi),
                       prims_mf.array(mfi));
    }
    BL_PROFILE_VAR_STOP(prof_cons2prims);
  }
  {
    // Sync wall-motion time before GP reconstruction so that
    // compute_surfIB() sees the correct wall velocity.
#ifdef CNS_USE_FSI
    PROB::Motion::sim_time = cur_time;
#endif

    BL_PROFILE_VAR("IBM::computeAllGPs", prof_gp);
    if (IBM::ib.nsGPCellAverageRecoveryEnabled()) {
      IBM::ib.computeAllGPsCartesianConservativeAverage(
          prims_mf, conservative_flux_state, cls_d, level);
    } else if (IBM::ib.rzAnnularCellAverageEnabled(level)) {
      IBM::ib.computeAllGPsRZAnnular(
          prims_mf, conservative_flux_state, cls_d, level);
    } else {
      IBM::ib.computeAllGPs(prims_mf, cls_d, level);
    }
    BL_PROFILE_VAR_STOP(prof_gp);
  }
#endif

  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    const Array4<Real>& conservative_state = statemf.array(mfi);
#ifdef AMREX_USE_GPIBM
    const Array4<Real>& flux_state = conservative_flux_state.array(mfi);
#else
    const Array4<Real>& flux_state = conservative_state;
#endif

    const Box& cell_box = mfi.growntilebox(0);
    const Box& ghost_cell_box = mfi.growntilebox(cls_h.NGHOST);
#ifdef CNS_USE_EB
    // The experimental EB correction reads one additional face-flux layer.
    const Box& flux_stencil_box =
        mfi.growntilebox(cls_h.NGHOST + 1);
#else
    const Box& flux_stencil_box = mfi.growntilebox(cls_h.NGHOST);
#endif

    // Primitive variables are level-wide for GP-IBM and FAB-local otherwise.
#ifdef AMREX_USE_GPIBM
    // The level-wide array was converted and GP-corrected before this loop.
    Array4<Real> const& prims = prims_mf.array(mfi);
#else
    FArrayBox primitive_scratch(
        ghost_cell_box, cls_h.NPRIM, The_Async_Arena());
    Array4<Real> const& prims = primitive_scratch.array();
#endif

    std::unique_ptr<FArrayBox> rz_center_pressure_owner;
    Array4<const Real> rz_center_pressure;
#ifndef AMREX_USE_GPIBM
    if (rz_euler_cell_average_deconvolution) {
      rz_center_pressure_owner = std::make_unique<FArrayBox>(
          ghost_cell_box, 1, The_Async_Arena());
      const Array4<Real> pressure = rz_center_pressure_owner->array();
      const auto rz_prob_lo = geom.ProbLoArray();
      const auto rz_cell_size = geom.CellSizeArray();
      const int radial_data_lo = ghost_cell_box.smallEnd(0);
      const int radial_data_hi = ghost_cell_box.bigEnd(0);
      const Real radial_origin_over_dr =
          rz_prob_lo[0] / rz_cell_size[0];
      ParallelFor(
          ghost_cell_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pressure(i, j, k) =
                cerisse::rz_fv::annular_average_to_radial_center_pressure(
                    i, j, k, radial_data_lo, radial_data_hi,
                    radial_origin_over_dr, flux_state, cls_d);
          });
      rz_center_pressure = rz_center_pressure_owner->const_array();
    }
#endif
#ifdef CNS_USE_EB
    // Scratch data used only by the isolated experimental EB redistribution.
    FArrayBox eb_divergence_scratch(
        ghost_cell_box, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& eb_divergence = eb_divergence_scratch.array();
    FArrayBox eb_state_scratch(
        ghost_cell_box, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& eb_conservative_state = eb_state_scratch.array();
    ParallelFor(
        ghost_cell_box, cls_h.NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          eb_conservative_state(i, j, k, n) =
              conservative_state(i, j, k, n);
        });
#endif

    // Face-centred flux storage, one FArrayBox for each coordinate direction.
    std::array<FArrayBox, AMREX_SPACEDIM> face_flux_storage;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      face_flux_storage[dir].resize(
          amrex::surroundingNodes(flux_stencil_box, dir), cls_h.NCONS,
          The_Async_Arena());
      face_flux_storage[dir].setVal<RunOn::Device>(Real(0.0));
    }
    const FaceFluxArray face_fluxes{
        AMREX_D_DECL(
            &face_flux_storage[0], &face_flux_storage[1],
            &face_flux_storage[2])};

    std::array<std::unique_ptr<BaseFab<int>>, AMREX_SPACEDIM>
        llf_fallback_mask_storage;
    FaceFallbackMaskArray llf_fallback_mask_views{};
    if (capture_llf_fallback) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        llf_fallback_mask_storage[dir] = std::make_unique<BaseFab<int>>(
            amrex::surroundingNodes(flux_stencil_box, dir), 1,
            The_Async_Arena());
        llf_fallback_mask_storage[dir]->setVal<RunOn::Device>(0);
        llf_fallback_mask_views[dir] =
            llf_fallback_mask_storage[dir]->array();
      }
    }

    // The fixed all-fluid RZ LLF-WENO operator also exposes its paired
    // pressure flux as stage-local scratch for the split divergence below.
    // Its contribution is already folded into the final radial-momentum face
    // flux. That final shared flux is captured and entered in flux registers.
    std::unique_ptr<FArrayBox> rz_pressure_face_owner;
    Array4<Real> rz_pressure_face_flux;
    if (use_rz_paired_pressure) {
      rz_pressure_face_owner = std::make_unique<FArrayBox>(
          amrex::surroundingNodes(flux_stencil_box, 0), 1,
          The_Async_Arena());
      rz_pressure_face_owner->setVal<RunOn::Device>(Real(0.0));
      rz_pressure_face_flux = rz_pressure_face_owner->array();
    }

    // Convert the local conservative state before evaluating the face fluxes.
#ifndef AMREX_USE_GPIBM
    cls_h.cons2prims(mfi, conservative_state, prims);
#endif

    // The companion pressure must use radial-centre values with the same
    // stage-local semantics as the flux stencil.  The pure shared-GP annular
    // path has already published marker-safe fluid/unique-GP centre states to
    // prims; the generic all-fluid annular recovery may use its separate
    // centre-pressure scratch.
    Array4<const Real> pressure_reconstruction_states = prims;
    int pressure_reconstruction_component = cls_h.QPRES;
#ifndef AMREX_USE_GPIBM
    if (rz_euler_cell_average_deconvolution) {
      pressure_reconstruction_states = rz_center_pressure;
      pressure_reconstruction_component = 0;
    }
#endif

    // Geometry markers are present only in the GP-IBM and experimental EB
    // builds. The GP path aliases the existing marker MultiFab.
#if (AMREX_USE_GPIBM && !CNS_USE_EB)
    auto& marker_mf = *IBM::ib.bmf_a[level];
    auto const& geoMarkers = marker_mf.array(mfi);
#elif (CNS_USE_EB && !AMREX_USE_GPIBM)
    // Convert the EB bool marker to the uint8_t marker interface.
    auto& eb_marker_mf = *EBM::eb.bmf_a[level];
    auto const& eb_bool_markers = eb_marker_mf.array(mfi);
    BaseFab<uint8_t> marker_scratch(
        ghost_cell_box, 2, The_Async_Arena());
    auto const& geoMarkers = marker_scratch.array();
    ParallelFor(
        ghost_cell_box, 2,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          geoMarkers(i, j, k, n) =
              static_cast<uint8_t>(eb_bool_markers(i, j, k, n));
        });
#endif

    // Compute the inviscid numerical face fluxes. flux_state contains the
    // conservative variables used by the selected Euler operator.
    {
      BL_PROFILE_VAR("CNS::compute_rhs::eflux", prof_eflux);
#if defined(AMREX_USE_GPIBM)
      compute_shared_gp_euler_fluxes(
          prob_rhs, geom, mfi, prims, face_fluxes, flux_state, cls_d,
          geoMarkers, pressure_reconstruction_states,
          pressure_reconstruction_component, rz_pressure_face_flux,
          capture_llf_fallback, llf_fallback_mask_views);
#elif defined(CNS_USE_EB)
      prob_rhs.eflux_ibm(
          geom, mfi, prims, face_fluxes, flux_state, cls_d, geoMarkers);
#else
      compute_all_fluid_euler_fluxes(
          prob_rhs, geom, mfi, prims, face_fluxes, flux_state, cls_d,
          pressure_reconstruction_states,
          pressure_reconstruction_component,
          rz_pressure_face_flux, capture_llf_fallback,
          llf_fallback_mask_views);
#endif
      BL_PROFILE_VAR_STOP(prof_eflux);
    }
    if (capture_llf_fallback) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const Box face_box = amrex::surroundingNodes(cell_box, dir);
        const auto src = llf_fallback_mask_storage[dir]->const_array();
        const auto dst = captured_llf_fallback_mask[dir]->array(mfi);
        ParallelFor(
            face_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              dst(i, j, k, 0) = src(i, j, k, 0);
            });
      }
    }
    {
      BL_PROFILE_VAR("CNS::compute_rhs::dflux", prof_dflux);
#if defined(AMREX_USE_GPIBM) || defined(CNS_USE_EB)
      prob_rhs.dflux_ibm(
          geom, mfi, prims, face_fluxes, flux_state, cls_d, geoMarkers);
#else
      prob_rhs.dflux(
          geom, mfi, prims, face_fluxes, flux_state, cls_d);
#endif
      BL_PROFILE_VAR_STOP(prof_dflux);
    }
    // Preserve the combined Euler+diffusive face flux before the conservative
    // state storage is reused for the RHS. MFIter is untiled here, so every
    // valid face of every level box is copied exactly once into its owning FAB.
    if (retain_level_face_fluxes) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const Box face_box = amrex::surroundingNodes(cell_box, dir);
        auto const src = face_flux_storage[dir].const_array();
        auto const dst = level_face_fluxes[dir]->array(mfi);
        ParallelFor(face_box, cls_h.NCONS,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                      dst(i, j, k, n) = src(i, j, k, n);
                    });
      }
    }
    if (capture_rz_pressure) {
      const Box radial_face_box =
          amrex::surroundingNodes(cell_box, 0);
      const auto src = rz_pressure_face_owner->const_array();
      const auto dst = level_rz_pressure_face_flux->array(mfi);
      ParallelFor(
          radial_face_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            dst(i, j, k, 0) = src(i, j, k, 0);
          });
    }

    // The input conservative state is no longer needed after both face-flux
    // operators finish. Reuse its storage for the stage RHS.
    const Array4<Real>& rhs = conservative_state;
    const int ncons = cls_h.NCONS;
    ParallelFor(
        ghost_cell_box, ncons,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) = Real(0.0);
        });

    // Convert the finalized face fluxes into the conservative flux-difference
    // contribution to dU/dt. The helper selects Cartesian or axisymmetric RZ
    // metrics from geom without changing the face fluxes.
    assemble_face_flux_divergence<PROB::ProbRHS>(
        geom, cell_box, prims, face_fluxes, rz_center_pressure,
        rz_euler_cell_average_deconvolution, rz_pressure_face_flux,
        use_rz_paired_pressure, rhs, ncons);

    // Add RZ viscous geometric terms that are not part of the face-flux
    // difference. Cartesian operators return immediately.
    if (geom.IsRZ()) {
      prob_rhs.rz_geometric_source(geom, mfi, prims, rhs, cls_d);
    }
#if CNS_USE_EB
    // Experimental EB-only correction. This block is not compiled into the
    // pure shared-GP production executable.
    const Box& eb_box = mfi.growntilebox(0);
    const auto& flag = (*EBM::eb.ebflags_a[level])[mfi];
    const FabType fab_type = flag.getType(eb_box);
    const bool fab_with_eb = (FabType::singlevalued == fab_type);
    if (fab_with_eb) {
      EBM::eb.ebflux(
          geom, mfi, prims, face_fluxes, rhs, cls_d, level);
    }

    ParallelFor(
        ghost_cell_box, ncons,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          eb_divergence(i, j, k, n) = rhs(i, j, k, n);
        });

    if (eb_redistribution && fab_with_eb) {
      EBM::eb.redist(
          geom, mfi, eb_conservative_state, eb_divergence, face_fluxes,
          rhs, cls_d, level, dt, h_phys_bc);
    }
#endif

    // A problem-specific NSCBC hook may replace boundary-cell dU/dt after the
    // conservative face-flux difference has been assembled.
    try_rhs_nscbc(
        0, geom, mfi, prims, rhs, cls_d, pparm_d, dt, cur_time,
        capture_stage_flux);

    // Add problem-defined source terms.
#if defined(AMREX_USE_GPIBM) || defined(CNS_USE_EB)
    prob_rhs.src(
        geom, mfi, prims, rhs, cls_d, dt, cur_time, geoMarkers);
#else
    prob_rhs.src(geom, mfi, prims, rhs, cls_d, dt, cur_time);
#endif

    // Solid cells supply ghost states to crossing stencils but are not
    // advanced. Fluid cells retain the full-Cartesian flux-difference RHS.
#if defined(AMREX_USE_GPIBM) || defined(CNS_USE_EB)
    ParallelFor(
        ghost_cell_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (geoMarkers(i, j, k, 0) == 0) return;
          for (int n = 0; n < ncons; ++n) {
            rhs(i, j, k, n) = Real(0.0);
          }
        });
#endif
  }

  // The production positivity limiter needs the conservative face-flux
  // contribution used by this stage RHS. Copy it only when the caller
  // provides stage-local capture storage.
  if (capture_stage_flux) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      MultiFab::Copy(*captured_face_flux[dir], *level_face_fluxes[dir], 0, 0,
                     cls_h.NCONS, 0);
      captured_face_flux[dir]->OverrideSync(geom.periodicity());
    }
  }
  if (capture_rz_pressure) {
    MultiFab::Copy(*captured_rz_pressure_face_flux,
                   *level_rz_pressure_face_flux, 0, 0, 1, 0);
    captured_rz_pressure_face_flux->OverrideSync(geom.periodicity());
  }
  if (capture_llf_fallback) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      captured_llf_fallback_mask[dir]->OverrideSync(geom.periodicity());
    }
  }
  if (register_fluxes) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (fr_as_crse != nullptr) {
        fr_as_crse->CrseInit(
            *level_face_fluxes[dir], *level_face_areas[dir], dir, 0, 0,
            cls_h.NCONS,
            -reflux_dt, FluxRegister::ADD);
      }
      if (fr_as_fine != nullptr) {
        fr_as_fine->FineAdd(
            *level_face_fluxes[dir], *level_face_areas[dir], dir, 0, 0,
            cls_h.NCONS,
            reflux_dt);
      }
    }
    // The face MultiFabs use The_Async_Arena(), whose deallocation is ordered
    // after work already submitted to the active stream.  An extra stream
    // wait here only serializes the next RK stage.
  }

}
