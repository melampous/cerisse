#ifndef AFD_H_
#define AFD_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <Riemann.h>
#include <IBMSharedGPFluxUtils.h>
#include <Weno.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

// Point-value reconstructions used by the alternative finite-difference (AFD)
// flux.  These weights differ from the finite-volume/flux-splitting weights in
// ReconScheme: the optimal point-interpolation weights are (1, 10, 5).
namespace AfdReconScheme {
namespace detail {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
pow6(const amrex::Real x) noexcept
{
  const amrex::Real x2 = x * x;
  return x2 * x2 * x2;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
smoothness(const amrex::Real s[5], amrex::Real beta[3]) noexcept
{
  using amrex::Real;
  const Real d0 = s[0] - Real(2.0) * s[1] + s[2];
  const Real d1 = s[1] - Real(2.0) * s[2] + s[3];
  const Real d2 = s[2] - Real(2.0) * s[3] + s[4];
  beta[0] = Real(13.0 / 12.0) * d0 * d0 +
            Real(0.25) * (s[0] - Real(4.0) * s[1] + Real(3.0) * s[2]) *
                (s[0] - Real(4.0) * s[1] + Real(3.0) * s[2]);
  beta[1] = Real(13.0 / 12.0) * d1 * d1 +
            Real(0.25) * (s[1] - s[3]) * (s[1] - s[3]);
  beta[2] = Real(13.0 / 12.0) * d2 * d2 +
            Real(0.25) * (Real(3.0) * s[2] - Real(4.0) * s[3] + s[4]) *
                (Real(3.0) * s[2] - Real(4.0) * s[3] + s[4]);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
point_polynomials(const amrex::Real s[5], amrex::Real q[3]) noexcept
{
  using amrex::Real;
  q[0] = (Real(3.0) * s[0] - Real(10.0) * s[1] + Real(15.0) * s[2]) /
         Real(8.0);
  q[1] = (-s[1] + Real(6.0) * s[2] + Real(3.0) * s[3]) / Real(8.0);
  q[2] = (Real(3.0) * s[2] + Real(6.0) * s[3] - s[4]) / Real(8.0);
}

}  // namespace detail

struct WenoZ5 {
  using FluxSplitScheme = ReconScheme::WenoZ5;
  static constexpr const char* scheme_name = "afd_hllc_weno_z5";
  static constexpr int default_shock_hllc_muscl = 1;

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  right(const amrex::Real s[5], const amrex::Real /*cutoff*/) noexcept
  {
    using amrex::Real;
    constexpr Real eps = (std::numeric_limits<Real>::digits >= 53)
                             ? Real(1.0e-40)
                             : std::numeric_limits<Real>::epsilon();
    Real beta[3];
    Real q[3];
    detail::smoothness(s, beta);
    detail::point_polynomials(s, q);

    const Real tau = std::abs(beta[0] - beta[2]);
    Real alpha[3];
    const Real r0 = tau / (eps + beta[0]);
    const Real r1 = tau / (eps + beta[1]);
    const Real r2 = tau / (eps + beta[2]);
    alpha[0] = Real(1.0) * (Real(1.0) + r0 * r0);
    alpha[1] = Real(10.0) * (Real(1.0) + r1 * r1);
    alpha[2] = Real(5.0) * (Real(1.0) + r2 * r2);
    const Real denom = alpha[0] + alpha[1] + alpha[2];
    return (alpha[0] * q[0] + alpha[1] * q[1] + alpha[2] * q[2]) /
           denom;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  right_masked(const amrex::Real s[5], const bool valid[3],
               const amrex::Real /*cutoff*/) noexcept
  {
    using amrex::Real;
    constexpr Real eps = (std::numeric_limits<Real>::digits >= 53)
                             ? Real(1.0e-40)
                             : std::numeric_limits<Real>::epsilon();
    Real beta[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real q[3] = {Real(0.0), Real(0.0), Real(0.0)};
    if (valid[0]) {
      const Real d = s[0] - Real(2.0) * s[1] + s[2];
      beta[0] = Real(13.0 / 12.0) * d * d +
                Real(0.25) *
                    (s[0] - Real(4.0) * s[1] + Real(3.0) * s[2]) *
                    (s[0] - Real(4.0) * s[1] + Real(3.0) * s[2]);
      q[0] = (Real(3.0) * s[0] - Real(10.0) * s[1] +
              Real(15.0) * s[2]) / Real(8.0);
    }
    if (valid[1]) {
      const Real d = s[1] - Real(2.0) * s[2] + s[3];
      beta[1] = Real(13.0 / 12.0) * d * d +
                Real(0.25) * (s[1] - s[3]) * (s[1] - s[3]);
      q[1] = (-s[1] + Real(6.0) * s[2] + Real(3.0) * s[3]) /
             Real(8.0);
    }
    if (valid[2]) {
      const Real d = s[2] - Real(2.0) * s[3] + s[4];
      beta[2] = Real(13.0 / 12.0) * d * d +
                Real(0.25) *
                    (Real(3.0) * s[2] - Real(4.0) * s[3] + s[4]) *
                    (Real(3.0) * s[2] - Real(4.0) * s[3] + s[4]);
      q[2] = (Real(3.0) * s[2] + Real(6.0) * s[3] - s[4]) /
             Real(8.0);
    }

    Real tau = Real(0.0);
    if (valid[0] && valid[2]) tau = std::abs(beta[0] - beta[2]);
    else if (valid[0] && valid[1]) tau = std::abs(beta[0] - beta[1]);
    else if (valid[1] && valid[2]) tau = std::abs(beta[1] - beta[2]);

    Real alpha[3] = {Real(0.0), Real(0.0), Real(0.0)};
    const Real optimal[3] = {Real(1.0), Real(10.0), Real(5.0)};
    for (int n = 0; n < 3; ++n) {
      if (!valid[n]) continue;
      const Real ratio = tau / (eps + beta[n]);
      alpha[n] = optimal[n] * (Real(1.0) + ratio * ratio);
    }
    const Real denom = alpha[0] + alpha[1] + alpha[2];
    if (!(denom > std::numeric_limits<Real>::min())) {
      if (valid[1]) return q[1];
      if (valid[2]) return q[2];
      return q[0];
    }
    return (alpha[0] * q[0] + alpha[1] * q[1] + alpha[2] * q[2]) /
           denom;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  left_masked(const amrex::Real s[5], const bool valid[3],
              const amrex::Real cutoff) noexcept
  {
    const amrex::Real reversed[5] = {s[4], s[3], s[2], s[1], s[0]};
    const bool reversed_valid[3] = {valid[2], valid[1], valid[0]};
    return right_masked(reversed, reversed_valid, cutoff);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  left(const amrex::Real s[5], const amrex::Real cutoff) noexcept
  {
    const amrex::Real reversed[5] = {s[4], s[3], s[2], s[1], s[0]};
    return right(reversed, cutoff);
  }

};

struct Teno5 {
  using FluxSplitScheme = ReconScheme::Teno5;
  static constexpr const char* scheme_name = "afd_hllc_teno5";
  static constexpr int default_shock_hllc_muscl = 1;

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  right(const amrex::Real s[5], const amrex::Real cutoff) noexcept
  {
    using amrex::Real;
    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    Real beta[3];
    Real q[3];
    detail::smoothness(s, beta);
    detail::point_polynomials(s, q);

    const Real tau = std::abs(
        std::abs(beta[0] - beta[2]) -
        (beta[0] + Real(4.0) * beta[1] + beta[2]) / Real(6.0));
    Real gamma[3];
    for (int n = 0; n < 3; ++n) {
      gamma[n] = detail::pow6(Real(1.0) + tau / (eps + beta[n]));
    }
    const Real inv_gamma = Real(1.0) / (gamma[0] + gamma[1] + gamma[2]);
    const Real alpha[3] = {
        gamma[0] * inv_gamma >= cutoff ? Real(1.0) : Real(0.0),
        gamma[1] * inv_gamma >= cutoff ? Real(10.0) : Real(0.0),
        gamma[2] * inv_gamma >= cutoff ? Real(5.0) : Real(0.0),
    };
    const Real denom = alpha[0] + alpha[1] + alpha[2];
    return (alpha[0] * q[0] + alpha[1] * q[1] + alpha[2] * q[2]) /
           denom;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  right_masked(const amrex::Real s[5], const bool valid[3],
               const amrex::Real cutoff) noexcept
  {
    using amrex::Real;
    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    Real beta[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real q[3] = {Real(0.0), Real(0.0), Real(0.0)};
    if (valid[0]) {
      const Real d = s[0] - Real(2.0) * s[1] + s[2];
      beta[0] = Real(13.0 / 12.0) * d * d +
                Real(0.25) *
                    (s[0] - Real(4.0) * s[1] + Real(3.0) * s[2]) *
                    (s[0] - Real(4.0) * s[1] + Real(3.0) * s[2]);
      q[0] = (Real(3.0) * s[0] - Real(10.0) * s[1] +
              Real(15.0) * s[2]) / Real(8.0);
    }
    if (valid[1]) {
      const Real d = s[1] - Real(2.0) * s[2] + s[3];
      beta[1] = Real(13.0 / 12.0) * d * d +
                Real(0.25) * (s[1] - s[3]) * (s[1] - s[3]);
      q[1] = (-s[1] + Real(6.0) * s[2] + Real(3.0) * s[3]) /
             Real(8.0);
    }
    if (valid[2]) {
      const Real d = s[2] - Real(2.0) * s[3] + s[4];
      beta[2] = Real(13.0 / 12.0) * d * d +
                Real(0.25) *
                    (Real(3.0) * s[2] - Real(4.0) * s[3] + s[4]) *
                    (Real(3.0) * s[2] - Real(4.0) * s[3] + s[4]);
      q[2] = (Real(3.0) * s[2] + Real(6.0) * s[3] - s[4]) /
             Real(8.0);
    }

    Real tau = Real(0.0);
    if (valid[0] && valid[1] && valid[2]) {
      tau = std::abs(std::abs(beta[0] - beta[2]) -
                     (beta[0] + Real(4.0) * beta[1] + beta[2]) /
                         Real(6.0));
    } else if (valid[0] && valid[2]) {
      tau = std::abs(beta[0] - beta[2]);
    } else if (valid[0] && valid[1]) {
      tau = std::abs(beta[0] - beta[1]);
    } else if (valid[1] && valid[2]) {
      tau = std::abs(beta[1] - beta[2]);
    }

    Real gamma[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real gamma_sum = Real(0.0);
    for (int n = 0; n < 3; ++n) {
      if (!valid[n]) continue;
      gamma[n] = detail::pow6(Real(1.0) + tau / (eps + beta[n]));
      gamma_sum += gamma[n];
    }
    if (!(gamma_sum > std::numeric_limits<Real>::min())) {
      return WenoZ5::right_masked(s, valid, cutoff);
    }
    const Real inv_gamma = Real(1.0) / gamma_sum;
    const Real optimal[3] = {Real(1.0), Real(10.0), Real(5.0)};
    Real alpha[3] = {Real(0.0), Real(0.0), Real(0.0)};
    for (int n = 0; n < 3; ++n) {
      if (valid[n] && gamma[n] * inv_gamma >= cutoff) {
        alpha[n] = optimal[n];
      }
    }
    const Real denom = alpha[0] + alpha[1] + alpha[2];
    if (!(denom > std::numeric_limits<Real>::min())) {
      return WenoZ5::right_masked(s, valid, cutoff);
    }
    return (alpha[0] * q[0] + alpha[1] * q[1] + alpha[2] * q[2]) /
           denom;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  left_masked(const amrex::Real s[5], const bool valid[3],
              const amrex::Real cutoff) noexcept
  {
    const amrex::Real reversed[5] = {s[4], s[3], s[2], s[1], s[0]};
    const bool reversed_valid[3] = {valid[2], valid[1], valid[0]};
    return right_masked(reversed, reversed_valid, cutoff);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  left(const amrex::Real s[5], const amrex::Real cutoff) noexcept
  {
    const amrex::Real reversed[5] = {s[4], s[3], s[2], s[1], s[0]};
    return right(reversed, cutoff);
  }

};

}  // namespace AfdReconScheme

// Fifth-order point-value WENO/TENO states + HLLC interface flux + AFD flux
// corrections for the five-equation ideal-gas system.  In R-Z coordinates,
// the radial correction is applied to the metric flux rF while the stored
// face flux retains its per-unit-area F semantics.  In GPIBM builds the full
// six-point AFD formula is retained on all-fluid stencils, while immersed-
// boundary crossings use the standard marker-aware shared-GP reconstruction.
// The reconstruction requires three primitive ghost cells.
template <typename PointScheme, typename cls_t>
class afd_hllc_t : public riemann_t<false, cls_t> {
 public:
  static constexpr int ng = 3;
  static constexpr const char* scheme_name = PointScheme::scheme_name;
  // The radial-axis face stores the AFD auxiliary metric flux \hat{rF}.
  // Non-axis face entries retain the usual per-unit-area flux convention.
  static constexpr bool rz_radial_axis_face_flux_is_metric = true;
  static constexpr bool rz_paired_pressure_flux_capable = true;
  static constexpr bool rz_paired_pressure_flux_required = true;
  static constexpr bool rz_annular_pressure_consistency_capable = false;
  static constexpr bool rz_euler_geometric_source_active = true;

  static_assert(NUM_SPECIES == 1,
                "AFD-HLLC currently supports one species only");
  static_assert(cls_t::NCONS == 5,
                "AFD-HLLC expects the five-equation ideal-gas system");
  static_assert(cls_t::NGHOST >= ng,
                "AFD-HLLC requires at least three primitive ghost cells");

  static void print_local_llf_fallback_manifest()
  {
    const RuntimeOptions opts = runtime_options();
    amrex::Print()
        << "[Numerics] local_first_order_llf_fallback="
        << (opts.local_llf_fallback != 0 ? "enabled" : "disabled")
        << " density_jump=" << opts.fallback_density_jump
        << " density_curvature=" << opts.fallback_density_curvature
        << " scope=cartesian_rz_all_fluid_shared_gp"
        << " first_interior_rz_gate=density_jump_only"
        << " axis_face=excluded\n";
  }

  void eflux(const amrex::Geometry& geom, const amrex::MFIter& mfi,
             const amrex::Array4<const amrex::Real>& prims,
             std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
             const amrex::Array4<amrex::Real>& /*rhs*/, const cls_t* cls)
  {
    amrex::Array4<amrex::Real> unused_pressure_face_flux;
    compute_face_fluxes<false>(
        geom, mfi, prims, flxt, cls, prims, cls_t::QPRES,
        unused_pressure_face_flux);
  }

  void eflux_with_rz_paired_pressure(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims,
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& /*rhs*/, const cls_t* cls,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux)
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.IsRZ(),
        "AFD-HLLC paired pressure flux requires cylindrical R-Z geometry");
    compute_face_fluxes<true>(
        geom, mfi, prims, flxt, cls, pressure_states,
        pressure_component, radial_pressure_face_flux);
  }

  template <bool ComputeRzPairedPressure>
  void compute_face_fluxes(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims,
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const cls_t* cls,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux)
  {
    const amrex::Box bx = mfi.tilebox();

    const RuntimeOptions opts = runtime_options();
    const int correction = opts.correction;
    const int shock_llf = opts.shock_llf;
    const int shock_hllc_muscl = opts.shock_hllc_muscl;
    const amrex::Real hllc_reconstructed_speed_ratio =
        opts.hllc_reconstructed_speed_ratio;
    const amrex::Real smoothness_threshold = opts.smoothness_threshold;
    const amrex::Real shock_pressure_jump = opts.shock_pressure_jump;
    const amrex::Real shock_compression = opts.shock_compression;
    const amrex::Real teno_cutoff = opts.teno_cutoff;
    const int local_llf_fallback = opts.local_llf_fallback;
    const amrex::Real fallback_density_jump =
        opts.fallback_density_jump;
    const amrex::Real fallback_density_curvature =
        opts.fallback_density_curvature;
#ifdef CERISSE_ENABLE_FACE_TRACE
    const int face_trace = opts.face_trace;
    const int face_trace_i = opts.face_trace_i;
    const int face_trace_j = opts.face_trace_j;
    const int face_trace_k = opts.face_trace_k;
    const int face_trace_dir = opts.face_trace_dir;
#endif
    const auto cell_size = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const int radial_axis_face_index = geom.Domain().smallEnd(0);
    const int radial_domain_high_index = geom.Domain().bigEnd(0);
    const bool is_rz = geom.IsRZ();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const bool rz_radial_direction = is_rz && dir == 0;
      const amrex::Box flxbx = amrex::surroundingNodes(bx, dir);
      const auto flx = flxt[dir]->array();
      amrex::ParallelFor<128>(
          flxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // NVCC 12 rejects a variable first captured inside if constexpr.
            amrex::ignore_unused(
                radial_domain_high_index, radial_pressure_face_flux);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const bool trace_this_face =
#ifdef CERISSE_ENABLE_FACE_TRACE
                face_trace != 0 && dir == face_trace_dir &&
                i == face_trace_i && j == face_trace_j &&
                k == face_trace_k;
#else
                false;
#endif
            amrex::Real pressure_face = this->face_flux_bulk(
                iv, dir, prims, flx, *cls, correction, shock_llf,
                shock_hllc_muscl,
                hllc_reconstructed_speed_ratio,
                smoothness_threshold, shock_pressure_jump,
                shock_compression, teno_cutoff, cell_size,
                local_llf_fallback, fallback_density_jump,
                fallback_density_curvature,
                rz_radial_direction, prob_lo[0], cell_size[0],
                radial_axis_face_index,
                pressure_states, pressure_component, trace_this_face);
            if constexpr (ComputeRzPairedPressure) {
              if (rz_radial_direction) {
                if (iv[0] == radial_axis_face_index &&
                    prob_lo[0] == amrex::Real(0.0)) {
                  const amrex::IntVect radial =
                      amrex::IntVect::TheDimensionVector(0);
                  const bool use_high_order_axis_pressure =
                      correction != 0 &&
                      this->smooth_stencil(
                          iv, radial, prims, smoothness_threshold);
                  pressure_face = this->axis_pressure_auxiliary_flux(
                      iv, pressure_states, pressure_component,
                      radial_domain_high_index,
                      use_high_order_axis_pressure);
                  // Preserve the AFD auxiliary metric h_(rF) for mass,
                  // axial momentum, and energy.  Those radial metric fluxes
                  // are even at the axis and their discrete h-function need
                  // not vanish.  Only r*(rho*u_r^2+p) is odd.
                  flx(iv, cls_t::UMX) = amrex::Real(0.0);
                }
                radial_pressure_face_flux(iv, 0) = pressure_face;
                if (!std::isfinite(pressure_face)) {
                  const amrex::Real invalid =
                      std::numeric_limits<amrex::Real>::quiet_NaN();
                  for (int n = 0; n < cls_t::NCONS; ++n) {
                    flx(iv, n) = invalid;
                  }
                }
              }
            }
          });
    }
  }

  // IBM-specific AFD-HLLC entry points.
  #include "AfdIBMDrivers.h"

 private:
  struct RuntimeOptions {
    int correction;
    int shock_llf;
    int shock_hllc_muscl;
    amrex::Real hllc_reconstructed_speed_ratio;
    amrex::Real smoothness_threshold;
    amrex::Real shock_pressure_jump;
    amrex::Real shock_compression;
    amrex::Real teno_cutoff;
    amrex::Real ibm_llf_threshold;
    int local_llf_fallback;
    amrex::Real fallback_density_jump;
    amrex::Real fallback_density_curvature;
    int face_trace;
    int face_trace_i;
    int face_trace_j;
    int face_trace_k;
    int face_trace_dir;
  };

  static RuntimeOptions runtime_options()
  {
    static const RuntimeOptions opts = [] {
      RuntimeOptions value{1, 1, PointScheme::default_shock_hllc_muscl,
                           amrex::Real(0.0), amrex::Real(0.08),
                           amrex::Real(1.0e-2), amrex::Real(1.0e-3),
                           amrex::Real(1.0e-4), amrex::Real(0.05),
                           1, amrex::Real(3.0e-3), amrex::Real(8.0e-2),
                           0, 0, 0, 0, 0};
      amrex::ParmParse pp("cns");
      pp.query("afd_correction", value.correction);
      pp.query("afd_shock_llf", value.shock_llf);
      pp.query("afd_shock_hllc_muscl", value.shock_hllc_muscl);
      pp.query("afd_hllc_reconstructed_speed_ratio",
               value.hllc_reconstructed_speed_ratio);
      pp.query("afd_smoothness_threshold", value.smoothness_threshold);
      pp.query("afd_shock_pressure_jump", value.shock_pressure_jump);
      pp.query("afd_shock_compression", value.shock_compression);
      pp.query("afd_teno_cutoff", value.teno_cutoff);
      pp.query("afd_ibm_llf_threshold", value.ibm_llf_threshold);
      pp.query("llf_first_order_fallback", value.local_llf_fallback);
      pp.query("llf_fallback_density_jump", value.fallback_density_jump);
      pp.query(
          "llf_fallback_density_curvature",
          value.fallback_density_curvature);
#ifdef CERISSE_ENABLE_FACE_TRACE
      pp.query("afd_face_trace", value.face_trace);
      pp.query("afd_face_trace_i", value.face_trace_i);
      pp.query("afd_face_trace_j", value.face_trace_j);
      pp.query("afd_face_trace_k", value.face_trace_k);
      pp.query("afd_face_trace_dir", value.face_trace_dir);
#endif
      if (value.correction != 0 && value.correction != 1) {
        amrex::Abort("cns.afd_correction must be 0 or 1");
      }
      if (value.shock_llf != 0 && value.shock_llf != 1) {
        amrex::Abort("cns.afd_shock_llf must be 0 or 1");
      }
      if (value.shock_hllc_muscl != 0 &&
          value.shock_hllc_muscl != 1) {
        amrex::Abort("cns.afd_shock_hllc_muscl must be 0 or 1");
      }
      if (!std::isfinite(value.hllc_reconstructed_speed_ratio) ||
          value.hllc_reconstructed_speed_ratio < amrex::Real(0.0) ||
          (value.hllc_reconstructed_speed_ratio > amrex::Real(0.0) &&
           value.hllc_reconstructed_speed_ratio < amrex::Real(1.0))) {
        amrex::Abort(
            "cns.afd_hllc_reconstructed_speed_ratio must be 0 or >= 1");
      }
      if (!std::isfinite(value.smoothness_threshold) ||
          value.smoothness_threshold < amrex::Real(0.0) ||
          value.smoothness_threshold > amrex::Real(1.0)) {
        amrex::Abort("cns.afd_smoothness_threshold must be in [0,1]");
      }
      if (!std::isfinite(value.shock_pressure_jump) ||
          value.shock_pressure_jump < amrex::Real(0.0) ||
          value.shock_pressure_jump > amrex::Real(1.0)) {
        amrex::Abort("cns.afd_shock_pressure_jump must be in [0,1]");
      }
      if (!std::isfinite(value.shock_compression) ||
          value.shock_compression < amrex::Real(0.0) ||
          value.shock_compression > amrex::Real(1.0)) {
        amrex::Abort("cns.afd_shock_compression must be in [0,1]");
      }
      if (!std::isfinite(value.teno_cutoff) ||
          value.teno_cutoff < amrex::Real(0.0) ||
          value.teno_cutoff >= amrex::Real(1.0 / 3.0)) {
        amrex::Abort("cns.afd_teno_cutoff must be in [0,1/3)");
      }
      if (!std::isfinite(value.ibm_llf_threshold) ||
          value.ibm_llf_threshold < amrex::Real(0.0) ||
          value.ibm_llf_threshold > amrex::Real(1.0)) {
        amrex::Abort("cns.afd_ibm_llf_threshold must be in [0,1]");
      }
      if (value.local_llf_fallback != 0 &&
          value.local_llf_fallback != 1) {
        amrex::Abort("cns.llf_first_order_fallback must be 0 or 1");
      }
      if (!std::isfinite(value.fallback_density_jump) ||
          value.fallback_density_jump < amrex::Real(0.0)) {
        amrex::Abort(
            "cns.llf_fallback_density_jump must be finite and non-negative");
      }
      if (!std::isfinite(value.fallback_density_curvature) ||
          value.fallback_density_curvature < amrex::Real(0.0) ||
          value.fallback_density_curvature > amrex::Real(1.0)) {
        amrex::Abort(
            "cns.llf_fallback_density_curvature must be in [0,1]");
      }
#ifdef CERISSE_ENABLE_FACE_TRACE
      if (value.face_trace != 0 && value.face_trace != 1) {
        amrex::Abort("cns.afd_face_trace must be 0 or 1");
      }
      if (value.face_trace_dir < 0 ||
          value.face_trace_dir >= AMREX_SPACEDIM) {
        amrex::Abort("cns.afd_face_trace_dir is outside AMREX_SPACEDIM");
      }
#endif
      return value;
    }();
    return opts;
  }

  // IBM-specific reconstruction and checked Riemann kernels.
  #include "AfdIBM.h"

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  shock_stencil(
      const amrex::IntVect& iv, const amrex::IntVect& ivd,
      const amrex::Array4<const amrex::Real>& q,
      const amrex::Real pressure_jump_threshold,
      const amrex::Real compression_threshold,
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cell_size) noexcept
  {
    using amrex::Real;
    const Real tiny = std::numeric_limits<Real>::min();

    Real maximum_pressure_jump = Real(0.0);
    for (int m = 0; m < 5; ++m) {
      const amrex::IntVect left = iv + (m - 3) * ivd;
      const amrex::IntVect right = left + ivd;
      const Real pressure_left = q(left, cls_t::QPRES);
      const Real pressure_right = q(right, cls_t::QPRES);
      const Real pressure_scale =
          amrex::max(pressure_left, pressure_right);
      const Real pressure_jump =
          std::abs(pressure_right - pressure_left) /
          amrex::max(pressure_scale, tiny);
      maximum_pressure_jump =
          amrex::max(maximum_pressure_jump, pressure_jump);
    }
    if (!(maximum_pressure_jump > pressure_jump_threshold)) {
      return false;
    }

    Real minimum_cell_size = cell_size[0];
    for (int d = 1; d < AMREX_SPACEDIM; ++d) {
      minimum_cell_size = amrex::min(minimum_cell_size, cell_size[d]);
    }

    // The compression gate uses negative velocity divergence.  This keeps
    // transverse shock coupling for a grid-aligned shock while excluding
    // density jumps and tangential-velocity curvature by themselves.
    Real maximum_compression = Real(0.0);
    for (int m = 1; m < 5; ++m) {
      const amrex::IntVect center = iv + (m - 3) * ivd;
      Real divergence = Real(0.0);
      Real sound_scale = q(center, cls_t::QC);
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const amrex::IntVect direction =
            amrex::IntVect::TheDimensionVector(d);
        const amrex::IntVect minus = center - direction;
        const amrex::IntVect plus = center + direction;
        divergence +=
            (q(plus, cls_t::QU + d) - q(minus, cls_t::QU + d)) /
            (Real(2.0) * cell_size[d]);
        sound_scale = amrex::max(
            sound_scale,
            amrex::max(q(minus, cls_t::QC), q(plus, cls_t::QC)));
      }
      const Real compression =
          amrex::max(-divergence, Real(0.0)) * minimum_cell_size /
          amrex::max(sound_scale, tiny);
      maximum_compression =
          amrex::max(maximum_compression, compression);
    }
    return maximum_compression > compression_threshold;
  }

  /**
   * Return true when a face must use the shared two-state LLF fallback.
   *
   * The detector is identical in every coordinate direction and for both
   * point-reconstruction schemes.  The availability mask prevents a
   * shared-GP stencil from reading an unreconstructed deep-solid state.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  local_first_order_llf_needed(
      const amrex::IntVect& iv, const amrex::IntVect& ivd,
      const amrex::Array4<const amrex::Real>& q,
      const bool state_available[2 * ng], const int enabled,
      const amrex::Real density_jump_threshold,
      const amrex::Real density_curvature_threshold,
      const bool rz_radial_direction,
      const amrex::Real radial_origin,
      const amrex::Real radial_spacing,
      const int radial_axis_face_index) noexcept
  {
    using amrex::Real;
    if (enabled == 0 || !state_available[ng - 1] ||
        !state_available[ng]) {
      return false;
    }

    const Real face_radius = rz_radial_direction
                                 ? radial_origin +
                                       Real(iv[0] - radial_axis_face_index) *
                                           radial_spacing
                                 : Real(1.0);
    if (rz_radial_direction && !(face_radius > Real(0.0))) {
      return false;
    }

    const Real density_left = q(iv - ivd, cls_t::QRHO);
    const Real density_right = q(iv, cls_t::QRHO);
    if (!(density_left > Real(0.0)) || !(density_right > Real(0.0)) ||
        !std::isfinite(density_left + density_right)) {
      return true;
    }
    const Real density_floor = amrex::max(
        amrex::min(density_left, density_right),
        std::numeric_limits<Real>::min());
    const Real relative_density_jump =
        std::abs(density_right - density_left) / density_floor;
    if (!(relative_density_jump > density_jump_threshold)) {
      return false;
    }

    Real density_scale = Real(0.0);
    for (int sample = 0; sample < 2 * ng; ++sample) {
      if (!state_available[sample]) {
        continue;
      }
      density_scale = amrex::max(
          density_scale,
          std::abs(q(iv + (sample - ng) * ivd, cls_t::QRHO)));
    }

    bool density_is_nonsmooth = false;
    for (int group = ng - 2; group <= ng - 1; ++group) {
      if (!state_available[group] || !state_available[group + 1] ||
          !state_available[group + 2]) {
        continue;
      }
      const amrex::IntVect first = iv + (group - ng) * ivd;
      const Real first_density = q(first, cls_t::QRHO);
      const Real second_density = q(first + ivd, cls_t::QRHO);
      const Real third_density = q(first + 2 * ivd, cls_t::QRHO);
      const Real denominator =
          std::abs(first_density) + Real(2.0) * std::abs(second_density) +
          std::abs(third_density) + density_scale +
          std::numeric_limits<Real>::min();
      const Real normalized_curvature =
          std::abs(first_density - Real(2.0) * second_density +
                   third_density) /
          denominator;
      density_is_nonsmooth = density_is_nonsmooth ||
          normalized_curvature > density_curvature_threshold;
    }

    const Real radius_tolerance =
        Real(8.0) * std::numeric_limits<Real>::epsilon() *
        amrex::max(Real(1.0), std::abs(radial_spacing));
    const bool first_interior_rz_radial_face =
        rz_radial_direction && radial_origin == Real(0.0) &&
        std::abs(face_radius - radial_spacing) <= radius_tolerance;
    return density_is_nonsmooth || first_interior_rz_radial_face;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  axis_pressure_auxiliary_flux(
      const amrex::IntVect& axis_face,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const int radial_domain_high_index,
      const bool use_high_order) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect radial =
        amrex::IntVect::TheDimensionVector(0);
    const Real pressure_ring1 =
        pressure_states(axis_face, pressure_component);
    if (!use_high_order ||
        axis_face[0] + 1 > radial_domain_high_index) {
      return pressure_ring1;
    }
    const Real pressure_ring2 =
        pressure_states(axis_face + radial, pressure_component);
    if (axis_face[0] + 2 > radial_domain_high_index) {
      return (Real(7.0) * pressure_ring1 - pressure_ring2) /
             Real(6.0);
    }
    const Real pressure_ring3 =
        pressure_states(axis_face + 2 * radial, pressure_component);
    // Parity-optimal H6 boundary h-function. Together with the AFD
    // nonmetric C[p] correction on the first positive-radius face, this is
    // exact for the even pressure moments 1, r^2, r^4 and r^6.
    return (Real(37.0) * pressure_ring1 -
            Real(8.0) * pressure_ring2 + pressure_ring3) /
           Real(30.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  adjacent_pressure_average(
      const amrex::IntVect& face, const amrex::IntVect& direction,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component) const noexcept
  {
    return amrex::Real(0.5) *
           (pressure_states(face - direction, pressure_component) +
            pressure_states(face, pressure_component));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  face_flux(const amrex::IntVect& iv, const int dir,
            const amrex::Array4<const amrex::Real>& q,
            const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
            const int correction, const int shock_llf,
            const int shock_hllc_muscl,
            const amrex::Real hllc_reconstructed_speed_ratio,
            const amrex::Real smoothness_threshold,
            const amrex::Real teno_cutoff,
            const int local_llf_fallback,
            const amrex::Real fallback_density_jump,
            const amrex::Real fallback_density_curvature,
            const bool rz_radial_direction,
            const amrex::Real radial_origin,
            const amrex::Real radial_spacing,
            const int radial_axis_face_index,
            const amrex::Array4<const amrex::Real>& pressure_states,
            const int pressure_component) const noexcept
  {
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> unused_cell_size{};
    return face_flux_impl(
        iv, dir, q, flx, cls, correction, shock_llf,
        shock_hllc_muscl, hllc_reconstructed_speed_ratio,
        smoothness_threshold, amrex::Real(0.0), amrex::Real(0.0),
        teno_cutoff, unused_cell_size, false, local_llf_fallback,
        fallback_density_jump, fallback_density_curvature,
        rz_radial_direction, radial_origin, radial_spacing,
        radial_axis_face_index, pressure_states, pressure_component,
        false);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  face_flux_bulk(
      const amrex::IntVect& iv, const int dir,
      const amrex::Array4<const amrex::Real>& q,
      const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
      const int correction, const int shock_llf,
      const int shock_hllc_muscl,
      const amrex::Real hllc_reconstructed_speed_ratio,
      const amrex::Real smoothness_threshold,
      const amrex::Real shock_pressure_jump,
      const amrex::Real shock_compression,
      const amrex::Real teno_cutoff,
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cell_size,
      const int local_llf_fallback,
      const amrex::Real fallback_density_jump,
      const amrex::Real fallback_density_curvature,
      const bool rz_radial_direction,
      const amrex::Real radial_origin,
      const amrex::Real radial_spacing,
      const int radial_axis_face_index,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const bool trace_this_face) const noexcept
  {
    return face_flux_impl(
        iv, dir, q, flx, cls, correction, shock_llf,
        shock_hllc_muscl, hllc_reconstructed_speed_ratio,
        smoothness_threshold, shock_pressure_jump, shock_compression,
        teno_cutoff, cell_size, true, local_llf_fallback,
        fallback_density_jump, fallback_density_curvature,
        rz_radial_direction, radial_origin, radial_spacing,
        radial_axis_face_index, pressure_states, pressure_component,
        trace_this_face);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  face_flux_impl(
      const amrex::IntVect& iv, const int dir,
      const amrex::Array4<const amrex::Real>& q,
      const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
      const int correction, const int shock_llf,
      const int shock_hllc_muscl,
      const amrex::Real hllc_reconstructed_speed_ratio,
      const amrex::Real smoothness_threshold,
      const amrex::Real shock_pressure_jump,
      const amrex::Real shock_compression,
      const amrex::Real teno_cutoff,
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cell_size,
      const bool use_pressure_compression_gate,
      const int local_llf_fallback,
      const amrex::Real fallback_density_jump,
      const amrex::Real fallback_density_curvature,
      const bool rz_radial_direction,
      const amrex::Real radial_origin,
      const amrex::Real radial_spacing,
      const int radial_axis_face_index,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const bool trace_this_face) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    const Real low_order_pressure =
        adjacent_pressure_average(
            iv, ivd, pressure_states, pressure_component);
    const int trace_k = AMREX_D_PICK(0, 0, iv[2]);
    const Real face_radius = rz_radial_direction
                                 ? radial_origin +
                                       Real(iv[0] - radial_axis_face_index) *
                                           radial_spacing
                                 : Real(1.0);
    const bool radial_axis_face =
        rz_radial_direction && !(face_radius > Real(0.0));

    if (!valid_stencil(iv, ivd, q)) {
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] face=(%d,%d,%d) dir=%d "
               "branch=invalid_stencil_cell_llf\n",
               iv[0], iv[1], trace_k, dir);
      }
      cell_centered_llf(iv, dir, q, flx, cls);
      if (radial_axis_face) {
        for (int n = 0; n < cls_t::NCONS; ++n) {
          flx(iv, n) = Real(0.0);
        }
      }
      return low_order_pressure;
    }

    const bool all_states_available[2 * ng] = {
        true, true, true, true, true, true};
    if (local_first_order_llf_needed(
            iv, ivd, q, all_states_available, local_llf_fallback,
            fallback_density_jump, fallback_density_curvature,
            rz_radial_direction, radial_origin, radial_spacing,
            radial_axis_face_index)) {
      cell_centered_llf(iv, dir, q, flx, cls);
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] branch=local_first_order_llf\n");
      }
      return low_order_pressure;
    }
    const bool smooth = smooth_stencil(iv, ivd, q, smoothness_threshold);
    bool shock_candidate = false;
    if (use_pressure_compression_gate &&
        ((((shock_llf != 0) || (shock_hllc_muscl != 0)) && !smooth) ||
         trace_this_face)) {
      shock_candidate = shock_stencil(iv, ivd, q, shock_pressure_jump,
                                      shock_compression, cell_size);
    }
    const bool troubled_shock_face =
        !smooth && (!use_pressure_compression_gate || shock_candidate);
    const bool use_hllc_muscl =
        shock_hllc_muscl != 0 && troubled_shock_face;
    const bool use_high_order_llf = shock_llf != 0 && troubled_shock_face;
    if (trace_this_face) {
      printf("[AFD-FACE-TRACE] face=(%d,%d,%d) dir=%d smooth=%d "
             "shock_candidate=%d shock_hllc_muscl=%d shock_llf=%d "
             "selected_hllc_muscl=%d selected_high_order_llf=%d "
             "llf_emergency_available=%d\n",
             iv[0], iv[1], trace_k, dir, int(smooth), int(shock_candidate),
             shock_hllc_muscl, shock_llf, int(use_hllc_muscl),
             int(use_high_order_llf && !use_hllc_muscl),
             int(use_high_order_llf));
    }
    if (use_hllc_muscl) {
      Real muscl_pressure = low_order_pressure;
      if (muscl_hllc_flux(
              iv, dir, q, flx, cls, trace_this_face,
              &muscl_pressure)) {
        if (radial_axis_face) {
          for (int n = 0; n < cls_t::NCONS; ++n) {
            flx(iv, n) = Real(0.0);
          }
        }
        return muscl_pressure;
      }
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] branch=invalid_hllc_muscl_candidate "
               "emergency=%s\n",
               use_high_order_llf ? "characteristic_llf_high_order"
                                  : "cell_centered_llf");
      }
      if (!use_high_order_llf) {
        cell_centered_llf(iv, dir, q, flx, cls);
        if (radial_axis_face) {
          for (int n = 0; n < cls_t::NCONS; ++n) {
            flx(iv, n) = Real(0.0);
          }
        }
        return low_order_pressure;
      }
    }
    if (use_high_order_llf) {
      characteristic_llf_high_order(iv, dir, q, flx, cls);
      if (radial_axis_face) {
        for (int n = 0; n < cls_t::NCONS; ++n) {
          flx(iv, n) = Real(0.0);
        }
      }
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] branch=characteristic_llf_high_order "
               "flux=(rho=%.17g,mn=%.17g,mt1=%.17g,mt2=%.17g,E=%.17g)\n",
               double(flx(iv, cls_t::URHO)),
               double(flx(iv, cls_t::UMX + dir)),
               double(flx(iv, cls_t::UMX + (dir == 0 ? 1 : 0))),
               double(flx(iv, cls_t::UMX + (dir == 2 ? 1 : 2))),
               double(flx(iv, cls_t::UET)));
      }
      return low_order_pressure;
    }

    const int t1 = dir == 0 ? 1 : 0;
    const int t2 = dir == 2 ? 1 : 2;
    Real wc[6][6];
    for (int m = 0; m < 6; ++m) {
      primitive_to_characteristic(iv + (m - 3) * ivd, dir, t1, t2, q,
                                  wc[m]);
    }

    Real wl[6];
    Real wr[6];
    for (int n = 0; n < 6; ++n) {
      Real sl[5];
      Real sr[5];
      for (int m = 0; m < 5; ++m) {
        sl[m] = wc[m][n];
        sr[m] = wc[m + 1][n];
      }
      wl[n] = PointScheme::right(sl, teno_cutoff);
      wr[n] = PointScheme::left(sr, teno_cutoff);
    }

    Real rl, ul, ut1l, ut2l, pl, gl;
    Real rr, ur, ut1r, ut2r, pr, gr;
    if (!characteristic_to_primitive(wl, rl, ul, ut1l, ut2l, pl, gl) ||
        !characteristic_to_primitive(wr, rr, ur, ut1r, ut2r, pr, gr)) {
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] branch=invalid_reconstructed_state_cell_llf\n");
      }
      cell_centered_llf(iv, dir, q, flx, cls);
      if (radial_axis_face) {
        for (int n = 0; n < cls_t::NCONS; ++n) {
          flx(iv, n) = Real(0.0);
        }
      }
      return low_order_pressure;
    }

    const Real yl[1] = {Real(1.0)};
    const Real yr[1] = {Real(1.0)};
    Real el;
    Real er;
    cls.RYP2E(rl, yl, pl, el);
    cls.RYP2E(rr, yr, pr, er);
    el += Real(0.5) * (ul * ul + ut1l * ut1l + ut2l * ut2l);
    er += Real(0.5) * (ur * ur + ut1r * ut1r + ut2r * ut2r);
    const Real cl = std::sqrt(gl * pl / rl);
    const Real cr = std::sqrt(gr * pr / rr);

    // A posteriori guard for a finite but nonphysical point-WENO/HLLC
    // candidate.  Each AFD acoustic variable currently has its own nonlinear
    // WENO weights, so cancellation during the inverse transform can leave a
    // small positive pressure and an enormous velocity.  Positivity checks on
    // rho and p alone do not reject that state.  Compare the reconstructed
    // total signal speed with the physical six-cell envelope and replace only
    // that shared face by the existing conservative HLLC-MUSCL candidate.
    // A value of zero disables the guard; smooth solutions have a ratio near
    // one and therefore retain the original fifth-order AFD flux bit-for-bit.
    if (hllc_reconstructed_speed_ratio > Real(0.0)) {
      Real stencil_signal_speed = std::numeric_limits<Real>::epsilon();
      for (int m = 0; m < 6; ++m) {
        const amrex::IntVect p = iv + (m - 3) * ivd;
        Real velocity_squared = Real(0.0);
        for (int d = 0; d < 3; ++d) {
          const Real velocity = q(p, cls_t::QU + d);
          velocity_squared += velocity * velocity;
        }
        stencil_signal_speed = amrex::max(
            stencil_signal_speed,
            std::sqrt(velocity_squared) + q(p, cls_t::QC));
      }
      const Real reconstructed_left_speed =
          std::sqrt(ul * ul + ut1l * ut1l + ut2l * ut2l) + cl;
      const Real reconstructed_right_speed =
          std::sqrt(ur * ur + ut1r * ut1r + ut2r * ut2r) + cr;
      const Real reconstructed_signal_speed =
          amrex::max(reconstructed_left_speed, reconstructed_right_speed);
      const Real reconstructed_speed_ratio =
          reconstructed_signal_speed / stencil_signal_speed;
      if (!std::isfinite(reconstructed_speed_ratio) ||
          reconstructed_speed_ratio > hllc_reconstructed_speed_ratio) {
        if (trace_this_face) {
          printf("[AFD-FACE-TRACE] branch=reconstructed_speed_guard "
                 "candidate_speed=%.17g stencil_speed=%.17g ratio=%.17g "
                 "limit=%.17g\n",
                 double(reconstructed_signal_speed),
                 double(stencil_signal_speed),
                 double(reconstructed_speed_ratio),
                 double(hllc_reconstructed_speed_ratio));
        }
        Real muscl_pressure = low_order_pressure;
        if (muscl_hllc_flux(
                iv, dir, q, flx, cls, trace_this_face,
                &muscl_pressure)) {
          if (radial_axis_face) {
            for (int n = 0; n < cls_t::NCONS; ++n) {
              flx(iv, n) = Real(0.0);
            }
          }
          return muscl_pressure;
        }
        cell_centered_llf(iv, dir, q, flx, cls);
        if (radial_axis_face) {
          for (int n = 0; n < cls_t::NCONS; ++n) {
            flx(iv, n) = Real(0.0);
          }
        }
        return low_order_pressure;
      }
    }

    if (trace_this_face) {
      Real trace_sl = amrex::min(ul - cl, ur - cr);
      Real trace_sr = amrex::max(ul + cl, ur + cr);
      const Real trace_rp = std::sqrt(rr / rl);
      const Real trace_uroe =
          (ul + ur * trace_rp) / (Real(1.0) + trace_rp);
      const Real trace_croe =
          (cl + cr * trace_rp) / (Real(1.0) + trace_rp);
      trace_sl = amrex::min(trace_sl, trace_uroe - trace_croe);
      trace_sr = amrex::max(trace_sr, trace_uroe + trace_croe);
      const Real trace_sstar_denom =
          rl * (trace_sl - ul) - rr * (trace_sr - ur);
      const Real trace_sstar_numer =
          pr - pl + rl * ul * (trace_sl - ul) -
          rr * ur * (trace_sr - ur);
      const Real trace_sstar = trace_sstar_numer / trace_sstar_denom;
      printf("[AFD-FACE-TRACE] reconstructed_left="
             "(rho=%.17g,un=%.17g,ut1=%.17g,ut2=%.17g,p=%.17g,gamma=%.17g,c=%.17g,EperMass=%.17g)\n",
             double(rl), double(ul), double(ut1l), double(ut2l), double(pl),
             double(gl), double(cl), double(el));
      printf("[AFD-FACE-TRACE] reconstructed_right="
             "(rho=%.17g,un=%.17g,ut1=%.17g,ut2=%.17g,p=%.17g,gamma=%.17g,c=%.17g,EperMass=%.17g)\n",
             double(rr), double(ur), double(ut1r), double(ut2r), double(pr),
             double(gr), double(cr), double(er));
      printf("[AFD-FACE-TRACE] hllc_speeds="
             "(sl=%.17g,sr=%.17g,sstar_numer=%.17g,sstar_denom=%.17g,sstar=%.17g,"
             "sl_minus_sstar=%.17g,sr_minus_sstar=%.17g)\n",
             double(trace_sl), double(trace_sr), double(trace_sstar_numer),
             double(trace_sstar_denom), double(trace_sstar),
             double(trace_sl - trace_sstar),
             double(trace_sr - trace_sstar));
    }

    Real fn = Real(0.0);
    Real ft1 = Real(0.0);
    Real ft2 = Real(0.0);
    Real fe = Real(0.0);
    Real fry[1] = {Real(0.0)};
    Real pressure_face = low_order_pressure;
    this->hllc(rl, ul, pl, ut1l, ut2l, el, yl, cl,
               rr, ur, pr, ut1r, ut2r, er, yr, cr,
               fn, ft1, ft2, fe, fry, &pressure_face);
    if (!std::isfinite(fn + ft1 + ft2 + fe + fry[0])) {
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] branch=nonfinite_hllc_cell_llf\n");
      }
      cell_centered_llf(iv, dir, q, flx, cls);
      if (radial_axis_face) {
        for (int n = 0; n < cls_t::NCONS; ++n) {
          flx(iv, n) = Real(0.0);
        }
      }
      return low_order_pressure;
    }

    flx(iv, cls_t::UMX + dir) = fn;
    flx(iv, cls_t::UMX + t1) = ft1;
    flx(iv, cls_t::UMX + t2) = ft2;
    flx(iv, cls_t::UET) = fe;
    flx(iv, cls_t::URHO) = fry[0];

    if (trace_this_face) {
      printf("[AFD-FACE-TRACE] base_hllc_flux="
             "(rho=%.17g,mn=%.17g,mt1=%.17g,mt2=%.17g,E=%.17g)\n",
             double(fry[0]), double(fn), double(ft1), double(ft2), double(fe));
    }

    if (correction == 0 || !smooth) {
      if (radial_axis_face) {
        for (int n = 0; n < cls_t::NCONS; ++n) {
          flx(iv, n) = Real(0.0);
        }
      }
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] final_branch=base_hllc "
               "afd_correction_applied=0 reason=%s\n",
               correction == 0 ? "disabled" : "nonsmooth_stencil");
      }
      return pressure_face;
    }

    // In the radial direction the AFD correction must act on the metric
    // advective flux r*(F-P), while pressure uses its ordinary Cartesian
    // h-function.  Correcting the complete r*F and only afterwards
    // subtracting the ordinary pressure companion leaves a pressure
    // commutator in F-P; at the first ring that already fails the r^2 even
    // moment.  The signed cell radii below also give the regular parity
    // extension needed by the axis stencil.
    Real pressure_samples[6] = {};
    if (rz_radial_direction) {
      for (int m = 0; m < 6; ++m) {
        pressure_samples[m] = pressure_states(
            iv + (m - 3) * ivd, pressure_component);
      }
    }
    Real f[6][cls_t::NCONS];
    for (int m = 0; m < 6; ++m) {
      const amrex::IntVect sample = iv + (m - 3) * ivd;
      cls.prims2flux(sample, dir, q, f[m]);
      if (rz_radial_direction) {
        const Real cell_radius =
            face_radius + (Real(m) - Real(2.5)) * radial_spacing;
        for (int n = 0; n < cls_t::NCONS; ++n) {
          if (n == cls_t::UMX) {
            f[m][n] -= pressure_samples[m];
          }
          f[m][n] *= cell_radius;
        }
      }
    }

    // The companion pressure approximates the ordinary derivative dp/dr,
    // so its AFD correction is deliberately nonmetric C[p].  The radial
    // momentum entry in f stores r*(F-P), whereas every other entry stores
    // its complete metric flux.  Combining them below gives
    // D_metric(F-P)+D_cartesian(P), and the C[p] linear limit pairs exactly
    // with the H6/H4/H2 axis closure.
    Real pressure_correction = Real(0.0);
    if (rz_radial_direction) {
      pressure_correction =
          -Real(1.0 / 1152.0) *
              (-Real(5.0) * pressure_samples[0] +
               Real(39.0) * pressure_samples[1] -
               Real(34.0) * pressure_samples[2] -
               Real(34.0) * pressure_samples[3] +
               Real(39.0) * pressure_samples[4] -
               Real(5.0) * pressure_samples[5]) +
          Real(7.0 / 11520.0) *
              (pressure_samples[0] - Real(3.0) * pressure_samples[1] +
               Real(2.0) * pressure_samples[2] +
               Real(2.0) * pressure_samples[3] -
               Real(3.0) * pressure_samples[4] + pressure_samples[5]);
    }
    for (int n = 0; n < cls_t::NCONS; ++n) {
      const Real correction_flux =
          -Real(1.0 / 1152.0) *
              (-Real(5.0) * f[0][n] + Real(39.0) * f[1][n] -
               Real(34.0) * f[2][n] - Real(34.0) * f[3][n] +
               Real(39.0) * f[4][n] - Real(5.0) * f[5][n]) +
          Real(7.0 / 11520.0) *
              (f[0][n] - Real(3.0) * f[1][n] + Real(2.0) * f[2][n] +
               Real(2.0) * f[3][n] - Real(3.0) * f[4][n] + f[5][n]);

      if (rz_radial_direction) {
        if (n == cls_t::UMX) {
          const Real pressure_auxiliary_flux =
              pressure_face + pressure_correction;
          const Real metric_advective_flux =
              face_radius * (flx(iv, n) - pressure_face) +
              correction_flux;
          flx(iv, n) = face_radius > Real(0.0)
                           ? metric_advective_flux / face_radius +
                                 pressure_auxiliary_flux
                           : metric_advective_flux;
        } else {
          const Real metric_flux =
              face_radius * flx(iv, n) + correction_flux;
          flx(iv, n) = face_radius > Real(0.0)
                           ? metric_flux / face_radius
                           : metric_flux;
        }
      } else {
        flx(iv, n) += correction_flux;
      }
      if (trace_this_face) {
        printf("[AFD-FACE-TRACE] correction component=%d value=%.17g "
               "final_flux=%.17g\n",
               n, double(correction_flux), double(flx(iv, n)));
      }
    }
    return pressure_face + pressure_correction;
  }
};

template <typename cls_t>
using afd_hllc_teno5_t = afd_hllc_t<AfdReconScheme::Teno5, cls_t>;

template <typename cls_t>
using afd_hllc_wenoz5_t = afd_hllc_t<AfdReconScheme::WenoZ5, cls_t>;

#endif
