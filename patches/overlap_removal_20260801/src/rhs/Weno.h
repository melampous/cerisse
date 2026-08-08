#ifndef CERISSE_WENO_H_
#define CERISSE_WENO_H_

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

/*
 * Compact characteristic local-LLF reconstruction family.
 *
 * The common driver retains the conservative finite-difference sequence used
 * by the historical CERISSE implementation:
 *
 *   local scalar LLF splitting
 *   face-local Roe characteristic projection
 *   nonlinear reconstruction of both split fluxes
 *   inverse characteristic projection
 *
 * WENO-Z5, TENO5, and TENO6 are reconstruction policies.  The only immersed
 * boundary adapter is the pure shared-GP stencil policy.  No cut-cell,
 * embedded-boundary, face-local wall-flux, or direct-wall-flux path is
 * implemented here.
 */

namespace ReconScheme {

using amrex::Real;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Real
square(const Real value) noexcept
{
  return value * value;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Real
power6(const Real value) noexcept
{
  const Real value2 = value * value;
  return value2 * value2 * value2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Real
weno_epsilon() noexcept
{
  return std::numeric_limits<Real>::digits >= 53
             ? Real(1.0e-40)
             : std::numeric_limits<Real>::epsilon();
}

template <unsigned int NCAND>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
select_teno_candidates(
    const bool active[NCAND], const Real beta[NCAND], const Real tau,
    const Real cutoff, bool selected[NCAND]) noexcept
{
  constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
  Real separation[NCAND] = {};
  Real maximum = Real(0.0);
  bool has_infinite = false;

  if (!(tau >= Real(0.0)) || !std::isfinite(tau)) {
    return false;
  }

  for (int candidate = 0; candidate < NCAND; ++candidate) {
    selected[candidate] = false;
    if (!active[candidate]) {
      continue;
    }
    if (!(beta[candidate] >= Real(0.0)) ||
        !std::isfinite(beta[candidate])) {
      return false;
    }

    separation[candidate] =
        Real(1.0) + tau / (epsilon + beta[candidate]);
    if (std::isinf(separation[candidate])) {
      has_infinite = true;
    } else if (!std::isfinite(separation[candidate]) ||
               !(separation[candidate] > Real(0.0))) {
      return false;
    } else {
      maximum = amrex::max(maximum, separation[candidate]);
    }
  }

  Real gamma[NCAND] = {};
  Real gamma_sum = Real(0.0);
  for (int candidate = 0; candidate < NCAND; ++candidate) {
    if (!active[candidate]) {
      continue;
    }
    if (has_infinite) {
      gamma[candidate] =
          std::isinf(separation[candidate]) ? Real(1.0) : Real(0.0);
    } else {
      if (!(maximum > Real(0.0))) {
        return false;
      }
      gamma[candidate] = power6(separation[candidate] / maximum);
    }
    gamma_sum += gamma[candidate];
  }

  if (!(gamma_sum > std::numeric_limits<Real>::min()) ||
      !std::isfinite(gamma_sum)) {
    return false;
  }

  const Real inverse_gamma_sum = Real(1.0) / gamma_sum;
  bool any_selected = false;
  for (int candidate = 0; candidate < NCAND; ++candidate) {
    selected[candidate] =
        active[candidate] &&
        gamma[candidate] * inverse_gamma_sum >= cutoff;
    any_selected = any_selected || selected[candidate];
  }
  return any_selected;
}

template <unsigned int NCAND>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
linear_candidate_fallback(
    const bool active[NCAND], const Real candidate_value[NCAND],
    const amrex::GpuArray<Real, NCAND>& optimal_weight) noexcept
{
  Real weighted_value = Real(0.0);
  Real weight_sum = Real(0.0);
  for (int candidate = 0; candidate < NCAND; ++candidate) {
    if (!active[candidate] ||
        !std::isfinite(candidate_value[candidate])) {
      continue;
    }
    weighted_value +=
        optimal_weight[candidate] * candidate_value[candidate];
    weight_sum += optimal_weight[candidate];
  }
  if (weight_sum > std::numeric_limits<Real>::min() &&
      std::isfinite(weighted_value + weight_sum)) {
    return weighted_value / (Real(6.0) * weight_sum);
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

struct WenoZ5 {
  static constexpr int ng = 3;
  static constexpr int stencil_size = 5;
  static constexpr int ncand = 3;
  static constexpr bool uses_cutoff = false;
  static constexpr Real default_cutoff = Real(0.0);
  static constexpr const char* cutoff_key = "";
  static constexpr const char* scheme_name =
      "characteristic_llf_weno_z5";

  using WeightArray = amrex::GpuArray<Real, ncand>;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr WeightArray
  optimal_weights() noexcept
  {
    return {Real(3.0), Real(6.0), Real(1.0)};
  }

  /*
   * The point indices refer to the six LLF samples stored by the common
   * driver.  Candidate numbering follows beta[0], beta[1], and beta[2].
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
  candidate_uses_point(
      const bool positive_split, const int candidate,
      const int point) noexcept
  {
    if (positive_split) {
      if (candidate == 0) return point >= 2 && point <= 4;
      if (candidate == 1) return point >= 1 && point <= 3;
      return point >= 0 && point <= 2;
    }
    if (candidate == 0) return point >= 1 && point <= 3;
    if (candidate == 1) return point >= 2 && point <= 4;
    return point >= 3 && point <= 5;
  }

  template <std::size_t NCONS>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  left_stencil(
      const int component, const Real split_flux[2 * ng][NCONS],
      Real stencil[stencil_size]) noexcept
  {
    for (int point = 1; point < 2 * ng; ++point) {
      stencil[point - 1] =
          split_flux[2 * ng - 1 - point][component];
    }
  }

  template <std::size_t NCONS>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  right_stencil(
      const int component, const Real split_flux[2 * ng][NCONS],
      Real stencil[stencil_size]) noexcept
  {
    for (int point = 1; point < 2 * ng; ++point) {
      stencil[point - 1] = split_flux[point][component];
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  smoothness_indicator(
      const Real stencil[stencil_size], Real beta[ncand]) noexcept
  {
    beta[2] =
        Real(13.0 / 12.0) *
            square(stencil[4] - Real(2.0) * stencil[3] + stencil[2]) +
        Real(0.25) *
            square(stencil[4] - Real(4.0) * stencil[3] +
                   Real(3.0) * stencil[2]);
    beta[1] =
        Real(13.0 / 12.0) *
            square(stencil[3] - Real(2.0) * stencil[2] + stencil[1]) +
        Real(0.25) * square(stencil[3] - stencil[1]);
    beta[0] =
        Real(13.0 / 12.0) *
            square(stencil[2] - Real(2.0) * stencil[1] + stencil[0]) +
        Real(0.25) *
            square(Real(3.0) * stencil[2] -
                   Real(4.0) * stencil[1] + stencil[0]);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  linear_polynomial_recon(
      const Real stencil[stencil_size],
      Real candidate_value[ncand]) noexcept
  {
    candidate_value[2] =
        Real(11.0) * stencil[2] - Real(7.0) * stencil[3] +
        Real(2.0) * stencil[4];
    candidate_value[1] =
        -stencil[3] + Real(5.0) * stencil[2] +
        Real(2.0) * stencil[1];
    candidate_value[0] =
        Real(2.0) * stencil[2] + Real(5.0) * stencil[1] -
        stencil[0];
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  evaluate_candidate(
      const int candidate, const Real stencil[stencil_size],
      Real& beta, Real& candidate_value) noexcept
  {
    if (candidate == 2) {
      beta =
          Real(13.0 / 12.0) *
              square(stencil[4] - Real(2.0) * stencil[3] + stencil[2]) +
          Real(0.25) *
              square(stencil[4] - Real(4.0) * stencil[3] +
                     Real(3.0) * stencil[2]);
      candidate_value =
          Real(11.0) * stencil[2] - Real(7.0) * stencil[3] +
          Real(2.0) * stencil[4];
    } else if (candidate == 1) {
      beta =
          Real(13.0 / 12.0) *
              square(stencil[3] - Real(2.0) * stencil[2] + stencil[1]) +
          Real(0.25) * square(stencil[3] - stencil[1]);
      candidate_value =
          -stencil[3] + Real(5.0) * stencil[2] +
          Real(2.0) * stencil[1];
    } else {
      beta =
          Real(13.0 / 12.0) *
              square(stencil[2] - Real(2.0) * stencil[1] + stencil[0]) +
          Real(0.25) *
              square(Real(3.0) * stencil[2] -
                     Real(4.0) * stencil[1] + stencil[0]);
      candidate_value =
          Real(2.0) * stencil[2] + Real(5.0) * stencil[1] -
          stencil[0];
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(
      const Real stencil[stencil_size],
      const Real /*unused_cutoff*/) noexcept
  {
    Real beta[ncand];
    Real candidate_value[ncand];
    smoothness_indicator(stencil, beta);
    linear_polynomial_recon(stencil, candidate_value);

    const Real tau = std::abs(beta[2] - beta[0]);
    const WeightArray linear_weight = optimal_weights();
    Real nonlinear_weight[ncand];
    Real weight_sum = Real(0.0);
    for (int candidate = 0; candidate < ncand; ++candidate) {
      nonlinear_weight[candidate] =
          (Real(1.0) +
           tau / (weno_epsilon() + beta[candidate])) *
          linear_weight[candidate];
      weight_sum += nonlinear_weight[candidate];
    }

    if (!(weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(weight_sum)) {
      const bool active[ncand] = {true, true, true};
      return linear_candidate_fallback<ncand>(
          active, candidate_value, linear_weight);
    }

    return (nonlinear_weight[0] * candidate_value[0] +
            nonlinear_weight[1] * candidate_value[1] +
            nonlinear_weight[2] * candidate_value[2]) /
           (Real(6.0) * weight_sum);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(const Real stencil[stencil_size]) noexcept
  {
    return recon(stencil, default_cutoff);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(
      const Real stencil[stencil_size],
      const WeightArray& optimal_weight,
      const Real /*unused_cutoff*/) noexcept
  {
    bool active[ncand] = {};
    Real beta[ncand] = {};
    Real candidate_value[ncand] = {};
    int active_count = 0;
    for (int candidate = 0; candidate < ncand; ++candidate) {
      active[candidate] = optimal_weight[candidate] > Real(0.0);
      if (!active[candidate]) {
        continue;
      }
      evaluate_candidate(
          candidate, stencil, beta[candidate],
          candidate_value[candidate]);
      ++active_count;
    }

    if (active_count == 0) {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    Real tau = Real(0.0);
    if (active[0] && active[1] && active[2]) {
      tau = std::abs(beta[2] - beta[0]);
    } else {
      Real beta_min = std::numeric_limits<Real>::max();
      Real beta_max = Real(0.0);
      for (int candidate = 0; candidate < ncand; ++candidate) {
        if (!active[candidate]) continue;
        beta_min = amrex::min(beta_min, beta[candidate]);
        beta_max = amrex::max(beta_max, beta[candidate]);
      }
      tau = active_count > 1 ? beta_max - beta_min : Real(0.0);
    }

    Real nonlinear_weight[ncand] = {};
    Real weight_sum = Real(0.0);
    for (int candidate = 0; candidate < ncand; ++candidate) {
      if (!active[candidate]) continue;
      nonlinear_weight[candidate] =
          (Real(1.0) +
           tau / (weno_epsilon() + beta[candidate])) *
          optimal_weight[candidate];
      weight_sum += nonlinear_weight[candidate];
    }

    if (!(weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(weight_sum)) {
      return linear_candidate_fallback<ncand>(
          active, candidate_value, optimal_weight);
    }

    Real reconstructed = Real(0.0);
    for (int candidate = 0; candidate < ncand; ++candidate) {
      reconstructed +=
          nonlinear_weight[candidate] * candidate_value[candidate];
    }
    return reconstructed / (Real(6.0) * weight_sum);
  }
};

/*
 * This policy retains the historical CERISSE TENO5 global indicator and its
 * pre-2026 cutoff.  The indicator is a CERISSE variant rather than the
 * canonical tau5 = |beta_2-beta_0| definition.  Keeping it here permits a
 * controlled comparison without silently changing both the driver and the
 * reconstruction at once.
 */
struct Teno5 : public WenoZ5 {
  static constexpr bool uses_cutoff = true;
  static constexpr Real default_cutoff = Real(1.0e-4);
  static constexpr const char* cutoff_key = "llf_teno_cutoff";
  static constexpr const char* scheme_name =
      "characteristic_llf_teno5";

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(const Real stencil[stencil_size], const Real cutoff) noexcept
  {
    return recon(stencil, optimal_weights(), cutoff);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(const Real stencil[stencil_size]) noexcept
  {
    return recon(stencil, default_cutoff);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(
      const Real stencil[stencil_size],
      const WeightArray& optimal_weight, const Real cutoff) noexcept
  {
    bool active[ncand] = {};
    Real beta[ncand] = {};
    Real candidate_value[ncand] = {};
    int active_count = 0;
    for (int candidate = 0; candidate < ncand; ++candidate) {
      active[candidate] = optimal_weight[candidate] > Real(0.0);
      if (!active[candidate]) {
        continue;
      }
      WenoZ5::evaluate_candidate(
          candidate, stencil, beta[candidate],
          candidate_value[candidate]);
      ++active_count;
    }

    if (active_count == 0) {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    Real tau = Real(0.0);
    if (active[0] && active[1] && active[2]) {
      tau = std::abs(
          std::abs(beta[2] - beta[0]) -
          (beta[2] + Real(4.0) * beta[1] + beta[0]) /
              Real(6.0));
    } else {
      Real beta_min = std::numeric_limits<Real>::max();
      Real beta_max = Real(0.0);
      for (int candidate = 0; candidate < ncand; ++candidate) {
        if (!active[candidate]) continue;
        beta_min = amrex::min(beta_min, beta[candidate]);
        beta_max = amrex::max(beta_max, beta[candidate]);
      }
      tau = active_count > 1 ? beta_max - beta_min : Real(0.0);
    }

    bool selected[ncand] = {};
    if (!select_teno_candidates<ncand>(
            active, beta, tau, cutoff, selected)) {
      return WenoZ5::recon(
          stencil, optimal_weight, Real(0.0));
    }

    Real selected_weight_sum = Real(0.0);
    Real reconstructed = Real(0.0);
    for (int candidate = 0; candidate < ncand; ++candidate) {
      if (!selected[candidate]) continue;
      selected_weight_sum += optimal_weight[candidate];
      reconstructed +=
          optimal_weight[candidate] * candidate_value[candidate];
    }
    if (!(selected_weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(reconstructed + selected_weight_sum)) {
      return WenoZ5::recon(
          stencil, optimal_weight, Real(0.0));
    }
    return reconstructed / (Real(6.0) * selected_weight_sum);
  }
};

/*
 * The full-stencil policy follows the six-point TENO reconstruction.  Its
 * three-point indicators are aligned with the candidate supports.  This
 * corrects a one-sample shift in the historical CERISSE implementation.
 */
struct Teno6 {
  static constexpr int ng = 3;
  static constexpr int stencil_size = 6;
  static constexpr int ncand = 4;
  static constexpr bool uses_cutoff = true;
  static constexpr Real default_cutoff = Real(1.0e-5);
  static constexpr const char* cutoff_key = "llf_teno6_cutoff";
  static constexpr const char* scheme_name =
      "characteristic_llf_teno6";

  using WeightArray = amrex::GpuArray<Real, ncand>;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr WeightArray
  optimal_weights() noexcept
  {
    /*
     * Candidate order is downwind three-point, centred three-point,
     * upwind three-point, and four-point downwind.  The normalized weights
     * are 6:9:1:4.  Their smooth combination gives the sixth-order central
     * finite-difference flux [1,-8,37,37,-8,1]/60.
     */
    return {Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
  candidate_uses_point(
      const bool positive_split, const int candidate,
      const int point) noexcept
  {
    if (positive_split) {
      if (candidate == 0) return point >= 2 && point <= 4;
      if (candidate == 1) return point >= 1 && point <= 3;
      if (candidate == 2) return point >= 0 && point <= 2;
      return point >= 2 && point <= 5;
    }
    if (candidate == 0) return point >= 1 && point <= 3;
    if (candidate == 1) return point >= 2 && point <= 4;
    if (candidate == 2) return point >= 3 && point <= 5;
    return point >= 0 && point <= 3;
  }

  template <std::size_t NCONS>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  left_stencil(
      const int component, const Real split_flux[2 * ng][NCONS],
      Real stencil[stencil_size]) noexcept
  {
    for (int point = 0; point < 2 * ng; ++point) {
      stencil[point] =
          split_flux[2 * ng - 1 - point][component];
    }
  }

  template <std::size_t NCONS>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  right_stencil(
      const int component, const Real split_flux[2 * ng][NCONS],
      Real stencil[stencil_size]) noexcept
  {
    for (int point = 0; point < 2 * ng; ++point) {
      stencil[point] = split_flux[point][component];
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  evaluate_candidate(
      const int candidate, const Real stencil[stencil_size],
      Real& beta, Real& candidate_value) noexcept
  {
    /*
     * Each smoothness indicator uses exactly the support of its candidate.
     * This is the Fu--Hu--Adams TENO6 stencil alignment.  The historical
     * CERISSE Teno6 code shifted the three three-point indicators by one
     * sample relative to their candidate polynomials.
     */
    if (candidate == 0) {
      beta =
          Real(13.0 / 12.0) *
              square(stencil[1] - Real(2.0) * stencil[2] + stencil[3]) +
          Real(0.25) *
              square(stencil[1] - Real(4.0) * stencil[2] +
                     Real(3.0) * stencil[3]);
      candidate_value =
          Real(2.0) * stencil[3] + Real(5.0) * stencil[2] -
          stencil[1];
    } else if (candidate == 1) {
      beta =
          Real(13.0 / 12.0) *
              square(stencil[2] - Real(2.0) * stencil[3] + stencil[4]) +
          Real(0.25) * square(stencil[2] - stencil[4]);
      candidate_value =
          -stencil[4] + Real(5.0) * stencil[3] +
          Real(2.0) * stencil[2];
    } else if (candidate == 2) {
      beta =
          Real(13.0 / 12.0) *
              square(stencil[3] - Real(2.0) * stencil[4] + stencil[5]) +
          Real(0.25) *
              square(Real(3.0) * stencil[3] -
                     Real(4.0) * stencil[4] + stencil[5]);
      candidate_value =
          Real(11.0) * stencil[3] - Real(7.0) * stencil[4] +
          Real(2.0) * stencil[5];
    } else {
      beta =
          Real(1.0 / 36.0) *
              square(
                  -Real(11.0) * stencil[3] +
                  Real(18.0) * stencil[2] -
                  Real(9.0) * stencil[1] +
                  Real(2.0) * stencil[0]) +
          Real(13.0 / 12.0) *
              square(
                  Real(2.0) * stencil[3] -
                  Real(5.0) * stencil[2] +
                  Real(4.0) * stencil[1] - stencil[0]) +
          Real(781.0 / 720.0) *
              square(
                  -stencil[3] + Real(3.0) * stencil[2] -
                  Real(3.0) * stencil[1] + stencil[0]);
      candidate_value =
          Real(0.5) *
          (stencil[0] - Real(5.0) * stencil[1] +
           Real(13.0) * stencil[2] + Real(3.0) * stencil[3]);
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(const Real stencil[stencil_size], const Real cutoff) noexcept
  {
    return recon(stencil, optimal_weights(), cutoff);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(const Real stencil[stencil_size]) noexcept
  {
    return recon(stencil, default_cutoff);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  recon(
      const Real stencil[stencil_size],
      const WeightArray& optimal_weight, const Real cutoff) noexcept
  {
    bool active[ncand] = {};
    Real beta[ncand] = {};
    Real candidate_value[ncand] = {};
    int active_count = 0;
    for (int candidate = 0; candidate < ncand; ++candidate) {
      active[candidate] = optimal_weight[candidate] > Real(0.0);
      if (!active[candidate]) {
        continue;
      }
      evaluate_candidate(
          candidate, stencil, beta[candidate],
          candidate_value[candidate]);
      ++active_count;
    }

    if (active_count == 0) {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    Real tau = Real(0.0);
    if (active[0] && active[1] && active[2] && active[3]) {
      tau = std::abs(
          beta[3] -
          (beta[2] + Real(4.0) * beta[1] + beta[0]) /
              Real(6.0));
    } else {
      /*
       * The full TENO6 reference indicator requires all four candidates.
       * Near a shared GP, use the range of the available indicators only for
       * ENO-like selection.  It never reads a disabled solid-side candidate.
       */
      Real beta_min = std::numeric_limits<Real>::max();
      Real beta_max = Real(0.0);
      for (int candidate = 0; candidate < ncand; ++candidate) {
        if (!active[candidate]) continue;
        beta_min = amrex::min(beta_min, beta[candidate]);
        beta_max = amrex::max(beta_max, beta[candidate]);
      }
      tau = active_count > 1 ? beta_max - beta_min : Real(0.0);
    }

    bool selected[ncand] = {};
    if (!select_teno_candidates<ncand>(
            active, beta, tau, cutoff, selected)) {
      return linear_candidate_fallback<ncand>(
          active, candidate_value, optimal_weight);
    }

    Real selected_weight_sum = Real(0.0);
    Real reconstructed = Real(0.0);
    for (int candidate = 0; candidate < ncand; ++candidate) {
      if (!selected[candidate]) continue;
      selected_weight_sum += optimal_weight[candidate];
      reconstructed +=
          optimal_weight[candidate] * candidate_value[candidate];
    }
    if (!(selected_weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(reconstructed + selected_weight_sum)) {
      return linear_candidate_fallback<ncand>(
          active, candidate_value, optimal_weight);
    }
    return reconstructed / (Real(6.0) * selected_weight_sum);
  }
};

template <typename Scheme>
struct CandidateMask {
  using WeightArray = typename Scheme::WeightArray;

  WeightArray left_weight{};
  WeightArray right_weight{};
  bool point_needed[2 * Scheme::ng] = {};
  bool high_order_available = false;
  bool full_stencil = false;
};

template <typename Scheme>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE CandidateMask<Scheme>
build_candidate_mask(const bool usable[2 * Scheme::ng]) noexcept
{
  CandidateMask<Scheme> mask;
  mask.left_weight = Scheme::optimal_weights();
  mask.right_weight = Scheme::optimal_weights();

  bool has_left = false;
  bool has_right = false;
  bool all_left = true;
  bool all_right = true;

  for (int candidate = 0; candidate < Scheme::ncand; ++candidate) {
    bool left_active = true;
    bool right_active = true;
    for (int point = 0; point < 2 * Scheme::ng; ++point) {
      if (Scheme::candidate_uses_point(true, candidate, point) &&
          !usable[point]) {
        left_active = false;
      }
      if (Scheme::candidate_uses_point(false, candidate, point) &&
          !usable[point]) {
        right_active = false;
      }
    }

    if (!left_active) {
      mask.left_weight[candidate] = Real(0.0);
      all_left = false;
    } else {
      has_left = true;
      for (int point = 0; point < 2 * Scheme::ng; ++point) {
        if (Scheme::candidate_uses_point(true, candidate, point)) {
          mask.point_needed[point] = true;
        }
      }
    }

    if (!right_active) {
      mask.right_weight[candidate] = Real(0.0);
      all_right = false;
    } else {
      has_right = true;
      for (int point = 0; point < 2 * Scheme::ng; ++point) {
        if (Scheme::candidate_uses_point(false, candidate, point)) {
          mask.point_needed[point] = true;
        }
      }
    }
  }

  mask.high_order_available = has_left && has_right;
  mask.full_stencil = all_left && all_right;
  return mask;
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_solid(
    const amrex::IntVect& cell, const MarkerArray& marker) noexcept
{
  return marker(cell, 0) != std::uint8_t(0);
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_shared_gp_usable(
    const amrex::IntVect& cell, const MarkerArray& marker) noexcept
{
  return marker(cell, 0) == std::uint8_t(0) ||
         marker(cell, 1) != std::uint8_t(0);
}

template <typename PrimitiveArray, typename Closure>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
adjacent_states_admissible(
    const amrex::IntVect& face, const amrex::IntVect& direction,
    const PrimitiveArray& primitive) noexcept
{
  const amrex::IntVect left = face - direction;
  const Real finite_sum =
      primitive(left, Closure::QRHO) +
      primitive(face, Closure::QRHO) +
      primitive(left, Closure::QPRES) +
      primitive(face, Closure::QPRES) +
      primitive(left, Closure::QC) +
      primitive(face, Closure::QC);
  return primitive(left, Closure::QRHO) > Real(0.0) &&
         primitive(face, Closure::QRHO) > Real(0.0) &&
         primitive(left, Closure::QPRES) > Real(0.0) &&
         primitive(face, Closure::QPRES) > Real(0.0) &&
         primitive(left, Closure::QC) > Real(0.0) &&
         primitive(face, Closure::QC) > Real(0.0) &&
         std::isfinite(finite_sum);
}

template <typename FluxArray, typename Closure>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
set_zero_flux(
    const amrex::IntVect& face, const FluxArray& flux) noexcept
{
  for (int component = 0; component < Closure::NCONS; ++component) {
    flux(face, component) = Real(0.0);
  }
}

template <typename FluxArray, typename Closure>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
set_invalid_flux(
    const amrex::IntVect& face, const FluxArray& flux) noexcept
{
  const Real invalid = std::numeric_limits<Real>::quiet_NaN();
  for (int component = 0; component < Closure::NCONS; ++component) {
    flux(face, component) = invalid;
  }
}

}  // namespace ReconScheme

template <typename Scheme, typename cls_t>
class weno_t {
  using CandidateMask = ReconScheme::CandidateMask<Scheme>;

  static constexpr int ng = Scheme::ng;
  static_assert(ng <= cls_t::NGHOST);
  static_assert(Scheme::stencil_size <= 2 * ng);

  static amrex::Real configured_reconstruction_cutoff()
  {
    static const amrex::Real value = [] {
      using amrex::Real;
      if constexpr (Scheme::uses_cutoff) {
        Real configured = Scheme::default_cutoff;
        amrex::ParmParse parameters("cns");
        parameters.query(Scheme::cutoff_key, configured);
        if (!std::isfinite(configured) ||
            !(configured > Real(0.0)) || !(configured < Real(1.0))) {
          amrex::Abort(
              "LLF TENO cutoff must be finite and lie in (0,1)");
        }
        return configured;
      } else {
        return Real(0.0);
      }
    }();
    return value;
  }

 public:
  static constexpr const char* scheme_name = Scheme::scheme_name;

  static constexpr bool rz_annular_pressure_consistency_capable = false;
  static constexpr bool rz_euler_geometric_source_active = true;

  AMREX_GPU_HOST_DEVICE
  weno_t() = default;

  void eflux(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims_in,
      const std::array<amrex::FArrayBox*, AMREX_SPACEDIM>& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls)
  {
    eflux(geom, mfi.tilebox(), prims_in, flxt, rhs, cls);
  }

  void eflux(
      const amrex::Geometry& geom, const amrex::Box& bx,
      const amrex::Array4<const amrex::Real>& prims_in,
      const std::array<amrex::FArrayBox*, AMREX_SPACEDIM>& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
      const amrex::Box& skip_cells = amrex::Box())
  {
    amrex::Array4<std::uint8_t> unused_marker;
    compute_flux<false>(
        geom, bx, prims_in, flxt, rhs, cls, unused_marker, skip_cells);
  }

  void eflux_ibm(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims_in,
      const std::array<amrex::FArrayBox*, AMREX_SPACEDIM>& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
      const amrex::Array4<std::uint8_t>& marker)
  {
#ifdef AMREX_USE_GPIBM
    if constexpr (
        std::is_same_v<Scheme, ReconScheme::WenoZ5>) {
      compute_flux<true>(
          geom, mfi.tilebox(), prims_in, flxt, rhs, cls, marker,
          amrex::Box());
    } else {
      amrex::ignore_unused(
          geom, mfi, prims_in, flxt, rhs, cls, marker);
      amrex::Abort(
          "TENO schemes are all-fluid verification options; "
          "the frozen pure shared-GP production scheme is LLF-WENO-Z5");
    }
#else
    amrex::ignore_unused(
        geom, mfi, prims_in, flxt, rhs, cls, marker);
    amrex::Abort(
        "weno_t::eflux_ibm requires the pure shared-GP IBM build");
#endif
  }

  template <bool UseSharedGP, typename MarkerArray>
  static void compute_flux(
      const amrex::Geometry& geom, const amrex::Box& bx,
      const amrex::Array4<const amrex::Real>& prims_in,
      const std::array<amrex::FArrayBox*, AMREX_SPACEDIM>& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
      const MarkerArray& marker, const amrex::Box& skip_cells)
  {
    using amrex::IntVect;
    using amrex::Real;
    using namespace ReconScheme;

    amrex::ignore_unused(geom, rhs);

    const Real reconstruction_cutoff =
        configured_reconstruction_cutoff();

    const bool skip_valid = skip_cells.ok();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const amrex::Box flux_box = amrex::surroundingNodes(bx, dir);
      const auto flux = flxt[dir]->array();

      amrex::ParallelFor(
          flux_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const IntVect face(AMREX_D_DECL(i, j, k));
            const IntVect direction(
                IntVect::TheDimensionVector(dir));
            amrex::ignore_unused(marker, flux);

            if (skip_valid &&
                skip_cells.contains(face) &&
                skip_cells.contains(face - direction)) {
              return;
            }

            bool usable[2 * ng];
            for (int point = 0; point < 2 * ng; ++point) {
              usable[point] = true;
            }

            if constexpr (UseSharedGP) {
              if (is_solid(face - direction, marker) &&
                  is_solid(face, marker)) {
                set_zero_flux<decltype(flux), cls_t>(face, flux);
                return;
              }
              for (int point = 0; point < 2 * ng; ++point) {
                const IntVect cell =
                    face + (point - ng) * direction;
                usable[point] =
                    is_shared_gp_usable(cell, marker);
              }
            }

            const CandidateMask mask =
                build_candidate_mask<Scheme>(usable);
            const bool adjacent_states_usable =
                usable[ng - 1] && usable[ng];
            const bool adjacent_states_valid =
                adjacent_states_usable &&
                adjacent_states_admissible<
                    decltype(prims_in), cls_t>(
                    face, direction, prims_in);
            AMREX_ASSERT_WITH_MESSAGE(
                adjacent_states_valid,
                "LLF reconstruction found an invalid adjacent state");
            if (!adjacent_states_valid) {
              set_invalid_flux<decltype(flux), cls_t>(face, flux);
              return;
            }

            Real alpha = Real(0.0);
            if (mask.full_stencil) {
              alpha =
                  cls->max_char_speed(face, dir, ng, prims_in);
            } else {
              for (int point = 0; point < 2 * ng; ++point) {
                const bool read_point =
                    mask.high_order_available
                        ? mask.point_needed[point]
                        : (point == ng - 1 || point == ng);
                if (!read_point) continue;

                const IntVect cell =
                    face + (point - ng) * direction;
                const Real wave_speed =
                    std::abs(
                        prims_in(cell, cls_t::QU + dir)) +
                    prims_in(cell, cls_t::QC);
                alpha = amrex::max(alpha, wave_speed);
              }
            }

            AMREX_ASSERT_WITH_MESSAGE(
                alpha > Real(0.0),
                "LLF reconstruction found a non-positive wave speed");
            if (!(alpha > Real(0.0)) || !std::isfinite(alpha)) {
              set_invalid_flux<decltype(flux), cls_t>(face, flux);
              return;
            }

            if (!mask.high_order_available) {
              Real conservative_left[cls_t::NCONS];
              Real conservative_right[cls_t::NCONS];
              Real physical_flux_left[cls_t::NCONS];
              Real physical_flux_right[cls_t::NCONS];

              cls->prims2flux(
                  face - direction, dir, prims_in,
                  physical_flux_left);
              cls->prims2flux(
                  face, dir, prims_in, physical_flux_right);
              cls->prims2cons(
                  face - direction, prims_in,
                  conservative_left);
              cls->prims2cons(
                  face, prims_in, conservative_right);

              for (int component = 0;
                   component < cls_t::NCONS; ++component) {
                flux(face, component) =
                    Real(0.5) *
                        (physical_flux_left[component] +
                         physical_flux_right[component]) -
                    Real(0.5) * alpha *
                        (conservative_right[component] -
                         conservative_left[component]);
              }
              return;
            }

            const auto roe_average =
                cls->roe_avg_state(face, dir, prims_in);
            Real conservative[cls_t::NCONS];
            Real physical_flux[cls_t::NCONS];
            Real positive_split[2 * ng][cls_t::NCONS] = {};
            Real negative_split[2 * ng][cls_t::NCONS] = {};

            for (int point = 0; point < 2 * ng; ++point) {
              if (!mask.point_needed[point]) continue;

              const IntVect cell =
                  face + (point - ng) * direction;
              cls->prims2flux(
                  cell, dir, prims_in, physical_flux);
              cls->prims2cons(
                  cell, prims_in, conservative);

              for (int component = 0;
                   component < cls_t::NCONS; ++component) {
                positive_split[point][component] =
                    Real(0.5) *
                    (conservative[component] +
                     physical_flux[component] / alpha);
                negative_split[point][component] =
                    Real(0.5) *
                    (conservative[component] -
                     physical_flux[component] / alpha);
              }
              cls->cons2char(
                  roe_average, positive_split[point]);
              cls->cons2char(
                  roe_average, negative_split[point]);
            }

            Real positive_face[cls_t::NCONS];
            Real negative_face[cls_t::NCONS];
            Real stencil[Scheme::stencil_size];
            for (int component = 0;
                 component < cls_t::NCONS; ++component) {
              Scheme::left_stencil(
                  component, positive_split, stencil);
              positive_face[component] =
                  mask.full_stencil
                      ? Scheme::recon(
                            stencil, reconstruction_cutoff)
                      : Scheme::recon(
                            stencil, mask.left_weight,
                            reconstruction_cutoff);

              Scheme::right_stencil(
                  component, negative_split, stencil);
              negative_face[component] =
                  mask.full_stencil
                      ? Scheme::recon(
                            stencil, reconstruction_cutoff)
                      : Scheme::recon(
                            stencil, mask.right_weight,
                            reconstruction_cutoff);
            }

            cls->char2cons(roe_average, positive_face);
            cls->char2cons(roe_average, negative_face);

            for (int component = 0;
                 component < cls_t::NCONS; ++component) {
              flux(face, component) =
                  alpha *
                  (positive_face[component] -
                   negative_face[component]);
            }
          });
    }
  }
};

template <typename cls_t>
using llf_wenoz5_t =
    weno_t<ReconScheme::WenoZ5, cls_t>;

template <typename cls_t>
using llf_teno5_t =
    weno_t<ReconScheme::Teno5, cls_t>;

template <typename cls_t>
using llf_teno6_t =
    weno_t<ReconScheme::Teno6, cls_t>;

#endif
