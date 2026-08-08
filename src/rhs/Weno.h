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

/*
 * Compact characteristic local Lax--Friedrichs (LLF) reconstruction family.
 *
 * The common driver retains the conservative finite-difference sequence used
 * by the historical CERISSE implementation:
 *
 *   local scalar LLF splitting
 *   face-local Roe characteristic projection
 *   nonlinear reconstruction of both split fluxes
 *   inverse characteristic projection
 *
 * WENO-Z5, TENO5, and TENO6 are reconstruction methods.  The immersed
 * boundary adapter uses only reconstructed ghost-point states.  No cut-cell,
 * embedded-boundary, face-local wall-flux, or direct-wall-flux path is
 * implemented here.
 */

/**
 * \brief Reconstruction methods used by the characteristic local LLF driver.
 *
 * A method supplies only the one-dimensional nonlinear reconstruction.  It
 * does not own flow data and does not compute the multidimensional flux
 * divergence.  The driver below supplies LLF splitting, characteristic
 * projection, face traversal, and the pure ghost-point IBM stencil rules.
 */
namespace Reconstruction {

using amrex::Real;

/// Return the square of a scalar without invoking a general power routine.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Real
square(const Real value) noexcept
{
  return value * value;
}

/// Return the sixth power used by the TENO scale-separation sensor.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Real
sixth_power(const Real value) noexcept
{
  const Real value_squared = value * value;
  return value_squared * value_squared * value_squared;
}

/// Return the regularisation used in the WENO smoothness ratios.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Real
smoothness_epsilon() noexcept
{
  return std::numeric_limits<Real>::digits >= 53
             ? Real(1.0e-40)
             : std::numeric_limits<Real>::epsilon();
}

// WENO-Z nonlinear-weight exponent q.  The production configuration uses
// q = 1 (Enson-compatible; bit-identical default).  Building with
// -DCERISSE_WENOZ_POWER=2 selects the squared smoothness ratio, which
// restores fifth order at first-order critical points (Borges et al. 2008)
// at the cost of stronger candidate suppression near discontinuities.
#ifndef CERISSE_WENOZ_POWER
#define CERISSE_WENOZ_POWER 1
#endif

// TENO5 candidate-selection cutoff.  The default preserves the production
// configuration; smooth-flow verification builds may override it explicitly.
#ifndef CERISSE_TENO5_SELECTION_CUTOFF
#define CERISSE_TENO5_SELECTION_CUTOFF 1.0e-4
#endif

// Reference indicator for TENO5: 0 preserves the historical CERISSE
// expression; 1 selects the standard tau_5 = abs(beta_0-beta_2).
#ifndef CERISSE_TENO5_STANDARD_TAU5
#define CERISSE_TENO5_STANDARD_TAU5 0
#endif

/// WENO-Z alpha factor 1 + (tau/(eps+beta))^q for the configured exponent.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
wenoz_alpha_factor(const Real global_smoothness_indicator,
                   const Real smoothness_indicator,
                   const Real epsilon) noexcept
{
#if (CERISSE_WENOZ_POWER == 2)
  const Real smoothness_ratio =
      global_smoothness_indicator /
      (epsilon + smoothness_indicator);
  return Real(1.0) + smoothness_ratio * smoothness_ratio;
#else
  return Real(1.0) +
         global_smoothness_indicator /
             (epsilon + smoothness_indicator);
#endif
}

/**
 * \brief Select the TENO candidates whose normalized sensor exceeds a cutoff.
 *
 * \tparam CandidateCount Number of candidate substencils.
 * \param candidate_available True when geometry permits a candidate.
 * \param smoothness_indicators Jiang--Shu indicator for each candidate.
 * \param global_smoothness_indicator TENO reference indicator.
 * \param selection_cutoff Fixed TENO acceptance threshold.
 * \param candidate_selected Output flag for every retained candidate.
 * \return True when at least one candidate is retained.
 */
template <unsigned int CandidateCount>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
select_teno_smooth_candidates(
    const bool candidate_available[CandidateCount],
    const Real smoothness_indicators[CandidateCount],
    const Real global_smoothness_indicator,
    const Real selection_cutoff,
    bool candidate_selected[CandidateCount]) noexcept
{
  // Machine epsilon prevents division by zero in the TENO sensor.
  constexpr Real roundoff_epsilon = std::numeric_limits<Real>::epsilon();
  // smoothness_ratio stores 1 + tau/(beta_k + epsilon).
  Real smoothness_ratio[CandidateCount] = {};
  // maximum_smoothness_ratio is used to evaluate the sixth power safely.
  Real maximum_smoothness_ratio = Real(0.0);
  // has_infinite_ratio records a mathematically dominant candidate.
  bool has_infinite_ratio = false;

  if (!(global_smoothness_indicator >= Real(0.0)) ||
      !std::isfinite(global_smoothness_indicator)) {
    return false;
  }

  for (int candidate_index = 0;
       candidate_index < CandidateCount;
       ++candidate_index) {
    candidate_selected[candidate_index] = false;
    if (!candidate_available[candidate_index]) {
      continue;
    }
    if (!(smoothness_indicators[candidate_index] >= Real(0.0)) ||
        !std::isfinite(smoothness_indicators[candidate_index])) {
      return false;
    }

    smoothness_ratio[candidate_index] =
        Real(1.0) +
        global_smoothness_indicator /
            (roundoff_epsilon + smoothness_indicators[candidate_index]);
    if (std::isinf(smoothness_ratio[candidate_index])) {
      has_infinite_ratio = true;
    } else if (!std::isfinite(smoothness_ratio[candidate_index]) ||
               !(smoothness_ratio[candidate_index] > Real(0.0))) {
      return false;
    } else {
      maximum_smoothness_ratio = amrex::max(
          maximum_smoothness_ratio, smoothness_ratio[candidate_index]);
    }
  }

  // teno_sensor is the unnormalized TENO value gamma_k.
  Real teno_sensor[CandidateCount] = {};
  // teno_sensor_sum normalizes the candidate sensors.
  Real teno_sensor_sum = Real(0.0);
  for (int candidate_index = 0;
       candidate_index < CandidateCount;
       ++candidate_index) {
    if (!candidate_available[candidate_index]) {
      continue;
    }
    if (has_infinite_ratio) {
      teno_sensor[candidate_index] =
          std::isinf(smoothness_ratio[candidate_index])
              ? Real(1.0)
              : Real(0.0);
    } else {
      if (!(maximum_smoothness_ratio > Real(0.0))) {
        return false;
      }
      teno_sensor[candidate_index] = sixth_power(
          smoothness_ratio[candidate_index] /
          maximum_smoothness_ratio);
    }
    teno_sensor_sum += teno_sensor[candidate_index];
  }

  if (!(teno_sensor_sum > std::numeric_limits<Real>::min()) ||
      !std::isfinite(teno_sensor_sum)) {
    return false;
  }

  // inverse_teno_sensor_sum avoids one division per candidate.
  const Real inverse_teno_sensor_sum = Real(1.0) / teno_sensor_sum;
  // has_selected_candidate is the success flag returned to the method.
  bool has_selected_candidate = false;
  for (int candidate_index = 0;
       candidate_index < CandidateCount;
       ++candidate_index) {
    candidate_selected[candidate_index] =
        candidate_available[candidate_index] &&
        teno_sensor[candidate_index] * inverse_teno_sensor_sum >=
            selection_cutoff;
    has_selected_candidate =
        has_selected_candidate || candidate_selected[candidate_index];
  }
  return has_selected_candidate;
}

/**
 * \brief Form the optimal linear combination of all available candidates.
 *
 * Candidate polynomial values in this file are stored with a factor of six.
 * The final division by six restores the interface value.  This routine is a
 * deterministic fallback for a failed nonlinear normalization.
 */
template <unsigned int CandidateCount>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
reconstruct_from_available_linear_candidates(
    const bool candidate_available[CandidateCount],
    const Real scaled_candidate_values[CandidateCount],
    const amrex::GpuArray<Real, CandidateCount>& linear_weights) noexcept
{
  // weighted_scaled_value accumulates weight_k times six times q_k.
  Real weighted_scaled_value = Real(0.0);
  // available_weight_sum normalizes only the available candidates.
  Real available_weight_sum = Real(0.0);
  for (int candidate_index = 0;
       candidate_index < CandidateCount;
       ++candidate_index) {
    if (!candidate_available[candidate_index] ||
        !std::isfinite(scaled_candidate_values[candidate_index])) {
      continue;
    }
    weighted_scaled_value +=
        linear_weights[candidate_index] *
        scaled_candidate_values[candidate_index];
    available_weight_sum += linear_weights[candidate_index];
  }
  if (available_weight_sum > std::numeric_limits<Real>::min() &&
      std::isfinite(weighted_scaled_value + available_weight_sum)) {
    return weighted_scaled_value /
           (Real(6.0) * available_weight_sum);
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

/** \brief Fifth-order WENO-Z reconstruction method. */
struct WenoZ5 {
  /// Number of cell layers required on each side of a face.
  static constexpr int required_ghost_cells = 3;
  /// Number of ordered samples used by one split-flux reconstruction.
  static constexpr int reconstruction_sample_count = 5;
  /// Number of three-point candidate substencils.
  static constexpr int candidate_count = 3;
  /// True because WENO-Z nonlinear weights use the supplied epsilon.
  static constexpr bool uses_weno_regularisation_epsilon = true;
  /// The complete WENO-Z weights may be frozen and reused by the paired RZ
  /// pressure reconstruction.  The same interface is implemented by TENO5,
  /// TENO6, and the marker-reduced pure shared-GP stencil plans below.
  static constexpr bool supports_rz_paired_pressure_flux = true;
  /// Method-manifest name exposed by the characteristic LLF driver.
  static constexpr const char* scheme_name =
      "characteristic_llf_weno_z5";

  /// Fixed-size storage for the candidate linear weights.
  using CandidateWeightArray = amrex::GpuArray<Real, candidate_count>;

  /**
   * \brief Unnormalized nonlinear weights frozen from one complete stencil.
   *
   * Keeping the common denominator preserves the arithmetic used by the
   * established WENO-Z reconstruction while allowing the same nonlinear
   * operator to be applied to the paired pressure field.
   */
  struct FrozenWeightSet {
    CandidateWeightArray candidate_weights{};
    Real weight_sum = Real(0.0);
  };

  /// Return the unnormalized WENO-Z linear weights 3:6:1.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr CandidateWeightArray
  optimal_linear_weights() noexcept
  {
    return {Real(3.0), Real(6.0), Real(1.0)};
  }

  /**
   * \brief Report whether a positive-branch candidate uses a driver sample.
   *
   * Sample indices refer to offsets -3,-2,-1,0,1,2 from the face index.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
  positive_candidate_uses_sample(
      const int candidate_index,
      const int sample_index) noexcept
  {
    if (candidate_index == 0) {
      return sample_index >= 2 && sample_index <= 4;
    }
    if (candidate_index == 1) {
      return sample_index >= 1 && sample_index <= 3;
    }
    return sample_index >= 0 && sample_index <= 2;
  }

  /**
   * \brief Report whether a negative-branch candidate uses a driver sample.
   *
   * Sample indices refer to offsets -3,-2,-1,0,1,2 from the face index.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
  negative_candidate_uses_sample(
      const int candidate_index,
      const int sample_index) noexcept
  {
    if (candidate_index == 0) {
      return sample_index >= 1 && sample_index <= 3;
    }
    if (candidate_index == 1) {
      return sample_index >= 2 && sample_index <= 4;
    }
    return sample_index >= 3 && sample_index <= 5;
  }

  /**
   * \brief Gather the five ordered samples for the positive LLF branch.
   *
   * \tparam ConservativeComponentCount Number of conserved flux components.
   */
  template <std::size_t ConservativeComponentCount>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  gather_positive_split_stencil(
      const int component_index,
      const Real split_samples[2 * required_ghost_cells]
                              [ConservativeComponentCount],
      Real stencil_values[reconstruction_sample_count]) noexcept
  {
    for (int sample_index = 1;
         sample_index < 2 * required_ghost_cells;
         ++sample_index) {
      stencil_values[sample_index - 1] =
          split_samples[2 * required_ghost_cells - 1 - sample_index]
                       [component_index];
    }
  }

  /**
   * \brief Gather the five ordered samples for the negative LLF branch.
   *
   * \tparam ConservativeComponentCount Number of conserved flux components.
   */
  template <std::size_t ConservativeComponentCount>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  gather_negative_split_stencil(
      const int component_index,
      const Real split_samples[2 * required_ghost_cells]
                              [ConservativeComponentCount],
      Real stencil_values[reconstruction_sample_count]) noexcept
  {
    for (int sample_index = 1;
         sample_index < 2 * required_ghost_cells;
         ++sample_index) {
      stencil_values[sample_index - 1] =
          split_samples[sample_index][component_index];
    }
  }

  /// Compute the three Jiang--Shu smoothness indicators.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  compute_smoothness_indicators(
      const Real stencil_values[reconstruction_sample_count],
      Real smoothness_indicators[candidate_count]) noexcept
  {
    smoothness_indicators[2] =
        Real(13.0 / 12.0) *
            square(stencil_values[4] - Real(2.0) * stencil_values[3] +
                   stencil_values[2]) +
        Real(0.25) *
            square(stencil_values[4] - Real(4.0) * stencil_values[3] +
                   Real(3.0) * stencil_values[2]);
    smoothness_indicators[1] =
        Real(13.0 / 12.0) *
            square(stencil_values[3] - Real(2.0) * stencil_values[2] +
                   stencil_values[1]) +
        Real(0.25) * square(stencil_values[3] - stencil_values[1]);
    smoothness_indicators[0] =
        Real(13.0 / 12.0) *
            square(stencil_values[2] - Real(2.0) * stencil_values[1] +
                   stencil_values[0]) +
        Real(0.25) *
            square(Real(3.0) * stencil_values[2] -
                   Real(4.0) * stencil_values[1] + stencil_values[0]);
  }

  /**
   * \brief Evaluate the three candidate polynomials with their factor of six.
   */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  evaluate_candidate_polynomials(
      const Real stencil_values[reconstruction_sample_count],
      Real scaled_candidate_values[candidate_count]) noexcept
  {
    scaled_candidate_values[2] =
        Real(11.0) * stencil_values[2] -
        Real(7.0) * stencil_values[3] +
        Real(2.0) * stencil_values[4];
    scaled_candidate_values[1] =
        -stencil_values[3] + Real(5.0) * stencil_values[2] +
        Real(2.0) * stencil_values[1];
    scaled_candidate_values[0] =
        Real(2.0) * stencil_values[2] +
        Real(5.0) * stencil_values[1] - stencil_values[0];
  }

  /// Evaluate one candidate without reading disabled stencil samples.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  evaluate_candidate_stencil(
      const int candidate_index,
      const Real stencil_values[reconstruction_sample_count],
      Real& smoothness_indicator,
      Real& scaled_candidate_value) noexcept
  {
    if (candidate_index == 2) {
      smoothness_indicator =
          Real(13.0 / 12.0) *
              square(stencil_values[4] - Real(2.0) * stencil_values[3] +
                     stencil_values[2]) +
          Real(0.25) *
              square(stencil_values[4] - Real(4.0) * stencil_values[3] +
                     Real(3.0) * stencil_values[2]);
      scaled_candidate_value =
          Real(11.0) * stencil_values[2] -
          Real(7.0) * stencil_values[3] +
          Real(2.0) * stencil_values[4];
    } else if (candidate_index == 1) {
      smoothness_indicator =
          Real(13.0 / 12.0) *
              square(stencil_values[3] - Real(2.0) * stencil_values[2] +
                     stencil_values[1]) +
          Real(0.25) * square(stencil_values[3] - stencil_values[1]);
      scaled_candidate_value =
          -stencil_values[3] + Real(5.0) * stencil_values[2] +
          Real(2.0) * stencil_values[1];
    } else {
      smoothness_indicator =
          Real(13.0 / 12.0) *
              square(stencil_values[2] - Real(2.0) * stencil_values[1] +
                     stencil_values[0]) +
          Real(0.25) *
              square(Real(3.0) * stencil_values[2] -
                     Real(4.0) * stencil_values[1] + stencil_values[0]);
      scaled_candidate_value =
          Real(2.0) * stencil_values[2] +
          Real(5.0) * stencil_values[1] - stencil_values[0];
    }
  }

  /// Freeze the WENO-Z weights formed from one complete five-point stencil.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static FrozenWeightSet
  freeze_full_stencil_weights(
      const Real stencil_values[reconstruction_sample_count],
      const Real epsilon = smoothness_epsilon()) noexcept
  {
    Real smoothness_indicators[candidate_count];
    compute_smoothness_indicators(stencil_values, smoothness_indicators);

    const Real global_smoothness_indicator = std::abs(
        smoothness_indicators[2] - smoothness_indicators[0]);
    const CandidateWeightArray linear_weights = optimal_linear_weights();

    FrozenWeightSet frozen_weights;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      frozen_weights.candidate_weights[candidate_index] =
          wenoz_alpha_factor(
              global_smoothness_indicator,
              smoothness_indicators[candidate_index], epsilon) *
          linear_weights[candidate_index];
      frozen_weights.weight_sum +=
          frozen_weights.candidate_weights[candidate_index];
    }

    if (!(frozen_weights.weight_sum >
          std::numeric_limits<Real>::min()) ||
        !std::isfinite(frozen_weights.weight_sum)) {
      // This is the existing high-order optimal-weight recovery.  It is not
      // a first-order flux fallback.
      frozen_weights.candidate_weights = linear_weights;
      frozen_weights.weight_sum = Real(10.0);
    }
    return frozen_weights;
  }

  /** Freeze the marker-permitted WENO-Z candidate weights at one face. */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static FrozenWeightSet
  freeze_available_stencil_weights(
      const Real stencil_values[reconstruction_sample_count],
      const CandidateWeightArray& available_linear_weights,
      const Real epsilon = smoothness_epsilon()) noexcept
  {
    bool candidate_available[candidate_count] = {};
    Real smoothness_indicators[candidate_count] = {};
    int available_candidate_count = 0;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      candidate_available[candidate_index] =
          available_linear_weights[candidate_index] > Real(0.0);
      if (!candidate_available[candidate_index]) {
        continue;
      }
      Real unused_scaled_candidate = Real(0.0);
      evaluate_candidate_stencil(
          candidate_index, stencil_values,
          smoothness_indicators[candidate_index],
          unused_scaled_candidate);
      ++available_candidate_count;
    }

    FrozenWeightSet frozen_weights;
    if (available_candidate_count == 0) {
      return frozen_weights;
    }

    Real global_smoothness_indicator = Real(0.0);
    if (candidate_available[0] && candidate_available[1] &&
        candidate_available[2]) {
      global_smoothness_indicator = std::abs(
          smoothness_indicators[2] - smoothness_indicators[0]);
    } else {
      Real minimum_smoothness_indicator =
          std::numeric_limits<Real>::max();
      Real maximum_smoothness_indicator = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        if (!candidate_available[candidate_index]) {
          continue;
        }
        minimum_smoothness_indicator = amrex::min(
            minimum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
        maximum_smoothness_indicator = amrex::max(
            maximum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
      }
      global_smoothness_indicator =
          available_candidate_count > 1
              ? maximum_smoothness_indicator - minimum_smoothness_indicator
              : Real(0.0);
    }

    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (!candidate_available[candidate_index]) {
        continue;
      }
      frozen_weights.candidate_weights[candidate_index] =
          wenoz_alpha_factor(
              global_smoothness_indicator,
              smoothness_indicators[candidate_index], epsilon) *
          available_linear_weights[candidate_index];
      frozen_weights.weight_sum +=
          frozen_weights.candidate_weights[candidate_index];
    }

    if (!(frozen_weights.weight_sum >
          std::numeric_limits<Real>::min()) ||
        !std::isfinite(frozen_weights.weight_sum)) {
      frozen_weights.candidate_weights = available_linear_weights;
      frozen_weights.weight_sum = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        frozen_weights.weight_sum +=
            available_linear_weights[candidate_index];
      }
    }
    return frozen_weights;
  }

  /// Apply previously frozen complete-stencil WENO-Z weights to new data.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_with_frozen_weights(
      const Real stencil_values[reconstruction_sample_count],
      const FrozenWeightSet& frozen_weights) noexcept
  {
    Real scaled_candidate_values[candidate_count];
    evaluate_candidate_polynomials(stencil_values, scaled_candidate_values);
    return (
        frozen_weights.candidate_weights[0] *
            scaled_candidate_values[0] +
        frozen_weights.candidate_weights[1] *
            scaled_candidate_values[1] +
        frozen_weights.candidate_weights[2] *
            scaled_candidate_values[2]) /
        (Real(6.0) * frozen_weights.weight_sum);
  }

  /// Reconstruct one interface value from the complete five-point stencil.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_interface_value(
      const Real stencil_values[reconstruction_sample_count],
      const Real epsilon = smoothness_epsilon()) noexcept
  {
    const FrozenWeightSet frozen_weights =
        freeze_full_stencil_weights(stencil_values, epsilon);
    return reconstruct_with_frozen_weights(
        stencil_values, frozen_weights);
  }

  /// Reconstruct from the candidate subset permitted by the IBM stencil plan.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_interface_value(
      const Real stencil_values[reconstruction_sample_count],
      const CandidateWeightArray& available_linear_weights,
      const Real epsilon = smoothness_epsilon()) noexcept
  {
    // candidate_available mirrors the nonzero weights in the stencil plan.
    bool candidate_available[candidate_count] = {};
    // smoothness_indicators stores only indicators safe to evaluate.
    Real smoothness_indicators[candidate_count] = {};
    // scaled_candidate_values stores only safe candidate polynomial values.
    Real scaled_candidate_values[candidate_count] = {};
    // available_candidate_count determines the reduced tau definition.
    int available_candidate_count = 0;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      candidate_available[candidate_index] =
          available_linear_weights[candidate_index] > Real(0.0);
      if (!candidate_available[candidate_index]) {
        continue;
      }
      evaluate_candidate_stencil(
          candidate_index, stencil_values,
          smoothness_indicators[candidate_index],
          scaled_candidate_values[candidate_index]);
      ++available_candidate_count;
    }

    if (available_candidate_count == 0) {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    // Use tau_5 only when the complete candidate set is available.
    Real global_smoothness_indicator = Real(0.0);
    if (candidate_available[0] && candidate_available[1] &&
        candidate_available[2]) {
      global_smoothness_indicator = std::abs(
          smoothness_indicators[2] - smoothness_indicators[0]);
    } else {
      // The available-indicator range avoids reading a disabled solid sample.
      Real minimum_smoothness_indicator =
          std::numeric_limits<Real>::max();
      Real maximum_smoothness_indicator = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        if (!candidate_available[candidate_index]) {
          continue;
        }
        minimum_smoothness_indicator = amrex::min(
            minimum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
        maximum_smoothness_indicator = amrex::max(
            maximum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
      }
      global_smoothness_indicator =
          available_candidate_count > 1
              ? maximum_smoothness_indicator - minimum_smoothness_indicator
              : Real(0.0);
    }

    // nonlinear_weights contains WENO-Z alpha_k for available candidates.
    Real nonlinear_weights[candidate_count] = {};
    // nonlinear_weight_sum normalizes the reduced set.
    Real nonlinear_weight_sum = Real(0.0);
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (!candidate_available[candidate_index]) {
        continue;
      }
      nonlinear_weights[candidate_index] =
          wenoz_alpha_factor(
              global_smoothness_indicator,
              smoothness_indicators[candidate_index],
              epsilon) *
          available_linear_weights[candidate_index];
      nonlinear_weight_sum += nonlinear_weights[candidate_index];
    }

    if (!(nonlinear_weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(nonlinear_weight_sum)) {
      return reconstruct_from_available_linear_candidates<candidate_count>(
          candidate_available, scaled_candidate_values,
          available_linear_weights);
    }

    // weighted_scaled_value accumulates the nonlinear candidate combination.
    Real weighted_scaled_value = Real(0.0);
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      weighted_scaled_value +=
          nonlinear_weights[candidate_index] *
          scaled_candidate_values[candidate_index];
    }
    return weighted_scaled_value /
           (Real(6.0) * nonlinear_weight_sum);
  }
};

/**
 * \brief Fifth-order CERISSE TENO reconstruction method.
 *
 * This method reuses the WENO-Z5 stencil geometry and candidate polynomials.
 * It replaces continuous WENO weighting with TENO candidate selection.  The
 * global smoothness indicator is the historical CERISSE expression rather
 * than the canonical tau_5 = abs(beta_2-beta_0) definition.
 */
struct Teno5 : public WenoZ5 {
  /// Frozen TENO selections can be reused by the R-Z pressure companion.
  static constexpr bool supports_rz_paired_pressure_flux = true;
  /// TENO selection uses its own machine-scale regularisation.
  static constexpr bool uses_weno_regularisation_epsilon = false;
  /// Fixed threshold used to remove a nonsmooth candidate.
  static constexpr Real selection_cutoff =
      Real(CERISSE_TENO5_SELECTION_CUTOFF);
  /// Method-manifest name exposed by the characteristic LLF driver.
  static constexpr const char* scheme_name =
      "characteristic_llf_teno5";

  /** Freeze the binary TENO5 candidate decision for companion fields. */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static FrozenWeightSet
  freeze_available_stencil_weights(
      const Real stencil_values[reconstruction_sample_count],
      const CandidateWeightArray& available_linear_weights,
      const Real /*epsilon*/ = smoothness_epsilon()) noexcept
  {
    bool candidate_available[candidate_count] = {};
    Real smoothness_indicators[candidate_count] = {};
    int available_candidate_count = 0;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      candidate_available[candidate_index] =
          available_linear_weights[candidate_index] > Real(0.0);
      if (!candidate_available[candidate_index]) {
        continue;
      }
      Real unused_scaled_candidate = Real(0.0);
      WenoZ5::evaluate_candidate_stencil(
          candidate_index, stencil_values,
          smoothness_indicators[candidate_index],
          unused_scaled_candidate);
      ++available_candidate_count;
    }

    FrozenWeightSet frozen_weights;
    if (available_candidate_count == 0) {
      return frozen_weights;
    }

    Real global_smoothness_indicator = Real(0.0);
    if (candidate_available[0] && candidate_available[1] &&
        candidate_available[2]) {
#if (CERISSE_TENO5_STANDARD_TAU5 == 1)
      global_smoothness_indicator = std::abs(
          smoothness_indicators[2] - smoothness_indicators[0]);
#else
      global_smoothness_indicator = std::abs(
          std::abs(
              smoothness_indicators[2] - smoothness_indicators[0]) -
          (smoothness_indicators[2] +
           Real(4.0) * smoothness_indicators[1] +
           smoothness_indicators[0]) /
              Real(6.0));
#endif
    } else {
      Real minimum_smoothness_indicator =
          std::numeric_limits<Real>::max();
      Real maximum_smoothness_indicator = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        if (!candidate_available[candidate_index]) {
          continue;
        }
        minimum_smoothness_indicator = amrex::min(
            minimum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
        maximum_smoothness_indicator = amrex::max(
            maximum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
      }
      global_smoothness_indicator =
          available_candidate_count > 1
              ? maximum_smoothness_indicator - minimum_smoothness_indicator
              : Real(0.0);
    }

    bool candidate_selected[candidate_count] = {};
    const bool selection_succeeded =
        select_teno_smooth_candidates<candidate_count>(
            candidate_available, smoothness_indicators,
            global_smoothness_indicator, selection_cutoff,
            candidate_selected);
    if (!selection_succeeded) {
      // Preserve the established TENO5 failure semantics: its nonlinear
      // selection falls back to reduced WENO-Z, not to linear TENO weights.
      return WenoZ5::freeze_available_stencil_weights(
          stencil_values, available_linear_weights);
    }
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (candidate_selected[candidate_index]) {
        frozen_weights.candidate_weights[candidate_index] =
            available_linear_weights[candidate_index];
      }
      frozen_weights.weight_sum +=
          frozen_weights.candidate_weights[candidate_index];
    }
    return frozen_weights;
  }

  /// Reconstruct one interface value from the complete five-point stencil.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_interface_value(
      const Real stencil_values[reconstruction_sample_count],
      const Real /*epsilon*/ = smoothness_epsilon()) noexcept
  {
    return reconstruct_interface_value(
        stencil_values, optimal_linear_weights());
  }

  /// Reconstruct from the candidate subset permitted by the stencil plan.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_interface_value(
      const Real stencil_values[reconstruction_sample_count],
      const CandidateWeightArray& available_linear_weights,
      const Real /*epsilon*/ = smoothness_epsilon()) noexcept
  {
    // candidate_available identifies geometrically permitted candidates.
    bool candidate_available[candidate_count] = {};
    // smoothness_indicators stores beta_k only for permitted candidates.
    Real smoothness_indicators[candidate_count] = {};
    // scaled_candidate_values stores six times each candidate value.
    Real scaled_candidate_values[candidate_count] = {};
    // available_candidate_count selects the full or reduced tau definition.
    int available_candidate_count = 0;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      candidate_available[candidate_index] =
          available_linear_weights[candidate_index] > Real(0.0);
      if (!candidate_available[candidate_index]) {
        continue;
      }
      WenoZ5::evaluate_candidate_stencil(
          candidate_index, stencil_values,
          smoothness_indicators[candidate_index],
          scaled_candidate_values[candidate_index]);
      ++available_candidate_count;
    }

    if (available_candidate_count == 0) {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    // global_smoothness_indicator supplies the TENO scale separation.
    Real global_smoothness_indicator = Real(0.0);
    if (candidate_available[0] && candidate_available[1] &&
        candidate_available[2]) {
#if (CERISSE_TENO5_STANDARD_TAU5 == 1)
      global_smoothness_indicator = std::abs(
          smoothness_indicators[2] - smoothness_indicators[0]);
#else
      global_smoothness_indicator = std::abs(
          std::abs(
              smoothness_indicators[2] - smoothness_indicators[0]) -
          (smoothness_indicators[2] +
           Real(4.0) * smoothness_indicators[1] +
           smoothness_indicators[0]) /
              Real(6.0));
#endif
    } else {
      // The available-indicator range avoids reading a disabled solid sample.
      Real minimum_smoothness_indicator =
          std::numeric_limits<Real>::max();
      Real maximum_smoothness_indicator = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        if (!candidate_available[candidate_index]) {
          continue;
        }
        minimum_smoothness_indicator = amrex::min(
            minimum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
        maximum_smoothness_indicator = amrex::max(
            maximum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
      }
      global_smoothness_indicator =
          available_candidate_count > 1
              ? maximum_smoothness_indicator - minimum_smoothness_indicator
              : Real(0.0);
    }

    // candidate_selected is the binary TENO decision for each candidate.
    bool candidate_selected[candidate_count] = {};
    if (!select_teno_smooth_candidates<candidate_count>(
            candidate_available, smoothness_indicators,
            global_smoothness_indicator, selection_cutoff,
            candidate_selected)) {
      return WenoZ5::reconstruct_interface_value(
          stencil_values, available_linear_weights);
    }

    // selected_linear_weight_sum normalizes the retained candidates.
    Real selected_linear_weight_sum = Real(0.0);
    // weighted_scaled_value accumulates weight_k times six times q_k.
    Real weighted_scaled_value = Real(0.0);
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (!candidate_selected[candidate_index]) {
        continue;
      }
      selected_linear_weight_sum +=
          available_linear_weights[candidate_index];
      weighted_scaled_value +=
          available_linear_weights[candidate_index] *
          scaled_candidate_values[candidate_index];
    }
    if (!(selected_linear_weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(
            weighted_scaled_value + selected_linear_weight_sum)) {
      return WenoZ5::reconstruct_interface_value(
          stencil_values, available_linear_weights);
    }
    return weighted_scaled_value /
           (Real(6.0) * selected_linear_weight_sum);
  }
};

/**
 * \brief Sixth-order TENO reconstruction method on a six-point stencil.
 *
 * The three-point smoothness indicators use the same samples as their
 * candidate polynomials.  The fourth candidate uses four points.  This is
 * the corrected alignment relative to the historical CERISSE implementation.
 */
struct Teno6 {
  /// Frozen TENO6 selections can be reused by the R-Z pressure companion.
  static constexpr bool supports_rz_paired_pressure_flux = true;
  /// Number of cell layers required on each side of a face.
  static constexpr int required_ghost_cells = 3;
  /// Number of ordered samples used by each split-flux reconstruction.
  static constexpr int reconstruction_sample_count = 6;
  /// Three three-point candidates plus one four-point candidate.
  static constexpr int candidate_count = 4;
  /// TENO selection does not use the WENO-Z epsilon argument.
  static constexpr bool uses_weno_regularisation_epsilon = false;
  /// Fixed threshold used to remove a nonsmooth candidate.
  static constexpr Real selection_cutoff = Real(1.0e-5);
  /// Method-manifest name exposed by the characteristic LLF driver.
  static constexpr const char* scheme_name =
      "characteristic_llf_teno6";

  /// Fixed-size storage for the four candidate linear weights.
  using CandidateWeightArray = amrex::GpuArray<Real, candidate_count>;

  /// Binary-selected linear weights frozen from one TENO6 stencil.
  struct FrozenWeightSet {
    CandidateWeightArray candidate_weights{};
    Real weight_sum = Real(0.0);
  };

  /// Return the unnormalized TENO6 linear weights 6:9:1:4.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr CandidateWeightArray
  optimal_linear_weights() noexcept
  {
    /*
     * Candidate order is downwind three-point, centred three-point,
     * upwind three-point, and four-point downwind.  The normalized weights
     * are 6:9:1:4.  Their smooth combination gives the sixth-order central
     * finite-difference flux [1,-8,37,37,-8,1]/60.
     */
    return {Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
  }

  /// Report whether a positive-branch candidate uses a driver sample.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
  positive_candidate_uses_sample(
      const int candidate_index,
      const int sample_index) noexcept
  {
    if (candidate_index == 0) {
      return sample_index >= 2 && sample_index <= 4;
    }
    if (candidate_index == 1) {
      return sample_index >= 1 && sample_index <= 3;
    }
    if (candidate_index == 2) {
      return sample_index >= 0 && sample_index <= 2;
    }
    return sample_index >= 2 && sample_index <= 5;
  }

  /// Report whether a negative-branch candidate uses a driver sample.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
  negative_candidate_uses_sample(
      const int candidate_index,
      const int sample_index) noexcept
  {
    if (candidate_index == 0) {
      return sample_index >= 1 && sample_index <= 3;
    }
    if (candidate_index == 1) {
      return sample_index >= 2 && sample_index <= 4;
    }
    if (candidate_index == 2) {
      return sample_index >= 3 && sample_index <= 5;
    }
    return sample_index >= 0 && sample_index <= 3;
  }

  /**
   * \brief Gather the six ordered samples for the positive LLF branch.
   *
   * \tparam ConservativeComponentCount Number of conserved flux components.
   */
  template <std::size_t ConservativeComponentCount>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  gather_positive_split_stencil(
      const int component_index,
      const Real split_samples[2 * required_ghost_cells]
                              [ConservativeComponentCount],
      Real stencil_values[reconstruction_sample_count]) noexcept
  {
    for (int sample_index = 0;
         sample_index < 2 * required_ghost_cells;
         ++sample_index) {
      stencil_values[sample_index] =
          split_samples[2 * required_ghost_cells - 1 - sample_index]
                       [component_index];
    }
  }

  /**
   * \brief Gather the six ordered samples for the negative LLF branch.
   *
   * \tparam ConservativeComponentCount Number of conserved flux components.
   */
  template <std::size_t ConservativeComponentCount>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  gather_negative_split_stencil(
      const int component_index,
      const Real split_samples[2 * required_ghost_cells]
                              [ConservativeComponentCount],
      Real stencil_values[reconstruction_sample_count]) noexcept
  {
    for (int sample_index = 0;
         sample_index < 2 * required_ghost_cells;
         ++sample_index) {
      stencil_values[sample_index] =
          split_samples[sample_index][component_index];
    }
  }

  /// Evaluate one candidate smoothness indicator and polynomial value.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
  evaluate_candidate_stencil(
      const int candidate_index,
      const Real stencil_values[reconstruction_sample_count],
      Real& smoothness_indicator,
      Real& scaled_candidate_value) noexcept
  {
    /*
     * Each smoothness indicator uses exactly the support of its candidate.
     * This is the Fu--Hu--Adams TENO6 stencil alignment.  The historical
     * CERISSE Teno6 code shifted the three three-point indicators by one
     * sample relative to their candidate polynomials.
     */
    if (candidate_index == 0) {
      smoothness_indicator =
          Real(13.0 / 12.0) *
              square(stencil_values[1] - Real(2.0) * stencil_values[2] +
                     stencil_values[3]) +
          Real(0.25) *
              square(stencil_values[1] - Real(4.0) * stencil_values[2] +
                     Real(3.0) * stencil_values[3]);
      scaled_candidate_value =
          Real(2.0) * stencil_values[3] +
          Real(5.0) * stencil_values[2] - stencil_values[1];
    } else if (candidate_index == 1) {
      smoothness_indicator =
          Real(13.0 / 12.0) *
              square(stencil_values[2] - Real(2.0) * stencil_values[3] +
                     stencil_values[4]) +
          Real(0.25) * square(stencil_values[2] - stencil_values[4]);
      scaled_candidate_value =
          -stencil_values[4] + Real(5.0) * stencil_values[3] +
          Real(2.0) * stencil_values[2];
    } else if (candidate_index == 2) {
      smoothness_indicator =
          Real(13.0 / 12.0) *
              square(stencil_values[3] - Real(2.0) * stencil_values[4] +
                     stencil_values[5]) +
          Real(0.25) *
              square(Real(3.0) * stencil_values[3] -
                     Real(4.0) * stencil_values[4] + stencil_values[5]);
      scaled_candidate_value =
          Real(11.0) * stencil_values[3] -
          Real(7.0) * stencil_values[4] +
          Real(2.0) * stencil_values[5];
    } else {
      smoothness_indicator =
          Real(1.0 / 36.0) *
              square(
                  -Real(11.0) * stencil_values[3] +
                  Real(18.0) * stencil_values[2] -
                  Real(9.0) * stencil_values[1] +
                  Real(2.0) * stencil_values[0]) +
          Real(13.0 / 12.0) *
              square(
                  Real(2.0) * stencil_values[3] -
                  Real(5.0) * stencil_values[2] +
                  Real(4.0) * stencil_values[1] - stencil_values[0]) +
          Real(781.0 / 720.0) *
              square(
                  -stencil_values[3] + Real(3.0) * stencil_values[2] -
                  Real(3.0) * stencil_values[1] + stencil_values[0]);
      scaled_candidate_value =
          Real(0.5) *
          (stencil_values[0] - Real(5.0) * stencil_values[1] +
           Real(13.0) * stencil_values[2] +
           Real(3.0) * stencil_values[3]);
    }
  }

  /** Freeze the binary TENO6 candidate decision for companion fields. */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static FrozenWeightSet
  freeze_available_stencil_weights(
      const Real stencil_values[reconstruction_sample_count],
      const CandidateWeightArray& available_linear_weights,
      const Real /*epsilon*/ = smoothness_epsilon()) noexcept
  {
    bool candidate_available[candidate_count] = {};
    Real smoothness_indicators[candidate_count] = {};
    int available_candidate_count = 0;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      candidate_available[candidate_index] =
          available_linear_weights[candidate_index] > Real(0.0);
      if (!candidate_available[candidate_index]) {
        continue;
      }
      Real unused_scaled_candidate = Real(0.0);
      evaluate_candidate_stencil(
          candidate_index, stencil_values,
          smoothness_indicators[candidate_index],
          unused_scaled_candidate);
      ++available_candidate_count;
    }

    FrozenWeightSet frozen_weights;
    if (available_candidate_count == 0) {
      return frozen_weights;
    }

    Real global_smoothness_indicator = Real(0.0);
    if (candidate_available[0] && candidate_available[1] &&
        candidate_available[2] && candidate_available[3]) {
      global_smoothness_indicator = std::abs(
          smoothness_indicators[3] -
          (smoothness_indicators[2] +
           Real(4.0) * smoothness_indicators[1] +
           smoothness_indicators[0]) /
              Real(6.0));
    } else {
      Real minimum_smoothness_indicator =
          std::numeric_limits<Real>::max();
      Real maximum_smoothness_indicator = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        if (!candidate_available[candidate_index]) {
          continue;
        }
        minimum_smoothness_indicator = amrex::min(
            minimum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
        maximum_smoothness_indicator = amrex::max(
            maximum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
      }
      global_smoothness_indicator =
          available_candidate_count > 1
              ? maximum_smoothness_indicator - minimum_smoothness_indicator
              : Real(0.0);
    }

    bool candidate_selected[candidate_count] = {};
    const bool selection_succeeded =
        select_teno_smooth_candidates<candidate_count>(
            candidate_available, smoothness_indicators,
            global_smoothness_indicator, selection_cutoff,
            candidate_selected);
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (selection_succeeded && candidate_selected[candidate_index]) {
        frozen_weights.candidate_weights[candidate_index] =
            available_linear_weights[candidate_index];
      } else if (!selection_succeeded &&
                 candidate_available[candidate_index]) {
        frozen_weights.candidate_weights[candidate_index] =
            available_linear_weights[candidate_index];
      }
      frozen_weights.weight_sum +=
          frozen_weights.candidate_weights[candidate_index];
    }
    return frozen_weights;
  }

  /// Apply one frozen TENO6 candidate set to a companion stencil.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_with_frozen_weights(
      const Real stencil_values[reconstruction_sample_count],
      const FrozenWeightSet& frozen_weights) noexcept
  {
    Real weighted_scaled_value = Real(0.0);
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (!(frozen_weights.candidate_weights[candidate_index] >
            Real(0.0))) {
        continue;
      }
      Real unused_smoothness = Real(0.0);
      Real scaled_candidate_value = Real(0.0);
      evaluate_candidate_stencil(
          candidate_index, stencil_values,
          unused_smoothness, scaled_candidate_value);
      weighted_scaled_value +=
          frozen_weights.candidate_weights[candidate_index] *
          scaled_candidate_value;
    }
    return weighted_scaled_value /
           (Real(6.0) * frozen_weights.weight_sum);
  }

  /// Reconstruct one interface value from the complete six-point stencil.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_interface_value(
      const Real stencil_values[reconstruction_sample_count],
      const Real /*epsilon*/ = smoothness_epsilon()) noexcept
  {
    return reconstruct_interface_value(
        stencil_values, optimal_linear_weights());
  }

  /// Reconstruct from the candidate subset permitted by the stencil plan.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static Real
  reconstruct_interface_value(
      const Real stencil_values[reconstruction_sample_count],
      const CandidateWeightArray& available_linear_weights,
      const Real /*epsilon*/ = smoothness_epsilon()) noexcept
  {
    // candidate_available identifies geometrically permitted candidates.
    bool candidate_available[candidate_count] = {};
    // smoothness_indicators stores beta_k only for permitted candidates.
    Real smoothness_indicators[candidate_count] = {};
    // scaled_candidate_values stores six times each candidate value.
    Real scaled_candidate_values[candidate_count] = {};
    // available_candidate_count selects the full or reduced tau definition.
    int available_candidate_count = 0;
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      candidate_available[candidate_index] =
          available_linear_weights[candidate_index] > Real(0.0);
      if (!candidate_available[candidate_index]) {
        continue;
      }
      evaluate_candidate_stencil(
          candidate_index, stencil_values,
          smoothness_indicators[candidate_index],
          scaled_candidate_values[candidate_index]);
      ++available_candidate_count;
    }

    if (available_candidate_count == 0) {
      return std::numeric_limits<Real>::quiet_NaN();
    }

    // global_smoothness_indicator is the full TENO6 reference tau.
    Real global_smoothness_indicator = Real(0.0);
    if (candidate_available[0] && candidate_available[1] &&
        candidate_available[2] && candidate_available[3]) {
      global_smoothness_indicator = std::abs(
          smoothness_indicators[3] -
          (smoothness_indicators[2] +
           Real(4.0) * smoothness_indicators[1] +
           smoothness_indicators[0]) /
              Real(6.0));
    } else {
      /*
       * The full TENO6 reference indicator requires all four candidates.
       * Near a shared GP, use the range of the available indicators only for
       * ENO-like selection.  It never reads a disabled solid-side candidate.
       */
      Real minimum_smoothness_indicator =
          std::numeric_limits<Real>::max();
      Real maximum_smoothness_indicator = Real(0.0);
      for (int candidate_index = 0;
           candidate_index < candidate_count;
           ++candidate_index) {
        if (!candidate_available[candidate_index]) {
          continue;
        }
        minimum_smoothness_indicator = amrex::min(
            minimum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
        maximum_smoothness_indicator = amrex::max(
            maximum_smoothness_indicator,
            smoothness_indicators[candidate_index]);
      }
      global_smoothness_indicator =
          available_candidate_count > 1
              ? maximum_smoothness_indicator - minimum_smoothness_indicator
              : Real(0.0);
    }

    // candidate_selected is the binary TENO decision for each candidate.
    bool candidate_selected[candidate_count] = {};
    if (!select_teno_smooth_candidates<candidate_count>(
            candidate_available, smoothness_indicators,
            global_smoothness_indicator, selection_cutoff,
            candidate_selected)) {
      return reconstruct_from_available_linear_candidates<candidate_count>(
          candidate_available, scaled_candidate_values,
          available_linear_weights);
    }

    // selected_linear_weight_sum normalizes the retained candidates.
    Real selected_linear_weight_sum = Real(0.0);
    // weighted_scaled_value accumulates weight_k times six times q_k.
    Real weighted_scaled_value = Real(0.0);
    for (int candidate_index = 0;
         candidate_index < candidate_count;
         ++candidate_index) {
      if (!candidate_selected[candidate_index]) {
        continue;
      }
      selected_linear_weight_sum +=
          available_linear_weights[candidate_index];
      weighted_scaled_value +=
          available_linear_weights[candidate_index] *
          scaled_candidate_values[candidate_index];
    }
    if (!(selected_linear_weight_sum > std::numeric_limits<Real>::min()) ||
        !std::isfinite(
            weighted_scaled_value + selected_linear_weight_sum)) {
      return reconstruct_from_available_linear_candidates<candidate_count>(
          candidate_available, scaled_candidate_values,
          available_linear_weights);
    }
    return weighted_scaled_value /
           (Real(6.0) * selected_linear_weight_sum);
  }
};

/**
 * \brief Candidate weights and sample requirements for one Cartesian face.
 *
 * A plan disables only candidate stencils that require an unavailable deep
 * solid cell.  Fluid cells and reconstructed ghost-point cells remain valid.
 * No cut-cell geometry, area fraction, or volume fraction enters this plan.
 */
template <typename ReconstructionMethod>
struct FaceReconstructionPlan {
  /// Candidate-weight storage supplied by the selected reconstruction method.
  using CandidateWeightArray =
      typename ReconstructionMethod::CandidateWeightArray;

  /// Linear weights retained for the positive LLF branch.
  CandidateWeightArray positive_linear_weights{};
  /// Linear weights retained for the negative LLF branch.
  CandidateWeightArray negative_linear_weights{};
  /// Samples that must be read by at least one retained candidate.
  bool sample_state_required[
      2 * ReconstructionMethod::required_ghost_cells] = {};
  /// True when both LLF branches retain at least one high-order candidate.
  bool high_order_reconstruction_available = false;
  /// True when no candidate has been disabled on either LLF branch.
  bool full_stencil_available = false;
};

/**
 * \brief Build the reconstruction plan from the state availability flags.
 *
 * \tparam ReconstructionMethod Selected WENO or TENO reconstruction method.
 * \param sample_state_available True for each stencil sample that can be read.
 * \return Candidate weights and required samples for both LLF branches.
 */
template <typename ReconstructionMethod>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
FaceReconstructionPlan<ReconstructionMethod>
build_face_reconstruction_plan(
    const bool sample_state_available[
        2 * ReconstructionMethod::required_ghost_cells]) noexcept
{
  // plan is the complete reconstruction decision returned to the face loop.
  FaceReconstructionPlan<ReconstructionMethod> plan;
  plan.positive_linear_weights =
      ReconstructionMethod::optimal_linear_weights();
  plan.negative_linear_weights =
      ReconstructionMethod::optimal_linear_weights();

  // These flags distinguish partial high-order support from a full stencil.
  bool has_positive_candidate = false;
  bool has_negative_candidate = false;
  bool all_positive_candidates_available = true;
  bool all_negative_candidates_available = true;

  for (int candidate_index = 0;
       candidate_index < ReconstructionMethod::candidate_count;
       ++candidate_index) {
    // A candidate remains available only when all of its samples have states.
    bool positive_candidate_available = true;
    bool negative_candidate_available = true;
    for (int sample_index = 0;
         sample_index <
             2 * ReconstructionMethod::required_ghost_cells;
         ++sample_index) {
      if (ReconstructionMethod::positive_candidate_uses_sample(
              candidate_index, sample_index) &&
          !sample_state_available[sample_index]) {
        positive_candidate_available = false;
      }
      if (ReconstructionMethod::negative_candidate_uses_sample(
              candidate_index, sample_index) &&
          !sample_state_available[sample_index]) {
        negative_candidate_available = false;
      }
    }

    if (!positive_candidate_available) {
      plan.positive_linear_weights[candidate_index] = Real(0.0);
      all_positive_candidates_available = false;
    } else {
      has_positive_candidate = true;
      for (int sample_index = 0;
           sample_index <
               2 * ReconstructionMethod::required_ghost_cells;
           ++sample_index) {
        if (ReconstructionMethod::positive_candidate_uses_sample(
                candidate_index, sample_index)) {
          plan.sample_state_required[sample_index] = true;
        }
      }
    }

    if (!negative_candidate_available) {
      plan.negative_linear_weights[candidate_index] = Real(0.0);
      all_negative_candidates_available = false;
    } else {
      has_negative_candidate = true;
      for (int sample_index = 0;
           sample_index <
               2 * ReconstructionMethod::required_ghost_cells;
           ++sample_index) {
        if (ReconstructionMethod::negative_candidate_uses_sample(
                candidate_index, sample_index)) {
          plan.sample_state_required[sample_index] = true;
        }
      }
    }
  }

  plan.high_order_reconstruction_available =
      has_positive_candidate && has_negative_candidate;
  plan.full_stencil_available =
      all_positive_candidates_available &&
      all_negative_candidates_available;
  return plan;
}

/// Marker component 0. Zero denotes fluid and a nonzero value denotes solid.
inline constexpr int solid_body_marker_component = 0;
/// Marker component 1. A nonzero value identifies a ghost-point cell.
inline constexpr int ghost_point_marker_component = 1;

/// Return true when the cell lies in the active fluid region.
template <typename IbmMarkerView>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_fluid_cell(
    const amrex::IntVect& cell_index,
    const IbmMarkerView& ibm_markers) noexcept
{
  return ibm_markers(cell_index, solid_body_marker_component) ==
         std::uint8_t(0);
}

/// Return true when the cell lies inside an immersed solid geometry.
template <typename IbmMarkerView>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_solid_cell(
    const amrex::IntVect& cell_index,
    const IbmMarkerView& ibm_markers) noexcept
{
  return !is_fluid_cell(cell_index, ibm_markers);
}

/// Return true when a solid cell owns a ghost-point boundary state.
template <typename IbmMarkerView>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_ghost_point_cell(
    const amrex::IntVect& cell_index,
    const IbmMarkerView& ibm_markers) noexcept
{
  return ibm_markers(cell_index, ghost_point_marker_component) !=
         std::uint8_t(0);
}

/// Return true when the flux stencil may read a state from this cell.
template <typename IbmMarkerView>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_flux_stencil_state_available(
    const amrex::IntVect& cell_index,
    const IbmMarkerView& ibm_markers) noexcept
{
  return is_fluid_cell(cell_index, ibm_markers) ||
         is_ghost_point_cell(cell_index, ibm_markers);
}

/// Check density, pressure, and sound speed next to a face.
template <typename PrimitiveStateView, typename Closure>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
face_adjacent_thermodynamic_values_are_valid(
    const amrex::IntVect& face_index,
    const amrex::IntVect& direction_offset,
    const PrimitiveStateView& primitive_states) noexcept
{
  // left_cell_index is the cell immediately below the face index.
  const amrex::IntVect left_cell_index = face_index - direction_offset;
  // finite_thermodynamic_sum checks all six scalar values at once.
  const Real finite_thermodynamic_sum =
      primitive_states(left_cell_index, Closure::QRHO) +
      primitive_states(face_index, Closure::QRHO) +
      primitive_states(left_cell_index, Closure::QPRES) +
      primitive_states(face_index, Closure::QPRES) +
      primitive_states(left_cell_index, Closure::QC) +
      primitive_states(face_index, Closure::QC);
  return primitive_states(left_cell_index, Closure::QRHO) > Real(0.0) &&
         primitive_states(face_index, Closure::QRHO) > Real(0.0) &&
         primitive_states(left_cell_index, Closure::QPRES) > Real(0.0) &&
         primitive_states(face_index, Closure::QPRES) > Real(0.0) &&
         primitive_states(left_cell_index, Closure::QC) > Real(0.0) &&
         primitive_states(face_index, Closure::QC) > Real(0.0) &&
         std::isfinite(finite_thermodynamic_sum);
}

/// Set every conservative component of an unused solid face to zero.
template <typename FaceFluxView, typename Closure>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
set_face_flux_to_zero(
    const amrex::IntVect& face_index,
    const FaceFluxView& face_flux) noexcept
{
  for (int component_index = 0;
       component_index < Closure::NCONS;
       ++component_index) {
    face_flux(face_index, component_index) = Real(0.0);
  }
}

/// Fill a failed face flux with NaNs so the invalid state cannot be hidden.
template <typename FaceFluxView, typename Closure>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
mark_face_flux_invalid(
    const amrex::IntVect& face_index,
    const FaceFluxView& face_flux) noexcept
{
  // invalid_value deliberately propagates through the subsequent RHS update.
  const Real invalid_value = std::numeric_limits<Real>::quiet_NaN();
  for (int component_index = 0;
       component_index < Closure::NCONS;
       ++component_index) {
    face_flux(face_index, component_index) = invalid_value;
  }
}

}  // namespace Reconstruction

// Existing prob.h files use ReconScheme; retain it as a public namespace alias.
namespace ReconScheme = Reconstruction;

/**
 * \brief Conservative finite-difference characteristic local LLF flux driver.
 *
 * \tparam ReconstructionMethod WENO or TENO interface reconstruction method.
 * \tparam Closure Thermodynamic conversion and characteristic operations.
 */
template <typename ReconstructionMethod, typename Closure>
class characteristic_llf_flux_t {
 public:
  /// Select whether face-stencil availability comes from IBM markers.
  /// Public: nvcc's generated kernel stubs reference this template argument
  /// from namespace scope, which fails if the nested type is private.
  enum class FaceStencilMode {
    all_fluid,          ///< Every stencil sample is assumed to have a state.
    ghost_point_aware  ///< Fluid and reconstructed GP states may be read.
  };

  // Optional, default-off diagnostic output.  Component zero is one when the
  // face used the two-state first-order LLF branch and zero otherwise.
  using DirectionalFaceFallbackMaskViews =
      std::array<amrex::Array4<int>, AMREX_SPACEDIM>;
  static constexpr bool local_llf_fallback_mask_capable = true;

 private:

  /**
   * \brief Runtime controls for the shared face-local LLF fallback.
   *
   * The same controls are used by WENO-Z5, TENO5, and TENO6 in Cartesian
   * and RZ coordinates, with either all-fluid or shared-GP stencils.
   */
  struct LocalLlfFallbackOptions {
    int enabled;
    amrex::Real density_jump_threshold;
    amrex::Real density_curvature_threshold;
    /// Curvature threshold used only on the axial faces of the axis-adjacent
    /// cell row in RZ.  The discrete RZ divergence amplifies the outer radial
    /// flux there by 2/dr, so the protection has to tighten as dr shrinks
    /// while every other face keeps density_curvature_threshold.
    amrex::Real axis_row_curvature_threshold;
  };

  /// Read and validate the face-local fallback controls once per executable.
  static const LocalLlfFallbackOptions& local_llf_fallback_options()
  {
    static const LocalLlfFallbackOptions options = [] {
      LocalLlfFallbackOptions value{
          1, amrex::Real(3.0e-3), amrex::Real(8.0e-2), amrex::Real(1.0e-2)};
      amrex::ParmParse parameters("cns");
      parameters.query("llf_first_order_fallback", value.enabled);
      parameters.query(
          "llf_fallback_density_jump", value.density_jump_threshold);
      parameters.query(
          "llf_fallback_density_curvature",
          value.density_curvature_threshold);
      parameters.query(
          "llf_fallback_axis_row_curvature",
          value.axis_row_curvature_threshold);

      if (value.enabled != 0 && value.enabled != 1) {
        amrex::Abort("cns.llf_first_order_fallback must be 0 or 1");
      }
      if (!std::isfinite(value.density_jump_threshold) ||
          value.density_jump_threshold < amrex::Real(0.0)) {
        amrex::Abort(
            "cns.llf_fallback_density_jump must be finite and non-negative");
      }
      if (!std::isfinite(value.density_curvature_threshold) ||
          value.density_curvature_threshold < amrex::Real(0.0) ||
          value.density_curvature_threshold > amrex::Real(1.0)) {
        amrex::Abort(
            "cns.llf_fallback_density_curvature must be in [0,1]");
      }
      if (!std::isfinite(value.axis_row_curvature_threshold) ||
          value.axis_row_curvature_threshold < amrex::Real(0.0) ||
          value.axis_row_curvature_threshold > amrex::Real(1.0)) {
        amrex::Abort(
            "cns.llf_fallback_axis_row_curvature must be in [0,1]");
      }
      return value;
    }();
    return options;
  }

  /// One writable AMReX face-flux array for each spatial direction.
  using DirectionalFaceFluxArrays =
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM>;

  /// Per-face plan for the candidates and samples that may be evaluated.
  using ReconstructionPlan =
      Reconstruction::FaceReconstructionPlan<ReconstructionMethod>;

  /// Number of samples required on either side of a Cartesian face.
  static constexpr int required_ghost_cells =
      ReconstructionMethod::required_ghost_cells;
  static_assert(required_ghost_cells <= Closure::NGHOST);
  static_assert(
      ReconstructionMethod::reconstruction_sample_count <=
      2 * required_ghost_cells);

 public:
  /// Name printed in the runtime numerical-method manifest.
  static constexpr const char* scheme_name =
      ReconstructionMethod::scheme_name;

  /// Every characteristic LLF selector exposes a frozen candidate decision
  /// for the paired R-Z radial-pressure operator.
  static constexpr bool rz_paired_pressure_flux_capable =
      ReconstructionMethod::supports_rz_paired_pressure_flux;

  /// Characteristic LLF WENO/TENO must never silently fall back to the legacy
  /// direct-p/r R-Z assembly.  rhs_dt propagates this independently of the
  /// diffusive operator so an incompatible viscous/LES composition fails
  /// closed instead of changing the inviscid discretization.
  static constexpr bool rz_paired_pressure_flux_required = true;

  /// This driver does not perform RZ annular pressure recovery internally.
  static constexpr bool rz_annular_pressure_consistency_capable = false;
  /// The matching Euler geometric source remains active in RZ coordinates.
  static constexpr bool rz_euler_geometric_source_active = true;
  /// The radial-axis face stores the reconstructed metric h-flux directly.
  static constexpr bool rz_radial_axis_face_flux_is_metric = true;

  /// Print the resolved fallback configuration in the method manifest.
  static void print_local_llf_fallback_manifest()
  {
    const LocalLlfFallbackOptions& options =
        local_llf_fallback_options();
    amrex::Print()
        << "[Numerics] local_first_order_llf_fallback="
        << (options.enabled != 0 ? "enabled" : "disabled")
        << " density_jump=" << options.density_jump_threshold
        << " density_curvature=" << options.density_curvature_threshold
        << " axis_row_curvature=" << options.axis_row_curvature_threshold
        << " scope=cartesian_rz_all_fluid_shared_gp"
        << " first_interior_rz_gate=density_jump_only"
        << " rz_axis_row_axial_faces=tightened_curvature"
        << " axis_face=excluded\n";
  }

  /// Construct a stateless characteristic LLF flux operator.
  AMREX_GPU_HOST_DEVICE
  characteristic_llf_flux_t() = default;

  /**
   * \brief Compute all-fluid inviscid face fluxes for one AMReX grid box.
   *
   * The method name is retained because rhs_dt calls every inviscid operator
   * through the common eflux interface.
   */
  void eflux(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure)
  {
    // unused_ibm_markers is never read by the all-fluid compile-time branch.
    amrex::Array4<std::uint8_t> unused_ibm_markers;
    amrex::Array4<amrex::Real> unused_paired_pressure_flux;
    DirectionalFaceFallbackMaskViews unused_fallback_mask;
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::all_fluid, false, false>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        primitive_states, Closure::QPRES,
        unused_ibm_markers, unused_paired_pressure_flux,
        unused_fallback_mask);
  }

  /// All-fluid Cartesian flux with an exact face-local fallback mask.
  void eflux_with_local_llf_fallback_mask(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const DirectionalFaceFallbackMaskViews& fallback_mask)
  {
    amrex::Array4<std::uint8_t> unused_ibm_markers;
    amrex::Array4<amrex::Real> unused_paired_pressure_flux;
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::all_fluid, false, true>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        primitive_states, Closure::QPRES,
        unused_ibm_markers, unused_paired_pressure_flux, fallback_mask);
  }

  /**
   * \brief Compute all-fluid RZ fluxes with a paired pressure face flux.
   *
   * The pressure reconstruction view contains stage-local radial-centre
   * pressures.  It aliases primitive pressure for point-state runs and may
   * name a separately recovered centre-pressure field for all-fluid annular
   * cell averages.
   */
  void eflux_with_rz_paired_pressure(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<const amrex::Real>&
          pressure_reconstruction_states,
      const int pressure_reconstruction_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux)
  {
    static_assert(
        ReconstructionMethod::supports_rz_paired_pressure_flux,
        "Paired RZ pressure flux requires a frozen LLF reconstruction");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geometry.IsRZ(),
        "Paired RZ pressure flux requires cylindrical RZ geometry");
    amrex::Array4<std::uint8_t> unused_ibm_markers;
    DirectionalFaceFallbackMaskViews unused_fallback_mask;
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::all_fluid, true, false>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        pressure_reconstruction_states, pressure_reconstruction_component,
        unused_ibm_markers, radial_pressure_face_flux,
        unused_fallback_mask);
  }

  /// All-fluid paired R-Z flux with an exact face-local fallback mask.
  void eflux_with_rz_paired_pressure_and_local_llf_fallback_mask(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<const amrex::Real>& pressure_reconstruction_states,
      const int pressure_reconstruction_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux,
      const DirectionalFaceFallbackMaskViews& fallback_mask)
  {
    static_assert(
        ReconstructionMethod::supports_rz_paired_pressure_flux,
        "Paired RZ pressure flux requires a frozen LLF reconstruction");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geometry.IsRZ(),
        "Paired RZ pressure flux requires cylindrical RZ geometry");
    amrex::Array4<std::uint8_t> unused_ibm_markers;
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::all_fluid, true, true>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        pressure_reconstruction_states, pressure_reconstruction_component,
        unused_ibm_markers, radial_pressure_face_flux, fallback_mask);
  }

  /**
   * \brief Compute pure ghost-point IBM inviscid face fluxes for one grid box.
   *
   * Every reconstruction method using this driver shares the same marker-aware
   * ghost-point stencil planning and characteristic LLF flux construction.
   */
  void eflux_ibm(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<std::uint8_t>& ibm_markers)
  {
#ifdef AMREX_USE_GPIBM
    amrex::Array4<amrex::Real> unused_paired_pressure_flux;
    DirectionalFaceFallbackMaskViews unused_fallback_mask;
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::ghost_point_aware, false, false>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        primitive_states, Closure::QPRES,
        ibm_markers, unused_paired_pressure_flux, unused_fallback_mask);
#else
    amrex::ignore_unused(
        geometry, grid_iterator, primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        ibm_markers);
    amrex::Abort(
        "weno_t::eflux_ibm requires the pure shared-GP IBM build");
#endif
  }

  /// Pure shared-GP Cartesian flux with an exact fallback mask.
  void eflux_ibm_with_local_llf_fallback_mask(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<std::uint8_t>& ibm_markers,
      const DirectionalFaceFallbackMaskViews& fallback_mask)
  {
#ifdef AMREX_USE_GPIBM
    amrex::Array4<amrex::Real> unused_paired_pressure_flux;
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::ghost_point_aware, false, true>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        primitive_states, Closure::QPRES,
        ibm_markers, unused_paired_pressure_flux, fallback_mask);
#else
    amrex::ignore_unused(
        geometry, grid_iterator, primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        ibm_markers, fallback_mask);
    amrex::Abort(
        "WENO fallback-mask capture requires the pure shared-GP build");
#endif
  }

  /** Compute pure shared-GP R-Z fluxes with the paired pressure operator. */
  void eflux_ibm_with_rz_paired_pressure(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<std::uint8_t>& ibm_markers,
      const amrex::Array4<const amrex::Real>&
          pressure_reconstruction_states,
      const int pressure_reconstruction_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux)
  {
#ifdef AMREX_USE_GPIBM
    static_assert(
        ReconstructionMethod::supports_rz_paired_pressure_flux,
        "GP paired RZ pressure requires a frozen LLF reconstruction");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geometry.IsRZ(),
        "GP paired RZ pressure requires cylindrical RZ geometry");
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::ghost_point_aware, true, false>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        pressure_reconstruction_states, pressure_reconstruction_component,
        ibm_markers, radial_pressure_face_flux,
        DirectionalFaceFallbackMaskViews{});
#else
    amrex::ignore_unused(
        geometry, grid_iterator, primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        ibm_markers, pressure_reconstruction_states,
        pressure_reconstruction_component, radial_pressure_face_flux);
    amrex::Abort(
        "weno_t::eflux_ibm_with_rz_paired_pressure requires pure shared-GP");
#endif
  }

  /// Pure shared-GP paired R-Z flux with an exact fallback mask.
  void eflux_ibm_with_rz_paired_pressure_and_local_llf_fallback_mask(
      const amrex::Geometry& geometry,
      const amrex::MFIter& grid_iterator,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<std::uint8_t>& ibm_markers,
      const amrex::Array4<const amrex::Real>& pressure_reconstruction_states,
      const int pressure_reconstruction_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux,
      const DirectionalFaceFallbackMaskViews& fallback_mask)
  {
#ifdef AMREX_USE_GPIBM
    static_assert(
        ReconstructionMethod::supports_rz_paired_pressure_flux,
        "GP paired RZ pressure requires a frozen LLF reconstruction");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geometry.IsRZ(),
        "GP paired RZ pressure requires cylindrical RZ geometry");
    compute_characteristic_llf_face_fluxes<
        FaceStencilMode::ghost_point_aware, true, true>(
        geometry, grid_iterator.tilebox(), primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        pressure_reconstruction_states, pressure_reconstruction_component,
        ibm_markers, radial_pressure_face_flux, fallback_mask);
#else
    amrex::ignore_unused(
        geometry, grid_iterator, primitive_states,
        directional_face_flux_arrays, unused_rhs_state, closure,
        ibm_markers, pressure_reconstruction_states,
        pressure_reconstruction_component, radial_pressure_face_flux,
        fallback_mask);
    amrex::Abort(
        "paired R-Z fallback-mask capture requires pure shared-GP");
#endif
  }

  /*
   * This implementation function must remain public.  NVCC does not permit
   * an extended device lambda inside a private or protected member function.
   */
  /**
   * \brief Evaluate characteristic LLF fluxes on every face of a cell box.
   *
   * \tparam SelectedFaceStencilMode All-fluid or marker-aware GP operation.
   * \tparam IbmMarkerView Array view exposing the two IBM marker components.
   */
  template <FaceStencilMode SelectedFaceStencilMode,
            bool ComputeRzPairedPressure,
            bool RecordLocalLlfFallbackMask,
            typename IbmMarkerView>
  static void compute_characteristic_llf_face_fluxes(
      const amrex::Geometry& geometry,
      const amrex::Box& cell_box,
      const amrex::Array4<const amrex::Real>& primitive_states,
      const DirectionalFaceFluxArrays& directional_face_flux_arrays,
      const amrex::Array4<amrex::Real>& unused_rhs_state,
      const Closure* closure,
      const amrex::Array4<const amrex::Real>&
          pressure_reconstruction_states,
      const int pressure_reconstruction_component,
      const IbmMarkerView& ibm_markers,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux,
      const DirectionalFaceFallbackMaskViews& fallback_mask_views)
  {
    using amrex::IntVect;
    using amrex::Real;

    const auto cell_size = geometry.CellSizeArray();
    const auto prob_lo = geometry.ProbLoArray();
    const bool is_rz = geometry.IsRZ();
    const int radial_axis_face_index = geometry.Domain().smallEnd(0);
    const int radial_domain_high_index = geometry.Domain().bigEnd(0);
    amrex::ignore_unused(unused_rhs_state);
    if constexpr (ComputeRzPairedPressure) {
      static_assert(
          ReconstructionMethod::supports_rz_paired_pressure_flux,
          "Paired RZ pressure requires reusable frozen LLF weights");
    } else {
      amrex::ignore_unused(
          pressure_reconstruction_states,
          pressure_reconstruction_component,
          radial_pressure_face_flux);
    }

    // Optional face-scaled WENO-Z regularisation:
    // epsilon_f = C * S_f^2, with S_f the maximum absolute characteristic
    // split-branch sample at the face. Cartesian reconstruction retains the
    // fixed smoothness_epsilon() unless the global option is set.
    static const Real llf_weno_epsilon_relative = [] {
      amrex::ParmParse parmparse("cns");
      amrex::Real value = 0.0;
      parmparse.query("llf_weno_epsilon_relative", value);
      if (!std::isfinite(value) || value < Real(0.0) ||
          value >= Real(1.0)) {
        amrex::Abort(
            "cns.llf_weno_epsilon_relative must be finite and in [0,1)");
      }
      return Real(value);
    }();
    // R-Z fields contain regular critical points at the axis, and the radial
    // metric split variables r(U +/- F/alpha) additionally vanish with
    // radius.  A fixed dimensional epsilon does not preserve the uniform
    // critical-point scaling of the WENO-Z weights.  The R-Z coefficient is
    // used in both coordinate directions, as in the established R-Z path.
    // The coefficient is applied uniformly at every reconstructed face.
    static const Real rz_weno_epsilon_relative = [] {
      amrex::ParmParse parmparse("cns");
      amrex::Real value = 1.0e-6;
      parmparse.query("rz_weno_epsilon_relative", value);
      if (!std::isfinite(value) || value < Real(0.0) ||
          value >= Real(1.0)) {
        amrex::Abort(
            "cns.rz_weno_epsilon_relative must be finite and in [0,1)");
      }
      return Real(value);
    }();
    // Device lambdas cannot reference function-local statics; capture the
    // resolved values by value instead.
    const Real llf_weno_epsilon_relative_v = llf_weno_epsilon_relative;
    const Real rz_weno_epsilon_relative_v = rz_weno_epsilon_relative;
    const LocalLlfFallbackOptions fallback_options =
        local_llf_fallback_options();
    const int local_llf_fallback_enabled = fallback_options.enabled;
    const Real fallback_density_jump_threshold =
        fallback_options.density_jump_threshold;
    const Real fallback_density_curvature_threshold =
        fallback_options.density_curvature_threshold;
    const Real fallback_axis_row_curvature_threshold =
        fallback_options.axis_row_curvature_threshold;

    // Each direction owns an independent face-centred flux array.
    for (int direction_index = 0;
         direction_index < AMREX_SPACEDIM;
         ++direction_index) {
      // face_box contains N+1 faces in the current coordinate direction.
      const amrex::Box face_box =
          amrex::surroundingNodes(cell_box, direction_index);
      // face_flux is the writable face-centred conservative flux view.
      const auto face_flux =
          directional_face_flux_arrays[direction_index]->array();
      const auto fallback_mask =
          fallback_mask_views[direction_index];
      const bool reconstruct_radial_metric_flux =
          is_rz && direction_index == 0;
      const Real radial_origin = prob_lo[0];
      const Real radial_spacing = cell_size[0];

      amrex::ParallelFor(
          face_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // NVCC requires variables used below an if constexpr branch to
            // be captured before that branch is entered.  These no-op uses
            // keep the CPU and GPU instantiations on the same code path.
            amrex::ignore_unused(
                is_rz, llf_weno_epsilon_relative_v,
                rz_weno_epsilon_relative_v, radial_pressure_face_flux,
                local_llf_fallback_enabled,
                fallback_density_jump_threshold,
                fallback_density_curvature_threshold,
                fallback_axis_row_curvature_threshold, fallback_mask,
                pressure_reconstruction_states,
                pressure_reconstruction_component,
                radial_domain_high_index);
            // face_index identifies the current Cartesian face.
            const IntVect face_index(AMREX_D_DECL(i, j, k));
            // Subtracting direction_offset gives the lower-index face cell.
            const IntVect direction_offset(
                IntVect::TheDimensionVector(direction_index));
            amrex::ignore_unused(ibm_markers, face_flux);
            if constexpr (RecordLocalLlfFallbackMask) {
              fallback_mask(face_index, 0) = 0;
            }

            // In R-Z, reconstruct the metric flux r*F_r rather than
            // reconstructing F_r and multiplying the completed face flux by
            // r.  The face array retains its per-unit-area convention by
            // storing the metric reconstruction divided by r_face.
            Real radial_face = Real(1.0);
            if (reconstruct_radial_metric_flux) {
              radial_face =
                  radial_origin +
                  Real(face_index[0] - radial_axis_face_index) *
                      radial_spacing;
            }
            // The lower index is an axis only when the physical radial origin
            // is exactly zero.  An annular R-Z domain may also have this mesh
            // index and must not receive an axis closure.
            const bool radial_axis_face =
                reconstruct_radial_metric_flux &&
                face_index[0] == radial_axis_face_index &&
                radial_origin == Real(0.0);

            // stencil_state_available marks fluid and reconstructed GP states.
            bool stencil_state_available[2 * required_ghost_cells];
            for (int sample_index = 0;
                 sample_index < 2 * required_ghost_cells;
                 ++sample_index) {
              stencil_state_available[sample_index] = true;
            }

            if constexpr (
                SelectedFaceStencilMode == FaceStencilMode::ghost_point_aware) {
              // Faces with two geometrically solid neighbours do not update
              // any active fluid cell and therefore carry zero flux.
              IntVect lower_marker_index = face_index - direction_offset;
              IntVect upper_marker_index = face_index;
              if (radial_axis_face) {
                lower_marker_index[0] =
                    2 * radial_axis_face_index - 1 - lower_marker_index[0];
              }
              if (Reconstruction::is_solid_cell(
                      lower_marker_index, ibm_markers) &&
                  Reconstruction::is_solid_cell(
                      upper_marker_index, ibm_markers)) {
                Reconstruction::set_face_flux_to_zero<
                    decltype(face_flux), Closure>(face_index, face_flux);
                if constexpr (ComputeRzPairedPressure) {
                  if (reconstruct_radial_metric_flux) {
                    radial_pressure_face_flux(face_index, 0) = Real(0.0);
                  }
                }
                return;
              }

              // The six driver samples span offsets -3 through +2 from the
              // face index when required_ghost_cells equals three.
              for (int sample_index = 0;
                   sample_index < 2 * required_ghost_cells;
                   ++sample_index) {
                const IntVect stencil_cell_index =
                    face_index +
                    (sample_index - required_ghost_cells) *
                        direction_offset;
                IntVect marker_cell_index = stencil_cell_index;
                if (reconstruct_radial_metric_flux &&
                    radial_origin == Real(0.0) &&
                    marker_cell_index[0] < radial_axis_face_index) {
                  // A negative-radius physical ghost is the parity image of
                  // one positive-radius background cell.  Its availability
                  // must therefore follow the mirrored fluid/shared-GP
                  // marker.  Treating the geometry halo's default marker as
                  // fluid could otherwise expose a stale deep-solid value.
                  marker_cell_index[0] =
                      2 * radial_axis_face_index - 1 -
                      marker_cell_index[0];
                }
                stencil_state_available[sample_index] =
                    Reconstruction::is_flux_stencil_state_available(
                        marker_cell_index, ibm_markers);
              }
            }

            // reconstruction_plan disables candidates that cross an
            // unavailable deep-solid state.
            const ReconstructionPlan reconstruction_plan =
                Reconstruction::build_face_reconstruction_plan<
                    ReconstructionMethod>(stencil_state_available);

            // The adjacent cells have sample indices radius-1 and radius.
            const bool adjacent_states_available =
                stencil_state_available[required_ghost_cells - 1] &&
                stencil_state_available[required_ghost_cells];
            // Marker availability and thermodynamic validity are separate.
            // This prevents a missing GP state from being read silently.
            const bool adjacent_thermodynamic_values_are_valid =
                adjacent_states_available &&
                Reconstruction::face_adjacent_thermodynamic_values_are_valid<
                    decltype(primitive_states), Closure>(
                    face_index, direction_offset, primitive_states);
            AMREX_ASSERT_WITH_MESSAGE(
                adjacent_thermodynamic_values_are_valid,
                "LLF reconstruction found an invalid adjacent state");
            if (!adjacent_thermodynamic_values_are_valid) {
              Reconstruction::mark_face_flux_invalid<
                  decltype(face_flux), Closure>(face_index, face_flux);
              if constexpr (ComputeRzPairedPressure) {
                if (reconstruct_radial_metric_flux) {
                  radial_pressure_face_flux(face_index, 0) =
                      std::numeric_limits<Real>::quiet_NaN();
                }
              }
              return;
            }

            // llf_max_wave_speed is one scalar alpha shared by all conserved
            // components at this face.
            Real llf_max_wave_speed = Real(0.0);
            // Use every readable fluid/GP state in the six-point driver
            // stencil. Deep-solid states have no reconstructed primitive state
            // and are skipped.
            for (int sample_index = 0;
                 sample_index < 2 * required_ghost_cells;
                 ++sample_index) {
              if (!stencil_state_available[sample_index]) {
                continue;
              }

              const IntVect stencil_cell_index =
                  face_index +
                  (sample_index - required_ghost_cells) * direction_offset;
              const Real local_wave_speed =
                  std::abs(primitive_states(
                      stencil_cell_index,
                      Closure::QU + direction_index)) +
                  primitive_states(stencil_cell_index, Closure::QC);
              llf_max_wave_speed =
                  amrex::max(llf_max_wave_speed, local_wave_speed);
            }

            AMREX_ASSERT_WITH_MESSAGE(
                llf_max_wave_speed > Real(0.0),
                "LLF reconstruction found a non-positive wave speed");
            if (!(llf_max_wave_speed > Real(0.0)) ||
                !std::isfinite(llf_max_wave_speed)) {
              Reconstruction::mark_face_flux_invalid<
                  decltype(face_flux), Closure>(face_index, face_flux);
              if constexpr (ComputeRzPairedPressure) {
                if (reconstruct_radial_metric_flux) {
                  radial_pressure_face_flux(face_index, 0) =
                      std::numeric_limits<Real>::quiet_NaN();
                }
              }
              return;
            }

            /*
             * A relative density jump alone is O(dx) in a smooth Cartesian
             * field and would activate on the practical MMS grids.  The
             * general gate therefore also requires a normalized density
             * second difference.  The first interior RZ radial face retains
             * the independently verified density-jump gate because even
             * parity makes its smooth jump O(dr^2), while the near-axis
             * contact instability is not accompanied by a pressure jump.
             */
            bool use_local_first_order_llf = false;
            if (local_llf_fallback_enabled != 0 &&
                !radial_axis_face) {
              const IntVect left_cell_index =
                  face_index - direction_offset;
              const IntVect right_cell_index = face_index;
              const Real density_left =
                  primitive_states(left_cell_index, Closure::QRHO);
              const Real density_right =
                  primitive_states(right_cell_index, Closure::QRHO);
              const Real density_floor =
                  amrex::max(
                      amrex::min(density_left, density_right),
                      std::numeric_limits<Real>::min());
              const Real relative_density_jump =
                  std::abs(density_right - density_left) / density_floor;

              if (relative_density_jump >
                  fallback_density_jump_threshold) {
                Real density_scale = Real(0.0);
                for (int sample_index = 0;
                     sample_index < 2 * required_ghost_cells;
                     ++sample_index) {
                  if (!stencil_state_available[sample_index]) {
                    continue;
                  }
                  const IntVect stencil_cell_index =
                      face_index +
                      (sample_index - required_ghost_cells) *
                          direction_offset;
                  density_scale = amrex::max(
                      density_scale,
                      std::abs(primitive_states(
                          stencil_cell_index, Closure::QRHO)));
                }

                /*
                 * Axial faces of the axis-adjacent cell row in RZ.  The
                 * discrete divergence there is
                 *   (r_{i+1/2} F_{i+1/2} - r_{i-1/2} F_{i-1/2}) / (r_i dr)
                 * with r_{i-1/2} = 0 and r_i = dr/2, so the outer radial flux
                 * is amplified by 2/dr.  That amplification doubles when dr is
                 * halved, while the generic curvature threshold is a
                 * grid-independent constant; a shock sweeping axially through
                 * this row therefore escapes the fallback on fine grids.  The
                 * existing collar only covers the first interior RADIAL face,
                 * so tighten the curvature test here and only here.
                 */
                const bool axis_row_axial_face =
                    is_rz && direction_index != 0 &&
                    face_index[0] == radial_axis_face_index;
                const Real effective_curvature_threshold =
                    axis_row_axial_face
                        ? fallback_axis_row_curvature_threshold
                        : fallback_density_curvature_threshold;

                bool density_is_nonsmooth = false;
                // Only these two groups contain both face-adjacent cells.
                // A remote discontinuity must not create a broad first-order
                // collar by itself.
                for (int group_index = required_ghost_cells - 2;
                     group_index <= required_ghost_cells - 1;
                     ++group_index) {
                  if (!stencil_state_available[group_index] ||
                      !stencil_state_available[group_index + 1] ||
                      !stencil_state_available[group_index + 2]) {
                    continue;
                  }
                  const IntVect first_cell =
                      face_index +
                      (group_index - required_ghost_cells) *
                          direction_offset;
                  const IntVect second_cell =
                      first_cell + direction_offset;
                  const IntVect third_cell =
                      second_cell + direction_offset;
                  const Real first_density =
                      primitive_states(first_cell, Closure::QRHO);
                  const Real second_density =
                      primitive_states(second_cell, Closure::QRHO);
                  const Real third_density =
                      primitive_states(third_cell, Closure::QRHO);
                  const Real curvature_denominator =
                      std::abs(first_density) +
                      Real(2.0) * std::abs(second_density) +
                      std::abs(third_density) + density_scale +
                      std::numeric_limits<Real>::min();
                  const Real normalized_density_curvature =
                      std::abs(
                          first_density - Real(2.0) * second_density +
                          third_density) /
                      curvature_denominator;
                  density_is_nonsmooth =
                      density_is_nonsmooth ||
                      normalized_density_curvature >
                          effective_curvature_threshold;
                }

                const bool first_interior_rz_radial_face =
                    reconstruct_radial_metric_flux &&
                    face_index[0] == radial_axis_face_index + 1;
                use_local_first_order_llf =
                    density_is_nonsmooth ||
                    first_interior_rz_radial_face;
              }
            }

            if (use_local_first_order_llf) {
              if constexpr (RecordLocalLlfFallbackMask) {
                fallback_mask(face_index, 0) = 1;
              }
              const IntVect left_cell_index =
                  face_index - direction_offset;
              const IntVect right_cell_index = face_index;
              Real conservative_left[Closure::NCONS];
              Real conservative_right[Closure::NCONS];
              Real physical_flux_left[Closure::NCONS];
              Real physical_flux_right[Closure::NCONS];
              closure->prims2cons(
                  left_cell_index, primitive_states, conservative_left);
              closure->prims2cons(
                  right_cell_index, primitive_states, conservative_right);
              closure->prims2flux(
                  left_cell_index, direction_index, primitive_states,
                  physical_flux_left);
              closure->prims2flux(
                  right_cell_index, direction_index, primitive_states,
                  physical_flux_right);

              Real pressure_face = Real(0.0);
              if constexpr (ComputeRzPairedPressure) {
                if (reconstruct_radial_metric_flux) {
                  const Real pressure_left =
                      pressure_reconstruction_states(
                          left_cell_index,
                          pressure_reconstruction_component);
                  const Real pressure_right =
                      pressure_reconstruction_states(
                          right_cell_index,
                          pressure_reconstruction_component);
                  physical_flux_left[Closure::UMX] -= pressure_left;
                  physical_flux_right[Closure::UMX] -= pressure_right;
                  pressure_face =
                      Real(0.5) * (pressure_left + pressure_right);
                }
              }

              if (reconstruct_radial_metric_flux) {
                const Real left_radius =
                    radial_origin +
                    (Real(left_cell_index[0] - radial_axis_face_index) +
                     Real(0.5)) *
                        radial_spacing;
                const Real right_radius =
                    radial_origin +
                    (Real(right_cell_index[0] - radial_axis_face_index) +
                     Real(0.5)) *
                        radial_spacing;
                const Real inverse_face_radius = Real(1.0) / radial_face;
                for (int component_index = 0;
                     component_index < Closure::NCONS;
                     ++component_index) {
                  face_flux(face_index, component_index) =
                      Real(0.5) * inverse_face_radius *
                          (left_radius * physical_flux_left[component_index] +
                           right_radius *
                               physical_flux_right[component_index]) -
                      Real(0.5) * llf_max_wave_speed *
                          (conservative_right[component_index] -
                           conservative_left[component_index]);
                }
              } else {
                for (int component_index = 0;
                     component_index < Closure::NCONS;
                     ++component_index) {
                  face_flux(face_index, component_index) =
                      Real(0.5) *
                          (physical_flux_left[component_index] +
                           physical_flux_right[component_index]) -
                      Real(0.5) * llf_max_wave_speed *
                          (conservative_right[component_index] -
                           conservative_left[component_index]);
                }
              }

              if constexpr (ComputeRzPairedPressure) {
                if (reconstruct_radial_metric_flux) {
                  face_flux(face_index, Closure::UMX) += pressure_face;
                  radial_pressure_face_flux(face_index, 0) = pressure_face;
                }
              }
              return;
            }

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                reconstruction_plan.high_order_reconstruction_available,
                "LLF-WENO reconstruction has no complete three-point "
                "candidate on both flux-split branches");
            if (!reconstruction_plan.high_order_reconstruction_available) {
              // Device assertions may be disabled in optimized GPU builds.
              // Preserve fail-closed behavior and never apply a two-state flux.
              Reconstruction::mark_face_flux_invalid<
                  decltype(face_flux), Closure>(face_index, face_flux);
              if constexpr (ComputeRzPairedPressure) {
                if (reconstruct_radial_metric_flux) {
                  radial_pressure_face_flux(face_index, 0) =
                      std::numeric_limits<Real>::quiet_NaN();
                }
              }
              return;
            }

            // face_roe_state defines one characteristic basis for both LLF
            // branches at the current face.
            const auto face_roe_state = closure->roe_avg_state(
                face_index, direction_index, primitive_states);
            // Temporary conservative state and physical flux for one sample.
            Real conservative_state[Closure::NCONS];
            Real physical_flux[Closure::NCONS];
            // These arrays store 0.5*(U+F/alpha) and 0.5*(U-F/alpha).
            // The final flux is alpha times the positive-minus-negative result.
            Real positive_split_samples[2 * required_ghost_cells]
                                       [Closure::NCONS] = {};
            Real negative_split_samples[2 * required_ghost_cells]
                                       [Closure::NCONS] = {};
            // The RZ pressure pair reuses the metric multiplier applied to
            // each cell-centred sample.  Cartesian and axial values are one.
            Real sample_metric_scales[2 * required_ghost_cells];
            for (int sample_index = 0;
                 sample_index < 2 * required_ghost_cells;
                 ++sample_index) {
              sample_metric_scales[sample_index] = Real(1.0);
            }

            for (int sample_index = 0;
                 sample_index < 2 * required_ghost_cells;
                 ++sample_index) {
              if (!reconstruction_plan.sample_state_required[sample_index]) {
                continue;
              }

              const IntVect stencil_cell_index =
                  face_index +
                  (sample_index - required_ghost_cells) *
                      direction_offset;
              closure->prims2flux(
                  stencil_cell_index, direction_index,
                  primitive_states, physical_flux);
              closure->prims2cons(
                  stencil_cell_index, primitive_states,
                  conservative_state);

              Real radial_metric_scale = Real(1.0);
              if (reconstruct_radial_metric_flux) {
                const Real radial_cell =
                    radial_origin +
                    (Real(stencil_cell_index[0] - radial_axis_face_index) +
                     Real(0.5)) *
                        radial_spacing;
                // Away from the axis the face array stores h_(rF)/r_face.
                // At r_face=0 it stores h_(rF) itself: unlike the physical
                // metric flux, this finite-difference h-flux need not vanish.
                radial_metric_scale = radial_axis_face
                    ? radial_cell
                    : radial_cell / radial_face;
              }
              sample_metric_scales[sample_index] = radial_metric_scale;

              for (int component_index = 0;
                   component_index < Closure::NCONS;
                   ++component_index) {
                // Divide F by the scalar LLF wave speed for both branches.
                const Real physical_flux_divided_by_wave_speed =
                    physical_flux[component_index] /
                    llf_max_wave_speed;
                if (reconstruct_radial_metric_flux) {
                  positive_split_samples[sample_index][component_index] =
                      Real(0.5) * radial_metric_scale *
                      (conservative_state[component_index] +
                       physical_flux_divided_by_wave_speed);
                  negative_split_samples[sample_index][component_index] =
                      Real(0.5) * radial_metric_scale *
                      (conservative_state[component_index] -
                       physical_flux_divided_by_wave_speed);
                } else {
                  positive_split_samples[sample_index][component_index] =
                      Real(0.5) *
                      (conservative_state[component_index] +
                       physical_flux_divided_by_wave_speed);
                  negative_split_samples[sample_index][component_index] =
                      Real(0.5) *
                      (conservative_state[component_index] -
                       physical_flux_divided_by_wave_speed);
                }
              }
              closure->cons2char(
                  face_roe_state, positive_split_samples[sample_index]);
              closure->cons2char(
                  face_roe_state, negative_split_samples[sample_index]);
            }

            // The scale-aware value is one face-wide epsilon formed from
            // every characteristic family and both LLF branches.
            Real face_epsilon = ReconScheme::smoothness_epsilon();
            if constexpr (
                ReconstructionMethod::uses_weno_regularisation_epsilon) {
              const Real direction_relative_epsilon =
                  is_rz ? rz_weno_epsilon_relative_v
                        : llf_weno_epsilon_relative_v;
              const bool use_scaled_epsilon =
                  direction_relative_epsilon > Real(0.0);
              if (use_scaled_epsilon) {
                Real face_scale = Real(0.0);
                for (int sample_index = 0;
                     sample_index < 2 * required_ghost_cells;
                     ++sample_index) {
                  if (!reconstruction_plan
                           .sample_state_required[sample_index]) {
                    continue;
                  }
                  for (int component_index = 0;
                       component_index < Closure::NCONS;
                       ++component_index) {
                    face_scale = amrex::max(
                        face_scale,
                        std::abs(positive_split_samples[sample_index]
                                                       [component_index]));
                    face_scale = amrex::max(
                        face_scale,
                        std::abs(negative_split_samples[sample_index]
                                                       [component_index]));
                  }
                }
                face_epsilon = amrex::max(
                    std::numeric_limits<Real>::min(),
                    direction_relative_epsilon * face_scale * face_scale);
              }
            }

            if constexpr (ComputeRzPairedPressure) {
              if (reconstruct_radial_metric_flux) {
                using FrozenWeightSet =
                    typename ReconstructionMethod::FrozenWeightSet;
                // These three arrays hold completed characteristic fluxes.
                // Processing one characteristic family at a time avoids
                // retaining 2*NCONS sets of nonlinear weights in every GPU
                // thread.
                Real reconstructed_full_flux[Closure::NCONS];
                Real reconstructed_advective_flux[Closure::NCONS];
                Real reconstructed_pressure_flux[Closure::NCONS];
                Real stencil_values[
                    ReconstructionMethod::reconstruction_sample_count];

                // L(P)/(2 alpha), where P contains only radial pressure in
                // the radial-momentum flux component.
                Real pressure_characteristic_samples[
                    2 * required_ghost_cells][Closure::NCONS] = {};
                for (int sample_index = 0;
                     sample_index < 2 * required_ghost_cells;
                     ++sample_index) {
                  if (!reconstruction_plan
                           .sample_state_required[sample_index]) {
                    continue;
                  }
                  const IntVect stencil_cell_index =
                      face_index +
                      (sample_index - required_ghost_cells) *
                          direction_offset;
                  pressure_characteristic_samples[sample_index]
                                                     [Closure::UMX] =
                      Real(0.5) *
                      pressure_reconstruction_states(
                          stencil_cell_index,
                          pressure_reconstruction_component) /
                      llf_max_wave_speed;
                  closure->cons2char(
                      face_roe_state,
                      pressure_characteristic_samples[sample_index]);
                }

                // Freeze and reuse the positive and negative weights one
                // characteristic family at a time.  The reconstruction is
                // linear once these weights have been fixed.
                for (int component_index = 0;
                     component_index < Closure::NCONS;
                     ++component_index) {
                  ReconstructionMethod::gather_positive_split_stencil(
                      component_index, positive_split_samples,
                      stencil_values);
                  const FrozenWeightSet positive_weights =
                      ReconstructionMethod::freeze_available_stencil_weights(
                          stencil_values,
                          reconstruction_plan.positive_linear_weights,
                          face_epsilon);
                  const Real reconstructed_full_positive =
                      ReconstructionMethod::reconstruct_with_frozen_weights(
                          stencil_values, positive_weights);

                  ReconstructionMethod::gather_negative_split_stencil(
                      component_index, negative_split_samples,
                      stencil_values);
                  const FrozenWeightSet negative_weights =
                      ReconstructionMethod::freeze_available_stencil_weights(
                          stencil_values,
                          reconstruction_plan.negative_linear_weights,
                          face_epsilon);
                  const Real reconstructed_full_negative =
                      ReconstructionMethod::reconstruct_with_frozen_weights(
                          stencil_values, negative_weights);

                  // Remove metric pressure while retaining the full LLF U
                  // dissipation in both split branches.
                  for (int sample_index = 0;
                       sample_index < 2 * required_ghost_cells;
                       ++sample_index) {
                    const Real metric_pressure_sample =
                        sample_metric_scales[sample_index] *
                        pressure_characteristic_samples[sample_index]
                                                           [component_index];
                    positive_split_samples[sample_index][component_index] -=
                        metric_pressure_sample;
                    negative_split_samples[sample_index][component_index] +=
                        metric_pressure_sample;
                  }

                  ReconstructionMethod::gather_positive_split_stencil(
                      component_index, positive_split_samples,
                      stencil_values);
                  const Real reconstructed_advective_positive =
                      ReconstructionMethod::reconstruct_with_frozen_weights(
                          stencil_values, positive_weights);

                  ReconstructionMethod::gather_negative_split_stencil(
                      component_index, negative_split_samples,
                      stencil_values);
                  const Real reconstructed_advective_negative =
                      ReconstructionMethod::reconstruct_with_frozen_weights(
                          stencil_values, negative_weights);

                  // Apply the same weights to the nonmetric pressure-only
                  // split.  Its ordinary face difference is the pressure
                  // derivative used by the paired RZ momentum operator.
                  for (int sample_index = 0;
                       sample_index < 2 * required_ghost_cells;
                       ++sample_index) {
                    const Real pressure_sample =
                        pressure_characteristic_samples[sample_index]
                                                           [component_index];
                    positive_split_samples[sample_index][component_index] =
                        pressure_sample;
                    negative_split_samples[sample_index][component_index] =
                        -pressure_sample;
                  }

                  ReconstructionMethod::gather_positive_split_stencil(
                      component_index, positive_split_samples,
                      stencil_values);
                  const Real reconstructed_pressure_positive =
                      ReconstructionMethod::reconstruct_with_frozen_weights(
                          stencil_values, positive_weights);

                  ReconstructionMethod::gather_negative_split_stencil(
                      component_index, negative_split_samples,
                      stencil_values);
                  const Real reconstructed_pressure_negative =
                      ReconstructionMethod::reconstruct_with_frozen_weights(
                          stencil_values, negative_weights);

                  reconstructed_full_flux[component_index] =
                      llf_max_wave_speed *
                      (reconstructed_full_positive -
                       reconstructed_full_negative);
                  reconstructed_advective_flux[component_index] =
                      llf_max_wave_speed *
                      (reconstructed_advective_positive -
                       reconstructed_advective_negative);
                  reconstructed_pressure_flux[component_index] =
                      llf_max_wave_speed *
                      (reconstructed_pressure_positive -
                       reconstructed_pressure_negative);
                }

                closure->char2cons(
                    face_roe_state, reconstructed_full_flux);
                closure->char2cons(
                    face_roe_state, reconstructed_advective_flux);
                closure->char2cons(
                    face_roe_state, reconstructed_pressure_flux);

                for (int component_index = 0;
                     component_index < Closure::NCONS;
                     ++component_index) {
                  face_flux(face_index, component_index) =
                      reconstructed_full_flux[component_index];
                }

                Real pressure_face =
                    reconstructed_pressure_flux[Closure::UMX];
                if (radial_axis_face) {
                  // This is the parity-optimal finite-difference auxiliary
                  // pressure flux h_p(0), not the physical point value p(0).
                  // For a complete mirrored stencil it is the common linear
                  // limit of WENO-Z5, TENO5, and TENO6.  A marker-reduced GP
                  // stencil degrades only as far as its available positive-r
                  // centre states permit; no deep-solid value is read.
                  const int ring1_sample = required_ghost_cells;
                  const int ring2_sample = required_ghost_cells + 1;
                  const int ring3_sample = required_ghost_cells + 2;
                  const bool ring1_available =
                      face_index[0] <= radial_domain_high_index &&
                      stencil_state_available[ring1_sample];
                  const bool ring2_available =
                      face_index[0] + 1 <= radial_domain_high_index &&
                      stencil_state_available[ring2_sample];
                  const bool ring3_available =
                      face_index[0] + 2 <= radial_domain_high_index &&
                      stencil_state_available[ring3_sample];
                  AMREX_ASSERT_WITH_MESSAGE(
                      ring1_available,
                      "Active R-Z axis face lacks its first-ring state");
                  if (!ring1_available) {
                    Reconstruction::mark_face_flux_invalid<
                        decltype(face_flux), Closure>(face_index, face_flux);
                    radial_pressure_face_flux(face_index, 0) =
                        std::numeric_limits<Real>::quiet_NaN();
                    return;
                  }
                  const Real pressure_ring1 =
                      pressure_reconstruction_states(
                          face_index, pressure_reconstruction_component);
                  if (!std::isfinite(pressure_ring1)) {
                    Reconstruction::mark_face_flux_invalid<
                        decltype(face_flux), Closure>(face_index, face_flux);
                    radial_pressure_face_flux(face_index, 0) =
                        std::numeric_limits<Real>::quiet_NaN();
                    return;
                  }
                  pressure_face = pressure_ring1;
                  if (ring2_available) {
                    const Real pressure_ring2 =
                        pressure_reconstruction_states(
                            face_index + direction_offset,
                            pressure_reconstruction_component);
                    if (!std::isfinite(pressure_ring2)) {
                      Reconstruction::mark_face_flux_invalid<
                          decltype(face_flux), Closure>(face_index, face_flux);
                      radial_pressure_face_flux(face_index, 0) =
                          std::numeric_limits<Real>::quiet_NaN();
                      return;
                    }
                    pressure_face =
                        (Real(7.0) * pressure_ring1 - pressure_ring2) /
                        Real(6.0);
                    if (ring3_available) {
                      const Real pressure_ring3 =
                          pressure_reconstruction_states(
                              face_index + 2 * direction_offset,
                              pressure_reconstruction_component);
                      if (!std::isfinite(pressure_ring3)) {
                        Reconstruction::mark_face_flux_invalid<
                            decltype(face_flux), Closure>(face_index, face_flux);
                        radial_pressure_face_flux(face_index, 0) =
                            std::numeric_limits<Real>::quiet_NaN();
                        return;
                      }
                      pressure_face =
                          (Real(37.0) * pressure_ring1 -
                           Real(8.0) * pressure_ring2 + pressure_ring3) /
                          Real(30.0);
                    }
                  }
                }
                const Real metric_advective_flux_per_radius =
                    reconstructed_advective_flux[Closure::UMX];
                radial_pressure_face_flux(face_index, 0) = pressure_face;

                // At non-axis faces the standard array remains a complete
                // per-area Euler flux.  The axis slot instead stores the
                // metric flux, whose regular radial-momentum limit is zero.
                face_flux(face_index, Closure::UMX) =
                    radial_axis_face
                        ? Real(0.0)
                        : metric_advective_flux_per_radius + pressure_face;
                return;
              }
            }

            // Reconstructed split states remain in characteristic space until
            // all characteristic components have been processed.
            Real reconstructed_positive_split[Closure::NCONS];
            Real reconstructed_negative_split[Closure::NCONS];
            // stencil_values is the ordered one-dimensional input expected by
            // the selected WENO or TENO reconstruction method.
            Real stencil_values[
                ReconstructionMethod::reconstruction_sample_count];
            for (int component_index = 0;
                 component_index < Closure::NCONS;
                 ++component_index) {
              ReconstructionMethod::gather_positive_split_stencil(
                  component_index, positive_split_samples, stencil_values);
              reconstructed_positive_split[component_index] =
                  reconstruction_plan.full_stencil_available
                      ? ReconstructionMethod::reconstruct_interface_value(
                            stencil_values, face_epsilon)
                      : ReconstructionMethod::reconstruct_interface_value(
                            stencil_values,
                            reconstruction_plan.positive_linear_weights,
                            face_epsilon);

              ReconstructionMethod::gather_negative_split_stencil(
                  component_index, negative_split_samples, stencil_values);
              reconstructed_negative_split[component_index] =
                  reconstruction_plan.full_stencil_available
                      ? ReconstructionMethod::reconstruct_interface_value(
                            stencil_values, face_epsilon)
                      : ReconstructionMethod::reconstruct_interface_value(
                            stencil_values,
                            reconstruction_plan.negative_linear_weights,
                            face_epsilon);
            }

            closure->char2cons(
                face_roe_state, reconstructed_positive_split);
            closure->char2cons(
                face_roe_state, reconstructed_negative_split);

            // Combine the two reconstructed branches into the numerical flux.
            for (int component_index = 0;
                 component_index < Closure::NCONS;
                 ++component_index) {
              face_flux(face_index, component_index) =
                  llf_max_wave_speed *
                  (reconstructed_positive_split[component_index] -
                   reconstructed_negative_split[component_index]);
            }

            // Every non-paired LLF-WENO/TENO R-Z path, including reduced
            // shared-GP stencils, receives the same exact radial-momentum
            // parity projection.  Other components retain their finite-
            // difference auxiliary metric h-flux and must not be blanked.
            if (radial_axis_face) {
              face_flux(face_index, Closure::UMX) = Real(0.0);
            }

          });
    }
  }
};

/**
 * \brief Public compatibility name used by existing prob.h files.
 *
 * The implementation class name states that this is a characteristic local
 * LLF flux driver.  The legacy weno_t spelling remains source-compatible.
 */
template <typename ReconstructionMethod, typename Closure>
using weno_t =
    characteristic_llf_flux_t<ReconstructionMethod, Closure>;

/// Characteristic LLF driver with fifth-order WENO-Z reconstruction.
template <typename Closure>
using llf_wenoz5_t =
    characteristic_llf_flux_t<Reconstruction::WenoZ5, Closure>;

/// Characteristic LLF driver with fifth-order TENO reconstruction.
template <typename Closure>
using llf_teno5_t =
    characteristic_llf_flux_t<Reconstruction::Teno5, Closure>;

/// Characteristic LLF driver with sixth-order TENO reconstruction.
template <typename Closure>
using llf_teno6_t =
    characteristic_llf_flux_t<Reconstruction::Teno6, Closure>;

#endif
