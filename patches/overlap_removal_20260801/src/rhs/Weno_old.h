#ifndef CERISSE_WENO_OLD_H_
#define CERISSE_WENO_OLD_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include "Closures.h"
#include "IBMSharedGPFluxUtils.h"
#include "RZFiniteVolume.h"
#ifdef CNS_IBM_VALIDATION
#include "rz_shock_face_diagnostic.h"
#endif

#define POWER2(x) ((x) * (x))
#define POWER6(x) ((x) * (x) * (x) * (x) * (x) * (x))

// Choices are WenoZ5, Teno5. (Teno6 is slightly unstable)
namespace OldReconScheme {
struct WenoZ5 {
  static constexpr int ng = 3;
  static constexpr int ncand = 3;
  static constexpr const char* scheme_name =
      "characteristic_llf_weno_z5_old";

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  nonlinear_factor(const amrex::Real tau,
                   const amrex::Real beta,
                   const amrex::Real epsilon,
                   const bool squared_ratio) noexcept {
    using amrex::Real;
    const Real ratio = tau / (epsilon + beta);
    return Real(1.0) +
           (squared_ratio ? ratio * ratio : ratio);
  }

  static constexpr amrex::Real legacy_epsilon() noexcept {
    using amrex::Real;
    return (std::numeric_limits<Real>::digits >= 53)
               ? Real(1e-40)
               : std::numeric_limits<Real>::epsilon();
  }

  /**
   * @brief From a 2*ng x NCONS array, extract the left stencil.
   * @param[in] n  Index of the component.
   * @param[in] fp Flux array.
   * @param[out] s Stencil array [i+1, i, i-1, i-2, i-3].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void left_stencil(
      int n, amrex::Real const fp[2 * ng][NCONS], amrex::Real s[5]) noexcept {
    for (int m = 1; m < 2 * ng; ++m) s[m - 1] = fp[2 * ng - 1 - m][n];
  }

  /**
   * @brief From a 2*ng x NCONS array, extract the right stencil.
   * @param[in] n  Index of the component.
   * @param[in] fm Flux array.
   * @param[out] s Stencil array [i-2, i-1, i, i+1, i+2].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void right_stencil(
      int n, amrex::Real const fm[2 * ng][NCONS], amrex::Real s[5]) noexcept {
    for (int m = 1; m < 2 * ng; ++m) s[m - 1] = fm[m][n];
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void smoothness_indicator(
      const amrex::Real s[5], amrex::Real beta[3]) noexcept {
    // constexpr Real w13o12 = 13.0 / 12.0;
    // constexpr Real w1o4 = 0.25;
    // beta[2] = w13o12 * POWER2(s[4] - 2.0 * s[3] + s[2]) +
    //           w1o4 * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
    // beta[1] =
    //     w13o12 * POWER2(s[3] - 2.0 * s[2] + s[1]) + w1o4 * POWER2(s[3] - s[1]);
    // beta[0] = w13o12 * POWER2(s[2] - 2.0 * s[1] + s[0]) +
    //           w1o4 * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);

    beta[2] = Real(13. / 12.) * POWER2(s[4] - 2.0 * s[3] + s[2]) +
              Real(0.25) * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
    beta[1] = Real(13. / 12.) * POWER2(s[3] - 2.0 * s[2] + s[1]) +
              Real(0.25) * POWER2(s[3] - s[1]);
    beta[0] = Real(13. / 12.) * POWER2(s[2] - 2.0 * s[1] + s[0]) +
              Real(0.25) * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void linear_polynomial_recon(
      const amrex::Real s[5], amrex::Real vr[3]) noexcept {
    vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
    vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
    vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];
  }

  /**
   * \brief WENO-Z5 reconstruction of characteristic split fluxes for the
   * conservative finite-difference operator.
   * Ref https://doi.org/10.1016/j.jcp.2010.11.028.
   * \param[in] s Stencil values at cells [i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_masked_impl(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight,
      const amrex::Real epsilon,
      const bool squared_ratio) noexcept {
    using amrex::Real;
    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);

    Real vr[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real beta[3] = {Real(0.0), Real(0.0), Real(0.0)};

    // Do not even evaluate a disabled candidate.  A masked stencil can
    // contain NaN/garbage in an interior solid cell; evaluating that beta or
    // polynomial and multiplying by a zero weight still produces NaN.
    if (v2) {
      beta[2] = Real(13. / 12.) * POWER2(s[4] - 2.0 * s[3] + s[2]) +
                Real(0.25) * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
      vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
    }
    if (v1) {
      beta[1] = Real(13. / 12.) * POWER2(s[3] - 2.0 * s[2] + s[1]) +
                Real(0.25) * POWER2(s[3] - s[1]);
      vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
    }
    if (v0) {
      beta[0] = Real(13. / 12.) * POWER2(s[2] - 2.0 * s[1] + s[0]) +
                Real(0.25) * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);
      vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];
    }

    Real tau;
    if (v0 && v1 && v2) {
      tau = std::abs(beta[2] - beta[0]);
    } else if (v0 && v1) {
      tau = std::abs(beta[1] - beta[0]);
    } else if (v1 && v2) {
      tau = std::abs(beta[2] - beta[1]);
    } else {
      tau = Real(0.0);
    }

    if (v2) {
      beta[2] =
          nonlinear_factor(tau, beta[2], epsilon, squared_ratio) *
          optimal_weight[2];
    }
    if (v1) {
      beta[1] =
          nonlinear_factor(tau, beta[1], epsilon, squared_ratio) *
          optimal_weight[1];
    }
    if (v0) {
      beta[0] =
          nonlinear_factor(tau, beta[0], epsilon, squared_ratio) *
          optimal_weight[0];
    }

    const Real denom = beta[2] + beta[1] + beta[0];
    if (!(denom > std::numeric_limits<Real>::min())) {
      if (v1) return vr[1] / Real(6.0);
      if (v0) return vr[0] / Real(6.0);
      if (v2) return vr[2] / Real(6.0);
      return vr[1] / Real(6.0);
    }
    const Real inv_denom = 1.0 / denom;

    return inv_denom / Real(6.0) *
           (beta[2] * vr[2] + beta[1] * vr[1] + beta[0] * vr[0]);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_masked(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight) noexcept {
    return recon_masked_impl(
        s, optimal_weight, legacy_epsilon(), false);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  recon_masked_scaled(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight,
      const amrex::Real epsilon) noexcept {
    // Preserve the q=1 weighting used by Enson's characteristic LLF branch.
    // This diagnostic changes only the scale of epsilon, so its effect can be
    // separated from a change in the WENO-Z exponent.
    return recon_masked_impl(s, optimal_weight, epsilon, false);
  }

  // WENO-Z interpolation of pointwise split fluxes to the physical face.
  // The ordinary finite-difference numerical flux contains a face-constant
  // O(dx^2) term which cancels in Cartesian flux differences but not after an
  // R-Z face-area multiplication. These quadratic candidates instead return
  // the actual half-cell point value; the full optimal weights are 1:10:5 in
  // upwind-to-downwind order.
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  recon_point_masked(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight) noexcept {
    using amrex::Real;
    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);
    Real candidate[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real beta[3] = {Real(0.0), Real(0.0), Real(0.0)};
    if (v2) {
      beta[2] = Real(13.0 / 12.0) *
                    POWER2(s[4] - Real(2.0) * s[3] + s[2]) +
                Real(0.25) *
                    POWER2(s[4] - Real(4.0) * s[3] + Real(3.0) * s[2]);
      candidate[2] =
          (Real(15.0) * s[2] - Real(10.0) * s[3] + Real(3.0) * s[4]) /
          Real(8.0);
    }
    if (v1) {
      beta[1] = Real(13.0 / 12.0) *
                    POWER2(s[3] - Real(2.0) * s[2] + s[1]) +
                Real(0.25) * POWER2(s[3] - s[1]);
      candidate[1] =
          (-s[3] + Real(6.0) * s[2] + Real(3.0) * s[1]) /
          Real(8.0);
    }
    if (v0) {
      beta[0] = Real(13.0 / 12.0) *
                    POWER2(s[2] - Real(2.0) * s[1] + s[0]) +
                Real(0.25) *
                    POWER2(Real(3.0) * s[2] - Real(4.0) * s[1] + s[0]);
      candidate[0] =
          (Real(3.0) * s[2] + Real(6.0) * s[1] - s[0]) /
          Real(8.0);
    }

    Real tau = Real(0.0);
    if (v0 && v1 && v2) {
      tau = std::abs(beta[2] - beta[0]);
    } else if (v0 && v1) {
      tau = std::abs(beta[1] - beta[0]);
    } else if (v1 && v2) {
      tau = std::abs(beta[2] - beta[1]);
    }
    if (v2) {
      beta[2] =
          nonlinear_factor(tau, beta[2], legacy_epsilon(), false) *
          optimal_weight[2];
    }
    if (v1) {
      beta[1] =
          nonlinear_factor(tau, beta[1], legacy_epsilon(), false) *
          optimal_weight[1];
    }
    if (v0) {
      beta[0] =
          nonlinear_factor(tau, beta[0], legacy_epsilon(), false) *
          optimal_weight[0];
    }
    const Real denominator = beta[0] + beta[1] + beta[2];
    if (!(denominator > std::numeric_limits<Real>::min())) {
      if (v1) return candidate[1];
      if (v2) return candidate[2];
      if (v0) return candidate[0];
      return s[2];
    }
    return (beta[0] * candidate[0] + beta[1] * candidate[1] +
            beta[2] * candidate[2]) /
           denominator;
  }

#ifdef CNS_IBM_VALIDATION
  // Apply WENO-Z weights computed from a frozen reference stencil to candidate
  // polynomials evaluated from a second stencil. This is an operator audit,
  // not a production reconstruction: it separates ghost-state truncation from
  // the nonlinear response of the smoothness indicators.
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  recon_masked_frozen(
      const amrex::Real value[5], const amrex::Real indicator[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight) noexcept {
    using amrex::Real;
    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);
    Real candidate[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real beta[3] = {Real(0.0), Real(0.0), Real(0.0)};
    if (v2) {
      beta[2] = Real(13.0 / 12.0) *
                    POWER2(indicator[4] - Real(2.0) * indicator[3] +
                           indicator[2]) +
                Real(0.25) *
                    POWER2(indicator[4] - Real(4.0) * indicator[3] +
                           Real(3.0) * indicator[2]);
      candidate[2] = Real(11.0) * value[2] - Real(7.0) * value[3] +
                     Real(2.0) * value[4];
    }
    if (v1) {
      beta[1] = Real(13.0 / 12.0) *
                    POWER2(indicator[3] - Real(2.0) * indicator[2] +
                           indicator[1]) +
                Real(0.25) * POWER2(indicator[3] - indicator[1]);
      candidate[1] = -value[3] + Real(5.0) * value[2] +
                     Real(2.0) * value[1];
    }
    if (v0) {
      beta[0] = Real(13.0 / 12.0) *
                    POWER2(indicator[2] - Real(2.0) * indicator[1] +
                           indicator[0]) +
                Real(0.25) *
                    POWER2(Real(3.0) * indicator[2] -
                           Real(4.0) * indicator[1] + indicator[0]);
      candidate[0] = Real(2.0) * value[2] + Real(5.0) * value[1] - value[0];
    }
    Real tau = Real(0.0);
    if (v0 && v1 && v2) {
      tau = std::abs(beta[2] - beta[0]);
    } else if (v0 && v1) {
      tau = std::abs(beta[1] - beta[0]);
    } else if (v1 && v2) {
      tau = std::abs(beta[2] - beta[1]);
    }
    if (v2) {
      beta[2] =
          nonlinear_factor(tau, beta[2], legacy_epsilon(), false) *
          optimal_weight[2];
    }
    if (v1) {
      beta[1] =
          nonlinear_factor(tau, beta[1], legacy_epsilon(), false) *
          optimal_weight[1];
    }
    if (v0) {
      beta[0] =
          nonlinear_factor(tau, beta[0], legacy_epsilon(), false) *
          optimal_weight[0];
    }
    const Real denominator = beta[0] + beta[1] + beta[2];
    if (!(denominator > std::numeric_limits<Real>::min())) {
      if (v1) return candidate[1] / Real(6.0);
      if (v0) return candidate[0] / Real(6.0);
      if (v2) return candidate[2] / Real(6.0);
      return candidate[1] / Real(6.0);
    }
    return (beta[0] * candidate[0] + beta[1] * candidate[1] +
            beta[2] * candidate[2]) /
           (Real(6.0) * denominator);
  }
#endif

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon(
    const amrex::Real s[5], const int gl, const int gr) noexcept {
    using amrex::Real;
    amrex::GpuArray<Real, 3> optimal_weight{Real(3.0), Real(6.0), Real(1.0)};

    // Near-wall masking via linear (optimal) weights: disable stencils crossing IB
    if (gr == 1) {
      optimal_weight[0] = Real(1.0); // keep S0
      optimal_weight[1] = Real(1.0); // keep S1
      optimal_weight[2] = Real(0.0); // drop S2 (uses i+2)
    }
    if (gr == 0) {
      optimal_weight[0] = Real(1.0); // keep S0
      optimal_weight[1] = Real(0.0); // drop S1
      optimal_weight[2] = Real(0.0); // drop S2
    }
    if (gl == 0) {
      optimal_weight[0] = Real(0.0); // drop S0 (uses i-2)
      optimal_weight[1] = Real(3.0); // favor S1
      optimal_weight[2] = Real(1.0); // keep S2
    }

    return recon_masked(s, optimal_weight);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_scaled(
      const amrex::Real s[5], const int gl, const int gr,
      const amrex::Real epsilon) noexcept {
    using amrex::Real;
    amrex::GpuArray<Real, 3> optimal_weight{
        Real(3.0), Real(6.0), Real(1.0)};

    if (gr == 1) {
      optimal_weight[0] = Real(1.0);
      optimal_weight[1] = Real(1.0);
      optimal_weight[2] = Real(0.0);
    }
    if (gr == 0) {
      optimal_weight[0] = Real(1.0);
      optimal_weight[1] = Real(0.0);
      optimal_weight[2] = Real(0.0);
    }
    if (gl == 0) {
      optimal_weight[0] = Real(0.0);
      optimal_weight[1] = Real(3.0);
      optimal_weight[2] = Real(1.0);
    }

    return recon_masked_scaled(s, optimal_weight, epsilon);
  }
};  // struct WenoZ5

struct Teno5 : public WenoZ5 {
  static constexpr const char* scheme_name =
      "characteristic_llf_teno5_old";
  /**
   * \brief TENO-5 reconstruction of characteristic split fluxes for the
   * conservative finite-difference operator.
   * Ref https://doi.org/10.1016/j.jcp.2015.10.037.
   * \param[in] s Stencil values at cells [i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon(
      const amrex::Real s[5], const int gl, const int gr) noexcept {
    return recon_with_cutoff(s, gl, gr, amrex::Real(1.0e-3));
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  recon_with_cutoff(const amrex::Real s[5], const int gl, const int gr,
                    const amrex::Real cutoff) noexcept {
    using amrex::Real;
    amrex::GpuArray<Real, 3> optimal_weight{Real(3.0), Real(6.0), Real(1.0)};

    // Near-wall masking via linear (optimal) weights: disable stencils crossing IB
    if (gr == 1) {
      optimal_weight[0] = Real(1.0); // keep S0
      optimal_weight[1] = Real(1.0); // keep S1
      optimal_weight[2] = Real(0.0); // drop S2 (uses i+2)
    }
    if (gr == 0) {
      optimal_weight[0] = Real(1.0); // keep S0
      optimal_weight[1] = Real(0.0); // drop S1
      optimal_weight[2] = Real(0.0); // drop S2
    }
    if (gl == 0) {
      optimal_weight[0] = Real(0.0); // drop S0 (uses i-2)
      optimal_weight[1] = Real(3.0); // favor S1
      optimal_weight[2] = Real(1.0); // keep S2
    }

    return recon_masked_with_cutoff(s, optimal_weight, cutoff);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_masked(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight) noexcept {
    return recon_masked_with_cutoff(
        s, optimal_weight, amrex::Real(1.0e-3));
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  recon_masked_with_cutoff(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight,
      const amrex::Real cutoff) noexcept {
    using amrex::Real;

    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);

    Real vr[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real beta[3] = {Real(0.0), Real(0.0), Real(0.0)};

    if (v2) {
      beta[2] = Real(13. / 12.) * POWER2(s[4] - 2.0 * s[3] + s[2]) +
                Real(0.25) * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
      vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
    }
    if (v1) {
      beta[1] = Real(13. / 12.) * POWER2(s[3] - 2.0 * s[2] + s[1]) +
                Real(0.25) * POWER2(s[3] - s[1]);
      vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
    }
    if (v0) {
      beta[0] = Real(13. / 12.) * POWER2(s[2] - 2.0 * s[1] + s[0]) +
                Real(0.25) * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);
      vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];
    }

    Real tau;
    if (v0 && v1 && v2) {
      tau = std::abs(std::abs(beta[2] - beta[0]) -
                     (beta[2] + Real(4.0) * beta[1] + beta[0]) /
                         Real(6.0));
    } else if (v0 && v1) {
      tau = std::abs(beta[1] - beta[0]);
    } else if (v1 && v2) {
      tau = std::abs(beta[2] - beta[1]);
    } else {
      tau = Real(0.0);
    }

    Real gamma[3] = {Real(0.0), Real(0.0), Real(0.0)};
    if (v2) gamma[2] = POWER6(1.0 + tau / (eps + beta[2]));
    if (v1) gamma[1] = POWER6(1.0 + tau / (eps + beta[1]));
    if (v0) gamma[0] = POWER6(1.0 + tau / (eps + beta[0]));
    const Real gamma_sum = gamma[2] + gamma[1] + gamma[0];
    if (!(gamma_sum > std::numeric_limits<Real>::min())) {
      return WenoZ5::recon_masked(s, optimal_weight);
    }
    const Real inv_gamma_sum = Real(1.0) / gamma_sum;
    beta[2] = v2 && gamma[2] * inv_gamma_sum >= cutoff ? optimal_weight[2] : Real(0.0);
    beta[1] = v1 && gamma[1] * inv_gamma_sum >= cutoff ? optimal_weight[1] : Real(0.0);
    beta[0] = v0 && gamma[0] * inv_gamma_sum >= cutoff ? optimal_weight[0] : Real(0.0);
    const Real denom = beta[2] + beta[1] + beta[0];
    if (!(denom > std::numeric_limits<Real>::min())) {
      return WenoZ5::recon_masked(s, optimal_weight);
    }
    const Real inv_denom = 1.0 / denom;

    return inv_denom / Real(6.0) *
           (beta[2] * vr[2] + beta[1] * vr[1] + beta[0] * vr[0]);
  }

};  // struct Teno5

struct Teno6 {
  static constexpr int ng = 3;
  static constexpr int ncand = 4;
  static constexpr const char* scheme_name =
      "characteristic_llf_teno6_old";

  /**
   * @brief From a 2*ng x NCONS array, extract the left stencil.
   * @param[in] n  Index of the component.
   * @param[in] fp Flux array.
   * @param[out] s Stencil array [i+2, i+1, i, i-1, i-2, i-3].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void left_stencil(
      int n, amrex::Real const fp[2 * ng][NCONS], amrex::Real s[6]) noexcept {
    for (int m = 0; m < 2 * ng; ++m) s[m] = fp[2 * ng - 1 - m][n];
  }

  /**
   * @brief From a 2*ng x NCONS array, extract the left stencil.
   * @param[in] n  Index of the component.
   * @param[in] fm Flux array.
   * @param[out] s Stencil array [i-3, i-2, i-1, i, i+1, i+2].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void right_stencil(
      int n, amrex::Real const fm[2 * ng][NCONS], amrex::Real s[6]) noexcept {
    for (int m = 0; m < 2 * ng; ++m) s[m] = fm[m][n];
  }

  /**
   * \brief TENO-6 reconstruction of characteristic split fluxes for the
   * conservative finite-difference operator.
   * Ref https://doi.org/10.1016/j.jcp.2015.10.037.
   * \param[in] s Stencil values at cells [i-3, i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_masked(
      const amrex::Real s[6],
      const amrex::GpuArray<amrex::Real, 4>& optimal_weight) noexcept {
    using amrex::Real;

    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    constexpr Real cutoff = 1e-5;
    Real vr[4], beta[4], tmp;

    constexpr Real w13o12 = 13.0 / 12.0;
    constexpr Real w1o4 = 0.25;
    beta[3] =
        Real(1. / 36.) *
            POWER2(-11.0 * s[3] + 18.0 * s[2] - 9.0 * s[1] + 2.0 * s[0]) +
        w13o12 * POWER2(2.0 * s[3] - 5.0 * s[2] + 4.0 * s[1] - s[0]) +
        Real(781. / 720.) * POWER2(-s[3] + 3.0 * s[2] - 3.0 * s[1] + s[0]);
    beta[2] = w13o12 * POWER2(s[4] - 2.0 * s[3] + s[2]) +
              w1o4 * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
    beta[1] =
        w13o12 * POWER2(s[3] - 2.0 * s[2] + s[1]) + w1o4 * POWER2(s[3] - s[1]);
    beta[0] = w13o12 * POWER2(s[2] - 2.0 * s[1] + s[0]) +
              w1o4 * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);

    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);
    const bool v3 = optimal_weight[3] > Real(0.0);
    if (v0 && v1 && v2 && v3) {
      tmp = std::abs(beta[3] - (beta[2] + Real(4.0) * beta[1] + beta[0]) / Real(6.0));
    } else {
      Real beta_min = std::numeric_limits<Real>::max();
      Real beta_max = Real(0.0);
      int nvalid = 0;
      if (v0) { beta_min = amrex::min(beta_min, beta[0]); beta_max = amrex::max(beta_max, beta[0]); ++nvalid; }
      if (v1) { beta_min = amrex::min(beta_min, beta[1]); beta_max = amrex::max(beta_max, beta[1]); ++nvalid; }
      if (v2) { beta_min = amrex::min(beta_min, beta[2]); beta_max = amrex::max(beta_max, beta[2]); ++nvalid; }
      if (v3) { beta_min = amrex::min(beta_min, beta[3]); beta_max = amrex::max(beta_max, beta[3]); ++nvalid; }
      tmp = (nvalid > 1) ? (beta_max - beta_min) : Real(0.0);
    }
    beta[3] = POWER6(1.0 + tmp / (eps + beta[3]));
    beta[2] = POWER6(1.0 + tmp / (eps + beta[2]));
    beta[1] = POWER6(1.0 + tmp / (eps + beta[1]));
    beta[0] = POWER6(1.0 + tmp / (eps + beta[0]));
    tmp = 1.0 / (beta[3] + beta[2] + beta[1] + beta[0]);
    beta[3] = beta[3] * tmp < cutoff ? Real(0.0) : optimal_weight[3];
    beta[2] = beta[2] * tmp < cutoff ? Real(0.0) : optimal_weight[2];
    beta[1] = beta[1] * tmp < cutoff ? Real(0.0) : optimal_weight[1];
    beta[0] = beta[0] * tmp < cutoff ? Real(0.0) : optimal_weight[0];

    vr[3] = 0.5 * (s[0] - 5.0 * s[1] + 13.0 * s[2] + 3.0 * s[3]);
    vr[2] = 11.0 * s[3] - 7.0 * s[4] + 2.0 * s[5];
    vr[1] = -s[4] + 5.0 * s[3] + 2.0 * s[2];
    vr[0] = 2.0 * s[3] + 5.0 * s[2] - s[1];

    const Real denom = beta[3] + beta[2] + beta[1] + beta[0];
    if (!(denom > std::numeric_limits<Real>::min())) {
      if (v1) return vr[1] / Real(6.0);
      if (v0) return vr[0] / Real(6.0);
      if (v2) return vr[2] / Real(6.0);
      if (v3) return vr[3] / Real(6.0);
      return vr[1] / Real(6.0);
    }
    tmp = 1.0 / denom;

    return tmp / Real(6.0) *
           (beta[3] * vr[3] + beta[2] * vr[2] + beta[1] * vr[1] +
            beta[0] * vr[0]);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon(
      const amrex::Real s[6], int /*gl*/, int /*gr*/) noexcept {
    using amrex::Real;
    amrex::GpuArray<Real, 4> optimal_weight{Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
    return recon_masked(s, optimal_weight);
  }
};  // struct Teno6
};  // namespace OldReconScheme

// TODO: this class name is not very accurate. llf? high_order_fd?
template <typename Scheme, typename cls_t>
class weno_old_t {
  static constexpr int ng = Scheme::ng;  // number of ghost cells on each side
  static_assert(ng <= cls_t::NGHOST);    // ensure enough ghost cells

 public:
  static constexpr const char* scheme_name = Scheme::scheme_name;
  // Only the WENO/TENO eflux implements the RZ pressure-split (it removes the
  // radial pressure flux to pair with compute_rhs's -dp/dr form). Other flux
  // schemes do NOT, so enabling cns.rz_pressure_split with them would
  // double-count pressure. compute_rhs guards on this trait.
  static constexpr bool rz_pressure_split_capable = true;
  // Pressure-only near-axis well balancing reconstructs a bounded physical
  // pressure at radial faces.  It deliberately does not reuse the conservative
  // momentum flux, which also contains advection and LLF dissipation.
  static constexpr bool rz_pressure_face_capable =
      Scheme::ng == 3 && Scheme::ncand == 3;
  static constexpr bool rz_annular_pressure_consistency_capable = true;

  template <typename PrimArray>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  rz_annular_pressure_from_cell_average(
      const amrex::IntVect& cell, const PrimArray& prims,
      const amrex::Real radial, const amrex::Real dr) noexcept
  {
    using amrex::Real;
    const Real radial2 = radial * radial;
    const Real dr2 = dr * dr;
    const Real midpoint_ur =
        prims(cell, cls_t::QU) * radial2 /
        (radial2 + dr2 / Real(12.0));
    const Real radial_slope = midpoint_ur / radial;
    const Real radial_variance =
        dr2 / Real(12.0) - dr2 * dr2 /
                                  (Real(144.0) * radial2);
    const Real unresolved_radial_ke =
        Real(0.5) * prims(cell, cls_t::QRHO) * radial_slope *
        radial_slope * amrex::max(radial_variance, Real(0.0));
    return prims(cell, cls_t::QPRES) -
           (prims(cell, cls_t::QG) - Real(1.0)) *
               unresolved_radial_ke;
  }

  template <typename PrimArray>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  rz_pressure_face(const amrex::IntVect& face, const int dir,
                   const PrimArray& prims,
                   const int pressure_component) noexcept
  {
    static_assert(Scheme::ng == 3 && Scheme::ncand == 3,
                  "RZ pressure-face reconstruction requires WENO-Z5/TENO5");
    using amrex::Real;
    const amrex::IntVect ivd =
        amrex::IntVect::TheDimensionVector(dir);

    // Reverse the left-state stencil so Scheme::recon evaluates its right
    // face.  The ordinary right-state ordering evaluates its left face.
    Real left_stencil[5];
    Real right_stencil[5];
    for (int m = 0; m < 5; ++m) {
      left_stencil[m] =
          prims(face + (1 - m) * ivd, pressure_component);
      right_stencil[m] =
          prims(face + (m - 2) * ivd, pressure_component);
    }
    const Real left_trace = Scheme::recon(left_stencil, ng, ng);
    const Real right_trace = Scheme::recon(right_stencil, ng, ng);

    const Real p_left = prims(face - ivd, pressure_component);
    const Real p_right = prims(face, pressure_component);
    const Real fallback = Real(0.5) * (p_left + p_right);
    const Real reconstructed = Real(0.5) * (left_trace + right_trace);
    if (!std::isfinite(reconstructed)) return fallback;

    // A geometric source must not introduce a new pressure extremum.  At r=0
    // reflected adjacent states are identical, so this also returns p(cell 0)
    // exactly on the symmetry face.
    const Real pressure_min = amrex::min(p_left, p_right);
    const Real pressure_max = amrex::max(p_left, p_right);
    return amrex::min(amrex::max(reconstructed, pressure_min), pressure_max);
  }

  AMREX_GPU_HOST_DEVICE
  weno_old_t() {}

  AMREX_GPU_HOST_DEVICE
  ~weno_old_t() {}

#if (AMREX_USE_GPIBM || CNS_USE_EB )
  void inline eflux_ibm(const amrex::Geometry& geom, const amrex::MFIter& mfi,
                        const amrex::Array4<const amrex::Real>& prims_in,
                        std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                        const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
                        //const amrex::Array4<const bool>& ibMarkers)
                        const amrex::Array4<uint8_t>& ibMarkers
#ifdef CNS_IBM_VALIDATION
                        , const std::array<FArrayBox*, AMREX_SPACEDIM>*
                              split_positive_audit = nullptr
                        , const std::array<FArrayBox*, AMREX_SPACEDIM>*
                              split_negative_audit = nullptr
                        , const amrex::Array4<const amrex::Real>*
                              frozen_coefficient_prims = nullptr
                        , const amrex::Array4<amrex::Real>*
                              rz_shock_face_audit = nullptr
#endif
                        )
  {
    const amrex::Box bx_region = mfi.tilebox();
    const amrex::Box skip_cells;  // invalid: no face skipping in the IBM path
#else
  void inline eflux(const amrex::Geometry& geom, const amrex::MFIter& mfi,
                    const amrex::Array4<const amrex::Real>& prims_in,
                    std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const amrex::Array4<amrex::Real>& rhs, const cls_t* cls)
  {
    eflux(geom, mfi.tilebox(), prims_in, flxt, rhs, cls);
  }

  // Region-parameterized overload (comm/comp overlap): computes fluxes only
  // on the faces of surroundingNodes(bx_region, dir). Faces whose BOTH
  // adjacent cells lie inside the optional 'skip_cells' box are skipped
  // (used by the overlap shell pass so interior faces already handled in
  // Pass 1 are not recomputed). Per-face numerics are identical to the
  // MFIter version, so splitting a tilebox into interior + shell regions
  // yields bitwise-identical fluxes at every face (seam faces are
  // recomputed from the same prims).
  void inline eflux(const amrex::Geometry& geom, const amrex::Box& bx_region,
                    const amrex::Array4<const amrex::Real>& prims_in,
                    std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
                    const amrex::Box& skip_cells = amrex::Box())
  {
#endif
    using amrex::Array4, amrex::Box, amrex::Dim3, amrex::IntVect, amrex::Real;

#if AMREX_USE_GPIBM
    static_assert(
        Scheme::ncand == 3,
        "IBM-safe stencil masking is implemented only for WenoZ5/Teno5; "
        "Teno6 must not be used with GPIBM until its four-candidate masking "
        "is made solid-state safe.");
#endif

    const Box& bx = bx_region;
    const Box skipbox = skip_cells;
    const bool skip_ok = skipbox.ok();
#if defined(CNS_IBM_VALIDATION) && (AMREX_USE_GPIBM || CNS_USE_EB)
    const bool capture_split_flux =
        split_positive_audit != nullptr && split_negative_audit != nullptr;
    const bool freeze_coefficients = frozen_coefficient_prims != nullptr;
    Array4<const Real> coefficient_prims;
    if (freeze_coefficients) coefficient_prims = *frozen_coefficient_prims;
    if ((split_positive_audit == nullptr) !=
        (split_negative_audit == nullptr)) {
      amrex::Abort("IBM WENO split-flux audit requires both F+ and F- arrays");
    }
#else
    constexpr bool freeze_coefficients = false;
#ifdef CNS_IBM_VALIDATION
    Array4<const Real> coefficient_prims;
#endif
#endif
#if defined(CNS_IBM_VALIDATION) && AMREX_USE_GPIBM
    const bool capture_rz_shock_face_audit =
        rz_shock_face_audit != nullptr;
    Array4<Real> rz_shock_face_diagnostic;
    if (capture_rz_shock_face_audit) {
      rz_shock_face_diagnostic = *rz_shock_face_audit;
    }
#else
    constexpr bool capture_rz_shock_face_audit = false;
#endif
    // Experimental RZ pressure split.  This must be paired with the matching
    // compute_rhs.cpp source form, and is only implemented for WENO/TENO.
    static const int s_rz_pressure_split = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_pressure_split", v);
      return v;
    }();
#if !AMREX_USE_GPIBM && !CNS_USE_EB
    // Opt-in scale-aware WENO-Z regularisation for the standard Cartesian
    // all-fluid LLF path. The scale is formed after characteristic projection
    // and is shared by every component and both split branches at a face.
    static const Real s_llf_weno_epsilon_relative = [] {
      Real value = Real(0.0);
      amrex::ParmParse pp("cns");
      pp.query("llf_weno_epsilon_relative", value);
      if (!std::isfinite(value) || value < Real(0.0) ||
          value > Real(1.0)) {
        amrex::Abort(
            "cns.llf_weno_epsilon_relative must be finite and in [0,1]");
      }
      return value;
    }();
    static const Real s_llf_teno_cutoff = [] {
      Real value = Real(1.0e-3);
      amrex::ParmParse pp("cns");
      pp.query("llf_teno_cutoff", value);
      if (!std::isfinite(value) || value <= Real(0.0) ||
          value >= Real(1.0)) {
        amrex::Abort("cns.llf_teno_cutoff must be finite and in (0,1)");
      }
      return value;
    }();
    if (geom.IsRZ() && s_llf_weno_epsilon_relative > Real(0.0)) {
      amrex::Abort(
          "cns.llf_weno_epsilon_relative is currently qualified only for "
          "Cartesian all-fluid reconstruction");
    }
    const Real llf_weno_epsilon_relative =
        s_llf_weno_epsilon_relative;
    const Real llf_teno_cutoff = s_llf_teno_cutoff;
#else
    constexpr Real llf_weno_epsilon_relative = Real(0.0);
    constexpr Real llf_teno_cutoff = Real(1.0e-3);
#endif

    // Opt-in composite R-Z Euler paths.  The original hybrid is retained for
    // exact reproduction of the frozen H2 campaign.  Version 2 keeps the
    // proven annular-average correction everywhere, but restricts point-flux
    // reconstruction to locally smooth radial stencils and applies the
    // grid-aligned-shock LLF only to transverse (radial) faces.  This avoids
    // sending a curved positive-flow bow shock through the experimental point
    // reconstruction and avoids replacing the shock-normal Mach-disk flux by
    // first order.  Every default remains unchanged.
    static const int s_rz_euler_consistent_shock_hybrid = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_consistent_shock_hybrid", v);
      return v;
    }();
    static const int s_rz_euler_consistent_shock_hybrid_v2 = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_consistent_shock_hybrid_v2", v);
      return v;
    }();
    const bool rz_euler_consistent_shock_hybrid =
        geom.IsRZ() && (s_rz_euler_consistent_shock_hybrid != 0);
    const bool rz_euler_consistent_shock_hybrid_v2 =
        geom.IsRZ() && (s_rz_euler_consistent_shock_hybrid_v2 != 0);
    if (rz_euler_consistent_shock_hybrid &&
        rz_euler_consistent_shock_hybrid_v2) {
      amrex::Abort(
          "select only one of cns.rz_euler_consistent_shock_hybrid and "
          "cns.rz_euler_consistent_shock_hybrid_v2");
    }
    const bool rz_euler_any_consistent_shock_hybrid =
        rz_euler_consistent_shock_hybrid ||
        rz_euler_consistent_shock_hybrid_v2;

    // In an R-Z finite-volume cell, the stored odd radial velocity is an
    // annular average.  Treating it as a Cartesian midpoint value gives an
    // O(1) continuity defect at the axis even for u_r=a*r.  Recover the
    // leading regular midpoint value for radial WENO stencils only.  This is
    // opt-in while the production SRP comparison is being audited.
    static const int s_rz_euler_annular_average = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_annular_average", v);
      return v;
    }();
    const bool rz_euler_annular_average =
        geom.IsRZ() && ((s_rz_euler_annular_average != 0) ||
                       rz_euler_any_consistent_shock_hybrid);
    static const int s_rz_euler_point_flux = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_point_flux", v);
      return v;
    }();
    const bool rz_euler_point_flux =
        geom.IsRZ() && ((s_rz_euler_point_flux != 0) ||
                       rz_euler_any_consistent_shock_hybrid);
    static const int s_rz_euler_annular_pressure_consistency = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_annular_pressure_consistency", v);
      return v;
    }();
    const bool rz_euler_annular_pressure_consistency =
        geom.IsRZ() &&
        (s_rz_euler_annular_pressure_consistency != 0);
    static const int s_rz_euler_cell_average_deconvolution = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_cell_average_deconvolution", v);
      return v;
    }();
    const bool rz_euler_cell_average_deconvolution =
        geom.IsRZ() && (s_rz_euler_cell_average_deconvolution != 0);
    static const int s_rz_euler_axial_annular_flux = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_axial_annular_flux", v);
      return v;
    }();
    const bool rz_euler_axial_annular_flux =
        geom.IsRZ() && (s_rz_euler_axial_annular_flux != 0);
    static const int s_rz_gp_annular_bic = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_gp_annular_bic", v);
      return v;
    }();
    const bool rz_gp_annular_bic =
        geom.IsRZ() && (s_rz_gp_annular_bic != 0);
    if (rz_euler_point_flux && !rz_euler_annular_average) {
      amrex::Abort(
          "cns.rz_euler_point_flux requires "
          "cns.rz_euler_annular_average=1");
    }
    if (rz_euler_annular_pressure_consistency &&
        (!rz_euler_point_flux || rz_euler_any_consistent_shock_hybrid)) {
      amrex::Abort(
          "cns.rz_euler_annular_pressure_consistency requires an ungated "
          "R-Z annular point-flux path");
    }
    if (rz_euler_cell_average_deconvolution &&
        (!rz_euler_annular_average || !rz_euler_point_flux ||
         !rz_euler_annular_pressure_consistency ||
         rz_euler_any_consistent_shock_hybrid)) {
      amrex::Abort(
          "cns.rz_euler_cell_average_deconvolution requires the ungated "
          "annular-average, point-flux, and pressure-consistency R-Z path");
    }
    if (rz_euler_axial_annular_flux &&
        !rz_euler_cell_average_deconvolution) {
      amrex::Abort(
          "cns.rz_euler_axial_annular_flux requires "
          "cns.rz_euler_cell_average_deconvolution=1");
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
          "R-Z GP-IBM requires centre-state recovery and annular BI-CWLS "
          "to be enabled together");
    }
#else
    if (rz_gp_annular_bic) {
      amrex::Abort("cns.rz_gp_annular_bic requires AMREX_USE_GPIBM");
    }
#endif
    // Conservative transverse coupling for an axis-spanning, grid-aligned
    // shock. A strong axial pressure jump switches only nearby radial faces
    // to first-order LLF; identical radial states retain the exact same flux.
    static const int s_rz_grid_aligned_shock_filter = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_grid_aligned_shock_filter", v);
      return v;
    }();
    static const Real s_rz_grid_aligned_shock_threshold = [] {
      Real v = Real(0.2);
      amrex::ParmParse pp("cns");
      pp.query("rz_grid_aligned_shock_threshold", v);
      return v;
    }();
    static const Real s_rz_grid_aligned_shock_ratio = [] {
      Real v = Real(2.0);
      amrex::ParmParse pp("cns");
      pp.query("rz_grid_aligned_shock_ratio", v);
      return v;
    }();
    static const int s_rz_grid_aligned_shock_flow_sign = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_grid_aligned_shock_flow_sign", v);
      return v;
    }();
    const bool rz_grid_aligned_shock_filter =
        geom.IsRZ() && ((s_rz_grid_aligned_shock_filter != 0) ||
                       rz_euler_any_consistent_shock_hybrid);
    const Real rz_grid_aligned_shock_threshold =
        amrex::max(Real(0.0), s_rz_grid_aligned_shock_threshold);
    const Real rz_grid_aligned_shock_ratio =
        amrex::max(Real(0.0), s_rz_grid_aligned_shock_ratio);
    const int rz_grid_aligned_shock_flow_sign =
        s_rz_grid_aligned_shock_flow_sign < 0
            ? -1
            : (s_rz_grid_aligned_shock_flow_sign > 0 ? 1 : 0);
    // Opt-in H-correction-type transverse stabilization.  Unlike the legacy
    // grid-aligned filter, this retains the complete high-order WENO flux and
    // adds only a conservative Rusanov dissipation term on radial faces:
    //   dF_H = -0.5 * weight_H * alpha_r * (U_R - U_L).
    // The axial H-stencil supplies weight_H, while identical radial states
    // make the correction exactly zero.
    static const int s_rz_euler_h_correction = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_h_correction", v);
      return v;
    }();
    static const Real s_rz_euler_h_correction_strength = [] {
      Real v = Real(0.125);
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_h_correction_strength", v);
      return v;
    }();
    const bool rz_euler_h_correction =
        geom.IsRZ() && (s_rz_euler_h_correction != 0);
    const Real rz_euler_h_correction_strength =
        s_rz_euler_h_correction_strength;
#if (AMREX_SPACEDIM != 2)
    if (s_rz_euler_h_correction != 0) {
      amrex::Abort("cns.rz_euler_h_correction requires a 2-D R-Z build");
    }
#endif
    if (s_rz_euler_h_correction != 0 && !geom.IsRZ()) {
      amrex::Abort("cns.rz_euler_h_correction requires R-Z geometry");
    }
    if (rz_euler_h_correction &&
        (!(rz_euler_h_correction_strength > Real(0.0)) ||
         !(rz_euler_h_correction_strength <= Real(1.0)))) {
      amrex::Abort(
          "cns.rz_euler_h_correction_strength must lie in (0,1]");
    }
    if (rz_euler_h_correction &&
        (!(s_rz_grid_aligned_shock_threshold >= Real(0.0)) ||
         !(s_rz_grid_aligned_shock_threshold < Real(1.0)))) {
      amrex::Abort(
          "H-correction requires cns.rz_grid_aligned_shock_threshold in "
          "[0,1)");
    }
    if (rz_euler_h_correction &&
        rz_grid_aligned_shock_flow_sign == 0) {
      amrex::Abort(
          "H-correction requires an explicit nonzero "
          "cns.rz_grid_aligned_shock_flow_sign");
    }
    if (rz_euler_h_correction && rz_grid_aligned_shock_filter) {
      amrex::Abort(
          "select either cns.rz_euler_h_correction or the legacy full-LLF "
          "grid-aligned shock filter/hybrid, not both");
    }
    static const int s_rz_euler_point_flux_shock_gate = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_euler_point_flux_shock_gate", v);
      return v;
    }();
    const bool rz_euler_point_flux_shock_gate =
        geom.IsRZ() && ((s_rz_euler_point_flux_shock_gate != 0) ||
                       rz_euler_consistent_shock_hybrid_v2);
    static const int s_rz_grid_aligned_shock_transverse_only = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_grid_aligned_shock_transverse_only", v);
      return v;
    }();
    const bool rz_grid_aligned_shock_transverse_only =
        geom.IsRZ() &&
        ((s_rz_grid_aligned_shock_transverse_only != 0) ||
         rz_euler_consistent_shock_hybrid_v2);
    if (rz_euler_any_consistent_shock_hybrid &&
        rz_grid_aligned_shock_flow_sign == 0) {
      amrex::Abort(
          "the R-Z consistent-shock hybrid requires an explicit nonzero "
          "cns.rz_grid_aligned_shock_flow_sign; unrestricted filtering is "
          "known to smear shocks moving in the opposite direction");
    }

    std::unique_ptr<FArrayBox> rz_annular_prims_owner;
    Array4<const Real> rz_annular_prims = prims_in;
    if (rz_euler_annular_average) {
      const Box prim_box = amrex::grow(bx, ng);
      rz_annular_prims_owner = std::make_unique<FArrayBox>(
          prim_box, cls_t::NPRIM, amrex::The_Async_Arena());
      const Array4<Real> corrected = rz_annular_prims_owner->array();
      const auto prob_lo = geom.ProbLoArray();
      const auto cell_size = geom.CellSizeArray();
      if (rz_euler_cell_average_deconvolution) {
        const int radial_data_lo = prim_box.smallEnd(0);
        const int radial_data_hi = prim_box.bigEnd(0);
        const Real radial_origin_over_dr = prob_lo[0] / cell_size[0];
        ParallelFor(
            prim_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real point_state[cls_t::NCONS];
              Real point_primitives[cls_t::NPRIM];
              for (int n = 0; n < cls_t::NPRIM; ++n) {
                point_primitives[n] = prims_in(i, j, k, n);
              }
              cerisse::rz_fv::annular_average_to_radial_center<
                  cls_t::NCONS>(
                  i, j, k, radial_data_lo, radial_data_hi,
                  radial_origin_over_dr, rhs, point_state);
              cls->cons2prims_point(point_state, point_primitives);
              for (int n = 0; n < cls_t::NPRIM; ++n) {
                corrected(i, j, k, n) = point_primitives[n];
              }
            });
      } else {
        ParallelFor(
            prim_box, cls_t::NPRIM,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
              Real value = prims_in(i, j, k, n);
              const Real dr = cell_size[0];
              const Real r =
                  prob_lo[0] + (Real(i) + Real(0.5)) * dr;
              if (n == cls_t::QU) {
                const Real r2 = r * r;
                value *= r2 / (r2 + dr * dr / Real(12.0));
              }
              if (rz_euler_annular_pressure_consistency &&
                  (n == cls_t::QPRES || n == cls_t::QT ||
                   n == cls_t::QC || n == cls_t::QEINT)) {
                const IntVect cell(AMREX_D_DECL(i, j, k));
                const Real corrected_pressure =
                    rz_annular_pressure_from_cell_average(
                        cell, prims_in, r, dr);
                const Real pressure_ratio =
                    corrected_pressure / prims_in(cell, cls_t::QPRES);
                if (n == cls_t::QPRES) {
                  value = corrected_pressure;
                } else if (n == cls_t::QC) {
                  value *= std::sqrt(pressure_ratio);
                } else {
                  value *= pressure_ratio;
                }
              }
              corrected(i, j, k, n) = value;
            });
      }
      rz_annular_prims = rz_annular_prims_owner->const_array();
    }

    // The radial point-flux path above needs meridional centre states.  The
    // axial finite-volume flux has different data semantics: its radial
    // momentum support must remain the r-weighted annular cell average.  In
    // the first axis-adjacent ring, replacing that average by its centre value
    // changes a regular odd field m_r=a*r by the fixed factor 3/4 and leaves
    // an O(dr) radial-momentum residual.  Reconstruct Q(Ubar_RZ) from the
    // conservative storage for axial WENO while retaining the centre-state
    // field for radial faces and BI-CWLS.
    std::unique_ptr<FArrayBox> rz_axial_annular_prims_owner;
    Array4<const Real> rz_axial_annular_prims = prims_in;
    if (rz_euler_axial_annular_flux) {
      const Box prim_box = amrex::grow(bx, ng);
      rz_axial_annular_prims_owner = std::make_unique<FArrayBox>(
          prim_box, cls_t::NPRIM, amrex::The_Async_Arena());
      const Array4<Real> annular_prims =
          rz_axial_annular_prims_owner->array();
      ParallelFor(
          prim_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real cell_average[cls_t::NCONS];
            Real cell_primitives[cls_t::NPRIM] = {Real(0.0)};
            for (int n = 0; n < cls_t::NCONS; ++n) {
              cell_average[n] = rhs(i, j, k, n);
            }
            cls->cons2prims_point(cell_average, cell_primitives);
            for (int n = 0; n < cls_t::NPRIM; ++n) {
              annular_prims(i, j, k, n) = cell_primitives[n];
            }
          });
      rz_axial_annular_prims =
          rz_axial_annular_prims_owner->const_array();
    }

    // for each direction
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {

      const Box& flxbx = amrex::surroundingNodes(bx, dir);

      auto const& flx = flxt[dir]->array(); // snm
#if defined(CNS_IBM_VALIDATION) && (AMREX_USE_GPIBM || CNS_USE_EB)
      Array4<Real> split_positive;
      Array4<Real> split_negative;
      if (capture_split_flux) {
        AMREX_ALWAYS_ASSERT((*split_positive_audit)[dir] != nullptr);
        AMREX_ALWAYS_ASSERT((*split_negative_audit)[dir] != nullptr);
        split_positive = (*split_positive_audit)[dir]->array();
        split_negative = (*split_negative_audit)[dir]->array();
      }
#endif
      const bool rz_pressure_split = geom.IsRZ() && dir == 0 &&
                                     (s_rz_pressure_split != 0);

      ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept { // [=] only: do not capture 'this' in GPU lambda
        // NVCC must see these captures before the scheme-dependent
        // constexpr branches below.
        (void)llf_weno_epsilon_relative;
        (void)llf_teno_cutoff;

        IntVect iv(AMREX_D_DECL(i, j, k));
        IntVect ivd(IntVect::TheDimensionVector(dir));

        // comm/comp overlap shell pass: skip faces interior to skip_cells
        // (both adjacent cells inside) — they were computed in Pass 1.
        if (skip_ok && skipbox.contains(iv) && skipbox.contains(iv - ivd)) {
          return;
        }

        int gl = ng, gr = ng;  // ghost point position on left and right

        // modify stencils near IBM:  ibMarkers(iv,0) true means solid, false means fluid
#if AMREX_USE_GPIBM
        for (int mm = 0; mm < ng; ++mm) {
          if (ibMarkers(iv + mm * ivd, 0))       gr = amrex::min(gr, mm);
          if (ibMarkers(iv - (mm + 1) * ivd, 0)) gl = amrex::min(gl, mm);
        }
        if (gl == 0 && gr == 0) {
          return;  // skip solid cells
        }
        if ((gl == 0 && gr == 1) || (gl == 1 && gr == 0)) {
          AMREX_ASSERT_WITH_MESSAGE(false, "Cell is fluid but both neighbors are solid: no valid stencil");
        }
#endif
        // modify stencils near EBM  ibMarkers(iv,0) true means solid, false means fluid
#if CNS_USE_EB
        // TODO
        if (ibMarkers(iv +  2*ivd, 0 )) gr = 1;
        if (ibMarkers(iv +  ivd, 0)   ) gr = 0;
        if (ibMarkers(iv -  2*ivd, 0) ) gl = 0;

        if (gl == 0 && gr == 0) {
          return;  // skip solid cells
        }
#endif

#if AMREX_USE_GPIBM
        amrex::GpuArray<Real, 3> fp_point_weight3{
            Real(5.0), Real(10.0), Real(1.0)};
        amrex::GpuArray<Real, 3> fm_point_weight3{
            Real(5.0), Real(10.0), Real(1.0)};
        // Build the complete IBM usability/candidate mask before any stencil
        // state is read.  Fluid and reconstructed GP cells are usable;
        // interior solid cells are not.  point_needed is the union of the
        // surviving left/right candidate supports and therefore defines the
        // only cells that LLF splitting and its wave-speed estimate may read.
        bool usable[2 * ng];
        bool point_needed[2 * ng] = {};
        for (int m = 0; m < 2 * ng; ++m) {
          const IntVect c = iv + (m - ng) * ivd;
          usable[m] = ibm_flux::is_usable(c, ibMarkers);
        }

        amrex::GpuArray<Real, 3> fp_weight3{Real(3.0), Real(6.0), Real(1.0)};
        amrex::GpuArray<Real, 3> fm_weight3{Real(3.0), Real(6.0), Real(1.0)};
        amrex::GpuArray<Real, 4> fp_weight4{Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
        amrex::GpuArray<Real, 4> fm_weight4{Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
        bool has_fp = true;
        bool has_fm = true;

        if constexpr (Scheme::ncand == 3) {
          const bool fp_ok0 = usable[2] && usable[3] && usable[4];
          const bool fp_ok1 = usable[1] && usable[2] && usable[3];
          const bool fp_ok2 = usable[0] && usable[1] && usable[2];
          const bool fm_ok0 = usable[1] && usable[2] && usable[3];
          const bool fm_ok1 = usable[2] && usable[3] && usable[4];
          const bool fm_ok2 = usable[3] && usable[4] && usable[5];

          // Keep the nominal WENO weights and only remove invalid candidates.
          // At an immersed boundary the surviving polynomials can straddle a
          // reflected ghost extension; re-optimizing them as one smooth
          // polynomial family destroys the wall-normal momentum coupling.
          if (!fp_ok0) {
            fp_weight3[0] = Real(0.0);
            fp_point_weight3[0] = Real(0.0);
          }
          if (!fp_ok1) {
            fp_weight3[1] = Real(0.0);
            fp_point_weight3[1] = Real(0.0);
          }
          if (!fp_ok2) {
            fp_weight3[2] = Real(0.0);
            fp_point_weight3[2] = Real(0.0);
          }
          if (!fm_ok0) {
            fm_weight3[0] = Real(0.0);
            fm_point_weight3[0] = Real(0.0);
          }
          if (!fm_ok1) {
            fm_weight3[1] = Real(0.0);
            fm_point_weight3[1] = Real(0.0);
          }
          if (!fm_ok2) {
            fm_weight3[2] = Real(0.0);
            fm_point_weight3[2] = Real(0.0);
          }

          if (fp_ok0) for (int m = 2; m <= 4; ++m) point_needed[m] = true;
          if (fp_ok1) for (int m = 1; m <= 3; ++m) point_needed[m] = true;
          if (fp_ok2) for (int m = 0; m <= 2; ++m) point_needed[m] = true;
          if (fm_ok0) for (int m = 1; m <= 3; ++m) point_needed[m] = true;
          if (fm_ok1) for (int m = 2; m <= 4; ++m) point_needed[m] = true;
          if (fm_ok2) for (int m = 3; m <= 5; ++m) point_needed[m] = true;

          has_fp = fp_ok0 || fp_ok1 || fp_ok2;
          has_fm = fm_ok0 || fm_ok1 || fm_ok2;
        } else if constexpr (Scheme::ncand == 4) {
          const bool fp_ok0 = usable[2] && usable[3] && usable[4];
          const bool fp_ok1 = usable[1] && usable[2] && usable[3];
          const bool fp_ok2 = usable[0] && usable[1] && usable[2];
          const bool fp_ok3 = usable[2] && usable[3] && usable[4] && usable[5];
          const bool fm_ok0 = usable[1] && usable[2] && usable[3];
          const bool fm_ok1 = usable[2] && usable[3] && usable[4];
          const bool fm_ok2 = usable[3] && usable[4] && usable[5];
          const bool fm_ok3 = usable[0] && usable[1] && usable[2] && usable[3];

          if (!fp_ok0) fp_weight4[0] = Real(0.0);
          if (!fp_ok1) fp_weight4[1] = Real(0.0);
          if (!fp_ok2) fp_weight4[2] = Real(0.0);
          if (!fp_ok3) fp_weight4[3] = Real(0.0);
          if (!fm_ok0) fm_weight4[0] = Real(0.0);
          if (!fm_ok1) fm_weight4[1] = Real(0.0);
          if (!fm_ok2) fm_weight4[2] = Real(0.0);
          if (!fm_ok3) fm_weight4[3] = Real(0.0);

          if (fp_ok0) for (int m = 2; m <= 4; ++m) point_needed[m] = true;
          if (fp_ok1) for (int m = 1; m <= 3; ++m) point_needed[m] = true;
          if (fp_ok2) for (int m = 0; m <= 2; ++m) point_needed[m] = true;
          if (fp_ok3) for (int m = 2; m <= 5; ++m) point_needed[m] = true;
          if (fm_ok0) for (int m = 1; m <= 3; ++m) point_needed[m] = true;
          if (fm_ok1) for (int m = 2; m <= 4; ++m) point_needed[m] = true;
          if (fm_ok2) for (int m = 3; m <= 5; ++m) point_needed[m] = true;
          if (fm_ok3) for (int m = 0; m <= 3; ++m) point_needed[m] = true;

          has_fp = fp_ok0 || fp_ok1 || fp_ok2 || fp_ok3;
          has_fm = fm_ok0 || fm_ok1 || fm_ok2 || fm_ok3;
        }
        const bool force_first_order = !has_fp || !has_fm;

        // A face adjacent to fluid must have a valid reconstructed GP on its
        // solid side.  P0 initialisation fails fast if this invariant cannot
        // be met; retain a device assertion here as a last line of defence.
        if (force_first_order) {
          AMREX_ASSERT(usable[ng - 1] && usable[ng]);
        }
#else
        const bool force_first_order = false;
        const amrex::GpuArray<Real, 3> fp_point_weight3{
            Real(5.0), Real(10.0), Real(1.0)};
        const amrex::GpuArray<Real, 3> fm_point_weight3{
            Real(5.0), Real(10.0), Real(1.0)};
#endif

        bool use_rz_annular_prims =
            rz_euler_annular_average && dir == 0;
#if AMREX_USE_GPIBM
        // Reconstructed IBM ghost states are point-reflected wall data, not
        // annular cell averages.  Keep the established closure whenever the
        // active WENO stencil touches an immersed solid.
        use_rz_annular_prims =
            use_rz_annular_prims &&
            ibm_flux::stencil_all_fluid(
                iv, ivd, -ng, 2 * ng, ibMarkers);
#elif CNS_USE_EB
        use_rz_annular_prims = false;
#endif
        const bool use_rz_axial_annular_prims =
            rz_euler_axial_annular_flux && dir == 1;
        const Array4<const Real> flux_prims =
            use_rz_axial_annular_prims
                ? rz_axial_annular_prims
                : (use_rz_annular_prims ? rz_annular_prims : prims_in);
        // The annular BI-CWLS path already supplies centre-state primitives at
        // real shared GPs. Use point-flux candidates there, while retaining
        // prims_in: rebuilding a GP centre state from annular neighbours would
        // require deeper solid cells that a one-layer GP closure does not own.
        bool use_rz_point_flux =
            rz_euler_point_flux && dir == 0 &&
            (use_rz_annular_prims || rz_gp_annular_bic);
        const bool rz_point_flux_requested = use_rz_point_flux;
        bool rz_point_gate_would_block = false;
        Real rz_point_gate_sensor = Real(0.0);

#if (AMREX_SPACEDIM == 2)
        if (use_rz_point_flux &&
            (rz_euler_point_flux_shock_gate ||
             capture_rz_shock_face_audit)) {
          bool gate_usable = true;
#if AMREX_USE_GPIBM
          // The radial reconstruction stencil was already proven all-fluid.
          // The extra axial probes close the two-dimensional shock sensor; if
          // one touches solid, conservatively use the established FV-WENO
          // reconstruction rather than reading a reflected point state.
          const IntVect gate_axial_vector =
              IntVect::TheDimensionVector(1);
          for (int radial_side = -1; radial_side <= 0; ++radial_side) {
            const IntVect radial_cell = iv + radial_side * ivd;
            for (int axial_offset = -1; axial_offset <= 1;
                 ++axial_offset) {
              gate_usable =
                  gate_usable &&
                  !ibm_flux::is_solid(
                      radial_cell + axial_offset * gate_axial_vector,
                      ibMarkers);
            }
          }
#endif
          if (!gate_usable) {
            rz_point_gate_would_block = true;
            if (rz_euler_point_flux_shock_gate) {
              use_rz_point_flux = false;
            }
          } else {
            Real gate_pressure_min = std::numeric_limits<Real>::max();
            Real gate_pressure_max = Real(0.0);
            // Scan the exact six-cell radial reconstruction support.  This
            // catches an oblique shock crossing the point-flux stencil.
            for (int m = 0; m < 2 * ng; ++m) {
              const IntVect c = iv + (m - ng) * ivd;
              const Real pressure = prims_in(c, cls_t::QPRES);
              gate_pressure_min = amrex::min(gate_pressure_min, pressure);
              gate_pressure_max = amrex::max(gate_pressure_max, pressure);
            }
            // Also sample the central radial pair in the axial direction.
            // This catches an axis-normal bow/Mach shock whose radial states
            // are locally identical but whose neighbouring axial cells are
            // discontinuous.
            const IntVect ivz = IntVect::TheDimensionVector(1);
            for (int radial_side = -1; radial_side <= 0; ++radial_side) {
              const IntVect radial_cell = iv + radial_side * ivd;
              for (int axial_offset = -1; axial_offset <= 1;
                   ++axial_offset) {
                const Real pressure = prims_in(
                    radial_cell + axial_offset * ivz, cls_t::QPRES);
                gate_pressure_min = amrex::min(gate_pressure_min, pressure);
                gate_pressure_max = amrex::max(gate_pressure_max, pressure);
              }
            }
            const Real gate_pressure_scale = amrex::max(
                gate_pressure_max + gate_pressure_min,
                std::numeric_limits<Real>::min());
            const Real gate_pressure_sensor =
                (gate_pressure_max - gate_pressure_min) /
                gate_pressure_scale;
            rz_point_gate_sensor = gate_pressure_sensor;
            rz_point_gate_would_block =
                gate_pressure_sensor >= rz_grid_aligned_shock_threshold;
            if (rz_euler_point_flux_shock_gate) {
              use_rz_point_flux = !rz_point_gate_would_block;
            }
          }
        }
#endif

        bool axis_shock_trigger = false;
        bool axis_shock_sensor_usable = false;
        bool axis_incident_flow_matches =
            rz_grid_aligned_shock_flow_sign == 0;
        Real axis_shock_pressure_sensor = Real(0.0);
#if (AMREX_SPACEDIM == 2)
        const bool evaluate_axis_shock_sensor =
            (rz_euler_h_correction && dir == 0) ||
            (rz_grid_aligned_shock_filter &&
             (!rz_grid_aligned_shock_transverse_only || dir == 0)) ||
            (capture_rz_shock_face_audit && dir == 0);
        if (evaluate_axis_shock_sensor) {
          bool sensor_usable = true;
#if AMREX_USE_GPIBM
          sensor_usable = ibm_flux::stencil_all_fluid(
              iv, ivd, -ng, 2 * ng, ibMarkers);
          if (sensor_usable && dir == 0) {
            const IntVect filter_axial_vector =
                IntVect::TheDimensionVector(1);
            for (int radial_side = -1; radial_side <= 0; ++radial_side) {
              const IntVect radial_cell = iv + radial_side * ivd;
              for (int axial_offset = -1; axial_offset <= 1;
                   ++axial_offset) {
                sensor_usable =
                    sensor_usable &&
                    !ibm_flux::is_solid(
                        radial_cell +
                            axial_offset * filter_axial_vector,
                        ibMarkers);
              }
            }
          }
#elif CNS_USE_EB
          sensor_usable = false;
#endif
          axis_shock_sensor_usable = sensor_usable;
          if (sensor_usable) {
            Real pressure_min = std::numeric_limits<Real>::max();
            Real pressure_max = Real(0.0);
            bool incident_flow_matches =
                rz_grid_aligned_shock_flow_sign == 0;
            if (dir == 0) {
              const IntVect ivz = IntVect::TheDimensionVector(1);
              for (int radial_side = -1; radial_side <= 0; ++radial_side) {
                const IntVect radial_cell = iv + radial_side * ivd;
                for (int axial_offset = -1; axial_offset <= 1;
                     ++axial_offset) {
                  const Real pressure = flux_prims(
                      radial_cell + axial_offset * ivz, cls_t::QPRES);
                  pressure_min = amrex::min(pressure_min, pressure);
                  pressure_max = amrex::max(pressure_max, pressure);
                }
                if (rz_grid_aligned_shock_flow_sign != 0) {
                  for (int axial_offset = -1; axial_offset <= 0;
                       ++axial_offset) {
                    const IntVect left_cell =
                        radial_cell + axial_offset * ivz;
                    const IntVect right_cell = left_cell + ivz;
                    const Real pressure_left =
                        flux_prims(left_cell, cls_t::QPRES);
                    const Real pressure_right =
                        flux_prims(right_cell, cls_t::QPRES);
                    const Real pair_scale = amrex::max(
                        pressure_left + pressure_right,
                        std::numeric_limits<Real>::min());
                    const Real pair_sensor =
                        std::abs(pressure_right - pressure_left) / pair_scale;
                    const IntVect upstream_cell =
                        pressure_left <= pressure_right ? left_cell
                                                        : right_cell;
                    const Real upstream_axial_velocity =
                        flux_prims(upstream_cell, cls_t::QU + 1);
                    incident_flow_matches =
                        incident_flow_matches ||
                        (pair_sensor >= rz_grid_aligned_shock_threshold &&
                         Real(rz_grid_aligned_shock_flow_sign) *
                                 upstream_axial_velocity >
                             Real(0.0));
                  }
                }
              }
            } else {
              const Real pressure_left =
                  flux_prims(iv - ivd, cls_t::QPRES);
              const Real pressure_right = flux_prims(iv, cls_t::QPRES);
              pressure_min = amrex::min(pressure_left, pressure_right);
              pressure_max = amrex::max(pressure_left, pressure_right);
              if (rz_grid_aligned_shock_flow_sign != 0) {
                const IntVect upstream_cell =
                    pressure_left <= pressure_right ? iv - ivd : iv;
                const Real upstream_axial_velocity =
                    flux_prims(upstream_cell, cls_t::QU + 1);
                incident_flow_matches =
                    Real(rz_grid_aligned_shock_flow_sign) *
                        upstream_axial_velocity >
                    Real(0.0);
              }
            }
            const Real pressure_scale = amrex::max(
                pressure_max + pressure_min,
                std::numeric_limits<Real>::min());
            const Real pressure_sensor =
                (pressure_max - pressure_min) / pressure_scale;
            axis_shock_pressure_sensor = pressure_sensor;
            axis_incident_flow_matches = incident_flow_matches;
            if (dir == 0) {
              // A radial face supplies the conservative transverse coupling.
              // It is active wherever either radial line contains a strong
              // axial shock, including an already displaced shock front.
              axis_shock_trigger =
                  pressure_sensor >= rz_grid_aligned_shock_threshold &&
                  incident_flow_matches;
            } else {
              // Reduce the normal reconstruction only when the pressure jump
              // is predominantly axial. This avoids imposing a radial band
              // boundary and naturally turns the filter off on oblique shocks.
              const IntVect ivr = IntVect::TheDimensionVector(0);
              const IntVect left_cell = iv - ivd;
              const IntVect right_cell = iv;
              const Real transverse_left = Real(0.5) * std::abs(
                  flux_prims(left_cell + ivr, cls_t::QPRES) -
                  flux_prims(left_cell - ivr, cls_t::QPRES));
              const Real transverse_right = Real(0.5) * std::abs(
                  flux_prims(right_cell + ivr, cls_t::QPRES) -
                  flux_prims(right_cell - ivr, cls_t::QPRES));
              const Real transverse_jump =
                  amrex::max(transverse_left, transverse_right);
              const Real normal_jump = pressure_max - pressure_min;
              axis_shock_trigger =
                  pressure_sensor >= rz_grid_aligned_shock_threshold &&
                  incident_flow_matches &&
                  normal_jump >=
                      rz_grid_aligned_shock_ratio * transverse_jump;
            }
          }
        }
#endif

        Real alpha_raw = Real(0.0);
#if AMREX_USE_GPIBM
        for (int m = 0; m < 2 * ng; ++m) {
          const bool read_point = force_first_order
                                      ? (m == ng - 1 || m == ng)
                                      : point_needed[m];
          if (!read_point) continue;
          const IntVect c = iv + (m - ng) * ivd;
#ifdef CNS_IBM_VALIDATION
          // Idempotent validation flag. Marker component 2 exists only in a
          // CNS_IBM_VALIDATION GPIBM build and records the exact GP copies
          // entering the active LLF split/reconstruction stencil.
          if (ibm_flux::is_solid(c, ibMarkers) &&
              ibm_flux::is_usable(c, ibMarkers)) {
            ibMarkers(c, 2) = uint8_t(1);
          }
#endif
          Real wavespeed = std::abs(flux_prims(c, cls_t::QU + dir)) +
                           flux_prims(c, cls_t::QC);
#ifdef CNS_IBM_VALIDATION
          if (freeze_coefficients) {
            wavespeed =
                std::abs(coefficient_prims(c, cls_t::QU + dir)) +
                coefficient_prims(c, cls_t::QC);
          }
#endif
          if (wavespeed > alpha_raw) alpha_raw = wavespeed;
        }
#else
        alpha_raw = cls->max_char_speed(iv, dir, ng, flux_prims);
#endif
        AMREX_ASSERT_WITH_MESSAGE(alpha_raw > Real(0.0), "Non-positive LLF alpha in WENO/TENO flux split");
        const Real alpha = (alpha_raw > Real(0.0))
                               ? alpha_raw
                               : std::numeric_limits<Real>::epsilon();

        // --- Fix C: positivity-preserving first-order fallback --------------
        // The high-order WENO/TENO reconstruction of the LLF-split conservative
        // fluxes is not positivity-preserving: in the under-expanded jet near-
        // field it overshoots the energy component and the update yields
        // rho*e < 0 (negative internal energy) even though LLF keeps rho > 0.
        // Detect a rarefaction / near-vacuum pocket (deep local minimum of
        // density OR pressure across this face, along this direction) and use
        // the first-order LLF (Rusanov) flux there, which is positivity-robust
        // under the usual CFL. A shock is a monotone jump, not a local minimum,
        // so the bow shock is not smeared. Mirrors the HLLC rarefaction_pocket
        // sensor in Riemann.h. No-op in smooth flow (a 10x local drop is needed).
        bool first_order_flux =
            force_first_order ||
            (rz_grid_aligned_shock_filter && axis_shock_trigger);
        if (!first_order_flux) {
          constexpr Real RAREFY_RATIO = Real(0.1);
          bool pocket = false;
#if AMREX_USE_GPIBM
          // The sensor uses cells m=1..4.  Disable it near an IBM surface
          // unless every one is usable; it is a robustness heuristic and
          // must never reintroduce reads from masked interior-solid cells.
          const bool sensor_usable = usable[1] && usable[2] &&
                                     usable[3] && usable[4];
#else
          const bool sensor_usable = true;
#endif
          for (int side = -1; side <= 0 && sensor_usable && !pocket; ++side) {
            const IntVect c = iv + side * ivd;
            Real rc = flux_prims(c, cls_t::QRHO);
            Real pc = flux_prims(c, cls_t::QPRES);
            Real rminus = flux_prims(c - ivd, cls_t::QRHO);
            Real rplus = flux_prims(c + ivd, cls_t::QRHO);
            Real pminus = flux_prims(c - ivd, cls_t::QPRES);
            Real pplus = flux_prims(c + ivd, cls_t::QPRES);
#ifdef CNS_IBM_VALIDATION
            if (freeze_coefficients) {
              rc = coefficient_prims(c, cls_t::QRHO);
              pc = coefficient_prims(c, cls_t::QPRES);
              rminus = coefficient_prims(c - ivd, cls_t::QRHO);
              rplus = coefficient_prims(c + ivd, cls_t::QRHO);
              pminus = coefficient_prims(c - ivd, cls_t::QPRES);
              pplus = coefficient_prims(c + ivd, cls_t::QPRES);
            }
#endif
            pocket =
                (rc < RAREFY_RATIO * amrex::max(rminus, rplus)) ||
                (pc < RAREFY_RATIO * amrex::max(pminus, pplus));
          }
          first_order_flux = pocket;
        }
        if (first_order_flux) {
          Real fL[cls_t::NCONS], fR[cls_t::NCONS];
          Real uL[cls_t::NCONS], uR[cls_t::NCONS];
          cls->prims2flux(iv - ivd, dir, flux_prims, fL);
          cls->prims2cons(iv - ivd, flux_prims, uL);
          cls->prims2flux(iv, dir, flux_prims, fR);
          cls->prims2cons(iv, flux_prims, uR);
          if (rz_pressure_split) {
            fL[cls_t::UMX] -= flux_prims(iv - ivd, cls_t::QPRES);
            fR[cls_t::UMX] -= flux_prims(iv, cls_t::QPRES);
          }
          for (int n = 0; n < cls_t::NCONS; ++n) {
            const Real positive = Real(0.5) * (fL[n] + alpha * uL[n]);
            const Real negative = Real(0.5) * (fR[n] - alpha * uR[n]);
            // Preserve the production arithmetic exactly; the split values
            // are diagnostic and may differ from this sum by roundoff.
            flx(iv, n) = Real(0.5) * (fL[n] + fR[n]) -
                         Real(0.5) * alpha * (uR[n] - uL[n]);
#if defined(CNS_IBM_VALIDATION) && (AMREX_USE_GPIBM || CNS_USE_EB)
            if (capture_split_flux) {
              split_positive(iv, n) = positive;
              split_negative(iv, n) = negative;
            }
#endif
          }
          return;
        }
        // -------------------------------------------------------------------

#ifdef CNS_IBM_VALIDATION
        const auto roe_avg = freeze_coefficients
                                 ? cls->roe_avg_state(
                                       iv, dir, coefficient_prims)
                                 : cls->roe_avg_state(iv, dir, flux_prims);
#else
        const auto roe_avg = cls->roe_avg_state(iv, dir, flux_prims);
#endif

        Real cons[cls_t::NCONS], f[cls_t::NCONS], fp[2 * ng][cls_t::NCONS],
            fm[2 * ng][cls_t::NCONS];
#ifdef CNS_IBM_VALIDATION
        Real coefficient_cons[cls_t::NCONS];
        Real coefficient_flux[cls_t::NCONS];
        Real coefficient_fp[2 * ng][cls_t::NCONS];
        Real coefficient_fm[2 * ng][cls_t::NCONS];
#endif
        for (int m = 0; m < 2 * ng; ++m) {
#if AMREX_USE_GPIBM
          if (!point_needed[m]) {
            for (int n = 0; n < cls_t::NCONS; ++n) {
              fp[m][n] = Real(0.0);
              fm[m][n] = Real(0.0);
#ifdef CNS_IBM_VALIDATION
              coefficient_fp[m][n] = Real(0.0);
              coefficient_fm[m][n] = Real(0.0);
#endif
            }
            continue;
          }
#endif
          // LLF splitting into left- and right-running fluxes
          const IntVect c = iv + (m - ng) * ivd;
          cls->prims2flux(c, dir, flux_prims, f);
          cls->prims2cons(c, flux_prims, cons);
          if (rz_pressure_split) {
            f[cls_t::UMX] -= flux_prims(c, cls_t::QPRES);
          }

          for (int n = 0; n < cls_t::NCONS; ++n) {
            fp[m][n] = 0.5 * (cons[n] + f[n] / alpha);
            fm[m][n] = 0.5 * (cons[n] - f[n] / alpha);
          }

          // Convert into characteristic variables
          cls->cons2char(roe_avg, fp[m]);
          cls->cons2char(roe_avg, fm[m]);
#ifdef CNS_IBM_VALIDATION
          if (freeze_coefficients) {
            cls->prims2flux(c, dir, coefficient_prims, coefficient_flux);
            cls->prims2cons(c, coefficient_prims, coefficient_cons);
            if (rz_pressure_split) {
              coefficient_flux[cls_t::UMX] -=
                  coefficient_prims(c, cls_t::QPRES);
            }
            for (int n = 0; n < cls_t::NCONS; ++n) {
              coefficient_fp[m][n] =
                  Real(0.5) * (coefficient_cons[n] +
                               coefficient_flux[n] / alpha);
              coefficient_fm[m][n] =
                  Real(0.5) * (coefficient_cons[n] -
                               coefficient_flux[n] / alpha);
            }
            cls->cons2char(roe_avg, coefficient_fp[m]);
            cls->cons2char(roe_avg, coefficient_fm[m]);
          }
#endif
        }

        Real face_weno_epsilon = OldReconScheme::WenoZ5::legacy_epsilon();
        if (llf_weno_epsilon_relative > Real(0.0)) {
          Real face_characteristic_scale = Real(0.0);
          for (int m = 0; m < 2 * ng; ++m) {
            for (int n = 0; n < cls_t::NCONS; ++n) {
              face_characteristic_scale =
                  amrex::max(face_characteristic_scale, std::abs(fp[m][n]));
              face_characteristic_scale =
                  amrex::max(face_characteristic_scale, std::abs(fm[m][n]));
            }
          }
          face_weno_epsilon = amrex::max(
              std::numeric_limits<Real>::min(),
              llf_weno_epsilon_relative * face_characteristic_scale *
                  face_characteristic_scale);
        }

        // Reconstruct with upwind stencils
        Real fpL[cls_t::NCONS], fmR[cls_t::NCONS], s[2 * ng];
#ifdef CNS_IBM_VALIDATION
        Real coefficient_stencil[2 * ng];
#endif
        for (int n = 0; n < cls_t::NCONS; ++n) {
          Scheme::left_stencil(n, fp, s);
#if AMREX_USE_GPIBM
          if constexpr (Scheme::ncand == 3) {
#ifdef CNS_IBM_VALIDATION
            if (freeze_coefficients) {
              Scheme::left_stencil(n, coefficient_fp, coefficient_stencil);
              fpL[n] = Scheme::recon_masked_frozen(
                  s, coefficient_stencil, fp_weight3);
            } else
#endif
            {
              fpL[n] = use_rz_point_flux
                           ? Scheme::recon_point_masked(s, fp_point_weight3)
                           : Scheme::recon_masked(s, fp_weight3);
            }
          } else if constexpr (Scheme::ncand == 4) {
            fpL[n] = Scheme::recon_masked(s, fp_weight4);
          } else {
            fpL[n] = Scheme::recon(s, gr, gl);
          }
#else
          if constexpr (Scheme::ncand == 3) {
            if (use_rz_point_flux) {
              fpL[n] =
                  Scheme::recon_point_masked(s, fp_point_weight3);
            } else if constexpr (
                std::is_same_v<Scheme, OldReconScheme::WenoZ5>) {
              fpL[n] =
                  llf_weno_epsilon_relative > Real(0.0)
                      ? Scheme::recon_scaled(
                            s, gr, gl, face_weno_epsilon)
                      : Scheme::recon(s, gr, gl);
            } else if constexpr (
                std::is_same_v<Scheme, OldReconScheme::Teno5>) {
              fpL[n] = Scheme::recon_with_cutoff(
                  s, gr, gl, llf_teno_cutoff);
            } else {
              fpL[n] = Scheme::recon(s, gr, gl);
            }
          } else {
            fpL[n] = Scheme::recon(s, gr, gl);
          }
#endif

          Scheme::right_stencil(n, fm, s);
#if AMREX_USE_GPIBM
          if constexpr (Scheme::ncand == 3) {
#ifdef CNS_IBM_VALIDATION
            if (freeze_coefficients) {
              Scheme::right_stencil(n, coefficient_fm, coefficient_stencil);
              fmR[n] = Scheme::recon_masked_frozen(
                  s, coefficient_stencil, fm_weight3);
            } else
#endif
            {
              fmR[n] = use_rz_point_flux
                           ? Scheme::recon_point_masked(s, fm_point_weight3)
                           : Scheme::recon_masked(s, fm_weight3);
            }
          } else if constexpr (Scheme::ncand == 4) {
            fmR[n] = Scheme::recon_masked(s, fm_weight4);
          } else {
            fmR[n] = Scheme::recon(s, gl, gr);
          }
#else
          if constexpr (Scheme::ncand == 3) {
            if (use_rz_point_flux) {
              fmR[n] =
                  Scheme::recon_point_masked(s, fm_point_weight3);
            } else if constexpr (
                std::is_same_v<Scheme, OldReconScheme::WenoZ5>) {
              fmR[n] =
                  llf_weno_epsilon_relative > Real(0.0)
                      ? Scheme::recon_scaled(
                            s, gl, gr, face_weno_epsilon)
                      : Scheme::recon(s, gl, gr);
            } else if constexpr (
                std::is_same_v<Scheme, OldReconScheme::Teno5>) {
              fmR[n] = Scheme::recon_with_cutoff(
                  s, gl, gr, llf_teno_cutoff);
            } else {
              fmR[n] = Scheme::recon(s, gl, gr);
            }
          } else {
            fmR[n] = Scheme::recon(s, gl, gr);
          }
#endif
        }

        // Convert back to conservative variables
        cls->char2cons(roe_avg, fpL);
        cls->char2cons(roe_avg, fmR);

        Real base_weno_flux[cls_t::NCONS];
        Real rz_h_correction_flux[cls_t::NCONS];
        for (int n = 0; n < cls_t::NCONS; ++n) {
          const Real positive = alpha * fpL[n];
          const Real negative = -alpha * fmR[n];
          // Preserve the production arithmetic exactly; the split values
          // are diagnostic and may differ from this sum by roundoff.
          flx(iv, n) = alpha * (fpL[n] - fmR[n]);
          base_weno_flux[n] = flx(iv, n);
          rz_h_correction_flux[n] = Real(0.0);
#if defined(CNS_IBM_VALIDATION) && (AMREX_USE_GPIBM || CNS_USE_EB)
          if (capture_split_flux) {
            split_positive(iv, n) = positive;
            split_negative(iv, n) = negative;
          }
#endif

        }
        Real rz_h_correction_weight = Real(0.0);
#if (AMREX_SPACEDIM == 2)
        if (rz_euler_h_correction && dir == 0 && axis_shock_trigger) {
          const Real sensor_denominator = amrex::max(
              Real(1.0) - rz_grid_aligned_shock_threshold,
              std::numeric_limits<Real>::epsilon());
          const Real sensor_coordinate = amrex::min(
              Real(1.0),
              amrex::max(
                  Real(0.0),
                  (axis_shock_pressure_sensor -
                   rz_grid_aligned_shock_threshold) /
                      sensor_denominator));
          const Real smooth_sensor =
              sensor_coordinate * sensor_coordinate *
              (Real(3.0) - Real(2.0) * sensor_coordinate);
          rz_h_correction_weight =
              rz_euler_h_correction_strength * smooth_sensor;

          Real left_state[cls_t::NCONS];
          Real right_state[cls_t::NCONS];
          cls->prims2cons(iv - ivd, flux_prims, left_state);
          cls->prims2cons(iv, flux_prims, right_state);
          for (int n = 0; n < cls_t::NCONS; ++n) {
            rz_h_correction_flux[n] =
                -Real(0.5) * rz_h_correction_weight * alpha *
                (right_state[n] - left_state[n]);
            flx(iv, n) += rz_h_correction_flux[n];
          }
        }
#endif
#if AMREX_USE_GPIBM
#ifdef CNS_IBM_VALIDATION
        if (capture_rz_shock_face_audit && dir == 0 &&
            bx.contains(iv)) {
          using namespace cerisse::rz_shock_face_diag;
          Real fL[cls_t::NCONS], fR[cls_t::NCONS];
          Real uL[cls_t::NCONS], uR[cls_t::NCONS];
          const IntVect left_cell = iv - ivd;
          cls->prims2flux(left_cell, dir, flux_prims, fL);
          cls->prims2cons(left_cell, flux_prims, uL);
          cls->prims2flux(iv, dir, flux_prims, fR);
          cls->prims2cons(iv, flux_prims, uR);
          if (rz_pressure_split) {
            fL[cls_t::UMX] -= flux_prims(left_cell, cls_t::QPRES);
            fR[cls_t::UMX] -= flux_prims(iv, cls_t::QPRES);
          }

          rz_shock_face_diagnostic(iv, valid) = Real(1.0);
          rz_shock_face_diagnostic(iv, point_flux_requested) =
              rz_point_flux_requested ? Real(1.0) : Real(0.0);
          rz_shock_face_diagnostic(iv, point_flux_used) =
              use_rz_point_flux ? Real(1.0) : Real(0.0);
          rz_shock_face_diagnostic(iv, point_gate_would_block) =
              rz_point_gate_would_block ? Real(1.0) : Real(0.0);
          rz_shock_face_diagnostic(iv, point_gate_sensor) =
              rz_point_gate_sensor;
          rz_shock_face_diagnostic(iv, shock_sensor_usable) =
              axis_shock_sensor_usable ? Real(1.0) : Real(0.0);
          rz_shock_face_diagnostic(iv, shock_pressure_sensor) =
              axis_shock_pressure_sensor;
          rz_shock_face_diagnostic(iv, incident_flow_matches) =
              axis_incident_flow_matches ? Real(1.0) : Real(0.0);
          rz_shock_face_diagnostic(iv, shock_trigger) =
              axis_shock_trigger ? Real(1.0) : Real(0.0);
          rz_shock_face_diagnostic(iv, h_correction_active) =
              rz_h_correction_weight > Real(0.0) ? Real(1.0)
                                                  : Real(0.0);
          rz_shock_face_diagnostic(iv, h_correction_weight) =
              rz_h_correction_weight;
          rz_shock_face_diagnostic(iv, llf_alpha) = alpha;
          rz_shock_face_diagnostic(iv, pair_density_min) =
              amrex::min(flux_prims(left_cell, cls_t::QRHO),
                         flux_prims(iv, cls_t::QRHO));
          rz_shock_face_diagnostic(iv, pair_pressure_min) =
              amrex::min(flux_prims(left_cell, cls_t::QPRES),
                         flux_prims(iv, cls_t::QPRES));
          for (int n = 0; n < cls_t::NCONS; ++n) {
            const Real llf =
                Real(0.5) * (fL[n] + fR[n]) -
                Real(0.5) * alpha * (uR[n] - uL[n]);
            rz_shock_face_diagnostic(iv, weno_flux(cls_t::NCONS) + n) =
                base_weno_flux[n];
            rz_shock_face_diagnostic(iv, llf_flux(cls_t::NCONS) + n) =
                llf;
            rz_shock_face_diagnostic(
                iv, llf_minus_weno(cls_t::NCONS) + n) =
                llf - base_weno_flux[n];
            rz_shock_face_diagnostic(iv, h_correction(cls_t::NCONS) + n) =
                rz_h_correction_flux[n];
            rz_shock_face_diagnostic(iv, selected_flux(cls_t::NCONS) + n) =
                flx(iv, n);
          }
        }
#endif
#endif
      });
    }  // end of for each direction
  }

#if defined(CNS_IBM_VALIDATION) && (AMREX_USE_GPIBM || CNS_USE_EB)
  void eflux_ibm_frozen_coefficients(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims_in,
      const amrex::Array4<const amrex::Real>& coefficient_prims,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
      const amrex::Array4<uint8_t>& ibMarkers)
  {
    eflux_ibm(geom, mfi, prims_in, flxt, rhs, cls, ibMarkers,
              nullptr, nullptr, &coefficient_prims);
  }

  void eflux_ibm_split_audit(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims_in,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
      const amrex::Array4<uint8_t>& ibMarkers,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& split_positive,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& split_negative)
  {
    eflux_ibm(geom, mfi, prims_in, flxt, rhs, cls, ibMarkers,
              &split_positive, &split_negative);
  }
#endif

#if defined(CNS_IBM_VALIDATION) && AMREX_USE_GPIBM
  void eflux_ibm_rz_shock_face_audit(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims_in,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
      const amrex::Array4<uint8_t>& ibMarkers,
      const amrex::Array4<amrex::Real>& diagnostic)
  {
    eflux_ibm(
        geom, mfi, prims_in, flxt, rhs, cls, ibMarkers,
        nullptr, nullptr, nullptr, &diagnostic);
  }
#endif

};

template <typename cls_t>
using llf_wenoz5_old_t =
    weno_old_t<OldReconScheme::WenoZ5, cls_t>;

template <typename cls_t>
using llf_teno5_old_t =
    weno_old_t<OldReconScheme::Teno5, cls_t>;

template <typename cls_t>
using llf_teno6_old_t =
    weno_old_t<OldReconScheme::Teno6, cls_t>;

#endif
