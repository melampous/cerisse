#ifndef Weno_H
#define Weno_H

#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <limits>

#include "Closures.h"

#define POWER2(x) ((x) * (x))
#define POWER6(x) ((x) * (x) * (x) * (x) * (x) * (x))

// Choices are WenoZ5, Teno5. (Teno6 is slightly unstable)
namespace ReconScheme {
struct WenoZ5 {
  static constexpr int ng = 3;
  static constexpr int ncand = 3;

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
   * \brief (FV)WENO-Z5. Ref https://doi.org/10.1016/j.jcp.2010.11.028.
   * \param[in] s Stencil values at cells [i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_masked(
    const amrex::Real s[5],
    const amrex::GpuArray<amrex::Real, 3>& optimal_weight) noexcept {
    using amrex::Real;
    constexpr Real eps = (std::numeric_limits<Real>::digits >= 53)
                             ? Real(1e-40)
                             : std::numeric_limits<Real>::epsilon();
    Real vr[3], beta[3], tmp;
    
    smoothness_indicator(s, beta);
    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);
    if (v0 && v1 && v2) {
      tmp = std::abs(beta[2] - beta[0]);
    } else if (v0 && v1) {
      tmp = std::abs(beta[1] - beta[0]);
    } else if (v1 && v2) {
      tmp = std::abs(beta[2] - beta[1]);
    } else {
      tmp = Real(0.0);
    }

    beta[2] = (1.0 + tmp / (eps + beta[2])) * optimal_weight[2];
    beta[1] = (1.0 + tmp / (eps + beta[1])) * optimal_weight[1];
    beta[0] = (1.0 + tmp / (eps + beta[0])) * optimal_weight[0];

    linear_polynomial_recon(s, vr);

    const Real denom = beta[2] + beta[1] + beta[0];
    if (!(denom > std::numeric_limits<Real>::min())) {
      if (v1) return vr[1] / Real(6.0);
      if (v0) return vr[0] / Real(6.0);
      if (v2) return vr[2] / Real(6.0);
      return vr[1] / Real(6.0);
    }
    tmp = 1.0 / denom;

    return tmp / Real(6.0) *
           (beta[2] * vr[2] + beta[1] * vr[1] + beta[0] * vr[0]);
  }

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
};  // struct WenoZ5

struct Teno5 : public WenoZ5 {
  /**
   * \brief (FV)TENO-5. Ref https://doi.org/10.1016/j.jcp.2015.10.037.
   * \param[in] s Stencil values at cells [i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
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

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon_masked(
      const amrex::Real s[5],
      const amrex::GpuArray<amrex::Real, 3>& optimal_weight) noexcept {
    using amrex::Real;

    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    // TENO cutoff CT: stencils with normalized smoothness below CT are
    // discarded (weight 0). Larger CT -> more stencils dropped near strong
    // gradients -> more upwind/dissipative -> more robust against the density
    // undershoot that crashed L3 NPR=3 (raised 1e-5 -> 1e-3 for robustness).
    constexpr Real cutoff = 1e-3;
    Real vr[3], beta[3], tmp;
    
    smoothness_indicator(s, beta);
    const bool v0 = optimal_weight[0] > Real(0.0);
    const bool v1 = optimal_weight[1] > Real(0.0);
    const bool v2 = optimal_weight[2] > Real(0.0);
    if (v0 && v1 && v2) {
      tmp = std::abs(std::abs(beta[2] - beta[0]) -
                     (beta[2] + 4.0 * beta[1] + beta[0]) / 6.0);
    } else if (v0 && v1) {
      tmp = std::abs(beta[1] - beta[0]);
    } else if (v1 && v2) {
      tmp = std::abs(beta[2] - beta[1]);
    } else {
      tmp = Real(0.0);
    }

    beta[2] = POWER6(1.0 + tmp / (eps + beta[2]));
    beta[1] = POWER6(1.0 + tmp / (eps + beta[1]));
    beta[0] = POWER6(1.0 + tmp / (eps + beta[0]));
    tmp = 1.0 / (beta[2] + beta[1] + beta[0]);
    beta[2] = beta[2] * tmp < cutoff ? 0.0 : optimal_weight[2];
    beta[1] = beta[1] * tmp < cutoff ? 0.0 : optimal_weight[1];
    beta[0] = beta[0] * tmp < cutoff ? 0.0 : optimal_weight[0];
    const Real denom = beta[2] + beta[1] + beta[0];
    if (!(denom > std::numeric_limits<Real>::min())) {
      return WenoZ5::recon_masked(s, optimal_weight);
    }
    tmp = 1.0 / denom;

    linear_polynomial_recon(s, vr);

    return tmp / Real(6.0) *
           (beta[2] * vr[2] + beta[1] * vr[1] + beta[0] * vr[0]);
  }
};  // struct Teno5

struct Teno6 {
  static constexpr int ng = 3;
  static constexpr int ncand = 4;

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
   * \brief (FV)TENO-6. Ref https://doi.org/10.1016/j.jcp.2015.10.037.
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
};  // namespace ReconScheme

// TODO: this class name is not very accurate. llf? high_order_fd?
template <typename Scheme, typename cls_t>
class weno_t {
  static constexpr int ng = Scheme::ng;  // number of ghost cells on each side
  static_assert(ng <= cls_t::NGHOST);    // ensure enough ghost cells

 public:
  // Only the WENO/TENO eflux implements the RZ pressure-split (it removes the
  // radial pressure flux to pair with compute_rhs's -dp/dr form). Other flux
  // schemes do NOT, so enabling cns.rz_pressure_split with them would
  // double-count pressure. compute_rhs guards on this trait.
  static constexpr bool rz_pressure_split_capable = true;

  AMREX_GPU_HOST_DEVICE
  weno_t() {}

  AMREX_GPU_HOST_DEVICE
  ~weno_t() {}

#if (AMREX_USE_GPIBM || CNS_USE_EB )
  void inline eflux_ibm(const amrex::Geometry& geom, const amrex::MFIter& mfi,
                        const amrex::Array4<const amrex::Real>& prims_in,
                        std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                        const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
                        //const amrex::Array4<const bool>& ibMarkers)
                        const amrex::Array4<uint8_t>& ibMarkers)
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

    const Box& bx = bx_region;
    const Box skipbox = skip_cells;
    const bool skip_ok = skipbox.ok();
    // Experimental RZ pressure split.  This must be paired with the matching
    // compute_rhs.cpp source form, and is only implemented for WENO/TENO.
    static const int s_rz_pressure_split = [] {
      int v = 0;
      amrex::ParmParse pp("cns");
      pp.query("rz_pressure_split", v);
      return v;
    }();

    // for each direction
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {

      const Box& flxbx = amrex::surroundingNodes(bx, dir);

      auto const& flx = flxt[dir]->array(); // snm
      const bool rz_pressure_split = geom.IsRZ() && dir == 0 &&
                                     (s_rz_pressure_split != 0);

      ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept { // [=] only: do not capture 'this' in GPU lambda
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
        

        const Real alpha_raw = cls->max_char_speed(iv, dir, ng, prims_in);
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
        {
          constexpr Real RAREFY_RATIO = Real(0.1);
          bool pocket = false;
          for (int side = -1; side <= 0 && !pocket; ++side) {  // cells iv-ivd, iv
            const IntVect c = iv + side * ivd;
            const Real rc = prims_in(c, cls_t::QRHO);
            const Real pc = prims_in(c, cls_t::QPRES);
            pocket = (rc < RAREFY_RATIO * amrex::max(prims_in(c - ivd, cls_t::QRHO),
                                                     prims_in(c + ivd, cls_t::QRHO))) ||
                     (pc < RAREFY_RATIO * amrex::max(prims_in(c - ivd, cls_t::QPRES),
                                                     prims_in(c + ivd, cls_t::QPRES)));
          }
          if (pocket) {
            Real fL[cls_t::NCONS], fR[cls_t::NCONS], uL[cls_t::NCONS], uR[cls_t::NCONS];
            cls->prims2flux(iv - ivd, dir, prims_in, fL);
            cls->prims2cons(iv - ivd, prims_in, uL);
            cls->prims2flux(iv, dir, prims_in, fR);
            cls->prims2cons(iv, prims_in, uR);
            if (rz_pressure_split) {
              fL[cls_t::UMX] -= prims_in(iv - ivd, cls_t::QPRES);
              fR[cls_t::UMX] -= prims_in(iv, cls_t::QPRES);
            }
            for (int n = 0; n < cls_t::NCONS; ++n)
              flx(iv, n) = Real(0.5) * (fL[n] + fR[n]) -
                           Real(0.5) * alpha * (uR[n] - uL[n]);
            return;
          }
        }
        // -------------------------------------------------------------------

        const auto roe_avg = cls->roe_avg_state(iv, dir, prims_in);

        Real cons[cls_t::NCONS], f[cls_t::NCONS], fp[2 * ng][cls_t::NCONS],
            fm[2 * ng][cls_t::NCONS];
        for (int m = 0; m < 2 * ng; ++m) {
          // LLF splitting into left- and right-running fluxes
          const IntVect c = iv + (m - ng) * ivd;
          cls->prims2flux(c, dir, prims_in, f);
          cls->prims2cons(c, prims_in, cons);
          if (rz_pressure_split) {
            f[cls_t::UMX] -= prims_in(c, cls_t::QPRES);
          }

          for (int n = 0; n < cls_t::NCONS; ++n) {
            fp[m][n] = 0.5 * (cons[n] + f[n] / alpha);
            fm[m][n] = 0.5 * (cons[n] - f[n] / alpha);
          }

          // Convert into characteristic variables
          cls->cons2char(roe_avg, fp[m]);
          cls->cons2char(roe_avg, fm[m]);
        }

#if AMREX_USE_GPIBM
        amrex::GpuArray<Real, 3> fp_weight3{Real(3.0), Real(6.0), Real(1.0)};
        amrex::GpuArray<Real, 3> fm_weight3{Real(3.0), Real(6.0), Real(1.0)};
        amrex::GpuArray<Real, 4> fp_weight4{Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
        amrex::GpuArray<Real, 4> fm_weight4{Real(6.0), Real(9.0), Real(1.0), Real(4.0)};
        bool has_fp = true;
        bool has_fm = true;
        if constexpr (Scheme::ncand == 3) {
          bool usable[2 * ng];
          for (int m = 0; m < 2 * ng; ++m) {
            const IntVect c = iv + (m - ng) * ivd;
            usable[m] = (!ibMarkers(c, 0)) || (ibMarkers(c, 1) != 0);
          }

          // Candidate maps for the WENO5/TENO5 stencils built above. A
          // reconstructed IBM ghost point is usable; an interior solid cell is
          // not. This masks only the candidate that actually reads bad data.
          const bool fp_ok0 = usable[2] && usable[3] && usable[4];
          const bool fp_ok1 = usable[1] && usable[2] && usable[3];
          const bool fp_ok2 = usable[0] && usable[1] && usable[2];
          const bool fm_ok0 = usable[1] && usable[2] && usable[3];
          const bool fm_ok1 = usable[2] && usable[3] && usable[4];
          const bool fm_ok2 = usable[3] && usable[4] && usable[5];

          if (!fp_ok0) fp_weight3[0] = Real(0.0);
          if (!fp_ok1) fp_weight3[1] = Real(0.0);
          if (!fp_ok2) fp_weight3[2] = Real(0.0);
          if (!fm_ok0) fm_weight3[0] = Real(0.0);
          if (!fm_ok1) fm_weight3[1] = Real(0.0);
          if (!fm_ok2) fm_weight3[2] = Real(0.0);

          has_fp = fp_ok0 || fp_ok1 || fp_ok2;
          has_fm = fm_ok0 || fm_ok1 || fm_ok2;
        } else if constexpr (Scheme::ncand == 4) {
          bool usable[2 * ng];
          for (int m = 0; m < 2 * ng; ++m) {
            const IntVect c = iv + (m - ng) * ivd;
            usable[m] = (!ibMarkers(c, 0)) || (ibMarkers(c, 1) != 0);
          }

          // Candidate maps for TENO6.  The first three candidates use three
          // points, and the extra candidate uses four points.  A reconstructed
          // IBM ghost is usable; an interior solid cell is not.
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

          has_fp = fp_ok0 || fp_ok1 || fp_ok2 || fp_ok3;
          has_fm = fm_ok0 || fm_ok1 || fm_ok2 || fm_ok3;
        }
        if (!has_fp || !has_fm) {
          Real fL[cls_t::NCONS], fR[cls_t::NCONS], uL[cls_t::NCONS], uR[cls_t::NCONS];
          cls->prims2flux(iv - ivd, dir, prims_in, fL);
          cls->prims2cons(iv - ivd, prims_in, uL);
          cls->prims2flux(iv, dir, prims_in, fR);
          cls->prims2cons(iv, prims_in, uR);
          if (rz_pressure_split) {
            fL[cls_t::UMX] -= prims_in(iv - ivd, cls_t::QPRES);
            fR[cls_t::UMX] -= prims_in(iv, cls_t::QPRES);
          }
          for (int n = 0; n < cls_t::NCONS; ++n)
            flx(iv, n) = Real(0.5) * (fL[n] + fR[n]) -
                         Real(0.5) * alpha * (uR[n] - uL[n]);
          return;
        }
#endif

        // Reconstruct with upwind stencils
        Real fpL[cls_t::NCONS], fmR[cls_t::NCONS], s[2 * ng];
        for (int n = 0; n < cls_t::NCONS; ++n) {
          Scheme::left_stencil(n, fp, s);
#if AMREX_USE_GPIBM
          if constexpr (Scheme::ncand == 3) {
            fpL[n] = Scheme::recon_masked(s, fp_weight3);
          } else if constexpr (Scheme::ncand == 4) {
            fpL[n] = Scheme::recon_masked(s, fp_weight4);
          } else {
            fpL[n] = Scheme::recon(s, gr, gl);
          }
#else
          fpL[n] = Scheme::recon(s, gr, gl);
#endif

          Scheme::right_stencil(n, fm, s);
#if AMREX_USE_GPIBM
          if constexpr (Scheme::ncand == 3) {
            fmR[n] = Scheme::recon_masked(s, fm_weight3);
          } else if constexpr (Scheme::ncand == 4) {
            fmR[n] = Scheme::recon_masked(s, fm_weight4);
          } else {
            fmR[n] = Scheme::recon(s, gl, gr);
          }
#else
          fmR[n] = Scheme::recon(s, gl, gr);
#endif
        }

        // Convert back to conservative variables
        cls->char2cons(roe_avg, fpL);
        cls->char2cons(roe_avg, fmR);

        for (int n = 0; n < cls_t::NCONS; ++n) {
                
          flx(iv, n) = alpha * (fpL[n] - fmR[n]);

        }
      });
    }  // end of for each direction
  }

///// OBSOLETE //////// // ........
#if (AMREX_USE_GPIBM || CNS_USE_EB )
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool fill_solid_prims(
      amrex::IntVect iv, amrex::IntVect ivd, int cdir,
      amrex::Array4<const amrex::Real> const &prims_in,
      amrex::Array4<amrex::Real> const &prims,
      amrex::Array4<const bool> const &ib_mask) {
    // Find the first solid point, which will be the ghost point
    int gl = ng, gr = ng;  // ghost point position on left and right
    for (int m = 0; m < ng; ++m) {
      if (ib_mask(iv + m * ivd,0)) {  // snm
        gr = std::min(gr, m);
      }
      if (ib_mask(iv - (m + 1) * ivd,0)) { //snm
        gl = std::min(gl, m);
      }
    }
    if (gl == 0 && gr == 0) return false;  // skip solid cells


  
    // 3. fill using 1st order BC (adiabatic no-slip)
    for (int m = 0; m < 2 * ng; ++m) {
      for (int n = 0; n < cls_t::NPRIM; ++n) {        
        prims(iv + (m - ng) * ivd, n) = prims_in(iv + (m - ng) * ivd, n);
      }
    }


    for (int m = 0; m < ng; ++m) {
   
      if (m >= gl) {
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(iv - (m + 1) * ivd, n) = prims(iv - (2 * gl - m) * ivd, n);
        }
        // prims(iv - (m + 1) * ivd, cls_t::QU + cdir) *= -1; // this is non-permeable
        for (int c = 0; c < amrex::SpaceDim; ++c) 
	    prims(iv - (m + 1) * ivd, cls_t::QU + c) *= -1; // this is no-slip non-permeable
      }
      if (m >= gr) {
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(iv + m * ivd, n) = prims(iv + (2 * gr - m - 1) * ivd, n);
        }
        // prims(iv + m * ivd, cls_t::QU + cdir) *= -1; // this is non-permeable
	for (int c = 0; c < amrex::SpaceDim; ++c)
            prims(iv + m * ivd, cls_t::QU + c) *= -1; // this is no-slip non-permeable
      }
    }

    return true;
  } // ........
#endif
};

#endif
