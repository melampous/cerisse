#ifndef AFD_H_
#define AFD_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <Riemann.h>
#include <IbmFluxUtils.h>

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
// corrections.  This production path is restricted to the Cartesian,
// five-equation ideal-gas system.  In GPIBM builds the full six-point AFD
// formula is retained on all-fluid stencils; wall-intersecting stencils use
// marker-masked point reconstruction with HLLC and a two-point LLF fallback.
// It needs three primitive ghost cells and supplies the Box overload required
// by communication overlap in no-IBM builds.
template <typename PointScheme, typename cls_t>
class afd_hllc_t : public riemann_t<false, cls_t> {
 public:
  static constexpr bool rz_pressure_split_capable = false;
  static constexpr int ng = 3;

  static_assert(NUM_SPECIES == 1,
                "AFD-HLLC currently supports one species only");
  static_assert(cls_t::NCONS == 5,
                "AFD-HLLC expects the five-equation ideal-gas system");
  static_assert(cls_t::NGHOST >= ng,
                "AFD-HLLC requires at least three primitive ghost cells");

  void eflux(const amrex::Geometry& geom, const amrex::MFIter& mfi,
             const amrex::Array4<const amrex::Real>& prims,
             std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
             const amrex::Array4<amrex::Real>& rhs, const cls_t* cls)
  {
    eflux(geom, mfi.tilebox(), prims, flxt, rhs, cls);
  }

  void eflux(const amrex::Geometry& geom, const amrex::Box& bx,
             const amrex::Array4<const amrex::Real>& prims,
             std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
             const amrex::Array4<amrex::Real>& /*rhs*/, const cls_t* cls,
             const amrex::Box& skip_cells = amrex::Box())
  {
    if (geom.IsRZ()) {
      amrex::Abort("AFD-HLLC is Cartesian-only; use LLF-WENO/TENO for RZ");
    }

    const RuntimeOptions opts = runtime_options();
    const int correction = opts.correction;
    const amrex::Real smoothness_threshold = opts.smoothness_threshold;
    const amrex::Real teno_cutoff = opts.teno_cutoff;
    const amrex::Box skipbox = skip_cells;
    const bool skip_ok = skipbox.ok();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const amrex::Box flxbx = amrex::surroundingNodes(bx, dir);
      const auto flx = flxt[dir]->array();
      const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
      amrex::ParallelFor<128>(
          flxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (skip_ok && skipbox.contains(iv) &&
                skipbox.contains(iv - ivd)) {
              return;
            }
            this->face_flux(iv, dir, prims, flx, *cls, correction,
                            smoothness_threshold, teno_cutoff);
          });
    }
  }

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void eflux_ibm(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<amrex::Real>& prims,
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& /*rhs*/, const cls_t* cls,
      const amrex::Array4<std::uint8_t>& ibMarkers)
  {
    if (geom.IsRZ()) {
      amrex::Abort("AFD-HLLC is Cartesian-only; use LLF-WENO/TENO for RZ");
    }

    const RuntimeOptions opts = runtime_options();
    const int correction = opts.correction;
    const amrex::Real smoothness_threshold = opts.smoothness_threshold;
    const amrex::Real teno_cutoff = opts.teno_cutoff;
    const amrex::Box bx = mfi.tilebox();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const amrex::Box flxbx = amrex::surroundingNodes(bx, dir);
      const auto flx = flxt[dir]->array();
      const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
      amrex::ParallelFor<128>(
          flxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (ibm_flux::is_solid_solid_face(iv, ivd, ibMarkers)) {
              ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
              return;
            }
            if (!ibm_flux::stencil_all_fluid(iv, ivd, -ng, 2 * ng,
                                             ibMarkers)) {
              this->face_flux_ibm_boundary(iv, dir, prims, flx, *cls,
                                           ibMarkers, teno_cutoff);
              return;
            }
            this->face_flux(iv, dir, prims, flx, *cls, correction,
                            smoothness_threshold, teno_cutoff);
          });
    }
  }
#endif

 private:
  struct RuntimeOptions {
    int correction;
    amrex::Real smoothness_threshold;
    amrex::Real teno_cutoff;
  };

  static RuntimeOptions runtime_options()
  {
    static const RuntimeOptions opts = [] {
      RuntimeOptions value{1, amrex::Real(0.08), amrex::Real(1.0e-4)};
      amrex::ParmParse pp("cns");
      pp.query("afd_correction", value.correction);
      pp.query("afd_smoothness_threshold", value.smoothness_threshold);
      pp.query("afd_teno_cutoff", value.teno_cutoff);
      if (value.correction != 0 && value.correction != 1) {
        amrex::Abort("cns.afd_correction must be 0 or 1");
      }
      if (!std::isfinite(value.smoothness_threshold) ||
          value.smoothness_threshold < amrex::Real(0.0) ||
          value.smoothness_threshold > amrex::Real(1.0)) {
        amrex::Abort("cns.afd_smoothness_threshold must be in [0,1]");
      }
      if (!std::isfinite(value.teno_cutoff) ||
          value.teno_cutoff < amrex::Real(0.0) ||
          value.teno_cutoff >= amrex::Real(1.0 / 3.0)) {
        amrex::Abort("cns.afd_teno_cutoff must be in [0,1/3)");
      }
      return value;
    }();
    return opts;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  valid_stencil(const amrex::IntVect& iv, const amrex::IntVect& ivd,
                const amrex::Array4<const amrex::Real>& q) noexcept
  {
    using amrex::Real;
    for (int m = 0; m < 6; ++m) {
      const amrex::IntVect p = iv + (m - 3) * ivd;
      const Real rho = q(p, cls_t::QRHO);
      const Real pressure = q(p, cls_t::QPRES);
      const Real gamma = q(p, cls_t::QG);
      const Real sound = q(p, cls_t::QC);
      const Real finite_sum = rho + pressure + gamma + sound +
                              q(p, cls_t::QU) + q(p, cls_t::QU + 1) +
                              q(p, cls_t::QU + 2);
      if (!(rho > Real(0.0)) || !(pressure > Real(0.0)) ||
          !(gamma > Real(1.0)) || !(sound > Real(0.0)) ||
          !std::isfinite(finite_sum)) {
        return false;
      }
    }
    return true;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  smooth_stencil(const amrex::IntVect& iv, const amrex::IntVect& ivd,
                 const amrex::Array4<const amrex::Real>& q,
                 const amrex::Real threshold) noexcept
  {
    using amrex::Real;
    if (!(threshold > Real(0.0))) {
      return true;  // threshold=0 reproduces the historical unconditional AFD
    }

    Real sound_scale = Real(0.0);
    for (int m = 0; m < 6; ++m) {
      sound_scale = amrex::max(
          sound_scale, std::abs(q(iv + (m - 3) * ivd, cls_t::QC)));
    }
    for (int field = 0; field < 5; ++field) {
      const int comp = field == 0 ? cls_t::QRHO
                       : field == 1 ? cls_t::QPRES
                                    : cls_t::QU + field - 2;
      Real field_scale = Real(0.0);
      for (int m = 0; m < 6; ++m) {
        field_scale = amrex::max(
            field_scale, std::abs(q(iv + (m - 3) * ivd, comp)));
      }
      const Real floor = field >= 2 ? sound_scale : field_scale;
      for (int m = 0; m < 4; ++m) {
        const Real a = q(iv + (m - 3) * ivd, comp);
        const Real b = q(iv + (m - 2) * ivd, comp);
        const Real c = q(iv + (m - 1) * ivd, comp);
        const Real denom = std::abs(a) + Real(2.0) * std::abs(b) +
                           std::abs(c) + floor +
                           std::numeric_limits<Real>::min();
        if (std::abs(a - Real(2.0) * b + c) / denom > threshold) {
          return false;
        }
      }
    }
    return true;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  rarefaction_stencil(const amrex::IntVect& iv,
                      const amrex::IntVect& ivd,
                      const amrex::Array4<const amrex::Real>& q) const noexcept
  {
    for (int m = 1; m < 5; ++m) {
      const amrex::IntVect p = iv + (m - 3) * ivd;
      if (this->rarefaction_pocket(
              q(p, cls_t::QRHO), q(p - ivd, cls_t::QRHO),
              q(p + ivd, cls_t::QRHO), q(p, cls_t::QPRES),
              q(p - ivd, cls_t::QPRES), q(p + ivd, cls_t::QPRES))) {
        return true;
      }
    }
    return false;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  primitive_to_characteristic(
      const amrex::IntVect& iv, const int dir, const int t1, const int t2,
      const amrex::Array4<const amrex::Real>& q,
      amrex::Real w[6]) noexcept
  {
    using amrex::Real;
    const Real rho = q(iv, cls_t::QRHO);
    const Real pressure = q(iv, cls_t::QPRES);
    const Real gamma = q(iv, cls_t::QG);
    const Real acoustic_scale = std::sqrt(gamma * rho * pressure);
    w[0] = rho * (Real(1.0) - Real(1.0) / gamma);
    w[1] = Real(0.5) *
           (pressure + acoustic_scale * q(iv, cls_t::QU + dir));
    w[2] = Real(0.5) *
           (pressure - acoustic_scale * q(iv, cls_t::QU + dir));
    w[3] = q(iv, cls_t::QU + t1);
    w[4] = q(iv, cls_t::QU + t2);
    w[5] = gamma;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  characteristic_to_primitive(
      const amrex::Real w[6], amrex::Real& rho, amrex::Real& un,
      amrex::Real& ut1, amrex::Real& ut2, amrex::Real& pressure,
      amrex::Real& gamma) noexcept
  {
    using amrex::Real;
    gamma = w[5];
    pressure = w[1] + w[2];
    if (!(gamma > Real(1.0)) || !(pressure > Real(0.0)) ||
        !std::isfinite(gamma + pressure)) {
      return false;
    }
    rho = w[0] / (Real(1.0) - Real(1.0) / gamma);
    if (!(rho > Real(0.0)) || !std::isfinite(rho)) {
      return false;
    }
    un = (w[1] - w[2]) / std::sqrt(gamma * rho * pressure);
    ut1 = w[3];
    ut2 = w[4];
    return std::isfinite(un + ut1 + ut2);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  cell_centered_llf(const amrex::IntVect& iv, const int dir,
                    const amrex::Array4<const amrex::Real>& q,
                    const amrex::Array4<amrex::Real>& flx,
                    const cls_t& cls) noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    const amrex::IntVect il = iv - ivd;
    Real ul[cls_t::NCONS];
    Real ur[cls_t::NCONS];
    Real fl[cls_t::NCONS];
    Real fr[cls_t::NCONS];
    cls.prims2cons(il, q, ul);
    cls.prims2cons(iv, q, ur);
    cls.prims2flux(il, dir, q, fl);
    cls.prims2flux(iv, dir, q, fr);
    const Real alpha = amrex::max(
        std::abs(q(il, cls_t::QU + dir)) + q(il, cls_t::QC),
        std::abs(q(iv, cls_t::QU + dir)) + q(iv, cls_t::QC));
    for (int n = 0; n < cls_t::NCONS; ++n) {
      flx(iv, n) = Real(0.5) * (fl[n] + fr[n] - alpha * (ur[n] - ul[n]));
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  face_flux_ibm_boundary(
      const amrex::IntVect& iv, const int dir,
      const amrex::Array4<const amrex::Real>& q,
      const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
      const amrex::Array4<std::uint8_t>& marker,
      const amrex::Real teno_cutoff) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    bool usable[6];
    for (int m = 0; m < 6; ++m) {
      usable[m] = ibm_flux::is_usable(iv + (m - 3) * ivd, marker);
    }

    const bool valid_l[3] = {
        usable[0] && usable[1] && usable[2],
        usable[1] && usable[2] && usable[3],
        usable[2] && usable[3] && usable[4]};
    const bool valid_r[3] = {
        usable[1] && usable[2] && usable[3],
        usable[2] && usable[3] && usable[4],
        usable[3] && usable[4] && usable[5]};
    if (!(valid_l[0] || valid_l[1] || valid_l[2]) ||
        !(valid_r[0] || valid_r[1] || valid_r[2])) {
      if (usable[2] && usable[3]) {
        ibm_flux::llf_flux(iv, dir, q, flx, cls);
      } else {
        ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
      }
      return;
    }

    const int t1 = dir == 0 ? 1 : 0;
    const int t2 = dir == 2 ? 1 : 2;
    Real wc[6][6] = {};
    for (int m = 0; m < 6; ++m) {
      if (usable[m]) {
        primitive_to_characteristic(iv + (m - 3) * ivd, dir, t1, t2, q,
                                    wc[m]);
      }
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
      wl[n] = PointScheme::right_masked(sl, valid_l, teno_cutoff);
      wr[n] = PointScheme::left_masked(sr, valid_r, teno_cutoff);
    }

    Real rl, ul, ut1l, ut2l, pl, gl;
    Real rr, ur, ut1r, ut2r, pr, gr;
    if (!characteristic_to_primitive(wl, rl, ul, ut1l, ut2l, pl, gl) ||
        !characteristic_to_primitive(wr, rr, ur, ut1r, ut2r, pr, gr)) {
      ibm_flux::llf_flux(iv, dir, q, flx, cls);
      return;
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

    Real fn = Real(0.0);
    Real ft1 = Real(0.0);
    Real ft2 = Real(0.0);
    Real fe = Real(0.0);
    Real fry[1] = {Real(0.0)};
    this->hllc(rl, ul, pl, ut1l, ut2l, el, yl, cl,
               rr, ur, pr, ut1r, ut2r, er, yr, cr,
               fn, ft1, ft2, fe, fry);
    if (!std::isfinite(fn + ft1 + ft2 + fe + fry[0])) {
      ibm_flux::llf_flux(iv, dir, q, flx, cls);
      return;
    }

    flx(iv, cls_t::UMX + dir) = fn;
    flx(iv, cls_t::UMX + t1) = ft1;
    flx(iv, cls_t::UMX + t2) = ft2;
    flx(iv, cls_t::UET) = fe;
    flx(iv, cls_t::URHO) = fry[0];
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  face_flux(const amrex::IntVect& iv, const int dir,
            const amrex::Array4<const amrex::Real>& q,
            const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
            const int correction, const amrex::Real smoothness_threshold,
            const amrex::Real teno_cutoff) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    if (!valid_stencil(iv, ivd, q) || rarefaction_stencil(iv, ivd, q)) {
      cell_centered_llf(iv, dir, q, flx, cls);
      return;
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
      cell_centered_llf(iv, dir, q, flx, cls);
      return;
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

    Real fn = Real(0.0);
    Real ft1 = Real(0.0);
    Real ft2 = Real(0.0);
    Real fe = Real(0.0);
    Real fry[1] = {Real(0.0)};
    this->hllc(rl, ul, pl, ut1l, ut2l, el, yl, cl,
               rr, ur, pr, ut1r, ut2r, er, yr, cr,
               fn, ft1, ft2, fe, fry);
    if (!std::isfinite(fn + ft1 + ft2 + fe + fry[0])) {
      cell_centered_llf(iv, dir, q, flx, cls);
      return;
    }

    flx(iv, cls_t::UMX + dir) = fn;
    flx(iv, cls_t::UMX + t1) = ft1;
    flx(iv, cls_t::UMX + t2) = ft2;
    flx(iv, cls_t::UET) = fe;
    flx(iv, cls_t::URHO) = fry[0];

    if (correction == 0 ||
        !smooth_stencil(iv, ivd, q, smoothness_threshold)) {
      return;
    }

    Real f[6][cls_t::NCONS];
    for (int m = 0; m < 6; ++m) {
      cls.prims2flux(iv + (m - 3) * ivd, dir, q, f[m]);
    }
    for (int n = 0; n < cls_t::NCONS; ++n) {
      flx(iv, n) +=
          -Real(1.0 / 1152.0) *
              (-Real(5.0) * f[0][n] + Real(39.0) * f[1][n] -
               Real(34.0) * f[2][n] - Real(34.0) * f[3][n] +
               Real(39.0) * f[4][n] - Real(5.0) * f[5][n]) +
          Real(7.0 / 11520.0) *
              (f[0][n] - Real(3.0) * f[1][n] + Real(2.0) * f[2][n] +
               Real(2.0) * f[3][n] - Real(3.0) * f[4][n] + f[5][n]);
    }
  }
};

template <typename cls_t>
using afd_hllc_teno5_t = afd_hllc_t<AfdReconScheme::Teno5, cls_t>;

template <typename cls_t>
using afd_hllc_wenoz5_t = afd_hllc_t<AfdReconScheme::WenoZ5, cls_t>;

#endif
