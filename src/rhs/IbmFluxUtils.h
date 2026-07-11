#ifndef IBM_FLUX_UTILS_H_
#define IBM_FLUX_UTILS_H_

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>

#include <cmath>
#include <cstdint>

namespace ibm_flux {

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_solid(const amrex::IntVect& iv, const MarkerArray& marker) noexcept
{
  return marker(iv, 0) != 0;
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_usable(const amrex::IntVect& iv, const MarkerArray& marker) noexcept
{
#if AMREX_USE_GPIBM
  return marker(iv, 0) == 0 || marker(iv, 1) != 0;
#else
  // The EB marker adapter does not currently expose the GPIBM
  // "successfully reconstructed GP" channel semantics.
  return marker(iv, 0) == 0;
#endif
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
is_solid_solid_face(const amrex::IntVect& iv,
                    const amrex::IntVect& ivd,
                    const MarkerArray& marker) noexcept
{
  return is_solid(iv, marker) && is_solid(iv - ivd, marker);
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
stencil_usable(const amrex::IntVect& face, const amrex::IntVect& ivd,
               const int first_offset, const int count,
               const MarkerArray& marker) noexcept
{
  for (int n = 0; n < count; ++n) {
    if (!is_usable(face + (first_offset + n) * ivd, marker)) {
      return false;
    }
  }
  return true;
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
stencil_all_fluid(const amrex::IntVect& face, const amrex::IntVect& ivd,
                  const int first_offset, const int count,
                  const MarkerArray& marker) noexcept
{
  for (int n = 0; n < count; ++n) {
    if (is_solid(face + (first_offset + n) * ivd, marker)) {
      return false;
    }
  }
  return true;
}

template <typename PrimArray, typename MarkerArray, typename Cls>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
adjacent_states_valid(const amrex::IntVect& face,
                      const amrex::IntVect& ivd,
                      const PrimArray& prims,
                      const MarkerArray& marker) noexcept
{
  using amrex::Real;
  const amrex::IntVect left = face - ivd;
  if (!is_usable(left, marker) || !is_usable(face, marker)) {
    return false;
  }
  const Real finite_sum = prims(left, Cls::QRHO) +
                          prims(face, Cls::QRHO) +
                          prims(left, Cls::QPRES) +
                          prims(face, Cls::QPRES) +
                          prims(left, Cls::QC) + prims(face, Cls::QC);
  return prims(left, Cls::QRHO) > Real(0.0) &&
         prims(face, Cls::QRHO) > Real(0.0) &&
         prims(left, Cls::QPRES) > Real(0.0) &&
         prims(face, Cls::QPRES) > Real(0.0) &&
         prims(left, Cls::QC) > Real(0.0) &&
         prims(face, Cls::QC) > Real(0.0) && std::isfinite(finite_sum);
}

template <typename PrimArray, typename FluxArray, typename Cls>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
two_point_central_flux(const amrex::IntVect& face, const int dir,
                       const PrimArray& prims, const FluxArray& flux,
                       const Cls& cls) noexcept
{
  using amrex::Real;
  const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
  Real left[Cls::NCONS];
  Real right[Cls::NCONS];
  cls.prims2flux(face - ivd, dir, prims, left);
  cls.prims2flux(face, dir, prims, right);
  for (int n = 0; n < Cls::NCONS; ++n) {
    flux(face, n) = Real(0.5) * (left[n] + right[n]);
  }
}

template <typename MarkerArray>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
usable_block(const amrex::IntVect& face, const amrex::IntVect& ivd,
             const int start, const int count,
             const MarkerArray& marker) noexcept
{
  return stencil_usable(face, ivd, start, count, marker);
}

// Interpolate the physical flux to a face using a one-sided block that still
// contains both face-adjacent cells.  Four points give a fourth-order face
// value; the three-point fallback gives a third-order face value.  With the
// quadratic (eorder=2) GP extension this is the minimum closure needed to keep
// the face-to-face divergence second-order where a high-order bulk stencil
// changes orientation near the immersed surface.
template <typename PrimArray, typename FluxArray, typename MarkerArray,
          typename Cls>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
one_sided_polynomial_flux(const amrex::IntVect& face, const int dir,
                          const PrimArray& prims, const FluxArray& flux,
                          const MarkerArray& marker,
                          const Cls& cls) noexcept
{
  using amrex::Real;
  const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);

  int start = 99;
  int count = 0;
  if (usable_block(face, ivd, -2, 4, marker)) {
    start = -2;
    count = 4;
  } else if (usable_block(face, ivd, -1, 4, marker)) {
    start = -1;
    count = 4;
  } else if (usable_block(face, ivd, -3, 4, marker)) {
    start = -3;
    count = 4;
  } else if (usable_block(face, ivd, -1, 3, marker)) {
    start = -1;
    count = 3;
  } else if (usable_block(face, ivd, -2, 3, marker)) {
    start = -2;
    count = 3;
  } else {
    return false;
  }

  Real weight[4] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};
  if (count == 4 && start == -2) {
    weight[0] = Real(-1.0 / 16.0);
    weight[1] = Real( 9.0 / 16.0);
    weight[2] = Real( 9.0 / 16.0);
    weight[3] = Real(-1.0 / 16.0);
  } else if (count == 4 && start == -1) {
    weight[0] = Real( 5.0 / 16.0);
    weight[1] = Real(15.0 / 16.0);
    weight[2] = Real(-5.0 / 16.0);
    weight[3] = Real( 1.0 / 16.0);
  } else if (count == 4) {
    weight[0] = Real( 1.0 / 16.0);
    weight[1] = Real(-5.0 / 16.0);
    weight[2] = Real(15.0 / 16.0);
    weight[3] = Real( 5.0 / 16.0);
  } else if (start == -1) {
    weight[0] = Real( 3.0 / 8.0);
    weight[1] = Real( 3.0 / 4.0);
    weight[2] = Real(-1.0 / 8.0);
  } else {
    weight[0] = Real(-1.0 / 8.0);
    weight[1] = Real( 3.0 / 4.0);
    weight[2] = Real( 3.0 / 8.0);
  }

  Real value[Cls::NCONS] = {};
  Real cell_flux[Cls::NCONS];
  for (int m = 0; m < count; ++m) {
    cls.prims2flux(face + (start + m) * ivd, dir, prims, cell_flux);
    for (int n = 0; n < Cls::NCONS; ++n) {
      value[n] += weight[m] * cell_flux[n];
    }
  }
  for (int n = 0; n < Cls::NCONS; ++n) {
    flux(face, n) = value[n];
  }
  return true;
}

template <typename PrimArray, typename FluxArray, typename Cls>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
llf_flux(const amrex::IntVect& face, const int dir,
         const PrimArray& prims, const FluxArray& flux,
         const Cls& cls) noexcept
{
  using amrex::Real;
  const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
  const amrex::IntVect left_iv = face - ivd;
  Real left_cons[Cls::NCONS];
  Real right_cons[Cls::NCONS];
  Real left_flux[Cls::NCONS];
  Real right_flux[Cls::NCONS];
  cls.prims2cons(left_iv, prims, left_cons);
  cls.prims2cons(face, prims, right_cons);
  cls.prims2flux(left_iv, dir, prims, left_flux);
  cls.prims2flux(face, dir, prims, right_flux);
  const Real alpha = amrex::max(
      std::abs(prims(left_iv, Cls::QU + dir)) + prims(left_iv, Cls::QC),
      std::abs(prims(face, Cls::QU + dir)) + prims(face, Cls::QC));
  for (int n = 0; n < Cls::NCONS; ++n) {
    flux(face, n) = Real(0.5) *
                    (left_flux[n] + right_flux[n] -
                     alpha * (right_cons[n] - left_cons[n]));
  }
}

template <typename FluxArray, typename Cls>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
zero_flux(const amrex::IntVect& face, const FluxArray& flux) noexcept
{
  for (int n = 0; n < Cls::NCONS; ++n) {
    flux(face, n) = amrex::Real(0.0);
  }
}

}  // namespace ibm_flux

#endif
