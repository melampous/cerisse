#ifndef DIFFUSION_CD_H_
#define DIFFUSION_CD_H_

#include <AMReX_FArrayBox.H>

#include <array>
#include <cstdint>

#include "IBMSharedGPFluxUtils.h"
#include "diff_ops.H"

// Heat conduction only.  This class contributes the conductive energy flux
// -lambda grad(T); compute_rhs owns the subsequent flux divergence.
template <typename param, typename cls_t>
class diffusiveheat_t {
 public:
  static_assert(param::order == 2 || param::order == 4 || param::order == 6,
                "diffusiveheat_t supports orders 2, 4, and 6 only");

  static constexpr int halfsten = param::order / 2;
  // The conductive metric flux vanishes at the axis; do not overwrite an
  // inviscid metric h-flux stored in the shared radial face slot.
  static constexpr bool rz_radial_axis_face_addition_is_zero = true;

  AMREX_GPU_HOST_DEVICE constexpr diffusiveheat_t() = default;
  AMREX_GPU_HOST_DEVICE ~diffusiveheat_t() = default;

  AMREX_GPU_HOST void init_coeffs() {}

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void dflux_ibm(
      const Geometry& geom, const MFIter& mfi, const Array4<Real>& prims,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const Array4<Real>& /*cons*/, const cls_t* cls,
      const Array4<uint8_t>& ibMarkers) {
    compute_fluxes(geom, mfi, prims, flxt, cls, ibMarkers);
  }
#else
  void dflux(const Geometry& geom, const MFIter& mfi,
             const Array4<Real>& prims,
             std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
             const Array4<Real>& /*cons*/, const cls_t* cls) {
    compute_fluxes(geom, mfi, prims, flxt, cls);
  }
#endif

  // Heat conduction has no additional cylindrical hoop-stress source.
  void rz_geometric_source(const Geometry& /*geom*/, const MFIter& /*mfi*/,
                           const Array4<Real>& /*prims*/,
                           const Array4<Real>& /*state*/,
                           const cls_t* /*cls*/) {}

 public:
  // NVCC requires a member enclosing an extended device lambda to be public.
  // Keep this as the shared implementation for dflux/dflux_ibm.
#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void compute_fluxes(
      const Geometry& geom, const MFIter& mfi, const Array4<Real>& prims,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const cls_t* cls, const Array4<uint8_t>& ibMarkers) {
#else
  void compute_fluxes(
      const Geometry& geom, const MFIter& mfi, const Array4<Real>& prims,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const cls_t* cls) {
#endif
    const auto dxinv = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const bool is_rz = geom.IsRZ();
    const Box bxg = mfi.growntilebox(cls->NGHOST);

    // Conductivity is reused by every directional face stencil.  Cache only
    // lambda (the legacy implementation also computed an unused viscosity).
    FArrayBox conductivity(bxg, 1, The_Async_Arena());
    const auto lam = conductivity.array();
    amrex::ParallelFor(
        bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if (AMREX_USE_GPIBM || CNS_USE_EB)
          const IntVect iv{AMREX_D_DECL(i, j, k)};
          if (!ibm_flux::is_usable(iv, ibMarkers)) {
            lam(i, j, k) = Real(0.0);
            return;
          }
#endif
          lam(i, j, k) = cls->cond(prims(i, j, k, cls_t::QT));
        });

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const auto flx = flxt[dir]->array();
      const Box bxface = mfi.grownnodaltilebox(dir, 0);
      const IntVect ivd = IntVect::TheDimensionVector(dir);

      amrex::ParallelFor(
          bxface, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const IntVect iv{AMREX_D_DECL(i, j, k)};
            if (is_rz && dir == 0 &&
                !(prob_lo[0] + Real(i) * dx[0] > Real(0.0))) {
              return;
            }
#if (AMREX_USE_GPIBM || CNS_USE_EB)
            if (ibm_flux::is_solid_solid_face(iv, ivd, ibMarkers)) {
              return;
            }

            if (!ibm_flux::stencil_usable(iv, ivd, -halfsten,
                                          param::order, ibMarkers)) {
              if (!ibm_flux::is_usable(iv - ivd, ibMarkers) ||
                  !ibm_flux::is_usable(iv, ibMarkers)) {
                return;
              }
              const Real lamf = Real(0.5) *
                  (lam(iv - ivd) + lam(iv));
              const Real dTdn =
                  (prims(iv, cls_t::QT) - prims(iv - ivd, cls_t::QT)) *
                  dxinv[dir];
              flx(iv, cls_t::UET) -= lamf * dTdn;
              return;
            }
#endif
            const Real lamf = interp<param::order>(iv, dir, 0, lam);
            const Real dTdn =
                normal_diff<param::order>(iv, dir, cls_t::QT, prims, dxinv);
            flx(iv, cls_t::UET) -= lamf * dTdn;
          });
    }
  }
};

#endif
