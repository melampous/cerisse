#ifndef CERISSE_EULER_ADMISSIBLE_PROLONGATION_H_
#define CERISSE_EULER_ADMISSIBLE_PROLONGATION_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Interpolater.H>

#include <cmath>
#include <limits>

namespace cerisse::amr {

struct EulerStateLayout {
  int momentum_x = -1;
  int momentum_y = -1;
  int momentum_z = -1;
  int total_energy = -1;
  int density = -1;
  amrex::Real density_floor = amrex::Real(0.0);
  amrex::Real internal_energy_density_floor = amrex::Real(0.0);
  amrex::Real specific_internal_energy_floor = amrex::Real(0.0);
};

// AMReX's limited conservative linear interpolater bounds each conservative
// component, but that does not bound rho*e = rho*E-|rho*u|^2/(2*rho).  This
// mapper first constructs the ordinary limited-linear children.  If any child
// leaves the configured admissible set, it applies one common theta to every
// conservative child deviation in that parent block.  The set always enforces
// positive density and internal-energy density and can additionally enforce a
// specific-internal-energy floor used by the thermodynamic closure.  The
// common theta keeps the coarse-cell average unchanged in Cartesian and
// metric-weighted R-Z interpolation.  No state is clipped independently.
class EulerAdmissibleConservativeLinear final
    : public amrex::CellConservativeLinear {
 public:
  explicit EulerAdmissibleConservativeLinear(EulerStateLayout layout)
      : amrex::CellConservativeLinear(true), layout_(layout) {
    if (layout_.momentum_x < 0 || layout_.momentum_y < 0 ||
        layout_.momentum_z < 0 || layout_.total_energy < 0 ||
        layout_.density < 0 || !(layout_.density_floor > amrex::Real(0.0)) ||
        !(layout_.internal_energy_density_floor > amrex::Real(0.0)) ||
        !(layout_.specific_internal_energy_floor >= amrex::Real(0.0)) ||
        !std::isfinite(layout_.density_floor) ||
        !std::isfinite(layout_.internal_energy_density_floor) ||
        !std::isfinite(layout_.specific_internal_energy_floor)) {
      amrex::Abort(
          "EulerAdmissibleConservativeLinear received an invalid state "
          "layout or admissibility floor");
    }
  }

  void interp(const amrex::FArrayBox& crse, int crse_comp,
              amrex::FArrayBox& fine, int fine_comp, int ncomp,
              const amrex::Box& fine_region, const amrex::IntVect& ratio,
              const amrex::Geometry& crse_geom,
              const amrex::Geometry& fine_geom,
              amrex::Vector<amrex::BCRec> const& bcr, int actual_comp,
              int actual_state, amrex::RunOn runon) override {
    using namespace amrex;

    const int required_components[] = {
        layout_.momentum_x, layout_.momentum_y, layout_.momentum_z,
        layout_.total_energy, layout_.density};
    for (const int component : required_components) {
      if (component < actual_comp || component >= actual_comp + ncomp) {
        amrex::Abort(
            "Euler-admissible AMR prolongation requires momentum, total "
            "energy, and density in one interpolation component group");
      }
    }
    if (!ratio.allGE(IntVect::TheUnitVector())) {
      amrex::Abort(
          "Euler-admissible AMR prolongation requires positive refinement "
          "ratios");
    }

    const Box target = fine_region & fine.box();
    if (!target.ok()) return;
    const Box parent_region = amrex::coarsen(target, ratio);
    const Box complete_children = amrex::refine(parent_region, ratio);

    FArrayBox high_order(complete_children, ncomp, The_Async_Arena());
    CellConservativeLinear::interp(
        crse, crse_comp, high_order, 0, ncomp, complete_children, ratio,
        crse_geom, fine_geom, bcr, actual_comp, actual_state, runon);
    FArrayBox parent_theta(parent_region, 1, The_Async_Arena());

    const auto coarse = crse.const_array();
    const auto high = high_order.const_array();
    const auto theta = parent_theta.array();
    const int rho_n = layout_.density - actual_comp;
    const int mx_n = layout_.momentum_x - actual_comp;
    const int my_n = layout_.momentum_y - actual_comp;
    const int mz_n = layout_.momentum_z - actual_comp;
    const int energy_n = layout_.total_energy - actual_comp;
    const Real rho_floor = layout_.density_floor;
    const Real rhoe_floor = layout_.internal_energy_density_floor;
    const Real ei_floor = layout_.specific_internal_energy_floor;
    const int ratio_x = ratio[0];
#if (AMREX_SPACEDIM >= 2)
    const int ratio_y = ratio[1];
#else
    const int ratio_y = 1;
#endif
#if (AMREX_SPACEDIM == 3)
    const int ratio_z = ratio[2];
#else
    const int ratio_z = 1;
#endif

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(
        runon, parent_region, ic, jc, kc, {
          const Real parent_rho =
              coarse(ic, jc, kc, crse_comp + rho_n);
          const Real parent_mx = coarse(ic, jc, kc, crse_comp + mx_n);
          const Real parent_my = coarse(ic, jc, kc, crse_comp + my_n);
          const Real parent_mz = coarse(ic, jc, kc, crse_comp + mz_n);
          const Real parent_energy =
              coarse(ic, jc, kc, crse_comp + energy_n);
          const bool parent_finite =
              std::isfinite(parent_rho) && std::isfinite(parent_mx) &&
              std::isfinite(parent_my) && std::isfinite(parent_mz) &&
              std::isfinite(parent_energy);
          Real parent_rhoe = -std::numeric_limits<Real>::infinity();
          if (parent_finite && parent_rho > rho_floor) {
            parent_rhoe =
                parent_energy -
                Real(0.5) *
                    (parent_mx * parent_mx + parent_my * parent_my +
                     parent_mz * parent_mz) /
                    parent_rho;
          }
          const Real parent_rhoe_floor =
              amrex::max(rhoe_floor, parent_rho * ei_floor);
          if (!std::isfinite(parent_rhoe) ||
              !(parent_rhoe > parent_rhoe_floor)) {
            // An invalid accepted coarse state is outside this mapper's
            // repair contract.  Emit NaNs so the next state audit fails
            // closed; do not hide it by manufacturing a replacement state.
            theta(ic, jc, kc, 0) =
                std::numeric_limits<Real>::quiet_NaN();
            return;
          }

          Real block_theta = Real(1.0);
          for (int oz = 0; oz < ratio_z; ++oz) {
            for (int oy = 0; oy < ratio_y; ++oy) {
              for (int ox = 0; ox < ratio_x; ++ox) {
                const int i = ic * ratio_x + ox;
                const int j = jc * ratio_y + oy;
                const int k = kc * ratio_z + oz;
                const Real child_rho = high(i, j, k, rho_n);
                const Real child_mx = high(i, j, k, mx_n);
                const Real child_my = high(i, j, k, my_n);
                const Real child_mz = high(i, j, k, mz_n);
                const Real child_energy = high(i, j, k, energy_n);
                const bool child_finite =
                    std::isfinite(child_rho) && std::isfinite(child_mx) &&
                    std::isfinite(child_my) && std::isfinite(child_mz) &&
                    std::isfinite(child_energy);
                Real child_rhoe =
                    -std::numeric_limits<Real>::infinity();
                if (child_finite && child_rho > rho_floor) {
                  child_rhoe =
                      child_energy -
                      Real(0.5) *
                          (child_mx * child_mx + child_my * child_my +
                           child_mz * child_mz) /
                          child_rho;
                }
                const Real child_rhoe_floor =
                    amrex::max(rhoe_floor, child_rho * ei_floor);
                if (std::isfinite(child_rhoe) &&
                    child_rhoe > child_rhoe_floor) {
                  continue;
                }

                if (!child_finite) {
                  block_theta = Real(0.0);
                  continue;
                }
                Real lower = Real(0.0);
                Real upper = Real(1.0);
                for (int iteration = 0; iteration < 64; ++iteration) {
                  const Real midpoint = Real(0.5) * (lower + upper);
                  const Real rho_mid =
                      parent_rho + midpoint * (child_rho - parent_rho);
                  const Real mx_mid =
                      parent_mx + midpoint * (child_mx - parent_mx);
                  const Real my_mid =
                      parent_my + midpoint * (child_my - parent_my);
                  const Real mz_mid =
                      parent_mz + midpoint * (child_mz - parent_mz);
                  const Real energy_mid = parent_energy +
                      midpoint * (child_energy - parent_energy);
                  Real rhoe_mid =
                      -std::numeric_limits<Real>::infinity();
                  if (rho_mid > rho_floor) {
                    rhoe_mid =
                        energy_mid -
                        Real(0.5) *
                            (mx_mid * mx_mid + my_mid * my_mid +
                             mz_mid * mz_mid) /
                            rho_mid;
                  }
                  const Real rhoe_mid_floor =
                      amrex::max(rhoe_floor, rho_mid * ei_floor);
                  if (std::isfinite(rhoe_mid) &&
                      rhoe_mid > rhoe_mid_floor) {
                    lower = midpoint;
                  } else {
                    upper = midpoint;
                  }
                }
                block_theta = amrex::min(block_theta, lower);
              }
            }
          }
          if (block_theta < Real(1.0)) {
            block_theta *= Real(1.0) - Real(1.0e-10);
          }
          theta(ic, jc, kc, 0) = block_theta;
        });

    const auto destination = fine.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        runon, target, ncomp, i, j, k, n, {
          const int ic = amrex::coarsen(i, ratio_x);
#if (AMREX_SPACEDIM >= 2)
          const int jc = amrex::coarsen(j, ratio_y);
#else
          const int jc = 0;
#endif
#if (AMREX_SPACEDIM == 3)
          const int kc = amrex::coarsen(k, ratio_z);
#else
          const int kc = 0;
#endif
          const Real block_theta = theta(ic, jc, kc, 0);
          const Real parent = coarse(ic, jc, kc, crse_comp + n);
          if (block_theta == Real(1.0)) {
            destination(i, j, k, fine_comp + n) = high(i, j, k, n);
          } else if (block_theta == Real(0.0)) {
            destination(i, j, k, fine_comp + n) = parent;
          } else {
            destination(i, j, k, fine_comp + n) =
                parent + block_theta * (high(i, j, k, n) - parent);
          }
        });
  }

 private:
  EulerStateLayout layout_;
};

}  // namespace cerisse::amr

#endif
