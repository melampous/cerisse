#ifndef CERISSE_IBM_POSITIVITY_LIMITER_H_
#define CERISSE_IBM_POSITIVITY_LIMITER_H_

// Full-cell, shared-Cartesian-face convex limiting for the frozen pure GP-IBM
// production path. Keep validation-only traces and analytic oracles outside
// this implementation.

#include <AMReX_FabArrayUtility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_iMultiFab.H>
#include <CNSconstants.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <prob.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <type_traits>
#include <vector>

namespace IBM::positivity {

struct StageAdmissibilityThresholds {
  amrex::Real rho = CNSConstants::smallr;
  amrex::Real pressure_rhoe = CNSConstants::min_press();
  amrex::Real temperature_ei = amrex::Real(0.0);

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  rhoe_min(const amrex::Real density) const noexcept {
    return amrex::max(pressure_rhoe,
                      density * temperature_ei);
  }
};

struct StageAdmissibilityStats {
  amrex::Long active = 0;
  amrex::Long nonfinite = 0;
  amrex::Long invalid_rho = 0;
  amrex::Long invalid_rhoe = 0;
  amrex::Real min_rho = std::numeric_limits<amrex::Real>::max();
  amrex::Real min_rhoe = std::numeric_limits<amrex::Real>::max();

  [[nodiscard]] bool admissible() const noexcept {
    return nonfinite == 0 && invalid_rho == 0 && invalid_rhoe == 0;
  }
};

template <typename Closure>
inline StageAdmissibilityThresholds stageAdmissibilityThresholds(
    const Closure& closure) {
  if (!(closure.gamma_m1 > amrex::Real(0.0)) ||
      !std::isfinite(closure.gamma_m1)) {
    amrex::Abort(
        "IBM positivity limiter requires a constant-gamma ideal gas");
  }
  const amrex::Real pressure_rhoe =
      CNSConstants::min_press() / closure.gamma_m1;
  // The Euler admissible set is defined by positive density and positive
  // internal-energy density. CLIP_TEMPERATURE_MIN is a primitive-conversion
  // safeguard; treating it as an invariant-domain constraint activates the
  // conservative flux limiter for otherwise admissible low-temperature
  // states and changes the shock solution. Temperature clipping must not
  // enter the positivity proof or the face-limiting coefficient.
  return {CNSConstants::smallr, pressure_rhoe, amrex::Real(0.0)};
}

template <typename Closure>
inline StageAdmissibilityStats collectStageAdmissibility(
    const amrex::MultiFab& state, int level, const Closure& closure) {
  using namespace amrex;
  const auto thresholds = stageAdmissibilityThresholds(closure);
  const auto& marker_mf = *IBM::ib.bmf_a[level];

  ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,
            ReduceOpMin, ReduceOpMin>
      reduce_ops;
  ReduceData<Long, Long, Long, Long, Real, Real> reduce_data(reduce_ops);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (MFIter mfi(state, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const auto state_a = state.const_array(mfi);
    const auto marker = marker_mf.const_array(mfi);
    reduce_ops.eval(
        bx, reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
          if (marker(i, j, k, 0) != 0) {
            return {Long(0), Long(0), Long(0), Long(0),
                    std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::max()};
          }

          bool finite = true;
          for (int n = 0; n < Closure::NCONS; ++n) {
            finite = finite && std::isfinite(state_a(i, j, k, n));
          }
          const Real rho = state_a(i, j, k, Closure::URHO);
          const Real mx = state_a(i, j, k, Closure::UMX);
          const Real my = state_a(i, j, k, Closure::UMY);
#if (AMREX_SPACEDIM == 3)
          const Real mz = state_a(i, j, k, Closure::UMZ);
#else
          const Real mz = Real(0.0);
#endif
          const Real energy = state_a(i, j, k, Closure::UET);
          Real rhoe = std::numeric_limits<Real>::quiet_NaN();
          if (finite && rho > thresholds.rho) {
            rhoe = energy -
                Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
          }
          const bool finite_rhoe = std::isfinite(rhoe);
          return {
              Long(1),
              finite && finite_rhoe ? Long(0) : Long(1),
              finite && rho > thresholds.rho ? Long(0) : Long(1),
              finite && finite_rhoe && rhoe > thresholds.rhoe_min(rho)
                  ? Long(0)
                  : Long(1),
              std::isfinite(rho) ? rho
                                 : std::numeric_limits<Real>::max(),
              finite_rhoe ? rhoe : std::numeric_limits<Real>::max()};
        });
  }

  const auto reduced = reduce_data.value(reduce_ops);
  StageAdmissibilityStats stats;
  stats.active = amrex::get<0>(reduced);
  stats.nonfinite = amrex::get<1>(reduced);
  stats.invalid_rho = amrex::get<2>(reduced);
  stats.invalid_rhoe = amrex::get<3>(reduced);
  stats.min_rho = amrex::get<4>(reduced);
  stats.min_rhoe = amrex::get<5>(reduced);
  ParallelDescriptor::ReduceLongSum(stats.active);
  ParallelDescriptor::ReduceLongSum(stats.nonfinite);
  ParallelDescriptor::ReduceLongSum(stats.invalid_rho);
  ParallelDescriptor::ReduceLongSum(stats.invalid_rhoe);
  ParallelDescriptor::ReduceRealMin(stats.min_rho);
  ParallelDescriptor::ReduceRealMin(stats.min_rhoe);
  return stats;
}

inline void printStageAdmissibility(
    const char* label, const StageAdmissibilityStats& stats) {
  amrex::Print() << std::setprecision(17)
                 << "[IBM-Stage-Admissibility-Oracle] " << label
                 << " active=" << stats.active
                 << " nonfinite=" << stats.nonfinite
                 << " invalid_rho=" << stats.invalid_rho
                 << " invalid_rhoe=" << stats.invalid_rhoe
                 << " min_rho=" << stats.min_rho
                 << " min_rhoe=" << stats.min_rhoe << '\n';
}

// Failure-only, read-only detail for an inadmissible SSPRK Forward-Euler
// input.  This deliberately consumes the same active-fluid definition and
// thresholds as collectStageAdmissibility().  It neither fills ghosts nor
// changes a conservative state, marker, flux, or AMR data structure.
//
// ``finer_coverage`` is the next level's BoxArray coarsened to this level.  The
// optional regrid context labels valid fine cells as old-fine overlap or
// coarse-prolongated and preserves the exact parent source stencil.  These are
// read-only diagnostic copies; they never enter FillPatch, the GP extension,
// flux reconstruction, or an accepted state update.
template <typename Closure>
inline void printStageInputBadCells(
    const amrex::MultiFab& state, const int level, const Closure& closure,
    const amrex::Geometry& geom, const amrex::BoxArray* finer_coverage,
    const amrex::iMultiFab* regrid_new_from_coarse,
    const amrex::MultiFab* regrid_coarse_source,
    const amrex::IntVect& regrid_ref_ratio,
    const bool new_grid_context, const char* scheme, const int stage_index,
    const int stage_count, const amrex::Real stage_time,
    const int maximum_to_print = 64) {
  using namespace amrex;
  if (maximum_to_print <= 0) return;

  const auto thresholds = stageAdmissibilityThresholds(closure);
  const auto& marker_mf = *IBM::ib.bmf_a[level];
  if (state.boxArray() != marker_mf.boxArray() ||
      state.DistributionMap() != marker_mf.DistributionMap()) {
    amrex::Abort(
        "FE-input bad-cell diagnostic requires state/marker layout identity");
  }
  if (regrid_new_from_coarse != nullptr &&
      (state.boxArray() != regrid_new_from_coarse->boxArray() ||
       state.DistributionMap() !=
           regrid_new_from_coarse->DistributionMap())) {
    amrex::Abort(
        "regrid provenance diagnostic requires state/mask layout identity");
  }
  if ((regrid_new_from_coarse == nullptr) !=
      (regrid_coarse_source == nullptr)) {
    amrex::Abort(
        "regrid provenance diagnostic requires both mask and coarse source");
  }
  if (regrid_coarse_source != nullptr) {
    BoxArray expected_coarse_boxes = state.boxArray();
    expected_coarse_boxes.coarsen(regrid_ref_ratio);
    if (expected_coarse_boxes != regrid_coarse_source->boxArray() ||
        state.DistributionMap() !=
            regrid_coarse_source->DistributionMap() ||
        regrid_coarse_source->nGrow() < 1) {
      amrex::Abort(
          "regrid coarse-source diagnostic layout/ref-ratio mismatch");
    }
  }

  // Failure reporting is host-side.  Explicit copies keep this diagnostic
  // valid when the normal arenas are device-only rather than managed.
  MultiFab state_host(
      state.boxArray(), state.DistributionMap(), state.nComp(),
      state.nGrowVect(), MFInfo().SetArena(The_Pinned_Arena()));
  IBMultiFab<std::uint8_t> marker_host(
      marker_mf.boxArray(), marker_mf.DistributionMap(), marker_mf.nComp(),
      marker_mf.nGrow(), MFInfo().SetArena(The_Pinned_Arena()));
  std::unique_ptr<iMultiFab> provenance_host;
  std::unique_ptr<MultiFab> coarse_source_host;
  if (regrid_new_from_coarse != nullptr) {
    provenance_host = std::make_unique<iMultiFab>(
        regrid_new_from_coarse->boxArray(),
        regrid_new_from_coarse->DistributionMap(), 1, 0,
        MFInfo().SetArena(The_Pinned_Arena()));
    coarse_source_host = std::make_unique<MultiFab>(
        regrid_coarse_source->boxArray(),
        regrid_coarse_source->DistributionMap(),
        regrid_coarse_source->nComp(), regrid_coarse_source->nGrowVect(),
        MFInfo().SetArena(The_Pinned_Arena()));
  }
  amrex::dtoh_memcpy(state_host, state);
  amrex::dtoh_memcpy(marker_host, marker_mf);
  if (provenance_host != nullptr) {
    amrex::dtoh_memcpy(*provenance_host, *regrid_new_from_coarse);
    amrex::dtoh_memcpy(*coarse_source_host, *regrid_coarse_source);
  }
  Gpu::streamSynchronize();

  const auto is_bad_active_fluid =
      [&](const auto& state_a, const auto& marker, const int i, const int j,
          const int k) {
        if (marker(i, j, k, 0) != 0) return false;
        bool finite = true;
        for (int n = 0; n < Closure::NCONS; ++n) {
          finite = finite && std::isfinite(state_a(i, j, k, n));
        }
        const Real rho = state_a(i, j, k, Closure::URHO);
        const Real mx = state_a(i, j, k, Closure::UMX);
        const Real my = state_a(i, j, k, Closure::UMY);
#if (AMREX_SPACEDIM == 3)
        const Real mz = state_a(i, j, k, Closure::UMZ);
#else
        const Real mz = Real(0.0);
#endif
        Real oracle_rhoe = std::numeric_limits<Real>::quiet_NaN();
        if (finite && rho > thresholds.rho) {
          oracle_rhoe = state_a(i, j, k, Closure::UET) -
                        Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
        }
        return !finite || !(rho > thresholds.rho) ||
               !std::isfinite(oracle_rhoe) ||
               !(oracle_rhoe > thresholds.rhoe_min(rho));
      };

  // Count every bad cell and select a bounded set for detailed printing.
  // Rank order only affects which cells are printed if the global cap is hit.
  int local_candidates = 0;
  Long local_bad_total = 0;
  Long local_bad_new_from_coarse = 0;
  Long local_bad_old_overlap = 0;
  Long local_bad_unknown_origin = 0;
  for (MFIter mfi(state_host, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    const auto state_a = state_host.const_array(mfi);
    const auto marker = marker_host.const_array(mfi);
    const bool have_provenance = provenance_host != nullptr;
    Array4<const int> provenance;
    if (have_provenance) {
      provenance = provenance_host->const_array(mfi);
    }
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          if (is_bad_active_fluid(state_a, marker, i, j, k)) {
            ++local_bad_total;
            if (!have_provenance) {
              ++local_bad_unknown_origin;
            } else if (provenance(i, j, k, 0) == 1) {
              ++local_bad_new_from_coarse;
            } else if (provenance(i, j, k, 0) == 0) {
              ++local_bad_old_overlap;
            } else {
              ++local_bad_unknown_origin;
            }
            if (local_candidates < maximum_to_print) {
              ++local_candidates;
            }
          }
        }
      }
    }
  }

  const int my_rank = ParallelDescriptor::MyProc();
  const int rank_count = ParallelDescriptor::NProcs();
  std::vector<int> candidates_by_rank(rank_count, 0);
  candidates_by_rank[my_rank] = local_candidates;
  ParallelDescriptor::ReduceIntSum(candidates_by_rank.data(), rank_count);
  int global_offset = 0;
  for (int rank = 0; rank < my_rank; ++rank) {
    global_offset += candidates_by_rank[rank];
  }
  const int local_allowance =
      amrex::max(0, amrex::min(local_candidates,
                               maximum_to_print - global_offset));

  Long bad_total = local_bad_total;
  Long bad_new_from_coarse = local_bad_new_from_coarse;
  Long bad_old_overlap = local_bad_old_overlap;
  Long bad_unknown_origin = local_bad_unknown_origin;
  ParallelDescriptor::ReduceLongSum(bad_total);
  ParallelDescriptor::ReduceLongSum(bad_new_from_coarse);
  ParallelDescriptor::ReduceLongSum(bad_old_overlap);
  ParallelDescriptor::ReduceLongSum(bad_unknown_origin);

  const Box& domain = geom.Domain();
  const BoxArray& level_boxes = state.boxArray();
  const auto dx = geom.CellSizeArray();
  const auto prob_lo = geom.ProbLoArray();

  const auto marker_class = [](const int solid_marker,
                               const int gp_marker) noexcept {
    if (solid_marker == 0) return "fluid";
    return gp_marker != 0 ? "GP" : "solid";
  };
  const auto covered_by_finer = [&](const IntVect& iv) {
    return finer_coverage != nullptr && finer_coverage->contains(iv);
  };
  const auto coarse_fine_adjacent = [&](const IntVect& iv) {
    const bool center_covered = covered_by_finer(iv);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      for (int side : {-1, 1}) {
        IntVect neighbor = iv;
        neighbor[dir] += side;
        if (!domain.contains(neighbor)) continue;
        // On a fine level, a same-level hole is filled from the coarse level.
        if (!level_boxes.contains(neighbor)) return true;
        // On a coarse level, a change in fine coverage is the other side of
        // the same coarse-fine interface.
        if (finer_coverage != nullptr &&
            covered_by_finer(neighbor) != center_covered) {
          return true;
        }
      }
    }
    return false;
  };

  amrex::Print()
      << "[IBM-FE-INPUT-BAD-CELLS] scheme=" << scheme
      << " FE_bracket=" << stage_index << '/' << stage_count
      << " level=" << level << " t=" << std::setprecision(17) << stage_time
      << " global_print_cap=" << maximum_to_print
      << " new_grid_context=" << (new_grid_context ? "YES" : "no")
      << " provenance_retained="
      << (provenance_host != nullptr ? "YES" : "no")
      << " bad_total=" << bad_total
      << " bad_new_from_coarse=" << bad_new_from_coarse
      << " bad_old_fine_overlap=" << bad_old_overlap
      << " bad_unknown_origin=" << bad_unknown_origin << '\n';

  int local_printed = 0;
  for (MFIter mfi(state_host, false);
       mfi.isValid() && local_printed < local_allowance; ++mfi) {
    const Box& bx = mfi.validbox();
    const Box& state_fab_box = state_host[mfi].box();
    const Box& marker_fab_box = marker_host[mfi].box();
    const auto state_a = state_host.const_array(mfi);
    const auto marker = marker_host.const_array(mfi);
    const bool have_regrid_context = provenance_host != nullptr;
    Array4<const int> provenance;
    Array4<const Real> coarse_source;
    const Box* provenance_fab_box = nullptr;
    const Box* coarse_source_fab_box = nullptr;
    if (have_regrid_context) {
      provenance = provenance_host->const_array(mfi);
      coarse_source = coarse_source_host->const_array(mfi);
      provenance_fab_box = &(*provenance_host)[mfi].box();
      coarse_source_fab_box = &(*coarse_source_host)[mfi].box();
    }
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z && local_printed < local_allowance; ++k) {
      for (int j = lo.y; j <= hi.y && local_printed < local_allowance; ++j) {
        for (int i = lo.x; i <= hi.x && local_printed < local_allowance; ++i) {
          if (!is_bad_active_fluid(state_a, marker, i, j, k)) continue;

          const IntVect center(AMREX_D_DECL(i, j, k));
          IntVect parent_index = center;
          parent_index.coarsen(regrid_ref_ratio);
          const int marker0 = static_cast<int>(marker(i, j, k, 0));
          const int marker1 = static_cast<int>(marker(i, j, k, 1));
          const int provenance_value =
              have_regrid_context ? provenance(i, j, k, 0) : -1;
          const char* provenance_name =
              provenance_value == 1
                  ? "new-from-coarse"
                  : (provenance_value == 0 ? "old-fine-overlap" : "unknown");
          std::ostringstream report;
          report << std::setprecision(17)
                 << "[IBM-FE-INPUT-BAD-CELL #"
                 << global_offset + local_printed << "] rank=" << my_rank
                 << " level=" << level << " index=" << center << " coord=("
                 << prob_lo[0] + (Real(i) + Real(0.5)) * dx[0]
#if (AMREX_SPACEDIM >= 2)
                 << ',' << prob_lo[1] + (Real(j) + Real(0.5)) * dx[1]
#endif
#if (AMREX_SPACEDIM == 3)
                 << ',' << prob_lo[2] + (Real(k) + Real(0.5)) * dx[2]
#endif
                 << ") coord_system=" << (geom.IsRZ() ? "RZ" : "Cartesian")
                 << " marker=(" << marker0 << ',' << marker1 << ") class="
                 << marker_class(marker0, marker1)
                 << " coarse_fine="
                 << (coarse_fine_adjacent(center) ? "YES" : "no")
                 << " covered_by_finer="
                 << (covered_by_finer(center) ? "YES" : "no")
                 << " new_grid_context="
                 << (new_grid_context ? "YES" : "no")
                 << " provenance=" << provenance_name
                 << " parent_index=" << parent_index << '\n';

          const auto append_state = [&](const auto& values, const int ii,
                                        const int jj, const int kk) {
            report << " U=[";
            for (int n = 0; n < Closure::NCONS; ++n) {
              if (n != 0) report << ',';
              report << values(ii, jj, kk, n);
            }
            report << ']';
            const Real rho = values(ii, jj, kk, Closure::URHO);
            const Real mx = values(ii, jj, kk, Closure::UMX);
            const Real my = values(ii, jj, kk, Closure::UMY);
#if (AMREX_SPACEDIM == 3)
            const Real mz = values(ii, jj, kk, Closure::UMZ);
#else
            const Real mz = Real(0.0);
#endif
            const Real energy = values(ii, jj, kk, Closure::UET);
            Real rhoe = std::numeric_limits<Real>::quiet_NaN();
            if (std::isfinite(rho) && std::isfinite(mx) &&
                std::isfinite(my) && std::isfinite(mz) &&
                std::isfinite(energy) && rho != Real(0.0)) {
              rhoe = energy -
                     Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
            }
            const Real pressure = closure.gamma_m1 * rhoe;
            const Real temperature = rhoe / (rho * closure.cv);
            report << " rho=" << rho << " rhoe=" << rhoe
                   << " p=" << pressure << " T=" << temperature;
          };

          report << "  center:";
          append_state(state_a, i, j, k);
          report << '\n';

          if (have_regrid_context &&
              coarse_source_fab_box->contains(parent_index)) {
            const int pi = parent_index[0];
#if (AMREX_SPACEDIM >= 2)
            const int pj = parent_index[1];
#else
            const int pj = 0;
#endif
#if (AMREX_SPACEDIM == 3)
            const int pk = parent_index[2];
#else
            const int pk = 0;
#endif
            report << "  coarse_parent index=" << parent_index << ':';
            append_state(coarse_source, pi, pj, pk);
            report << '\n';
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              for (int side : {-1, 1}) {
                IntVect coarse_neighbor = parent_index;
                coarse_neighbor[dir] += side;
                report << "  coarse_neighbor[d" << dir
                       << (side < 0 ? "-" : "+")
                       << "] index=" << coarse_neighbor;
                if (!coarse_source_fab_box->contains(coarse_neighbor)) {
                  report << " data=unavailable\n";
                  continue;
                }
                const int ci = coarse_neighbor[0];
#if (AMREX_SPACEDIM >= 2)
                const int cj = coarse_neighbor[1];
#else
                const int cj = 0;
#endif
#if (AMREX_SPACEDIM == 3)
                const int ck = coarse_neighbor[2];
#else
                const int ck = 0;
#endif
                append_state(coarse_source, ci, cj, ck);
                report << '\n';
              }
            }

            const int rr_x = regrid_ref_ratio[0];
#if (AMREX_SPACEDIM >= 2)
            const int rr_y = regrid_ref_ratio[1];
#else
            const int rr_y = 1;
#endif
#if (AMREX_SPACEDIM == 3)
            const int rr_z = regrid_ref_ratio[2];
#else
            const int rr_z = 1;
#endif
            for (int oz = 0; oz < rr_z; ++oz) {
              for (int oy = 0; oy < rr_y; ++oy) {
                for (int ox = 0; ox < rr_x; ++ox) {
                  const IntVect child(AMREX_D_DECL(
                      pi * rr_x + ox, pj * rr_y + oy, pk * rr_z + oz));
                  report << "  fine_child[offset=(" << ox;
#if (AMREX_SPACEDIM >= 2)
                  report << ',' << oy;
#endif
#if (AMREX_SPACEDIM == 3)
                  report << ',' << oz;
#endif
                  report << ")] index=" << child;
                  const bool state_available = state_fab_box.contains(child);
                  const bool provenance_available =
                      provenance_fab_box->contains(child);
                  if (provenance_available) {
                    const int child_origin = provenance(
                        child[0],
#if (AMREX_SPACEDIM >= 2)
                        child[1],
#else
                        0,
#endif
#if (AMREX_SPACEDIM == 3)
                        child[2],
#else
                        0,
#endif
                        0);
                    report << " provenance="
                           << (child_origin == 1
                                   ? "new-from-coarse"
                                   : (child_origin == 0 ? "old-fine-overlap"
                                                        : "unknown"));
                  } else {
                    report << " provenance=unavailable";
                  }
                  if (!state_available) {
                    report << " state=unavailable\n";
                    continue;
                  }
                  append_state(state_a, child[0],
#if (AMREX_SPACEDIM >= 2)
                               child[1],
#else
                               0,
#endif
#if (AMREX_SPACEDIM == 3)
                               child[2]
#else
                               0
#endif
                  );
                  report << '\n';
                }
              }
            }
          } else {
            report << "  regrid_parent_and_children=data-unavailable\n";
          }

          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int side : {-1, 1}) {
              IntVect neighbor = center;
              neighbor[dir] += side;
              const bool in_domain = domain.contains(neighbor);
              const bool same_level_valid = level_boxes.contains(neighbor);
              const bool data_available =
                  state_fab_box.contains(neighbor) &&
                  marker_fab_box.contains(neighbor);
              report << "  neighbor[d" << dir
                     << (side < 0 ? "-" : "+") << "] index=" << neighbor
                     << " topology="
                     << (!in_domain
                             ? "physical-boundary-ghost"
                             : (same_level_valid ? "same-level-valid"
                                                 : "coarse-fine-FillPatch"))
                     << " covered_by_finer="
                     << (covered_by_finer(neighbor) ? "YES" : "no");
              if (!data_available) {
                report << " data=unavailable\n";
                continue;
              }
              const int ni = neighbor[0];
#if (AMREX_SPACEDIM >= 2)
              const int nj = neighbor[1];
#else
              const int nj = 0;
#endif
#if (AMREX_SPACEDIM == 3)
              const int nk = neighbor[2];
#else
              const int nk = 0;
#endif
              const int neighbor_marker0 =
                  static_cast<int>(marker(ni, nj, nk, 0));
              const int neighbor_marker1 =
                  static_cast<int>(marker(ni, nj, nk, 1));
              report << " marker=(" << neighbor_marker0 << ','
                     << neighbor_marker1 << ") class="
                     << marker_class(neighbor_marker0, neighbor_marker1);
              append_state(state_a, ni, nj, nk);
              report << '\n';
            }
          }
          amrex::AllPrint() << report.str() << std::flush;
          ++local_printed;
        }
      }
    }
  }
  ParallelDescriptor::Barrier();
}

// Default-off failure diagnostic for the face-local WENO/LLF decision.  It
// reports the exact six masks and high/low conservative fluxes surrounding the
// worst high-order cell.  It never changes a state or a face flux.
template <typename Closure>
inline void printWorstCellFluxComparison(
    const amrex::MultiFab& high_candidate,
    const amrex::MultiFab& low_candidate,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& high_face_flux,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& low_face_flux,
    const std::array<amrex::iMultiFab*, AMREX_SPACEDIM>& fallback_mask,
    const amrex::Geometry& geom, const int level,
    const Closure& closure, const char* scheme, const int stage_index,
    const int stage_count, const amrex::Real stage_time,
    const amrex::Real stage_dt) {
  using namespace amrex;
  constexpr int face_count = 2 * AMREX_SPACEDIM;
  constexpr int state_value_count = 2 * Closure::NCONS;
  constexpr int face_family_value_count = face_count * Closure::NCONS;
  constexpr int real_value_count =
      state_value_count + 2 * face_family_value_count;
  constexpr int int_value_count = 1 + face_count;

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (high_face_flux[dir] == nullptr || low_face_flux[dir] == nullptr ||
        fallback_mask[dir] == nullptr) {
      amrex::Abort(
          "worst-cell LLF diagnostic requires high/low fluxes and masks");
    }
  }

  const auto& marker_mf = *IBM::ib.bmf_a[level];
  MultiFab high_rhoe(
      high_candidate.boxArray(), high_candidate.DistributionMap(), 1, 0,
      MFInfo().SetArena(The_Async_Arena()));
  for (MFIter mfi(high_candidate, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    const auto high = high_candidate.const_array(mfi);
    const auto marker = marker_mf.const_array(mfi);
    const auto rhoe = high_rhoe.array(mfi);
    ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (marker(i, j, k, 0) != 0) {
            rhoe(i, j, k, 0) = std::numeric_limits<Real>::max();
            return;
          }
          const Real density = high(i, j, k, Closure::URHO);
          const Real mx = high(i, j, k, Closure::UMX);
          const Real my = high(i, j, k, Closure::UMY);
#if (AMREX_SPACEDIM == 3)
          const Real mz = high(i, j, k, Closure::UMZ);
#else
          const Real mz = Real(0.0);
#endif
          const Real energy = high(i, j, k, Closure::UET);
          rhoe(i, j, k, 0) =
              density > Real(0.0)
                  ? energy - Real(0.5) * (mx * mx + my * my + mz * mz) /
                                 density
                  : -std::numeric_limits<Real>::max();
        });
  }
  Gpu::streamSynchronize();
  Real worst_rhoe = high_rhoe.min(0, 0);
  ParallelDescriptor::ReduceRealMin(worst_rhoe);
  const IntVect worst_cell = high_rhoe.minIndex(0);

  Gpu::DeviceVector<Real> device_values(real_value_count);
  Gpu::fillAsync(
      device_values.begin(), device_values.end(),
      [] AMREX_GPU_HOST_DEVICE(Real& value, Long) noexcept {
        value = std::numeric_limits<Real>::quiet_NaN();
      });
  Gpu::DeviceVector<int> device_ints(int_value_count);
  Gpu::fillAsync(
      device_ints.begin(), device_ints.end(),
      [] AMREX_GPU_HOST_DEVICE(int& value, Long) noexcept { value = -1; });
  Real* values = device_values.data();
  int* ints = device_ints.data();

  for (MFIter mfi(high_candidate, false); mfi.isValid(); ++mfi) {
    if (!mfi.validbox().contains(worst_cell)) continue;
    const auto high = high_candidate.const_array(mfi);
    const auto low = low_candidate.const_array(mfi);
    GpuArray<Array4<const Real>, AMREX_SPACEDIM> high_flux_views;
    GpuArray<Array4<const Real>, AMREX_SPACEDIM> low_flux_views;
    GpuArray<Array4<const int>, AMREX_SPACEDIM> mask_views;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      high_flux_views[dir] = high_face_flux[dir]->const_array(mfi);
      low_flux_views[dir] = low_face_flux[dir]->const_array(mfi);
      mask_views[dir] = fallback_mask[dir]->const_array(mfi);
    }
    const Box point_box(worst_cell, worst_cell);
    ParallelFor(
        point_box,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          ints[0] = 1;
          for (int n = 0; n < Closure::NCONS; ++n) {
            values[n] = high(i, j, k, n);
            values[Closure::NCONS + n] = low(i, j, k, n);
          }
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int side = 0; side < 2; ++side) {
              const int face_ordinal = 2 * dir + side;
              const int fi = i + ((dir == 0 && side == 1) ? 1 : 0);
#if (AMREX_SPACEDIM >= 2)
              const int fj = j + ((dir == 1 && side == 1) ? 1 : 0);
#else
              const int fj = j;
#endif
#if (AMREX_SPACEDIM == 3)
              const int fk = k + ((dir == 2 && side == 1) ? 1 : 0);
#else
              const int fk = k;
#endif
              ints[1 + face_ordinal] =
                  mask_views[dir](fi, fj, fk, 0);
              const int high_base =
                  state_value_count + face_ordinal * Closure::NCONS;
              const int low_base = state_value_count +
                  face_family_value_count +
                  face_ordinal * Closure::NCONS;
              for (int n = 0; n < Closure::NCONS; ++n) {
                values[high_base + n] =
                    high_flux_views[dir](fi, fj, fk, n);
                values[low_base + n] =
                    low_flux_views[dir](fi, fj, fk, n);
              }
            }
          }
        });
  }

  Gpu::HostVector<Real> host_values(real_value_count);
  Gpu::HostVector<int> host_ints(int_value_count);
  Gpu::copy(
      Gpu::deviceToHost, device_values.begin(), device_values.end(),
      host_values.begin());
  Gpu::copy(
      Gpu::deviceToHost, device_ints.begin(), device_ints.end(),
      host_ints.begin());
  Gpu::streamSynchronize();
  if (host_ints[0] != 1) return;

  const auto prob_lo = geom.ProbLoArray();
  const auto dx = geom.CellSizeArray();
  amrex::AllPrint()
      << std::setprecision(17)
      << "[LLF-FALLBACK-WORST-CELL] scheme=" << scheme
      << " stage=" << stage_index << '/' << stage_count
      << " level=" << level << " t=" << stage_time
      << " dt_FE=" << stage_dt << " cell=" << worst_cell
      << " xyz=(" << prob_lo[0] + (Real(worst_cell[0]) + Real(0.5)) * dx[0]
#if (AMREX_SPACEDIM >= 2)
      << ',' << prob_lo[1] + (Real(worst_cell[1]) + Real(0.5)) * dx[1]
#endif
#if (AMREX_SPACEDIM == 3)
      << ',' << prob_lo[2] + (Real(worst_cell[2]) + Real(0.5)) * dx[2]
#endif
      << ") min_high_rhoe=" << worst_rhoe << '\n';

  amrex::AllPrint() << "  U_high=";
  for (int n = 0; n < Closure::NCONS; ++n) {
    amrex::AllPrint() << (n == 0 ? "(" : ",") << host_values[n];
  }
  amrex::AllPrint() << ") U_low=";
  for (int n = 0; n < Closure::NCONS; ++n) {
    amrex::AllPrint()
        << (n == 0 ? "(" : ",") << host_values[Closure::NCONS + n];
  }
  amrex::AllPrint() << ")\n";

  GpuArray<Real, Closure::NCONS> reconstructed_correction{};
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    for (int side = 0; side < 2; ++side) {
      const int face_ordinal = 2 * dir + side;
      const int high_base =
          state_value_count + face_ordinal * Closure::NCONS;
      const int low_base = state_value_count + face_family_value_count +
                           face_ordinal * Closure::NCONS;
      amrex::AllPrint()
          << "  face dir=" << dir << " side=" << (side == 0 ? "lo" : "hi")
          << " first_order_llf_mask=" << host_ints[1 + face_ordinal]
          << " dF_high_minus_low=(";
      for (int n = 0; n < Closure::NCONS; ++n) {
        const Real delta_flux =
            host_values[high_base + n] - host_values[low_base + n];
        amrex::AllPrint() << (n == 0 ? "" : ",") << delta_flux;
        if (!geom.IsRZ()) {
          reconstructed_correction[n] +=
              (side == 0 ? Real(1.0) : Real(-1.0)) *
              stage_dt * delta_flux / dx[dir];
        }
      }
      amrex::AllPrint() << ")\n";
    }
  }
  if (!geom.IsRZ()) {
    amrex::AllPrint() << "  Cartesian dt*div(F_high-F_low)=(";
    for (int n = 0; n < Closure::NCONS; ++n) {
      amrex::AllPrint()
          << (n == 0 ? "" : ",") << reconstructed_correction[n];
    }
    amrex::AllPrint() << ") observed U_high-U_low=(";
    for (int n = 0; n < Closure::NCONS; ++n) {
      amrex::AllPrint()
          << (n == 0 ? "" : ",")
          << host_values[n] - host_values[Closure::NCONS + n];
    }
    amrex::AllPrint() << ")\n";
  }
  amrex::AllPrint() << "[LLF-FALLBACK-WORST-CELL-END]\n";
  amrex::ignore_unused(closure);
}

// Construct the deliberately low-order positivity baseline from the same
// stage state used by the production RHS. The fluid and IBM ghost primitives
// are piecewise constant. Cartesian faces receive the usual local
// Lax-Friedrichs/Rusanov flux. R-Z radial faces use the matching metric
// advective average, unweighted physical-state LLF jump, and a separate
// centred pressure companion. Nodal overlap copies are synchronized before
// their divergence is taken, so both cells use one shared numerical face.
template <typename Rhs, typename Closure>
inline amrex::Real computePiecewiseConstantRusanovRhs(
    amrex::MultiFab& stage_state, amrex::MultiFab& low_rhs,
    const amrex::Geometry& geom, int level, const Closure* closure_device,
    const Closure& closure_host, amrex::Real stage_dt,
    std::array<amrex::MultiFab*, AMREX_SPACEDIM>
        captured_face_flux = {},
    amrex::MultiFab* captured_rz_pressure_face_flux = nullptr) {
  using namespace amrex;

  if constexpr (!std::is_base_of_v<no_diffusive_t, Rhs> ||
                !std::is_base_of_v<no_source_t, Rhs>) {
    amrex::Abort(
        "IBM positivity limiter currently supports source-free Euler only");
    return Real(0.0);
  } else {
    const bool is_rz = geom.IsRZ();
    if (is_rz) {
#if (AMREX_SPACEDIM == 2)
      if constexpr (!rhs_detail::rz_paired_pressure_flux_capable<Rhs>::value) {
        amrex::Abort(
            "R-Z IBM positivity limiting requires a paired-pressure Euler "
            "operator");
      }
      const BoxArray expected_pressure_boxes = amrex::convert(
          stage_state.boxArray(), IntVect::TheDimensionVector(0));
      if (captured_rz_pressure_face_flux == nullptr ||
          captured_rz_pressure_face_flux->boxArray() !=
              expected_pressure_boxes ||
          captured_rz_pressure_face_flux->DistributionMap() !=
              stage_state.DistributionMap() ||
          captured_rz_pressure_face_flux->nComp() != 1) {
        amrex::Abort(
            "R-Z low-order positivity baseline requires a one-component "
            "radial pressure-face capture");
      }
#else
      amrex::Abort("R-Z IBM positivity limiting requires a 2-D build");
#endif
    } else if (captured_rz_pressure_face_flux != nullptr) {
      amrex::Abort(
          "Cartesian IBM positivity limiting cannot capture an R-Z pressure "
          "companion");
    }
    if (stage_state.nComp() != Closure::NCONS ||
        low_rhs.nComp() != Closure::NCONS) {
      amrex::Abort("IBM positivity limiter conservative-component mismatch");
    }
    if (!(stage_dt > Real(0.0)) || !std::isfinite(stage_dt)) {
      amrex::Abort("IBM positivity limiter requires a positive stage dt");
    }

#if (AMREX_SPACEDIM < 3)
    stage_state.setVal(Real(0.0), Closure::UMZ, 1, stage_state.nGrow());
#endif

    MultiFab* flux_state = &stage_state;

    const int primitive_nghost = IBM::ib.volumeInterpolationNghost(level);
    MultiFab primitive(flux_state->boxArray(), flux_state->DistributionMap(),
                       Closure::NPRIM, primitive_nghost,
                       MFInfo().SetArena(The_Async_Arena()));
    for (MFIter mfi(*flux_state, false); mfi.isValid(); ++mfi) {
      closure_host.cons2prims(
          mfi, flux_state->array(mfi), primitive.array(mfi));
    }
    if (IBM::ib.nsGPCellAverageRecoveryEnabled()) {
      IBM::ib.computeAllGPsCartesianConservativeAverage(
          primitive, *flux_state, closure_device, level);
    } else if (IBM::ib.rzAnnularCellAverageEnabled(level)) {
      IBM::ib.computeAllGPsRZAnnular(
          primitive, *flux_state, closure_device, level);
    } else {
      IBM::ib.computeAllGPs(primitive, closure_device, level);
    }

    std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> face_flux;
    std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> face_speed;
    std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> face_area;
    std::unique_ptr<MultiFab> radial_pressure_face_flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      BoxArray face_boxes = amrex::convert(
          stage_state.boxArray(), IntVect::TheDimensionVector(dir));
      face_flux[dir] = std::make_unique<MultiFab>(
          face_boxes, stage_state.DistributionMap(), Closure::NCONS, 0,
          MFInfo().SetArena(The_Async_Arena()));
      face_speed[dir] = std::make_unique<MultiFab>(
          face_boxes, stage_state.DistributionMap(), 1, 0,
          MFInfo().SetArena(The_Async_Arena()));
      face_area[dir] = std::make_unique<MultiFab>(
          face_boxes, stage_state.DistributionMap(), 1, 0,
          MFInfo().SetArena(The_Async_Arena()));
      geom.GetFaceArea(*face_area[dir], dir);
      face_flux[dir]->setVal(Real(0.0));
      face_speed[dir]->setVal(Real(0.0));
    }
    if (is_rz) {
      const BoxArray radial_face_boxes = amrex::convert(
          stage_state.boxArray(), IntVect::TheDimensionVector(0));
      radial_pressure_face_flux = std::make_unique<MultiFab>(
          radial_face_boxes, stage_state.DistributionMap(), 1, 0,
          MFInfo().SetArena(The_Async_Arena()));
      radial_pressure_face_flux->setVal(Real(0.0));
    }
    MultiFab cell_volume(stage_state.boxArray(),
                         stage_state.DistributionMap(), 1, 0,
                         MFInfo().SetArena(The_Async_Arena()));
    geom.GetVolume(cell_volume);

    const auto& marker_mf = *IBM::ib.bmf_a[level];
    const auto prob_lo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const int radial_axis_face_index = geom.Domain().smallEnd(0);
    for (MFIter mfi(stage_state, false); mfi.isValid(); ++mfi) {
      const auto prim = primitive.const_array(mfi);
      const auto marker = marker_mf.const_array(mfi);
      Array4<Real> pressure_flux;
      if (is_rz) {
        pressure_flux = radial_pressure_face_flux->array(mfi);
      }
      const Box& cell_box = mfi.tilebox();
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const Box face_box = amrex::surroundingNodes(cell_box, dir);
        const auto flux = face_flux[dir]->array(mfi);
        const auto speed = face_speed[dir]->array(mfi);
        ParallelFor(
            face_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const IntVect right(AMREX_D_DECL(i, j, k));
              const IntVect left =
                  right - IntVect::TheDimensionVector(dir);
              IntVect left_marker = left;
              IntVect right_marker = right;
#if (AMREX_SPACEDIM == 2)
              const bool radial_axis_face =
                  is_rz && dir == 0 && prob_lo[0] == Real(0.0) &&
                  right[0] == radial_axis_face_index;
              if (radial_axis_face &&
                  left_marker[0] < radial_axis_face_index) {
                left_marker[0] =
                    2 * radial_axis_face_index - 1 - left_marker[0];
              }
#else
              const bool radial_axis_face = false;
#endif
              const bool left_solid = marker(left_marker, 0) != 0;
              const bool right_solid = marker(right_marker, 0) != 0;
              if (left_solid && right_solid) {
                speed(i, j, k) = Real(0.0);
                for (int n = 0; n < Closure::NCONS; ++n) {
                  flux(i, j, k, n) = Real(0.0);
                }
                if (is_rz && dir == 0) {
                  pressure_flux(i, j, k, 0) = Real(0.0);
                }
                return;
              }

              const Real rho_left = prim(left, Closure::QRHO);
              const Real rho_right = prim(right, Closure::QRHO);
              const Real pressure_left = prim(left, Closure::QPRES);
              const Real pressure_right = prim(right, Closure::QPRES);
              const Real sound_left = prim(left, Closure::QC);
              const Real sound_right = prim(right, Closure::QC);
              const Real velocity_left = prim(left, Closure::QU + dir);
              const Real velocity_right = prim(right, Closure::QU + dir);
              const bool side_states_admissible =
                  rho_left > CNSConstants::smallr &&
                  rho_right > CNSConstants::smallr &&
                  pressure_left > CNSConstants::min_press() &&
                  pressure_right > CNSConstants::min_press() &&
                  sound_left >= Real(0.0) && sound_right >= Real(0.0) &&
                  std::isfinite(rho_left + rho_right + pressure_left +
                                pressure_right + sound_left + sound_right +
                                velocity_left + velocity_right);
              if (!side_states_admissible) {
                const Real nan =
                    std::numeric_limits<Real>::quiet_NaN();
                speed(i, j, k) = nan;
                for (int n = 0; n < Closure::NCONS; ++n) {
                  flux(i, j, k, n) = nan;
                }
                if (is_rz && dir == 0) {
                  pressure_flux(i, j, k, 0) = nan;
                }
                return;
              }

              Real conservative_left[Closure::NCONS];
              Real conservative_right[Closure::NCONS];
              Real physical_flux_left[Closure::NCONS];
              Real physical_flux_right[Closure::NCONS];
              closure_device->prims2cons(
                  left, prim, conservative_left);
              closure_device->prims2cons(
                  right, prim, conservative_right);
              closure_device->prims2flux(
                  left, dir, prim, physical_flux_left);
              closure_device->prims2flux(
                  right, dir, prim, physical_flux_right);
              const Real alpha = amrex::max(
                  amrex::Math::abs(velocity_left) + sound_left,
                  amrex::Math::abs(velocity_right) + sound_right);
              speed(i, j, k) = alpha;
#if (AMREX_SPACEDIM == 2)
              if (is_rz && dir == 0) {
                // The low-order member uses the same paired R-Z operator as
                // the high-order WENO member.  Away from the axis its central
                // advective flux is metric weighted, while the LLF
                // dissipation acts on the unweighted physical state jump so a
                // uniform field has exactly zero dissipation.
                if (radial_axis_face) {
                  for (int n = 0; n < Closure::NCONS; ++n) {
                    // At the axis this slot is a metric h-flux rather than a
                    // per-area flux.  The monotone low-order baseline uses the
                    // regular zero metric flux; the high-order axis auxiliary
                    // enters the antidiffusive correction explicitly.
                    flux(i, j, k, n) = Real(0.0);
                  }
                  pressure_flux(i, j, k, 0) = pressure_right;
                  return;
                }

                const Real radial_face =
                    prob_lo[0] +
                    Real(right[0] - radial_axis_face_index) * dx[0];
                const Real radial_left = radial_face - Real(0.5) * dx[0];
                const Real radial_right = radial_face + Real(0.5) * dx[0];
                if (!(radial_face > Real(0.0))) {
                  const Real nan =
                      std::numeric_limits<Real>::quiet_NaN();
                  pressure_flux(i, j, k, 0) = nan;
                  for (int n = 0; n < Closure::NCONS; ++n) {
                    flux(i, j, k, n) = nan;
                  }
                  return;
                }
                const Real pressure_face =
                    Real(0.5) * (pressure_left + pressure_right);
                physical_flux_left[Closure::UMX] -= pressure_left;
                physical_flux_right[Closure::UMX] -= pressure_right;
                for (int n = 0; n < Closure::NCONS; ++n) {
                  flux(i, j, k, n) =
                      Real(0.5) *
                          (radial_left * physical_flux_left[n] +
                           radial_right * physical_flux_right[n]) /
                          radial_face -
                      Real(0.5) * alpha *
                          (conservative_right[n] - conservative_left[n]);
                }
                flux(i, j, k, Closure::UMX) += pressure_face;
                pressure_flux(i, j, k, 0) = pressure_face;
                return;
              }
#endif
              for (int n = 0; n < Closure::NCONS; ++n) {
                flux(i, j, k, n) = Real(0.5) *
                    (physical_flux_left[n] + physical_flux_right[n] -
                     alpha *
                         (conservative_right[n] - conservative_left[n]));
              }
            });
      }
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      face_flux[dir]->OverrideSync(geom.periodicity());
      face_speed[dir]->OverrideSync(geom.periodicity());
      if (captured_face_flux[dir] != nullptr) {
        if (captured_face_flux[dir]->boxArray() !=
                face_flux[dir]->boxArray() ||
            captured_face_flux[dir]->DistributionMap() !=
                face_flux[dir]->DistributionMap() ||
            captured_face_flux[dir]->nComp() != Closure::NCONS) {
          amrex::Abort(
              "low-order face-flux capture layout mismatch");
        }
        MultiFab::Copy(*captured_face_flux[dir], *face_flux[dir], 0, 0,
                       Closure::NCONS, 0);
      }
    }
    if (is_rz) {
      radial_pressure_face_flux->OverrideSync(geom.periodicity());
      MultiFab::Copy(*captured_rz_pressure_face_flux,
                     *radial_pressure_face_flux, 0, 0, 1, 0);
      captured_rz_pressure_face_flux->OverrideSync(geom.periodicity());
    }

    low_rhs.setVal(Real(0.0));
    ReduceOps<ReduceOpMax> cfl_reduce;
    ReduceData<Real> cfl_data(cfl_reduce);
    using CflTuple = typename decltype(cfl_data)::Type;
    for (MFIter mfi(low_rhs, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const auto rhs = low_rhs.array(mfi);
      const auto marker = marker_mf.const_array(mfi);
      const auto volume = cell_volume.const_array(mfi);
      const auto fx = face_flux[0]->const_array(mfi);
      const auto ax = face_speed[0]->const_array(mfi);
      const auto area_x = face_area[0]->const_array(mfi);
      Array4<const Real> pressure_x;
      if (is_rz) {
        pressure_x = radial_pressure_face_flux->const_array(mfi);
      }
#if (AMREX_SPACEDIM >= 2)
      const auto fy = face_flux[1]->const_array(mfi);
      const auto ay = face_speed[1]->const_array(mfi);
      const auto area_y = face_area[1]->const_array(mfi);
#endif
#if (AMREX_SPACEDIM == 3)
      const auto fz = face_flux[2]->const_array(mfi);
      const auto az = face_speed[2]->const_array(mfi);
      const auto area_z = face_area[2]->const_array(mfi);
#endif
      ParallelFor(
          bx, Closure::NCONS,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (marker(i, j, k, 0) != 0) {
              rhs(i, j, k, n) = Real(0.0);
              return;
            }
            const Real inverse_volume = Real(1.0) / volume(i, j, k);
            Real value;
#if (AMREX_SPACEDIM == 2)
            if (is_rz && n == Closure::UMX) {
              const Real pressure_lo = pressure_x(i, j, k, 0);
              const Real pressure_hi = pressure_x(i + 1, j, k, 0);
              value = inverse_volume *
                      (area_x(i, j, k) *
                           (fx(i, j, k, n) - pressure_lo) -
                       area_x(i + 1, j, k) *
                           (fx(i + 1, j, k, n) - pressure_hi)) +
                  (pressure_lo - pressure_hi) / dx[0];
            } else
#endif
            {
              value = inverse_volume *
                  (area_x(i, j, k) * fx(i, j, k, n) -
                   area_x(i + 1, j, k) * fx(i + 1, j, k, n));
            }
#if (AMREX_SPACEDIM >= 2)
            value += inverse_volume *
                (area_y(i, j, k) * fy(i, j, k, n) -
                 area_y(i, j + 1, k) * fy(i, j + 1, k, n));
#endif
#if (AMREX_SPACEDIM == 3)
            value += inverse_volume *
                (area_z(i, j, k) * fz(i, j, k, n) -
                 area_z(i, j, k + 1) * fz(i, j, k + 1, n));
#endif
            rhs(i, j, k, n) = value;
          });
      cfl_reduce.eval(
          bx, cfl_data,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) -> CflTuple {
            if (marker(i, j, k, 0) != 0) return {Real(0.0)};
            const Real inverse_volume = Real(1.0) / volume(i, j, k);
            Real cfl = Real(0.5) * stage_dt * inverse_volume *
                (area_x(i, j, k) * ax(i, j, k) +
                 area_x(i + 1, j, k) * ax(i + 1, j, k));
#if (AMREX_SPACEDIM >= 2)
            cfl += Real(0.5) * stage_dt * inverse_volume *
                (area_y(i, j, k) * ay(i, j, k) +
                 area_y(i, j + 1, k) * ay(i, j + 1, k));
#endif
#if (AMREX_SPACEDIM == 3)
            cfl += Real(0.5) * stage_dt * inverse_volume *
                (area_z(i, j, k) * az(i, j, k) +
                 area_z(i, j, k + 1) * az(i, j, k + 1));
#endif
            return {cfl};
          });
    }
    Real maximum_cfl = amrex::get<0>(cfl_data.value(cfl_reduce));
    ParallelDescriptor::ReduceRealMax(maximum_cfl);
    return maximum_cfl;
  }
}

template <typename Closure>
inline amrex::Real maximumGlobalAdmissibleTheta(
    const amrex::MultiFab& low_state, const amrex::MultiFab& high_state,
    int level, const Closure& closure) {
  using namespace amrex;
  const auto thresholds = stageAdmissibilityThresholds(closure);
  const auto& marker_mf = *IBM::ib.bmf_a[level];
  ReduceOps<ReduceOpMin> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (MFIter mfi(low_state, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const auto low = low_state.const_array(mfi);
    const auto high = high_state.const_array(mfi);
    const auto marker = marker_mf.const_array(mfi);
    reduce_op.eval(
        bx, reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
          if (marker(i, j, k, 0) != 0) return {Real(1.0)};

          auto admissible_at = [=] AMREX_GPU_DEVICE(Real theta) noexcept {
            Real state[Closure::NCONS];
            bool finite = true;
            for (int n = 0; n < Closure::NCONS; ++n) {
              state[n] = low(i, j, k, n);
              if (theta > Real(0.0)) {
                state[n] +=
                    theta * (high(i, j, k, n) - low(i, j, k, n));
              }
              finite = finite && std::isfinite(state[n]);
            }
            const Real rho = state[Closure::URHO];
            if (!finite || !(rho > thresholds.rho)) return false;
            const Real mx = state[Closure::UMX];
            const Real my = state[Closure::UMY];
#if (AMREX_SPACEDIM == 3)
            const Real mz = state[Closure::UMZ];
#else
            const Real mz = Real(0.0);
#endif
            const Real rhoe = state[Closure::UET] -
                Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
            return std::isfinite(rhoe) &&
                   rhoe > thresholds.rhoe_min(rho);
          };

          if (!admissible_at(Real(0.0))) return {Real(0.0)};
          if (admissible_at(Real(1.0))) return {Real(1.0)};
          Real lower = Real(0.0);
          Real upper = Real(1.0);
          for (int iteration = 0; iteration < 64; ++iteration) {
            const Real midpoint = Real(0.5) * (lower + upper);
            if (admissible_at(midpoint)) {
              lower = midpoint;
            } else {
              upper = midpoint;
            }
          }
          return {lower};
        });
  }

  Real theta = amrex::get<0>(reduce_data.value(reduce_op));
  ParallelDescriptor::ReduceRealMin(theta);
  return theta;
}

struct LocalSharedFaceLimiterStats {
  bool high_order_admissible = false;
  bool low_order_admissible = false;
  bool limited_state_admissible = false;
  amrex::Long active_faces = 0;
  amrex::Long limited_faces = 0;
  amrex::Long limited_fluid_faces = 0;
  amrex::Long limited_crossing_faces = 0;
  amrex::Real minimum_theta = amrex::Real(1.0);
  amrex::Real conservation_residual_linf = amrex::Real(0.0);
  amrex::Real conservation_relative_linf = amrex::Real(0.0);
};

// Sufficient local convex limiter for a source-free Cartesian/RZ Euler stage.
// Write the antidiffusive correction of cell i as
//
//   U_i = U_i^L + sum_f theta_f C_if
//       = (1/N_f) sum_f (U_i^L + N_f theta_f C_if),
//
// where N_f=2*dim. In paired R-Z radial momentum, C_if contains both the
// metric complete-flux correction and the ordinary-gradient pressure
// correction.  The same shared theta is applied to F and P, including the
// zero-area axis face. A per-cell/per-face bound makes every state in the
// final average admissible. A shared face then receives the
// minimum of the two bounds supplied by its adjacent active cells. This differs from applying
// min(theta_i,theta_j) to cellwise high/low line segments: the convex
// decomposition above remains valid when different faces receive different
// theta values, while the naive cellwise construction does not.
template <typename Closure>
inline LocalSharedFaceLimiterStats applyLocalSharedFaceConvexLimiter(
    amrex::MultiFab& high_candidate, const amrex::MultiFab& low_candidate,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& high_face_flux,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& low_face_flux,
    amrex::MultiFab* high_rz_pressure_face_flux,
    const amrex::MultiFab* low_rz_pressure_face_flux,
    const amrex::Geometry& geom, int level, const Closure& closure,
    const char* scheme, int stage_index, int stage_count,
    amrex::Real stage_time, amrex::Real stage_dt,
    amrex::Real multidimensional_cfl, bool abort_on_failure = true,
    bool replace_high_order_flux = false) {
  using namespace amrex;
  constexpr int face_count = 2 * AMREX_SPACEDIM;
  LocalSharedFaceLimiterStats result;
  const auto thresholds = stageAdmissibilityThresholds(closure);
  const auto high_stats =
      collectStageAdmissibility(high_candidate, level, closure);
  const auto low_stats =
      collectStageAdmissibility(low_candidate, level, closure);
  result.high_order_admissible = high_stats.admissible();
  result.low_order_admissible = low_stats.admissible();
  if (result.high_order_admissible) {
    result.limited_state_admissible = true;
    return result;
  }

  amrex::Print()
      << "\n============= IBM LOCAL SHARED-FACE CONVEX LIMITER =============\n"
      << "  scheme=" << scheme << " stage=" << stage_index << '/'
      << stage_count << " t=" << std::setprecision(17) << stage_time
      << " stage_dt=" << stage_dt
      << " multidimensional_cfl=" << multidimensional_cfl << '\n';
  printStageAdmissibility("high-order", high_stats);
  printStageAdmissibility("all-low-order", low_stats);
  if (!result.low_order_admissible) {
    amrex::Print()
        << "  verdict=ALL_LOW_ORDER_INADMISSIBLE; local flux limiting is not "
           "permitted.\n"
        << "=================================================================\n";
    if (abort_on_failure) {
      amrex::Abort(
          "local shared-face limiter: all-low-order update is inadmissible");
    }
    return result;
  }

  if (geom.isAnyPeriodic()) {
    amrex::Abort(
        "local shared-face limiter currently requires non-periodic "
        "boundaries");
  }
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (high_face_flux[dir] == nullptr || low_face_flux[dir] == nullptr) {
      amrex::Abort("local shared-face limiter requires both flux families");
    }
  }
  const bool is_rz = geom.IsRZ();
  if (is_rz) {
#if (AMREX_SPACEDIM == 2)
    const BoxArray expected_pressure_boxes = amrex::convert(
        low_candidate.boxArray(), IntVect::TheDimensionVector(0));
    if (high_rz_pressure_face_flux == nullptr ||
        low_rz_pressure_face_flux == nullptr ||
        high_rz_pressure_face_flux->boxArray() != expected_pressure_boxes ||
        low_rz_pressure_face_flux->boxArray() != expected_pressure_boxes ||
        high_rz_pressure_face_flux->DistributionMap() !=
            low_candidate.DistributionMap() ||
        low_rz_pressure_face_flux->DistributionMap() !=
            low_candidate.DistributionMap() ||
        high_rz_pressure_face_flux->nComp() != 1 ||
        low_rz_pressure_face_flux->nComp() != 1) {
      amrex::Abort(
          "paired R-Z local limiter pressure-face layout mismatch");
    }
#else
    amrex::Abort("paired R-Z local limiter requires a 2-D build");
#endif
  } else if (high_rz_pressure_face_flux != nullptr ||
             low_rz_pressure_face_flux != nullptr) {
    amrex::Abort(
        "Cartesian local limiter cannot receive R-Z pressure companions");
  }

  const auto& marker_mf = *IBM::ib.bmf_a[level];
  const auto cell_size = geom.CellSizeArray();
  const auto prob_lo = geom.ProbLoArray();
  const int radial_axis_face_index = geom.Domain().smallEnd(0);
  MultiFab cell_face_bound(
      low_candidate.boxArray(), low_candidate.DistributionMap(), face_count, 1,
      MFInfo().SetArena(The_Async_Arena()));
  cell_face_bound.setVal(Real(1.0));
  // A level BoxArray need not cover the full problem domain.  In particular,
  // a fine-grid face at a coarse-fine interface has only one same-level valid
  // neighbour.  Geometry markers exist in ghost cells too, so they cannot by
  // themselves distinguish that interface from a true active/active face.
  // Keep an explicit same-level mask for shared-theta ownership and closure
  // diagnostics.  This is a full-background-cell mask; it contains no cut
  // volume, open-area fraction, or aggregate topology.
  iMultiFab valid_on_level(
      low_candidate.boxArray(), low_candidate.DistributionMap(), 1, 1,
      MFInfo().SetArena(The_Async_Arena()));
  valid_on_level.setVal(0);
  for (MFIter mfi(valid_on_level, false); mfi.isValid(); ++mfi) {
    const Box& valid_box = mfi.validbox();
    const auto valid = valid_on_level.array(mfi);
    ParallelFor(
        valid_box,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          valid(i, j, k) = 1;
        });
  }
  valid_on_level.FillBoundary(geom.periodicity());
  MultiFab cell_volume(low_candidate.boxArray(),
                       low_candidate.DistributionMap(), 1, 0,
                       MFInfo().SetArena(The_Async_Arena()));
  geom.GetVolume(cell_volume);
  std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> metric_face_area;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    BoxArray face_boxes = amrex::convert(
        low_candidate.boxArray(), IntVect::TheDimensionVector(dir));
    metric_face_area[dir] = std::make_unique<MultiFab>(
        face_boxes, low_candidate.DistributionMap(), 1, 0,
        MFInfo().SetArena(The_Async_Arena()));
    geom.GetFaceArea(*metric_face_area[dir], dir);
  }

  for (MFIter mfi(low_candidate, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const auto low_state = low_candidate.const_array(mfi);
    const auto marker = marker_mf.const_array(mfi);
    const auto bound = cell_face_bound.array(mfi);
    const auto volume = cell_volume.const_array(mfi);
    const auto high_fx = high_face_flux[0]->const_array(mfi);
    const auto low_fx = low_face_flux[0]->const_array(mfi);
    const auto area_x = metric_face_area[0]->const_array(mfi);
    Array4<const Real> high_pressure_x;
    Array4<const Real> low_pressure_x;
    if (is_rz) {
      high_pressure_x = high_rz_pressure_face_flux->const_array(mfi);
      low_pressure_x = low_rz_pressure_face_flux->const_array(mfi);
    }
#if (AMREX_SPACEDIM >= 2)
    const auto high_fy = high_face_flux[1]->const_array(mfi);
    const auto low_fy = low_face_flux[1]->const_array(mfi);
    const auto area_y = metric_face_area[1]->const_array(mfi);
#endif
#if (AMREX_SPACEDIM == 3)
    const auto high_fz = high_face_flux[2]->const_array(mfi);
    const auto low_fz = low_face_flux[2]->const_array(mfi);
    const auto area_z = metric_face_area[2]->const_array(mfi);
#endif
    ParallelFor(
        bx, face_count,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int face_slot) noexcept {
          if (marker(i, j, k, 0) != 0) {
            bound(i, j, k, face_slot) = Real(1.0);
            return;
          }
          const int dir = face_slot / 2;
          const int side = face_slot - 2 * dir;
          int fi = i;
          int fj = j;
          int fk = k;
          if (side != 0) {
            if (dir == 0) ++fi;
#if (AMREX_SPACEDIM >= 2)
            if (dir == 1) ++fj;
#endif
#if (AMREX_SPACEDIM == 3)
            if (dir == 2) ++fk;
#endif
          }
          const Real orientation = side == 0 ? Real(1.0) : Real(-1.0);
          Real selected_face_area = area_x(fi, fj, fk);
#if (AMREX_SPACEDIM >= 2)
          if (dir == 1) selected_face_area = area_y(fi, fj, fk);
#endif
#if (AMREX_SPACEDIM == 3)
          if (dir == 2) selected_face_area = area_z(fi, fj, fk);
#endif
          const Real inverse_volume = Real(1.0) / volume(i, j, k);
          const Real scale = orientation * stage_dt * selected_face_area *
                             inverse_volume;

          auto correction = [=] AMREX_GPU_DEVICE(int n) noexcept {
            if (dir == 0) {
              const Real flux_difference =
                  high_fx(fi, fj, fk, n) - low_fx(fi, fj, fk, n);
#if (AMREX_SPACEDIM == 2)
              if (is_rz) {
                const bool axis_face =
                    prob_lo[0] == Real(0.0) &&
                    fi == radial_axis_face_index;
                if (n == Closure::UMX) {
                  const Real pressure_difference =
                      high_pressure_x(fi, fj, fk, 0) -
                      low_pressure_x(fi, fj, fk, 0);
                  // D_pair(F,P)=D_metric(F-P)+D_cart(P).  This form keeps
                  // the ordinary pressure derivative visible at the zero-area
                  // axis face and applies one theta to both channels.
                  return scale * flux_difference +
                      orientation * stage_dt *
                          (Real(1.0) / cell_size[0] -
                           selected_face_area * inverse_volume) *
                          pressure_difference;
                }
                if (axis_face) {
                  // WENO stores the undivided metric h-flux at r=0 for the
                  // non-radial-momentum components.  AMReX's physical axis
                  // area is zero, so its finite-difference correction needs
                  // this explicit coefficient.
                  return orientation * stage_dt * Real(2.0) /
                         (cell_size[0] * cell_size[0]) * flux_difference;
                }
              }
#endif
              return scale * flux_difference;
            }
#if (AMREX_SPACEDIM >= 2)
            if (dir == 1) {
              return scale *
                  (high_fy(fi, fj, fk, n) - low_fy(fi, fj, fk, n));
            }
#endif
#if (AMREX_SPACEDIM == 3)
            return scale *
                (high_fz(fi, fj, fk, n) - low_fz(fi, fj, fk, n));
#else
            return Real(0.0);
#endif
          };

          auto admissible_at = [=] AMREX_GPU_DEVICE(Real theta) noexcept {
            Real state[Closure::NCONS];
            bool finite = true;
            for (int n = 0; n < Closure::NCONS; ++n) {
              state[n] = low_state(i, j, k, n);
              if (theta > Real(0.0)) {
                state[n] += Real(face_count) * theta * correction(n);
              }
              finite = finite && std::isfinite(state[n]);
            }
            const Real rho = state[Closure::URHO];
            if (!finite || !(rho > thresholds.rho)) return false;
            const Real mx = state[Closure::UMX];
            const Real my = state[Closure::UMY];
#if (AMREX_SPACEDIM == 3)
            const Real mz = state[Closure::UMZ];
#else
            const Real mz = Real(0.0);
#endif
            const Real rhoe = state[Closure::UET] -
                Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
            return std::isfinite(rhoe) &&
                   rhoe > thresholds.rhoe_min(rho);
          };

          if (!admissible_at(Real(0.0))) {
            bound(i, j, k, face_slot) = Real(0.0);
            return;
          }
          if (admissible_at(Real(1.0))) {
            bound(i, j, k, face_slot) = Real(1.0);
            return;
          }
          Real lower = Real(0.0);
          Real upper = Real(1.0);
          for (int iteration = 0; iteration < 64; ++iteration) {
            const Real midpoint = Real(0.5) * (lower + upper);
            if (admissible_at(midpoint)) {
              lower = midpoint;
            } else {
              upper = midpoint;
            }
          }
          bound(i, j, k, face_slot) =
              lower * (Real(1.0) - Real(1.0e-12));
        });
  }
  cell_face_bound.FillBoundary(geom.periodicity());

  std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> face_theta;
  const Box domain = geom.Domain();
  const int ilo = domain.smallEnd(0);
  const int ihi = domain.bigEnd(0);
#if (AMREX_SPACEDIM >= 2)
  const int jlo = domain.smallEnd(1);
  const int jhi = domain.bigEnd(1);
#endif
#if (AMREX_SPACEDIM == 3)
  const int klo = domain.smallEnd(2);
  const int khi = domain.bigEnd(2);
#endif
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    BoxArray face_boxes = amrex::convert(
        low_candidate.boxArray(), IntVect::TheDimensionVector(dir));
    face_theta[dir] = std::make_unique<MultiFab>(
        face_boxes, low_candidate.DistributionMap(), 1, 0,
        MFInfo().SetArena(The_Async_Arena()));
    face_theta[dir]->setVal(Real(1.0));
  }

  for (MFIter mfi(low_candidate, false); mfi.isValid(); ++mfi) {
    const auto marker = marker_mf.const_array(mfi);
    const auto valid = valid_on_level.const_array(mfi);
    const auto bound = cell_face_bound.const_array(mfi);
    const Box& cell_box = mfi.tilebox();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const Box face_box = amrex::surroundingNodes(cell_box, dir);
      const auto theta = face_theta[dir]->array(mfi);
      ParallelFor(
          face_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int li = i;
            int lj = j;
            int lk = k;
            if (dir == 0) --li;
#if (AMREX_SPACEDIM >= 2)
            if (dir == 1) --lj;
#endif
#if (AMREX_SPACEDIM == 3)
            if (dir == 2) --lk;
#endif
            const bool left_inside = li >= ilo && li <= ihi
#if (AMREX_SPACEDIM >= 2)
                && lj >= jlo && lj <= jhi
#endif
#if (AMREX_SPACEDIM == 3)
                && lk >= klo && lk <= khi
#endif
                ;
            const bool right_inside = i >= ilo && i <= ihi
#if (AMREX_SPACEDIM >= 2)
                && j >= jlo && j <= jhi
#endif
#if (AMREX_SPACEDIM == 3)
                && k >= klo && k <= khi
#endif
                ;
            Real value = Real(1.0);
            if (left_inside && valid(li, lj, lk) != 0 &&
                marker(li, lj, lk, 0) == 0) {
              value = amrex::min(value, bound(li, lj, lk, 2 * dir + 1));
            }
            if (right_inside && valid(i, j, k) != 0 &&
                marker(i, j, k, 0) == 0) {
              value = amrex::min(value, bound(i, j, k, 2 * dir));
            }
            theta(i, j, k) = value;
          });
    }
  }
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    face_theta[dir]->OverrideSync(geom.periodicity());
  }

  Long invalid_face_theta = 0;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    ReduceOps<ReduceOpSum> theta_reduce_ops;
    ReduceData<Long> theta_reduce_data(theta_reduce_ops);
    using ThetaReduceTuple = typename decltype(theta_reduce_data)::Type;
    for (MFIter mfi(*face_theta[dir], false); mfi.isValid(); ++mfi) {
      const Box& face_box = mfi.validbox();
      const auto theta = face_theta[dir]->const_array(mfi);
      theta_reduce_ops.eval(
          face_box, theta_reduce_data,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ThetaReduceTuple {
            const Real value = theta(i, j, k);
            return {std::isfinite(value) && value >= Real(0.0) &&
                            value <= Real(1.0)
                        ? Long(0)
                        : Long(1)};
          });
    }
    invalid_face_theta += amrex::get<0>(
        theta_reduce_data.value(theta_reduce_ops));
  }
  ParallelDescriptor::ReduceLongSum(invalid_face_theta);
  if (invalid_face_theta != 0) {
    amrex::Abort(
        "local shared-face limiter produced a nonfinite or out-of-range "
        "face theta");
  }

  MultiFab limited_candidate(
      low_candidate.boxArray(), low_candidate.DistributionMap(),
      Closure::NCONS, 0, MFInfo().SetArena(The_Async_Arena()));
  MultiFab::Copy(limited_candidate, low_candidate, 0, 0, Closure::NCONS, 0);
  for (MFIter mfi(limited_candidate, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const auto limited = limited_candidate.array(mfi);
    const auto marker = marker_mf.const_array(mfi);
    const auto volume = cell_volume.const_array(mfi);
    const auto high_fx = high_face_flux[0]->const_array(mfi);
    const auto low_fx = low_face_flux[0]->const_array(mfi);
    const auto theta_x = face_theta[0]->const_array(mfi);
    const auto area_x = metric_face_area[0]->const_array(mfi);
    Array4<const Real> high_pressure_x;
    Array4<const Real> low_pressure_x;
    if (is_rz) {
      high_pressure_x = high_rz_pressure_face_flux->const_array(mfi);
      low_pressure_x = low_rz_pressure_face_flux->const_array(mfi);
    }
#if (AMREX_SPACEDIM >= 2)
    const auto high_fy = high_face_flux[1]->const_array(mfi);
    const auto low_fy = low_face_flux[1]->const_array(mfi);
    const auto theta_y = face_theta[1]->const_array(mfi);
    const auto area_y = metric_face_area[1]->const_array(mfi);
#endif
#if (AMREX_SPACEDIM == 3)
    const auto high_fz = high_face_flux[2]->const_array(mfi);
    const auto low_fz = low_face_flux[2]->const_array(mfi);
    const auto theta_z = face_theta[2]->const_array(mfi);
    const auto area_z = metric_face_area[2]->const_array(mfi);
#endif
    ParallelFor(
        bx, Closure::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          if (marker(i, j, k, 0) != 0) return;
          auto limited_difference = [] AMREX_GPU_DEVICE(
                                        Real theta, Real high,
                                        Real low) noexcept {
            return theta > Real(0.0) ? theta * (high - low) : Real(0.0);
          };
          const Real inverse_volume = Real(1.0) / volume(i, j, k);
          auto radial_face_correction =
              [=] AMREX_GPU_DEVICE(int face_i, Real orientation) noexcept {
                const Real theta = theta_x(face_i, j, k);
                const Real flux_difference = limited_difference(
                    theta, high_fx(face_i, j, k, n),
                    low_fx(face_i, j, k, n));
                const Real area = area_x(face_i, j, k);
#if (AMREX_SPACEDIM == 2)
                if (is_rz) {
                  const bool axis_face =
                      prob_lo[0] == Real(0.0) &&
                      face_i == radial_axis_face_index;
                  if (n == Closure::UMX) {
                    const Real pressure_difference = limited_difference(
                        theta, high_pressure_x(face_i, j, k, 0),
                        low_pressure_x(face_i, j, k, 0));
                    return orientation * stage_dt *
                        (area * inverse_volume * flux_difference +
                         (Real(1.0) / cell_size[0] -
                          area * inverse_volume) * pressure_difference);
                  }
                  if (axis_face) {
                    return orientation * stage_dt * Real(2.0) /
                           (cell_size[0] * cell_size[0]) * flux_difference;
                  }
                }
#endif
                return orientation * stage_dt * area * inverse_volume *
                       flux_difference;
              };
          Real correction =
              radial_face_correction(i, Real(1.0)) +
              radial_face_correction(i + 1, Real(-1.0));
#if (AMREX_SPACEDIM >= 2)
          correction += stage_dt * inverse_volume *
              (area_y(i, j, k) *
                   limited_difference(theta_y(i, j, k),
                                      high_fy(i, j, k, n),
                                      low_fy(i, j, k, n)) -
               area_y(i, j + 1, k) *
                   limited_difference(theta_y(i, j + 1, k),
                                      high_fy(i, j + 1, k, n),
                                      low_fy(i, j + 1, k, n)));
#endif
#if (AMREX_SPACEDIM == 3)
          correction += stage_dt * inverse_volume *
              (area_z(i, j, k) *
                   limited_difference(theta_z(i, j, k),
                                      high_fz(i, j, k, n),
                                      low_fz(i, j, k, n)) -
               area_z(i, j, k + 1) *
                   limited_difference(theta_z(i, j, k + 1),
                                      high_fz(i, j, k + 1, n),
                                      low_fz(i, j, k + 1, n)));
#endif
          limited(i, j, k, n) += correction;
        });
  }

  // Verify conservation over the active metric control volumes. Internal
  // active/active faces must cancel because both cells read the same theta_f.
  // The remaining state change must equal the limited antidiffusive complete
  // flux through outer-domain/active-solid faces plus the explicitly reported
  // paired-pressure (and axis metric-h) auxiliary.  The latter is source-like
  // in the cylindrical balance and must not be mistaken for a conservative
  // physical-face flux.
  for (int n = 0; n < Closure::NCONS; ++n) {
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum>
        cell_reduce_ops;
    ReduceData<Real, Real, Real, Real, Real> cell_reduce_data(cell_reduce_ops);
    using CellReduceTuple = typename decltype(cell_reduce_data)::Type;
    for (MFIter mfi(limited_candidate, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      const auto limited = limited_candidate.const_array(mfi);
      const auto low = low_candidate.const_array(mfi);
      const auto marker = marker_mf.const_array(mfi);
      const auto volume = cell_volume.const_array(mfi);
      const auto theta_x = face_theta[0]->const_array(mfi);
      const auto area_x = metric_face_area[0]->const_array(mfi);
      const auto high_fx = high_face_flux[0]->const_array(mfi);
      const auto low_fx = low_face_flux[0]->const_array(mfi);
      Array4<const Real> high_pressure_x;
      Array4<const Real> low_pressure_x;
      if (is_rz) {
        high_pressure_x = high_rz_pressure_face_flux->const_array(mfi);
        low_pressure_x = low_rz_pressure_face_flux->const_array(mfi);
      }
      cell_reduce_ops.eval(
          bx, cell_reduce_data,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) -> CellReduceTuple {
            if (marker(i, j, k, 0) != 0) {
              return {Real(0.0), Real(0.0), Real(0.0), Real(0.0),
                      Real(0.0)};
            }
            const Real delta = volume(i, j, k) *
                (limited(i, j, k, n) - low(i, j, k, n));
            Real auxiliary = Real(0.0);
#if (AMREX_SPACEDIM == 2)
            if (is_rz) {
              auto selected_difference = [] AMREX_GPU_DEVICE(
                                             Real theta, Real high,
                                             Real low) noexcept {
                return theta > Real(0.0) ? theta * (high - low) : Real(0.0);
              };
              if (n == Closure::UMX) {
                const Real pressure_lo = selected_difference(
                    theta_x(i, j, k), high_pressure_x(i, j, k, 0),
                    low_pressure_x(i, j, k, 0));
                const Real pressure_hi = selected_difference(
                    theta_x(i + 1, j, k),
                    high_pressure_x(i + 1, j, k, 0),
                    low_pressure_x(i + 1, j, k, 0));
                auxiliary = stage_dt *
                    ((volume(i, j, k) / cell_size[0] -
                      area_x(i, j, k)) * pressure_lo -
                     (volume(i, j, k) / cell_size[0] -
                      area_x(i + 1, j, k)) * pressure_hi);
              } else if (prob_lo[0] == Real(0.0) &&
                         i == radial_axis_face_index) {
                const Real axis_metric_flux = selected_difference(
                    theta_x(i, j, k), high_fx(i, j, k, n),
                    low_fx(i, j, k, n));
                auxiliary = volume(i, j, k) * stage_dt * Real(2.0) /
                            (cell_size[0] * cell_size[0]) *
                            axis_metric_flux;
              }
            }
#endif
            return {delta, amrex::Math::abs(delta), auxiliary,
                    amrex::Math::abs(auxiliary),
                    volume(i, j, k) * amrex::Math::abs(low(i, j, k, n))};
          });
    }
    const auto cell_reduced = cell_reduce_data.value(cell_reduce_ops);
    Real cell_integral = amrex::get<0>(cell_reduced);
    Real cell_absolute_integral = amrex::get<1>(cell_reduced);
    Real auxiliary_integral = amrex::get<2>(cell_reduced);
    Real auxiliary_absolute_integral = amrex::get<3>(cell_reduced);
    Real base_absolute_integral = amrex::get<4>(cell_reduced);
    ParallelDescriptor::ReduceRealSum(cell_integral);
    ParallelDescriptor::ReduceRealSum(cell_absolute_integral);
    ParallelDescriptor::ReduceRealSum(auxiliary_integral);
    ParallelDescriptor::ReduceRealSum(auxiliary_absolute_integral);
    ParallelDescriptor::ReduceRealSum(base_absolute_integral);

    Real exposed_face_integral = Real(0.0);
    Real exposed_face_absolute_integral = Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      auto owner = face_theta[dir]->OwnerMask(geom.periodicity());
      ReduceOps<ReduceOpSum, ReduceOpSum> face_reduce_ops;
      ReduceData<Real, Real> face_reduce_data(face_reduce_ops);
      using FaceReduceTuple = typename decltype(face_reduce_data)::Type;
      for (MFIter mfi(*face_theta[dir], false); mfi.isValid(); ++mfi) {
        const Box& face_box = mfi.validbox();
        const auto theta = face_theta[dir]->const_array(mfi);
        const auto high_flux = high_face_flux[dir]->const_array(mfi);
        const auto low_flux = low_face_flux[dir]->const_array(mfi);
        const auto owned = owner->const_array(mfi);
        const auto marker = marker_mf.const_array(mfi);
        const auto valid = valid_on_level.const_array(mfi);
        const auto area = metric_face_area[dir]->const_array(mfi);
        face_reduce_ops.eval(
            face_box, face_reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> FaceReduceTuple {
              if (owned(i, j, k) == 0) {
                return {Real(0.0), Real(0.0)};
              }
              int li = i;
              int lj = j;
              int lk = k;
              if (dir == 0) --li;
#if (AMREX_SPACEDIM >= 2)
              if (dir == 1) --lj;
#endif
#if (AMREX_SPACEDIM == 3)
              if (dir == 2) --lk;
#endif
              const bool left_inside = li >= ilo && li <= ihi
#if (AMREX_SPACEDIM >= 2)
                  && lj >= jlo && lj <= jhi
#endif
#if (AMREX_SPACEDIM == 3)
                  && lk >= klo && lk <= khi
#endif
                  ;
              const bool right_inside = i >= ilo && i <= ihi
#if (AMREX_SPACEDIM >= 2)
                  && j >= jlo && j <= jhi
#endif
#if (AMREX_SPACEDIM == 3)
                  && k >= klo && k <= khi
#endif
                  ;
              const bool left_active = left_inside &&
                  valid(li, lj, lk) != 0 && marker(li, lj, lk, 0) == 0;
              const bool right_active = right_inside &&
                  valid(i, j, k) != 0 && marker(i, j, k, 0) == 0;
              if (left_active == right_active) {
                return {Real(0.0), Real(0.0)};
              }
              const Real orientation =
                  (right_active ? Real(1.0) : Real(0.0)) -
                  (left_active ? Real(1.0) : Real(0.0));
              const Real theta_value = theta(i, j, k);
              const Real correction = theta_value > Real(0.0)
                  ? theta_value *
                        (high_flux(i, j, k, n) - low_flux(i, j, k, n))
                  : Real(0.0);
              const Real contribution =
                  orientation * stage_dt * area(i, j, k) * correction;
              return {contribution, amrex::Math::abs(contribution)};
            });
      }
      const auto face_reduced = face_reduce_data.value(face_reduce_ops);
      exposed_face_integral += amrex::get<0>(face_reduced);
      exposed_face_absolute_integral += amrex::get<1>(face_reduced);
    }
    ParallelDescriptor::ReduceRealSum(exposed_face_integral);
    ParallelDescriptor::ReduceRealSum(exposed_face_absolute_integral);

    // Do not feed a non-finite reduction into amrex::max: depending on
    // argument order, a NaN can otherwise leave an earlier finite maximum
    // unchanged and silently pass the post-loop gate.
    if (!std::isfinite(cell_integral) ||
        !std::isfinite(cell_absolute_integral) ||
        !std::isfinite(auxiliary_integral) ||
        !std::isfinite(auxiliary_absolute_integral) ||
        !std::isfinite(exposed_face_integral) ||
        !std::isfinite(exposed_face_absolute_integral) ||
        !std::isfinite(base_absolute_integral)) {
      amrex::Print()
          << "[IBM-Positivity-Limiter] non-finite shared-face closure "
          << "reduction for component=" << n << '\n';
      amrex::Abort(
          "local shared-face limiter produced a non-finite full-cell "
          "operator-closure reduction");
    }

    const Real residual =
        cell_integral - exposed_face_integral - auxiliary_integral;
    const Real scale = cell_absolute_integral +
        exposed_face_absolute_integral + auxiliary_absolute_integral;
    const Real absolute_residual = amrex::Math::abs(residual);
    Real relative = std::numeric_limits<Real>::max();
    if (absolute_residual == Real(0.0)) {
      relative = Real(0.0);
    } else if (scale > Real(0.0) &&
               scale >= absolute_residual /
                            std::numeric_limits<Real>::max()) {
      relative = absolute_residual / scale;
    }
    const Real conservation_relative_tolerance = amrex::max(
        Real(1.0e-10),
        Real(4096.0) * std::numeric_limits<Real>::epsilon());
    // When the antidiffusive correction is tiny, its global sum is obtained
    // by subtracting two nearly identical stage states.  Admit only the
    // roundoff implied by that subtraction, measured against the full low
    // state; this is an absolute allowance, not a relaxation of operator
    // closure at finite correction magnitude.
    // The production qualification uses double precision.  Keep a much
    // tighter fail-closed allowance if a single-precision executable is
    // nevertheless built, instead of silently accepting O(1e-3) base-state
    // closure errors from a large epsilon multiplier.
    constexpr Real base_roundoff_multiplier =
        std::numeric_limits<Real>::digits >= 53 ? Real(16384.0) : Real(64.0);
    const Real roundoff_allowance = base_roundoff_multiplier *
        std::numeric_limits<Real>::epsilon() * base_absolute_integral;
    const Real allowed_absolute_residual =
        conservation_relative_tolerance * scale + roundoff_allowance;
    if (!std::isfinite(residual) || !std::isfinite(scale) ||
        !std::isfinite(allowed_absolute_residual) ||
        absolute_residual > allowed_absolute_residual) {
      amrex::Print()
          << "[IBM-Positivity-Limiter] shared-face operator closure failed: "
          << "component=" << n
          << " absolute=" << absolute_residual
          << " correction_scale=" << scale
          << " base_scale=" << base_absolute_integral
          << " relative=" << relative
          << " allowed_absolute=" << allowed_absolute_residual << '\n';
      amrex::Abort(
          "local shared-face limiter failed its full-cell conservative/"
          "paired-pressure closure gate");
    }
    result.conservation_residual_linf = amrex::max(
        result.conservation_residual_linf, absolute_residual);
    result.conservation_relative_linf = amrex::max(
        result.conservation_relative_linf, relative);
  }

  if (!std::isfinite(result.conservation_residual_linf) ||
      !std::isfinite(result.conservation_relative_linf)) {
    amrex::Print()
        << "[IBM-Positivity-Limiter] non-finite shared-face closure summary: "
        << "absolute=" << result.conservation_residual_linf
        << " relative=" << result.conservation_relative_linf << '\n';
    amrex::Abort(
        "local shared-face limiter produced a non-finite closure summary");
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    auto owner = face_theta[dir]->OwnerMask(geom.periodicity());
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMin>
        reduce_ops;
    ReduceData<Long, Long, Long, Long, Real> reduce_data(reduce_ops);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(*face_theta[dir], false); mfi.isValid(); ++mfi) {
      const Box& face_box = mfi.validbox();
      const auto theta = face_theta[dir]->const_array(mfi);
      const auto owned = owner->const_array(mfi);
      const auto marker = marker_mf.const_array(mfi);
      const auto valid = valid_on_level.const_array(mfi);
      reduce_ops.eval(
          face_box, reduce_data,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            if (owned(i, j, k) == 0) {
              return {Long(0), Long(0), Long(0), Long(0), Real(1.0)};
            }
            int li = i;
            int lj = j;
            int lk = k;
            if (dir == 0) --li;
#if (AMREX_SPACEDIM >= 2)
            if (dir == 1) --lj;
#endif
#if (AMREX_SPACEDIM == 3)
            if (dir == 2) --lk;
#endif
            const bool left_inside = li >= ilo && li <= ihi
#if (AMREX_SPACEDIM >= 2)
                && lj >= jlo && lj <= jhi
#endif
#if (AMREX_SPACEDIM == 3)
                && lk >= klo && lk <= khi
#endif
                ;
            const bool right_inside = i >= ilo && i <= ihi
#if (AMREX_SPACEDIM >= 2)
                && j >= jlo && j <= jhi
#endif
#if (AMREX_SPACEDIM == 3)
                && k >= klo && k <= khi
#endif
                ;
            const bool left_on_level =
                left_inside && valid(li, lj, lk) != 0;
            const bool right_on_level =
                right_inside && valid(i, j, k) != 0;
            const bool left_active =
                left_on_level && marker(li, lj, lk, 0) == 0;
            const bool right_active =
                right_on_level && marker(i, j, k, 0) == 0;
            if (!left_active && !right_active) {
              return {Long(0), Long(0), Long(0), Long(0), Real(1.0)};
            }
            const Real value = theta(i, j, k);
            const bool limited_face = value < Real(1.0) - Real(1.0e-13);
            const bool fluid_face = left_active && right_active;
            const bool crossing_face = left_on_level && right_on_level &&
                (left_active != right_active);
            return {Long(1), limited_face ? Long(1) : Long(0),
                    limited_face && fluid_face ? Long(1) : Long(0),
                    limited_face && crossing_face ? Long(1) : Long(0),
                    value};
          });
    }
    const auto reduced = reduce_data.value(reduce_ops);
    result.active_faces += amrex::get<0>(reduced);
    result.limited_faces += amrex::get<1>(reduced);
    result.limited_fluid_faces += amrex::get<2>(reduced);
    result.limited_crossing_faces += amrex::get<3>(reduced);
    result.minimum_theta =
        amrex::min(result.minimum_theta, amrex::get<4>(reduced));
  }
  ParallelDescriptor::ReduceLongSum(result.active_faces);
  ParallelDescriptor::ReduceLongSum(result.limited_faces);
  ParallelDescriptor::ReduceLongSum(result.limited_fluid_faces);
  ParallelDescriptor::ReduceLongSum(result.limited_crossing_faces);
  ParallelDescriptor::ReduceRealMin(result.minimum_theta);

  const auto limited_stats =
      collectStageAdmissibility(limited_candidate, level, closure);
  result.limited_state_admissible = limited_stats.admissible();
  printStageAdmissibility("local-shared-face", limited_stats);
  amrex::Print()
      << "  active_faces=" << result.active_faces
      << " limited_faces=" << result.limited_faces
      << " limited_fluid_faces=" << result.limited_fluid_faces
      << " limited_crossing_faces=" << result.limited_crossing_faces
      << " minimum_theta=" << result.minimum_theta
      << " conservation_residual_linf="
      << result.conservation_residual_linf
      << " conservation_relative_linf="
      << result.conservation_relative_linf
      << " verdict="
      << (limited_stats.admissible() ? "LOCAL_SHARED_FACE_PASS"
                                     : "LOCAL_SHARED_FACE_FAIL")
      << '\n'
      << "=================================================================\n";
  if (!result.limited_state_admissible) {
    if (abort_on_failure) {
      amrex::Abort(
          "local shared-face limiter produced an inadmissible stage state");
    }
    return result;
  }
  MultiFab::Copy(high_candidate, limited_candidate, 0, 0, Closure::NCONS, 0);
  if (replace_high_order_flux) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      for (MFIter mfi(*high_face_flux[dir], false); mfi.isValid(); ++mfi) {
        const Box& face_box = mfi.validbox();
        const auto high_flux = high_face_flux[dir]->array(mfi);
        const auto low_flux = low_face_flux[dir]->const_array(mfi);
        const auto theta = face_theta[dir]->const_array(mfi);
        ParallelFor(
            face_box, Closure::NCONS,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
              const Real theta_value = theta(i, j, k);
              const Real low_value = low_flux(i, j, k, n);
              const Real high_value = high_flux(i, j, k, n);
              high_flux(i, j, k, n) =
                  theta_value > Real(0.0) && std::isfinite(high_value)
                  ? low_value + theta_value * (high_value - low_value)
                  : low_value;
            });
      }
      high_face_flux[dir]->OverrideSync(geom.periodicity());
    }
    if (is_rz) {
      for (MFIter mfi(*high_rz_pressure_face_flux, false); mfi.isValid();
           ++mfi) {
        const Box& face_box = mfi.validbox();
        const auto high_pressure =
            high_rz_pressure_face_flux->array(mfi);
        const auto low_pressure =
            low_rz_pressure_face_flux->const_array(mfi);
        const auto theta = face_theta[0]->const_array(mfi);
        ParallelFor(
            face_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const Real theta_value = theta(i, j, k);
              const Real low_value = low_pressure(i, j, k, 0);
              const Real high_value = high_pressure(i, j, k, 0);
              high_pressure(i, j, k, 0) =
                  theta_value > Real(0.0) && std::isfinite(high_value)
                      ? low_value + theta_value * (high_value - low_value)
                      : low_value;
            });
      }
      high_rz_pressure_face_flux->OverrideSync(geom.periodicity());
    }
  }
  return result;
}

template <typename Closure>
inline void replaceFaceFluxWithGlobalBlend(
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& high_face_flux,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& low_face_flux,
    amrex::MultiFab* high_rz_pressure_face_flux,
    const amrex::MultiFab* low_rz_pressure_face_flux,
    const amrex::Geometry& geom, amrex::Real theta) {
  using namespace amrex;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (high_face_flux[dir] == nullptr || low_face_flux[dir] == nullptr) {
      amrex::Abort("global flux blend requires both flux families");
    }
    for (MFIter mfi(*high_face_flux[dir], false); mfi.isValid(); ++mfi) {
      const Box& face_box = mfi.validbox();
      const auto high = high_face_flux[dir]->array(mfi);
      const auto low = low_face_flux[dir]->const_array(mfi);
      ParallelFor(
          face_box, Closure::NCONS,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            const Real low_value = low(i, j, k, n);
            const Real high_value = high(i, j, k, n);
            high(i, j, k, n) = theta > Real(0.0) &&
                                      std::isfinite(high_value)
                ? low_value + theta * (high_value - low_value)
                : low_value;
          });
    }
    high_face_flux[dir]->OverrideSync(geom.periodicity());
  }
  if (geom.IsRZ()) {
#if (AMREX_SPACEDIM == 2)
    if (high_rz_pressure_face_flux == nullptr ||
        low_rz_pressure_face_flux == nullptr) {
      amrex::Abort(
          "paired R-Z global flux blend requires both pressure families");
    }
    for (MFIter mfi(*high_rz_pressure_face_flux, false); mfi.isValid(); ++mfi) {
      const Box& face_box = mfi.validbox();
      const auto high = high_rz_pressure_face_flux->array(mfi);
      const auto low = low_rz_pressure_face_flux->const_array(mfi);
      ParallelFor(
          face_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real low_value = low(i, j, k, 0);
            const Real high_value = high(i, j, k, 0);
            high(i, j, k, 0) = theta > Real(0.0) &&
                                       std::isfinite(high_value)
                ? low_value + theta * (high_value - low_value)
                : low_value;
          });
    }
    high_rz_pressure_face_flux->OverrideSync(geom.periodicity());
#else
    amrex::Abort("paired R-Z global flux blend requires a 2-D build");
#endif
  } else if (high_rz_pressure_face_flux != nullptr ||
             low_rz_pressure_face_flux != nullptr) {
    amrex::Abort(
        "Cartesian global flux blend cannot receive R-Z pressure companions");
  }
}

// Fallback from the identical high/low candidates and face fluxes used by the
// local limiter. The result is first assembled in scratch storage, so a failed
// global blend cannot contaminate the subsequent all-low-order fallback.
template <typename Closure>
inline bool tryGlobalThetaFluxFallback(
    amrex::MultiFab& high_candidate, const amrex::MultiFab& low_candidate,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& high_face_flux,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& low_face_flux,
    amrex::MultiFab* high_rz_pressure_face_flux,
    const amrex::MultiFab* low_rz_pressure_face_flux,
    const amrex::Geometry& geom, int level, const Closure& closure,
    amrex::Real& theta_used) {
  using namespace amrex;
  const Real theta_max = maximumGlobalAdmissibleTheta(
      low_candidate, high_candidate, level, closure);
  if (!(theta_max > Real(0.0)) || !std::isfinite(theta_max)) {
    theta_used = Real(0.0);
    return false;
  }
  theta_used = amrex::min(
      Real(1.0), theta_max * (Real(1.0) - Real(1.0e-12)));

  MultiFab candidate(
      low_candidate.boxArray(), low_candidate.DistributionMap(),
      Closure::NCONS, 0, MFInfo().SetArena(The_Async_Arena()));
  const auto& marker_mf = *IBM::ib.bmf_a[level];
  MultiFab::Copy(candidate, low_candidate, 0, 0, Closure::NCONS, 0);
  for (MFIter mfi(candidate, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const auto candidate_a = candidate.array(mfi);
    const auto high = high_candidate.const_array(mfi);
    const auto low = low_candidate.const_array(mfi);
    const auto marker = marker_mf.const_array(mfi);
    ParallelFor(
        bx, Closure::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          if (marker(i, j, k, 0) != 0) return;
          const Real low_value = low(i, j, k, n);
          const Real high_value = high(i, j, k, n);
          candidate_a(i, j, k, n) = std::isfinite(high_value)
              ? low_value + theta_used * (high_value - low_value)
              : low_value;
        });
  }
  const auto stats = collectStageAdmissibility(candidate, level, closure);
  if (!stats.admissible()) return false;

  MultiFab::Copy(high_candidate, candidate, 0, 0, Closure::NCONS, 0);
  replaceFaceFluxWithGlobalBlend<Closure>(
      high_face_flux, low_face_flux, high_rz_pressure_face_flux,
      low_rz_pressure_face_flux, geom, theta_used);
  return true;
}

template <typename Closure>
inline void replaceWithAllLowOrderFluxUpdate(
    amrex::MultiFab& high_candidate, const amrex::MultiFab& low_candidate,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& high_face_flux,
    const std::array<amrex::MultiFab*, AMREX_SPACEDIM>& low_face_flux,
    amrex::MultiFab* high_rz_pressure_face_flux,
    const amrex::MultiFab* low_rz_pressure_face_flux,
    const amrex::Geometry& geom) {
  amrex::MultiFab::Copy(
      high_candidate, low_candidate, 0, 0, Closure::NCONS, 0);
  replaceFaceFluxWithGlobalBlend<Closure>(
      high_face_flux, low_face_flux, high_rz_pressure_face_flux,
      low_rz_pressure_face_flux, geom, amrex::Real(0.0));
}

}  // namespace IBM::positivity

#endif
