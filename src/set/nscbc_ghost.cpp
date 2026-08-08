#include <CNS.h>
#include <nscbc.h>

#include <AMReX_GpuAtomic.H>

using namespace amrex;

namespace {

Box physical_ghost_slab(const Box& domain, int dir, int side, int ng)
{
  Box slab(domain);
  if (side < 0) {
    slab.setSmall(dir, domain.smallEnd(dir) - ng);
    slab.setBig(dir, domain.smallEnd(dir) - 1);
  } else {
    slab.setSmall(dir, domain.bigEnd(dir) + 1);
    slab.setBig(dir, domain.bigEnd(dir) + ng);
  }
  return slab;
}

} // namespace

void CNS::initialize_nscbc_ghost_state(Real time)
{
  if (!use_nscbc) return;

  if (nscbc_parm.ghost_update_model !=
      nscbc::GHOST_UPDATE_PERSISTENT_ODE) {
    nscbc_ghost_state.reset();
    nscbc_ghost_initialized = false;
    return;
  }

  const int ncons = h_prob_closures->NCONS;
  const int ng = h_prob_closures->NGHOST;

  // Restart constructs AmrLevels through the default constructor, while a
  // regrid may replace grids/dmap in an existing level.  Allocate here, after
  // both paths have installed their final layout, so every lifecycle path uses
  // a correctly shaped persistent state.
  nscbc_ghost_state = std::make_unique<MultiFab>(
    grids, dmap, ncons, ng, MFInfo(), Factory());
  nscbc_ghost_state->setVal(Real(0.0));

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if ((nscbc_lo[dir] > 0 || nscbc_hi[dir] > 0) && Geom().isPeriodic(dir)) {
      amrex::Abort("NSCBC cannot be enabled on a periodic direction");
    }
  }

  MultiFab initial(grids, dmap, ncons, ng, MFInfo(), Factory());
  FillPatch(*this, initial, ng, time, State_Type, 0, ncons);
#if (AMREX_SPACEDIM < 3)
  initial.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, ng);
#endif
  MultiFab::Copy(*nscbc_ghost_state, initial, 0, 0, ncons, ng);
  nscbc_ghost_initialized = true;
}

void CNS::copy_nscbc_ghost_to_state(MultiFab& state,
                                     MultiFab const& ghost) const
{
  if (!use_nscbc) return;
  const int ncons = h_prob_closures->NCONS;
  const int ng = h_prob_closures->NGHOST;
  const Box domain = Geom().Domain();

  for (MFIter mfi(state, false); mfi.isValid(); ++mfi) {
    auto const& dst = state.array(mfi);
    auto const& src = ghost.const_array(mfi);
    const Box fabbox = mfi.growntilebox(ng);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (nscbc_lo[dir] > 0) {
        const Box b = fabbox & physical_ghost_slab(domain, dir, -1, ng);
        if (b.ok()) {
          ParallelFor(b, ncons,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
              dst(i, j, k, n) = src(i, j, k, n);
            });
        }
      }
      if (nscbc_hi[dir] > 0) {
        const Box b = fabbox & physical_ghost_slab(domain, dir, +1, ng);
        if (b.ok()) {
          ParallelFor(b, ncons,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
              dst(i, j, k, n) = src(i, j, k, n);
            });
        }
      }
    }
  }
}

void CNS::compute_nscbc_ghost_rhs(MultiFab& state, MultiFab& ghost_rhs) const
{
  if (!use_nscbc) return;
  const PROB::ProbClosures& cls_h = *h_prob_closures;
  const int ncons = cls_h.NCONS;
  const int ng = cls_h.NGHOST;
  const Box domain = Geom().Domain();
  const auto dxinv = Geom().InvCellSizeArray();
  const nscbc::Parm parm = nscbc_parm;

#if (AMREX_SPACEDIM < 3)
  // Thermodynamic conversion reads all three momenta in every dimensional
  // build.  Keep the inactive component deterministic before constructing q.
  state.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, state.nGrow());
#endif
  ghost_rhs.setVal(Real(0.0), 0, ncons, ng);

  for (MFIter mfi(state, false); mfi.isValid(); ++mfi) {
    const Box fabbox = mfi.growntilebox(ng);
    FArrayBox primfab(fabbox, cls_h.NPRIM, The_Async_Arena());
    auto const& cons = state.array(mfi);
    auto const& q = primfab.array();
    auto const& rhs = ghost_rhs.array(mfi);
    cls_h.cons2prims(mfi, cons, q);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (nscbc_lo[dir] > 0) {
        const Box b = fabbox & physical_ghost_slab(domain, dir, -1, ng);
        const int type = nscbc_lo[dir];
        if (b.ok()) {
          ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            nscbc::add_lodi_rhs_to_cons<PROB::ProbClosures>(
              IntVect(AMREX_D_DECL(i, j, k)), dir, +1,
              domain.smallEnd(dir), dxinv, q, rhs, type, parm);
          });
        }
      }
      if (nscbc_hi[dir] > 0) {
        const Box b = fabbox & physical_ghost_slab(domain, dir, +1, ng);
        const int type = nscbc_hi[dir];
        if (b.ok()) {
          ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            nscbc::add_lodi_rhs_to_cons<PROB::ProbClosures>(
              IntVect(AMREX_D_DECL(i, j, k)), dir, -1,
              domain.bigEnd(dir), dxinv, q, rhs, type, parm);
          });
        }
      }
    }
  }
}

void CNS::fill_nscbc_gc_ghosts(MultiFab& state) const
{
  if (!use_nscbc ||
      nscbc_parm.ghost_update_model !=
        nscbc::GHOST_UPDATE_STAGE_LOCAL_GC) {
    return;
  }

#if (AMREX_SPACEDIM != 2)
  amrex::Abort("stage-local GC-NSCBC is currently restricted to 2-D");
#else
  if (Geom().IsRZ()) {
    amrex::Abort(
      "stage-local GC-NSCBC is currently restricted to Cartesian geometry");
  }
  if (level != 0) {
    amrex::Abort(
      "stage-local GC-NSCBC has not been certified on AMR refined levels");
  }
  if constexpr (NUM_SPECIES != 1) {
    amrex::Abort(
      "stage-local GC-NSCBC is currently restricted to a single-species "
      "calorically perfect gas");
  }

  const Box domain = Geom().Domain();
  const int ng = h_prob_closures->NGHOST;
  if (ng < 1 || ng > 3) {
    amrex::Abort("stage-local GC-NSCBC supports one to three ghost layers");
  }

  int active_count = 0;
  int active_dir = -1;
  int active_side = 0;
  int active_type = 0;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (nscbc_lo[dir] > 0) {
      ++active_count;
      active_dir = dir;
      active_side = -1;
      active_type = nscbc_lo[dir];
    }
    if (nscbc_hi[dir] > 0) {
      ++active_count;
      active_dir = dir;
      active_side = +1;
      active_type = nscbc_hi[dir];
    }
  }
  if (active_count != 1 || (active_type != 2 && active_type != 3)) {
    amrex::Abort(
      "stage-local GC-NSCBC currently requires exactly one type-2 or type-3 "
      "outflow boundary");
  }
  const int tangent_dir = 1 - active_dir;
  if (!Geom().isPeriodic(tangent_dir)) {
    amrex::Abort(
      "stage-local GC-NSCBC currently requires a periodic tangential "
      "direction; physical-corner ownership has not yet been certified");
  }
  if (domain.length(active_dir) < 3) {
    amrex::Abort(
      "stage-local GC-NSCBC requires at least three interior cells in the "
      "boundary-normal direction");
  }

#if (AMREX_SPACEDIM < 3)
  state.setVal(Real(0.0), PROB::ProbClosures::UMZ, 1, state.nGrow());
#endif

  const auto dx = Geom().CellSizeArray();
  const nscbc::Parm parm = nscbc_parm;
  const PROB::ProbClosures* cls_d = d_prob_closures;
  Gpu::DeviceScalar<int> failure(0);
  int* failure_ptr = failure.dataPtr();

  for (MFIter mfi(state, false); mfi.isValid(); ++mfi) {
    const Box valid = mfi.validbox();
    const int boundary_index =
      (active_side < 0) ? domain.smallEnd(active_dir)
                        : domain.bigEnd(active_dir);
    if (boundary_index < valid.smallEnd(active_dir) ||
        boundary_index > valid.bigEnd(active_dir)) {
      continue;
    }

    Box boundary_box(valid);
    boundary_box.setSmall(active_dir, boundary_index);
    boundary_box.setBig(active_dir, boundary_index);

    const Box fabbox = mfi.growntilebox(ng);
    FArrayBox primfab(fabbox, h_prob_closures->NPRIM,
                      The_Async_Arena());
    auto const& cons = state.array(mfi);
    auto const& q = primfab.array();
    h_prob_closures->cons2prims(fabbox, cons, q);

    const int dir = active_dir;
    const int tdir = tangent_dir;
    const int side = active_side;
    const int type = active_type;
    const Real h = dx[dir];
    const Real inv_2h = Real(0.5) / h;
    const Real inv_2ht = Real(0.5) / dx[tdir];
    const Real normal_sign = (side < 0) ? Real(-1.0) : Real(1.0);
    const auto en = IntVect::TheDimensionVector(dir);
    const auto et = IntVect::TheDimensionVector(tdir);

    ParallelFor(boundary_box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        using PC = PROB::ProbClosures;
        const IntVect ib(AMREX_D_DECL(i, j, k));
        const IntVect im1 = ib - static_cast<int>(normal_sign) * en;
        const IntVect im2 = ib - 2 * static_cast<int>(normal_sign) * en;
        const IntVect tp = ib + et;
        const IntVect tm = ib - et;

        const Real rho = q(ib, PC::QRHO);
        const Real p = q(ib, PC::QPRES);
        const Real c = q(ib, PC::QC);
        const Real un = normal_sign * q(ib, nscbc::qvel<PC>(dir));
        const Real ut = q(ib, nscbc::qvel<PC>(tdir));

        const auto outward_derivative =
          [=] AMREX_GPU_DEVICE (int component, Real component_sign) noexcept {
            const Real qb = component_sign * q(ib, component);
            const Real q1 = component_sign * q(im1, component);
            const Real q2 = component_sign * q(im2, component);
            return (Real(3.0) * qb - Real(4.0) * q1 + q2) * inv_2h;
          };

        const Real rho_n_i = outward_derivative(PC::QRHO, Real(1.0));
        const Real p_n_i = outward_derivative(PC::QPRES, Real(1.0));
        const Real un_n_i =
          outward_derivative(nscbc::qvel<PC>(dir), normal_sign);
        const Real ut_n_i =
          outward_derivative(nscbc::qvel<PC>(tdir), Real(1.0));
        const Real tangent_divergence =
          (q(tp, nscbc::qvel<PC>(tdir)) -
           q(tm, nscbc::qvel<PC>(tdir))) * inv_2ht;

        const bool finite_boundary = std::isfinite(rho) &&
          std::isfinite(p) && std::isfinite(c) && std::isfinite(un) &&
          std::isfinite(ut) && rho > parm.min_rho && p > parm.min_p &&
          c > Real(0.0);
        if (!finite_boundary || un <= Real(0.0)) {
          amrex::Gpu::Atomic::Max(failure_ptr, 1);
          return;
        }

        Real rho_n = rho_n_i;
        Real un_n = un_n_i;
        Real ut_n = ut_n_i;
        Real p_n = p_n_i;
        if (un < c) {
          Real pressure_relaxation = Real(0.0);
          if (type == 3) {
            pressure_relaxation = parm.sigma * c *
              amrex::max(Real(1.0) - parm.Mmax * parm.Mmax, Real(0.0)) *
              (p - parm.Ptarget) / parm.Lchar;
          }
          nscbc::giles_second_order_normal_derivatives(
            rho, c, un, rho_n_i, un_n_i, ut_n_i, p_n_i,
            tangent_divergence, pressure_relaxation,
            rho_n, un_n, ut_n, p_n);
        }

        const Real rho_in = q(im1, PC::QRHO);
        const Real p_in = q(im1, PC::QPRES);
        const Real un_in = normal_sign * q(im1, nscbc::qvel<PC>(dir));
        const Real ut_in = q(im1, nscbc::qvel<PC>(tdir));

        Real rho_g[3] = {Real(0.0), Real(0.0), Real(0.0)};
        Real p_g[3] = {Real(0.0), Real(0.0), Real(0.0)};
        Real un_g[3] = {Real(0.0), Real(0.0), Real(0.0)};
        Real ut_g[3] = {Real(0.0), Real(0.0), Real(0.0)};

        for (int layer = 1; layer <= ng; ++layer) {
          rho_g[layer - 1] = nscbc::gc_nscbc_value(
            layer, rho, rho_in, h * rho_n,
            rho_g[0], rho_g[1]);
          p_g[layer - 1] = nscbc::gc_nscbc_value(
            layer, p, p_in, h * p_n, p_g[0], p_g[1]);
          un_g[layer - 1] = nscbc::gc_nscbc_value(
            layer, un, un_in, h * un_n,
            un_g[0], un_g[1]);
          ut_g[layer - 1] = nscbc::gc_nscbc_value(
            layer, ut, ut_in, h * ut_n,
            ut_g[0], ut_g[1]);

          const Real rg = rho_g[layer - 1];
          const Real pg = p_g[layer - 1];
          const Real ung = un_g[layer - 1];
          const Real utg = ut_g[layer - 1];
          if (!std::isfinite(rg) || !std::isfinite(pg) ||
              !std::isfinite(ung) || !std::isfinite(utg) ||
              rg <= parm.min_rho || pg <= parm.min_p) {
            amrex::Gpu::Atomic::Max(failure_ptr, 2);
            return;
          }

          const IntVect ig = ib + layer * static_cast<int>(normal_sign) * en;
          Real velocity[3] = {Real(0.0), Real(0.0), Real(0.0)};
          velocity[dir] = normal_sign * ung;
          velocity[tdir] = utg;
          Real Y[NUM_SPECIES] = {Real(1.0)};
          Real internal_energy = Real(0.0);
          cls_d->RYP2E(rg, Y, pg, internal_energy);
          if (!std::isfinite(internal_energy) ||
              internal_energy <= Real(0.0)) {
            amrex::Gpu::Atomic::Max(failure_ptr, 3);
            return;
          }

          cons(ig, PC::UMX) = rg * velocity[0];
          cons(ig, PC::UMY) = rg * velocity[1];
          cons(ig, PC::UMZ) = Real(0.0);
          cons(ig, PC::UET) = rg * internal_energy + Real(0.5) * rg *
            (velocity[0] * velocity[0] + velocity[1] * velocity[1]);
          cons(ig, PC::URHO) = rg;
        }
      });
  }

  Gpu::streamSynchronize();
  const int failure_code = failure.dataValue();
  if (failure_code != 0) {
    amrex::Abort(
      "stage-local GC-NSCBC rejected a boundary stage: code=" +
      std::to_string(failure_code) +
      " (1=invalid/reverse boundary state, 2=invalid ghost primitive, "
      "3=invalid ghost internal energy)");
  }
#endif
}
