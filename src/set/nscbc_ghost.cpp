#include <CNS.h>
#include <nscbc.h>

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
