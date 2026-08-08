#include <CNS.h>
#include <prob.h>

using namespace amrex;

#if NUM_SPECIES > 1
void CNS::clip_species_state(MultiFab& statemf)
{
  BL_PROFILE("CNS::clip_species_state()");

  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;

  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    const Array4<Real> cons = statemf.array(mfi);
    const Box bx = mfi.growntilebox(0);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      Real rho_species = Real(0.0);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        rho_species += amrex::max(cons(i, j, k, cls_h.UFS + n), Real(0.0));
      }

      const Real rho_species_inv = Real(1.0) / rho_species;
      Real mass_fraction[NUM_SPECIES];
      Real mass_fraction_sum = Real(0.0);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        mass_fraction[n] =
            amrex::max(cons(i, j, k, cls_h.UFS + n), Real(0.0)) *
            rho_species_inv;
        mass_fraction_sum += mass_fraction[n];
      }
      for (int n = 0; n < NUM_SPECIES; ++n) {
        mass_fraction[n] /= mass_fraction_sum;
        cons(i, j, k, cls_h.UFS + n) =
            rho_species * mass_fraction[n];
      }
    });
  }
}
#endif
