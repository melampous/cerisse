#ifndef IBM_TRIPWIRE_H
#define IBM_TRIPWIRE_H

// Lightweight IB/AMR pipeline tripwire for diagnosing state corruption.
//
// Scans all fluid cells (marker[0] == 0) in `S` and prints one summary line
// with: min(rho), min(p), max(|u|+c), n_bad (rho<=0 OR p<=0 OR NaN), n_total.
// Only fires when (step_lo <= step <= step_hi) so it's zero-cost outside.
//
// Conservative layout is hardcoded (indicies_t):
//   0=UMX  1=UMY  2=UMZ  3=UET  4=URHO
// Uses Gpu::Atomic scalars instead of ReduceOps to avoid interactions with
// nvcc's SFINAE in ibm_containers.

#include <AMReX_MultiFab.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <cmath>

template<class IBMFabArr>
inline void
ib_tripwire(const amrex::MultiFab& S,
            const IBMFabArr& ib_mf,
            const char* stage,
            int level, int step,
            int step_lo, int step_hi)
{
    if (step < step_lo || step > step_hi) return;

    using namespace amrex;
    constexpr int IMX   = 0;
    constexpr int IMY   = 1;
    constexpr int IET   = 3;
    constexpr int IRHO  = 4;
    constexpr Real GAMMA    = Real(1.4);
    constexpr Real GAMMA_M1 = Real(0.4);

    Gpu::DeviceScalar<Real> d_min_rho(Real( 1.0e30));
    Gpu::DeviceScalar<Real> d_min_p  (Real( 1.0e30));
    Gpu::DeviceScalar<Real> d_max_sig(Real(-1.0e30));
    Gpu::DeviceScalar<int>  d_n_bad  (0);
    Gpu::DeviceScalar<int>  d_n_tot  (0);

    Real* p_min_rho = d_min_rho.dataPtr();
    Real* p_min_p   = d_min_p.dataPtr();
    Real* p_max_sig = d_max_sig.dataPtr();
    int*  p_n_bad   = d_n_bad.dataPtr();
    int*  p_n_tot   = d_n_tot.dataPtr();

    for (MFIter mfi(S, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& state = S.const_array(mfi);
        auto const& mk    = ib_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mk(i,j,k,0) != 0) return;     // skip solid cells

            const Real rho = state(i,j,k, IRHO);
            const Real mx  = state(i,j,k, IMX);
            const Real my  = state(i,j,k, IMY);
            const Real et  = state(i,j,k, IET);

            const Real rho_safe = (rho > Real(1.0e-30)) ? rho : Real(1.0e-30);
            const Real ke = Real(0.5) * (mx*mx + my*my) / rho_safe;
            const Real p  = GAMMA_M1 * (et - ke);
            const Real p_safe = (p > Real(1.0e-30)) ? p : Real(1.0e-30);
            const Real c  = std::sqrt(GAMMA * p_safe / rho_safe);
            const Real umag = std::sqrt((mx*mx + my*my)) / rho_safe;
            const Real sig  = umag + c;

            const bool is_nan = !(rho == rho) || !(p == p) || !(et == et);
            const bool bad    = (rho <= Real(0)) || (p <= Real(0)) || is_nan;

            Gpu::Atomic::Min(p_min_rho, rho);
            Gpu::Atomic::Min(p_min_p,   p);
            Gpu::Atomic::Max(p_max_sig, sig);
            Gpu::Atomic::Add(p_n_tot,   1);
            if (bad) Gpu::Atomic::Add(p_n_bad, 1);
        });
    }

    Gpu::streamSynchronize();

    Real min_rho = d_min_rho.dataValue();
    Real min_p   = d_min_p.dataValue();
    Real max_sig = d_max_sig.dataValue();
    int  n_bad   = d_n_bad.dataValue();
    int  n_total = d_n_tot.dataValue();

    ParallelDescriptor::ReduceRealMin(min_rho);
    ParallelDescriptor::ReduceRealMin(min_p);
    ParallelDescriptor::ReduceRealMax(max_sig);
    ParallelDescriptor::ReduceIntSum(n_bad);
    ParallelDescriptor::ReduceIntSum(n_total);

    const char* tag = (n_bad > 0) ? " !!BAD" : "";
    amrex::Print() << "[TRIPWIRE step=" << step
                   << " lev="   << level
                   << " stage=" << stage << "]"
                   << " n="       << n_total
                   << " min_rho=" << min_rho
                   << " min_p="   << min_p
                   << " max_sig=" << max_sig
                   << " n_bad="   << n_bad
                   << tag << "\n";
}

#endif // IBM_TRIPWIRE_H
