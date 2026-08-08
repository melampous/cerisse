#ifndef IBM_RZ_MOMENTS_H_
#define IBM_RZ_MOMENTS_H_

#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace cerisse::ibm {

using amrex::GpuArray;
using amrex::Real;

// Physical R-Z cell bounds. The radial coordinate is non-negative; the
// omitted 2*pi factor cancels from every normalized annular moment.
struct RZCellBounds {
    Real r_lo;
    Real r_hi;
    Real z_lo;
    Real z_hi;
};

enum RZQuadraticMoment : int {
    rz_one = 0,
    rz_z,
    rz_r,
    rz_z2,
    rz_zr,
    rz_r2,
    rz_quadratic_moment_count
};

// Return exact annular cell averages of {1,z,r,z^2,zr,r^2}:
//
//   <q>_RZ = integral(q r dz dr) / integral(r dz dr).
//
// The centroid and variance formulas avoid subtracting nearly equal powers
// far from the axis. In the first ring they give <r>=2*dr/3 and
// Var(r)=dr^2/18, rather than the Cartesian midpoint values.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool rz_annular_quadratic_moments(
    const RZCellBounds& cell,
    GpuArray<Real, rz_quadratic_moment_count>& moments)
{
    for (Real& value : moments) value = Real(0.0);

    const Real dr = cell.r_hi - cell.r_lo;
    const Real dz = cell.z_hi - cell.z_lo;
    const Real rc = Real(0.5) * (cell.r_lo + cell.r_hi);
    const Real zc = Real(0.5) * (cell.z_lo + cell.z_hi);
    if (!(cell.r_lo >= Real(0.0)) || !(dr > Real(0.0)) ||
        !(dz > Real(0.0)) || !(rc > Real(0.0)) ||
        !amrex::Math::isfinite(cell.r_lo + cell.r_hi + cell.z_lo +
                              cell.z_hi)) {
        return false;
    }

    const Real mean_r = rc + dr * dr / (Real(12.0) * rc);
    const Real mean_z = zc;
    const Real mean_r2 = rc * rc + dr * dr / Real(4.0);
    const Real mean_z2 = zc * zc + dz * dz / Real(12.0);

    moments[rz_one] = Real(1.0);
    moments[rz_z] = mean_z;
    moments[rz_r] = mean_r;
    moments[rz_z2] = mean_z2;
    moments[rz_zr] = mean_z * mean_r;
    moments[rz_r2] = mean_r2;
    return true;
}

// Return exact annular averages of the local BI-CWLS basis
// {1,s,t,s^2,st,t^2}. The local coordinates are normalized by length_scale,
// with frame rows ordered as (normal,tangent) in physical (r,z) components.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool rz_annular_local_quadratic_moments(
    const RZCellBounds& cell, const GpuArray<Real, 2>& wall_point,
    const GpuArray<Real, 2>& normal_in,
    const GpuArray<Real, 2>& tangent_in, Real length_scale,
    GpuArray<Real, rz_quadratic_moment_count>& moments)
{
    for (Real& value : moments) value = Real(0.0);
    if (!(length_scale > Real(0.0)) ||
        !amrex::Math::isfinite(length_scale)) {
        return false;
    }

    GpuArray<Real, 2> normal = normal_in;
    GpuArray<Real, 2> tangent = tangent_in;
    const Real normal_norm2 =
        normal[0] * normal[0] + normal[1] * normal[1];
    const Real tangent_norm2 =
        tangent[0] * tangent[0] + tangent[1] * tangent[1];
    if (!(normal_norm2 > Real(0.0)) ||
        !(tangent_norm2 > Real(0.0))) {
        return false;
    }
    const Real inv_normal = Real(1.0) / std::sqrt(normal_norm2);
    const Real inv_tangent = Real(1.0) / std::sqrt(tangent_norm2);
    for (int d = 0; d < 2; ++d) {
        normal[d] *= inv_normal;
        tangent[d] *= inv_tangent;
    }
    if (amrex::Math::abs(normal[0] * tangent[0] +
                         normal[1] * tangent[1]) > Real(1.0e-10)) {
        return false;
    }

    GpuArray<Real, rz_quadratic_moment_count> global{};
    if (!rz_annular_quadratic_moments(cell, global)) return false;

    const Real mean_r = global[rz_r];
    const Real mean_z = global[rz_z];
    const Real dr = cell.r_hi - cell.r_lo;
    const Real dz = cell.z_hi - cell.z_lo;
    const Real rc = Real(0.5) * (cell.r_lo + cell.r_hi);

    // Stable central moments for the separable measure r dr dz.
    const Real variance_r =
        dr * dr / Real(12.0) -
        dr * dr * dr * dr / (Real(144.0) * rc * rc);
    const Real variance_z = dz * dz / Real(12.0);
    const Real inverse_scale = Real(1.0) / length_scale;
    const Real inverse_scale2 = inverse_scale * inverse_scale;

    const Real mean_s =
        (normal[0] * (mean_r - wall_point[0]) +
         normal[1] * (mean_z - wall_point[1])) * inverse_scale;
    const Real mean_t =
        (tangent[0] * (mean_r - wall_point[0]) +
         tangent[1] * (mean_z - wall_point[1])) * inverse_scale;
    const Real covariance_ss =
        (normal[0] * normal[0] * variance_r +
         normal[1] * normal[1] * variance_z) * inverse_scale2;
    const Real covariance_st =
        (normal[0] * tangent[0] * variance_r +
         normal[1] * tangent[1] * variance_z) * inverse_scale2;
    const Real covariance_tt =
        (tangent[0] * tangent[0] * variance_r +
         tangent[1] * tangent[1] * variance_z) * inverse_scale2;

    moments[rz_one] = Real(1.0);
    moments[rz_z] = mean_s;
    moments[rz_r] = mean_t;
    moments[rz_z2] = mean_s * mean_s + covariance_ss;
    moments[rz_zr] = mean_s * mean_t + covariance_st;
    moments[rz_r2] = mean_t * mean_t + covariance_tt;
    return amrex::Math::isfinite(
        moments[rz_z] + moments[rz_r] + moments[rz_z2] +
        moments[rz_zr] + moments[rz_r2]);
}

} // namespace cerisse::ibm

#endif
