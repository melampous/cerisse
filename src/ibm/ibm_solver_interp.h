#ifndef IBM_SOLVER_INTERP_H_
#define IBM_SOLVER_INTERP_H_
// ============================================================================
// ibm_solver_interp.h — Interpolation, extrapolation, and coordinate-transform helpers
//
// This file contains in-class definitions for private member
// functions of ibm_solver_t.  It is included inside the class body in ibm_solver.h
// and must NOT be included independently.
//
// Contents:
//   1. valid_mirror              — Count fluid cells in an interpolation stencil
//   2. search_optimal_image_point— Best first image-point placement (multiple attempts)
//   3. search_image_point        — Subsequent image-point placement (single step)
//   4. computeIPweights          — Bi-/tri-linear interpolation weights
//   5. interpolateIMs            — Interpolate primitives at image points
//   6. extrapolate               — Lagrangian extrapolation to ghost point
//   7. global2local / local2global — Velocity coordinate transforms
//   8. check_interpolation_stencil — Validate stencil containment in box
// ============================================================================
// ============================================================================
// 1. valid_mirror — count fluid cells in the (iorder_t+1)^D interpolation block
//    whose lower corner is (i,j,k).  iorder_t=1: bilinear 2^D corners;
//    iorder_t=2: the 3^D WLS neighbourhood.
// ============================================================================
template <int iorder_t>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
int valid_mirror(
    int i, int j, int k,
    const Array4<const uint8_t>& ibMarkers)
{
    int fluid_count = 0;
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(k);
    for (int di = 0; di <= iorder_t; ++di) {
      for (int dj = 0; dj <= iorder_t; ++dj) {
          if (ibMarkers(i + di, j + dj, 0, 0) == 0) {
              ++fluid_count;
          }
      }
    }
    return fluid_count;
#else
    for (int di = 0; di <= iorder_t; ++di) {
      for (int dj = 0; dj <= iorder_t; ++dj) {
        for (int dk = 0; dk <= iorder_t; ++dk) {
          if (ibMarkers(i + di, j + dj, k + dk, 0) == 0) {
            ++fluid_count;
          }
        }
      }
    }
    return fluid_count;
#endif
}
// ============================================================================
// 1b. ibm_interp_base_index — lower-corner cell index of the interpolation
//     block for a point at physical coordinate x.
//     iorder_t=1: the 2-cell bracket (cell centres straddle x);
//     iorder_t=2: the 3-cell block centred on the cell containing x.
// ============================================================================
template <int iorder_t>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_interp_base_index(Real x, Real prob_lo_d, Real dx_d)
{
    if constexpr (iorder_t == 1) {
        return int(amrex::Math::floor((x - prob_lo_d) / dx_d - Real(0.5)));
    } else {
        return int(amrex::Math::floor((x - prob_lo_d) / dx_d)) - (iorder_t / 2);
    }
}

template <int NB>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_interp_matrix_resolvable(const Real matrix[NB][NB], int basis_size)
{
    Real factor[NB][NB] = {};
    Real maxdiag = matrix[0][0];
    for (int a = 1; a < basis_size; ++a) {
        maxdiag = amrex::max(maxdiag, matrix[a][a]);
    }
    const Real tolerance =
        Real(1.0e-11) * amrex::max(maxdiag, Real(1.0e-300));
    for (int a = 0; a < basis_size; ++a) {
        Real diagonal = matrix[a][a];
        for (int c = 0; c < a; ++c) {
            diagonal -= factor[a][c] * factor[a][c];
        }
        if (!(diagonal > tolerance) ||
            !amrex::Math::isfinite(diagonal)) {
            return false;
        }
        factor[a][a] = std::sqrt(diagonal);
        for (int b = a + 1; b < basis_size; ++b) {
            Real value = matrix[b][a];
            for (int c = 0; c < a; ++c) {
                value -= factor[b][c] * factor[a][c];
            }
            factor[b][a] = value / factor[a][a];
        }
    }
    return true;
}

// Return the highest polynomial order whose normal matrix is numerically
// resolvable for one candidate interpolation block.  The image-point search
// must use this information, not fluid-point count alone: at sharp 3-D STL
// features a stencil can contain many fluid cells that are nevertheless
// coplanar for the requested polynomial basis.
template <int iorder_t>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_interp_candidate_fit_order(
    const Array1D<int, 0, AMREX_SPACEDIM - 1>& base_ijk,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& image_xyz,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    const Array4<const uint8_t>& ibMarkers)
{
    if constexpr (iorder_t == 1) {
        amrex::ignore_unused(image_xyz, prob_lo, dxyz);
#if (AMREX_SPACEDIM == 2)
        return valid_mirror<iorder_t>(
                   base_ijk(0), base_ijk(1), 0, ibMarkers) > 0
                   ? 1
                   : 0;
#else
        return valid_mirror<iorder_t>(
                   base_ijk(0), base_ijk(1), base_ijk(2), ibMarkers) > 0
                   ? 1
                   : 0;
#endif
    } else {
        static_assert(iorder_t == 2,
                      "IBM interpolation order must be one or two");
        constexpr int NBQ = (AMREX_SPACEDIM == 2) ? 6 : 10;
        constexpr int NBL = AMREX_SPACEDIM + 1;
        constexpr int NIP = ipow(iorder_t + 1, AMREX_SPACEDIM);

        Real matrix[NBQ][NBQ] = {};
        int fluid_points = 0;
        for (int point = 0; point < NIP; ++point) {
            int rem = point;
            int ijk[AMREX_SPACEDIM];
            Real coordinate[AMREX_SPACEDIM];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                const int offset = rem % (iorder_t + 1);
                rem /= (iorder_t + 1);
                ijk[d] = base_ijk(d) + offset;
                const Real cell_center =
                    prob_lo[d] + (Real(ijk[d]) + Real(0.5)) * dxyz[d];
                coordinate[d] = (cell_center - image_xyz(d)) / dxyz[d];
            }
#if (AMREX_SPACEDIM == 3)
            const int kk = ijk[2];
#else
            const int kk = 0;
#endif
            if (ibMarkers(ijk[0], ijk[1], kk, 0) != 0) continue;
            ++fluid_points;

            Real basis[NBQ];
            basis[0] = Real(1.0);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                basis[1 + d] = coordinate[d];
            }
#if (AMREX_SPACEDIM == 2)
            basis[3] = coordinate[0] * coordinate[0];
            basis[4] = coordinate[0] * coordinate[1];
            basis[5] = coordinate[1] * coordinate[1];
#else
            basis[4] = coordinate[0] * coordinate[0];
            basis[5] = coordinate[1] * coordinate[1];
            basis[6] = coordinate[2] * coordinate[2];
            basis[7] = coordinate[0] * coordinate[1];
            basis[8] = coordinate[0] * coordinate[2];
            basis[9] = coordinate[1] * coordinate[2];
#endif
            for (int a = 0; a < NBQ; ++a) {
                for (int b = 0; b <= a; ++b) {
                    matrix[a][b] += basis[a] * basis[b];
                }
            }
        }
        for (int a = 0; a < NBQ; ++a) {
            for (int b = a + 1; b < NBQ; ++b) matrix[a][b] = matrix[b][a];
        }

        // Only the leading linear block is needed for the demoted fit because
        // the basis stores {1,x,y[,z]} first.
        if (fluid_points >= NBQ &&
            ibm_interp_matrix_resolvable<NBQ>(matrix, NBQ)) return 2;
        if (fluid_points >= NBL &&
            ibm_interp_matrix_resolvable<NBQ>(matrix, NBL)) return 1;
        return 0;
    }
}

// ============================================================================
// 2. search_optimal_image_point
// ============================================================================
template <int eorder_t, int iorder_t, typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
int search_optimal_image_point(
    const Point& cp_start,
    const LocalFrame& localframe,
    int lev,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_lev,
    Real di,
    const Box& bxg,
    const Array4<uint8_t const>& ibMarkers,
    const IPDATA& ipData,
    int f_idx,
    Array2D<Real, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
    Array2D< int, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
    Array1D<Real, 0, eorder_t - 1>& disIM,
    Array1D< int, 0, eorder_t - 1>& imp_ninterp)
{
    // Status bitmask returned to the caller so the GPU init kernel can
    // aggregate failures into device counters (device printf is lossy and
    // AMREX_ASSERT is stripped in release):
    //   bit 0: first-attempt stencil out of box (CPU aborts here)
    //   bit 1: best stencil below INTERP_THRESHOLD (placement failure)
    int status = 0;
    int best_fluid = 0;
    int best_fit_order = 0;
    Array1D<Real, 0, AMREX_SPACEDIM - 1> candi_xyz;
    Array1D< int, 0, AMREX_SPACEDIM - 1> candi_ijk;
    constexpr int N_ATTEMPTS = (GP_OR_SURF ? N_ATTEMPTS_GP : N_ATTEMPTS_SURF);
    // Image-point placement distance multipliers (multiples of di = alpha * cell_diagonal).
    // The first attempt places the image point at 1.0*di; if the stencil there is
    // insufficient (too many solid neighbours), we fall back to 1.5*di, 2.0*di, etc.
    // Local constexpr arrays are required for device-side address validity.
    constexpr Real IMP_FACTOR_GP_local[]   = {1.0, 1.5, 2.0};
    constexpr Real IMP_FACTOR_SURF_local[] = {1.0, 1.5, 2.0, 2.5, 3.0};
    const Real* IMP_FACTOR = (GP_OR_SURF ? IMP_FACTOR_GP_local : IMP_FACTOR_SURF_local);
    
    for (int attempt = 0; attempt < N_ATTEMPTS; ++attempt) {
    
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          candi_xyz(d) = cp_start[d] + IMP_FACTOR[attempt] * di * localframe.normal[d];
          candi_ijk(d) = ibm_interp_base_index<iorder_t>(candi_xyz(d), prob_lo[d], dx_lev[d]);
      }
    
#if (AMREX_SPACEDIM == 2)
      bool in_box = check_interpolation_stencil<IPDATA, iorder_t>(candi_ijk(0), candi_ijk(1), 0, 
                                                bxg, lev,
                                                ipData, f_idx,
                                                (attempt == 0) ? CheckMode::Abort : CheckMode::Silent);
      int n_fluid = (in_box) ? valid_mirror<iorder_t>(candi_ijk(0), candi_ijk(1), 0, ibMarkers) : -1;
#else
      bool in_box = check_interpolation_stencil<IPDATA, iorder_t>(candi_ijk(0), candi_ijk(1), candi_ijk(2),
                                                bxg, lev,
                                                ipData, f_idx,
                                                (attempt == 0) ? CheckMode::Abort : CheckMode::Silent);
      int n_fluid = (in_box) ? valid_mirror<iorder_t>(candi_ijk(0), candi_ijk(1), candi_ijk(2), ibMarkers) : -1;
#endif
      if (attempt == 0 && !in_box) { status |= 1; }
      int fit_order = 0;
      if (in_box) {
        fit_order = ibm_interp_candidate_fit_order<iorder_t>(
            candi_ijk, candi_xyz, prob_lo, dx_lev, ibMarkers);
      }
      if (fit_order > best_fit_order ||
          (fit_order == best_fit_order && n_fluid > best_fluid)) {
        best_fit_order = fit_order;
        best_fluid = n_fluid;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            imp_xyz(0, d) = candi_xyz(d);
            imp_ijk(0, d) = candi_ijk(d);
        }
        imp_ninterp(0) = n_fluid; 
        disIM(0) = IMP_FACTOR[attempt] * di;
      }
      constexpr int IDEAL_NINTERP = ipow(iorder_t + 1, AMREX_SPACEDIM);
      constexpr int IDEAL_FIT_ORDER = (iorder_t == 2) ? 2 : 1;
      if (best_fit_order == IDEAL_FIT_ORDER &&
          best_fluid == IDEAL_NINTERP) {
        break;
      }
    } // end loop on attempt
    constexpr int INTERP_THRESHOLD = (GP_OR_SURF ? INTERP_THRESHOLD_GP : INTERP_THRESHOLD_SURF);
    if (best_fluid < INTERP_THRESHOLD) {
        status |= 2;
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("Not enough valid interpolation points for first image point! "
                            "lev=%d f_idx=%d best_fluid=%d threshold=%d\n",
                            lev, f_idx, best_fluid, INTERP_THRESHOLD);
        AMREX_ASSERT(false);
#else
        int current_geom = -1;
        int current_elem = ipData.elemIdx[f_idx];
        const char* point_label;
        Real p_x = 0.0, p_y = 0.0, p_z = 0.0;
        if constexpr (GP_OR_SURF == 1) {
            current_geom = ipData.geomIdx[f_idx];
            point_label = "Ghost Point";
            
            int gp_i = ipData.gp_ijk[f_idx](0);
            int gp_j = ipData.gp_ijk[f_idx](1);
            p_x = prob_lo[0] + (0.5_rt + gp_i) * dx_lev[0];
            p_y = prob_lo[1] + (0.5_rt + gp_j) * dx_lev[1];
#if (AMREX_SPACEDIM == 3)
            int gp_k = ipData.gp_ijk[f_idx](2);
            p_z = prob_lo[2] + (0.5_rt + gp_k) * dx_lev[2];
#endif
        } else {
            point_label = "Surface Point";
            
            p_x = cp_start[0];
            p_y = cp_start[1];
#if (AMREX_SPACEDIM == 3)
            p_z = cp_start[2];
#endif
        }
        
#if (AMREX_SPACEDIM == 3)
        std::printf("Not enough valid interpolation points found for the first image point!\n"
                    "  Level: %d\n"
                    "  %s: (%f, %f, %f)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Best Fluid Points Found: %d (Threshold: %d)\n",
                    lev, point_label, p_x, p_y, p_z, 
                    current_geom, current_elem,
                    best_fluid, INTERP_THRESHOLD);
#else
        std::printf("Not enough valid interpolation points found for the first image point!\n"
                    "  Level: %d\n"
                    "  %s: (%f, %f)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Best Fluid Points Found: %d (Threshold: %d)\n",
                    lev, point_label, p_x, p_y,
                    current_geom, current_elem,
                    best_fluid, INTERP_THRESHOLD);
#endif
        std::fflush(stdout);
#endif // AMREX_DEVICE_COMPILE
    }
    return status;
}
// ============================================================================
// 3. search_image_point
// ============================================================================
template <int order_t, int iorder_t, typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void search_image_point(
    int jj,
    const Point& cp_start,
    const LocalFrame& localframe,
    int lev,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_lev,
    Real di,
    const Box& bxg,
    const Array4<uint8_t const>& ibMarkers,
    const IPDATA& ipData,
    int f_idx,
    Array2D<Real, 0, order_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
    Array2D< int, 0, order_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
    Array1D<Real, 0, order_t - 1>& disIM,
    Array1D< int, 0, order_t - 1>& imp_ninterp)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        imp_xyz(jj, d) = cp_start[d] + di * localframe.normal[d];
        imp_ijk(jj, d) = ibm_interp_base_index<iorder_t>(imp_xyz(jj, d), prob_lo[d], dx_lev[d]);
    }
#if (AMREX_SPACEDIM == 2)
    bool in_box = check_interpolation_stencil<IPDATA, iorder_t>(imp_ijk(jj, 0), imp_ijk(jj, 1), 0,
                                              bxg, lev,
                                              ipData, f_idx,
                                              CheckMode::Silent);
    int fluid = (in_box) ? valid_mirror<iorder_t>(imp_ijk(jj, 0), imp_ijk(jj, 1), 0, ibMarkers) : -1;
#else
    bool in_box = check_interpolation_stencil<IPDATA, iorder_t>(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2),
                                              bxg, lev,
                                              ipData, f_idx,
                                              CheckMode::Silent);
    int fluid = (in_box) ? valid_mirror<iorder_t>(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2), ibMarkers) : -1;
#endif
    disIM(jj) = (jj > 0) ? disIM(jj - 1) + di : di;
    imp_ninterp(jj) = fluid;
}
// ============================================================================
// 3b. place_image_points — chained walk along the surface normal
//
// Drives the eorder_t image-point sequence:
//   IP_0 ← search_optimal_image_point starting from the IB point on surface
//   IP_jj ← search_image_point starting from IP_{jj-1}, fixed step di
//
// Used identically by both the GPU and CPU paths in initialiseGPs.  The
// `cp_start` slide between iterations is the only place the two backend
// Point types diverge (CGAL Point uses ctor, BVH Point uses subscript) —
// kept in this helper so the call site stays single-line.
// ============================================================================
template <int eorder_t, int iorder_t, typename IPDATA>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
int place_image_points(
    const Point& cp,
    const LocalFrame& localframe,
    int lev,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_lev,
    Real di,
    const Box& bxg,
    const Array4<uint8_t const>& ibMarkers,
    const IPDATA& ipData,
    int gidx,
    Array2D<Real, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
    Array2D< int, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
    Array1D<Real, 0, eorder_t - 1>& disIM,
    Array1D< int, 0, eorder_t - 1>& imp_ninterp)
{
    Point cp_start = cp;
    int status = 0;
    for (int jj = 0; jj < eorder_t; ++jj) {
        if (jj == 0) {
            status = search_optimal_image_point<eorder_t, iorder_t>(
                cp_start, localframe,
                lev, prob_lo, dx_lev, di,
                bxg, ibMarkers,
                ipData, gidx,
                imp_xyz, imp_ijk, disIM, imp_ninterp);
        } else {
            search_image_point<eorder_t, iorder_t>(
                jj, cp_start, localframe,
                lev, prob_lo, dx_lev, di,
                bxg, ibMarkers,
                ipData, gidx,
                imp_xyz, imp_ijk, disIM, imp_ninterp);
        }
        // Slide cp_start to the IP just placed, ready for jj+1.
#if defined(AMREX_USE_CGAL)
#if (AMREX_SPACEDIM == 2)
        cp_start = Point(imp_xyz(jj, 0), imp_xyz(jj, 1));
#else
        cp_start = Point(imp_xyz(jj, 0), imp_xyz(jj, 1), imp_xyz(jj, 2));
#endif
#else
        for (int d = 0; d < AMREX_SPACEDIM; ++d) cp_start[d] = imp_xyz(jj, d);
#endif
    }
    return status;
}
// ============================================================================
// 3c. ibm_spd_solve_e0 — solve M y = e0 for a small SPD system via Cholesky
//     with a relative pivot guard.  Used by the iorder=2 WLS interpolation:
//     with the polynomial basis centred at the image point, the interpolated
//     value is the constant coefficient, whose normal-equation solution only
//     needs M^{-1} e0.  Returns false when the (masked) stencil is too
//     ill-conditioned for this basis — the caller then demotes the basis.
// ============================================================================
template <int NB>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_spd_solve_e0(const Real M[NB][NB], Real y[NB])
{
    Real maxdiag = M[0][0];
    for (int a = 1; a < NB; ++a) maxdiag = amrex::max(maxdiag, M[a][a]);
    const Real tol = Real(1.0e-11) * amrex::max(maxdiag, Real(1.0e-300));

    Real L[NB][NB] = {};
    for (int a = 0; a < NB; ++a) {
        Real sdiag = M[a][a];
        for (int c = 0; c < a; ++c) sdiag -= L[a][c] * L[a][c];
        if (!(sdiag > tol)) return false;
        L[a][a] = std::sqrt(sdiag);
        for (int b = a + 1; b < NB; ++b) {
            Real t = M[b][a];
            for (int c = 0; c < a; ++c) t -= L[b][c] * L[a][c];
            L[b][a] = t / L[a][a];
        }
    }
    Real z[NB];
    for (int a = 0; a < NB; ++a) {
        Real t = (a == 0) ? Real(1.0) : Real(0.0);
        for (int c = 0; c < a; ++c) t -= L[a][c] * z[c];
        z[a] = t / L[a][a];
    }
    for (int a = NB - 1; a >= 0; --a) {
        Real t = z[a];
        for (int c = a + 1; c < NB; ++c) t -= L[c][a] * y[c];
        y[a] = t / L[a][a];
    }
    return true;
}

template <int NB>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_spd_factor(const Real M[NB][NB], Real L[NB][NB])
{
    Real maxdiag = M[0][0];
    for (int a = 1; a < NB; ++a) maxdiag = amrex::max(maxdiag, M[a][a]);
    const Real tol = Real(1.0e-11) * amrex::max(maxdiag, Real(1.0e-300));
    for (int a = 0; a < NB; ++a) {
        Real sdiag = M[a][a];
        for (int c = 0; c < a; ++c) sdiag -= L[a][c] * L[a][c];
        if (!(sdiag > tol) || !amrex::Math::isfinite(sdiag)) return false;
        L[a][a] = std::sqrt(sdiag);
        for (int b = a + 1; b < NB; ++b) {
            Real value = M[b][a];
            for (int c = 0; c < a; ++c) value -= L[b][c] * L[a][c];
            L[b][a] = value / L[a][a];
        }
    }
    return true;
}

template <int NB>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ibm_spd_solve_factored(
    const Real L[NB][NB], const Real rhs[NB], Real solution[NB])
{
    Real work[NB];
    for (int a = 0; a < NB; ++a) {
        Real value = rhs[a];
        for (int c = 0; c < a; ++c) value -= L[a][c] * work[c];
        work[a] = value / L[a][a];
    }
    for (int a = NB - 1; a >= 0; --a) {
        Real value = work[a];
        for (int c = a + 1; c < NB; ++c) value -= L[c][a] * solution[c];
        solution[a] = value / L[a][a];
    }
}

template <int NB>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_matrix_numerical_rank(const Real matrix[NB][NB])
{
    Real work[NB][NB];
    Real scale = Real(0.0);
    for (int row = 0; row < NB; ++row) {
        for (int col = 0; col < NB; ++col) {
            work[row][col] = matrix[row][col];
            scale = amrex::max(scale, amrex::Math::abs(matrix[row][col]));
        }
    }
    if (!(scale > Real(0.0)) || !amrex::Math::isfinite(scale)) return 0;
    const Real tolerance = Real(1.0e-11) * scale;

    int rank = 0;
    for (int col = 0; col < NB && rank < NB; ++col) {
        int pivot = rank;
        Real pivot_abs = amrex::Math::abs(work[pivot][col]);
        for (int row = rank + 1; row < NB; ++row) {
            const Real candidate = amrex::Math::abs(work[row][col]);
            if (candidate > pivot_abs) {
                pivot = row;
                pivot_abs = candidate;
            }
        }
        if (!(pivot_abs > tolerance)) continue;
        if (pivot != rank) {
            for (int c = col; c < NB; ++c) {
                const Real temporary = work[rank][c];
                work[rank][c] = work[pivot][c];
                work[pivot][c] = temporary;
            }
        }
        const Real inverse = Real(1.0) / work[rank][col];
        for (int row = rank + 1; row < NB; ++row) {
            const Real factor = work[row][col] * inverse;
            for (int c = col; c < NB; ++c) {
                work[row][c] -= factor * work[rank][c];
            }
        }
        ++rank;
    }
    return rank;
}

/// Exact 1-norm condition number of a small SPD matrix, evaluated from its
/// Cholesky factor. Returns +inf when the same guarded factorisation used by
/// the production WLS rejects the matrix.
template <int NB>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ibm_spd_condition_1(const Real matrix[NB][NB])
{
    Real factor[NB][NB] = {};
    if (!ibm_spd_factor<NB>(matrix, factor)) {
        return std::numeric_limits<Real>::infinity();
    }

    Real matrix_norm = Real(0.0);
    for (int col = 0; col < NB; ++col) {
        Real column_sum = Real(0.0);
        for (int row = 0; row < NB; ++row) {
            column_sum += amrex::Math::abs(matrix[row][col]);
        }
        matrix_norm = amrex::max(matrix_norm, column_sum);
    }

    Real inverse_norm = Real(0.0);
    for (int col = 0; col < NB; ++col) {
        Real rhs[NB] = {};
        Real solution[NB] = {};
        rhs[col] = Real(1.0);
        ibm_spd_solve_factored<NB>(factor, rhs, solution);
        Real column_sum = Real(0.0);
        for (int row = 0; row < NB; ++row) {
            column_sum += amrex::Math::abs(solution[row]);
        }
        inverse_norm = amrex::max(inverse_norm, column_sum);
    }
    const Real condition = matrix_norm * inverse_norm;
    return amrex::Math::isfinite(condition)
               ? condition
               : std::numeric_limits<Real>::infinity();
}

// Cartesian static no-slip compatibility from an explicitly selected,
// visibility-filtered fluid support:
//
//   dp/dn = n_i d_j tau_ij,
//   tau_ij = mu (u_i,j + u_j,i - 2/3 theta delta_ij)
//            + xi theta delta_ij.
//
// Velocity, mu(T), and xi(T) are fitted with a quadratic WLS basis about the
// requested evaluation point. Returning valid=0 requests the wall closure's
// homogeneous-Neumann fallback.
template <int NIP, typename IjkArr, typename MaskArr, typename XYZArr,
          typename PrimArr, typename NormArr, typename Closure>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ibm_pressure_compatibility_t
ibm_viscous_pressure_compatibility_from_support(
    const IjkArr& support_ijk,
    const MaskArr& support_mask,
    const XYZArr& evaluation_xyz,
    const PrimArr& prims,
    const NormArr& normal,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    const Closure* cls,
    bool cartesian_geometry,
    Real maximum_condition)
{
    ibm_pressure_compatibility_t result{};
#ifdef USE_PELEPHYSICS
    // The lightweight closure object exposes only placeholder transport values
    // in PelePhysics builds; the true mixture coefficients live in viscous.h.
    // Using those placeholders here would be physically inconsistent.
    amrex::ignore_unused(support_ijk, support_mask, evaluation_xyz, prims,
                         normal, prob_lo, dxyz, cls, cartesian_geometry,
                         maximum_condition);
    return result;
#else
    static_assert(AMREX_SPACEDIM == 2 || AMREX_SPACEDIM == 3,
                  "viscous pressure compatibility supports 2-D/3-D Cartesian grids");
    constexpr int NBQ = (AMREX_SPACEDIM == 2) ? 6 : 10;
    if (!cartesian_geometry) return result;

    Real matrix[NBQ][NBQ] = {};
    Real velocity_rhs[AMREX_SPACEDIM][NBQ] = {};
    Real mu_rhs[NBQ] = {};
    Real xi_rhs[NBQ] = {};
    int support_count = 0;

    for (int point = 0; point < NIP; ++point) {
        if (support_mask(point) == 0) continue;
        ++support_count;
        int ijk[AMREX_SPACEDIM];
        Real coordinate[AMREX_SPACEDIM];
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ijk[d] = support_ijk(point, d);
            const Real cell_center =
                prob_lo[d] + (Real(ijk[d]) + Real(0.5)) * dxyz[d];
            coordinate[d] =
                (cell_center - evaluation_xyz(d)) / dxyz[d];
        }
#if (AMREX_SPACEDIM == 3)
        const int kk = ijk[2];
#else
        const int kk = 0;
#endif
        const int ii = ijk[0];
        const int jj = ijk[1];

        Real basis[NBQ];
        basis[0] = Real(1.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) basis[1 + d] = coordinate[d];
#if (AMREX_SPACEDIM == 2)
        basis[3] = coordinate[0] * coordinate[0];
        basis[4] = coordinate[0] * coordinate[1];
        basis[5] = coordinate[1] * coordinate[1];
#else
        basis[4] = coordinate[0] * coordinate[0];
        basis[5] = coordinate[1] * coordinate[1];
        basis[6] = coordinate[2] * coordinate[2];
        basis[7] = coordinate[0] * coordinate[1];
        basis[8] = coordinate[0] * coordinate[2];
        basis[9] = coordinate[1] * coordinate[2];
#endif

        const Real temperature = prims(ii, jj, kk, cls_t::QT);
        if (!(temperature > Real(0.0)) || !amrex::Math::isfinite(temperature)) {
            return result;
        }
        Real transport_temperature = temperature;
        const Real mu_value = cls->visc(transport_temperature);
        const Real xi_value = cls->xi(transport_temperature);
        if (!(mu_value >= Real(0.0)) || !amrex::Math::isfinite(mu_value) ||
            !amrex::Math::isfinite(xi_value)) {
            return result;
        }

        for (int a = 0; a < NBQ; ++a) {
            for (int b = 0; b <= a; ++b) {
                matrix[a][b] += basis[a] * basis[b];
            }
            for (int component = 0; component < AMREX_SPACEDIM; ++component) {
                const Real velocity = prims(ii, jj, kk, cls_t::QU + component);
                if (!amrex::Math::isfinite(velocity)) return result;
                velocity_rhs[component][a] += basis[a] * velocity;
            }
            mu_rhs[a] += basis[a] * mu_value;
            xi_rhs[a] += basis[a] * xi_value;
        }
    }
    if (support_count < NBQ) return result;
    for (int a = 0; a < NBQ; ++a) {
        for (int b = a + 1; b < NBQ; ++b) matrix[a][b] = matrix[b][a];
    }

    if (amrex::Math::isfinite(maximum_condition)) {
        const Real condition = ibm_spd_condition_1<NBQ>(matrix);
        if (!amrex::Math::isfinite(condition) ||
            condition > maximum_condition) {
            return result;
        }
    }

    Real factor[NBQ][NBQ] = {};
    if (!ibm_spd_factor<NBQ>(matrix, factor)) return result;

    Real velocity_coef[AMREX_SPACEDIM][NBQ] = {};
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        ibm_spd_solve_factored<NBQ>(
            factor, velocity_rhs[component], velocity_coef[component]);
    }
    Real mu_coef[NBQ] = {};
    Real xi_coef[NBQ] = {};
    ibm_spd_solve_factored<NBQ>(factor, mu_rhs, mu_coef);
    ibm_spd_solve_factored<NBQ>(factor, xi_rhs, xi_coef);

    Real grad_velocity[AMREX_SPACEDIM][AMREX_SPACEDIM] = {};
    Real hess_velocity[AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM] = {};
    Real grad_mu[AMREX_SPACEDIM] = {};
    Real grad_xi[AMREX_SPACEDIM] = {};
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            grad_velocity[component][d] =
                velocity_coef[component][1 + d] / dxyz[d];
        }
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        grad_mu[d] = mu_coef[1 + d] / dxyz[d];
        grad_xi[d] = xi_coef[1 + d] / dxyz[d];
    }

#if (AMREX_SPACEDIM == 2)
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        hess_velocity[component][0][0] =
            Real(2.0) * velocity_coef[component][3] / (dxyz[0] * dxyz[0]);
        hess_velocity[component][0][1] = hess_velocity[component][1][0] =
            velocity_coef[component][4] / (dxyz[0] * dxyz[1]);
        hess_velocity[component][1][1] =
            Real(2.0) * velocity_coef[component][5] / (dxyz[1] * dxyz[1]);
    }
#else
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        hess_velocity[component][0][0] =
            Real(2.0) * velocity_coef[component][4] / (dxyz[0] * dxyz[0]);
        hess_velocity[component][1][1] =
            Real(2.0) * velocity_coef[component][5] / (dxyz[1] * dxyz[1]);
        hess_velocity[component][2][2] =
            Real(2.0) * velocity_coef[component][6] / (dxyz[2] * dxyz[2]);
        hess_velocity[component][0][1] = hess_velocity[component][1][0] =
            velocity_coef[component][7] / (dxyz[0] * dxyz[1]);
        hess_velocity[component][0][2] = hess_velocity[component][2][0] =
            velocity_coef[component][8] / (dxyz[0] * dxyz[2]);
        hess_velocity[component][1][2] = hess_velocity[component][2][1] =
            velocity_coef[component][9] / (dxyz[1] * dxyz[2]);
    }
#endif

    const Real mu = mu_coef[0];
    const Real xi = xi_coef[0];
    if (!(mu >= Real(0.0)) || !amrex::Math::isfinite(mu) ||
        !amrex::Math::isfinite(xi)) {
        return result;
    }

    Real theta = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        theta += grad_velocity[d][d];
    }

    Real div_tau[AMREX_SPACEDIM] = {};
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        Real laplacian = Real(0.0);
        Real grad_theta = Real(0.0);
        Real variable_mu = Real(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            laplacian += hess_velocity[component][d][d];
            grad_theta += hess_velocity[d][d][component];
            variable_mu += grad_mu[d] *
                (grad_velocity[component][d] + grad_velocity[d][component]);
        }
        div_tau[component] =
            mu * laplacian + variable_mu +
            (xi + mu / Real(3.0)) * grad_theta +
            (grad_xi[component] - Real(2.0 / 3.0) * grad_mu[component]) * theta;
        if (!amrex::Math::isfinite(div_tau[component])) return result;
    }

    Real dpdn = Real(0.0);
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        dpdn += normal(component) * div_tau[component];
    }
    if (!amrex::Math::isfinite(dpdn)) return result;
    result.dpdn = dpdn;
    result.valid = 1;
    return result;
#endif
}

// Historical first-image-point adapter. A complete 3^D fluid block preserves
// the established pressure-compatibility path. The BI-CWLS caller may
// separately retry the same fit on its accepted visible support when this
// complete block is unavailable.
template <int iorder_t, int NIP, typename IjkArr, typename XYZArr,
          typename PrimArr, typename NormArr, typename Closure>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ibm_pressure_compatibility_t ibm_viscous_pressure_compatibility(
    const IjkArr& imp_ip_ijk,
    const XYZArr& imp_xyz,
    int first_ip_fluid_points,
    const PrimArr& prims,
    const NormArr& normal,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    const Closure* cls,
    bool cartesian_geometry)
{
    static_assert(iorder_t == 2,
                  "viscous pressure compatibility requires quadratic WLS");
    constexpr int FULL_STENCIL = ipow(3, AMREX_SPACEDIM);
    static_assert(NIP == FULL_STENCIL,
                  "quadratic IBM WLS must use the complete 3^D candidate block");

    ibm_pressure_compatibility_t result{};
    if (!cartesian_geometry || first_ip_fluid_points != FULL_STENCIL) {
        return result;
    }

    Array2D<int, 0, NIP - 1, 0, AMREX_SPACEDIM - 1> support_ijk{};
    Array1D<uint8_t, 0, NIP - 1> support_mask{};
    Array1D<Real, 0, AMREX_SPACEDIM - 1> evaluation_xyz{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        evaluation_xyz(d) = imp_xyz(0, d);
    }
    for (int point = 0; point < NIP; ++point) {
        support_mask(point) = uint8_t(1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            support_ijk(point, d) = imp_ip_ijk(0, point, d);
        }
    }
    return ibm_viscous_pressure_compatibility_from_support<NIP>(
        support_ijk, support_mask, evaluation_xyz, prims, normal, prob_lo,
        dxyz, cls, cartesian_geometry,
        std::numeric_limits<Real>::infinity());
}

// ============================================================================
// BI-constrained shared-GP functionals
// ============================================================================

template <int NB, int NMAX, int NIP, typename MaskArray,
          typename WeightArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_fit_target_functional(
    const Real (&support_basis)[NIP][NMAX],
    const MaskArray& support_mask, const Real (&target_basis)[NMAX],
    Real maximum_condition, WeightArray& weights, Real& condition)
{
    Real matrix[NB][NB] = {};
    int active_supports = 0;
    for (int support = 0; support < NIP; ++support) {
        if (support_mask(support) == 0) continue;
        ++active_supports;
        for (int row = 0; row < NB; ++row) {
            for (int col = 0; col <= row; ++col) {
                matrix[row][col] += support_basis[support][row] *
                                    support_basis[support][col];
            }
        }
    }
    if (active_supports < NB) return false;
    for (int row = 0; row < NB; ++row) {
        for (int col = row + 1; col < NB; ++col) {
            matrix[row][col] = matrix[col][row];
        }
    }
    if (ibm_matrix_numerical_rank<NB>(matrix) != NB) return false;
    condition = ibm_spd_condition_1<NB>(matrix);
    if (!amrex::Math::isfinite(condition) ||
        condition > maximum_condition) {
        return false;
    }

    Real factor[NB][NB] = {};
    if (!ibm_spd_factor<NB>(matrix, factor)) return false;
    Real rhs[NB];
    Real solution[NB] = {};
    for (int row = 0; row < NB; ++row) rhs[row] = target_basis[row];
    ibm_spd_solve_factored<NB>(factor, rhs, solution);

    for (int support = 0; support < NIP; ++support) {
        Real value = Real(0.0);
        if (support_mask(support) != 0) {
            for (int row = 0; row < NB; ++row) {
                value += support_basis[support][row] * solution[row];
            }
        }
        weights(support) = value;
    }
    return true;
}

/// Build a one-sided fluid trace at the boundary intercept from the same
/// visible support used by the BI-constrained ghost reconstruction.  Unlike
/// the Dirichlet/Neumann ghost functionals below, this fit imposes no wall
/// condition: it approximates the fluid-side limit required by compatibility
/// relations whose data depend on the current Runge--Kutta stage.
///
/// The quadratic basis is complete in the local wall frame.  A linear
/// fallback is retained for robustness and reported separately by the caller;
/// smooth-wall verification requires the quadratic path on every GP.
template <int NIP, typename IndexArray, typename MaskArray,
          typename PointArray, typename VectorArray, typename WeightArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_unconstrained_wall_trace_functional(
    const IndexArray& support_indices, const MaskArray& support_mask,
    const PointArray& wall_point, const VectorArray& wall_normal,
    const VectorArray& tangent1, const VectorArray& tangent2,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    Real length_scale, Real maximum_condition, WeightArray& trace_weights,
    Real& accepted_condition, Real& weight_l1)
{
    constexpr int NBL = AMREX_SPACEDIM + 1;
    constexpr int NBQ = (AMREX_SPACEDIM == 2) ? 6 : 10;

    for (int support = 0; support < NIP; ++support) {
        trace_weights(support) = Real(0.0);
    }
    accepted_condition = std::numeric_limits<Real>::infinity();
    weight_l1 = std::numeric_limits<Real>::infinity();
    if (!(length_scale > Real(0.0)) ||
        !amrex::Math::isfinite(length_scale)) {
        return -1;
    }

    Real normal[AMREX_SPACEDIM];
    Real tangent_a[AMREX_SPACEDIM];
    Real tangent_b[AMREX_SPACEDIM];
    Real normal_norm2 = Real(0.0);
    Real tangent_a_norm2 = Real(0.0);
    Real tangent_b_norm2 = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        normal[d] = wall_normal(d);
        tangent_a[d] = tangent1(d);
        tangent_b[d] = tangent2(d);
        normal_norm2 += normal[d] * normal[d];
        tangent_a_norm2 += tangent_a[d] * tangent_a[d];
        tangent_b_norm2 += tangent_b[d] * tangent_b[d];
    }
    if (!(normal_norm2 > Real(0.0)) ||
        !(tangent_a_norm2 > Real(0.0))) {
        return -1;
    }
    const Real inverse_normal = Real(1.0) / std::sqrt(normal_norm2);
    const Real inverse_tangent_a = Real(1.0) / std::sqrt(tangent_a_norm2);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        normal[d] *= inverse_normal;
        tangent_a[d] *= inverse_tangent_a;
    }
#if (AMREX_SPACEDIM == 3)
    if (!(tangent_b_norm2 > Real(0.0))) return -1;
    const Real inverse_tangent_b = Real(1.0) / std::sqrt(tangent_b_norm2);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        tangent_b[d] *= inverse_tangent_b;
    }
#else
    amrex::ignore_unused(tangent_b_norm2);
#endif

    Real support_basis[NIP][NBQ] = {};
    for (int support = 0; support < NIP; ++support) {
        if (support_mask(support) == 0) continue;
        Real s = Real(0.0);
        Real t = Real(0.0);
        Real r = Real(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int index = support_indices(support, d);
            const Real centre =
                prob_lo[d] + (Real(index) + Real(0.5)) * dxyz[d];
            const Real delta = (centre - wall_point(d)) / length_scale;
            s += delta * normal[d];
            t += delta * tangent_a[d];
#if (AMREX_SPACEDIM == 3)
            r += delta * tangent_b[d];
#endif
        }
#if (AMREX_SPACEDIM == 2)
        support_basis[support][0] = Real(1.0);
        support_basis[support][1] = s;
        support_basis[support][2] = t;
        support_basis[support][3] = s * s;
        support_basis[support][4] = s * t;
        support_basis[support][5] = t * t;
#else
        support_basis[support][0] = Real(1.0);
        support_basis[support][1] = s;
        support_basis[support][2] = t;
        support_basis[support][3] = r;
        support_basis[support][4] = s * s;
        support_basis[support][5] = s * t;
        support_basis[support][6] = s * r;
        support_basis[support][7] = t * t;
        support_basis[support][8] = t * r;
        support_basis[support][9] = r * r;
#endif
    }

    Real wall_target[NBQ] = {};
    wall_target[0] = Real(1.0);
    Array1D<Real, 0, NIP - 1> trial_weights{};
    Real condition = std::numeric_limits<Real>::infinity();
    int accepted_order = -1;
    if (ibm_fit_target_functional<NBQ, NBQ, NIP>(
            support_basis, support_mask, wall_target, maximum_condition,
            trial_weights, condition)) {
        accepted_order = 2;
    } else {
        condition = std::numeric_limits<Real>::infinity();
        if (!ibm_fit_target_functional<NBL, NBQ, NIP>(
                support_basis, support_mask, wall_target, maximum_condition,
                trial_weights, condition)) {
            return -1;
        }
        accepted_order = 1;
    }

    Real l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        trace_weights(support) = trial_weights(support);
        l1 += amrex::Math::abs(trial_weights(support));
    }
    if (!amrex::Math::isfinite(condition + l1)) return -1;
    accepted_condition = condition;
    weight_l1 = l1;
    return accepted_order;
}

// Exact annular target functional for the local 2-D BI frame. The wall
// constraint remains pointwise; only the stored solid-side ghost-cell target
// uses the cylindrical r dr dz measure. This helper is compiled out of every
// Cartesian production configuration.
template <typename PointArray, typename VectorArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ibm_rz_annular_ghost_target_moments(
    const PointArray& wall_point, const PointArray& ghost_point,
    const VectorArray& wall_normal, const VectorArray& tangent1,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz, Real length_scale,
    GpuArray<Real, cerisse::ibm::rz_quadratic_moment_count>& moments)
{
#if (AMREX_SPACEDIM == 2)
    const Real roundoff = Real(64.0) *
        std::numeric_limits<Real>::epsilon() *
        amrex::max(Real(1.0), amrex::Math::abs(ghost_point(0)));
    Real r_lo = ghost_point(0) - Real(0.5) * dxyz[0];
    if (amrex::Math::abs(r_lo) <= roundoff) r_lo = Real(0.0);
    const cerisse::ibm::RZCellBounds target{
        r_lo, ghost_point(0) + Real(0.5) * dxyz[0],
        ghost_point(1) - Real(0.5) * dxyz[1],
        ghost_point(1) + Real(0.5) * dxyz[1]};
    const GpuArray<Real, 2> wall{{wall_point(0), wall_point(1)}};
    const GpuArray<Real, 2> normal{{wall_normal(0), wall_normal(1)}};
    const GpuArray<Real, 2> tangent{{tangent1(0), tangent1(1)}};
    return cerisse::ibm::rz_annular_local_quadratic_moments(
        target, wall, normal, tangent, length_scale, moments);
#else
    amrex::ignore_unused(
        wall_point, ghost_point, wall_normal, tangent1, dxyz, length_scale,
        moments);
    return false;
#endif
}

/// Reconstruct a scalar Cartesian ghost-cell average from Cartesian support-
/// cell averages while imposing a pointwise Dirichlet value at the boundary
/// intercept.  The support and target rows contain exact cell averages of the
/// local quadratic monomials.  This is the finite-volume counterpart of the
/// point-evaluation Dirichlet functional used by
/// ibm_bi_constrained_gp_functionals().
///
/// For a stationary slip wall this functional is applied to the conserved
/// normal momentum, whose boundary value is exactly zero.  Reconstructing
/// \f$\overline{\rho u_n}\f$ avoids treating Q(Ubar) as a point sample of
/// primitive normal velocity.
template <int NIP,
          bool AnnularTarget = rz_annular_bic_cell_average,
          typename IndexArray, typename MaskArray,
          typename PointArray, typename VectorArray, typename WeightArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_bi_cell_average_dirichlet_functional(
    const IndexArray& support_indices, const MaskArray& support_mask,
    const PointArray& wall_point, const PointArray& ghost_point,
    const VectorArray& wall_normal, const VectorArray& tangent1,
    const VectorArray& tangent2,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    Real length_scale, Real maximum_condition, WeightArray& weights,
    Real& boundary_weight, Real& accepted_condition, Real& weight_l1)
{
    constexpr int NL = AMREX_SPACEDIM;
    constexpr int NQ = (AMREX_SPACEDIM == 2) ? 5 : 9;

    for (int support = 0; support < NIP; ++support) {
        weights(support) = Real(0.0);
    }
    boundary_weight = Real(0.0);
    accepted_condition = std::numeric_limits<Real>::infinity();
    weight_l1 = std::numeric_limits<Real>::infinity();
    if (!(length_scale > Real(0.0)) ||
        !amrex::Math::isfinite(length_scale)) {
        return -1;
    }

    Real frame[AMREX_SPACEDIM][AMREX_SPACEDIM] = {};
    Real frame_norm2[AMREX_SPACEDIM] = {};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        frame[0][d] = wall_normal(d);
        frame[1][d] = tangent1(d);
        frame_norm2[0] += frame[0][d] * frame[0][d];
        frame_norm2[1] += frame[1][d] * frame[1][d];
#if (AMREX_SPACEDIM == 3)
        frame[2][d] = tangent2(d);
        frame_norm2[2] += frame[2][d] * frame[2][d];
#endif
    }
    for (int a = 0; a < AMREX_SPACEDIM; ++a) {
        if (!(frame_norm2[a] > Real(0.0))) return -1;
        const Real inverse_norm = Real(1.0) / std::sqrt(frame_norm2[a]);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            frame[a][d] *= inverse_norm;
        }
    }

    // Cartesian mode consumes Cartesian support-cell averages. R-Z mode
    // consumes recovered centre point states, so its support covariance is
    // zero and only the target carries the exact annular moments.
    Real covariance[AMREX_SPACEDIM][AMREX_SPACEDIM] = {};
    if constexpr (!rz_annular_bic_cell_average) {
        const Real inverse_scale2 = Real(1.0) / (length_scale * length_scale);
        for (int a = 0; a < AMREX_SPACEDIM; ++a) {
            for (int b = 0; b < AMREX_SPACEDIM; ++b) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    covariance[a][b] +=
                        dxyz[d] * dxyz[d] * frame[a][d] * frame[b][d] *
                        inverse_scale2 / Real(12.0);
                }
            }
        }
    }

    Real support_basis[NIP][NQ] = {};
    for (int support = 0; support < NIP; ++support) {
        if (support_mask(support) == 0) continue;
        Real coordinate[AMREX_SPACEDIM] = {};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int index = support_indices(support, d);
            const Real centre =
                prob_lo[d] + (Real(index) + Real(0.5)) * dxyz[d];
            const Real delta = (centre - wall_point(d)) / length_scale;
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                coordinate[a] += delta * frame[a][d];
            }
        }
#if (AMREX_SPACEDIM == 2)
        const Real s = coordinate[0];
        const Real t = coordinate[1];
        support_basis[support][0] = s;
        support_basis[support][1] = t;
        support_basis[support][2] = s * s + covariance[0][0];
        support_basis[support][3] = s * t + covariance[0][1];
        support_basis[support][4] = t * t + covariance[1][1];
#else
        const Real s = coordinate[0];
        const Real t = coordinate[1];
        const Real r = coordinate[2];
        support_basis[support][0] = s;
        support_basis[support][1] = t;
        support_basis[support][2] = r;
        support_basis[support][3] = s * s + covariance[0][0];
        support_basis[support][4] = s * t + covariance[0][1];
        support_basis[support][5] = s * r + covariance[0][2];
        support_basis[support][6] = t * t + covariance[1][1];
        support_basis[support][7] = t * r + covariance[1][2];
        support_basis[support][8] = r * r + covariance[2][2];
#endif
    }

    Real ghost_coordinate[AMREX_SPACEDIM] = {};
    if constexpr (!AnnularTarget) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const Real delta = (ghost_point(d) - wall_point(d)) / length_scale;
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                ghost_coordinate[a] += delta * frame[a][d];
            }
        }
    }
    Real target_basis[NQ] = {};
#if (AMREX_SPACEDIM == 2)
    if constexpr (AnnularTarget) {
        GpuArray<Real, cerisse::ibm::rz_quadratic_moment_count> target{};
        if (!ibm_rz_annular_ghost_target_moments(
                wall_point, ghost_point, wall_normal, tangent1, dxyz,
                length_scale, target)) {
            return -1;
        }
        target_basis[0] = target[cerisse::ibm::rz_z];
        target_basis[1] = target[cerisse::ibm::rz_r];
        target_basis[2] = target[cerisse::ibm::rz_z2];
        target_basis[3] = target[cerisse::ibm::rz_zr];
        target_basis[4] = target[cerisse::ibm::rz_r2];
    } else {
        const Real ghost_s = ghost_coordinate[0];
        const Real ghost_t = ghost_coordinate[1];
        target_basis[0] = ghost_s;
        target_basis[1] = ghost_t;
        target_basis[2] = ghost_s * ghost_s + covariance[0][0];
        target_basis[3] = ghost_s * ghost_t + covariance[0][1];
        target_basis[4] = ghost_t * ghost_t + covariance[1][1];
    }
#else
    const Real ghost_s = ghost_coordinate[0];
    const Real ghost_t = ghost_coordinate[1];
    const Real ghost_r = ghost_coordinate[2];
    target_basis[0] = ghost_s;
    target_basis[1] = ghost_t;
    target_basis[2] = ghost_r;
    target_basis[3] = ghost_s * ghost_s + covariance[0][0];
    target_basis[4] = ghost_s * ghost_t + covariance[0][1];
    target_basis[5] = ghost_s * ghost_r + covariance[0][2];
    target_basis[6] = ghost_t * ghost_t + covariance[1][1];
    target_basis[7] = ghost_t * ghost_r + covariance[1][2];
    target_basis[8] = ghost_r * ghost_r + covariance[2][2];
#endif

    Array1D<Real, 0, NIP - 1> trial_weights{};
    Real condition = std::numeric_limits<Real>::infinity();
    int accepted_order = -1;
    if (ibm_fit_target_functional<NQ, NQ, NIP>(
            support_basis, support_mask, target_basis, maximum_condition,
            trial_weights, condition)) {
        accepted_order = 2;
    } else {
        condition = std::numeric_limits<Real>::infinity();
        if (!ibm_fit_target_functional<NL, NQ, NIP>(
                support_basis, support_mask, target_basis, maximum_condition,
                trial_weights, condition)) {
            return -1;
        }
        accepted_order = 1;
    }

    Real sum = Real(0.0);
    Real l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        weights(support) = trial_weights(support);
        sum += trial_weights(support);
        l1 += amrex::Math::abs(trial_weights(support));
    }
    boundary_weight = Real(1.0) - sum;
    l1 += amrex::Math::abs(boundary_weight);
    if (!amrex::Math::isfinite(boundary_weight + condition + l1)) return -1;
    accepted_condition = condition;
    weight_l1 = l1;
    return accepted_order;
}

/// Extrapolate a smooth auxiliary scalar from the visible fluid side to the
/// real ghost-cell centre without prescribing a wall-normal derivative.  The
/// full local quadratic polynomial is a one-sided fluid jet.  An embedded
/// linear jet is built from the identical support so the nonlinear wall
/// closure can reduce order without changing geometry or stencil ownership.
/// In the Euler-slip closure both are applied to
/// sigma(Q(Ubar)) = log(p)-gamma*log(rho); consequently the fitted samples and
/// target share the same finite-volume discrete-data semantics even though
/// sigma itself is not a conserved cell average.
template <int NIP,
          bool AnnularTarget = rz_annular_bic_cell_average,
          typename IndexArray, typename MaskArray,
          typename PointArray, typename VectorArray, typename WeightArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_one_sided_ghost_jet_functional(
    const IndexArray& support_indices, const MaskArray& support_mask,
    const PointArray& wall_point, const PointArray& ghost_point,
    const VectorArray& wall_normal, const VectorArray& tangent1,
    const VectorArray& tangent2,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    Real length_scale, Real maximum_condition, WeightArray& weights,
    WeightArray& linear_weights, Real& accepted_condition, Real& weight_l1,
    Real& linear_condition, Real& linear_weight_l1)
{
    constexpr int NL = AMREX_SPACEDIM + 1;
    constexpr int NQ = (AMREX_SPACEDIM == 2) ? 6 : 10;

    for (int support = 0; support < NIP; ++support) {
        weights(support) = Real(0.0);
        linear_weights(support) = Real(0.0);
    }
    accepted_condition = std::numeric_limits<Real>::infinity();
    weight_l1 = std::numeric_limits<Real>::infinity();
    linear_condition = std::numeric_limits<Real>::infinity();
    linear_weight_l1 = std::numeric_limits<Real>::infinity();
    if (!(length_scale > Real(0.0)) ||
        !amrex::Math::isfinite(length_scale)) {
        return -1;
    }

    Real frame[AMREX_SPACEDIM][AMREX_SPACEDIM] = {};
    Real frame_norm2[AMREX_SPACEDIM] = {};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        frame[0][d] = wall_normal(d);
        frame[1][d] = tangent1(d);
        frame_norm2[0] += frame[0][d] * frame[0][d];
        frame_norm2[1] += frame[1][d] * frame[1][d];
#if (AMREX_SPACEDIM == 3)
        frame[2][d] = tangent2(d);
        frame_norm2[2] += frame[2][d] * frame[2][d];
#endif
    }
    for (int a = 0; a < AMREX_SPACEDIM; ++a) {
        if (!(frame_norm2[a] > Real(0.0))) return -1;
        const Real inverse_norm = Real(1.0) / std::sqrt(frame_norm2[a]);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            frame[a][d] *= inverse_norm;
        }
    }

    Real support_basis[NIP][NQ] = {};
    for (int support = 0; support < NIP; ++support) {
        if (support_mask(support) == 0) continue;
        Real coordinate[AMREX_SPACEDIM] = {};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int index = support_indices(support, d);
            const Real centre =
                prob_lo[d] + (Real(index) + Real(0.5)) * dxyz[d];
            const Real delta = (centre - wall_point(d)) / length_scale;
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                coordinate[a] += delta * frame[a][d];
            }
        }
#if (AMREX_SPACEDIM == 2)
        const Real s = coordinate[0];
        const Real t = coordinate[1];
        support_basis[support][0] = Real(1.0);
        support_basis[support][1] = s;
        support_basis[support][2] = t;
        support_basis[support][3] = s * s;
        support_basis[support][4] = s * t;
        support_basis[support][5] = t * t;
#else
        const Real s = coordinate[0];
        const Real t = coordinate[1];
        const Real r = coordinate[2];
        support_basis[support][0] = Real(1.0);
        support_basis[support][1] = s;
        support_basis[support][2] = t;
        support_basis[support][3] = r;
        support_basis[support][4] = s * s;
        support_basis[support][5] = s * t;
        support_basis[support][6] = s * r;
        support_basis[support][7] = t * t;
        support_basis[support][8] = t * r;
        support_basis[support][9] = r * r;
#endif
    }

    Real coordinate[AMREX_SPACEDIM] = {};
    if constexpr (!AnnularTarget) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const Real delta = (ghost_point(d) - wall_point(d)) / length_scale;
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                coordinate[a] += delta * frame[a][d];
            }
        }
    }
    Real target_basis[NQ] = {};
#if (AMREX_SPACEDIM == 2)
    if constexpr (AnnularTarget) {
        GpuArray<Real, cerisse::ibm::rz_quadratic_moment_count> target{};
        if (!ibm_rz_annular_ghost_target_moments(
                wall_point, ghost_point, wall_normal, tangent1, dxyz,
                length_scale, target)) {
            return -1;
        }
        for (int m = 0; m < NQ; ++m) target_basis[m] = target[m];
    } else {
        const Real s = coordinate[0];
        const Real t = coordinate[1];
        target_basis[0] = Real(1.0);
        target_basis[1] = s;
        target_basis[2] = t;
        target_basis[3] = s * s;
        target_basis[4] = s * t;
        target_basis[5] = t * t;
    }
#else
    const Real s = coordinate[0];
    const Real t = coordinate[1];
    const Real r = coordinate[2];
    target_basis[0] = Real(1.0);
    target_basis[1] = s;
    target_basis[2] = t;
    target_basis[3] = r;
    target_basis[4] = s * s;
    target_basis[5] = s * t;
    target_basis[6] = s * r;
    target_basis[7] = t * t;
    target_basis[8] = t * r;
    target_basis[9] = r * r;
#endif

    Array1D<Real, 0, NIP - 1> linear_trial_weights{};
    if (!ibm_fit_target_functional<NL, NQ, NIP>(
            support_basis, support_mask, target_basis, maximum_condition,
            linear_trial_weights, linear_condition)) {
        return -1;
    }
    Real linear_l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        linear_weights(support) = linear_trial_weights(support);
        linear_l1 += amrex::Math::abs(linear_trial_weights(support));
    }
    if (!amrex::Math::isfinite(linear_condition + linear_l1)) return -1;
    linear_weight_l1 = linear_l1;

    Array1D<Real, 0, NIP - 1> trial_weights{};
    Real condition = std::numeric_limits<Real>::infinity();
    int accepted_order = -1;
    if (ibm_fit_target_functional<NQ, NQ, NIP>(
            support_basis, support_mask, target_basis, maximum_condition,
            trial_weights, condition)) {
        accepted_order = 2;
    } else {
        condition = linear_condition;
        for (int support = 0; support < NIP; ++support) {
            trial_weights(support) = linear_trial_weights(support);
        }
        accepted_order = 1;
    }

    Real l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        weights(support) = trial_weights(support);
        l1 += amrex::Math::abs(trial_weights(support));
    }
    if (!amrex::Math::isfinite(condition + l1)) return -1;
    accepted_condition = condition;
    weight_l1 = l1;
    return accepted_order;
}

/// One-sided entropy extension at an R-Z wall-axis junction.
///
/// A regular axisymmetric scalar is even in signed radius.  In the first
/// radial ring the unrestricted local quadratic basis can therefore fit an
/// unphysical O(r) component from one-sided r>=0 data.  This functional uses
/// the complete degree-two axis-regular scalar basis
///
///   {1, zeta, zeta^2, xi^2},  zeta=(z-z_B)/h, xi=(r-r_axis)/h,
///
/// and its embedded linear basis {1,zeta}.  The target may be either a point
/// or the exact r-weighted annular cell average.  No wall value or normal
/// derivative is imposed: this remains a one-sided fluid entropy jet, with
/// the independent physical regularity condition d(sigma)/dr=0 at r=0.
template <int NIP, bool AnnularTarget = rz_annular_bic_cell_average,
          typename IndexArray, typename MaskArray, typename PointArray,
          typename WeightArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_rz_axis_regular_ghost_jet_functional(
    const IndexArray& support_indices, const MaskArray& support_mask,
    const PointArray& wall_point, const PointArray& ghost_point,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz, Real length_scale,
    Real maximum_condition, WeightArray& weights,
    WeightArray& linear_weights, Real& accepted_condition, Real& weight_l1,
    Real& linear_condition, Real& linear_weight_l1)
{
#if (AMREX_SPACEDIM == 2)
    constexpr int NL = 2;
    constexpr int NQ = 4;
    for (int support = 0; support < NIP; ++support) {
        weights(support) = Real(0.0);
        linear_weights(support) = Real(0.0);
    }
    accepted_condition = std::numeric_limits<Real>::infinity();
    weight_l1 = std::numeric_limits<Real>::infinity();
    linear_condition = std::numeric_limits<Real>::infinity();
    linear_weight_l1 = std::numeric_limits<Real>::infinity();
    if (!(length_scale > Real(0.0)) ||
        !amrex::Math::isfinite(length_scale)) {
        return -1;
    }

    const Real inverse_scale = Real(1.0) / length_scale;
    const Real radial_axis = prob_lo[0];
    Real support_basis[NIP][NQ] = {};
    for (int support = 0; support < NIP; ++support) {
        if (support_mask(support) == 0) continue;
        const Real r =
            prob_lo[0] +
            (Real(support_indices(support, 0)) + Real(0.5)) * dxyz[0];
        const Real z =
            prob_lo[1] +
            (Real(support_indices(support, 1)) + Real(0.5)) * dxyz[1];
        const Real xi = (r - radial_axis) * inverse_scale;
        const Real zeta = (z - wall_point(1)) * inverse_scale;
        support_basis[support][0] = Real(1.0);
        support_basis[support][1] = zeta;
        support_basis[support][2] = zeta * zeta;
        support_basis[support][3] = xi * xi;
    }

    Real target_basis[NQ] = {};
    target_basis[0] = Real(1.0);
    if constexpr (AnnularTarget) {
        const Real roundoff =
            Real(64.0) * std::numeric_limits<Real>::epsilon() *
            amrex::max(Real(1.0), amrex::Math::abs(ghost_point(0)));
        Real r_lo = ghost_point(0) - Real(0.5) * dxyz[0];
        if (amrex::Math::abs(r_lo - radial_axis) <= roundoff) {
            r_lo = radial_axis;
        }
        const cerisse::ibm::RZCellBounds target{
            r_lo, ghost_point(0) + Real(0.5) * dxyz[0],
            ghost_point(1) - Real(0.5) * dxyz[1],
            ghost_point(1) + Real(0.5) * dxyz[1]};
        GpuArray<Real, cerisse::ibm::rz_quadratic_moment_count> moments{};
        if (!cerisse::ibm::rz_annular_quadratic_moments(target, moments)) {
            return -1;
        }
        const Real z_wall = wall_point(1);
        const Real inverse_scale2 = inverse_scale * inverse_scale;
        target_basis[1] =
            (moments[cerisse::ibm::rz_z] - z_wall) * inverse_scale;
        target_basis[2] =
            (moments[cerisse::ibm::rz_z2] -
             Real(2.0) * z_wall * moments[cerisse::ibm::rz_z] +
             z_wall * z_wall) *
            inverse_scale2;
        target_basis[3] =
            (moments[cerisse::ibm::rz_r2] -
             Real(2.0) * radial_axis * moments[cerisse::ibm::rz_r] +
             radial_axis * radial_axis) *
            inverse_scale2;
    } else {
        const Real xi =
            (ghost_point(0) - radial_axis) * inverse_scale;
        const Real zeta =
            (ghost_point(1) - wall_point(1)) * inverse_scale;
        target_basis[1] = zeta;
        target_basis[2] = zeta * zeta;
        target_basis[3] = xi * xi;
    }

    Array1D<Real, 0, NIP - 1> linear_trial_weights{};
    if (!ibm_fit_target_functional<NL, NQ, NIP>(
            support_basis, support_mask, target_basis, maximum_condition,
            linear_trial_weights, linear_condition)) {
        return -1;
    }
    Real linear_l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        linear_weights(support) = linear_trial_weights(support);
        linear_l1 += amrex::Math::abs(linear_trial_weights(support));
    }
    if (!amrex::Math::isfinite(linear_condition + linear_l1)) return -1;
    linear_weight_l1 = linear_l1;

    Array1D<Real, 0, NIP - 1> trial_weights{};
    Real condition = std::numeric_limits<Real>::infinity();
    int accepted_order = -1;
    if (ibm_fit_target_functional<NQ, NQ, NIP>(
            support_basis, support_mask, target_basis, maximum_condition,
            trial_weights, condition)) {
        accepted_order = 2;
    } else {
        condition = linear_condition;
        for (int support = 0; support < NIP; ++support) {
            trial_weights(support) = linear_trial_weights(support);
        }
        accepted_order = 1;
    }

    Real l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        weights(support) = trial_weights(support);
        l1 += amrex::Math::abs(trial_weights(support));
    }
    if (!amrex::Math::isfinite(condition + l1)) return -1;
    accepted_condition = condition;
    weight_l1 = l1;
    return accepted_order;
#else
    amrex::ignore_unused(
        support_indices, support_mask, wall_point, ghost_point, prob_lo,
        dxyz, length_scale, maximum_condition, weights, linear_weights,
        accepted_condition, weight_l1, linear_condition, linear_weight_l1);
    return -1;
#endif
}

/// Build geometry-only functionals for one single-valued shared ghost cell.
/// Coordinates are centred at the GP's own boundary intercept and resolved in
/// its local (normal,tangent) frame.  For a Dirichlet scalar,
///
///   q_G = sum_j w^D_j q_j + w_B q_B,
///
/// while a prescribed wall-normal derivative uses
///
///   q_G = sum_j w^N_j q_j + w_g (dq/dn)_B.
///
/// Both expressions evaluate the constrained polynomial directly at the real
/// ghost-cell centre.  No image-point primitive state or IP-to-GP extrapolation
/// appears in this operator.  The quadratic null spaces are
///
///   Dirichlet: {s,t[,r], all degree-2 monomials},
///   Neumann:   {1,t[,r], all degree-2 monomials},
///
/// where s is the wall-normal coordinate.  Their leading D terms form the
/// explicit linear fallback.
template <int NIP,
          bool AnnularTarget = rz_annular_bic_cell_average,
          typename IndexArray, typename MaskArray,
          typename PointArray, typename VectorArray,
          typename WeightArray>
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ibm_bi_constrained_gp_functionals(
    const IndexArray& support_indices, const MaskArray& support_mask,
    const PointArray& wall_point, const PointArray& ghost_point,
    const VectorArray& wall_normal, const VectorArray& tangent1,
    const VectorArray& tangent2,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxyz,
    Real length_scale, Real maximum_condition,
    WeightArray& dirichlet_weights, Real& dirichlet_boundary_weight,
    WeightArray& neumann_weights, Real& neumann_gradient_weight,
    int& support_count, Real& accepted_condition, Real& weight_l1)
{
    constexpr int NL = AMREX_SPACEDIM;
    constexpr int NQ =
        (AMREX_SPACEDIM == 2) ? 5 : 9;

    for (int support = 0; support < NIP; ++support) {
        dirichlet_weights(support) = Real(0.0);
        neumann_weights(support) = Real(0.0);
    }
    dirichlet_boundary_weight = Real(0.0);
    neumann_gradient_weight = Real(0.0);
    support_count = 0;
    accepted_condition = std::numeric_limits<Real>::infinity();
    weight_l1 = std::numeric_limits<Real>::infinity();
    if (!(length_scale > Real(0.0)) ||
        !amrex::Math::isfinite(length_scale)) {
        return -1;
    }

    Real normal[AMREX_SPACEDIM];
    Real tangent_a[AMREX_SPACEDIM];
    Real tangent_b[AMREX_SPACEDIM];
    Real normal_norm2 = Real(0.0);
    Real tangent_a_norm2 = Real(0.0);
    Real tangent_b_norm2 = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        normal[d] = wall_normal(d);
        tangent_a[d] = tangent1(d);
        tangent_b[d] = tangent2(d);
        normal_norm2 += normal[d] * normal[d];
        tangent_a_norm2 += tangent_a[d] * tangent_a[d];
        tangent_b_norm2 += tangent_b[d] * tangent_b[d];
    }
    if (!(normal_norm2 > Real(0.0)) ||
        !(tangent_a_norm2 > Real(0.0))) {
        return -1;
    }
    const Real inv_normal = Real(1.0) / std::sqrt(normal_norm2);
    const Real inv_tangent_a = Real(1.0) / std::sqrt(tangent_a_norm2);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        normal[d] *= inv_normal;
        tangent_a[d] *= inv_tangent_a;
    }
#if (AMREX_SPACEDIM == 3)
    if (!(tangent_b_norm2 > Real(0.0))) return -1;
    const Real inv_tangent_b = Real(1.0) / std::sqrt(tangent_b_norm2);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        tangent_b[d] *= inv_tangent_b;
    }
#else
    amrex::ignore_unused(tangent_b_norm2);
#endif

    Real dirichlet_basis[NIP][NQ] = {};
    Real neumann_basis[NIP][NQ] = {};
    Real support_normal_coordinate[NIP] = {};
    for (int support = 0; support < NIP; ++support) {
        if (support_mask(support) == 0) continue;
        ++support_count;
        Real s = Real(0.0);
        Real t = Real(0.0);
        Real r = Real(0.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int index = support_indices(support, d);
            const Real centre =
                prob_lo[d] + (Real(index) + Real(0.5)) * dxyz[d];
            const Real delta = centre - wall_point(d);
            s += delta * normal[d] / length_scale;
            t += delta * tangent_a[d] / length_scale;
#if (AMREX_SPACEDIM == 3)
            r += delta * tangent_b[d] / length_scale;
#endif
        }
        support_normal_coordinate[support] = s;
#if (AMREX_SPACEDIM == 2)
        dirichlet_basis[support][0] = s;
        dirichlet_basis[support][1] = t;
        dirichlet_basis[support][2] = s * s;
        dirichlet_basis[support][3] = s * t;
        dirichlet_basis[support][4] = t * t;
        neumann_basis[support][0] = Real(1.0);
        neumann_basis[support][1] = t;
        neumann_basis[support][2] = s * s;
        neumann_basis[support][3] = s * t;
        neumann_basis[support][4] = t * t;
#else
        dirichlet_basis[support][0] = s;
        dirichlet_basis[support][1] = t;
        dirichlet_basis[support][2] = r;
        dirichlet_basis[support][3] = s * s;
        dirichlet_basis[support][4] = s * t;
        dirichlet_basis[support][5] = s * r;
        dirichlet_basis[support][6] = t * t;
        dirichlet_basis[support][7] = t * r;
        dirichlet_basis[support][8] = r * r;
        neumann_basis[support][0] = Real(1.0);
        neumann_basis[support][1] = t;
        neumann_basis[support][2] = r;
        neumann_basis[support][3] = s * s;
        neumann_basis[support][4] = s * t;
        neumann_basis[support][5] = s * r;
        neumann_basis[support][6] = t * t;
        neumann_basis[support][7] = t * r;
        neumann_basis[support][8] = r * r;
#endif
    }

    Real ghost_s = Real(0.0);
    Real ghost_t = Real(0.0);
    Real ghost_r = Real(0.0);
    Real ghost_s2 = Real(0.0);
    Real ghost_st = Real(0.0);
    Real ghost_t2 = Real(0.0);
    if constexpr (AnnularTarget) {
        GpuArray<Real, cerisse::ibm::rz_quadratic_moment_count> target{};
        if (!ibm_rz_annular_ghost_target_moments(
                wall_point, ghost_point, wall_normal, tangent1, dxyz,
                length_scale, target)) {
            return -1;
        }
        ghost_s = target[cerisse::ibm::rz_z];
        ghost_t = target[cerisse::ibm::rz_r];
        ghost_s2 = target[cerisse::ibm::rz_z2];
        ghost_st = target[cerisse::ibm::rz_zr];
        ghost_t2 = target[cerisse::ibm::rz_r2];
    } else {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const Real delta = ghost_point(d) - wall_point(d);
            ghost_s += delta * normal[d] / length_scale;
            ghost_t += delta * tangent_a[d] / length_scale;
#if (AMREX_SPACEDIM == 3)
            ghost_r += delta * tangent_b[d] / length_scale;
#endif
        }
        ghost_s2 = ghost_s * ghost_s;
        ghost_st = ghost_s * ghost_t;
        ghost_t2 = ghost_t * ghost_t;
    }
    Real dirichlet_target[NQ] = {};
    Real neumann_target[NQ] = {};
#if (AMREX_SPACEDIM == 2)
    dirichlet_target[0] = ghost_s;
    dirichlet_target[1] = ghost_t;
    dirichlet_target[2] = ghost_s2;
    dirichlet_target[3] = ghost_st;
    dirichlet_target[4] = ghost_t2;
    neumann_target[0] = Real(1.0);
    neumann_target[1] = ghost_t;
    neumann_target[2] = ghost_s2;
    neumann_target[3] = ghost_st;
    neumann_target[4] = ghost_t2;
#else
    dirichlet_target[0] = ghost_s;
    dirichlet_target[1] = ghost_t;
    dirichlet_target[2] = ghost_r;
    dirichlet_target[3] = ghost_s * ghost_s;
    dirichlet_target[4] = ghost_s * ghost_t;
    dirichlet_target[5] = ghost_s * ghost_r;
    dirichlet_target[6] = ghost_t * ghost_t;
    dirichlet_target[7] = ghost_t * ghost_r;
    dirichlet_target[8] = ghost_r * ghost_r;
    neumann_target[0] = Real(1.0);
    neumann_target[1] = ghost_t;
    neumann_target[2] = ghost_r;
    neumann_target[3] = ghost_s * ghost_s;
    neumann_target[4] = ghost_s * ghost_t;
    neumann_target[5] = ghost_s * ghost_r;
    neumann_target[6] = ghost_t * ghost_t;
    neumann_target[7] = ghost_t * ghost_r;
    neumann_target[8] = ghost_r * ghost_r;
#endif

    Array1D<Real, 0, NIP - 1> dirichlet_trial;
    Array1D<Real, 0, NIP - 1> neumann_trial;
    Real dirichlet_condition = std::numeric_limits<Real>::infinity();
    Real neumann_condition = std::numeric_limits<Real>::infinity();
    int accepted_order = -1;
    const bool quadratic_dirichlet =
        ibm_fit_target_functional<NQ, NQ, NIP>(
            dirichlet_basis, support_mask, dirichlet_target,
            maximum_condition, dirichlet_trial, dirichlet_condition);
    const bool quadratic_neumann =
        ibm_fit_target_functional<NQ, NQ, NIP>(
            neumann_basis, support_mask, neumann_target,
            maximum_condition, neumann_trial, neumann_condition);
    if (quadratic_dirichlet && quadratic_neumann) {
        accepted_order = 2;
    } else {
        dirichlet_condition = std::numeric_limits<Real>::infinity();
        neumann_condition = std::numeric_limits<Real>::infinity();
        const bool linear_dirichlet =
            ibm_fit_target_functional<NL, NQ, NIP>(
                dirichlet_basis, support_mask, dirichlet_target,
                maximum_condition, dirichlet_trial, dirichlet_condition);
        const bool linear_neumann =
            ibm_fit_target_functional<NL, NQ, NIP>(
                neumann_basis, support_mask, neumann_target,
                maximum_condition, neumann_trial, neumann_condition);
        if (!(linear_dirichlet && linear_neumann)) return -1;
        accepted_order = 1;
    }

    Real dirichlet_sum = Real(0.0);
    Real neumann_normal_moment = Real(0.0);
    Real dirichlet_l1 = Real(0.0);
    Real neumann_l1 = Real(0.0);
    for (int support = 0; support < NIP; ++support) {
        dirichlet_weights(support) = dirichlet_trial(support);
        neumann_weights(support) = neumann_trial(support);
        dirichlet_sum += dirichlet_trial(support);
        neumann_normal_moment +=
            neumann_trial(support) * support_normal_coordinate[support];
        dirichlet_l1 += amrex::Math::abs(dirichlet_trial(support));
        neumann_l1 += amrex::Math::abs(neumann_trial(support));
    }
    dirichlet_boundary_weight = Real(1.0) - dirichlet_sum;
    neumann_gradient_weight =
        length_scale * (ghost_s - neumann_normal_moment);
    accepted_condition = amrex::max(dirichlet_condition, neumann_condition);
    weight_l1 = amrex::max(
        dirichlet_l1 + amrex::Math::abs(dirichlet_boundary_weight),
        neumann_l1);
    if (!amrex::Math::isfinite(
            dirichlet_boundary_weight + neumann_gradient_weight +
            accepted_condition + weight_l1)) {
        return -1;
    }
    return accepted_order;
}

// ============================================================================
// 4. computeIPweights
// ============================================================================
template <int eorder_t, int iorder_t, typename IPDATA, int N_InterP = ipow(iorder_t + 1, AMREX_SPACEDIM), int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void computeIPweights(
    Array2D<Real,0,eorder_t-1,0,N_InterP-1>&                     weights,
    Array3D< int,0,eorder_t-1,0,N_InterP-1,0,AMREX_SPACEDIM-1>&  ip_ijk,
    Array2D<Real,0,eorder_t-1,0,AMREX_SPACEDIM-1>&               imp_xyz,
    Array2D< int,0,eorder_t-1,0,AMREX_SPACEDIM-1>&               imp_ijk,
    Array1D< int,0,eorder_t-1>&                                  imp_ninterp,
    Array1D< int,0,eorder_t-1>&                                  imp_fit_order,
    const Array2D<uint8_t,0,eorder_t-1,0,N_InterP-1>&             support_visible,
    const GpuArray<Real, AMREX_SPACEDIM>&                        prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>&                        dxyz,
    const Array4<uint8_t const>&                                 ibFab,
    Real maximum_condition = std::numeric_limits<Real>::infinity(),
    int maximum_fit_order = iorder_t,
    bool support_values_are_cell_averages = false)
{
    constexpr int INTERP_THRESHOLD = (GP_OR_SURF ? INTERP_THRESHOLD_GP : INTERP_THRESHOLD_SURF);
    // An image point below the path-specific threshold is invalid.  For
    // iorder=2, a rank-deficient quadratic/linear stencil above that threshold
    // is retained as an explicit constant-fit fallback; this is needed at
    // non-smooth 3-D STL corners where no local first-order polynomial exists.
    constexpr int MIN_PTS = INTERP_THRESHOLD;
    for (int iim = 0; iim < eorder_t; ++iim) {
      if (imp_ninterp(iim) < MIN_PTS) {
        for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) = Real(0.0);
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = -99;
          }
        }
        imp_ninterp(iim) = 0;
        imp_fit_order(iim) = -1;
        continue;
      }
      int base_ijk[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {base_ijk[d] = imp_ijk(iim, d);}

      if constexpr (iorder_t == 1) {
      // ----------------------------------------------------------------------
      // Bilinear/trilinear path (legacy, bit-identical): tensor-product weights
      // on the 2^D corner cells, solid corners zeroed and renormalised.
      // ----------------------------------------------------------------------
      Real frac[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real lo = prob_lo[d] + Real(base_ijk[d] + 0.5_rt) * dxyz[d];
        frac[d] = (imp_xyz(iim, d) - lo) / dxyz[d];
      }
      int  sumfluid   = 0;
      Real sumweights = Real(0.0);
      for (int corner = 0; corner < N_InterP; ++corner) {
        int  ijk[AMREX_SPACEDIM];
        Real w = Real(1.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int bit = (corner >> d) & 1;
            ijk[d] = base_ijk[d] + bit;
            const Real fd = frac[d];
            w *= (bit ? fd : (Real(1.0) - fd));
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = ijk[d];
        }
        int ii = ijk[0];
        int jj = ijk[1];
#if (AMREX_SPACEDIM == 3)
        int kk = ijk[2];
#else
        int kk = 0;
#endif
        int fluid = !ibFab(ii, jj, kk, 0) && support_visible(iim, corner);
        weights(iim, corner) = w * Real(fluid);
        sumfluid   += fluid;
        sumweights += weights(iim, corner);
      }
      AMREX_ASSERT_WITH_MESSAGE(
          sumweights > Real(0.0),
          "computeIPweights: sum of raw weights is zero (unexpected numerical error).");
      AMREX_ASSERT_WITH_MESSAGE(
          sumfluid == imp_ninterp(iim),
          "computeIPweights: mismatch in fluid stencil count.");
      // Runtime guard: in release builds ASSERT may be stripped; protect against FPE.
      Real inv_sum = (sumweights > Real(1.0e-30))
                   ? Real(1.0) / sumweights
                   : Real(0.0);
      Real check_sum = Real(0.0);
      for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) *= inv_sum;
          check_sum += weights(iim, corner);
      }
      AMREX_ASSERT_WITH_MESSAGE(
          amrex::Math::abs(check_sum - Real(1.0)) < Real(1.0e-9),
          "Interpolation point weights do not sum to 1.0");
      imp_fit_order(iim) = 1;

      } else {
      // ----------------------------------------------------------------------
      // iorder >= 2: least-squares quadratic interpolation on the
      // (iorder+1)^D fluid cells of the block anchored at imp_ijk.  The
      // production GP path uses iorder=2.
      //
      // The monomial basis is centred AT the image point and scaled by dx, so
      // the interpolated point value equals the constant coefficient a0 and
      // the per-cell weights are lambda_i = b(xi_i) . M^{-1} e0 with
      // M = sum_fluid b b^T.
      //
      // For finite-volume support data, b is the exact full Cartesian
      // cell-average functional applied to each monomial. Thus
      // <xi_d^2> = xi_{d,c}^2 + 1/12, while mixed moments are unchanged.
      // This maps conservative cell averages directly to an image-point state
      // without treating Q(Ubar) as a point sample.
      // Constant and linear (and quadratic) reproduction hold by construction — masking solid
      // cells does NOT break partition of unity, unlike the renormalised
      // bilinear path.  Basis ordering keeps the linear terms first, so the
      // linear-basis fallback reuses the leading NBL x NBL block of M:
      //   2D: {1, x, y, x^2, xy, y^2}          NBQ = 6, NBL = 3
      //   3D: {1, x, y, z, x^2, y^2, z^2,
      //        xy, xz, yz}                     NBQ = 10, NBL = 4
      // Demotion ladder: quadratic (n_f >= NBQ and well-conditioned)
      //   -> linear LS (n_f >= NBL) -> constant LS (n_f >= threshold)
      //   -> IP invalidated.  The constant fallback is deliberately exposed
      // through imp_fit_order=0 so sharp-corner coverage is never mistaken for
      // a high-order reconstruction.
      // Register-heavy but init-time only (once per GP per regrid).
      // ----------------------------------------------------------------------
      static_assert(iorder_t >= 2,
                    "computeIPweights: nonlinear WLS needs iorder >= 2");
      constexpr int NBQ = (AMREX_SPACEDIM == 2) ? 6 : 10;
      constexpr int NBL = AMREX_SPACEDIM + 1;

      Real bmat[N_InterP][NBQ] = {};
      bool isfl[N_InterP] = {};
      Real Mq[NBQ][NBQ] = {};
      int  n_f = 0;

      for (int corner = 0; corner < N_InterP; ++corner) {
        int ijk[AMREX_SPACEDIM];
        Real xi[AMREX_SPACEDIM];
        int rem = corner;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int off = rem % (iorder_t + 1);
            rem /= (iorder_t + 1);
            ijk[d] = base_ijk[d] + off;
            ip_ijk(iim, corner, d) = ijk[d];
            const Real xc = prob_lo[d] + (Real(ijk[d]) + Real(0.5)) * dxyz[d];
            xi[d] = (xc - imp_xyz(iim, d)) / dxyz[d];
        }
        int ii = ijk[0];
        int jj = ijk[1];
#if (AMREX_SPACEDIM == 3)
        int kk = ijk[2];
#else
        int kk = 0;
#endif
        Real* b = bmat[corner];
        b[0] = Real(1.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) b[1 + d] = xi[d];
#if (AMREX_SPACEDIM == 2)
        const Real cell_second_moment =
            support_values_are_cell_averages ? Real(1.0 / 12.0) : Real(0.0);
        b[3] = xi[0] * xi[0] + cell_second_moment;
        b[4] = xi[0] * xi[1];
        b[5] = xi[1] * xi[1] + cell_second_moment;
#else
        const Real cell_second_moment =
            support_values_are_cell_averages ? Real(1.0 / 12.0) : Real(0.0);
        b[4] = xi[0] * xi[0] + cell_second_moment;
        b[5] = xi[1] * xi[1] + cell_second_moment;
        b[6] = xi[2] * xi[2] + cell_second_moment;
        b[7] = xi[0] * xi[1];
        b[8] = xi[0] * xi[2];
        b[9] = xi[1] * xi[2];
#endif
        isfl[corner] =
            (ibFab(ii, jj, kk, 0) == 0) && support_visible(iim, corner);
        if (isfl[corner]) {
            ++n_f;
            for (int a = 0; a < NBQ; ++a)
                for (int c = 0; c <= a; ++c)
                    Mq[a][c] += b[a] * b[c];
        }
      }
      if (n_f < MIN_PTS) {
          for (int corner = 0; corner < N_InterP; ++corner) {
              weights(iim, corner) = Real(0.0);
              for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                  ip_ijk(iim, corner, d) = -99;
              }
          }
          imp_ninterp(iim) = 0;
          imp_fit_order(iim) = -1;
          continue;
      }
      for (int a = 0; a < NBQ; ++a)
          for (int c = a + 1; c < NBQ; ++c)
              Mq[a][c] = Mq[c][a];

      int used_nb = 0;
      Real yq[NBQ];
      if (maximum_fit_order >= 2 && n_f >= NBQ) {
          if (ibm_spd_solve_e0<NBQ>(Mq, yq) &&
              (!amrex::Math::isfinite(maximum_condition) ||
               ibm_spd_condition_1<NBQ>(Mq) <= maximum_condition)) {
              used_nb = NBQ;
          }
      }
      if (used_nb == 0 && n_f >= NBL) {
          Real Ml[NBL][NBL];
          Real yl[NBL];
          for (int a = 0; a < NBL; ++a)
              for (int c = 0; c < NBL; ++c)
                  Ml[a][c] = Mq[a][c];
          if (ibm_spd_solve_e0<NBL>(Ml, yl) &&
              (!amrex::Math::isfinite(maximum_condition) ||
               ibm_spd_condition_1<NBL>(Ml) <= maximum_condition)) {
              for (int a = 0; a < NBL; ++a) yq[a] = yl[a];
              used_nb = NBL;
          }
      }

      if (used_nb == 0) {
          for (int corner = 0; corner < N_InterP; ++corner) {
              weights(iim, corner) = isfl[corner]
                  ? Real(1.0) / Real(n_f) : Real(0.0);
          }
          imp_ninterp(iim) = n_f;
          imp_fit_order(iim) = 0;
      } else {
          Real check_sum = Real(0.0);
          for (int corner = 0; corner < N_InterP; ++corner) {
              Real w = Real(0.0);
              if (isfl[corner]) {
                  for (int a = 0; a < used_nb; ++a)
                      w += bmat[corner][a] * yq[a];
              }
              weights(iim, corner) = w;
              check_sum += w;
          }
          imp_ninterp(iim) = n_f;
          imp_fit_order(iim) = (used_nb == NBQ) ? 2 : 1;
          AMREX_ASSERT_WITH_MESSAGE(
              amrex::Math::abs(check_sum - Real(1.0)) < Real(1.0e-8),
              "WLS interpolation weights do not reproduce constants");
          amrex::ignore_unused(check_sum);
      }
      } // end iorder branch
    } // end loop over image points

    if constexpr (iorder_t >= 2) {
      // A constant spatial fit is only first-order accurate for an off-grid
      // image value.  Keep IP1 when it is the only robust local estimate, but
      // do not combine it (or a later constant IP) with a higher-order normal
      // polynomial.  Truncating the contiguous IP sequence makes the existing
      // n_valid/extrapolation machinery expose the intended order reduction.
      bool truncate = (imp_fit_order(0) == 0);
      for (int iim = 1; iim < eorder_t; ++iim) {
        truncate = truncate || (imp_fit_order(iim) == 0);
        if (!truncate) continue;
        for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) = Real(0.0);
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = -99;
          }
        }
        imp_ninterp(iim) = 0;
        imp_fit_order(iim) = -1;
      }
    }
}
// ============================================================================
// 5. interpolateIMs
// ============================================================================
template <int eorder_t, int iorder_t, int N_InterP = ipow(iorder_t + 1, AMREX_SPACEDIM)>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static void interpolateIMs(
    const Array3D< int, 0, eorder_t - 1, 0, N_InterP - 1, 0, AMREX_SPACEDIM-1>&  imp_ip_ijk,
    const Array2D<Real, 0, eorder_t - 1, 0, N_InterP - 1>&                       imp_ipweights,
    const Array4<Real>&                                                          prims,
    Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM-1>&                           primsNormal) noexcept
{
    for (int iim = 0; iim < eorder_t; ++iim) {
        for (int iip = 0; iip < N_InterP; ++iip) {
            const Real w = imp_ipweights(iim, iip);
            if (w == 0.0) continue;
            const int ii = imp_ip_ijk(iim, iip, 0);
            const int jj = imp_ip_ijk(iim, iip, 1);
        #if (AMREX_SPACEDIM == 3)
            const int kk = imp_ip_ijk(iim, iip, 2);
        #else
            const int kk = 0;
        #endif
            for (int n = 0; n < cls_t::NPRIM; ++n) {
                primsNormal(iim + 2, n) += prims(ii, jj, kk, n) * w;
            }
        }
    }
}

template <int eorder_t, int iorder_t,
          int N_InterP = ipow(iorder_t + 1, AMREX_SPACEDIM)>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static bool interpolateConservativeCellAveragesToPointPrimitives(
    const Array3D<int, 0, eorder_t - 1, 0, N_InterP - 1,
                  0, AMREX_SPACEDIM - 1>& imp_ip_ijk,
    const Array2D<Real, 0, eorder_t - 1, 0, N_InterP - 1>&
        imp_ipweights,
    const Array1D<int, 0, eorder_t - 1>& imp_ninterp,
    const Array4<Real>& conservative,
    const cls_t* cls,
    Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM - 1>&
        primsNormal) noexcept
{
    for (int iim = 0; iim < eorder_t; ++iim) {
        if (imp_ninterp(iim) < INTERP_THRESHOLD_SURF) break;

        Real state[cls_t::NCONS]{};
        for (int iip = 0; iip < N_InterP; ++iip) {
            const Real weight = imp_ipweights(iim, iip);
            if (weight == Real(0.0)) continue;
            const int ii = imp_ip_ijk(iim, iip, 0);
            const int jj = imp_ip_ijk(iim, iip, 1);
#if (AMREX_SPACEDIM == 3)
            const int kk = imp_ip_ijk(iim, iip, 2);
#else
            const int kk = 0;
#endif
            for (int component = 0; component < cls_t::NCONS; ++component) {
                state[component] +=
                    weight * conservative(ii, jj, kk, component);
            }
        }

        bool finite = true;
        for (int component = 0; component < cls_t::NCONS; ++component) {
            finite = finite && amrex::Math::isfinite(state[component]);
        }
        const Real density = state[cls_t::URHO];
        if (!finite || !(density > Real(0.0))) return false;
        const Real kinetic =
            Real(0.5) *
            (state[cls_t::UMX] * state[cls_t::UMX] +
             state[cls_t::UMY] * state[cls_t::UMY] +
             state[cls_t::UMZ] * state[cls_t::UMZ]) /
            density;
        const Real internal_energy_density = state[cls_t::UET] - kinetic;
        if (!amrex::Math::isfinite(internal_energy_density) ||
            !(internal_energy_density > Real(0.0))) {
            return false;
        }

        Real point_primitives[cls_t::NPRIM]{};
        cls->cons2prims_point(state, point_primitives);
        for (int component = 0; component < cls_t::NPRIM; ++component) {
            if (!amrex::Math::isfinite(point_primitives[component])) {
                return false;
            }
            primsNormal(iim + 2, component) =
                point_primitives[component];
        }
    }
    return true;
}
// ============================================================================
// 6. extrapolate
// ============================================================================
template <int eorder_t>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static void extrapolate(
    Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM - 1>& prims, 
    const Array1D< int, 0, eorder_t - 1>& imp_ninterp,
    const Real disGP, const Array1D<Real, 0, eorder_t - 1>& disIM)
{
    // Determine effective order based on INTERP_THRESHOLD
    int eff_order = eorder_t;
    for (int k = 0; k < eorder_t; ++k) {
        if (imp_ninterp(k) < INTERP_THRESHOLD_GP) {
            eff_order = k;
            break;
        }
    }
    // Positivity-preserving floor for the ghost extrapolation of positive-definite
    // thermodynamic primitives (temperature, pressure). A steep wall-normal gradient
    // — e.g. a cold injection wall (SRP nozzle, T clamped to the cold throat value)
    // adjacent to a shock-/compression-heated image point — makes the linear/quadratic
    // extrapolation overshoot a positive quantity through zero. A negative ghost
    // temperature then reaches the EOS unclamped (when CLIP_MINTEMP is off) and yields
    // rho = P/(R T) < 0, which poisons the neighbouring fluid reconstruction (NaN).
    // We limit the ghost to at least this fraction of the (positive) surface value
    // rather than masking the result downstream.
    //
    // NOTE: this must be a *positivity rescue* only -- a small floor that catches
    // the genuine overshoot-through-zero, NOT a 10%-of-surface clamp. In an
    // under-expanded jet / nozzle expansion the near-wall ghost legitimately sits
    // well below the surface value (the gas is expanding), and the old 0.1 floor
    // wrongly raised those physical low pressures/temperatures, injecting a
    // spurious near-wall pressure jump and breaking the extrapolation's
    // normal-derivative continuity. 1e-3 only triggers on values heading through
    // zero; at the floor T and P scale together so rho=P/(RT) stays ~surface
    // density (bounded). Legitimate expansion (>0.1% of surface) is preserved.
    constexpr Real GP_POS_FLOOR_FRAC = Real(1.0e-3);
    // Reduce eff_order while the two outermost image points are nearly
    // coincident (vanishing Lagrange denominators). disIM is monotone by
    // construction (chained IP walk with step di), so checking the trailing
    // pair suffices; reducing by one reproduces the pre-existing
    // quadratic->linear fallback exactly.
    {
        constexpr Real eps_dist = Real(1.0e-12);
        while (eff_order >= 2) {
            const Real xa = disIM(eff_order - 2);
            const Real xb = disIM(eff_order - 1);
            if (amrex::Math::abs(xb - xa) > eps_dist * amrex::max(xa, xb)) break;
            --eff_order;
        }
    }

    // only extrapolate up to QLS (Last Species), skipping aux vars like QC, QG, QEINT.
    for (int n = 0; n <= cls_t::QLS; ++n) {

        // Optional mixed-order profile: cubic thermodynamic extension supplies
        // the extra value order lost by heat-flux differentiation, while the
        // velocity extension stays on the better-conditioned quadratic path.
        // Models without this member retain the historical uniform order.
        int component_order = eff_order;
        if constexpr (requires { param::cap_velocity_extrap_order_at_two; }) {
            if (param::cap_velocity_extrap_order_at_two &&
                n >= cls_t::QU && n < cls_t::QU + AMREX_SPACEDIM) {
                component_order = amrex::min(component_order, 2);
            }
        }

        if (component_order == 0) {
            prims(0, n) = prims(1, n);
        } else {
            // Degree-component_order Lagrange polynomial through the
            // wall-normal nodes {0, disIM(0), ..., disIM(component_order-1)}
            // and matching values {prims(1,n), ..., prims(component_order+1,n)},
            // evaluated at the ghost abscissa s = -disGP. Orders 1/2 reproduce
            // the previous linear/quadratic formulas identically; order >= 3
            // is available when enough image points were placed. Runtime
            // component_order <= eff_order <= eorder_t bounds all indices.
            const Real s = -disGP;
            Real acc = Real(0.0);
            for (int m = 0; m <= component_order; ++m) {
                const Real xm = (m == 0) ? Real(0.0) : disIM(m - 1);
                Real Lm = Real(1.0);
                for (int kk = 0; kk <= component_order; ++kk) {
                    if (kk == m) continue;
                    const Real xk = (kk == 0) ? Real(0.0) : disIM(kk - 1);
                    Lm *= (s - xk) / (xm - xk);
                }
                acc += Lm * prims(m + 1, n);
            }
            prims(0, n) = acc;
        }

        // Positivity floor for temperature and pressure (see note above). The
        // surface value prims(1,n) is positive for T and P, so the floor is a
        // positive lower bound that prevents the ghost from crossing zero.
        if (n == cls_t::QT || n == cls_t::QPRES) {
            const Real floor_val = GP_POS_FLOOR_FRAC * prims(1, n);
            if (!amrex::Math::isfinite(prims(0, n)) ||
                prims(0, n) < floor_val) {
                prims(0, n) = floor_val;
            }
        }
    }
}
// ============================================================================
// 7. global2local / local2global
// ============================================================================
template <int eorder_t>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static void global2local(
    int iip,
    Array2D<Real,0,eorder_t+1,0,cls_t::NPRIM-1>& primsNormal,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2)
{
    const Real ux = primsNormal(iip, cls_t::QU);
    const Real uy = primsNormal(iip, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    const Real uz = primsNormal(iip, cls_t::QW);
#endif
    primsNormal(iip, cls_t::QU) =
        ux * norm(0) + uy * norm(1)
#if (AMREX_SPACEDIM == 3)
        + uz * norm(2)
#endif
        ;
    primsNormal(iip, cls_t::QV) =
        ux * tan1(0) + uy * tan1(1)
#if (AMREX_SPACEDIM == 3)
        + uz * tan1(2)
#endif
        ;
#if (AMREX_SPACEDIM == 3)
    primsNormal(iip, cls_t::QW) =
        ux * tan2(0) + uy * tan2(1) + uz * tan2(2);
#endif
}
template <int eorder_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static void local2global(
    int jj,
    Array2D<Real,0,eorder_t+1,0,cls_t::NPRIM-1>& primsNormal,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2)
{
    const Real un  = primsNormal(jj, cls_t::QU);
    const Real ut1 = primsNormal(jj, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    const Real ut2 = primsNormal(jj, cls_t::QW);
#endif
    primsNormal(jj, cls_t::QU) =
        un  * norm(0) + ut1 * tan1(0)
#if (AMREX_SPACEDIM == 3)
      + ut2 * tan2(0)
#endif
      ;
    primsNormal(jj, cls_t::QV) =
        un  * norm(1) + ut1 * tan1(1)
#if (AMREX_SPACEDIM == 3)
      + ut2 * tan2(1)
#endif
      ;
#if (AMREX_SPACEDIM == 3)
    primsNormal(jj, cls_t::QW) =
        un * norm(2) + ut1 * tan1(2) + ut2 * tan2(2);
#endif
}
// ============================================================================
// 8. check_interpolation_stencil
// ============================================================================
template <typename IPDATA, int iorder_t = 1, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
bool check_interpolation_stencil(
    int i, int j, int k, 
    const amrex::Box& bx, 
    int lev,
    const IPDATA& ipData, 
    int f_idx,
    CheckMode mode)
{
#if (AMREX_SPACEDIM == 2)
    bool is_valid = bx.contains(amrex::IntVect(i, j)) &&
                    bx.contains(amrex::IntVect(i+iorder_t, j+iorder_t));
#else
    bool is_valid = bx.contains(amrex::IntVect(i, j, k)) &&
                    bx.contains(amrex::IntVect(i+iorder_t, j+iorder_t, k+iorder_t));
#endif
    if (!is_valid) {
        if (mode == CheckMode::Silent) {
            return false;
        }
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("Interpolation stencil out of box bounds! "
                            "lev=%d base=(%d,%d,%d) f_idx=%d\n",
                            lev, i, j, k, f_idx);
        AMREX_ASSERT(false);
#else
        int current_geom = -1;
        int current_elem = ipData.elemIdx[f_idx];
        if constexpr (GP_OR_SURF == 1) {
             current_geom = ipData.geomIdx[f_idx];
        }
#if (AMREX_SPACEDIM == 3)
        std::printf("Interpolation stencil out of box bounds!\n"
                    "  Level: %d\n"
                    "  Stencil Base: (%d, %d, %d)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n",
                    lev, i, j, k,
                    current_geom, current_elem);
#else
        std::printf("Interpolation stencil out of box bounds!\n"
                    "  Level: %d\n"
                    "  Stencil Base: (%d, %d)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n",
                    lev, i, j,
                    current_geom, current_elem);
#endif
        std::fflush(stdout);
        if (mode == CheckMode::Warn) {
            amrex::Warning("Interpolation stencil out of box bounds!");
        } else if (mode == CheckMode::Abort) {
            amrex::Abort("Interpolation stencil out of box bounds!");
        }
#endif // AMREX_DEVICE_COMPILE
    }
    return is_valid;
}
#endif // IBM_SOLVER_INTERP_H_
