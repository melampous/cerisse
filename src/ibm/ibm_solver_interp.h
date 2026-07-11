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
      if (n_fluid > best_fluid) {
        best_fluid = n_fluid;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            imp_xyz(0, d) = candi_xyz(d);
            imp_ijk(0, d) = candi_ijk(d);
        }
        imp_ninterp(0) = n_fluid; 
        disIM(0) = IMP_FACTOR[attempt] * di;
      }
      constexpr int IDEAL_NINTERP = ipow(iorder_t + 1, AMREX_SPACEDIM);
      if (best_fluid == IDEAL_NINTERP) {
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
    const GpuArray<Real, AMREX_SPACEDIM>&                        prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>&                        dxyz,
    const Array4<uint8_t const>&                                 ibFab)
{
    constexpr int INTERP_THRESHOLD = (GP_OR_SURF ? INTERP_THRESHOLD_GP : INTERP_THRESHOLD_SURF);
    // iorder=2 LS cannot fit anything below the linear-basis minimum D+1;
    // lift the usability gate accordingly.  At either interpolation order an
    // IP below the applicable threshold is invalid, so persist ninterp=0 rather
    // than leaving a nonzero marker paired with all-zero weights.
    constexpr int MIN_PTS = (iorder_t == 1)
        ? INTERP_THRESHOLD
        : ((INTERP_THRESHOLD > AMREX_SPACEDIM + 1) ? INTERP_THRESHOLD
                                                   : AMREX_SPACEDIM + 1);
    for (int iim = 0; iim < eorder_t; ++iim) {
      if (imp_ninterp(iim) < MIN_PTS) {
        for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) = Real(0.0);
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = -99;
          }
        }
        imp_ninterp(iim) = 0;
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
        int fluid = !ibFab(ii, jj, kk, 0);
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

      } else {
      // ----------------------------------------------------------------------
      // iorder == 2: weighted-least-squares quadratic interpolation on the 3^D
      // fluid cells of the block anchored at imp_ijk (containing cell - 1).
      //
      // The monomial basis is centred AT the image point and scaled by dx, so
      // the interpolated value equals the constant coefficient a0 and the
      // per-cell weights are lambda_i = b(xi_i) . M^{-1} e0 with
      // M = sum_fluid b b^T (plain LS, unit weights).  Constant and linear
      // (and quadratic) reproduction hold by construction — masking solid
      // cells does NOT break partition of unity, unlike the renormalised
      // bilinear path.  Basis ordering keeps the linear terms first, so the
      // linear-basis fallback reuses the leading NBL x NBL block of M:
      //   2D: {1, x, y, x^2, xy, y^2}          NBQ = 6, NBL = 3
      //   3D: {1, x, y, z, x^2, y^2, z^2,
      //        xy, xz, yz}                     NBQ = 10, NBL = 4
      // Demotion ladder: quadratic (n_f >= NBQ and well-conditioned)
      //   -> linear WLS (n_f >= NBL) -> IP invalidated (ninterp := 0), which
      // the existing INTERP_THRESHOLD / n_valid machinery then discards.
      // Register-heavy but init-time only (once per GP per regrid).
      // ----------------------------------------------------------------------
      static_assert(iorder_t == 2, "computeIPweights: iorder must be 1 or 2");
      constexpr int NBQ = (AMREX_SPACEDIM == 2) ? 6 : 10;
      constexpr int NBL = AMREX_SPACEDIM + 1;

      Real bmat[N_InterP][NBQ];
      bool isfl[N_InterP];
      Real Mq[NBQ][NBQ] = {};
      int  n_f = 0;

      for (int corner = 0; corner < N_InterP; ++corner) {
        int rem = corner;
        int ijk[AMREX_SPACEDIM];
        Real xi[AMREX_SPACEDIM];
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
        b[3] = xi[0] * xi[0];
        b[4] = xi[0] * xi[1];
        b[5] = xi[1] * xi[1];
#else
        b[4] = xi[0] * xi[0];
        b[5] = xi[1] * xi[1];
        b[6] = xi[2] * xi[2];
        b[7] = xi[0] * xi[1];
        b[8] = xi[0] * xi[2];
        b[9] = xi[1] * xi[2];
#endif
        isfl[corner] = (ibFab(ii, jj, kk, 0) == 0);
        if (isfl[corner]) {
            ++n_f;
            for (int a = 0; a < NBQ; ++a)
                for (int c = 0; c <= a; ++c)
                    Mq[a][c] += b[a] * b[c];
        }
      }
      for (int a = 0; a < NBQ; ++a)
          for (int c = a + 1; c < NBQ; ++c)
              Mq[a][c] = Mq[c][a];

      int used_nb = 0;
      Real yq[NBQ];
      if (n_f >= NBQ) {
          if (ibm_spd_solve_e0<NBQ>(Mq, yq)) used_nb = NBQ;
      }
      if (used_nb == 0 && n_f >= NBL) {
          Real Ml[NBL][NBL];
          Real yl[NBL];
          for (int a = 0; a < NBL; ++a)
              for (int c = 0; c < NBL; ++c)
                  Ml[a][c] = Mq[a][c];
          if (ibm_spd_solve_e0<NBL>(Ml, yl)) {
              for (int a = 0; a < NBL; ++a) yq[a] = yl[a];
              used_nb = NBL;
          }
      }

      if (used_nb == 0) {
          // Not enough resolvable fluid data even for a linear fit: invalidate
          // this image point so n_valid / eff_order treat it as missing.
          for (int corner = 0; corner < N_InterP; ++corner) {
              weights(iim, corner) = Real(0.0);
          }
          imp_ninterp(iim) = 0;
      } else {
          Real check_sum = Real(0.0);
          for (int corner = 0; corner < N_InterP; ++corner) {
              Real w = Real(0.0);
              if (isfl[corner]) {
                  for (int a = 0; a < used_nb; ++a) w += bmat[corner][a] * yq[a];
              }
              weights(iim, corner) = w;
              check_sum += w;
          }
          imp_ninterp(iim) = n_f;
          AMREX_ASSERT_WITH_MESSAGE(
              amrex::Math::abs(check_sum - Real(1.0)) < Real(1.0e-8),
              "WLS interpolation weights do not reproduce constants");
          amrex::ignore_unused(check_sum);
      }
      } // end iorder branch
    } // end loop over image points
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
            if (prims(0, n) < floor_val) prims(0, n) = floor_val;
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
#else
    const Real uz = Real(0.0);
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
#else
    const Real ut2 = Real(0.0);
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
