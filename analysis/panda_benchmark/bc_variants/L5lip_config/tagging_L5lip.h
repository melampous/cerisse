////////////////////////////////// tagging /////////////////////////////////////
// x_fc-focused refinement.  max_level = 5  (finest dx = 49.6 um = D_e/512).
//
//   L5 is FORCED ONLY, and only where the measurement needs it:
//        - the nozzle-lip band, widened downstream so it follows the initial
//          Prandtl-Meyer turning of the jet boundary (nu(M=1.42) = 9.6 deg,
//          hence the 0.17 slope);
//        - a thin strip on the axis, so x_fc is extracted at 49.6 um.
//     The density sensor is switched off at level 4, so it can NEVER build
//     L5 on the downstream shock cells.  In the production tagging that
//     sensor-driven L5/L6 carried 86% of the cost and none of the signal.
//
//   L4 (99.2 um) covers the whole first shock cell (z <= 2D, r <= 1D) plus
//     sensor-driven refinement inside tight caps, so the first cell is fully
//     resolved and the coarsening interface sits downstream of x_1.
//
// Justification for dropping L6: discretising every profile in the sweep onto
// the real cell centres and integrating the RZ mass flux gives a discharge-
// coefficient error <= 6.7e-5 at 49.6 um, which maps through dx_fc/dC_d = 0.504
// to a bias <= 3.4e-5 D_e -- 30x below the x_fc extraction resolution and
// 1500x below the s-ladder signal of 0.052 D_e.
template <typename TagFab, typename SDataFab, typename GeomData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, TagFab& tagfab,
                  const SDataFab& sdatafab, const GeomData& geomdata,
                  const ProbParm& pparm, int level)
{
  if (nt <= 0) return;

  using idx = indicies_t;
  const Real rho_c = sdatafab(i, j, k, idx::URHO);
  const Real drhor = amrex::Math::abs(sdatafab(i+1,j,k,idx::URHO)
                                    - sdatafab(i-1,j,k,idx::URHO))
                     / (Real(2.0) * rho_c);
  const Real drhoz = amrex::Math::abs(sdatafab(i,j+1,k,idx::URHO)
                                    - sdatafab(i,j-1,k,idx::URHO))
                     / (Real(2.0) * rho_c);
  const Real grad_rho = std::sqrt(drhor*drhor + drhoz*drhoz);

  const Real xr = geomdata.ProbLo()[0] + (Real(i)+Real(0.5))*geomdata.CellSize()[0];
  const Real xz = geomdata.ProbLo()[1] + (Real(j)+Real(0.5))*geomdata.CellSize()[1];
  const Real D  = Real(2.0) * pparm.r_jet;          // 25.4 mm

  // ---- per-level spatial caps: tagging at lev l creates level l+1.
  //                 lev:      0        1        2        3        4
  //                          7D       5D       3D     1.5D     0.90D
  const Real rcap[5] = { Real(0.17780), Real(0.12700), Real(0.07620),
                         Real(0.03810), Real(0.022860) };
  //                         24D      16D       8D       3D      1.5D
  const Real zcap[5] = { Real(0.60960), Real(0.40640), Real(0.20320),
                         Real(0.07620), Real(0.038100) };
  const int  lc = (level < 5) ? level : 4;
  const bool allowed = (xr <= rcap[lc]) && (xz <= zcap[lc]);

  // ---- density sensor: levels 0-3 only.  L5 is never sensor-driven.
  const Real thresh = (level < 3) ? Real(0.05) : Real(0.03);
  const bool sensor = (level < 4) && (grad_rho > thresh);

  // ---- forced L4: the whole first shock cell.
  const bool forceCore = (level < 4)
                      && (xz <= Real(2.0)*D) && (xr <= Real(1.0)*D);

  // ---- forced L5: lip band + axis strip, x <= 1.5 D_e (downstream of x_1).
  //      2.5 mm inside the lip covers the thickest profile in the matrix
  //      (wall tanh a = 679.5 um reaches 0.99 U_e at R - 1.9 mm).
  const Real rlo_lip = pparm.r_jet - Real(2.5e-3);
  const Real rhi_lip = pparm.r_jet + Real(2.5e-3) + Real(0.17) * xz;
  const bool inLip   = (xr >= rlo_lip) && (xr <= rhi_lip);
  const bool onAxis  = (xr <= Real(2.0e-3));
  const bool forceLip = (level == 4) && (xz <= Real(1.5)*D)
                     && (inLip || onAxis);

  const bool tag = (sensor || forceCore || forceLip) && allowed;
  tagfab(i, j, k) = tagfab(i, j, k) || tag;
}
