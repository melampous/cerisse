#ifndef IBM_SOLVER_H_
#define IBM_SOLVER_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Scan.H>
#include <AMReX_Reduce.H>
#include <AMReX_iMultiFab.H>

#include "ibm_containers.h"

#include <iomanip>  // std::setprecision, std::fixed
#include <sstream>

//===================================================================================
///-------------------------------- main class --------------------------------------
///
/// \brief ibm_solver_t is explicit geometry (triangulation based) immersed boundary method
/// class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
/// the geometry
///
template <typename wallmodel, typename param, typename cls_t>
class ibm_solver_t
{
public:
  // constant factor for image point
  static constexpr int  iorder_tparm = param::interp_order; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm = param::extrap_order; // number of image points used for ghost point extrapolation
  static_assert(param::interp_order == 1 || param::interp_order == 2,
      "IBM interp_order: 1 = bi/trilinear, 2 = WLS quadratic; no other value is implemented");
  static_assert(param::interp_order_surf == 1 || param::interp_order_surf == 2,
      "IBM interp_order_surf: 1 = bi/trilinear, 2 = WLS quadratic; no other value is implemented");
  static constexpr Real cim = param::alpha;
  
  static constexpr int  iorder_tparm_surf = param::interp_order_surf; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm_surf = param::extrap_order_surf; // number of image points used for surface reconstruction
  static constexpr Real cim_surf  = param::alpha_surf;

  // number of ghost layers needed for IB method
  static constexpr int  ghost_layers = param::ghost_layers;
  static_assert(ghost_layers == 1,
      "Volume IBM currently requires ghost_layers=1. ghost_layers=0 applies no "
      "ghost-cell wall condition, and deeper GP bands are not yet validated.");
  static_assert(ghost_layers <= cls_t::NGHOST,
      "IBM ghost_layers exceeds cls_t::NGHOST");
  
  // true: interior of closed geometry is solid; false: interior is fluid
  static constexpr bool interior_is_solid = param::interior_is_solid; 

  // ideal number of interpolation points for each image point(ghost point extrapolation and surface reconstruction)
  static constexpr int  N_InterP      = ipow(iorder_tparm + 1, AMREX_SPACEDIM);
  static constexpr int  N_InterP_surf = ipow(iorder_tparm_surf + 1, AMREX_SPACEDIM);

  using SURFIMP = surfImp_t<eorder_tparm_surf, iorder_tparm_surf>;
  using GPSTORE = GPStore<eorder_tparm, iorder_tparm>;
  using GPSTOREVIEW = GPStoreView<eorder_tparm, iorder_tparm>;

  // MultiFabs pointer to Amr class instance
  Amr* amr_p;
  // Per-level marker MultiFabs (uint8_t, 2 components: solid + ghost)
  Vector<IBMultiFab<uint8_t>*> bmf_a;
  // Level-wide flattened ghost-point storage (CSR indexed by local FAB)
  Vector<GPSTORE> gpstore_a;

  // parameters for cell size and refinement ratio
  Vector<IntVect> rratio_a;                                 // vector of refinement ratio per level in each direction
  Vector<GpuArray<Real, AMREX_SPACEDIM>> dx_a;              // vector of cell sizes per level in each direction
  Vector<Real> diag_a;                                      // vector of cell diagonal length per level
  Vector<Real> di_a;                                        // image point distance per level for ghost point extrapolation
  Vector<Real> di_a_surf;                                   // image point distance per level for surface reconstruction

  // geometry related data
  int ngeom = 0;                                            // number of geometries
  Vector<GeomType> geom_a;                                  // IB explicit geometry
#ifdef AMREX_USE_CGAL
  Vector<Tree> tree_a;                                      // CGAL AABB tree per geometry
  Vector<PrimitiveIndexMap> idxmap_a;                       // CGAL primitive-to-index map per geometry
#else
  Vector<BVH> bvh_a;                                        // BVH per geometry (replaces CGAL AABB tree)
#endif
  Vector<std::unique_ptr<inside_t>> inout_fa;                // in out testing function per geometry
  Vector<Bbox> bbox_a;                                      // bounding box per geometry (world-frame, for fast rejection)
  Vector<Bbox> bbox_body_a;                                 // bounding box per geometry (body-frame, constant after init)
  Vector<RigidTransform> transform_a;                       // rigid-body transform per geometry (body→world)
  Vector<std::string> geom_names;                           // geometry names stripped from input filenames

  Gpu::ManagedVector<LocalFrame> LocalFrame_a;              // local orthonormal frame matrix (flattened, body-frame)
  Gpu::ManagedVector<SurfElem> SurfElem_a;                  // surface element area and coordinates (flattened, body-frame)
  Gpu::ManagedVector<int> geom_offsets;                     // Start index for each geometry in flattened arrays
 
  // surface related data
  int ntotalfaces = 0;                                      // number of faces/edges across all geometries
  SURFIMP surfimp_soa;                                      // face/edge image point information (SoA structure)
  surfPhys_t surfphys_soa;                                  // face/edge physical data and identification (SoA structure)
  Vector<FaceCSR> faces_per_level;                          // faces integers per fab and per level [lev] (CSR format)

  /** 
   * \brief Destructor to release allocated memory
   */
  ~ibm_solver_t() noexcept
  {
    // Release per-level IBMultiFab pointers if any remain
    for (auto*& p : bmf_a) {
      if (p) { delete p; p = nullptr; }
    }

    // Release inside/outside testers (unique_ptr handles cleanup automatically)
    inout_fa.clear();
  }

  /**
   * \brief Explicitly release all GPU-managed memory before amrex::Finalize().
   *
   * Must be called while AMReX arenas are still alive. The global inline
   * IBM::ib outlives amrex::Finalize(), so its implicit destructor would
   * free arena memory after the arena is destroyed (static destruction
   * order problem).  Calling cleanup() first leaves the destructor with
   * nothing to free.
   */
  void cleanup() noexcept
  {
    // Delete raw-pointer members first (their internals use arena memory)
    for (auto*& p : bmf_a)    { if (p) { delete p; p = nullptr; } }
    inout_fa.clear();  // unique_ptr handles cleanup

    // Swap-with-empty idiom: guarantees capacity→0 and arena memory freed NOW,
    // so the post-Finalize destructor finds nothing to deallocate.
    { decltype(bmf_a)        tmp; tmp.swap(bmf_a);        }
    { decltype(inout_fa)     tmp; tmp.swap(inout_fa);     }
    { decltype(LocalFrame_a) tmp; tmp.swap(LocalFrame_a); }
    { decltype(SurfElem_a)   tmp; tmp.swap(SurfElem_a);   }
    { decltype(geom_offsets) tmp; tmp.swap(geom_offsets);  }

    // Compound types: destroying elements calls their ManagedVector destructors
    { decltype(geom_a) tmp; tmp.swap(geom_a); }
#ifdef AMREX_USE_CGAL
    { decltype(tree_a)   tmp; tmp.swap(tree_a);   }
    { decltype(idxmap_a) tmp; tmp.swap(idxmap_a); }
#else
    { decltype(bvh_a) tmp; tmp.swap(bvh_a); }
#endif
    { decltype(bbox_a)      tmp; tmp.swap(bbox_a);      }
    { decltype(bbox_body_a) tmp; tmp.swap(bbox_body_a); }
    { decltype(transform_a) tmp; tmp.swap(transform_a); }
    { decltype(geom_names)  tmp; tmp.swap(geom_names);  }
    { decltype(gpstore_a)  tmp; tmp.swap(gpstore_a);  }

    // These structs have their own clear() with shrink_to_fit()
    surfimp_soa.clear();
    surfphys_soa.clear();

    { decltype(faces_per_level) tmp; tmp.swap(faces_per_level); }

    amr_p = nullptr;
  }

  /**
   * \brief Initializes the Immersed Boundary (IB) method structures and geometry.
   *
   * This function sets up the AMR pointer, resizes internal data structures based on the maximum AMR level,
   * computes grid metrics (cell sizes, diagonals) for all levels, and loads the IB geometry from files.
   *
   * \param pointer_amr Pointer to the main Amr class instance.
   */
  void init(Amr* pointer_amr)
  {
    amr_p = pointer_amr;
    rratio_a = amr_p->refRatio();
    int lmax = amr_p->maxLevel();

    bmf_a.resize(lmax + 1);
    gpstore_a.resize(lmax + 1);
    faces_per_level.resize(lmax + 1);

    dx_a.resize(lmax + 1);
    dx_a[0] = amr_p->Geom(0).CellSizeArray();
    for (int i = 1; i <= lmax; i++) {
      for (int j = 0; j < AMREX_SPACEDIM; j++) {
        dx_a[i][j] = dx_a[i - 1][j] / rratio_a[i - 1][j];
      }
    }

    di_a.resize(lmax + 1);
    di_a_surf.resize(lmax + 1);
    diag_a.resize(lmax + 1);

    for (int i = 0; i <= lmax; i++) {
      diag_a[i] = std::sqrt(
        AMREX_D_TERM( std::pow(dx_a[i][0], 2),
                      + std::pow(dx_a[i][1], 2),
                      + std::pow(dx_a[i][2], 2)) );

      di_a[i] = cim * diag_a[i];
      di_a_surf[i] = cim_surf * diag_a[i];
    }

    // read geometry from file
    Real t_geom = amrex::second();
    read_geom();
    t_geom = amrex::second() - t_geom;
    amrex::Print() << "  IBM init: read_geom       = " << t_geom << " s\n";

    // One-time notice: at extrap_order>=2 a custom wall model that only
    // implements a legacy compute_surfIB signature never receives disIM /
    // n_valid, so its Neumann quantities (zero-grad P/T, slip u_t) silently
    // stay 1st-order while ghost extrapolation runs at 2nd. Make that loud.
    if constexpr (eorder_tparm >= 2) {
      using eo_tag  = std::integral_constant<int, eorder_tparm>;
      using Vec1D   = Array1D<Real, 0, AMREX_SPACEDIM - 1>;
      using DisArr  = Array1D<Real, 0, eorder_tparm - 1>;
      using Prims2D = Array2D<Real, 0, eorder_tparm + 1, 0, cls_t::NPRIM - 1>;
      constexpr bool wm_2nd = ibm_detail::is_detected_v<
          ibm_detail::compute_surfIB_expr, wallmodel,
          eo_tag, const Vec1D&, const Vec1D&, const Vec1D&, const Vec1D&,
          const DisArr&, int, Prims2D&, const int, const cls_t*>;
      if constexpr (!wm_2nd) {
        amrex::Print() << "  IBM note: extrap_order=" << eorder_tparm
                       << " but the wall model has only a legacy compute_surfIB signature;\n"
                       << "            wall Neumann BC stays 1st-order (see ibm_walltypes.h for the\n"
                       << "            disIM-aware template to adopt).\n";
      }
    }
  }

  /**
   * \brief create IBMultiFab at a level and store pointers to it
   * \param bxa BoxArray for the level
   * \param dm DistributionMapping for the level
   * \param lev The AMR level
   */
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    // Prevent memory leak if MF already exists
    if (lev < 0 || lev >= static_cast<int>(bmf_a.size())) {
        amrex::Abort("ibm_solver_t::build_mf: lev is out of bounds");
    }

    // Prevent memory leak if MF already exists
    if (bmf_a[lev] != nullptr) {
        destroy_mf(lev);
    }

    // Use default MFInfo (configured to use The_Managed_Arena in IBMultiFab.h)
    bmf_a[lev] = new IBMultiFab<uint8_t>(bxa, dm, 2, cls_t::NGHOST);
  }

  /**
   * \brief destroy IBMultiFab at a level
   * \param lev The AMR level of the IBMultiFab to be destroyed.
   */
  void destroy_mf(int lev)
  {
    if (lev >= 0 && lev < static_cast<int>(bmf_a.size())) { 
      // safe delete: delete nullptr is valid in C++
      delete bmf_a[lev]; 
      bmf_a[lev] = nullptr;  // Prevent dangling pointer
    }
    if (lev >= 0 && lev < static_cast<int>(gpstore_a.size())) {
      gpstore_a[lev].clear();
    }
  }

  /**
   * \brief Rebuild all geometry-level data from current vertex positions.
   *
   * After externally modifying vertex positions in geom_a (e.g. for moving
   * geometry), call this to reconstruct BVH trees, inside/outside testers,
   * bounding boxes, and the surface cache (LocalFrame_a, SurfElem_a).
   * Must be called BEFORE rebuildIBM().
   */
  /**
   * \brief Lightweight rigid-body transform update (replaces full BVH rebuild).
   *
   * For rigid-body FSI, the geometry shape never changes — only its position
   * and orientation.  Instead of moving all vertices and rebuilding the BVH
   * tree, InsideTester, LocalFrame, and SurfElem every time step, we store
   * the geometry in its reference (body) frame and only update the transform.
   *
   * All query functions (computeMarkers, initialiseGPs, computeAllGPs, etc.)
   * apply the inverse transform to query points before using the body-frame
   * BVH, then forward-transform the results back to the world frame.
   * The results are mathematically identical to a full rebuild.
   *
   * \param geomIdx  Index of the geometry to update.
   * \param T        New rigid-body transform (body → world).
   */
  void updateRigidTransform(int geomIdx, const RigidTransform& T)
  {
      transform_a[geomIdx] = T;
      // Update world-frame bounding box from body-frame bbox + new transform
      bbox_a[geomIdx] = T.transform_bbox(bbox_body_a[geomIdx]);
  }

  /**
   * \brief Sanity-check that surface normals point from solid to fluid.
   *
   * Samples N faces per geometry, computes a probe point at
   * centroid + eps*normal (body frame), and tests inside/outside.  If the
   * probe lands inside the body, the normal is pointing the wrong way,
   * which usually means the user's `interior_is_solid` parameter does not
   * match the mesh's actual orientation (e.g., STL exported with reversed
   * face winding, or polygon file in the wrong CCW/CW convention).
   *
   * Warns on any failure; aborts if the majority of samples fail (clear
   * miscofiguration vs. occasional boundary-grazing false positive).
   *
   * Called once at the end of rebuildGeometryData() — after convert_inout
   * has had its chance to flip normals based on `interior_is_solid`.
   */
  void verify_normal_orientation()
  {
    BL_PROFILE("IBM::verify_normal_orientation");
    constexpr int N_SAMPLES_PER_GEOM = 20;
    constexpr Real EPS_FACTOR        = Real(0.01);  // 1% of body diagonal

    for (int ii = 0; ii < ngeom; ++ii) {
      const int n_faces = geom_offsets[ii + 1] - geom_offsets[ii];
      if (n_faces == 0) continue;

      // Length scale: body-frame bounding box diagonal.
      Real diag2 = Real(0.0);
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real ext = bbox_max_d(bbox_body_a[ii], d)
                       - bbox_min_d(bbox_body_a[ii], d);
        diag2 += ext * ext;
      }
      const Real eps = EPS_FACTOR * std::sqrt(diag2);

      const int n_samples = std::min(N_SAMPLES_PER_GEOM, n_faces);
      const int step      = std::max(1, n_faces / n_samples);

      int n_failures = 0;
      int n_checked  = 0;
      for (int local_f = 0; local_f < n_faces; local_f += step) {
        const int f_idx = geom_offsets[ii] + local_f;
        const auto& se  = SurfElem_a[f_idx];
        const auto& lf  = LocalFrame_a[f_idx];

        Point probe;
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        probe = Point(se.centroid[0] + eps * lf.normal[0],
                      se.centroid[1] + eps * lf.normal[1]);
#else
        probe = Point(se.centroid[0] + eps * lf.normal[0],
                      se.centroid[1] + eps * lf.normal[1],
                      se.centroid[2] + eps * lf.normal[2]);
#endif
#else
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          probe[d] = se.centroid[d] + eps * lf.normal[d];
        }
#endif

        BoundedSide side = (*inout_fa[ii])(probe);
        if (side == BoundedSide::Inside) ++n_failures;
        ++n_checked;
      }

      if (n_failures > 0) {
        amrex::Print()
            << "[IBM] WARNING: normal-orientation sanity check failed for geometry " << ii
            << ": " << n_failures << "/" << n_checked
            << " sample probes (centroid + eps*normal) landed inside the body.\n"
            << "       eps = " << eps << " (body-frame, ~"
            << int(EPS_FACTOR * Real(100)) << "% of diag).\n"
            << "       Likely cause: the 'interior_is_solid' parameter does not match\n"
            << "       the mesh's actual outward-normal orientation. STL files with\n"
            << "       inverted face windings, or polygon files in the wrong CCW/CW\n"
            << "       convention, will trigger this. Visualise the mesh in ParaView\n"
            << "       and verify face normals point away from the solid interior.\n";

        if (n_failures > n_checked / 2) {
          amrex::Abort("[IBM] verify_normal_orientation: majority of samples failed; "
                       "check `interior_is_solid` parameter and mesh orientation.");
        }
      }
    }
  }

  /**
   * \brief Full geometry rebuild (for deformable bodies or re-initialization).
   *
   * This is the expensive path that rebuilds BVH, InsideTester, LocalFrame,
   * and SurfElem from scratch.  For rigid-body FSI, use updateRigidTransform()
   * instead — it is O(1) per geometry.
   */
  void rebuildGeometryData()
  {
    LocalFrame_a.clear();
    SurfElem_a.clear();

    for (int i = 0; i < ngeom; i++) {
#ifdef AMREX_USE_CGAL
      tree_a[i].clear();
#if (AMREX_SPACEDIM == 2)
      tree_a[i].insert(geom_a[i].edges_begin(), geom_a[i].edges_end());
#elif (AMREX_SPACEDIM == 3)
      tree_a[i].insert(faces(geom_a[i]).first, faces(geom_a[i]).second, geom_a[i]);
#endif
      tree_a[i].build();

      inout_fa[i] = std::make_unique<inside_t>(geom_a[i]);
#else
      bvh_a[i].build(geom_a[i]);

      inout_fa[i] = std::make_unique<inside_t>(geom_a[i], bvh_a[i]);
#endif

#if defined(AMREX_USE_CGAL) && (AMREX_SPACEDIM == 3)
      bbox_body_a[i] = PMP::bbox(geom_a[i]);
#else
      bbox_body_a[i] = geom_a[i].bbox();
#endif
      bbox_a[i] = transform_a[i].transform_bbox(bbox_body_a[i]);

      geom_offsets[i] = static_cast<int>(LocalFrame_a.size());
#ifdef AMREX_USE_CGAL
      build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                           idxmap_a[i], geom_offsets[i], i);
#elif defined(AMREX_USE_GPU)
      build_geometry_cache_gpu(geom_a[i], SurfElem_a, LocalFrame_a, i);
#else
      build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                           geom_offsets[i], i);
#endif
    }
    geom_offsets[ngeom] = static_cast<int>(LocalFrame_a.size());

    if (!interior_is_solid) {
      convert_inout(LocalFrame_a);
    }

    // Defensive check: detect mismatched `interior_is_solid` setting or
    // reversed mesh winding before downstream IBM machinery reads
    // LocalFrame_a / inout_fa.
    verify_normal_orientation();
  }

  /**
   * \brief Compute solid/fluid markers and identify ghost points on a level.
   *
   * Per cell, writes two uint8_t components into the IBMultiFab:
   *   comp 0 — solid marker:  0 = fluid, ii+1 = solid in geometry ii
   *   comp 1 — ghost marker:  0 = not a ghost point, otherwise carries the
   *                           geometry id of the owning solid cell
   *
   * Pipeline (per FAB):
   *   FAST-SKIP : reject FABs whose working region misses every body bbox.
   *   Step 1    : solid markers via fast bbox-reject + BVH inside test.
   *   Step 2    : ghost markers + ghost-point count (fused reduction).
   *
   * Per-FAB ghost-point counts feed the level-wide CSR GPStore allocation
   * at the end of the function.
   *
   * \param lev AMR level to process.
   */
  void computeMarkers(int lev)
  {
    BL_PROFILE("IBM::computeMarkers");

    // Step 2 reads markers in grow(bx, GP_BOX_EXTRA) that were written by
    // Step 1 in grow(bx, NGHOST). If GP_BOX_EXTRA > NGHOST, Step 2 would
    // read uninitialised marker memory — guard the invariant at compile time.
    static_assert(GP_BOX_EXTRA <= cls_t::NGHOST,
        "GP_BOX_EXTRA must not exceed cls_t::NGHOST: Step 2 would read "
        "marker cells that Step 1 never initialised.");

    // ---- Level-wide handles ----------------------------------------------
    auto& mfab        = *bmf_a[lev];
    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto dx_lev  = dx_a[lev];

    // ---- Per-FAB GP counts (drive level-wide CSR allocation below) -------
    const int nfabs_local = mfab.local_size();
    Vector<int> gp_counts(nfabs_local, 0);
    int ifab_local = 0;

#ifdef AMREX_USE_GPU
    // Batched per-FAB GP counters: Step 2 accumulates into one device array
    // and a single D2H copy after the MFIter loop replaces the previous
    // blocking ReduceData::value() (4 B D2H + stream sync) per FAB.
    // Pre-zeroed from gp_counts (all zeros) so fast-skipped FABs, whose
    // slot no kernel ever touches, read back 0.
    Gpu::DeviceVector<int> d_gp_counts(nfabs_local);
    Gpu::copyAsync(Gpu::hostToDevice, gp_counts.begin(), gp_counts.end(),
                   d_gp_counts.begin());
    int* const p_gp_counts = d_gp_counts.data();
#endif

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      const Box&  bx        = mfi.tilebox();
      const auto& ibMarkers = mfab.array(mfi);

      // ===================================================================
      // FAST-SKIP — short-circuit FABs that are entirely freestream.
      //
      // If the FAB's working region (valid box + the ghost cells consumed
      // downstream) is disjoint from every body's world-frame AABB, every
      // cell must be fluid; setVal(0) is mathematically equivalent to the
      // full kernel but skips the GPU launch and per-cell BVH query.
      //
      // Working region grows valid box by max(NGHOST, GP_BOX_EXTRA) + 2:
      //   NGHOST       — stencil width read by compute_rhs
      //   GP_BOX_EXTRA — extra growth for ghost-point detection (Step 2)
      //   +2           — safety margin for moving-geometry bbox inflation
      //
      // bbox_a[ii] is in the world frame: built once for static geometry,
      // refreshed by the caller every regrid for FSI before this runs.
      // ===================================================================
      {
        constexpr int check_grow = std::max(int(cls_t::NGHOST), GP_BOX_EXTRA) + 2;
        const Box bxcheck = amrex::grow(bx, check_grow);

        bool any_touches = false;
        for (int ii = 0; ii < ngeom; ++ii) {
          bool intersects = true;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const Real box_lo_d = prob_lo[d] +  bxcheck.smallEnd(d)      * dx_lev[d];
            const Real box_hi_d = prob_lo[d] + (bxcheck.bigEnd(d) + 1)   * dx_lev[d];
            if (box_hi_d < bbox_min_d(bbox_a[ii], d) ||
                box_lo_d > bbox_max_d(bbox_a[ii], d)) {
              intersects = false;
              break;
            }
          }
          if (intersects) { any_touches = true; break; }
        }

        if (!any_touches) {
          // Pure freestream: clear both marker components over the FAB's
          // full allocated region (valid + ghost), then move on.
          mfab.get(mfi).template setVal<amrex::RunOn::Device>(uint8_t(0));
          gp_counts[ifab_local] = 0;
          continue;
        }
      }

#ifdef AMREX_USE_GPU
      // ===================================================================
      // GPU path
      // ===================================================================
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
          "ngeom exceeds MAX_NGEOM; raise MAX_NGEOM in ibm_containers.h");
      const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);

      // Pre-pack trivially-copyable POD views for device capture.
      GpuArray<InsideTesterView, MAX_NGEOM> inout_views;
      GpuArray<AABB,             MAX_NGEOM> bboxes;
      GpuArray<RigidTransform,   MAX_NGEOM> transforms;
      for (int ii = 0; ii < ngeom_local; ++ii) {
        inout_views[ii] = inout_fa[ii]->view();
        bboxes[ii]      = bbox_a[ii];          // world-frame bbox
        transforms[ii]  = transform_a[ii];     // body→world transform
      }

      // --- Step 1: solid markers (comp 0) over grow(bx, NGHOST) ----------
      // Query points are in world frame: fast-reject vs world bbox, then
      // inverse-transform to body frame for the BVH inside/outside test.
      const Box& bxg = amrex::grow(bx, cls_t::NGHOST);
      amrex::ParallelFor(bxg,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        ibMarkers(i, j, k, 0) = uint8_t(0);
        ibMarkers(i, j, k, 1) = uint8_t(0);

        Point gridpoint = make_grid_point(prob_lo, dx_lev, i, j, k);
        for (int ii = 0; ii < ngeom_local; ++ii) {
          if (!bbox_contains(bboxes[ii], gridpoint)) continue;
          Point gp_body = transforms[ii].to_body(gridpoint);
          BoundedSide result = inout_views[ii](gp_body);
          if (result == BoundedSide::Inside || result == BoundedSide::OnBoundary) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      });

      // --- Step 2: ghost markers (comp 1) + count ------------------------
      // A solid cell is a ghost point iff at least one neighbour within
      // ±ghost_layers along an axis is fluid. Comp 1 stores the geometry id
      // (== comp 0) on ghost points; non-ghost cells keep the comp-1 zero
      // written by Step 1 (relies on GP_BOX_EXTRA <= NGHOST static_assert).
      // Count via per-FAB device atomic slot instead of a fused ReduceOps:
      // ReduceData::value() forced a blocking D2H + stream sync per FAB;
      // the slots are read back once for all FABs after the loop.
      {
        constexpr int gl   = ghost_layers;
        const Box&    bxgp = amrex::grow(bx, GP_BOX_EXTRA);
        int* const p_count_fab = p_gp_counts + ifab_local;

        amrex::ParallelFor(bxgp,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (!ibMarkers(i, j, k, 0)) return;  // fluid: not a candidate

          bool ghost = false;
          for (int d = 1; d <= gl; ++d) {
            ghost = ghost || (!ibMarkers(i-d, j,   k,   0))
                          || (!ibMarkers(i+d, j,   k,   0))
                          || (!ibMarkers(i,   j-d, k,   0))
                          || (!ibMarkers(i,   j+d, k,   0));
#if (AMREX_SPACEDIM == 3)
            ghost = ghost || (!ibMarkers(i,   j,   k-d, 0))
                          || (!ibMarkers(i,   j,   k+d, 0));
#endif
            if (ghost) break;
          }

          if (ghost) {
            ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
            Gpu::Atomic::Add(p_count_fab, 1);
          }
        });
      }

#else
      // ===================================================================
      // CPU path — semantically identical to GPU. CPU additionally emits
      // an OnBoundary diagnostic via IB_WarnOnBoundary (cannot print on GPU).
      // ===================================================================

      // --- Step 1: solid markers (comp 0) over grow(bx, NGHOST) ----------
      amrex::LoopOnCpu(amrex::grow(bx, cls_t::NGHOST),
        [&](int i, int j, int k)
      {
        ibMarkers(i, j, k, 0) = uint8_t(0);
        ibMarkers(i, j, k, 1) = uint8_t(0);

        Point gridpoint = make_grid_point(prob_lo, dx_lev, i, j, k);
        for (int ii = 0; ii < ngeom; ++ii) {
          if (!bbox_contains(bbox_a[ii], gridpoint)) continue;
          Point gp_body = transform_a[ii].to_body(gridpoint);
          inside_t& inside = *inout_fa[ii];
          BoundedSide result = inside(gp_body);
          IB_WarnOnBoundary(ii, lev, i, j, k, result, gridpoint);
          if (result == BoundedSide::Inside || result == BoundedSide::OnBoundary) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      });

      // --- Step 2: ghost markers (comp 1) + count ------------------------
      int ngps_fab = 0;
      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA),
        [&](int i, int j, int k)
      {
        if (!ibMarkers(i, j, k, 0)) return;  // fluid: not a candidate

        bool ghost = false;
        for (int d = 1; d <= ghost_layers; ++d) {
          ghost = ghost || (!ibMarkers(i-d, j,   k,   0))
                        || (!ibMarkers(i+d, j,   k,   0))
                        || (!ibMarkers(i,   j-d, k,   0))
                        || (!ibMarkers(i,   j+d, k,   0));
#if (AMREX_SPACEDIM == 3)
          ghost = ghost || (!ibMarkers(i,   j,   k-d, 0))
                        || (!ibMarkers(i,   j,   k+d, 0));
#endif
          if (ghost) break;
        }

        if (ghost) {
          ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
          ++ngps_fab;
        }
        // Non-ghost solid cells keep the comp-1 zero from Step 1 init
        // (relies on GP_BOX_EXTRA <= NGHOST static_assert).
      });
      gp_counts[ifab_local] = ngps_fab;

#endif // AMREX_USE_GPU

    } // end MFIter

#ifdef AMREX_USE_GPU
    // Single batched readback of all per-FAB GP counts (one sync for the
    // whole level; replaces one blocking readback per FAB).
    Gpu::copy(Gpu::deviceToHost, d_gp_counts.begin(), d_gp_counts.end(),
              gp_counts.begin());
#endif

    // ---- Build the level-wide flattened GPStore (CSR layout) -------------
    gpstore_a[lev].allocate(gp_counts);
  }

  /**
   * \brief Initialises geometric and interpolation data for Ghost Points (GPs).
   *
   * This function iterates over all ghost points identified in the `computeMarkers` step.
   * For each ghost point, it performs the following operations:
   * 1. Identifies the specific geometry (body) the ghost point belongs to.
   * 2. Finds the closest point on the surface (IB point) using BVH trees.
   * 3. Computes the normal distance from the ghost point to the surface.
   * 4. Constructs a local orthonormal frame (normal, tangent1, tangent2) at the IB point.
   * 5. Projects "Image Points" (IMs) into the fluid domain along the surface normal.
   * 6. Computes trilinear interpolation weights and indices for these Image Points.
   *
   * All computed data is stored in the level-wide `GPStore` (CSR layout) for use in boundary condition reconstruction.
   *
   * \param lev The current AMR level index.
   */
  void initialiseGPs(int lev) {
    BL_PROFILE("IBM::initialiseGPs");
    auto& mfab = *bmf_a[lev];
    auto& gpstore = gpstore_a[lev];
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

#ifdef AMREX_USE_GPU
    // ================================================================
    // GPU path: ParallelFor with atomic counter for GP index
    // ================================================================

    // ------------------------------------------------------------------
    // Pack per-geometry data into trivially-copyable POD GpuArrays so the
    // device lambda can capture-by-value.  Stack-allocated, fixed length
    // MAX_NGEOM (16) — multi-body limit; raise in ibm_containers.h if hit.
    // ------------------------------------------------------------------
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
        "ngeom exceeds MAX_NGEOM; increase MAX_NGEOM in ibm_containers.h");
    const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);

    // bqv[ii] : BVH4 query view for geometry ii (POD: raw ptrs to nodes/
    //           verts/faces). Used on device for closest-point queries in
    //           body frame. Owning BVH lives in bvh_a[ii].
    GpuArray<BVH4QueryView, MAX_NGEOM> bqv;

    // geom_off[ii] : level-wide flat-id offset for geometry ii's elements.
    //                Body-local prim_id + geom_off[ii] = global element id
    //                (used to index into LocalFrame_a / SurfElem_a).
    //                geom_off[ngeom] is the sentinel = total #elements.
    GpuArray<int, MAX_NGEOM + 1> geom_off;

    // transforms[ii] : rigid-body transform (body→world) for geometry ii.
    //                  Identity for static bodies; refreshed by FSI before
    //                  this kernel via updateRigidTransform().
    GpuArray<RigidTransform, MAX_NGEOM> transforms;

    for (int ii = 0; ii < ngeom_local; ii++) {
        bqv[ii]        = bvh_a[ii].query_view(geom_a[ii]);
        geom_off[ii]   = geom_offsets[ii];
        transforms[ii] = transform_a[ii];
    }
    geom_off[ngeom_local] = geom_offsets[ngeom_local];

    // ------------------------------------------------------------------
    // Level-scalar constants captured by value into the device lambda.
    // ------------------------------------------------------------------

    // Cell sizes on this level (dx, dy[, dz]).
    const auto dx_lev       = dx_a[lev];

    // Image-point spacing along the surface normal (= alpha * cell_diagonal).
    // search_*_image_point uses multiples of this to place IPs.
    const Real di_lev       = di_a[lev];

    // Cell diagonal length on this level — used as upper-bound sanity
    // check for ghost-point ↔ closest-surface-point distance.
    const Real diag_lev     = diag_a[lev];

    // Raw pointer to level-wide LocalFrame array (per-element body-frame
    // normal + tangents). Indexed by global f_idx = prim_id + geom_off[ii].
    const auto* lf_ptr      = LocalFrame_a.data();

    // GPU-capturable POD view of the CSR GPStore. Holds both const read
    // pointers and writable _w pointers used during this init pass only.
    auto gpview = gpstore.view();

    // All host-visible counters batched into ONE device array read back with
    // a single copy+sync after the MFIter loop (previously: one blocking
    // dataValue() readback + stream sync PER FAB for the count verification,
    // plus a reset kernel per FAB):
    //   slots [0, nfabs)  — per-FAB atomic GP counters (assign GP slots);
    //                       pre-zeroed here, so no per-FAB reset is needed.
    //   slot  nfabs       — first-IP stencil-out-of-box failure count
    //   slot  nfabs+1     — first-IP below-threshold placement failure count
    //   slot  nfabs+2     — first-IP interpolation fit failure count
    // (failure counters give GPU parity with the CPU path's abort/print
    // semantics: device printf is lossy and asserts are stripped in release,
    // so failures are aggregated and reported host-side).
    const int nfabs_local = mfab.local_size();

    // Host mirror of the managed CSR offsets. fab_offsets is managed memory,
    // and on devices with concurrentManagedAccess==0 (WSL) ANY host touch of
    // managed memory while device work is in flight faults. The old per-FAB
    // blocking readbacks serialized the device before every host read by
    // accident; with batched counters the launch loops below run ahead
    // asynchronously, so every host-side offset read must go through this
    // plain host copy instead (device kernels keep the managed copy via
    // gpview). The one sync here replaces this function's former 2-per-FAB.
    Gpu::streamSynchronize();
    Vector<int> h_fab_offsets(nfabs_local + 1);
    std::copy(gpstore.fab_offsets.begin(), gpstore.fab_offsets.end(),
              h_fab_offsets.begin());

    Vector<int> h_counts(nfabs_local + 3, 0);
    Gpu::DeviceVector<int> d_counts(nfabs_local + 3);
    Gpu::copyAsync(Gpu::hostToDevice, h_counts.begin(), h_counts.end(),
                   d_counts.begin());
    int* p_gp_count    = d_counts.data();                    // per-FAB slots
    int* p_err_stencil = d_counts.data() + nfabs_local;      // first-IP stencil out of box
    int* p_err_place   = d_counts.data() + nfabs_local + 1;  // first-IP below threshold
    int* p_err_fit     = d_counts.data() + nfabs_local + 2;  // first-IP LS/weight fit failure

    int ifab_local = 0;
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      const int gp_offset = h_fab_offsets[ifab_local];
      const int ngps_fab  = h_fab_offsets[ifab_local + 1] - gp_offset;

      if (ngps_fab == 0) continue;

      // bxg = valid box + NGHOST ghost layers. Used as the bounds container
      // for image-point stencil checks: the 2x2 (2D) or 2x2x2 (3D) bilinear
      // stencil corner cells must lie inside bxg, so that their markers
      // (written by computeMarkers in this same range) are valid reads.
      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);

      // This FAB's dedicated counter slot (pre-zeroed at allocation above,
      // so the old per-FAB single-thread reset kernel is gone).
      int* const p_fab_count = p_gp_count + ifab_local;

      const Box bxgp = amrex::grow(bx, GP_BOX_EXTRA);
      amrex::ParallelFor(bxgp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (!ibMarkers(i, j, k, 1)) return;

        // Atomically claim a slot in the flat GP array.
        // Note: local_idx assignment order is non-deterministic on GPU —
        // the GP at the same physical cell (i,j,k) may receive different
        // gidx across runs. This does not affect physics: downstream reads
        // gpview.gp_ijk[ii] to recover (i,j,k), and final ghost-cell
        // primitives are written back to those cells regardless of how
        // gidx is permuted.
        int local_idx = Gpu::Atomic::Add(p_fab_count, 1);
        const int gidx = gp_offset + local_idx;

        Point gp = make_grid_point(prob_lo, dx_lev, i, j, k);
        gpview.gp_ijk_w[gidx] = make_vec<int>(i, j, k);  // dim-aware: 2D ignores k

        int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;
        AMREX_ASSERT(geomIdx >= 0 && geomIdx < ngeom_local);
        gpview.geomIdx_w[gidx] = geomIdx;

        // Transform query point to body frame, then BVH closest-point query
        const auto& T = transforms[geomIdx];
        Point gp_body = T.to_body(gp);
        ClosestPointResult closest_elem = bqv[geomIdx].closest_point_query(gp_body);
        int f_idx = closest_elem.prim_id + geom_off[geomIdx];
        gpview.elemIdx_w[gidx] = f_idx;

        // Transform closest point back to world frame
        Point cp = T.to_world(closest_elem.point);

        // Distance (invariant under rigid transform, but computed in world frame)
        Real disGP_val = std::sqrt(point_distance_sq(gp, cp));

        // Sanity: ghost points are solid cells with at least one fluid
        // neighbour within ghost_layers, so the closest surface point must
        // be within roughly one cell diagonal. Failure here means one of:
        //   (a) BVH returned a wrong primitive (geometry data corruption)
        //   (b) marker comp 1 was set on a cell not adjacent to fluid
        //   (c) FAB layout / DistributionMapping inconsistency between
        //       computeMarkers and initialiseGPs (shouldn't happen, but
        //       would manifest here first if a regrid slipped between).
        AMREX_ASSERT(disGP_val < diag_lev);
        gpview.disGP_w[gidx] = disGP_val;
        gpview.ib_xyz_w[gidx] = make_vec<Real>(cp);

        // Rotate body-frame LocalFrame to world frame
        const LocalFrame& lf_body = lf_ptr[f_idx];
        LocalFrame localframe;
        T.rotate_to_world(lf_body.normal,   localframe.normal);
        T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
        T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif

        // Image points: stack-local SoA, written into GPStore once at the end.
        // Value-initialisation keeps diagnostics deterministic before the
        // aggregated host-side failure check aborts an invalid first IP.
        Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz{};
        Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk{};
        Array1D<Real, 0, eorder_tparm - 1> disIM{};
        Array1D< int, 0, eorder_tparm - 1> imp_ninterp{};

        // Chained walk along the outward normal — jj=0 from IB point,
        // jj>=1 from previous IP. Shared with CPU path.
        const int place_status = place_image_points<eorder_tparm, iorder_tparm>(
            cp, localframe,
            lev, prob_lo, dx_lev, di_lev,
            bxg, ibMarkers,
            gpview, gidx,
            imp_xyz, imp_ijk, disIM, imp_ninterp);
        if (place_status & 1) { Gpu::Atomic::Add(p_err_stencil, 1); }
        if (place_status & 2) { Gpu::Atomic::Add(p_err_place, 1); }

        gpview.imp_xyz_w[gidx]       = imp_xyz;
        gpview.imp_ijk_w[gidx]       = imp_ijk;
        gpview.disIM_w[gidx]         = disIM;

        Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1> imp_ipweights;
        Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;

        computeIPweights<eorder_tparm, iorder_tparm, GPSTOREVIEW>(
            imp_ipweights, imp_ip_ijk,
            imp_xyz, imp_ijk, imp_ninterp,
            prob_lo, dx_lev, ibMarkers);
        if (place_status == 0 && imp_ninterp(0) == 0) {
          Gpu::Atomic::Add(p_err_fit, 1);
        }

        // Store imp_ninterp AFTER computeIPweights: the iorder=2 WLS path
        // demotes an unfittable image point by setting imp_ninterp(iim)=0,
        // and n_valid/eff_order downstream must see that demotion.
        gpview.imp_ninterp_w[gidx]   = imp_ninterp;
        gpview.imp_ipweights_w[gidx] = imp_ipweights;
        gpview.imp_ip_ijk_w[gidx]    = imp_ip_ijk;

        ibMarkers(i, j, k, 1) = static_cast<uint8_t>(imp_ninterp(0));
      }); // end ParallelFor

    } // end MFIter

    // Single batched readback of per-FAB counts + error counters (one sync
    // for the whole level). Verifying after all FABs launched instead of
    // per-FAB does not add risk: any CSR overflow from a miscounted FAB
    // already happened inside that FAB's own kernel (the atomic claims
    // slots during the launch, before any readback could run), and in both
    // versions we abort before the GP data is consumed.
    Gpu::copy(Gpu::deviceToHost, d_counts.begin(), d_counts.end(),
              h_counts.begin());
    for (int f = 0; f < nfabs_local; ++f) {
      const int gp_offset = h_fab_offsets[f];
      const int ngps_fab  = h_fab_offsets[f + 1] - gp_offset;
      // Fast-skipped FABs (ngps_fab==0) never touch their pre-zeroed slot.
      if (h_counts[f] != ngps_fab) {
        amrex::Print() << "Error in initialiseGPs (GPU): GP count mismatch\n"
                       << "  level=" << lev
                       << "  FAB=" << f
                       << "  expected=" << ngps_fab
                       << "  got=" << h_counts[f]
                       << "  gp_offset=" << gp_offset << "\n";
        amrex::Abort("initialiseGPs: ghost point count mismatch");
      }
    }

    // Aggregated failure report (parity with the CPU path, which aborts on a
    // first-IP stencil escape and prints per-GP placement failures).
    // Counters were read back in the batched copy above — no extra sync.
    {
      const int n_stencil = h_counts[nfabs_local];
      const int n_place   = h_counts[nfabs_local + 1];
      const int n_fit     = h_counts[nfabs_local + 2];
      if (n_stencil > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_stencil
                       << " ghost point(s) on level " << lev
                       << " have their first image-point stencil outside the grown box.\n";
        amrex::Abort("initialiseGPs: interpolation stencil out of box (GPU)");
      }
      if (n_place > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_place
                       << " ghost point(s) on level " << lev
                       << " failed first image-point placement (below threshold).\n";
        amrex::Abort("initialiseGPs: invalid first image point (GPU)");
      }
      if (n_fit > 0) {
        amrex::Print() << "initialiseGPs (GPU): " << n_fit
                       << " ghost point(s) on level " << lev
                       << " have no resolvable interpolation fit at the first image point.\n";
        amrex::Abort("initialiseGPs: invalid first image-point fit (GPU)");
      }
    }

#else
    // ================================================================
    // CPU path: original LoopOnCpu implementation
    // ================================================================
    int ifab_local = 0;
    int n_stencil_fail = 0;
    int n_place_fail = 0;
    int n_fit_fail = 0;
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      // Get the pre-allocated range in the flattened GPStore for this FAB.
      // (CPU build: ManagedVector is plain host memory — direct read is safe;
      //  the h_fab_offsets host mirror exists only in the GPU branch above.)
      const int gp_offset = gpstore.fab_offsets[ifab_local];
      const int ngps_fab  = gpstore.fab_offsets[ifab_local + 1] - gp_offset;

      if (ngps_fab == 0) continue;

      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);  // uint8_t array

      // CPU loop: Write directly into gpstore at pre-allocated offsets (no push_back).
      int gp_count = 0; 
      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA), [&](int i, int j, int k) {
        // for each ghost point
        if (ibMarkers(i, j, k, 1)) {
          const int gidx = gp_offset + gp_count;  // global flat index in GPStore
          gp_count++;

          Point gp = make_grid_point(prob_lo, dx_a[lev], i, j, k);
          gpstore.gp_ijk[gidx] = make_vec<int>(i, j, k);
          
          int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;
          AMREX_ASSERT_WITH_MESSAGE(geomIdx >= 0 && geomIdx < ngeom, 
                  "Invalid geometry index in initialiseGPs");
          gpstore.geomIdx[gidx] = geomIdx;

          // Rigid-body transform: geometry trees/frames live in BODY frame;
          // query in body frame, return results to world frame (mirrors the
          // GPU path — previously missing here, which made CPU builds use
          // wrong closest points/normals once an FSI body moved).
          const auto& T = transform_a[geomIdx];
          Point gp_body = T.to_body(gp);

#ifdef AMREX_USE_CGAL
          ClosestPointResult closest_elem =
              cgal_closest_point_query(tree_a[geomIdx], idxmap_a[geomIdx],
                                       gp_body, this->geom_offsets[geomIdx]);
#else
          ClosestPointResult closest_elem =
              bvh_a[geomIdx].closest_point_query(gp_body, geom_a[geomIdx]);
#endif

          int f_idx = closest_elem.prim_id + this->geom_offsets[geomIdx];
          gpstore.elemIdx[gidx] = f_idx;

          Point cp = T.to_world(closest_elem.point);

          Real disGP = std::sqrt(point_distance_sq(gp, cp));

          AMREX_ASSERT_WITH_MESSAGE(
              disGP < diag_a[lev],  "Ghost point and IB point distance larger than mesh diagonal");
          gpstore.disGP[gidx] = disGP;

          gpstore.ib_xyz[gidx] = make_vec<Real>(cp);

          const LocalFrame& lf_body = LocalFrame_a[f_idx];
          LocalFrame localframe;
          T.rotate_to_world(lf_body.normal,   localframe.normal);
          T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
          T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif

          // Value-initialised ({}) so a failed search remains deterministic
          // until the aggregated host-side failure check below aborts.
          Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz{};
          Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk{};
          Array1D<Real, 0, eorder_tparm - 1> disIM{};
          Array1D< int, 0, eorder_tparm - 1> imp_ninterp{};

          // Chained walk along the outward normal — shared with GPU path.
          const int place_status = place_image_points<eorder_tparm, iorder_tparm>(
              cp, localframe,
              lev, prob_lo, dx_a[lev], di_a[lev],
              bxg, ibMarkers,
              gpstore, gidx,
              imp_xyz, imp_ijk, disIM, imp_ninterp);
          if (place_status & 1) { ++n_stencil_fail; }
          if (place_status & 2) { ++n_place_fail; }

          gpstore.imp_xyz[gidx] = imp_xyz;
          gpstore.imp_ijk[gidx] = imp_ijk;
          gpstore.disIM[gidx] = disIM;

          Array2D<Real, 0, eorder_tparm - 1 , 0, N_InterP -1 > imp_ipweights;
          Array3D< int, 0, eorder_tparm - 1 , 0, N_InterP -1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;
          
          computeIPweights<eorder_tparm, iorder_tparm, GPSTORE>(
              imp_ipweights, 
              imp_ip_ijk, 
              imp_xyz, 
              imp_ijk, 
              imp_ninterp,
              prob_lo, dx_a[lev], ibMarkers);
          if (place_status == 0 && imp_ninterp(0) == 0) {
            ++n_fit_fail;
          }
          
          // Store imp_ninterp AFTER computeIPweights (WLS demotion visibility).
          gpstore.imp_ninterp[gidx] = imp_ninterp;
          gpstore.imp_ipweights[gidx] = imp_ipweights;
          gpstore.imp_ip_ijk[gidx] = imp_ip_ijk;

          ibMarkers(i, j, k, 1) = static_cast<uint8_t>(imp_ninterp(0)); 

          } //end if (ibMarkers(i,j,k,1))
      });//end loop on bx
    
      if(gp_count != ngps_fab) {
        amrex::Abort("Error in initialiseGPs: mismatch in ghost point count");
      }
    } //end MFIter

    if (n_stencil_fail > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_stencil_fail
                     << " ghost point(s) on level " << lev
                     << " have their first image-point stencil outside the grown box.\n";
      amrex::Abort("initialiseGPs: interpolation stencil out of box (CPU)");
    }
    if (n_place_fail > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_place_fail
                     << " ghost point(s) on level " << lev
                     << " failed first image-point placement (below threshold).\n";
      amrex::Abort("initialiseGPs: invalid first image point (CPU)");
    }
    if (n_fit_fail > 0) {
      amrex::Print() << "initialiseGPs (CPU): " << n_fit_fail
                     << " ghost point(s) on level " << lev
                     << " have no resolvable interpolation fit at the first image point.\n";
      amrex::Abort("initialiseGPs: invalid first image-point fit (CPU)");
    }

#endif // AMREX_USE_GPU

    // initialiseGPs rewrites marker component 1 from the temporary geometry id
    // to the final first-IP interpolation count.  Only valid cells own that
    // value; same-level ghost copies must be refreshed before any flux stencil
    // inspects them.  Exchange both components so the solid and GP views have
    // one authoritative, box-layout-independent representation.
    mfab.FillBoundary(0, 2, amr_p->Geom(lev).periodicity());

    // NOTE: gpstore.shrink() disabled — shrink_to_fit() on PODVector
    // corrupts data on some platforms (observed with AMReX ManagedVector on CPU).

    bool ib_diag = false;
    bool gp_diag = false;
    {
      ParmParse pp("ib");
      pp.query("diag", ib_diag);
      gp_diag = ib_diag;
      pp.query("gp_diag", gp_diag);
    }
    if (gp_diag) reportGPDiagnostics(lev);
  }

  // ========================================================================
  // GP reconstruction quality diagnostics
  // ========================================================================

  /**
   * \brief Report statistics on ghost-point reconstruction quality.
   *
   * Scans the GPStore after initialiseGPs and prints:
   *   - Total GP count
   *   - Histogram of first image-point fluid stencil count (imp_ninterp[0])
   *   - Effective extrapolation order distribution (0th / 1st / 2nd)
   *   - disGP statistics (min, mean, max) in units of dx diagonal
   *   - Weight concentration: max single-point weight among all GPs
   *
   * \param lev  AMR level index.
   */
  void reportGPDiagnostics(int lev)
  {
    auto& gpstore = gpstore_a[lev];
    const int ngps = gpstore.total_ngps;
    const int istep = amr_p->levelSteps(0);

    // Ensure GPU data is accessible on host
    Gpu::streamSynchronize();

    constexpr int IDEAL = N_InterP;
    int hist[IDEAL + 1] = {};

    // Cross-tabulation: severity × effective order
    // severity bins: ideal (n_fluid==IDEAL), mild (IDEAL/2 < n < IDEAL),
    //               severe (n <= IDEAL/2)
    int sev_ideal = 0, sev_mild = 0, sev_severe = 0;
    int sev_mild_order[3] = {}, sev_severe_order[3] = {};  // indexed by eff_order

    int order_counts[3] = {};

    Real disGP_min = std::numeric_limits<Real>::max();
    Real disGP_max = Real(0.0);
    Real disGP_sum = Real(0.0);

    Real max_single_weight = Real(0.0);
    Real sum_weight_deficit = Real(0.0);
    int  n_low_stencil = 0;

    for (int ii = 0; ii < ngps; ++ii) {
      int n0 = gpstore.imp_ninterp[ii](0);
      int idx = amrex::max(0, amrex::min(IDEAL, n0));
      hist[idx]++;

      // Effective order
      int eff = eorder_tparm;
      for (int k = 0; k < eorder_tparm; ++k) {
        if (gpstore.imp_ninterp[ii](k) < INTERP_THRESHOLD_GP) {
          eff = k; break;
        }
      }
      order_counts[amrex::min(eff, 2)]++;

      // Severity classification
      if (n0 == IDEAL) {
        sev_ideal++;
      } else if (n0 > IDEAL / 2) {
        sev_mild++;
        sev_mild_order[amrex::min(eff, 2)]++;
      } else {
        sev_severe++;
        sev_severe_order[amrex::min(eff, 2)]++;
      }

      // disGP
      Real d = gpstore.disGP[ii];
      disGP_min = amrex::min(disGP_min, d);
      disGP_max = amrex::max(disGP_max, d);
      disGP_sum += d;

      // Weight concentration
      Real wmax = Real(0.0);
      for (int c = 0; c < N_InterP; ++c) {
        wmax = amrex::max(wmax, gpstore.imp_ipweights[ii](0, c));
      }
      max_single_weight = amrex::max(max_single_weight, wmax);

      if (n0 < IDEAL) {
        n_low_stencil++;
        sum_weight_deficit += wmax;
      }
    }

    int ngps_global = ngps;
    int sev_counts[3] = {sev_ideal, sev_mild, sev_severe};
    ParallelDescriptor::ReduceIntSum(ngps_global);
    ParallelDescriptor::ReduceIntSum(hist, IDEAL + 1);
    ParallelDescriptor::ReduceIntSum(sev_counts, 3);
    ParallelDescriptor::ReduceIntSum(sev_mild_order, 3);
    ParallelDescriptor::ReduceIntSum(sev_severe_order, 3);
    ParallelDescriptor::ReduceIntSum(order_counts, 3);
    ParallelDescriptor::ReduceIntSum(n_low_stencil);

    ParallelDescriptor::ReduceRealMin(disGP_min);
    ParallelDescriptor::ReduceRealMax(disGP_max);
    ParallelDescriptor::ReduceRealSum(disGP_sum);
    ParallelDescriptor::ReduceRealMax(max_single_weight);
    ParallelDescriptor::ReduceRealSum(sum_weight_deficit);

    sev_ideal  = sev_counts[0];
    sev_mild   = sev_counts[1];
    sev_severe = sev_counts[2];

    if (ngps_global == 0) {
      amrex::Print() << "[GP-Diag] Step " << istep << " Level " << lev << ": 0 ghost points\n";
      return;
    }

    Real diag = diag_a[lev];
    Real inv_diag = (diag > Real(1e-30)) ? Real(1.0) / diag : Real(1.0);

    // --- Console output ---
    amrex::Print()
      << "\n[GP-Diag] Step " << istep << " Level " << lev
      << "  |  eorder=" << eorder_tparm
      << "  iorder=" << iorder_tparm
      << "  ghost_layers=" << ghost_layers
      << "  alpha=" << cim
      << "\n";

    amrex::Print()
      << "[GP-Diag]   Total GPs: " << ngps_global << "\n";

    // Fluid stencil histogram (layered)
    amrex::Print() << "[GP-Diag]   Stencil quality (ideal=" << IDEAL << "):\n";
    amrex::Print() << "[GP-Diag]     IDEAL  (" << IDEAL << "/" << IDEAL << "): "
                   << sev_ideal << " GPs (" << std::fixed << std::setprecision(1)
                   << Real(100.0) * Real(sev_ideal) / Real(ngps_global) << "%)\n";
    if (sev_mild > 0) {
      amrex::Print() << "[GP-Diag]     MILD   (>" << IDEAL/2 << "/" << IDEAL << "): "
                     << sev_mild << " GPs (" << std::setprecision(1)
                     << Real(100.0) * Real(sev_mild) / Real(ngps_global) << "%)";
      // Detail per n_fluid
      for (int k = IDEAL - 1; k > IDEAL / 2; --k) {
        if (hist[k] > 0) amrex::Print() << "  [" << k << "/" << IDEAL << "]=" << hist[k];
      }
      amrex::Print() << "\n";
    }
    if (sev_severe > 0) {
      amrex::Print() << "[GP-Diag]     SEVERE (<=" << IDEAL/2 << "/" << IDEAL << "): "
                     << sev_severe << " GPs (" << std::setprecision(1)
                     << Real(100.0) * Real(sev_severe) / Real(ngps_global) << "%)";
      for (int k = IDEAL / 2; k >= 0; --k) {
        if (hist[k] > 0) amrex::Print() << "  [" << k << "/" << IDEAL << "]=" << hist[k];
      }
      amrex::Print() << "\n";
    }

    // Cross-tabulation: severity × order
    if (sev_mild + sev_severe > 0) {
      amrex::Print() << "[GP-Diag]   Severity x Order cross-tab:\n";
      if (sev_mild > 0) {
        amrex::Print() << "[GP-Diag]     MILD  : ";
        if (sev_mild_order[0]) amrex::Print() << "0th=" << sev_mild_order[0] << " ";
        if (sev_mild_order[1]) amrex::Print() << "1st=" << sev_mild_order[1] << " ";
        if (sev_mild_order[2]) amrex::Print() << "2nd=" << sev_mild_order[2] << " ";
        amrex::Print() << "\n";
      }
      if (sev_severe > 0) {
        amrex::Print() << "[GP-Diag]     SEVERE: ";
        if (sev_severe_order[0]) amrex::Print() << "0th=" << sev_severe_order[0] << " ";
        if (sev_severe_order[1]) amrex::Print() << "1st=" << sev_severe_order[1] << " ";
        if (sev_severe_order[2]) amrex::Print() << "2nd=" << sev_severe_order[2] << " ";
        amrex::Print() << "\n";
      }
      // Weight stats
      Real avg_wmax = (n_low_stencil > 0)
        ? sum_weight_deficit / Real(n_low_stencil) : Real(0.0);
      amrex::Print()
        << "[GP-Diag]   Weight: avg_max=" << std::setprecision(3) << avg_wmax
        << "  worst_max=" << std::setprecision(3) << max_single_weight
        << "\n";
    }

    // Effective order distribution
    amrex::Print() << "[GP-Diag]   Effective order:";
    if (order_counts[0] > 0) amrex::Print() << "  0th=" << order_counts[0];
    if (order_counts[1] > 0) amrex::Print() << "  1st=" << order_counts[1];
    if (order_counts[2] > 0) amrex::Print() << "  2nd+=" << order_counts[2];
    amrex::Print() << "\n";

    // disGP statistics
    amrex::Print()
      << "[GP-Diag]   disGP (diag units): min="
      << std::setprecision(4) << disGP_min * inv_diag
      << "  mean=" << std::setprecision(4) << (disGP_sum / Real(ngps_global)) * inv_diag
      << "  max=" << std::setprecision(4) << disGP_max * inv_diag
      << "  (diag=" << std::setprecision(6) << diag << ")\n";

    amrex::Print() << "\n";
  }

  /**
   * \brief Report statistics on surface image-point reconstruction quality.
   *
   * This is the surface counterpart of reportGPDiagnostics(). It measures how
   * many surface elements actually have enough leading valid image points to
   * use the requested extrap_order_surf branch.
   */
  void reportSurfDiagnostics(int lev)
  {
    const int istep = amr_p->levelSteps(0);
    Gpu::streamSynchronize();

    constexpr int IDEAL = N_InterP_surf;
    int hist[IDEAL + 1] = {};
    int sev_counts[3] = {};      // ideal, mild, severe
    int order_counts[3] = {};    // 0th, 1st, 2nd+
    int surf_faces = 0;
    int low_stencil_faces = 0;

    Real dis1_min = std::numeric_limits<Real>::max();
    Real dis1_max = Real(0.0);
    Real dis1_sum = Real(0.0);

    int ratio_count = 0;
    Real ratio_min = std::numeric_limits<Real>::max();
    Real ratio_max = Real(0.0);
    Real ratio_sum = Real(0.0);

    for (int f = 0; f < ntotalfaces; ++f) {
      if (surfphys_soa.elemfound[f] != 1) continue;
      if (surfphys_soa.lev[f] != lev) continue;

      surf_faces++;

      const int n0 = surfimp_soa.imp_ninterp[f](0);
      const int idx = amrex::max(0, amrex::min(IDEAL, n0));
      hist[idx]++;

      int eff = eorder_tparm_surf;
      for (int k = 0; k < eorder_tparm_surf; ++k) {
        if (surfimp_soa.imp_ninterp[f](k) < INTERP_THRESHOLD_SURF) {
          eff = k;
          break;
        }
      }
      order_counts[amrex::min(eff, 2)]++;

      if (n0 == IDEAL) {
        sev_counts[0]++;
      } else if (n0 > IDEAL / 2) {
        sev_counts[1]++;
      } else {
        sev_counts[2]++;
      }
      if (n0 < IDEAL) low_stencil_faces++;

      const Real d1 = surfimp_soa.disIM[f](0);
      dis1_min = amrex::min(dis1_min, d1);
      dis1_max = amrex::max(dis1_max, d1);
      dis1_sum += d1;

      if constexpr (eorder_tparm_surf >= 2) {
        if (eff >= 2 && d1 > Real(0.0)) {
          const Real r = surfimp_soa.disIM[f](1) / d1;
          ratio_min = amrex::min(ratio_min, r);
          ratio_max = amrex::max(ratio_max, r);
          ratio_sum += r;
          ratio_count++;
        }
      }
    }

    ParallelDescriptor::ReduceIntSum(surf_faces);
    ParallelDescriptor::ReduceIntSum(low_stencil_faces);
    ParallelDescriptor::ReduceIntSum(hist, IDEAL + 1);
    ParallelDescriptor::ReduceIntSum(sev_counts, 3);
    ParallelDescriptor::ReduceIntSum(order_counts, 3);
    ParallelDescriptor::ReduceIntSum(ratio_count);

    ParallelDescriptor::ReduceRealMin(dis1_min);
    ParallelDescriptor::ReduceRealMax(dis1_max);
    ParallelDescriptor::ReduceRealSum(dis1_sum);
    ParallelDescriptor::ReduceRealMin(ratio_min);
    ParallelDescriptor::ReduceRealMax(ratio_max);
    ParallelDescriptor::ReduceRealSum(ratio_sum);

    if (surf_faces == 0) {
      amrex::Print() << "[Surf-Diag] Step " << istep << " Level " << lev
                     << ": 0 owned surface faces\n";
      return;
    }

    const Real diag = diag_a[lev];
    const Real inv_diag = (diag > Real(1e-30)) ? Real(1.0) / diag : Real(1.0);

    amrex::Print()
      << "\n[Surf-Diag] Step " << istep << " Level " << lev
      << "  |  eorder_surf=" << eorder_tparm_surf
      << "  iorder_surf=" << iorder_tparm_surf
      << "  alpha_surf=" << cim_surf
      << "\n";
    amrex::Print() << "[Surf-Diag]   Owned faces: " << surf_faces << "\n";

    amrex::Print() << "[Surf-Diag]   Stencil quality (ideal=" << IDEAL << "):\n";
    amrex::Print() << "[Surf-Diag]     IDEAL  (" << IDEAL << "/" << IDEAL << "): "
                   << sev_counts[0] << " faces (" << std::fixed << std::setprecision(1)
                   << Real(100.0) * Real(sev_counts[0]) / Real(surf_faces) << "%)\n";
    if (sev_counts[1] > 0) {
      amrex::Print() << "[Surf-Diag]     MILD   (>" << IDEAL / 2 << "/" << IDEAL << "): "
                     << sev_counts[1] << " faces (" << std::setprecision(1)
                     << Real(100.0) * Real(sev_counts[1]) / Real(surf_faces) << "%)";
      for (int k = IDEAL - 1; k > IDEAL / 2; --k) {
        if (hist[k] > 0) amrex::Print() << "  [" << k << "/" << IDEAL << "]=" << hist[k];
      }
      amrex::Print() << "\n";
    }
    if (sev_counts[2] > 0) {
      amrex::Print() << "[Surf-Diag]     SEVERE (<=" << IDEAL / 2 << "/" << IDEAL << "): "
                     << sev_counts[2] << " faces (" << std::setprecision(1)
                     << Real(100.0) * Real(sev_counts[2]) / Real(surf_faces) << "%)";
      for (int k = IDEAL / 2; k >= 0; --k) {
        if (hist[k] > 0) amrex::Print() << "  [" << k << "/" << IDEAL << "]=" << hist[k];
      }
      amrex::Print() << "\n";
    }

    amrex::Print() << "[Surf-Diag]   Effective order:";
    if (order_counts[0] > 0) amrex::Print() << "  0th=" << order_counts[0];
    if (order_counts[1] > 0) amrex::Print() << "  1st=" << order_counts[1];
    if (order_counts[2] > 0) amrex::Print() << "  2nd+=" << order_counts[2];
    amrex::Print() << "\n";

    if constexpr (eorder_tparm_surf >= 2) {
      amrex::Print() << "[Surf-Diag]   2nd-order eligible: " << order_counts[2]
                     << " faces (" << std::fixed << std::setprecision(1)
                     << Real(100.0) * Real(order_counts[2]) / Real(surf_faces) << "%)\n";
    }

    amrex::Print()
      << "[Surf-Diag]   IP1 distance (diag units): min="
      << std::setprecision(4) << dis1_min * inv_diag
      << "  mean=" << std::setprecision(4) << (dis1_sum / Real(surf_faces)) * inv_diag
      << "  max=" << std::setprecision(4) << dis1_max * inv_diag
      << "  (diag=" << std::setprecision(6) << diag << ")\n";

    if constexpr (eorder_tparm_surf >= 2) {
      if (ratio_count > 0) {
        amrex::Print()
          << "[Surf-Diag]   IP spacing ratio x2/x1 on 2nd-order faces: min="
          << std::setprecision(4) << ratio_min
          << "  mean=" << std::setprecision(4) << (ratio_sum / Real(ratio_count))
          << "  max=" << std::setprecision(4) << ratio_max << "\n";
      }
    }

    amrex::Print() << "\n";
  }

  // ========================================================================
  // Ghost‐point reconstruction — single‐kernel, level‐wide
  // ========================================================================

  /**
   * \brief Level-wide single-kernel ghost-point reconstruction.
   *
   * Prerequisite: \c prims_mf must already be filled with primitives
   * (via cons2prims) for every local FAB on this rank.
   *
   * Builds device-accessible Array4 pointers (one per local FAB) and
   * launches a single ParallelFor over all ghost points at this level.
   *
   * \param prims_mf  MultiFab of primitive variables (same BA/DM as bmf_a[lev]).
   *                  Ghost cells will be overwritten with IB-reconstructed values.
   * \param cls       Pointer to the physics/closure class (device-accessible).
   * \param lev       Current AMR level index.
   */
  void computeAllGPs(MultiFab& prims_mf,
                     const cls_t* cls,
                     int lev)
  {
    BL_PROFILE("IBM::computeAllGPs");
    auto& gpstore = gpstore_a[lev];
    if (gpstore.total_ngps > 0) {

    auto gpview = gpstore.view();
    auto const* lf_ptr = LocalFrame_a.data();
    const int nfabs_local = gpview.nfabs;

    // Capture per-geometry transforms for rotating body-frame LocalFrame to world
    AMREX_ALWAYS_ASSERT(ngeom <= MAX_NGEOM);
    const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    for (int ii = 0; ii < ngeom_local; ii++) {
        transforms[ii] = transform_a[ii];
    }

    // Build device array of Array4 pointers (one per local FAB).
    // No snapshot needed: GPs write only to ghost cells (ibMarkers==1)
    // while image-point interpolation reads only from fluid cells (ibMarkers==0).
    auto& mfab = *bmf_a[lev];

    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
      }
      Gpu::copyAsync(Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      // No streamSynchronize: copyAsync and ParallelFor share the same stream.
    }

    auto* prims_arr = d_prims.data();

    // Single kernel launch over all ghost points on this level
    auto* copy = this;
    const int total_ngps = gpview.total_ngps;
    const int num_lframes = static_cast<int>(LocalFrame_a.size());

    ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE (int ii) noexcept
    {
      // Determine which FAB this GP belongs to
      int ifab = gpview.gp_fab[ii];
      auto prims = prims_arr[ifab];  // read image points & write ghost cells

      // 1) Reconstruct local orthonormal frame — rotate from body to world frame
      int elem_idx = gpview.elemIdx[ii];
      const auto& frame = lf_ptr[elem_idx];

      // Find which geometry this GP belongs to (for transform lookup)
      int geomIdx = gpview.geomIdx[ii];
      const auto& xform = transforms[geomIdx];

      int type_solid_bc = 0;

      // Rotate body-frame vectors to world frame
      amrex::Real n_w[AMREX_SPACEDIM], t1_w[AMREX_SPACEDIM], t2_w[AMREX_SPACEDIM];
      xform.rotate_to_world(frame.normal,   n_w);
      xform.rotate_to_world(frame.tangent1, t1_w);
#if (AMREX_SPACEDIM == 3)
      xform.rotate_to_world(frame.tangent2, t2_w);
#endif

      Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          nvec(d)  = n_w[d];
          t1vec(d) = t1_w[d];
#if (AMREX_SPACEDIM == 3)
          t2vec(d) = t2_w[d];
#else
          t2vec(d) = Real(0.0);
#endif
      }

      // 2) Storage for primitive variables along the normal
      Array2D<Real, 0, eorder_tparm + 1, 0, cls_t::NPRIM - 1> primsNormal;
      for (int p = 0; p <= eorder_tparm + 1; ++p) {
          for (int n = 0; n < cls_t::NPRIM; ++n) {
              primsNormal(p,n) = Real(0.0);
          }
      }

      // 3) Interpolate primitive variables at all image points
      copy->template interpolateIMs<eorder_tparm, iorder_tparm>(
          gpview.imp_ip_ijk[ii], gpview.imp_ipweights[ii], prims, primsNormal);

      // 4) Transform velocities at image points to local frame
      for (int iip = 2; iip < 2 + eorder_tparm; ++iip) {
          copy->template global2local<eorder_tparm>(iip, primsNormal, nvec, t1vec, t2vec);
      }

      // 5) Apply wall model at IB surface.
      //    n_valid_bc = number of leading valid image points. The wall model
      //    uses the 2nd-order one-sided form only when >=2 are available and
      //    extrap_order>=2; otherwise it degrades to the 1st-order single-IP
      //    form (identical to the previous behaviour at extrap_order=1).
      int n_valid_bc = 0;
      for (int kk = 0; kk < eorder_tparm; ++kk) {
          if (gpview.imp_ninterp[ii](kk) < INTERP_THRESHOLD_GP) break;
          n_valid_bc = kk + 1;
      }
      gpview.n_valid_w[ii] = n_valid_bc;
      ibm_detail::dispatch_compute_surfIB<eorder_tparm, wallmodel>(
          gpview.ib_xyz[ii], nvec, t1vec, t2vec, gpview.disIM[ii], n_valid_bc,
          primsNormal, type_solid_bc, cls);

      // 6) Extrapolate from surface/image points back to ghost point
      copy->template extrapolate<eorder_tparm>(
          primsNormal, gpview.imp_ninterp[ii], gpview.disGP[ii], gpview.disIM[ii]);

      // 7) Transform ghost-point velocity back to global coordinates
      int idx = 0;
      copy->template local2global<eorder_tparm>(idx, primsNormal, nvec, t1vec, t2vec);

      // 8) Extract primitive variables at ghost point (slot 0)
      Real P = primsNormal(0, cls_t::QPRES);
      Real T = primsNormal(0, cls_t::QT);

      Real Y[NUM_SPECIES] = { Real(0.0) };
#if NUM_SPECIES > 1
      for (int n = 0; n < NUM_SPECIES; ++n) {
          Y[n] = primsNormal(0, cls_t::QFS + n);
      }
#endif

      Real ux = primsNormal(0, cls_t::QU);
      Real uy = primsNormal(0, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
      Real uz = primsNormal(0, cls_t::QW);
#else
      Real uz = Real(0.0);
#endif

      // 9) Enforce thermodynamic consistency
      Real Q[cls_t::NPRIM];
      cls->ensurePTYfillq(P, T, Y, ux, uy, uz, Q);

      gpview.recon_prim_w[ii](0) = Q[cls_t::QRHO];
      gpview.recon_prim_w[ii](1) = ux;
      gpview.recon_prim_w[ii](2) = uy;
#if (AMREX_SPACEDIM == 3)
      gpview.recon_prim_w[ii](3) = uz;
#else
      gpview.recon_prim_w[ii](3) = Real(0.0);
#endif
      gpview.recon_prim_w[ii](4) = P;
      gpview.recon_prim_w[ii](5) = T;

      // 10) Write ghost-cell primitive variables back into prims
      int i = gpview.gp_ijk[ii](0);
      int j = gpview.gp_ijk[ii](1);
#if (AMREX_SPACEDIM == 3)
      int k = gpview.gp_ijk[ii](2);
#else
      int k = 0;
#endif

      for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(i,j,k,n) = Q[n];
      }
    }); // end ParallelFor over all ghost points
    }

    // GP reconstruction writes the valid owner cell only.  A neighboring FAB
    // can read that cell through its private ghost copy in WENO/TENO or viscous
    // stencils, so publish the reconstructed values before returning.  This is
    // intentionally unconditional: ranks with zero local GPs must still
    // participate in the same-level exchange initiated by ranks that own GPs.
    prims_mf.FillBoundary(amr_p->Geom(lev).periodicity());
  }

  /**
   * \brief Repair conservative state in cells exposed by moving geometry.
   *
   * A level-wide Jacobi front propagates data from cells that are fluid in
   * both topologies into cells that changed solid -> fluid.  Tags and state
   * ghosts are exchanged after every iteration, so the result is independent
   * of GPU scheduling, FAB decomposition and MPI rank boundaries.  Only a
   * one-cell halo is accessed; there are no unchecked expanding-ring reads.
   *
   * If an exposed component is disconnected from every old-fluid donor, the
   * state is physically undefined.  Abort instead of inventing a freestream
   * value or retaining stale solid data.
   */
  void fixExposedCells(const FabArray<BaseFab<uint8_t>>& old_markers,
                       MultiFab& state_mf,
                       int lev)
  {
    BL_PROFILE("IBM::fixExposedCells");

    auto& new_markers = *bmf_a[lev];
    AMREX_ALWAYS_ASSERT(old_markers.boxArray() == new_markers.boxArray());
    AMREX_ALWAYS_ASSERT(state_mf.boxArray() == new_markers.boxArray());
    AMREX_ALWAYS_ASSERT(state_mf.nGrowVect().allGE(IntVect(1)));

    constexpr int invalid = -1;
    constexpr int donor   = 0;
    constexpr int pending = 1;
    constexpr int fixed   = 2;

    // Two tag components are the old/new Jacobi buffers.  Ghosts start invalid
    // and become valid only through FillBoundary from an owned same-level cell.
    iMultiFab tags(new_markers.boxArray(), new_markers.DistributionMap(),
                   2, 1, MFInfo().SetArena(The_Async_Arena()));
    tags.setVal(invalid);

    Gpu::DeviceScalar<int> d_pending(0);
    int* const p_pending = d_pending.dataPtr();

    for (MFIter mfi(new_markers, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& old_mk = old_markers.const_array(mfi);
      auto const& new_mk = new_markers.const_array(mfi);
      auto const& tag = tags.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int t = invalid;
        if (new_mk(i,j,k,0) == 0) {
          if (old_mk(i,j,k,0) == 0) {
            t = donor;
          } else {
            t = pending;
            Gpu::Atomic::Add(p_pending, 1);
          }
        }
        tag(i,j,k,0) = t;
        tag(i,j,k,1) = t;
      });
    }

    int n_pending = d_pending.dataValue();
    ParallelDescriptor::ReduceIntSum(n_pending);
    if (n_pending == 0) return;

    const auto& periodicity = amr_p->Geom(lev).periodicity();
    state_mf.FillBoundary(periodicity);
    tags.FillBoundary(periodicity);

    int max_iters = 0;
    const Box& domain = amr_p->Geom(lev).Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      max_iters = amrex::max(max_iters, domain.length(d));
    }

    int previous_pending = n_pending;
    for (int iter = 0; iter < max_iters; ++iter) {
      const int told = iter & 1;
      const int tnew = 1 - told;
      const int zero = 0;
      Gpu::htod_memcpy(p_pending, &zero, sizeof(int));

      for (MFIter mfi(new_markers, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);
        auto const& state = state_mf.array(mfi);
        const int nc = cls_t::NCONS;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          const int old_tag = tag(i,j,k,told);
          if (old_tag != pending) {
            tag(i,j,k,tnew) = old_tag;
            return;
          }

          Real sum[cls_t::NCONS] = {};
          int count = 0;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            for (int s = -1; s <= 1; s += 2) {
              int ii = i;
              int jj = j;
              int kk = k;
              if (d == 0) ii += s;
              if (d == 1) jj += s;
#if (AMREX_SPACEDIM == 3)
              if (d == 2) kk += s;
#endif
              const int neighbour_tag = tag(ii,jj,kk,told);
              if (neighbour_tag != donor && neighbour_tag != fixed) continue;

              bool valid_state = state(ii,jj,kk,cls_t::URHO) > Real(0.0)
                              && state(ii,jj,kk,cls_t::UET)  > Real(0.0);
              for (int n = 0; n < nc; ++n) {
                valid_state = valid_state && std::isfinite(state(ii,jj,kk,n));
              }
              if (!valid_state) continue;

              for (int n = 0; n < nc; ++n) {
                sum[n] += state(ii,jj,kk,n);
              }
              ++count;
            }
          }

          if (count > 0) {
            const Real inv = Real(1.0) / Real(count);
            for (int n = 0; n < nc; ++n) {
              state(i,j,k,n) = sum[n] * inv;
            }
            tag(i,j,k,tnew) = fixed;
          } else {
            tag(i,j,k,tnew) = pending;
            Gpu::Atomic::Add(p_pending, 1);
          }
        });
      }

      // Publish newly fixed states and tags before the next Jacobi iteration.
      state_mf.FillBoundary(periodicity);
      tags.FillBoundary(tnew, 1, periodicity);

      n_pending = d_pending.dataValue();
      ParallelDescriptor::ReduceIntSum(n_pending);
      if (n_pending == 0) return;

      if (n_pending >= previous_pending) {
        amrex::Abort(
            "fixExposedCells: exposed fluid region is disconnected from all "
            "old-fluid donors on level " + std::to_string(lev) +
            " (" + std::to_string(n_pending) + " cell(s) remain)");
      }
      previous_pending = n_pending;
    }

    amrex::Abort("fixExposedCells: Jacobi propagation exceeded the level "
                 "domain extent on level " + std::to_string(lev));
  }

  /**
    * \brief Compute surface indices and interpolation data for all faces/edges at given level
    *
    * Algorithm:
    *  1. Build spatial lookup (global_fab_idx -> local_fab_idx)
    *  2. For each face: compute mirror point, find owning FAB, compute interpolation weights
    *  3. Build CSR structure for GPU-friendly access
    *
    * \param lev AMR level
    */
  void computeSurfIndices(int lev)
  {
    BL_PROFILE("IBM::computeSurfIndices");

    int myrank = amrex::ParallelDescriptor::MyProc();
    amrex::Print()  << "Compute Surface Index at LEVEL " << lev  << std::endl;

    auto& mfab = *bmf_a[lev];

    const BoxArray& ba_global     = mfab.boxArray();
    const DistributionMapping& dm = mfab.DistributionMap();

    const int nfab_local  = mfab.local_size();
    const int nfab_global = ba_global.size();

    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto& domain = amr_p->Geom(lev).Domain();

    // ========================================================================
    // Phase 0: Initialize surfimp_soa / surfphys_soa
    // ========================================================================
    if (lev == amr_p->finestLevel()) {
      if ((surfimp_soa.elemIdx.size() != ntotalfaces) || (surfphys_soa.elemIdx.size() != ntotalfaces)) {
        surfimp_soa.resize(ntotalfaces);
        surfphys_soa.resize(ntotalfaces);
      }
      surfphys_soa.reset();
    }

    // ========================================================================
    // Phase 1: Build per-FAB lookup structures
    // ========================================================================
    Vector<Array4<uint8_t const>> fab_markers(nfab_local);
    Vector<Box> fab_bxg(nfab_local);
    Vector<int> global_to_local_fab(nfab_global, -1);

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
        int lidx = mfi.LocalIndex();
        int gidx = mfi.index();
        fab_markers[lidx] = mfab.const_array(mfi);
        fab_bxg[lidx] = mfi.growntilebox(cls_t::NGHOST);
        global_to_local_fab[gidx] = lidx;
    }

    // ========================================================================
    // Phase 2: Fast locate — assign each face to a FAB (lightweight)
    //
    // This pass only does centroid→cell→FAB lookup and writes metadata.
    // The heavy image-point computation is deferred to Phase 3.
    // ========================================================================
    // Per-face FAB assignment: >=0 means locally owned, -1 means not local
    std::vector<int> face_local_fab(ntotalfaces, -1);

    int faces_found = 0;
    int faces_out_domain = 0;
    int faces_out_level  = 0;
    int faces_out_rank   = 0;
    int faces_in_finer   = 0;

    // Precompute inverse dx for faster centroid→cell conversion
    amrex::GpuArray<Real, AMREX_SPACEDIM> inv_dx;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        inv_dx[d] = Real(1.0) / dx_a[lev][d];
    }

    {
      std::vector<std::pair<int, Box>> isects;

      for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {

        // Skip faces already assigned to a finer level
        if (surfphys_soa.elemfound[f_idx] && surfphys_soa.lev[f_idx] > lev) {
            faces_in_finer++;
            continue;
        }

        // Convert body-frame centroid to world frame, then to cell index
        const SurfElem& surfelem = SurfElem_a[f_idx];
        // RZ fix: skip non-physical surface faces lying on the symmetry axis
        // (r~0).  Where a closed axisymmetric contour runs along r=0 through the
        // solid interior (e.g. a long centre body / plenum), those faces have
        // zero circumference (2*pi*r -> 0) and their image point falls inside
        // the solid, so the stencil search always fails and the failed-search
        // host printf floods stdout -- serializing one rank and stalling
        // refinement of long on-axis geometries.  They carry no valid surface
        // data, so skip them.  (RZ-gated: no effect on Cartesian.)
        //
        // FORCE-INTEGRAL NOTE (verified for the SRP center-nozzle cases): the
        // skipped faces are the contour's axis-closure segment (centroid r==0 to
        // ~0.5dx); for a body whose nose is a nozzle bore there is NO solid
        // on-axis apex, so no real surface is lost and the missing 2*pi*r force
        // contribution is negligible (-> 0 as r -> 0).  CAVEAT: a body with a
        // SOLID on-axis apex (e.g. a forebody-only sphere-cone, no centre
        // nozzle) WOULD have a real stagnation-apex surface element here; that
        // case needs the face KEPT with 2*pi*r weighting + a symmetry/mirror
        // fallback for the (in-solid) image point, NOT this blanket skip.
        if (amr_p->Geom(lev).IsRZ()
            && surfelem.centroid[0] < Real(0.5) * dx_a[lev][0]) {
            faces_out_domain++;
            continue;
        }
        int gIdx = getGeomIdx(f_idx);
        const auto& T = transform_a[gIdx];
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(surfelem.centroid[0], surfelem.centroid[1]);
#else
        Point c_body(surfelem.centroid[0], surfelem.centroid[1], surfelem.centroid[2]);
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = surfelem.centroid[dd];
#endif
        Point c_world = T.to_world(c_body);

        IntVect iv;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv[d] = static_cast<int>(std::floor(
                    (c_world[d] - prob_lo[d]) * inv_dx[d]));
        }

        if (!domain.contains(iv)) { faces_out_domain++; continue; }

        ba_global.intersections(Box(iv, iv), isects, true, IntVect::TheZeroVector());
        if (isects.empty()) { faces_out_level++; continue; }

        int gidx       = isects[0].first;
        int owner_rank = dm[gidx];
        int local_fab  = global_to_local_fab[gidx];

        surfphys_soa.ifab[f_idx] = local_fab;
        surfphys_soa.rank[f_idx] = owner_rank;
        surfphys_soa.lev[f_idx]  = lev;

        if (owner_rank != myrank) {
            surfphys_soa.elemfound[f_idx] = -1;
            faces_out_rank++;
            continue;
        }

        surfphys_soa.elemfound[f_idx] = 1;
        face_local_fab[f_idx] = local_fab;
        faces_found++;
      }
    } // isects freed here

    // ========================================================================
    // Phase 3: Group locally-owned faces by FAB for cache locality
    //
    // Faces in the same FAB access the same ibMarkers array.  Processing
    // them together keeps marker data in cache and avoids thrashing when
    // calling search_optimal_image_point / computeIPweights.
    // ========================================================================
    std::vector<std::vector<int>> fab_faces(nfab_local);
    for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {
        if (face_local_fab[f_idx] >= 0) {
            fab_faces[face_local_fab[f_idx]].push_back(f_idx);
        }
    }

    // ========================================================================
    // Phase 4: Compute image points and interpolation weights per FAB
    //
    // Each face writes exclusively to its own f_idx slot in surfimp_soa /
    // surfphys_soa — no cross-face data dependencies.  The inner helpers
    // (search_optimal_image_point, computeIPweights) are stateless and
    // thread-safe, so the outer FAB loop is safe to parallelize with OpenMP.
    // ========================================================================
    int surf_stencil_fail = 0;
    int surf_place_fail = 0;
    int surf_fit_fail = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+:surf_stencil_fail,surf_place_fail,surf_fit_fail)
#endif
    for (int lfab = 0; lfab < nfab_local; ++lfab) {
      const auto& flist = fab_faces[lfab];
      if (flist.empty()) continue;

      // Pin the marker array and box for this FAB — all faces in flist
      // will read from this same data, maximizing cache reuse.
      auto const ibMarkers = fab_markers[lfab];
      auto const bxg       = fab_bxg[lfab];

      for (int f_idx : flist) {
        surfimp_soa.elemIdx[f_idx]  = f_idx;
        surfphys_soa.elemIdx[f_idx] = f_idx;

        // Rotate body-frame LocalFrame to world frame
        const LocalFrame& lf_body = LocalFrame_a[f_idx];
        const SurfElem&   surfelem = SurfElem_a[f_idx];
        int gIdx = getGeomIdx(f_idx);
        const auto& T = transform_a[gIdx];

        LocalFrame localframe;
        T.rotate_to_world(lf_body.normal,   localframe.normal);
        T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
        T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif

        // Transform body-frame centroid to world frame
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(surfelem.centroid[0], surfelem.centroid[1]);
#else
        Point c_body(surfelem.centroid[0], surfelem.centroid[1], surfelem.centroid[2]);
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = surfelem.centroid[dd];
#endif
        Point surf_centroid = T.to_world(c_body);

        // Temporary storage for this face's image point data.  Keep it
        // deterministic so an invalid first IP can be diagnosed and rejected.
        Array2D<Real, 0, eorder_tparm_surf - 1, 0, IDIM> imp_xyz{};
        Array2D< int, 0, eorder_tparm_surf - 1, 0, IDIM> imp_ijk{};
        Array1D<Real, 0, eorder_tparm_surf - 1> disIM{};
        Array1D< int, 0, eorder_tparm_surf - 1> imp_ninterp{};
        int first_place_status = 0;

        for (int jj = 0; jj < eorder_tparm_surf; jj++) {
          Point cp_start;
          if (jj == 0) {
            cp_start = surf_centroid;
          } else {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
            cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1));
#else
            cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1), imp_xyz(jj - 1, 2));
#endif
#else
            for (int d = 0; d < AMREX_SPACEDIM; ++d) cp_start[d] = imp_xyz(jj - 1, d);
#endif
          }

          if (jj == 0) {
            first_place_status =
                search_optimal_image_point<eorder_tparm_surf, iorder_tparm_surf>(
                                          cp_start, localframe,
                                          lev, prob_lo, dx_a[lev], di_a_surf[lev],
                                          bxg, ibMarkers,
                                          surfimp_soa, f_idx,
                                          imp_xyz, imp_ijk, disIM, imp_ninterp);
            if (first_place_status & 1) { ++surf_stencil_fail; }
            if (first_place_status & 2) { ++surf_place_fail; }
            if (first_place_status != 0) { surfphys_soa.elemfound[f_idx] = 0; }
          }
          else {
            search_image_point<eorder_tparm_surf, iorder_tparm_surf>(
                                      jj, cp_start, localframe,
                                      lev, prob_lo, dx_a[lev], di_a_surf[lev],
                                      bxg, ibMarkers,
                                      surfimp_soa, f_idx,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);
          }
        } // end loop on image points

        surfimp_soa.imp_xyz[f_idx]     = imp_xyz;
        surfimp_soa.imp_ijk[f_idx]     = imp_ijk;
        surfimp_soa.disIM[f_idx]       = disIM;

        Array3D< int, 0, eorder_tparm_surf - 1, 0, N_InterP_surf - 1, 0, IDIM> imp_ip_ijk;
        Array2D<Real, 0, eorder_tparm_surf - 1, 0, N_InterP_surf - 1> imp_ipweights;

        computeIPweights<eorder_tparm_surf, iorder_tparm_surf, SURFIMP>(
            imp_ipweights, imp_ip_ijk,
            imp_xyz, imp_ijk, imp_ninterp,
            prob_lo, dx_a[lev], ibMarkers);
        if (first_place_status == 0 && imp_ninterp(0) == 0) {
          ++surf_fit_fail;
          surfphys_soa.elemfound[f_idx] = 0;
        }

        // Store imp_ninterp / ip_quality AFTER computeIPweights (WLS demotion).
        surfimp_soa.imp_ninterp[f_idx] = imp_ninterp;
        surfphys_soa.ip_quality[f_idx] = imp_ninterp(0);
        surfimp_soa.imp_ip_ijk[f_idx]    = imp_ip_ijk;
        surfimp_soa.imp_ipweights[f_idx] = imp_ipweights;
      } // end loop over faces in this FAB
    } // end loop over FABs

    ParallelDescriptor::ReduceIntSum(surf_stencil_fail);
    ParallelDescriptor::ReduceIntSum(surf_place_fail);
    ParallelDescriptor::ReduceIntSum(surf_fit_fail);
    if (surf_stencil_fail > 0) {
      amrex::Print() << "computeSurfIndices: " << surf_stencil_fail
                     << " surface first-IP stencil(s) on level " << lev
                     << " extend outside their owning FAB's grown box.\n";
      amrex::Abort("computeSurfIndices: surface interpolation stencil out of box");
    }
    if (surf_place_fail > 0) {
      amrex::Print() << "computeSurfIndices: " << surf_place_fail
                     << " surface element(s) on level " << lev
                     << " failed first image-point placement.\n";
      amrex::Abort("computeSurfIndices: invalid surface first image point");
    }
    if (surf_fit_fail > 0) {
      amrex::Print() << "computeSurfIndices: " << surf_fit_fail
                     << " surface element(s) on level " << lev
                     << " have no resolvable interpolation fit at the first image point.\n";
      amrex::Abort("computeSurfIndices: invalid surface first image-point fit");
    }

    // Build CSR at the coarsest level (surface is built from finest to coarsest)
    if (lev == 0) buildCSR();

    bool ib_diag = false;
    bool surf_diag = false;
    {
      ParmParse pp("ib");
      pp.query("diag", ib_diag);
      surf_diag = ib_diag;
      pp.query("surf_diag", surf_diag);
    }
    if (surf_diag) reportSurfDiagnostics(lev);
  }

  /**
   * \brief Computes surface properties (pressure, temperature, gradients) for each surface face.
   *
   * This function iterates over all surface element owned by the current process. For each face:
   * 1. Interpolates primitive variables at image points using the pre-computed weights.
   * 2. Applies the wall model to determine surface state (e.g., no-slip, adiabatic/isothermal).
   * 3. Computes gradients (e.g., dT/dn) at the surface.
   * 4. Stores the results back into the Surface Data SoA structure.
   *
   * \param stateprops MultiFab containing the fluid state properties.
   * \param cls        Pointer to the physics/closure class.
   * \param lev        Current AMR level.
   */
  void computeSURFs(MultiFab& prims_mf, const cls_t* cls, int lev) 
  {
    BL_PROFILE("IBM::computeSURFs");

    // Skip if no CSR data for this level
    if (lev >= static_cast<int>(faces_per_level.size())) return;
    auto& csr = faces_per_level[lev];
    const int nfaces_local = static_cast<int>(csr.face_indices.size());
    if (nfaces_local == 0) return;

    amrex::Print() << "computeSURFs at LEVEL " << lev
                   << " (" << nfaces_local << " faces)" << std::endl;

    auto& mfab = *bmf_a[lev];
    const int nfabs_local = mfab.local_size();

    // Build device array of Array4<Real> pointers (one per local FAB)
    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
      }
      Gpu::copyAsync(Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      // Sync the async H2D copy before the host-side reads of managed members
      // (geom_offsets) below: on WSL2, host access to managed memory while a
      // copy is in flight on the stream faults. (computeAllGPs reads these
      // BEFORE its copyAsync, so it is unaffected.)
      Gpu::streamSynchronize();
    }

    auto* prims_arr = d_prims.data();

    // Device-accessible SoA pointers
    const int* d_face_indices  = csr.face_indices.data();

    auto* sp_pressure    = surfphys_soa.pressure.data();
    auto* sp_temperature = surfphys_soa.temperature.data();
    auto* sp_dTdn        = surfphys_soa.dTdn.data();
    auto* sp_tau1        = surfphys_soa.tau1.data();
    auto* sp_tau2        = surfphys_soa.tau2.data();
    const int* sp_ifab   = surfphys_soa.ifab.data();

    const auto* si_imp_ip_ijk    = surfimp_soa.imp_ip_ijk.data();
    const auto* si_imp_ipweights = surfimp_soa.imp_ipweights.data();
    const auto* si_disIM         = surfimp_soa.disIM.data();
    const auto* si_imp_ninterp   = surfimp_soa.imp_ninterp.data();

    auto const* lf_ptr = LocalFrame_a.data();
    auto const* se_ptr = SurfElem_a.data();

    // Capture transforms and geom_offsets for body→world rotation
    AMREX_ALWAYS_ASSERT(ngeom <= MAX_NGEOM);
    const int ngeom_local = amrex::min(ngeom, MAX_NGEOM);
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    GpuArray<int, MAX_NGEOM + 1> geom_off;
    for (int ii = 0; ii < ngeom_local; ii++) {
        transforms[ii] = transform_a[ii];
        geom_off[ii]   = geom_offsets[ii];
    }
    geom_off[ngeom_local] = geom_offsets[ngeom_local];

    auto* copy = this;

    ParallelFor(nfaces_local, [=] AMREX_GPU_DEVICE (int ii) noexcept
    {
      // Global face index from CSR
      int f_idx = d_face_indices[ii];
      int ifab = sp_ifab[f_idx];
      auto prims = prims_arr[ifab];

      // Find geometry index for this face (to get the right transform)
      int gIdx = 0;
      for (int g = 0; g < ngeom_local; ++g) {
          if (f_idx >= geom_off[g] && f_idx < geom_off[g + 1]) { gIdx = g; break; }
      }
      const auto& T = transforms[gIdx];

      // 1) Reconstruct local orthonormal frame — rotate body→world
      const auto& frame = lf_ptr[f_idx];

      amrex::Real n_w[AMREX_SPACEDIM], t1_w[AMREX_SPACEDIM], t2_w[AMREX_SPACEDIM];
      T.rotate_to_world(frame.normal,   n_w);
      T.rotate_to_world(frame.tangent1, t1_w);
#if (AMREX_SPACEDIM == 3)
      T.rotate_to_world(frame.tangent2, t2_w);
#endif

      Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          nvec(d)  = n_w[d];
          t1vec(d) = t1_w[d];
#if (AMREX_SPACEDIM == 3)
          t2vec(d) = t2_w[d];
#else
          t2vec(d) = Real(0.0);
#endif
      }

      // Surface centroid coordinates — transform body→world
      Array1D<Real, 0, AMREX_SPACEDIM - 1> xyz;
      {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(se_ptr[f_idx].centroid[0], se_ptr[f_idx].centroid[1]);
#else
        Point c_body(se_ptr[f_idx].centroid[0], se_ptr[f_idx].centroid[1], se_ptr[f_idx].centroid[2]);
#endif
        Point c_world = T.to_world(c_body);
        xyz(0) = c_world.x(); xyz(1) = c_world.y();
#if (AMREX_SPACEDIM == 3)
        xyz(2) = c_world.z();
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = se_ptr[f_idx].centroid[dd];
        Point c_world = T.to_world(c_body);
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) xyz(dd) = c_world[dd];
#endif
      }

      // 2) Zero-initialize primsNormal
      Array2D<Real, 0, eorder_tparm_surf + 1, 0, cls_t::NPRIM - 1> primsNormal;
      for (int p = 0; p <= eorder_tparm_surf + 1; ++p) {
          for (int n = 0; n < cls_t::NPRIM; ++n) {
              primsNormal(p,n) = Real(0.0);
          }
      }

      // 3) Interpolate primitive variables at image points
      copy->template interpolateIMs<eorder_tparm_surf, iorder_tparm_surf>(
          si_imp_ip_ijk[f_idx], si_imp_ipweights[f_idx], prims, primsNormal);

      // 4) Transform velocities at image points to local frame
      for (int iip = 2; iip < 2 + eorder_tparm_surf; ++iip) {
          copy->template global2local<eorder_tparm_surf>(iip, primsNormal, nvec, t1vec, t2vec);
      }

      // 5) Apply wall model at IB surface (fills primsNormal slot 1).
      //    n_valid_surf = number of leading valid image points; the wall model
      //    and the wall-flux stencils below use the 2nd-order one-sided form
      //    only when >=2 are valid (extrap_order_surf>=2), else degrade to the
      //    1st-order single-image-point form.
      int n_valid_surf = 0;
      for (int kk = 0; kk < eorder_tparm_surf; ++kk) {
          if (si_imp_ninterp[f_idx](kk) < INTERP_THRESHOLD_SURF) break;
          n_valid_surf = kk + 1;
      }
      int type_solid_bc = 0;
      ibm_detail::dispatch_compute_surfIB<eorder_tparm_surf, wallmodel>(
          xyz, nvec, t1vec, t2vec, si_disIM[f_idx], n_valid_surf,
          primsNormal, type_solid_bc, cls);

      // 6) Extract surface quantities
      Real P_surf = primsNormal(1, cls_t::QPRES);
      Real T_surf = primsNormal(1, cls_t::QT);

      // Wall-normal gradients via 2nd-order one-sided stencil (surface + IP1 +
      // IP2, general spacing) when >=2 image points are valid; else 1st-order
      // one-sided (= the previous (IP1 - surface)/dis formula).
      //   dT/dn = grad(T)·n   (heat flux)
      //   tau   = mu * d(u_tangential)/dn   (local frame: QU=normal,
      //           QV=tangent1, QW=tangent2)
      Real dTdn = ibm_wall_normal_deriv<eorder_tparm_surf>(
          primsNormal, cls_t::QT, si_disIM[f_idx], n_valid_surf);

      Real mu = cls->visc(T_surf);
      Real tau1 = mu * ibm_wall_normal_deriv<eorder_tparm_surf>(
          primsNormal, cls_t::QV, si_disIM[f_idx], n_valid_surf);
#if (AMREX_SPACEDIM == 3)
      Real tau2 = mu * ibm_wall_normal_deriv<eorder_tparm_surf>(
          primsNormal, cls_t::QW, si_disIM[f_idx], n_valid_surf);
#else
      Real tau2 = Real(0.0);
#endif

      // 7) Store results
      sp_pressure[f_idx]    = P_surf;
      sp_temperature[f_idx] = T_surf;
      sp_dTdn[f_idx]        = dTdn;
      sp_tau1[f_idx]        = tau1;
      sp_tau2[f_idx]        = tau2;
    }); // end ParallelFor

    // Ensure GPU writes are visible to CPU before gatherSurfData / plotSURF
    Gpu::streamSynchronize();
  }


//============================================================================
///--------------------------- private functions -----------------------------
private:

  // Interpolation, extrapolation, coordinate-transform helpers
  #include "ibm_solver_interp.h"

  // Geometry I/O, VTK output, MPI gather, CSR builder
  #include "ibm_solver_io.h"

}; // end class ibm_solver_t

#endif // IBM_SOLVER_H_
