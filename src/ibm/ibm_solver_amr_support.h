// Member definitions included by ibm_solver_t.
//
// This fragment computes IBM interpolation halo requirements and audits
// coarse-fine support coverage. These checks protect the configured
// reconstruction order and remain part of production builds.

  int volumeInterpolationNghost(int lev) const
  {
    AMREX_ALWAYS_ASSERT(
        lev >= 0 && lev < static_cast<int>(volume_interp_nghost_a.size()));
    return volume_interp_nghost_a[lev];
  }

  int surfaceInterpolationNghost(int lev) const
  {
    AMREX_ALWAYS_ASSERT(
        lev >= 0 && lev < static_cast<int>(surface_interp_nghost_a.size()));
    return surface_interp_nghost_a[lev];
  }

  int interpolationMarkerNghost(int lev) const
  {
    return amrex::max(volumeInterpolationNghost(lev),
                      surfaceInterpolationNghost(lev));
  }

  static bool amrSupportBufferEnabled()
  {
    static const bool enabled = [] {
      int value = 1;
      ParmParse pp("ib");
      pp.query("amr_support_buffer", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1,
          "ib.amr_support_buffer must be 0 or 1");
      return value != 0;
    }();
    return enabled;
  }

  static bool amrSupportAuditEnabled()
  {
    static const bool enabled = [] {
      int value = 1;
      ParmParse pp("ib");
      pp.query("amr_support_audit", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1,
          "ib.amr_support_audit must be 0 or 1");
      return value != 0;
    }();
    return enabled;
  }

  static bool amrSupportStrict()
  {
    static const bool strict = [] {
      int value = 1;
      ParmParse pp("ib");
      pp.query("amr_support_strict", value);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          value == 0 || value == 1,
          "ib.amr_support_strict must be 0 or 1");
      return value != 0;
    }();
    return strict;
  }

  /**
   * Coarse-level tag radius required to keep the next level's complete IBM
   * image-point stencil on same-level valid cells.  Two coarse guard cells
   * cover the wall-to-GP offset and grid-clustering edge placement; the exact
   * nonzero support set is checked independently after the child level is
   * built.
   */
  IntVect amrSupportTagRadius(int lev) const
  {
    AMREX_ALWAYS_ASSERT(lev >= 0 && lev < amr_p->maxLevel());
    const IntVect ratio = amr_p->refRatio(lev);
    const int fine_halo = interpolationMarkerNghost(lev + 1);
    IntVect radius(0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      radius[d] = (fine_halo + ratio[d] - 1) / ratio[d] + 2;
    }
    return radius;
  }

  Box interpolationReadBox(const MFIter& mfi, int lev, int ngrow) const
  {
    Box result = mfi.growntilebox(ngrow);
    // Same-level exchange supplies the enlarged halo only from valid cells.
    // Periodic FillBoundary can populate the full IBM halo. At a non-periodic
    // physical boundary, FillPatch guarantees cls_t::NGHOST and no more, so
    // do not let an image-point stencil consume uninitialised data.
    const Geometry& level_geom = amr_p->Geom(lev);
    Box available = level_geom.Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      available.grow(d, level_geom.isPeriodic(d) ? ngrow : cls_t::NGHOST);
    }
    result &= available;
    return result;
  }

  struct AMRSupportCoverageSummary {
    Long targets = 0;
    Long active_supports = 0;
    Long physical_boundary_supports = 0;
    Long missing_targets = 0;
    Long missing_supports = 0;
  };

  /**
   * Map periodic support indices back into the level domain.  A support beyond
   * a non-periodic physical boundary is supplied by the physical BC and is not
   * a coarse-fine coverage failure.  Interior points must belong to the level's
   * global valid BoxArray; a marker ghost cell alone is not sufficient.
   */
  bool sameLevelSupportCovered(int lev, IntVect support,
                               bool& physical_boundary) const
  {
    const Geometry& level_geom = amr_p->Geom(lev);
    const Box& domain = level_geom.Domain();
    physical_boundary = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (support[d] >= domain.smallEnd(d) &&
          support[d] <= domain.bigEnd(d)) {
        continue;
      }
      if (!level_geom.isPeriodic(d)) {
        physical_boundary = true;
        return true;
      }
      const int period = domain.length(d);
      while (support[d] < domain.smallEnd(d)) support[d] += period;
      while (support[d] > domain.bigEnd(d)) support[d] -= period;
    }
    return bmf_a[lev]->boxArray().contains(support);
  }

  void reportAMRSupportCoverage(const char* target_kind, int lev,
                                AMRSupportCoverageSummary summary,
                                const std::string& representative) const
  {
    if (summary.missing_supports > 0) {
      amrex::AllPrint() << "[IBM-AMR-Support] rank "
                        << ParallelDescriptor::MyProc() << " "
                        << representative << "\n";
    }
    ParallelDescriptor::ReduceLongSum(summary.targets);
    ParallelDescriptor::ReduceLongSum(summary.active_supports);
    ParallelDescriptor::ReduceLongSum(summary.physical_boundary_supports);
    ParallelDescriptor::ReduceLongSum(summary.missing_targets);
    ParallelDescriptor::ReduceLongSum(summary.missing_supports);

    amrex::Print() << "[IBM-AMR-Support] " << target_kind
                   << " level=" << lev
                   << " targets=" << summary.targets
                   << " active-supports=" << summary.active_supports
                   << " physical-BC-supports="
                   << summary.physical_boundary_supports
                   << " coarse-fine-missing-targets="
                   << summary.missing_targets
                   << " coarse-fine-missing-supports="
                   << summary.missing_supports << "\n";

    if (amrSupportStrict() && summary.missing_supports > 0) {
      amrex::Abort(
          std::string("IBM AMR ") + target_kind +
          " interpolation support crosses a coarse-fine boundary. Keep "
          "ib.amr_support_buffer=1, increase AMR proper nesting if needed, "
          "and rebuild the hierarchy; coarse interpolation is not a valid "
          "replacement for the configured IBM reconstruction order.");
    }
  }

  void auditGPAMRSupportCoverage(int lev) const
  {
    if (!amrSupportAuditEnabled() || amr_p->maxLevel() == 0) return;
    Gpu::streamSynchronize();

    AMRSupportCoverageSummary summary;
    std::string representative;
    const auto& gpstore = gpstore_a[lev];
    for (int gp = 0; gp < gpstore.total_ngps; ++gp) {
      for (int image = 0; image < eorder_tparm; ++image) {
        bool active_target = false;
        bool missing_target = false;
        for (int support_id = 0; support_id < N_InterP; ++support_id) {
          if (gpstore.imp_ipweights[gp](image, support_id) == Real(0.0)) {
            continue;
          }
          active_target = true;
          ++summary.active_supports;
          IntVect support(AMREX_D_DECL(
              gpstore.imp_ip_ijk[gp](image, support_id, 0),
              gpstore.imp_ip_ijk[gp](image, support_id, 1),
              gpstore.imp_ip_ijk[gp](image, support_id, 2)));
          bool physical_boundary = false;
          if (sameLevelSupportCovered(lev, support, physical_boundary)) {
            if (physical_boundary) ++summary.physical_boundary_supports;
            continue;
          }
          ++summary.missing_supports;
          missing_target = true;
          if (representative.empty()) {
            std::ostringstream os;
            os << "GP representative: gp=" << gp
               << " image=" << image << " support=" << support;
            representative = os.str();
          }
        }
        if (active_target) ++summary.targets;
        if (missing_target) ++summary.missing_targets;
      }
    }
    reportAMRSupportCoverage("GP", lev, summary, representative);
  }

  void auditSurfaceAMRSupportCoverage(int lev) const
  {
    if (!amrSupportAuditEnabled() || amr_p->maxLevel() == 0) return;
    Gpu::streamSynchronize();

    AMRSupportCoverageSummary summary;
    std::string representative;
    const int my_rank = ParallelDescriptor::MyProc();
    for (int face = 0; face < ntotalfaces; ++face) {
      if (surfphys_soa.elemfound[face] != 1 ||
          surfphys_soa.lev[face] != lev ||
          surfphys_soa.rank[face] != my_rank) {
        continue;
      }
      for (int image = 0; image < eorder_tparm_surf; ++image) {
        bool active_target = false;
        bool missing_target = false;
        for (int support_id = 0; support_id < N_InterP_surf; ++support_id) {
          if (surfimp_soa.imp_ipweights[face](image, support_id) ==
              Real(0.0)) {
            continue;
          }
          active_target = true;
          ++summary.active_supports;
          IntVect support(AMREX_D_DECL(
              surfimp_soa.imp_ip_ijk[face](image, support_id, 0),
              surfimp_soa.imp_ip_ijk[face](image, support_id, 1),
              surfimp_soa.imp_ip_ijk[face](image, support_id, 2)));
          bool physical_boundary = false;
          if (sameLevelSupportCovered(lev, support, physical_boundary)) {
            if (physical_boundary) ++summary.physical_boundary_supports;
            continue;
          }
          ++summary.missing_supports;
          missing_target = true;
          if (representative.empty()) {
            std::ostringstream os;
            os << "surface representative: face=" << face
               << " image=" << image << " support=" << support;
            representative = os.str();
          }
        }
        if (active_target) ++summary.targets;
        if (missing_target) ++summary.missing_targets;
      }
    }
    reportAMRSupportCoverage("surface", lev, summary, representative);
  }
