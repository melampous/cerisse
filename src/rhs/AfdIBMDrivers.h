// Member definitions included by afd_hllc_t.

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void eflux_ibm(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims,
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& /*rhs*/, const cls_t* cls,
      const amrex::Array4<std::uint8_t>& ibMarkers)
  {
    amrex::Array4<amrex::Real> unused_pressure_face_flux;
    compute_ibm_face_fluxes<false>(
        geom, mfi, prims, flxt, cls, ibMarkers, prims, cls_t::QPRES,
        unused_pressure_face_flux);
  }

  void eflux_ibm_with_rz_paired_pressure(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims,
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const amrex::Array4<amrex::Real>& /*rhs*/, const cls_t* cls,
      const amrex::Array4<std::uint8_t>& ibMarkers,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux)
  {
#ifdef AMREX_USE_GPIBM
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.IsRZ(),
        "AFD-HLLC pure shared-GP paired pressure requires R-Z geometry");
    compute_ibm_face_fluxes<true>(
        geom, mfi, prims, flxt, cls, ibMarkers, pressure_states,
        pressure_component, radial_pressure_face_flux);
#else
    amrex::ignore_unused(
        geom, mfi, prims, flxt, cls, ibMarkers, pressure_states,
        pressure_component, radial_pressure_face_flux);
    amrex::Abort(
        "AFD-HLLC paired R-Z pressure is a pure shared-GP/full-cell path, "
        "not an embedded-boundary flow update");
#endif
  }

  template <bool ComputeRzPairedPressure>
  void compute_ibm_face_fluxes(
      const amrex::Geometry& geom, const amrex::MFIter& mfi,
      const amrex::Array4<const amrex::Real>& prims,
      std::array<amrex::FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const cls_t* cls,
      const amrex::Array4<std::uint8_t>& ibMarkers,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const amrex::Array4<amrex::Real>& radial_pressure_face_flux)
  {
    const RuntimeOptions opts = runtime_options();
    const int correction = opts.correction;
    const int shock_llf = opts.shock_llf;
    const int shock_hllc_muscl = opts.shock_hllc_muscl;
    const amrex::Real hllc_reconstructed_speed_ratio =
        opts.hllc_reconstructed_speed_ratio;
    const amrex::Real smoothness_threshold = opts.smoothness_threshold;
    const amrex::Real teno_cutoff = opts.teno_cutoff;
    const amrex::Real ibm_llf_threshold = opts.ibm_llf_threshold;
    const int local_llf_fallback = opts.local_llf_fallback;
    const amrex::Real fallback_density_jump =
        opts.fallback_density_jump;
    const amrex::Real fallback_density_curvature =
        opts.fallback_density_curvature;
    const amrex::Box bx = mfi.tilebox();
    const auto cell_size = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const int radial_axis_face_index = geom.Domain().smallEnd(0);
    const int radial_domain_high_index = geom.Domain().bigEnd(0);
    const bool is_rz = geom.IsRZ();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const bool rz_radial_direction = is_rz && dir == 0;
      const amrex::Box flxbx = amrex::surroundingNodes(bx, dir);
      const auto flx = flxt[dir]->array();
      const amrex::IntVect ivd =
          amrex::IntVect::TheDimensionVector(dir);
      amrex::ParallelFor<128>(
          flxbx,
          [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            using amrex::Real;
            // NVCC 12 rejects a variable first captured inside if constexpr.
            amrex::ignore_unused(
                radial_domain_high_index, radial_pressure_face_flux);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            amrex::IntVect lower_marker = iv - ivd;
            amrex::IntVect upper_marker = iv;
            if (rz_radial_direction && prob_lo[0] == Real(0.0) &&
                lower_marker[0] < radial_axis_face_index) {
              lower_marker[0] =
                  2 * radial_axis_face_index - 1 - lower_marker[0];
            }
            if (rz_radial_direction && prob_lo[0] == Real(0.0) &&
                upper_marker[0] < radial_axis_face_index) {
              upper_marker[0] =
                  2 * radial_axis_face_index - 1 - upper_marker[0];
            }

            if (ibm_flux::is_solid(lower_marker, ibMarkers) &&
                ibm_flux::is_solid(upper_marker, ibMarkers)) {
              ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
              if constexpr (ComputeRzPairedPressure) {
                if (rz_radial_direction) {
                  radial_pressure_face_flux(iv, 0) = Real(0.0);
                }
              }
              return;
            }

            bool stencil_is_all_fluid = true;
            for (int sample = 0; sample < 2 * ng; ++sample) {
              amrex::IntVect marker_cell =
                  iv + (sample - ng) * ivd;
              if (rz_radial_direction && prob_lo[0] == Real(0.0) &&
                  marker_cell[0] < radial_axis_face_index) {
                marker_cell[0] =
                    2 * radial_axis_face_index - 1 - marker_cell[0];
              }
              stencil_is_all_fluid =
                  stencil_is_all_fluid &&
                  !ibm_flux::is_solid(marker_cell, ibMarkers);
            }

            Real pressure_face;
            if (!stencil_is_all_fluid) {
              pressure_face = this->face_flux_ibm_boundary(
                  iv, dir, prims, flx, *cls, ibMarkers, teno_cutoff,
                  ibm_llf_threshold, local_llf_fallback,
                  fallback_density_jump, fallback_density_curvature,
                  rz_radial_direction, prob_lo[0], cell_size[0],
                  radial_axis_face_index, pressure_states,
                  pressure_component);
            } else {
              pressure_face = this->face_flux(
                  iv, dir, prims, flx, *cls, correction, shock_llf,
                  shock_hllc_muscl, hllc_reconstructed_speed_ratio,
                  smoothness_threshold, teno_cutoff,
                  local_llf_fallback, fallback_density_jump,
                  fallback_density_curvature, rz_radial_direction,
                  prob_lo[0], cell_size[0], radial_axis_face_index,
                  pressure_states, pressure_component);
            }

            if constexpr (ComputeRzPairedPressure) {
              if (rz_radial_direction) {
                const bool radial_axis_face =
                    iv[0] == radial_axis_face_index &&
                    prob_lo[0] == Real(0.0);
                if (radial_axis_face) {
                  // A complete all-fluid AFD stencil carries nonzero axis
                  // auxiliary h_(rF) values for the even mass, axial-momentum,
                  // and energy metric fluxes.  A reduced GP boundary stencil
                  // has no such AFD correction and safely degrades to the
                  // physical zero metric flux.
                  if (stencil_is_all_fluid) {
                    flx(iv, cls_t::UMX) = Real(0.0);
                  } else {
                    ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
                  }
                  const bool use_high_order_axis_pressure =
                      correction != 0 && stencil_is_all_fluid &&
                      this->smooth_stencil(
                          iv, ivd, prims, smoothness_threshold);
                  pressure_face =
                      this->axis_pressure_auxiliary_flux_ibm(
                          iv, pressure_states, pressure_component,
                          ibMarkers, radial_domain_high_index,
                          use_high_order_axis_pressure);
                }
                radial_pressure_face_flux(iv, 0) = pressure_face;
                if (!std::isfinite(pressure_face)) {
                  const Real invalid =
                      std::numeric_limits<Real>::quiet_NaN();
                  for (int n = 0; n < cls_t::NCONS; ++n) {
                    flx(iv, n) = invalid;
                  }
                }
              }
            }
          });
    }
  }
#endif
