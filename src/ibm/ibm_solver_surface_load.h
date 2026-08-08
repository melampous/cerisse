// Member definitions included by ibm_solver_t.
//
// Physical BI pressure/traction integration. Geometry I/O remains in
// ibm_solver_io.h; this fragment must be included inside the class body.

  /// Integrate force on all fixed IBM bodies and audit the oriented surface.
  /// The normal points from solid to fluid, so force on the body is
  ///   integral[-p n + tau_nn n + tau_n1 t_1 + tau_n2 t_2] dA.
  /// Each surface element contributes on its owning rank only; the returned
  /// values are MPI-global on every rank.
  ibm_surface_force_audit_t integrateSurfaceForce(bool include_viscous) const
  {
    BL_PROFILE("IBM::integrateSurfaceForce");
    if constexpr (!supports_complete_viscous_surface_traction) {
      if (include_viscous) {
        amrex::Abort(
            "complete IBM viscous force requires a wall model that declares "
            "stationary_no_slip=true");
      }
    }
    Gpu::streamSynchronize();
    ibm_surface_force_audit_t audit{};

#ifdef AMREX_USE_GPU
    // Surface fields and geometry are managed arrays written by device kernels.
    // Reducing them on the host would migrate every page CPU<->GPU each budget
    // step.  Keep the traversal on device and copy back only the audit scalars.
    constexpr int pressure_offset = 0;
    constexpr int viscous_offset = 3;
    constexpr int force_offset = 6;
    constexpr int normal_offset = 9;
    constexpr int measure_index = 12;
    constexpr int faces_index = 13;
    constexpr int nonfinite_index = 14;
    constexpr int num_sums = 15;

    // Pack managed geometry metadata before launching any device work.  Host
    // access to managed memory while a kernel is in flight can fault on WSL2.
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    GpuArray<int, MAX_NGEOM + 1> geometry_offsets;
    for (int geom_idx = 0; geom_idx < ngeom; ++geom_idx) {
      transforms[geom_idx] = transform_a[geom_idx];
      geometry_offsets[geom_idx] = geom_offsets[geom_idx];
    }
    geometry_offsets[ngeom] = geom_offsets[ngeom];
    const int geometry_count = ngeom;
    const int face_count = ntotalfaces;

    Gpu::DeviceVector<Real> device_sums(num_sums);
    Gpu::fillAsync(
        device_sums.begin(), device_sums.end(),
        [] AMREX_GPU_HOST_DEVICE(Real& value, Long) noexcept {
          value = Real(0.0);
        });

    const auto* elemfound = surfphys_soa.elemfound.data();
    const auto* pressure = surfphys_soa.pressure.data();
    const auto* tau_normal_values = surfphys_soa.tau_normal.data();
    const auto* tau1_values = surfphys_soa.tau1.data();
    const auto* tau2_values = surfphys_soa.tau2.data();
    const auto* surface_elements = SurfElem_a.data();
    const auto* local_frames = LocalFrame_a.data();
    auto* sums = device_sums.data();

    ParallelFor(face_count, [=] AMREX_GPU_DEVICE(int f_idx) noexcept {
      if (elemfound[f_idx] != 1) return;

      int geom_idx = -1;
      for (int candidate = 0; candidate < geometry_count; ++candidate) {
        if (f_idx >= geometry_offsets[candidate] &&
            f_idx < geometry_offsets[candidate + 1]) {
          geom_idx = candidate;
          break;
        }
      }
      if (geom_idx < 0) {
        Gpu::Atomic::Add(&sums[nonfinite_index], Real(1.0));
        return;
      }

      const Real p = pressure[f_idx];
      const Real tau_normal =
          include_viscous ? tau_normal_values[f_idx] : Real(0.0);
      const Real tau1 = include_viscous ? tau1_values[f_idx] : Real(0.0);
#if (AMREX_SPACEDIM == 3)
      const Real tau2 = include_viscous ? tau2_values[f_idx] : Real(0.0);
#else
      const Real tau2 = Real(0.0);
#endif
      const Real measure = surface_elements[f_idx].measure;
      if (!(measure > Real(0.0)) || !amrex::Math::isfinite(measure) ||
          !amrex::Math::isfinite(p) ||
          !amrex::Math::isfinite(tau_normal) ||
          !amrex::Math::isfinite(tau1) || !amrex::Math::isfinite(tau2)) {
        Gpu::Atomic::Add(&sums[nonfinite_index], Real(1.0));
        return;
      }

      Real normal[AMREX_SPACEDIM] = {};
      Real tangent1[AMREX_SPACEDIM] = {};
      Real tangent2[AMREX_SPACEDIM] = {};
      transforms[geom_idx].rotate_to_world(
          local_frames[f_idx].normal, normal);
      transforms[geom_idx].rotate_to_world(
          local_frames[f_idx].tangent1, tangent1);
#if (AMREX_SPACEDIM == 3)
      transforms[geom_idx].rotate_to_world(
          local_frames[f_idx].tangent2, tangent2);
#endif

      bool finite_frame = true;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        finite_frame = finite_frame && amrex::Math::isfinite(normal[d]) &&
                       amrex::Math::isfinite(tangent1[d]) &&
                       amrex::Math::isfinite(tangent2[d]);
      }
      if (!finite_frame) {
        Gpu::Atomic::Add(&sums[nonfinite_index], Real(1.0));
        return;
      }

      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real pressure_component = -p * normal[d] * measure;
        const Real viscous_component =
            (tau_normal * normal[d] + tau1 * tangent1[d] +
             tau2 * tangent2[d]) * measure;
        Gpu::Atomic::Add(
            &sums[pressure_offset + d], pressure_component);
        Gpu::Atomic::Add(
            &sums[viscous_offset + d], viscous_component);
        Gpu::Atomic::Add(
            &sums[force_offset + d], pressure_component + viscous_component);
        Gpu::Atomic::Add(&sums[normal_offset + d], normal[d] * measure);
      }
      Gpu::Atomic::Add(&sums[measure_index], measure);
      Gpu::Atomic::Add(&sums[faces_index], Real(1.0));
    });

    Gpu::HostVector<Real> host_sums(num_sums);
    Gpu::copy(
        Gpu::deviceToHost, device_sums.begin(), device_sums.end(),
        host_sums.begin());
    Gpu::streamSynchronize();
    for (int d = 0; d < 3; ++d) {
      audit.pressure_force[d] = host_sums[pressure_offset + d];
      audit.viscous_force[d] = host_sums[viscous_offset + d];
      audit.force[d] = host_sums[force_offset + d];
      audit.normal_closure[d] = host_sums[normal_offset + d];
    }
    audit.measure = host_sums[measure_index];
    audit.faces = static_cast<Long>(std::llround(host_sums[faces_index]));
    audit.nonfinite_faces =
        static_cast<Long>(std::llround(host_sums[nonfinite_index]));
#else
    for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {
      if (surfphys_soa.elemfound[f_idx] != 1) continue;
      const int geom_idx = getGeomIdx(f_idx);
      if (geom_idx < 0 || geom_idx >= ngeom) {
        ++audit.nonfinite_faces;
        continue;
      }

      const Real pressure = surfphys_soa.pressure[f_idx];
      const Real tau_normal =
          include_viscous ? surfphys_soa.tau_normal[f_idx] : Real(0.0);
      const Real tau1 = include_viscous ? surfphys_soa.tau1[f_idx] : Real(0.0);
#if (AMREX_SPACEDIM == 3)
      const Real tau2 = include_viscous ? surfphys_soa.tau2[f_idx] : Real(0.0);
#else
      const Real tau2 = Real(0.0);
#endif
      const Real measure = SurfElem_a[f_idx].measure;
      if (!(measure > Real(0.0)) || !amrex::Math::isfinite(measure) ||
          !amrex::Math::isfinite(pressure) ||
          !amrex::Math::isfinite(tau_normal) ||
          !amrex::Math::isfinite(tau1) || !amrex::Math::isfinite(tau2)) {
        ++audit.nonfinite_faces;
        continue;
      }

      const auto& transform = transform_a[geom_idx];
      Real normal[AMREX_SPACEDIM] = {};
      Real tangent1[AMREX_SPACEDIM] = {};
      Real tangent2[AMREX_SPACEDIM] = {};
      transform.rotate_to_world(LocalFrame_a[f_idx].normal, normal);
      transform.rotate_to_world(LocalFrame_a[f_idx].tangent1, tangent1);
#if (AMREX_SPACEDIM == 3)
      transform.rotate_to_world(LocalFrame_a[f_idx].tangent2, tangent2);
#endif

      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real pressure_component = -pressure * normal[d] * measure;
        const Real viscous_component =
            (tau_normal * normal[d] + tau1 * tangent1[d] +
             tau2 * tangent2[d]) * measure;
        audit.pressure_force[d] += pressure_component;
        audit.viscous_force[d] += viscous_component;
        audit.force[d] += pressure_component + viscous_component;
        audit.normal_closure[d] += normal[d] * measure;
      }
      audit.measure += measure;
      ++audit.faces;
    }
#endif

    ParallelDescriptor::ReduceRealSum(audit.force.data(), 3);
    ParallelDescriptor::ReduceRealSum(audit.pressure_force.data(), 3);
    ParallelDescriptor::ReduceRealSum(audit.viscous_force.data(), 3);
    ParallelDescriptor::ReduceRealSum(audit.normal_closure.data(), 3);
    ParallelDescriptor::ReduceRealSum(audit.measure);
    ParallelDescriptor::ReduceLongSum(audit.faces);
    ParallelDescriptor::ReduceLongSum(audit.nonfinite_faces);
    return audit;
  }
