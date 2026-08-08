// Standard marker-aware shared-GP boundary support included by afd_hllc_t.

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  valid_stencil(const amrex::IntVect& iv, const amrex::IntVect& ivd,
                const amrex::Array4<const amrex::Real>& q) noexcept
  {
    using amrex::Real;
    for (int m = 0; m < 6; ++m) {
      const amrex::IntVect p = iv + (m - 3) * ivd;
      const Real rho = q(p, cls_t::QRHO);
      const Real pressure = q(p, cls_t::QPRES);
      const Real gamma = q(p, cls_t::QG);
      const Real sound = q(p, cls_t::QC);
      const Real finite_sum = rho + pressure + gamma + sound +
                              q(p, cls_t::QU) + q(p, cls_t::QU + 1) +
                              q(p, cls_t::QU + 2);
      if (!(rho > Real(0.0)) || !(pressure > Real(0.0)) ||
          !(gamma > Real(1.0)) || !(sound > Real(0.0)) ||
          !std::isfinite(finite_sum)) {
        return false;
      }
    }
    return true;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  smooth_stencil(const amrex::IntVect& iv, const amrex::IntVect& ivd,
                 const amrex::Array4<const amrex::Real>& q,
                 const amrex::Real threshold) noexcept
  {
    using amrex::Real;
    if (!(threshold > Real(0.0))) {
      return true;  // threshold=0 reproduces the historical unconditional AFD
    }

    Real sound_scale = Real(0.0);
    for (int m = 0; m < 6; ++m) {
      sound_scale = amrex::max(
          sound_scale, std::abs(q(iv + (m - 3) * ivd, cls_t::QC)));
    }
    for (int field = 0; field < 5; ++field) {
      const int comp = field == 0 ? cls_t::QRHO
                       : field == 1 ? cls_t::QPRES
                                    : cls_t::QU + field - 2;
      Real field_scale = Real(0.0);
      for (int m = 0; m < 6; ++m) {
        field_scale = amrex::max(
            field_scale, std::abs(q(iv + (m - 3) * ivd, comp)));
      }
      const Real floor = field >= 2 ? sound_scale : field_scale;
      for (int m = 0; m < 4; ++m) {
        const Real a = q(iv + (m - 3) * ivd, comp);
        const Real b = q(iv + (m - 2) * ivd, comp);
        const Real c = q(iv + (m - 1) * ivd, comp);
        const Real denom = std::abs(a) + Real(2.0) * std::abs(b) +
                           std::abs(c) + floor +
                           std::numeric_limits<Real>::min();
        if (std::abs(a - Real(2.0) * b + c) / denom > threshold) {
          return false;
        }
      }
    }
    return true;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  primitive_to_characteristic(
      const amrex::IntVect& iv, const int dir, const int t1, const int t2,
      const amrex::Array4<const amrex::Real>& q,
      amrex::Real w[6]) noexcept
  {
    using amrex::Real;
    const Real rho = q(iv, cls_t::QRHO);
    const Real pressure = q(iv, cls_t::QPRES);
    const Real gamma = q(iv, cls_t::QG);
    const Real acoustic_scale = std::sqrt(gamma * rho * pressure);
    w[0] = rho * (Real(1.0) - Real(1.0) / gamma);
    w[1] = Real(0.5) *
           (pressure + acoustic_scale * q(iv, cls_t::QU + dir));
    w[2] = Real(0.5) *
           (pressure - acoustic_scale * q(iv, cls_t::QU + dir));
    w[3] = q(iv, cls_t::QU + t1);
    w[4] = q(iv, cls_t::QU + t2);
    w[5] = gamma;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  characteristic_to_primitive(
      const amrex::Real w[6], amrex::Real& rho, amrex::Real& un,
      amrex::Real& ut1, amrex::Real& ut2, amrex::Real& pressure,
      amrex::Real& gamma) noexcept
  {
    using amrex::Real;
    gamma = w[5];
    pressure = w[1] + w[2];
    if (!(gamma > Real(1.0)) || !(pressure > Real(0.0)) ||
        !std::isfinite(gamma + pressure)) {
      return false;
    }
    rho = w[0] / (Real(1.0) - Real(1.0) / gamma);
    if (!(rho > Real(0.0)) || !std::isfinite(rho)) {
      return false;
    }
    un = (w[1] - w[2]) / std::sqrt(gamma * rho * pressure);
    ut1 = w[3];
    ut2 = w[4];
    return std::isfinite(un + ut1 + ut2);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  cell_centered_llf(const amrex::IntVect& iv, const int dir,
                    const amrex::Array4<const amrex::Real>& q,
                    const amrex::Array4<amrex::Real>& flx,
                    const cls_t& cls) noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    const amrex::IntVect il = iv - ivd;
    Real ul[cls_t::NCONS];
    Real ur[cls_t::NCONS];
    Real fl[cls_t::NCONS];
    Real fr[cls_t::NCONS];
    cls.prims2cons(il, q, ul);
    cls.prims2cons(iv, q, ur);
    cls.prims2flux(il, dir, q, fl);
    cls.prims2flux(iv, dir, q, fr);
    const Real alpha = amrex::max(
        std::abs(q(il, cls_t::QU + dir)) + q(il, cls_t::QC),
        std::abs(q(iv, cls_t::QU + dir)) + q(iv, cls_t::QC));
    for (int n = 0; n < cls_t::NCONS; ++n) {
      flx(iv, n) = Real(0.5) * (fl[n] + fr[n] - alpha * (ur[n] - ul[n]));
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  muscl_characteristic_slope(
      const amrex::IntVect& center, const int dir,
      const amrex::Array4<const amrex::Real>& q,
      amrex::Real slope[6]) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    const amrex::IntVect minus = center - ivd;
    const amrex::IntVect plus = center + ivd;
    const int t1 = dir == 0 ? 1 : 0;
    const int t2 = dir == 2 ? 1 : 2;
    const Real rho = q(center, cls_t::QRHO);
    const Real sound = q(center, cls_t::QC) + Real(1.0e-40);
    const Real sound_squared = sound * sound;

    Real dlft =
        Real(0.5) *
            (q(center, cls_t::QPRES) - q(minus, cls_t::QPRES)) / sound -
        Real(0.5) * rho *
            (q(center, cls_t::QU + dir) - q(minus, cls_t::QU + dir));
    Real drgt =
        Real(0.5) *
            (q(plus, cls_t::QPRES) - q(center, cls_t::QPRES)) / sound -
        Real(0.5) * rho *
            (q(plus, cls_t::QU + dir) - q(center, cls_t::QU + dir));
    slope[0] = this->limiter(dlft, drgt);

    dlft = q(center, cls_t::QRHO) - q(minus, cls_t::QRHO) -
           (q(center, cls_t::QPRES) - q(minus, cls_t::QPRES)) /
               sound_squared;
    drgt = q(plus, cls_t::QRHO) - q(center, cls_t::QRHO) -
           (q(plus, cls_t::QPRES) - q(center, cls_t::QPRES)) /
               sound_squared;
    slope[1] = this->limiter(dlft, drgt);

    dlft =
        Real(0.5) *
            (q(center, cls_t::QPRES) - q(minus, cls_t::QPRES)) / sound +
        Real(0.5) * rho *
            (q(center, cls_t::QU + dir) - q(minus, cls_t::QU + dir));
    drgt =
        Real(0.5) *
            (q(plus, cls_t::QPRES) - q(center, cls_t::QPRES)) / sound +
        Real(0.5) * rho *
            (q(plus, cls_t::QU + dir) - q(center, cls_t::QU + dir));
    slope[2] = this->limiter(dlft, drgt);

    dlft = q(center, cls_t::QU + t1) - q(minus, cls_t::QU + t1);
    drgt = q(plus, cls_t::QU + t1) - q(center, cls_t::QU + t1);
    slope[3] = this->limiter(dlft, drgt);

    dlft = q(center, cls_t::QU + t2) - q(minus, cls_t::QU + t2);
    drgt = q(plus, cls_t::QU + t2) - q(center, cls_t::QU + t2);
    slope[4] = this->limiter(dlft, drgt);

    dlft = q(center, cls_t::QFS) - q(minus, cls_t::QFS);
    drgt = q(plus, cls_t::QFS) - q(center, cls_t::QFS);
    slope[5] = this->limiter(dlft, drgt);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
  muscl_hllc_flux(const amrex::IntVect& iv, const int dir,
                  const amrex::Array4<const amrex::Real>& q,
                  const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
                  const bool trace_this_face,
                  amrex::Real* numerical_pressure = nullptr) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    const amrex::IntVect left = iv - ivd;
    const amrex::IntVect right = iv;
    const int t1 = dir == 0 ? 1 : 0;
    const int t2 = dir == 2 ? 1 : 2;
    Real slope_left[6];
    Real slope_right[6];
    muscl_characteristic_slope(left, dir, q, slope_left);
    muscl_characteristic_slope(right, dir, q, slope_right);

    const Real sound_left = q(left, cls_t::QC) + Real(1.0e-40);
    const Real sound_right = q(right, cls_t::QC) + Real(1.0e-40);
    const Real rho_cell_left = q(left, cls_t::QRHO);
    const Real rho_cell_right = q(right, cls_t::QRHO);

    const Real rl =
        rho_cell_left +
        Real(0.5) * ((slope_left[0] + slope_left[2]) / sound_left +
                     slope_left[1]);
    const Real ul =
        q(left, cls_t::QU + dir) +
        Real(0.5) * (slope_left[2] - slope_left[0]) / rho_cell_left;
    const Real pl = q(left, cls_t::QPRES) +
                    Real(0.5) * (slope_left[0] + slope_left[2]) *
                        sound_left;
    const Real ut1l =
        q(left, cls_t::QU + t1) + Real(0.5) * slope_left[3];
    const Real ut2l =
        q(left, cls_t::QU + t2) + Real(0.5) * slope_left[4];
    Real yl[1] = {q(left, cls_t::QFS) + Real(0.5) * slope_left[5]};

    const Real rr =
        rho_cell_right -
        Real(0.5) * ((slope_right[0] + slope_right[2]) / sound_right +
                     slope_right[1]);
    const Real ur =
        q(right, cls_t::QU + dir) -
        Real(0.5) * (slope_right[2] - slope_right[0]) / rho_cell_right;
    const Real pr = q(right, cls_t::QPRES) -
                    Real(0.5) * (slope_right[0] + slope_right[2]) *
                        sound_right;
    const Real ut1r =
        q(right, cls_t::QU + t1) - Real(0.5) * slope_right[3];
    const Real ut2r =
        q(right, cls_t::QU + t2) - Real(0.5) * slope_right[4];
    Real yr[1] = {q(right, cls_t::QFS) - Real(0.5) * slope_right[5]};

    const Real primitive_sum =
        rl + ul + pl + ut1l + ut2l + yl[0] + sound_left + rr + ur + pr +
        ut1r + ut2r + yr[0] + sound_right;
    if (!(rl > Real(0.0)) || !(rr > Real(0.0)) || !(pl > Real(0.0)) ||
        !(pr > Real(0.0)) || !(sound_left > Real(0.0)) ||
        !(sound_right > Real(0.0)) || !std::isfinite(primitive_sum)) {
      return false;
    }

    Real el;
    Real er;
    cls.RYP2E(rl, yl, pl, el);
    cls.RYP2E(rr, yr, pr, er);
    el += Real(0.5) * (ul * ul + ut1l * ut1l + ut2l * ut2l);
    er += Real(0.5) * (ur * ur + ut1r * ut1r + ut2r * ut2r);
    if (!std::isfinite(el + er)) {
      return false;
    }

    Real fn = Real(0.0);
    Real ft1 = Real(0.0);
    Real ft2 = Real(0.0);
    Real fe = Real(0.0);
    Real fry[1] = {Real(0.0)};
    this->hllc(rl, ul, pl, ut1l, ut2l, el, yl, sound_left,
               rr, ur, pr, ut1r, ut2r, er, yr, sound_right,
               fn, ft1, ft2, fe, fry, numerical_pressure);
    if (!std::isfinite(fn + ft1 + ft2 + fe + fry[0])) {
      return false;
    }

    flx(iv, cls_t::UMX + dir) = fn;
    flx(iv, cls_t::UMX + t1) = ft1;
    flx(iv, cls_t::UMX + t2) = ft2;
    flx(iv, cls_t::UET) = fe;
    flx(iv, cls_t::URHO) = fry[0];
    if (trace_this_face) {
      printf("[AFD-FACE-TRACE] branch=hllc_muscl "
             "left=(rho=%.17g,un=%.17g,ut1=%.17g,ut2=%.17g,p=%.17g,c=%.17g,EperMass=%.17g)\n",
             double(rl), double(ul), double(ut1l), double(ut2l),
             double(pl), double(sound_left), double(el));
      printf("[AFD-FACE-TRACE] branch=hllc_muscl "
             "right=(rho=%.17g,un=%.17g,ut1=%.17g,ut2=%.17g,p=%.17g,c=%.17g,EperMass=%.17g)\n",
             double(rr), double(ur), double(ut1r), double(ut2r),
             double(pr), double(sound_right), double(er));
      printf("[AFD-FACE-TRACE] branch=hllc_muscl "
             "flux=(rho=%.17g,mn=%.17g,mt1=%.17g,mt2=%.17g,E=%.17g)\n",
             double(fry[0]), double(fn), double(ft1), double(ft2),
             double(fe));
    }
    return true;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  characteristic_llf_high_order(
      const amrex::IntVect& iv, const int dir,
      const amrex::Array4<const amrex::Real>& q,
      const amrex::Array4<amrex::Real>& flx, const cls_t& cls) noexcept
  {
    using amrex::Real;
    using FluxSplitScheme = typename PointScheme::FluxSplitScheme;
    constexpr int npoint = 6;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    Real alpha = Real(0.0);
    for (int m = 0; m < npoint; ++m) {
      const amrex::IntVect p = iv + (m - 3) * ivd;
      alpha = amrex::max(alpha,
                         std::abs(q(p, cls_t::QU + dir)) + q(p, cls_t::QC));
    }
    alpha = amrex::max(alpha, std::numeric_limits<Real>::epsilon());

    const auto roe_avg = cls.roe_avg_state(iv, dir, q);
    Real fp[npoint][cls_t::NCONS];
    Real fm[npoint][cls_t::NCONS];
    for (int m = 0; m < npoint; ++m) {
      const amrex::IntVect p = iv + (m - 3) * ivd;
      Real cons[cls_t::NCONS];
      Real flux[cls_t::NCONS];
      cls.prims2cons(p, q, cons);
      cls.prims2flux(p, dir, q, flux);
      for (int n = 0; n < cls_t::NCONS; ++n) {
        fp[m][n] = Real(0.5) * (cons[n] + flux[n] / alpha);
        fm[m][n] = Real(0.5) * (cons[n] - flux[n] / alpha);
      }
      cls.cons2char(roe_avg, fp[m]);
      cls.cons2char(roe_avg, fm[m]);
    }

    Real fp_left[cls_t::NCONS];
    Real fm_right[cls_t::NCONS];
    for (int n = 0; n < cls_t::NCONS; ++n) {
      Real stencil[5];
      FluxSplitScheme::gather_positive_split_stencil(n, fp, stencil);
      fp_left[n] =
          FluxSplitScheme::reconstruct_interface_value(stencil);
      FluxSplitScheme::gather_negative_split_stencil(n, fm, stencil);
      fm_right[n] =
          FluxSplitScheme::reconstruct_interface_value(stencil);
    }
    cls.char2cons(roe_avg, fp_left);
    cls.char2cons(roe_avg, fm_right);
    for (int n = 0; n < cls_t::NCONS; ++n) {
      flx(iv, n) = alpha * (fp_left[n] - fm_right[n]);
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  axis_pressure_auxiliary_flux_ibm(
      const amrex::IntVect& axis_face,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component,
      const amrex::Array4<std::uint8_t>& marker,
      const int radial_domain_high_index,
      const bool use_high_order) const noexcept
  {
    using amrex::Real;
    const Real invalid = std::numeric_limits<Real>::quiet_NaN();
    const amrex::IntVect radial =
        amrex::IntVect::TheDimensionVector(0);
    if (axis_face[0] > radial_domain_high_index ||
        !ibm_flux::is_usable(axis_face, marker)) {
      return invalid;
    }
    const Real pressure_ring1 =
        pressure_states(axis_face, pressure_component);
    if (!std::isfinite(pressure_ring1)) {
      return invalid;
    }
    if (!use_high_order ||
        axis_face[0] + 1 > radial_domain_high_index ||
        !ibm_flux::is_usable(axis_face + radial, marker)) {
      return pressure_ring1;
    }
    const Real pressure_ring2 =
        pressure_states(axis_face + radial, pressure_component);
    if (!std::isfinite(pressure_ring2)) {
      return invalid;
    }
    if (axis_face[0] + 2 > radial_domain_high_index ||
        !ibm_flux::is_usable(axis_face + 2 * radial, marker)) {
      return (Real(7.0) * pressure_ring1 - pressure_ring2) /
             Real(6.0);
    }
    const Real pressure_ring3 =
        pressure_states(axis_face + 2 * radial, pressure_component);
    if (!std::isfinite(pressure_ring3)) {
      return invalid;
    }
    return (Real(37.0) * pressure_ring1 -
            Real(8.0) * pressure_ring2 + pressure_ring3) /
           Real(30.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
  face_flux_ibm_boundary(
      const amrex::IntVect& iv, const int dir,
      const amrex::Array4<const amrex::Real>& q,
      const amrex::Array4<amrex::Real>& flx, const cls_t& cls,
      const amrex::Array4<std::uint8_t>& marker,
      const amrex::Real teno_cutoff,
      const amrex::Real ibm_llf_threshold,
      const int local_llf_fallback,
      const amrex::Real fallback_density_jump,
      const amrex::Real fallback_density_curvature,
      const bool rz_radial_direction,
      const amrex::Real radial_origin,
      const amrex::Real radial_spacing,
      const int radial_axis_face_index,
      const amrex::Array4<const amrex::Real>& pressure_states,
      const int pressure_component) const noexcept
  {
    using amrex::Real;
    const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
    const Real invalid = std::numeric_limits<Real>::quiet_NaN();
    bool usable[6];
    for (int m = 0; m < 6; ++m) {
      const amrex::IntVect point = iv + (m - 3) * ivd;
      amrex::IntVect marker_point = point;
      if (rz_radial_direction && radial_origin == Real(0.0) &&
          marker_point[0] < radial_axis_face_index) {
        marker_point[0] =
            2 * radial_axis_face_index - 1 - marker_point[0];
      }
      usable[m] = ibm_flux::is_usable(marker_point, marker);
    }

    if (local_first_order_llf_needed(
            iv, ivd, q, usable, local_llf_fallback,
            fallback_density_jump, fallback_density_curvature,
            rz_radial_direction, radial_origin, radial_spacing,
            radial_axis_face_index)) {
      if (usable[2] && usable[3]) {
        ibm_flux::llf_flux(iv, dir, q, flx, cls);
        return adjacent_pressure_average(
            iv, ivd, pressure_states, pressure_component);
      }
      ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
      return invalid;
    }

    // Reconstructing independent left/right states across a strong discrete
    // IBM jump can create HLLC energy anti-diffusion when the solid state is a
    // single reconstructed GP.  Preserve masked HLLC for smooth wall-compatible
    // data, and use two-point LLF only when the GP-fluid jump sensor fires.
    amrex::IntVect left_marker = iv - ivd;
    amrex::IntVect right_marker = iv;
    if (rz_radial_direction && radial_origin == Real(0.0) &&
        left_marker[0] < radial_axis_face_index) {
      left_marker[0] =
          2 * radial_axis_face_index - 1 - left_marker[0];
    }
    if (rz_radial_direction && radial_origin == Real(0.0) &&
        right_marker[0] < radial_axis_face_index) {
      right_marker[0] =
          2 * radial_axis_face_index - 1 - right_marker[0];
    }
    const bool crosses_ibm =
        ibm_flux::is_solid(right_marker, marker) !=
        ibm_flux::is_solid(left_marker, marker);
    if (crosses_ibm) {
      const Real rho_l = q(iv - ivd, cls_t::QRHO);
      const Real rho_r = q(iv, cls_t::QRHO);
      const Real p_l = q(iv - ivd, cls_t::QPRES);
      const Real p_r = q(iv, cls_t::QPRES);
      const Real c_l = q(iv - ivd, cls_t::QC);
      const Real c_r = q(iv, cls_t::QC);
      const Real finite_sum = rho_l + rho_r + p_l + p_r + c_l + c_r;
      const bool adjacent_valid =
          usable[2] && usable[3] && rho_l > Real(0.0) &&
          rho_r > Real(0.0) && p_l > Real(0.0) && p_r > Real(0.0) &&
          c_l > Real(0.0) && c_r > Real(0.0) &&
          std::isfinite(finite_sum);
      if (!adjacent_valid) {
        ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
        return invalid;
      }

      const Real tiny = std::numeric_limits<Real>::min();
      const Real rho_jump = std::abs(rho_r - rho_l) /
                            (std::abs(rho_r) + std::abs(rho_l) + tiny);
      const Real p_jump = std::abs(p_r - p_l) /
                          (std::abs(p_r) + std::abs(p_l) + tiny);
      Real velocity_jump = Real(0.0);
      for (int component = 0; component < 3; ++component) {
        const Real u_l = q(iv - ivd, cls_t::QU + component);
        const Real u_r = q(iv, cls_t::QU + component);
        velocity_jump = amrex::max(
            velocity_jump,
            std::abs(u_r - u_l) /
                (std::abs(u_r) + std::abs(u_l) + c_l + c_r + tiny));
      }
      const Real jump = amrex::max(rho_jump, amrex::max(p_jump, velocity_jump));
      if (jump > ibm_llf_threshold) {
        ibm_flux::llf_flux(iv, dir, q, flx, cls);
        return adjacent_pressure_average(
            iv, ivd, pressure_states, pressure_component);
      }
    }

    const bool valid_l[3] = {
        usable[0] && usable[1] && usable[2],
        usable[1] && usable[2] && usable[3],
        usable[2] && usable[3] && usable[4]};
    const bool valid_r[3] = {
        usable[1] && usable[2] && usable[3],
        usable[2] && usable[3] && usable[4],
        usable[3] && usable[4] && usable[5]};
    if (!(valid_l[0] || valid_l[1] || valid_l[2]) ||
        !(valid_r[0] || valid_r[1] || valid_r[2])) {
      if (usable[2] && usable[3]) {
        ibm_flux::llf_flux(iv, dir, q, flx, cls);
        return adjacent_pressure_average(
            iv, ivd, pressure_states, pressure_component);
      } else {
        ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
        return invalid;
      }
    }

    const int t1 = dir == 0 ? 1 : 0;
    const int t2 = dir == 2 ? 1 : 2;
    Real wc[6][6] = {};
    for (int m = 0; m < 6; ++m) {
      if (usable[m]) {
        primitive_to_characteristic(iv + (m - 3) * ivd, dir, t1, t2, q,
                                    wc[m]);
      }
    }

    Real wl[6];
    Real wr[6];
    for (int n = 0; n < 6; ++n) {
      Real sl[5];
      Real sr[5];
      for (int m = 0; m < 5; ++m) {
        sl[m] = wc[m][n];
        sr[m] = wc[m + 1][n];
      }
      wl[n] = PointScheme::right_masked(sl, valid_l, teno_cutoff);
      wr[n] = PointScheme::left_masked(sr, valid_r, teno_cutoff);
    }

    Real rl, ul, ut1l, ut2l, pl, gl;
    Real rr, ur, ut1r, ut2r, pr, gr;
    if (!characteristic_to_primitive(wl, rl, ul, ut1l, ut2l, pl, gl) ||
        !characteristic_to_primitive(wr, rr, ur, ut1r, ut2r, pr, gr)) {
      if (usable[2] && usable[3]) {
        ibm_flux::llf_flux(iv, dir, q, flx, cls);
        return adjacent_pressure_average(
            iv, ivd, pressure_states, pressure_component);
      }
      ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
      return invalid;
    }

    const Real yl[1] = {Real(1.0)};
    const Real yr[1] = {Real(1.0)};
    Real el;
    Real er;
    cls.RYP2E(rl, yl, pl, el);
    cls.RYP2E(rr, yr, pr, er);
    el += Real(0.5) * (ul * ul + ut1l * ut1l + ut2l * ut2l);
    er += Real(0.5) * (ur * ur + ut1r * ut1r + ut2r * ut2r);
    const Real cl = std::sqrt(gl * pl / rl);
    const Real cr = std::sqrt(gr * pr / rr);

    Real fn = Real(0.0);
    Real ft1 = Real(0.0);
    Real ft2 = Real(0.0);
    Real fe = Real(0.0);
    Real fry[1] = {Real(0.0)};
    Real pressure_face = invalid;
    this->hllc(rl, ul, pl, ut1l, ut2l, el, yl, cl,
               rr, ur, pr, ut1r, ut2r, er, yr, cr,
               fn, ft1, ft2, fe, fry, &pressure_face);
    if (!std::isfinite(fn + ft1 + ft2 + fe + fry[0] + pressure_face)) {
      if (usable[2] && usable[3]) {
        ibm_flux::llf_flux(iv, dir, q, flx, cls);
        return adjacent_pressure_average(
            iv, ivd, pressure_states, pressure_component);
      }
      ibm_flux::zero_flux<decltype(flx), cls_t>(iv, flx);
      return invalid;
    }

    flx(iv, cls_t::UMX + dir) = fn;
    flx(iv, cls_t::UMX + t1) = ft1;
    flx(iv, cls_t::UMX + t2) = ft2;
    flx(iv, cls_t::UET) = fe;
    flx(iv, cls_t::URHO) = fry[0];
    return pressure_face;
  }
