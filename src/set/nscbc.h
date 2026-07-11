#ifndef CNS_NSCBC_H_
#define CNS_NSCBC_H_

#include <AMReX_Array4.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

// Time-evolved ghost-cell NSCBC/LODI kernel.  The formulation follows the
// persistent ghost-state architecture on Salvador's nscbc branch, with a
// corrected entropy characteristic and an explicit projected transverse term
// in the incoming acoustic amplitude.
namespace nscbc {

static constexpr int LMINUS = 0;
static constexpr int LENT   = 1;
static constexpr int LTAN1  = 2;
static constexpr int LTAN2  = 3;
static constexpr int LPLUS  = 4;
static constexpr int LSP    = 5;

struct Parm {
  amrex::Real Lchar = amrex::Real(1.0);
  amrex::Real Mmax = amrex::Real(0.0);
  amrex::Real Ptarget = amrex::Real(1.0);
  amrex::Real sigma = amrex::Real(0.0);
  amrex::Real utarget = amrex::Real(0.0);
  amrex::Real vtarget = amrex::Real(0.0);
  amrex::Real wtarget = amrex::Real(0.0);
  amrex::Real Ttarget = amrex::Real(1.0);
  amrex::Real eta = amrex::Real(0.0);

  // beta=1 applies the multidimensional characteristic compatibility term;
  // beta=0 recovers a strictly local one-dimensional LODI boundary.
  // The native oblique-acoustic sweep shows that the full unweighted
  // transverse cancellation over-forces the incoming mode.  0.25 is the
  // current generic default; retain this as a runtime-tunable model constant.
  amrex::Real transverse_relax = amrex::Real(0.25);
  amrex::Real min_rho = amrex::Real(1.e-12);
  amrex::Real min_T = amrex::Real(1.e-12);
  amrex::Real min_p = amrex::Real(1.e-12);
  int derivative_order = 2;
  int use_transverse = 1;
};

template <typename cls_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int qvel(int dir)
{
  return cls_t::QU + dir;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void tangent_dirs(
  int dir, int& t1, int& t2)
{
#if (AMREX_SPACEDIM == 1)
  amrex::ignore_unused(dir);
  t1 = -1;
  t2 = -1;
#elif (AMREX_SPACEDIM == 2)
  t1 = 1 - dir;
  t2 = -1;
#else
  t1 = (dir + 1) % 3;
  t2 = (dir + 2) % 3;
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void one_sided_deriv_prim(
  amrex::IntVect const& iv, int dir, int side_sign, amrex::Real dxinv,
  amrex::Array4<const amrex::Real> const& q, amrex::Real& drho,
  amrex::Real& dun, amrex::Real& dut1, amrex::Real& dut2,
  amrex::Real& dT, amrex::Real* dY, int derivative_order)
{
  using amrex::Real;
  const auto e = amrex::IntVect::TheDimensionVector(dir) * side_sign;
  const Real sdx = Real(side_sign) * dxinv;
  int t1, t2;
  tangent_dirs(dir, t1, t2);

  auto deriv = [&] AMREX_GPU_DEVICE (int n) noexcept -> Real {
    if (derivative_order == 2) {
      return sdx * (-Real(1.5) * q(iv, n) + Real(2.0) * q(iv + e, n)
                    - Real(0.5) * q(iv + e + e, n));
    }
    return sdx * (-q(iv, n) + q(iv + e, n));
  };

  drho = deriv(cls_t::QRHO);
  dun = deriv(qvel<cls_t>(dir));
  dut1 = (t1 >= 0) ? deriv(qvel<cls_t>(t1)) : Real(0.0);
  dut2 = (t2 >= 0) ? deriv(qvel<cls_t>(t2)) : Real(0.0);
  dT = deriv(cls_t::QT);
#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) {
    dY[n] = deriv(cls_t::QFS + n);
  }
#else
  amrex::ignore_unused(dY);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_L_lodi(
  amrex::IntVect const& iv, int dir,
  amrex::Array4<const amrex::Real> const& q, amrex::Real drho,
  amrex::Real dun, amrex::Real dut1, amrex::Real dut2, amrex::Real dT,
  amrex::Real const* dY, amrex::Real* L, Parm const& parm)
{
  using amrex::Real;
  const Real rho = amrex::max(q(iv, cls_t::QRHO), parm.min_rho);
  const Real T = amrex::max(q(iv, cls_t::QT), parm.min_T);
  const Real p = amrex::max(q(iv, cls_t::QPRES), parm.min_p);
  const Real c = amrex::max(q(iv, cls_t::QC), Real(1.e-14));
  const Real un = q(iv, qvel<cls_t>(dir));
  const Real gamma = amrex::max(q(iv, cls_t::QG), Real(1.0) + Real(1.e-12));

  const Real lam_m = un - c;
  const Real lam_0 = un;
  const Real lam_p = un + c;

  L[LMINUS] = lam_m *
    (drho / (Real(2.0) * gamma) - rho * dun / (Real(2.0) * c)
     + rho * dT / (Real(2.0) * gamma * T));

  // The 1/gamma on the density term is required for the forward and inverse
  // transforms to recover rho_t + u rho_n + rho u_n = 0 exactly.
  L[LENT] = lam_0 *
    (-(gamma - Real(1.0)) * T * drho / (gamma * rho) + dT / gamma);
  L[LTAN1] = lam_0 * dut1;
  L[LTAN2] = lam_0 * dut2;
  L[LPLUS] = lam_p *
    (drho / (Real(2.0) * gamma) + rho * dun / (Real(2.0) * c)
     + rho * dT / (Real(2.0) * gamma * T));

#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) {
    L[LSP + n] = lam_0 * dY[n];
  }
#else
  amrex::ignore_unused(dY, p);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void transverse_rhs(
  amrex::IntVect const& iv, int dir,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
  amrex::Array4<const amrex::Real> const& q, amrex::Real& drhodt,
  amrex::Real& dudt, amrex::Real& dvdt, amrex::Real& dwdt,
  amrex::Real& dTdt, amrex::Real* dYdt, Parm const& parm)
{
  using amrex::Real;
  drhodt = Real(0.0);
  dudt = Real(0.0);
  dvdt = Real(0.0);
  dwdt = Real(0.0);
  dTdt = Real(0.0);
#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) dYdt[n] = Real(0.0);
#else
  amrex::ignore_unused(dYdt);
#endif

  if (!parm.use_transverse || AMREX_SPACEDIM == 1) return;

  const Real rho = amrex::max(q(iv, cls_t::QRHO), parm.min_rho);
  const Real T = amrex::max(q(iv, cls_t::QT), parm.min_T);
  const Real gamma = amrex::max(q(iv, cls_t::QG), Real(1.0) + Real(1.e-12));
  Real vel[3] = {q(iv, cls_t::QU), q(iv, cls_t::QV), q(iv, cls_t::QW)};

  auto dc = [&] AMREX_GPU_DEVICE (int tdir, int n) noexcept -> Real {
    const auto e = amrex::IntVect::TheDimensionVector(tdir);
    return Real(0.5) * dxinv[tdir] * (q(iv + e, n) - q(iv - e, n));
  };

  for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
    if (tdir == dir) continue;
    const Real ut = vel[tdir];
    const Real div_t = dc(tdir, qvel<cls_t>(tdir));

    drhodt -= ut * dc(tdir, cls_t::QRHO) + rho * div_t;
    dTdt -= ut * dc(tdir, cls_t::QT)
             + (gamma - Real(1.0)) * T * div_t;
    dudt -= ut * dc(tdir, cls_t::QU);
    dvdt -= ut * dc(tdir, cls_t::QV);
    dwdt -= ut * dc(tdir, cls_t::QW);

    if (tdir == 0) dudt -= dc(tdir, cls_t::QPRES) / rho;
    if (tdir == 1) dvdt -= dc(tdir, cls_t::QPRES) / rho;
    if (tdir == 2) dwdt -= dc(tdir, cls_t::QPRES) / rho;

#if (NUM_SPECIES > 1)
    for (int n = 0; n < NUM_SPECIES; ++n) {
      dYdt[n] -= ut * dc(tdir, cls_t::QFS + n);
    }
#endif
  }
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void transverse_acoustic_amplitudes(
  amrex::IntVect const& iv, int dir,
  amrex::Array4<const amrex::Real> const& q, amrex::Real drhodt,
  amrex::Real dudt, amrex::Real dvdt, amrex::Real dwdt,
  amrex::Real dTdt, amrex::Real& Ltminus, amrex::Real& Ltplus,
  Parm const& parm)
{
  using amrex::Real;
  const Real rho = amrex::max(q(iv, cls_t::QRHO), parm.min_rho);
  const Real T = amrex::max(q(iv, cls_t::QT), parm.min_T);
  const Real p = amrex::max(q(iv, cls_t::QPRES), parm.min_p);
  const Real c = amrex::max(q(iv, cls_t::QC), Real(1.e-14));
  const Real vel_rhs[3] = {dudt, dvdt, dwdt};
  const Real p_t = p * (drhodt / rho + dTdt / T);
  const Real un_t = vel_rhs[dir];

  // Equivalent characteristic amplitudes L_T defined by q_t,T = -S L_T.
  Ltminus = -p_t / (Real(2.0) * c * c) + rho * un_t / (Real(2.0) * c);
  Ltplus = -p_t / (Real(2.0) * c * c) - rho * un_t / (Real(2.0) * c);
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void apply_outflow(
  amrex::IntVect const& iv, int dir, int side_sign,
  amrex::Array4<const amrex::Real> const& q, amrex::Real* L,
  amrex::Real Ltminus, amrex::Real Ltplus, int type, Parm const& parm)
{
  using amrex::Real;
  const Real c = amrex::max(q(iv, cls_t::QC), Real(1.e-14));
  const Real un = q(iv, qvel<cls_t>(dir));
  const Real un_out = -Real(side_sign) * un;

  // All characteristics leave at a supersonic outflow.
  if (un_out >= c) return;

  const Real beta = parm.use_transverse ? parm.transverse_relax : Real(0.0);
  Real target = Real(0.0);
  if (type == 3) {
    const Real p = q(iv, cls_t::QPRES);
    target = parm.sigma * amrex::max(Real(1.0) - parm.Mmax * parm.Mmax,
                                    Real(0.0))
             * (p - parm.Ptarget) / (Real(2.0) * c * parm.Lchar);
  }

  if (side_sign > 0) {
    L[LPLUS] = target - beta * Ltplus;
  } else {
    L[LMINUS] = target - beta * Ltminus;
  }
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void apply_inflow(
  amrex::IntVect const& iv, int dir, int side_sign,
  amrex::Array4<const amrex::Real> const& q, amrex::Real* L,
  amrex::Real Ltminus, amrex::Real Ltplus, Parm const& parm)
{
  using amrex::Real;
  const Real rho = amrex::max(q(iv, cls_t::QRHO), parm.min_rho);
  const Real c = amrex::max(q(iv, cls_t::QC), Real(1.e-14));
  const Real T = q(iv, cls_t::QT);
  const Real target_vel[3] = {parm.utarget, parm.vtarget, parm.wtarget};
  int t1, t2;
  tangent_dirs(dir, t1, t2);

  L[LENT] = parm.eta * (T - parm.Ttarget);
  if (t1 >= 0) L[LTAN1] = parm.eta * (q(iv, qvel<cls_t>(t1)) - target_vel[t1]);
  if (t2 >= 0) L[LTAN2] = parm.eta * (q(iv, qvel<cls_t>(t2)) - target_vel[t2]);

  const Real acoustic_delta =
    rho * parm.eta * (q(iv, qvel<cls_t>(dir)) - target_vel[dir]) / c;
  const Real beta = parm.use_transverse ? parm.transverse_relax : Real(0.0);
  if (side_sign > 0) {
    L[LPLUS] = L[LMINUS] + acoustic_delta - beta * Ltplus;
  } else {
    L[LMINUS] = L[LPLUS] - acoustic_delta - beta * Ltminus;
  }

#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) L[LSP + n] = Real(0.0);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void primitive_rhs_from_L(
  amrex::IntVect const& iv, int dir,
  amrex::Array4<const amrex::Real> const& q, amrex::Real const* L,
  amrex::Real& drhodt, amrex::Real& dudt, amrex::Real& dvdt,
  amrex::Real& dwdt, amrex::Real& dTdt, amrex::Real* dYdt,
  Parm const& parm)
{
  using amrex::Real;
  const Real rho = amrex::max(q(iv, cls_t::QRHO), parm.min_rho);
  const Real T = amrex::max(q(iv, cls_t::QT), parm.min_T);
  const Real c = amrex::max(q(iv, cls_t::QC), Real(1.e-14));
  const Real gamma = amrex::max(q(iv, cls_t::QG), Real(1.0) + Real(1.e-12));

  drhodt = -L[LMINUS] - L[LPLUS] + rho * L[LENT] / T;
  dTdt = -(gamma - Real(1.0)) * T * (L[LMINUS] + L[LPLUS]) / rho
          - L[LENT];

  Real vel_rhs[3] = {Real(0.0), Real(0.0), Real(0.0)};
  vel_rhs[dir] = -c * (L[LPLUS] - L[LMINUS]) / rho;
  int t1, t2;
  tangent_dirs(dir, t1, t2);
  if (t1 >= 0) vel_rhs[t1] = -L[LTAN1];
  if (t2 >= 0) vel_rhs[t2] = -L[LTAN2];
  dudt = vel_rhs[0];
  dvdt = vel_rhs[1];
  dwdt = vel_rhs[2];

#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) dYdt[n] = -L[LSP + n];
#else
  amrex::ignore_unused(dYdt);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_lodi_rhs_to_cons(
  amrex::IntVect const& iv, int dir, int side_sign, int boundary_index,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& rhs, int type, Parm const& parm)
{
  using amrex::Real;
  Real drho, dun, dut1, dut2, dT;
  Real dY[NUM_SPECIES] = {Real(0.0)};
  one_sided_deriv_prim<cls_t>(iv, dir, side_sign, dxinv[dir], q,
                              drho, dun, dut1, dut2, dT, dY,
                              parm.derivative_order);

  Real L[5 + NUM_SPECIES] = {Real(0.0)};
  compute_L_lodi<cls_t>(iv, dir, q, drho, dun, dut1, dut2, dT, dY, L, parm);

  // Persistent physical ghosts are duplicated between tangentially adjacent
  // AMReX FABs; those duplicated values cannot be synchronized by the normal
  // MultiFab FillBoundary operation because they are outside every valid box.
  // Evaluate transverse derivatives on the adjacent valid boundary plane,
  // whose FAB ghosts are synchronized, and apply that compatibility forcing
  // to each normal ghost layer.  This removes box-decomposition dependence.
  amrex::IntVect transverse_iv = iv;
  transverse_iv[dir] = boundary_index;

  Real trho, tu, tv, tw, tT;
  Real tY[NUM_SPECIES] = {Real(0.0)};
  transverse_rhs<cls_t>(transverse_iv, dir, dxinv, q, trho, tu, tv, tw, tT,
                        tY, parm);
  Real Ltminus = Real(0.0);
  Real Ltplus = Real(0.0);
  transverse_acoustic_amplitudes<cls_t>(iv, dir, q, trho, tu, tv, tw, tT,
                                         Ltminus, Ltplus, parm);

  const Real un_out = -Real(side_sign) * q(iv, qvel<cls_t>(dir));
  if (type == 1 || un_out < Real(0.0)) {
    apply_inflow<cls_t>(iv, dir, side_sign, q, L, Ltminus, Ltplus, parm);
  } else {
    apply_outflow<cls_t>(iv, dir, side_sign, q, L, Ltminus, Ltplus,
                         type, parm);
  }

  Real rho_t, u_t, v_t, w_t, T_t;
  Real Y_t[NUM_SPECIES] = {Real(0.0)};
  primitive_rhs_from_L<cls_t>(iv, dir, q, L, rho_t, u_t, v_t, w_t, T_t,
                              Y_t, parm);
  rho_t += trho;
  u_t += tu;
  v_t += tv;
  w_t += tw;
  T_t += tT;
#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) Y_t[n] += tY[n];
#endif

  const Real rho = amrex::max(q(iv, cls_t::QRHO), parm.min_rho);
  const Real T = amrex::max(q(iv, cls_t::QT), parm.min_T);
  const Real p = amrex::max(q(iv, cls_t::QPRES), parm.min_p);
  const Real gamma = amrex::max(q(iv, cls_t::QG), Real(1.0) + Real(1.e-12));
  const Real u = q(iv, cls_t::QU);
  const Real v = q(iv, cls_t::QV);
  const Real w = q(iv, cls_t::QW);

#if (NUM_SPECIES > 1)
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const Real Y = q(iv, cls_t::QFS + n);
    rhs(iv, cls_t::UFS + n) = Y * rho_t + rho * Y_t[n];
  }
#else
  rhs(iv, cls_t::URHO) = rho_t;
#endif
  rhs(iv, cls_t::UMX) = u * rho_t + rho * u_t;
  rhs(iv, cls_t::UMY) = v * rho_t + rho * v_t;
  rhs(iv, cls_t::UMZ) = w * rho_t + rho * w_t;

  const Real p_t = p * (rho_t / rho + T_t / T);
  rhs(iv, cls_t::UET) = p_t / (gamma - Real(1.0))
    + Real(0.5) * rho_t * (u * u + v * v + w * w)
    + rho * (u * u_t + v * v_t + w * w_t);
}

} // namespace nscbc

#endif
