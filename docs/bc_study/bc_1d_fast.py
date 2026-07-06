#!/usr/bin/env python3
# =============================================================================
#  Vectorised 1-D far-field BC reflection benchmark (numpy, ~100x faster).
#  Physics identical to bc_1d_reflection.py:
#    1-D Euler, HLLC (Roe-averaged speeds, mirrors src/rhs/Riemann.h:145),
#    1st-order Godunov, fine grid.
#  BCs: FOE (foextrap), DIR (hard Dirichlet), NSC (Cerisse algebraic NSCBC,
#       bc_types.h:180-291), LODI (true Poinsot-Lele NSCBC, relaxation sigma).
#  Reflection coefficient is energy-based (Mach-robust).
# =============================================================================
import numpy as np
GAM = 1.4

def cons2prim(U):
    rho = U[0]; u = U[1] / rho
    p = (GAM - 1.0) * (U[2] - 0.5 * rho * u * u)
    return rho, u, p

def prim2cons_arr(rho, u, p):
    return np.array([rho, rho * u, p / (GAM - 1.0) + 0.5 * rho * u * u])

def csound(rho, p):
    return np.sqrt(GAM * p / rho)

def hllc_vec(UL, UR):
    """Vectorised HLLC over all faces. UL,UR shape (3, Nf). Roe signal speeds."""
    rL, uL, pL = cons2prim(UL); rR, uR, pR = cons2prim(UR)
    pL = np.maximum(pL, 1e-12); pR = np.maximum(pR, 1e-12)
    cL = csound(rL, pL); cR = csound(rR, pR)
    rp = np.sqrt(rR / rL)
    uroe = (uL + uR * rp) / (1.0 + rp)
    croe = (cL + cR * rp) / (1.0 + rp)
    SL = np.minimum(uL - cL, uroe - croe)
    SR = np.maximum(uR + cR, uroe + croe)
    EL = UL[2]; ER = UR[2]
    FL = np.array([rL * uL, rL * uL * uL + pL, uL * (EL + pL)])
    FR = np.array([rR * uR, rR * uR * uR + pR, uR * (ER + pR)])
    Sstar = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / \
            (rL * (SL - uL) - rR * (SR - uR))
    def star(U, S, u, r, p):
        E = U[2] / r
        fac = r * (S - u) / (S - Sstar)
        return fac * np.array([np.ones_like(Sstar), Sstar,
                               E + (Sstar - u) * (Sstar + p / (r * (S - u)))])
    UsL = star(UL, SL, uL, rL, pL)
    UsR = star(UR, SR, uR, rR, pR)
    FsL = FL + SL * (UsL - UL)
    FsR = FR + SR * (UsR - UR)
    F = np.empty_like(FL)
    for k in range(3):
        F[k] = np.where(SL >= 0, FL[k],
                np.where(Sstar >= 0, FsL[k],
                np.where(SR > 0, FsR[k], FR[k])))
    return F

# ----- ghost states (right boundary, outward normal +x) ----------------------
def ghost_nscbc(U_i, Uinf):
    rho_i, u_i, p_i = cons2prim(U_i.reshape(3, 1))
    rho_i = rho_i[0]; u_i = u_i[0]; p_i = max(p_i[0], 1e-12)
    rinf, uinf, pinf = cons2prim(Uinf.reshape(3, 1))
    rinf = rinf[0]; uinf = uinf[0]; pinf = pinf[0]
    ci = csound(rho_i, p_i); gm1 = GAM - 1.0; eps = 1e-14
    uni = u_i; unf = uinf; cf = csound(rinf, pinf)
    if uni <= -ci:
        return prim2cons_arr(np.array([rinf]), np.array([uinf]), np.array([pinf]))[:, 0]
    if uni >= ci:
        return U_i.copy()
    if uni < 0.0:
        jplus_i = uni + 2.0 * ci / gm1; jminus_f = unf - 2.0 * cf / gm1
        cb = max(0.25 * gm1 * (jplus_i - jminus_f), eps)
        sfar = pinf / max(rinf, eps) ** GAM
        rho = (cb * cb / (GAM * sfar)) ** (1.0 / gm1); p = rho * cb * cb / GAM
        return prim2cons_arr(np.array([rho]), np.array([uinf]), np.array([p]))[:, 0]
    jplus_i = uni + 2.0 * ci / gm1
    sint = p_i / max(rho_i, eps) ** GAM
    p = pinf; rho = (p / sint) ** (1.0 / GAM); cb = csound(rho, p)
    un = jplus_i - 2.0 * cb / gm1
    return prim2cons_arr(np.array([rho]), np.array([un]), np.array([p]))[:, 0]

def ghost_nscbc_fixed(U_i, Uinf):
    """Independent re-implementation of the OTHER AI's fixed bc_types.h:263-272.
    Subsonic outflow: J+ from interior, J- from freestream, entropy from interior
    (Thompson / Whitfield-Janus characteristic non-reflecting outflow)."""
    rho_i, u_i, p_i = cons2prim(U_i.reshape(3, 1))
    rho_i = max(rho_i[0], 1e-14); u_i = u_i[0]; p_i = max(p_i[0], 1e-14)
    rinf, uinf, pinf = cons2prim(Uinf.reshape(3, 1))
    rinf = max(rinf[0], 1e-14); uinf = uinf[0]; pinf = max(pinf[0], 1e-14)
    ci = csound(rho_i, p_i); cf = csound(rinf, pinf); gm1 = GAM - 1.0; eps = 1e-14
    uni = u_i; unf = uinf
    if uni <= -ci:                                       # supersonic inflow
        return prim2cons_arr(np.array([rinf]), np.array([uinf]), np.array([pinf]))[:, 0]
    if uni >= ci:                                        # supersonic outflow
        return U_i.copy()
    jplus_i = uni + 2.0 * ci / gm1
    jminus_f = unf - 2.0 * cf / gm1
    cb = max(0.25 * gm1 * (jplus_i - jminus_f), eps)
    if uni < 0.0:                                        # subsonic inflow
        sfar = pinf / rinf ** GAM
        un = 0.5 * (jplus_i + jminus_f)
        rho = (cb * cb / (GAM * sfar)) ** (1.0 / gm1); p = rho * cb * cb / GAM
    else:                                                # subsonic outflow (FIXED)
        sint = p_i / rho_i ** GAM
        un = 0.5 * (jplus_i + jminus_f)
        rho = (cb * cb / (GAM * sint)) ** (1.0 / gm1); p = rho * cb * cb / GAM
    return prim2cons_arr(np.array([rho]), np.array([un]), np.array([p]))[:, 0]

def lodi_rhs_right(U, dx, Uinf, sigma, Ldom):
    rho, u, p = cons2prim(U[:, -1].reshape(3, 1)); rho = rho[0]; u = u[0]; p = p[0]
    rho1, u1, p1 = cons2prim(U[:, -2].reshape(3, 1)); rho1 = rho1[0]; u1 = u1[0]; p1 = p1[0]
    c = csound(rho, p)
    dpdx = (p - p1) / dx; dudx = (u - u1) / dx; drdx = (rho - rho1) / dx
    rinf, uinf, pinf = cons2prim(Uinf.reshape(3, 1)); pinf = pinf[0]
    L5 = (u + c) * (dpdx + rho * c * dudx)
    L2 = u * (c * c * drdx - dpdx)
    M = abs(u) / c
    K = sigma * (1.0 - M * M) * c / Ldom
    L1 = K * (p - pinf)
    d_rho = -(L2 + 0.5 * (L5 + L1)) / (c * c)
    d_u = -(L5 - L1) / (2.0 * rho * c)
    d_p = -0.5 * (L5 + L1)
    dE = d_p / (GAM - 1.0) + 0.5 * d_rho * u * u + rho * u * d_u
    return np.array([d_rho, d_rho * u + rho * d_u, dE])

def step(U, dt, dx, bc, Uinf, sigma, Ldom):
    N = U.shape[1]
    if bc == "LODI":
        Uext = np.column_stack([U[:, 0], U])
        Fl = hllc_vec(Uext[:, :-1], Uext[:, 1:])   # left faces of cells 0..N-1
        Unew = U.copy()
        Fr = hllc_vec(U[:, :-1], U[:, 1:])          # interior faces (right of cells 0..N-2)
        Unew[:, :-1] = U[:, :-1] - dt / dx * (Fr - Fl[:, :-1])
        Unew[:, -1] = U[:, -1] + dt * lodi_rhs_right(U, dx, Uinf, sigma, Ldom)
        return Unew
    if bc == "FOE":
        Ug = U[:, -1].copy()
    elif bc == "DIR":
        Ug = Uinf.copy()
    elif bc == "NSC":
        Ug = ghost_nscbc(U[:, -1], Uinf)
    elif bc == "NSC_FIX":
        Ug = ghost_nscbc_fixed(U[:, -1], Uinf)
    Uext = np.column_stack([U[:, 0], U, Ug])
    F = hllc_vec(Uext[:, :-1], Uext[:, 1:])
    return U - dt / dx * (F[:, 1:] - F[:, :-1])

def acoustic_energy(U, dx):
    rho, u, p = cons2prim(U)
    r0 = rho[0]; p0 = p[0]; u0 = u[0]; c0 = csound(r0, p0)
    e = (p - p0) ** 2 / (2 * r0 * c0 ** 2) + 0.5 * r0 * (u - u0) ** 2
    return np.sum(e) * dx

def run_pulse(bc, M=0.0, amp=1e-3, sigma=0.25, L=10.0, N=1000, x0=3.0,
              width=0.25, t_end=13.0, sensor_x=6.0, record_xt=False, kind="acoustic"):
    r0 = 1.0; p0 = 1.0 / GAM; c0 = csound(r0, p0); u0 = M * c0
    x = np.linspace(0.5 * L / N, L - 0.5 * L / N, N); dx = L / N
    if kind == "acoustic":
        dp = amp * np.exp(-((x - x0) / width) ** 2)
        rho = r0 + dp / c0 ** 2; u = u0 + dp / (r0 * c0); p = p0 + dp
    else:  # entropy bump: density only
        drho = amp * np.exp(-((x - x0) / width) ** 2)
        rho = r0 + drho; u = np.full(N, u0); p = np.full(N, p0)
    U = prim2cons_arr(rho, u, p)
    Uinf = prim2cons_arr(np.array([r0]), np.array([u0]), np.array([p0]))[:, 0]
    dt = 0.4 * dx / (c0 + abs(u0))
    sensor_i = int(sensor_x / dx)
    ts, ps, xt, xtt = [], [], [], []
    E0 = max(acoustic_energy(U, dx), 1e-300)
    t_star = 1.3 * (L - x0) / (c0 + abs(u0))
    E_star = None
    t = 0.0; n = 0
    while t < t_end:
        if t + dt > t_end: dt = t_end - t
        U = step(U, dt, dx, bc, Uinf, sigma, L); t += dt; n += 1
        _, _, pc = cons2prim(U[:, sensor_i:sensor_i+1]); ps.append(pc[0] - p0); ts.append(t)
        if E_star is None and t >= t_star:
            E_star = acoustic_energy(U, dx)
        if record_xt and n % 4 == 0:
            _, _, pf = cons2prim(U); xt.append(pf - p0); xtt.append(t)
    if E_star is None: E_star = acoustic_energy(U, dx)
    out = {"t": np.array(ts), "p_sensor": np.array(ps), "x": x, "p0": p0, "amp": amp,
           "R_energy": float(np.sqrt(max(E_star, 0.0) / E0))}
    if record_xt: out["xt"] = np.array(xt); out["xt_t"] = np.array(xtt)
    return out

def reflection_coeff(res):
    return res["R_energy"]
