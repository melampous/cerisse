#!/usr/bin/env python3
# =============================================================================
#  1-D far-field boundary-condition reflection benchmark
#
#  Method follows the canonical acoustic-pulse test of
#      Poinsot & Lele, J. Comput. Phys. 101 (1992) 104-129
#      Rudy & Strikwerda, J. Comput. Phys. 36 (1980) 55-70
#
#  A right-running acoustic pulse is launched in a uniform mean flow and
#  propagates toward the outflow (right) boundary.  The amplitude returned
#  into the domain is the reflection coefficient
#          R = max |delta p_reflected| / max |delta p_incident|.
#
#  BC implementations under test (all in one finite-volume + HLLC code,
#  the same Riemann solver Cerisse uses, src/rhs/Riemann.h:145):
#    FOE  : foextrap                 (cns.lo_bc/hi_bc = 2)  ghost = interior
#    DIR  : hard Dirichlet freestream (CNS_DIRICHLET_FARFIELD) ghost = U_inf
#    NSC  : Cerisse "nscbc" algebraic (cns.lo_bc/hi_bc = 7)  bc_types.h:180-291
#    LODI : TRUE Poinsot-Lele NSCBC with relaxation sigma (reference method)
#
#  The LODI path is the academically-correct non-reflecting BC and is used to
#  prove that the algebraic "NSC" code is the sigma -> infinity (perfectly
#  reflecting) limit of a proper NSCBC.
# =============================================================================
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAM = 1.4

# -----------------------------------------------------------------------------
#  thermodynamics / Euler helpers
# -----------------------------------------------------------------------------
def cons2prim(U):
    rho = U[0]
    u = U[1] / rho
    E = U[2] / rho
    p = (GAM - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

def prim2cons(rho, u, p):
    E = p / ((GAM - 1.0) * rho) + 0.5 * u * u
    return np.array([rho, rho * u, rho * E])

def csound(rho, p):
    return np.sqrt(GAM * p / rho)

def euler_flux(U):
    rho, u, p = cons2prim(U)
    E = U[2]
    return np.array([rho * u, rho * u * u + p, u * (E + p)])

def hllc(UL, UR):
    """HLLC with Roe-averaged signal speeds — mirrors Cerisse hllc()."""
    rL, uL, pL = cons2prim(UL)
    rR, uR, pR = cons2prim(UR)
    cL = csound(rL, pL)
    cR = csound(rR, pR)
    rp = np.sqrt(rR / rL)
    uroe = (uL + uR * rp) / (1.0 + rp)
    croe = (cL + cR * rp) / (1.0 + rp)
    SL = min(uL - cL, uroe - croe)
    SR = max(uR + cR, uroe + croe)
    if SL >= 0.0:
        return euler_flux(UL)
    if SR <= 0.0:
        return euler_flux(UR)
    Sstar = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / \
            (rL * (SL - uL) - rR * (SR - uR))
    def star(U, S, u, r, p):
        E = U[2] / r
        fac = r * (S - u) / (S - Sstar)
        return fac * np.array([1.0, Sstar,
                               E + (Sstar - u) * (Sstar + p / (r * (S - u)))])
    if Sstar >= 0.0:
        return euler_flux(UL) + SL * (star(UL, SL, uL, rL, pL) - UL)
    else:
        return euler_flux(UR) + SR * (star(UR, SR, uR, rR, pR) - UR)

# -----------------------------------------------------------------------------
#  ghost-cell BC fillers : given first interior cell U_i, return ghost U_g
#  (right boundary, outward normal +x).  Faithful to src/set/bcs.cpp:38-88
#  which fills every ghost layer with the same bcnormal() output.
# -----------------------------------------------------------------------------
def ghost_foextrap(U_i, Uinf):
    return U_i.copy()

def ghost_dirichlet(U_i, Uinf):
    return Uinf.copy()

def ghost_nscbc_algebraic(U_i, Uinf):
    """Direct translation of GlobalBC::bc_nscbc_farfield, bc_types.h:180-291,
    right boundary nx=+1."""
    rho_i, u_i, p_i = cons2prim(U_i)
    rinf, uinf, pinf = cons2prim(Uinf)
    ci = csound(rho_i, p_i)
    gm1 = GAM - 1.0
    eps = 1e-14
    nx = 1.0
    uni = u_i * nx
    unf = uinf * nx
    cf = csound(rinf, pinf)
    if uni <= -ci:                                  # Case 1 supersonic inflow
        return prim2cons(rinf, uinf, pinf)
    if uni >= ci:                                   # Case 2 supersonic outflow
        return U_i.copy()
    if uni < 0.0:                                   # Case 3 subsonic inflow
        jplus_i = uni + 2.0 * ci / gm1
        jminus_f = unf - 2.0 * cf / gm1
        cb = max(0.25 * gm1 * (jplus_i - jminus_f), eps)
        sfar = pinf / max(rinf, eps) ** GAM
        un = 0.5 * (jplus_i + jminus_f)
        rho = (cb * cb / (GAM * sfar)) ** (1.0 / gm1)
        p = rho * cb * cb / GAM
        u = uinf                                    # tangential=0 in 1D
        return prim2cons(rho, u, p)
    # Case 4 subsonic outflow  (THE reflective branch: p hard-set to p_inf)
    jplus_i = uni + 2.0 * ci / gm1
    sint = p_i / max(rho_i, eps) ** GAM
    p = pinf                                        # <-- bc_types.h:266
    rho = (p / sint) ** (1.0 / GAM)
    cb = csound(rho, p)
    un = jplus_i - 2.0 * cb / gm1
    return prim2cons(rho, un, p)

GHOST = {"FOE": ghost_foextrap, "DIR": ghost_dirichlet, "NSC": ghost_nscbc_algebraic}

# -----------------------------------------------------------------------------
#  finite-volume driver (1st-order Godunov + HLLC, fine grid).
#  For ghost BCs we append one ghost cell on the right and run HLLC at the face.
#  For the LODI reference BC the right boundary cell RHS is replaced by the
#  characteristic (Poinsot-Lele) time-derivative with relaxation sigma.
# -----------------------------------------------------------------------------
def lodi_rhs_right(U, dx, Uinf, sigma, Ldom):
    """Return d/dt of conservative U in the right boundary cell using LODI.
    Outgoing waves (L5 acoustic, L2 entropy) from one-sided interior diff;
    incoming wave L1 set by relaxation toward p_inf."""
    rho, u, p = cons2prim(U[:, -1])
    rho1, u1, p1 = cons2prim(U[:, -2])
    c = csound(rho, p)
    # one-sided (backward) derivatives at boundary
    dpdx = (p - p1) / dx
    dudx = (u - u1) / dx
    drdx = (rho - rho1) / dx
    # characteristic wave amplitudes (Poinsot-Lele 1992, eq 5)
    L5 = (u + c) * (dpdx + rho * c * dudx)        # outgoing (u+c>0)
    L2 = u * (c * c * drdx - dpdx)                # entropy (u>0)
    rinf, uinf, pinf = cons2prim(Uinf)
    M = abs(u) / c
    K = sigma * (1.0 - M * M) * c / Ldom          # relaxation coefficient
    L1 = K * (p - pinf)                           # incoming, relaxed
    # LODI primitive time derivatives (Poinsot-Lele eq 37-40, 1D)
    d_rho = -(L2 + 0.5 * (L5 + L1)) / (c * c)
    d_u = -(L5 - L1) / (2.0 * rho * c)
    d_p = -0.5 * (L5 + L1)
    # convert to conservative d/dt
    dU = np.zeros(3)
    dU[0] = d_rho
    dU[1] = d_rho * u + rho * d_u
    E = p / ((GAM - 1.0)) + 0.5 * rho * u * u
    dE = d_p / (GAM - 1.0) + 0.5 * d_rho * u * u + rho * u * d_u
    dU[2] = dE
    return dU

def step(U, dt, dx, bc, Uinf, sigma=0.25, Ldom=10.0):
    """One forward-Euler Godunov step. Left boundary: foextrap (untested side)."""
    N = U.shape[1]
    if bc == "LODI":
        # interior faces 1..N-1 via HLLC with foextrap-extended left ghost
        Uext = np.column_stack([U[:, 0], U])      # left ghost = copy
        F = np.zeros((3, N))                       # faces 0..N-1 (left of each cell)
        for i in range(N):
            F[:, i] = hllc(Uext[:, i], Uext[:, i + 1])
        # right face of last cell handled by LODI (no HLLC flux there)
        Unew = U.copy()
        # cells 0..N-2 standard: dU = -(F_{i+1}-F_i)/dx ; need right face fluxes
        Fright = np.zeros((3, N))
        for i in range(N - 1):
            Fright[:, i] = hllc(U[:, i], U[:, i + 1])
        # last cell right face: use physical flux of boundary state (LODI sets RHS)
        for i in range(N - 1):
            Unew[:, i] = U[:, i] - dt / dx * (Fright[:, i] - F[:, i])
        # boundary cell via LODI
        dU = lodi_rhs_right(U, dx, Uinf, sigma, Ldom)
        Unew[:, -1] = U[:, -1] + dt * dU
        return Unew
    # ghost-cell BCs
    gfn = GHOST[bc]
    Ug = gfn(U[:, -1], Uinf)
    Uext = np.column_stack([U[:, 0], U, Ug])      # left copy + right ghost
    F = np.zeros((3, N + 1))
    for i in range(N + 1):
        F[:, i] = hllc(Uext[:, i], Uext[:, i + 1])
    return U - dt / dx * (F[:, 1:] - F[:, :-1])

# -----------------------------------------------------------------------------
#  acoustic-pulse experiment
# -----------------------------------------------------------------------------
def acoustic_energy(U, dx):
    """Domain-integrated linear acoustic energy E = sum[ dp^2/(2 r c^2) + r du^2 /2 ]dx.
    Reference state taken as the domain edges (assumed at freestream)."""
    rho, u, p = cons2prim(U)
    r0 = rho[0]; p0 = p[0]; u0 = u[0]
    c0 = csound(r0, p0)
    dp = p - p0
    du = u - u0
    e = dp**2 / (2.0 * r0 * c0**2) + 0.5 * r0 * du**2
    return np.sum(e) * dx

def run_pulse(bc, M=0.0, amp=1e-3, sigma=0.25,
              L=10.0, N=1000, x0=3.0, width=0.25, t_end=11.0,
              sensor_x=6.0, record_xt=False):
    r0 = 1.0
    p0 = 1.0 / GAM            # c0 = 1
    c0 = csound(r0, p0)
    u0 = M * c0
    x = np.linspace(0.5 * L / N, L - 0.5 * L / N, N)
    dx = L / N
    # pure right-running acoustic pulse: dp, du=dp/(r c), drho=dp/c^2
    dp = amp * np.exp(-((x - x0) / width) ** 2)
    rho = r0 + dp / c0 ** 2
    u = u0 + dp / (r0 * c0)
    p = p0 + dp
    U = np.array([prim2cons(rho[i], u[i], p[i]) for i in range(N)]).T
    Uinf = prim2cons(r0, u0, p0)
    dt = 0.4 * dx / (c0 + abs(u0))
    sensor_i = int(sensor_x / dx)
    ts, ps = [], []
    xt = [] if record_xt else None
    xt_t = [] if record_xt else None
    # energy-based reflection: sample E right after the incident has exited
    E0 = acoustic_energy(U, dx)
    t_star = 1.3 * (L - x0) / (c0 + abs(u0))   # incident gone, reflected still in domain
    E_star = None
    t = 0.0
    nstep = 0
    while t < t_end:
        if t + dt > t_end:
            dt = t_end - t
        U = step(U, dt, dx, bc, Uinf, sigma=sigma, Ldom=L)
        t += dt
        nstep += 1
        _, _, pcur = cons2prim(U[:, sensor_i])
        ts.append(t)
        ps.append(pcur - p0)
        if E_star is None and t >= t_star:
            E_star = acoustic_energy(U, dx)
        if record_xt and nstep % 4 == 0:
            _, _, pfield = cons2prim(U)
            xt.append(pfield - p0)
            xt_t.append(t)
    if E_star is None:
        E_star = acoustic_energy(U, dx)
    out = {"t": np.array(ts), "p_sensor": np.array(ps), "x": x, "p0": p0,
           "c0": c0, "u0": u0, "amp": amp,
           "R_energy": float(np.sqrt(max(E_star, 0.0) / E0)) if E0 > 0 else 0.0}
    if record_xt:
        out["xt"] = np.array(xt)
        out["xt_t"] = np.array(xt_t)
    return out

def reflection_coeff(res):
    """Energy-based reflection coefficient (Mach-robust). Falls back to the
    sensor-amplitude ratio if energy was not recorded."""
    if "R_energy" in res:
        return res["R_energy"]
    t = res["t"]
    p = res["p_sensor"]
    half = len(t) // 2
    inc = np.max(np.abs(p[:half]))
    ref = np.max(np.abs(p[half:]))
    return ref / inc if inc > 0 else 0.0
