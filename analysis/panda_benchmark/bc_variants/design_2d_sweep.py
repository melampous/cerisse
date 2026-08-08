#!/usr/bin/env python3
"""Design table for the 2D axisymmetric inlet-profile sweep whose single
objective is the law x_fc = F(inlet profile).

For every candidate value of delta_jet (and the two functional forms coded
in prob.h) this computes the inlet scalars that are candidate predictors:

  velocity-based (what a textbook shear layer quotes)
      delta_w  = dU / max|du/dr|              vorticity thickness
      theta_i  = int f (1-f) dr               incompressible momentum thickness
      dstar_i  = int (1-f) dr                 incompressible displacement thickness
  compressible / flux-based (what the solver actually sees, since the BC
  holds p = p_e uniform and T0 = T_amb, so rho varies across the layer)
      dstar_c  = int (1 - rho u / rho_e U_e) dr
      theta_c  = int (rho u / rho_e U_e)(1 - f) dr
  integral (exact in RZ, no thin-layer assumption)
      C_d      = mdot / mdot_ideal = A_eff / A       discharge coefficient
      C_J      = Jx  / Jx_ideal                      momentum-flux coefficient
      D_eff    = D sqrt(C_d)

and the resolution actually delivered by the production 2D grid
(finest dx = 24.8 um in the lip corridor x <= 1 D_e).
"""
import numpy as np

# NumPy 2.0 renamed ``trapz`` to ``trapezoid``.  Keep the analysis script
# reproducible with both the thesis environment and newer installations.
# getattr with np.trapz as the default fails on NumPy 2: the default is
# evaluated eagerly and trapz no longer exists there.
_trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

GAM, RGAS = 1.4, 8.31446261815324 / 28.96e-3
CP = GAM * RGAS / (GAM - 1.0)
T_AMB = 297.15
T_E = T_AMB / 1.2                      # 247.625 K sonic exit static T
R = 0.0127
D = 2.0 * R
U_E = np.sqrt(GAM * RGAS * T_E)        # sonic exit velocity
DX_LIP = 1.5875e-3 / 64                # 24.8 um : level 6 in the lip corridor


def fvel(r, a, form, s=3.0):
    """Velocity shape u/U_e exactly as coded in prob.h bcnormal.

    'shifted' is the lip-centred tanh displaced INTO the orifice by s*a,
    i.e. r0 = R - s*a.  s = 0 recovers the plain lip-centred tanh; the
    production cases use s = 3.  'wall' is the wall-attached tanh.
    """
    if a <= 0.0:
        return np.where(r <= R, 1.0, 0.0)
    if form == "shifted":
        r0 = R - s * a
        return 0.5 * (1.0 - np.tanh((r - r0) / a))
    if form == "wall":
        return np.tanh(np.maximum(R - r, 0.0) / a)
    raise ValueError(form)


def scalars(a, form, s=3.0):
    n = 400001
    r = np.linspace(0.0, R, n)
    f = fvel(r, a, form, s)
    u = U_E * f
    T = T_AMB - 0.5 * u * u / CP        # constant total temperature
    rho = 1.0 / (RGAS * T)              # p_e cancels in every ratio below
    rho_e = 1.0 / (RGAS * T_E)
    g = rho * u / (rho_e * U_E)         # normalised mass flux

    # vorticity thickness
    if a <= 0.0:
        dw = 0.0
    else:
        dw = 1.0 / np.max(np.abs(np.gradient(f, r)))

    dstar_i = _trapezoid(1.0 - f, r)
    theta_i = _trapezoid(f * (1.0 - f), r)
    dstar_c = _trapezoid(1.0 - g, r)
    theta_c = _trapezoid(g * (1.0 - f), r)

    # exact RZ integrals
    Cd = 2.0 * _trapezoid(g * r, r) / R**2
    Jn = rho * u * u
    Jd = rho_e * U_E**2
    CJ = 2.0 * _trapezoid((Jn / Jd) * r, r) / R**2
    leak = float(fvel(np.array([R - 1e-12]), a, form, s)[0])
    return dict(dw=dw, theta_i=theta_i, dstar_i=dstar_i,
                theta_c=theta_c, dstar_c=dstar_c, Cd=Cd, CJ=CJ,
                H=(dstar_i / theta_i if theta_i > 0 else np.nan),
                leak=leak, Deff=D * np.sqrt(Cd))


# ---------------------------------------------------------------- existing
EXIST = [("Baseline      (3D+2D)", 254.0, "shifted", 0.883),
         ("Thin shifted  (3D)   ", 100.0, "shifted", 0.922),
         ("Wall tanh     (3D)   ", 200.0, "wall",    0.930),
         ("Top hat       (3D)   ",   0.0, "shifted", 0.945)]

print("U_e = %.2f m/s   T_e = %.3f K   D_e = %.1f mm   dx_lip = %.1f um\n"
      % (U_E, T_E, D * 1e3, DX_LIP * 1e6))
hdr = ("%-22s %6s %7s %7s %7s %7s %7s %7s %7s %6s"
       % ("case", "a[um]", "dw[um]", "th_i", "d*_i", "th_c", "d*_c",
          "C_d", "Deff/D", "cells"))
print("EXISTING FOUR (x_fc measured in 3D)")
print(hdr)
rows = []
for lab, a, form, xfc in EXIST:
    s = scalars(a * 1e-6, form)
    cells = s["dw"] / DX_LIP
    rows.append((lab, a, form, xfc, s))
    print("%-22s %6.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.4f %7.4f %6.1f"
          % (lab, a, s["dw"] * 1e6, s["theta_i"] * 1e6, s["dstar_i"] * 1e6,
             s["theta_c"] * 1e6, s["dstar_c"] * 1e6, s["Cd"],
             s["Deff"] / D, cells))

# which predictor best explains the four measured x_fc
print("\nPREDICTOR RANKING on the existing four (x_fc vs P, linear, n=4)")
y = np.array([r[3] for r in rows])
for key, nm in [("dstar_c", "d*_c  (compressible)"), ("dstar_i", "d*_i"),
                ("theta_c", "th_c  (compressible)"), ("theta_i", "th_i"),
                ("dw", "delta_w"), ("Cd", "C_d"), ("CJ", "C_J")]:
    xv = np.array([r[4][key] for r in rows])
    xn = xv / D if key not in ("Cd", "CJ") else xv
    m, c = np.polyfit(xn, y, 1)
    r2 = 1 - np.sum((y - (m * xn + c))**2) / np.sum((y - y.mean())**2)
    print("   %-22s slope=%+9.4f  intercept=%.4f  R2=%.4f" % (nm, m, c, r2))

# the D_eff null hypothesis:  x_fc proportional to D_eff
xv = np.array([r[4]["Deff"] / D for r in rows])
m, c = np.polyfit(xv, y, 1)
r2 = 1 - np.sum((y - (m * xv + c))**2) / np.sum((y - y.mean())**2)
k0 = np.mean(y / xv)
print("   %-22s slope=%+9.4f  intercept=%.4f  R2=%.4f" % ("D_eff/D", m, c, r2))
print("   through-origin  x_fc = %.4f (D_eff/D)   spread %.4f"
      % (k0, np.ptp(y / xv)))

# fit used to predict the sweep
xd = np.array([r[4]["dstar_c"] / D for r in rows])
mA, cA = np.polyfit(xd, y, 1)
print("\n   working law used below:  x_fc = %.4f %+.4f (d*_c/D_e)" % (cA, mA))

# ---------------------------------------------------------------- proposal
ARM_A = [0.0, 100.0, 170.0, 254.0, 340.0, 425.0]
anchor_a = 170.0
sa = scalars(anchor_a * 1e-6, "shifted")
ARM_B = [("d*-matched",  sa["dstar_c"], "dstar_c"),
         ("dw-matched",  sa["dw"],      "dw"),
         ("th-matched",  sa["theta_c"], "theta_c")]


def solve_wall(target, key):
    lo, hi = 1e-6, 4e-3
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if scalars(mid, "wall")[key] < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print("\nARM A - shifted tanh ladder (form fixed, a is the only variable)")
print(hdr + "   x_fc_pred")
for a in ARM_A:
    s = scalars(a * 1e-6, "shifted")
    xp = cA + mA * s["dstar_c"] / D
    print("%-22s %6.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.4f %7.4f %6.1f   %.3f"
          % ("A: shifted a=%g" % a, a, s["dw"] * 1e6, s["theta_i"] * 1e6,
             s["dstar_i"] * 1e6, s["theta_c"] * 1e6, s["dstar_c"] * 1e6,
             s["Cd"], s["Deff"] / D, s["dw"] / DX_LIP, xp))

print("\nARM B - wall tanh matched to the ARM A anchor a=%g um, one scalar at "
      "a time" % anchor_a)
print(hdr + "   x_fc_pred")
for nm, tgt, key in ARM_B:
    aw = solve_wall(tgt, key)
    s = scalars(aw, "wall")
    xp = cA + mA * s["dstar_c"] / D
    print("%-22s %6.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.4f %7.4f %6.1f   %.3f"
          % ("B: wall %s" % nm, aw * 1e6, s["dw"] * 1e6, s["theta_i"] * 1e6,
             s["dstar_i"] * 1e6, s["theta_c"] * 1e6, s["dstar_c"] * 1e6,
             s["Cd"], s["Deff"] / D, s["dw"] / DX_LIP, xp))

print("\nARM C - 2D twin of the remaining 3D variant")
s = scalars(200e-6, "wall")
print("%-22s %6.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.4f %7.4f %6.1f   %.3f"
      % ("C: wall a=200", 200.0, s["dw"] * 1e6, s["theta_i"] * 1e6,
         s["dstar_i"] * 1e6, s["theta_c"] * 1e6, s["dstar_c"] * 1e6,
         s["Cd"], s["Deff"] / D, s["dw"] / DX_LIP,
         cA + mA * s["dstar_c"] / D))
