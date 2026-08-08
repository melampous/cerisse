#ifndef NOZZLEFUNCTIONS_H
#define NOZZLEFUNCTIONS_H

#include <cmath>

using namespace amrex;

// Isentropic-flow and nozzle geometry utilities for problem set-up.
// These functions are host-only helpers used to pre-compute initial/boundary
// conditions; they are not called inside GPU kernels.
// Note: the specific gas constant R = 287 J/(kg·K) assumes dry air (γ = 1.4).
//       Pass the appropriate Rgas for other working fluids.
namespace nozzle_functions {

  // Static temperature from stagnation temperature and Mach number (isentropic).
  static inline Real T_isen(Real T0, Real M, Real gamma)
  {
    return T0 / (1.0 + 0.5 * (gamma - 1.0) * M * M);
  }

  // Static pressure from stagnation pressure and Mach number (isentropic).
  static inline Real P_isen(Real P0, Real M, Real gamma)
  {
    const Real factor = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    return P0 * std::pow(factor, -gamma / (gamma - 1.0));
  }

  // Density from stagnation conditions and Mach number (ideal gas, R = 287 J/kg/K).
  static inline Real rho_isen(Real P0, Real T0, Real M, Real gamma,
                               Real Rgas = 287.0)
  {
    const Real T = T_isen(T0, M, gamma);
    const Real P = P_isen(P0, M, gamma);
    return P / (Rgas * T);
  }

  // Area ratio A/A* for isentropic flow at Mach M (choked throat reference).
  static inline Real nozzle_area_ratio(Real M, Real gamma)
  {
    const Real term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M * M);
    const Real exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    return std::pow(term, exponent) / M;
  }

  // Stagnation pressure from static pressure, Mach, and γ.
  static inline Real Pstag(Real P, Real M, Real gamma)
  {
    const Real term = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    return P * std::pow(term, gamma / (gamma - 1.0));
  }

  // Stagnation temperature from static temperature, Mach, and γ.
  static inline Real Tstag(Real T, Real M, Real gamma)
  {
    return T * (1.0 + 0.5 * (gamma - 1.0) * M * M);
  }

  // Choked (throat, M=1) pressure from stagnation pressure.
  static inline Real Pchok(Real P0, Real gamma)
  {
    return P0 * std::pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));
  }

  // Choked (throat, M=1) temperature from stagnation temperature.
  static inline Real Tchok(Real T0, Real gamma)
  {
    return T0 * (2.0 / (gamma + 1.0));
  }

  // Choked mass-flow rate [kg/s] per unit area (R = 287 J/kg/K for air).
  static inline Real masschok(Real T0, Real P0, Real gamma, Real Rgas = 287.0)
  {
    const Real gp1_o2 = 0.5 * (gamma + 1.0);
    return P0 * std::sqrt(gamma / (Rgas * T0))
               * std::pow(1.0 / gp1_o2, gp1_o2 / (gamma - 1.0));
  }

  // Static pressure at Mach M in a nozzle with stagnation pressure P0.
  static inline Real Pnozz(Real Mach, Real P0, Real gamma)
  {
    if (Mach < 0.0) return 0.0;
    const Real factor = 1.0 + 0.5 * (gamma - 1.0) * Mach * Mach;
    return P0 * std::pow(factor, -gamma / (gamma - 1.0));
  }

} // namespace nozzle_functions

#endif // NOZZLEFUNCTIONS_H
