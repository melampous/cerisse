#include <Weno.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

using amrex::Real;

bool close(const Real actual, const Real expected)
{
  const Real scale =
      std::max(Real(1.0), std::max(std::abs(actual), std::abs(expected)));
  return std::abs(actual - expected) <=
         Real(128.0) * std::numeric_limits<Real>::epsilon() * scale;
}

bool require_close(
    const char* name, const Real actual, const Real expected)
{
  if (std::isfinite(actual) && close(actual, expected)) {
    return true;
  }
  std::cerr << name << " failed: actual=" << actual
            << ", expected=" << expected << '\n';
  return false;
}

bool require_teno5_selection(
    const char* name, const std::array<Real, 5>& stencil,
    const std::array<bool, 3>& expected)
{
  Real beta[ReconScheme::WenoZ5::candidate_count] = {};
  ReconScheme::WenoZ5::compute_smoothness_indicators(
      stencil.data(), beta);
  const Real tau = std::abs(
      std::abs(beta[2] - beta[0]) -
      (beta[2] + Real(4.0) * beta[1] + beta[0]) / Real(6.0));
  const bool available[3] = {true, true, true};
  bool selected[3] = {};
  const bool valid = ReconScheme::select_teno_smooth_candidates<3>(
      available, beta, tau, Real(1.0e-9), selected);
  const bool passed =
      valid && selected[0] == expected[0] &&
      selected[1] == expected[1] && selected[2] == expected[2];
  if (!passed) {
    std::cerr << name << " failed: selected={"
              << selected[0] << ',' << selected[1] << ',' << selected[2]
              << "}, expected={" << expected[0] << ',' << expected[1]
              << ',' << expected[2] << "}\n";
  }
  return passed;
}

template <typename ReconstructionMethod>
Real periodic_sine_derivative_error(const int cell_count)
{
  const Real pi = std::acos(Real(-1.0));
  const Real dx = Real(1.0) / Real(cell_count);
  std::vector<Real> face_flux(cell_count, Real(0.0));

  const auto sample = [&](int cell) {
    cell %= cell_count;
    if (cell < 0) cell += cell_count;
    const Real x = (Real(cell) + Real(0.5)) * dx;
    return std::sin(Real(2.0) * pi * x);
  };

  for (int face = 0; face < cell_count; ++face) {
    Real stencil[ReconstructionMethod::reconstruction_sample_count];
    if constexpr (
        ReconstructionMethod::reconstruction_sample_count == 5) {
      for (int point = 0; point < 5; ++point) {
        stencil[point] = sample(face + 1 - point);
      }
    } else {
      for (int point = 0; point < 6; ++point) {
        stencil[point] = sample(face + 2 - point);
      }
    }
    face_flux[face] =
        ReconstructionMethod::reconstruct_interface_value(stencil);
  }

  Real error_sum = Real(0.0);
  for (int cell = 0; cell < cell_count; ++cell) {
    const Real derivative =
        (face_flux[(cell + 1) % cell_count] - face_flux[cell]) / dx;
    const Real x = (Real(cell) + Real(0.5)) * dx;
    const Real exact = Real(2.0) * pi * std::cos(Real(2.0) * pi * x);
    const Real error = derivative - exact;
    error_sum += error * error;
  }
  return std::sqrt(error_sum / Real(cell_count));
}

template <typename ReconstructionMethod>
bool require_spatial_order(
    const char* name, const Real minimum_order)
{
  const Real coarse_error =
      periodic_sine_derivative_error<ReconstructionMethod>(64);
  const Real fine_error =
      periodic_sine_derivative_error<ReconstructionMethod>(128);
  const Real order = std::log(coarse_error / fine_error) / std::log(Real(2.0));
  std::cout << name << ": E64=" << coarse_error
            << ", E128=" << fine_error << ", order=" << order << '\n';
  if (std::isfinite(order) && order >= minimum_order) {
    return true;
  }
  std::cerr << name << " failed: E64=" << coarse_error
            << ", E128=" << fine_error << ", order=" << order
            << ", minimum=" << minimum_order << '\n';
  return false;
}

}  // namespace

int main()
{
  using ReconScheme::WenoZ5;
  using ReconScheme::Teno5;
  using ReconScheme::Teno6;

  bool passed = true;

  {
    const Real stencil[5] = {
        Real(2.0), Real(2.0), Real(2.0), Real(2.0), Real(2.0)};
    passed &= require_close(
        "WENO-Z5 constant preservation",
        WenoZ5::reconstruct_interface_value(stencil), Real(2.0));
  }

  {
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real stencil[5] = {
        nan, nan, Real(1.0), Real(1.0), Real(1.0)};
    const WenoZ5::CandidateWeightArray only_upwind = {
        Real(0.0), Real(0.0), Real(1.0)};
    passed &= require_close(
        "WENO-Z5 inactive-candidate NaN isolation",
        WenoZ5::reconstruct_interface_value(stencil, only_upwind),
        Real(1.0));
  }

  {
    const Real stencil[5] = {
        Real(2.0), Real(2.0), Real(2.0), Real(2.0), Real(2.0)};
    passed &= require_close(
        "TENO5 constant preservation",
        Teno5::reconstruct_interface_value(stencil), Real(2.0));
  }

  passed &= require_teno5_selection(
      "TENO5 smooth-stencil retention at CT=1e-9",
      {Real(0.0), Real(1.0), Real(2.0), Real(3.0), Real(4.0)},
      {true, true, true});
  passed &= require_teno5_selection(
      "TENO5 rightmost jump rejection at CT=1e-9",
      {Real(0.0), Real(0.0), Real(0.0), Real(0.0), Real(1.0)},
      {true, true, false});
  passed &= require_teno5_selection(
      "TENO5 right jump rejection at CT=1e-9",
      {Real(0.0), Real(0.0), Real(0.0), Real(1.0), Real(1.0)},
      {true, false, false});
  passed &= require_teno5_selection(
      "TENO5 left jump rejection at CT=1e-9",
      {Real(0.0), Real(0.0), Real(1.0), Real(1.0), Real(1.0)},
      {false, false, true});
  passed &= require_teno5_selection(
      "TENO5 leftmost jump rejection at CT=1e-9",
      {Real(0.0), Real(1.0), Real(1.0), Real(1.0), Real(1.0)},
      {false, true, true});

  {
    const Real stencil[6] = {
        Real(0.0), Real(1.0), Real(2.0),
        Real(3.0), Real(4.0), Real(5.0)};
    passed &= require_close(
        "TENO6 smooth linear reconstruction",
        Teno6::reconstruct_interface_value(stencil), Real(2.5));
  }

  {
    const Real stencil[6] = {
        Real(0.0), Real(0.0), Real(0.0),
        Real(0.0), Real(0.0), Real(1.0)};
    passed &= require_close(
        "TENO6 sixth-point discontinuity detection",
        Teno6::reconstruct_interface_value(stencil), Real(0.0));
  }

  {
    const Real stencil[6] = {
        Real(0.0), Real(0.0), Real(0.0),
        Real(1.0), Real(1.0), Real(1.0)};
    passed &= require_close(
        "TENO6 interface discontinuity selection",
        Teno6::reconstruct_interface_value(stencil), Real(1.0));
  }

  {
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real stencil[6] = {
        nan, nan, nan, Real(1.0), Real(1.0), Real(1.0)};
    const Teno6::CandidateWeightArray only_upwind = {
        Real(0.0), Real(0.0), Real(1.0), Real(0.0)};
    passed &= require_close(
        "TENO6 inactive-candidate NaN isolation",
        Teno6::reconstruct_interface_value(stencil, only_upwind),
        Real(1.0));
  }

  {
    bool stencil_state_available[2 * Teno6::required_ghost_cells] = {
        true, true, true, true, true, false};
    const auto reconstruction_plan =
        ReconScheme::build_face_reconstruction_plan<Teno6>(
            stencil_state_available);
    const bool reconstruction_plan_passed =
        close(
            reconstruction_plan.positive_linear_weights[3],
            Real(0.0)) &&
        close(
            reconstruction_plan.negative_linear_weights[2],
            Real(0.0)) &&
        reconstruction_plan.high_order_reconstruction_available &&
        !reconstruction_plan.full_stencil_available;
    passed &= reconstruction_plan_passed;
    if (!reconstruction_plan_passed) {
      std::cerr << "TENO6 ghost-point reconstruction plan failed\n";
    }
  }

  passed &= require_spatial_order<WenoZ5>(
      "WENO-Z5 periodic smooth derivative order", Real(4.5));
  passed &= require_spatial_order<Teno5>(
      "TENO5 periodic smooth derivative order", Real(4.5));
  passed &= require_spatial_order<Teno6>(
      "TENO6 periodic smooth derivative order", Real(5.5));

  if (!passed) {
    return 1;
  }
  std::cout << "WENO reconstruction tests passed\n";
  return 0;
}
