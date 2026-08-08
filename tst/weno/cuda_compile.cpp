#include <Weno.h>
#include <Weno_old.h>

#include <type_traits>

namespace {

struct CompileClosure {
  static constexpr int NCONS = 5;
  static constexpr int NGHOST = 3;
  static constexpr int QRHO = 0;
  static constexpr int QPRES = 1;
  static constexpr int QC = 2;
  static constexpr int QU = 3;

  template <typename PrimitiveArray>
  AMREX_GPU_DEVICE amrex::Real max_char_speed(
      const amrex::IntVect& cell, int, int,
      const PrimitiveArray& primitive) const
  {
    return primitive(cell, QC);
  }

  template <typename PrimitiveArray>
  AMREX_GPU_DEVICE void prims2flux(
      const amrex::IntVect& cell, int,
      const PrimitiveArray& primitive, amrex::Real* output) const
  {
    for (int component = 0; component < NCONS; ++component) {
      output[component] = primitive(cell, component);
    }
  }

  template <typename PrimitiveArray>
  AMREX_GPU_DEVICE void prims2cons(
      const amrex::IntVect& cell, const PrimitiveArray& primitive,
      amrex::Real* output) const
  {
    for (int component = 0; component < NCONS; ++component) {
      output[component] = primitive(cell, component);
    }
  }

  template <typename PrimitiveArray>
  AMREX_GPU_DEVICE amrex::Real roe_avg_state(
      const amrex::IntVect&, int, const PrimitiveArray&) const
  {
    return amrex::Real(0.0);
  }

  AMREX_GPU_DEVICE void cons2char(
      amrex::Real, amrex::Real*) const {}

  AMREX_GPU_DEVICE void char2cons(
      amrex::Real, amrex::Real*) const {}
};

static_assert(
    !std::is_same_v<
        llf_wenoz5_t<CompileClosure>,
        llf_wenoz5_old_t<CompileClosure>>,
    "The formal and archived WENO drivers must remain distinct types");

}  // namespace

#ifdef AMREX_USE_CUDA
__global__ void weno_cuda_compile_kernel(amrex::Real* output)
{
  const amrex::Real stencil5[5] = {
      amrex::Real(0.0), amrex::Real(1.0), amrex::Real(2.0),
      amrex::Real(3.0), amrex::Real(4.0)};
  const amrex::Real stencil6[6] = {
      amrex::Real(0.0), amrex::Real(1.0), amrex::Real(2.0),
      amrex::Real(3.0), amrex::Real(4.0), amrex::Real(5.0)};
  output[0] =
      ReconScheme::WenoZ5::reconstruct_interface_value(stencil5);
  output[1] =
      ReconScheme::Teno5::reconstruct_interface_value(stencil5);
  output[2] =
      ReconScheme::Teno6::reconstruct_interface_value(stencil6);

  bool stencil_state_available[
      2 * ReconScheme::Teno6::required_ghost_cells] = {
      true, true, true, true, true, false};
  const auto reconstruction_plan =
      ReconScheme::build_face_reconstruction_plan<ReconScheme::Teno6>(
          stencil_state_available);
  output[3] =
      reconstruction_plan.positive_linear_weights[3] +
      reconstruction_plan.negative_linear_weights[2];
}
#else
void weno_cuda_compile_kernel(amrex::Real*) {}
#endif
