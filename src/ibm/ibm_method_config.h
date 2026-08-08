#ifndef IBM_METHOD_CONFIG_H_
#define IBM_METHOD_CONFIG_H_

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <array>
#include <string>

namespace IBM::method_config {

inline constexpr const char* production_method =
    "pure_shared_gp_full_cartesian";

inline void validate_and_print_runtime_contract()
{
  amrex::ParmParse ib("ib");
  amrex::ParmParse cns("cns");

  constexpr std::array<const char*, 7> forbidden_ib_keys{
      "cut_control_database",
      "cut_control_runtime_audit",
      "cut_control_conservation_tolerance",
      "cut_control_budget_file",
      "cut_control_cp_reference_pressure",
      "cut_control_cp_reference_dynamic_pressure",
      "cut_control_wall_pressure_output"};
  constexpr std::array<const char*, 3> forbidden_cns_keys{
      "cut_control_low_order_admissibility_audit",
      "cut_control_forced_global_theta",
      "cut_control_forced_local_theta"};

  for (const char* key : forbidden_ib_keys) {
    if (ib.contains(key)) {
      amrex::Abort(
          std::string("pure GP-IBM forbids ib.") + key +
          "; rebuild with IBM_EXPERIMENTAL_CUT_CONTROL=TRUE only for an "
          "explicitly isolated hybrid experiment");
    }
  }
  for (const char* key : forbidden_cns_keys) {
    if (cns.contains(key)) {
      amrex::Abort(
          std::string("pure GP-IBM forbids cns.") + key +
          "; cut-control and aggregate limiting are outside the production "
          "method");
    }
  }

  amrex::Print()
      << "[IBM-Method] family=" << production_method
      << " pure_gp_production=1 cut_control_compiled=0"
         " control_volume=full_cartesian_cell"
         " face_area=full_cartesian_face"
         " solid_extension=shared_ghost_point\n";
}

}  // namespace IBM::method_config

#endif
