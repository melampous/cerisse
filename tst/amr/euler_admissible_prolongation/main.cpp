#include <AMReX.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>

#include <EulerAdmissibleProlongation.H>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

using amrex::BCRec;
using amrex::Box;
using amrex::FArrayBox;
using amrex::Geometry;
using amrex::IntVect;
using amrex::Real;

constexpr int mx_comp = 0;
constexpr int my_comp = 1;
constexpr int mz_comp = 2;
constexpr int energy_comp = 3;
constexpr int rho_comp = 4;
constexpr int ncomp = 5;
constexpr Real rhoe_floor = Real(2.5e-8);

using State = std::array<Real, ncomp>;

void require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

Real internal_energy_density(const FArrayBox& state, const IntVect& iv) {
  const Real rho = state(iv, rho_comp);
  const Real mx = state(iv, mx_comp);
  const Real my = state(iv, my_comp);
  const Real mz = state(iv, mz_comp);
  return state(iv, energy_comp) -
         Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
}

Real specific_internal_energy(const FArrayBox& state, const IntVect& iv) {
  return internal_energy_density(state, iv) / state(iv, rho_comp);
}

void set_state(FArrayBox& fab, const IntVect& iv, const State& state) {
  for (int n = 0; n < ncomp; ++n) fab(iv, n) = state[n];
}

void interpolate(amrex::Interpolater& mapper, const FArrayBox& coarse,
                 FArrayBox& fine, const Box& fine_box,
                 const IntVect& ratio, const Geometry& coarse_geom,
                 const Geometry& fine_geom,
                 const amrex::Vector<BCRec>& bcs) {
#ifdef AMREX_USE_GPU
  constexpr amrex::RunOn runon = amrex::RunOn::Gpu;
#else
  constexpr amrex::RunOn runon = amrex::RunOn::Cpu;
#endif
  mapper.interp(coarse, 0, fine, 0, ncomp, fine_box, ratio, coarse_geom,
                fine_geom, bcs, 0, 0, runon);
  amrex::Gpu::streamSynchronize();
}

}  // namespace

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

#if (AMREX_SPACEDIM == 3)
  const IntVect coarse_lo(AMREX_D_DECL(-1, -1, -1));
  const IntVect coarse_hi(AMREX_D_DECL(1, 1, 1));
  const Box coarse_box(coarse_lo, coarse_hi);
  const IntVect ratio(AMREX_D_DECL(2, 2, 2));
  const Box fine_domain = amrex::refine(coarse_box, ratio);
  const amrex::RealBox physical_box(
      {AMREX_D_DECL(Real(-1.5), Real(-1.5), Real(-1.5))},
      {AMREX_D_DECL(Real(1.5), Real(1.5), Real(1.5))});
  const amrex::Vector<int> periodic(AMREX_SPACEDIM, 0);
  const Geometry coarse_geom(coarse_box, &physical_box, 0, periodic.data());
  const Geometry fine_geom(fine_domain, &physical_box, 0, periodic.data());
  const Box children(IntVect(AMREX_D_DECL(0, 0, 0)),
                     IntVect(AMREX_D_DECL(1, 1, 1)));

  amrex::Vector<BCRec> bcs(ncomp);
  for (auto& bc : bcs) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      bc.setLo(dir, amrex::BCType::int_dir);
      bc.setHi(dir, amrex::BCType::int_dir);
    }
  }

  const State parent = {33.495503118828978, -3.2320043225558757,
                        5.2191074520630023, 15235.465250857507,
                        0.053085159232775986};
  FArrayBox coarse(coarse_box, ncomp);
  for (amrex::BoxIterator iterator(coarse_box); iterator.ok(); ++iterator) {
    set_state(coarse, iterator(), parent);
  }
  set_state(coarse, IntVect(AMREX_D_DECL(-1, 0, 0)),
            {22.96435755216833, -0.28541561268893728, 0.40556702588223104,
             9959.7923848898117, 0.03170434393441865});
  set_state(coarse, IntVect(AMREX_D_DECL(1, 0, 0)),
            {44.795279616204432, -11.819455443160917, 19.125922207630193,
             26568.713646089782, 0.097527844278812295});
  set_state(coarse, IntVect(AMREX_D_DECL(0, -1, 0)),
            {29.359905981334975, -1.3092597970045508, 2.280906400315001,
             12484.796024345289, 0.043586432262205134});
  set_state(coarse, IntVect(AMREX_D_DECL(0, 1, 0)),
            {40.258108819133881, -8.5511528228091667, 14.013797206867203,
             22177.707042727296, 0.077222937608389414});
  set_state(coarse, IntVect(AMREX_D_DECL(0, 0, -1)),
            {43.524627604490902, -12.00421085165523, 18.697681135562718,
             26017.072367434495, 0.092633648566989357});
  set_state(coarse, IntVect(AMREX_D_DECL(0, 0, 1)),
            {24.036510472541202, -0.28297704582253869, 0.43833798581392946,
             10314.604821345292, 0.033008981656025511});

  FArrayBox legacy(children, ncomp);
  FArrayBox protected_state(children, ncomp);
  interpolate(amrex::lincc_interp, coarse, legacy, children, ratio,
              coarse_geom, fine_geom, bcs);
  cerisse::amr::EulerAdmissibleConservativeLinear protected_mapper(
      {mx_comp, my_comp, mz_comp, energy_comp, rho_comp, Real(1.0e-19),
       rhoe_floor});
  interpolate(protected_mapper, coarse, protected_state, children, ratio,
              coarse_geom, fine_geom, bcs);

  const IntVect logged_bad_child(AMREX_D_DECL(0, 0, 1));
  const Real legacy_bad_rhoe =
      internal_energy_density(legacy, logged_bad_child);
  require(std::abs(legacy_bad_rhoe - Real(-984.86910627266479)) < Real(1e-8),
          "legacy lincc_interp did not reproduce the frozen bad child");

  Real minimum_protected_rhoe = std::numeric_limits<Real>::max();
  State protected_average{};
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    const IntVect& iv = iterator();
    minimum_protected_rhoe = amrex::min(
        minimum_protected_rhoe, internal_energy_density(protected_state, iv));
    for (int n = 0; n < ncomp; ++n) {
      protected_average[n] += protected_state(iv, n) / Real(8.0);
    }
  }
  require(minimum_protected_rhoe > rhoe_floor,
          "protected child block is not Euler-admissible");
  Real conservation_residual = Real(0.0);
  for (int n = 0; n < ncomp; ++n) {
    conservation_residual = amrex::max(
        conservation_residual, std::abs(protected_average[n] - parent[n]));
  }
  require(conservation_residual < Real(5e-12),
          "protected child average does not recover the parent state");

  const Real observed_theta =
      (protected_state(logged_bad_child, mx_comp) - parent[mx_comp]) /
      (legacy(logged_bad_child, mx_comp) - parent[mx_comp]);
  const Real expected_theta =
      Real(0.8465959907344502) * (Real(1.0) - Real(1.0e-10));
  require(std::abs(observed_theta - expected_theta) < Real(2e-12),
          "parent-block theta differs from the frozen offline oracle");

  constexpr Real specific_ei_floor = Real(7000.0);
  FArrayBox thermodynamic_state(children, ncomp);
  cerisse::amr::EulerAdmissibleConservativeLinear thermodynamic_mapper(
      {mx_comp, my_comp, mz_comp, energy_comp, rho_comp, Real(1.0e-19),
       rhoe_floor, specific_ei_floor});
  interpolate(thermodynamic_mapper, coarse, thermodynamic_state, children,
              ratio, coarse_geom, fine_geom, bcs);
  Real minimum_thermodynamic_ei = std::numeric_limits<Real>::max();
  State thermodynamic_average{};
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    const IntVect& iv = iterator();
    minimum_thermodynamic_ei = amrex::min(
        minimum_thermodynamic_ei,
        specific_internal_energy(thermodynamic_state, iv));
    for (int n = 0; n < ncomp; ++n) {
      thermodynamic_average[n] += thermodynamic_state(iv, n) / Real(8.0);
    }
  }
  require(minimum_thermodynamic_ei > specific_ei_floor,
          "thermodynamic-floor child block violates its specific-e floor");
  Real thermodynamic_conservation_residual = Real(0.0);
  for (int n = 0; n < ncomp; ++n) {
    thermodynamic_conservation_residual = amrex::max(
        thermodynamic_conservation_residual,
        std::abs(thermodynamic_average[n] - parent[n]));
  }
  require(thermodynamic_conservation_residual < Real(5e-12),
          "thermodynamic-floor child average does not recover the parent");
  const Real thermodynamic_theta =
      (thermodynamic_state(logged_bad_child, mx_comp) - parent[mx_comp]) /
      (legacy(logged_bad_child, mx_comp) - parent[mx_comp]);
  require(thermodynamic_theta > Real(0.0) &&
              thermodynamic_theta < observed_theta,
          "specific-e floor did not strengthen the active parent theta");

  const State uniform = {2.0, -0.5, 0.25, 20000.0, 1.0};
  for (amrex::BoxIterator iterator(coarse_box); iterator.ok(); ++iterator) {
    set_state(coarse, iterator(), uniform);
  }
  interpolate(amrex::lincc_interp, coarse, legacy, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(protected_mapper, coarse, protected_state, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(thermodynamic_mapper, coarse, thermodynamic_state, children,
              ratio, coarse_geom, fine_geom, bcs);
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    for (int n = 0; n < ncomp; ++n) {
      require(protected_state(iterator(), n) == legacy(iterator(), n),
              "admissible uniform interpolation is not exactly transparent");
      require(thermodynamic_state(iterator(), n) == legacy(iterator(), n),
              "thermodynamic-floor uniform interpolation is not exactly transparent");
    }
  }

  for (amrex::BoxIterator iterator(coarse_box); iterator.ok(); ++iterator) {
    const IntVect& iv = iterator();
    set_state(coarse, iv,
              {Real(0.2) + Real(0.03) * iv[0],
               Real(-0.1) + Real(0.02) * iv[1],
               Real(0.05) - Real(0.01) * iv[2],
               Real(20000.0) + Real(0.4) * iv[0] - Real(0.2) * iv[1] +
                   Real(0.1) * iv[2],
               Real(1.0) + Real(0.02) * iv[0] + Real(0.01) * iv[1]});
  }
  interpolate(amrex::lincc_interp, coarse, legacy, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(protected_mapper, coarse, protected_state, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(thermodynamic_mapper, coarse, thermodynamic_state, children,
              ratio, coarse_geom, fine_geom, bcs);
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    for (int n = 0; n < ncomp; ++n) {
      require(protected_state(iterator(), n) == legacy(iterator(), n),
              "admissible smooth 3-D interpolation is not bitwise transparent");
      require(thermodynamic_state(iterator(), n) == legacy(iterator(), n),
              "thermodynamic-floor smooth 3-D interpolation is not bitwise transparent");
    }
  }

  std::cout << std::setprecision(17)
            << "PASS cartesian_3d legacy_bad_rhoe=" << legacy_bad_rhoe
            << " protected_min_rhoe=" << minimum_protected_rhoe
            << " theta=" << observed_theta
            << " conservation_linf=" << conservation_residual
            << " thermo_min_ei=" << minimum_thermodynamic_ei
            << " thermo_theta=" << thermodynamic_theta
            << " thermo_conservation_linf="
            << thermodynamic_conservation_residual << '\n';
#elif (AMREX_SPACEDIM == 2)
  const Box coarse_box(IntVect(AMREX_D_DECL(0, 0, 0)),
                       IntVect(AMREX_D_DECL(2, 2, 0)));
  const IntVect ratio(AMREX_D_DECL(2, 2, 1));
  const Box fine_domain = amrex::refine(coarse_box, ratio);
  const amrex::RealBox physical_box(
      {AMREX_D_DECL(Real(0.0), Real(-1.5), Real(0.0))},
      {AMREX_D_DECL(Real(3.0), Real(1.5), Real(1.0))});
  const amrex::Vector<int> periodic(AMREX_SPACEDIM, 0);
  const Geometry coarse_geom(coarse_box, &physical_box, 1, periodic.data());
  const Geometry fine_geom(fine_domain, &physical_box, 1, periodic.data());
  const Box children(IntVect(AMREX_D_DECL(2, 2, 0)),
                     IntVect(AMREX_D_DECL(3, 3, 0)));

  amrex::Vector<BCRec> bcs(ncomp);
  for (auto& bc : bcs) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      bc.setLo(dir, amrex::BCType::int_dir);
      bc.setHi(dir, amrex::BCType::int_dir);
    }
  }

  const State parent = {0.0, 0.0, 0.0, 2.0, 1.0};
  FArrayBox coarse(coarse_box, ncomp);
  for (amrex::BoxIterator iterator(coarse_box); iterator.ok(); ++iterator) {
    set_state(coarse, iterator(), parent);
  }
  set_state(coarse, IntVect(AMREX_D_DECL(0, 1, 0)),
            {-8.0, 0.0, 0.0, 40.0, 1.0});
  set_state(coarse, IntVect(AMREX_D_DECL(2, 1, 0)),
            {8.0, 0.0, 0.0, 40.0, 1.0});
  set_state(coarse, IntVect(AMREX_D_DECL(1, 0, 0)),
            {0.0, -8.0, 0.0, 40.0, 1.0});
  set_state(coarse, IntVect(AMREX_D_DECL(1, 2, 0)),
            {0.0, 8.0, 0.0, 40.0, 1.0});

  FArrayBox legacy(children, ncomp);
  FArrayBox protected_state(children, ncomp);
  FArrayBox thermodynamic_state(children, ncomp);
  interpolate(amrex::lincc_interp, coarse, legacy, children, ratio,
              coarse_geom, fine_geom, bcs);
  cerisse::amr::EulerAdmissibleConservativeLinear protected_mapper(
      {mx_comp, my_comp, mz_comp, energy_comp, rho_comp, Real(1.0e-19),
       rhoe_floor});
  interpolate(protected_mapper, coarse, protected_state, children, ratio,
              coarse_geom, fine_geom, bcs);

  Real minimum_legacy_rhoe = std::numeric_limits<Real>::max();
  Real minimum_protected_rhoe = std::numeric_limits<Real>::max();
  State protected_weighted_average{};
  Real total_weight = Real(0.0);
  const Real dr = fine_geom.CellSize(0);
  const Real radial_lo = fine_geom.ProbLo(0);
  Real common_theta = Real(-1.0);
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    const IntVect& iv = iterator();
    minimum_legacy_rhoe = amrex::min(
        minimum_legacy_rhoe, internal_energy_density(legacy, iv));
    minimum_protected_rhoe = amrex::min(
        minimum_protected_rhoe, internal_energy_density(protected_state, iv));

    const Real r_inner = radial_lo + Real(iv[0]) * dr;
    const Real r_outer = r_inner + dr;
    const Real annular_weight =
        r_outer * r_outer - r_inner * r_inner;
    total_weight += annular_weight;
    for (int n = 0; n < ncomp; ++n) {
      protected_weighted_average[n] +=
          annular_weight * protected_state(iv, n);
      const Real deviation = legacy(iv, n) - parent[n];
      if (std::abs(deviation) > Real(1.0e-12)) {
        const Real component_theta =
            (protected_state(iv, n) - parent[n]) / deviation;
        if (common_theta < Real(0.0)) common_theta = component_theta;
        require(std::abs(component_theta - common_theta) < Real(2.0e-12),
                "R-Z child components do not share one parent-block theta");
      }
    }
  }
  require(minimum_legacy_rhoe < Real(0.0),
          "R-Z legacy interpolation did not create the intended bad child");
  require(minimum_protected_rhoe > rhoe_floor,
          "R-Z protected child block is not Euler-admissible");
  require(common_theta > Real(0.0) && common_theta < Real(1.0),
          "R-Z regression did not activate a nontrivial theta limiter");

  Real conservation_residual = Real(0.0);
  for (int n = 0; n < ncomp; ++n) {
    protected_weighted_average[n] /= total_weight;
    conservation_residual = amrex::max(
        conservation_residual,
        std::abs(protected_weighted_average[n] - parent[n]));
  }
  require(conservation_residual < Real(5.0e-12),
          "R-Z protected annular average does not recover the parent state");

  constexpr Real specific_ei_floor = Real(1.0);
  cerisse::amr::EulerAdmissibleConservativeLinear thermodynamic_mapper(
      {mx_comp, my_comp, mz_comp, energy_comp, rho_comp, Real(1.0e-19),
       rhoe_floor, specific_ei_floor});
  interpolate(thermodynamic_mapper, coarse, thermodynamic_state, children,
              ratio, coarse_geom, fine_geom, bcs);
  Real minimum_thermodynamic_ei = std::numeric_limits<Real>::max();
  State thermodynamic_weighted_average{};
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    const IntVect& iv = iterator();
    minimum_thermodynamic_ei = amrex::min(
        minimum_thermodynamic_ei,
        specific_internal_energy(thermodynamic_state, iv));
    const Real r_inner = radial_lo + Real(iv[0]) * dr;
    const Real r_outer = r_inner + dr;
    const Real annular_weight =
        r_outer * r_outer - r_inner * r_inner;
    for (int n = 0; n < ncomp; ++n) {
      thermodynamic_weighted_average[n] +=
          annular_weight * thermodynamic_state(iv, n);
    }
  }
  require(minimum_thermodynamic_ei > specific_ei_floor,
          "R-Z thermodynamic-floor child block violates its specific-e floor");
  Real thermodynamic_conservation_residual = Real(0.0);
  for (int n = 0; n < ncomp; ++n) {
    thermodynamic_weighted_average[n] /= total_weight;
    thermodynamic_conservation_residual = amrex::max(
        thermodynamic_conservation_residual,
        std::abs(thermodynamic_weighted_average[n] - parent[n]));
  }
  require(thermodynamic_conservation_residual < Real(5.0e-12),
          "R-Z thermodynamic-floor annular average does not recover parent");

  const State uniform = {2.0, -0.5, 0.25, 20.0, 1.0};
  for (amrex::BoxIterator iterator(coarse_box); iterator.ok(); ++iterator) {
    set_state(coarse, iterator(), uniform);
  }
  interpolate(amrex::lincc_interp, coarse, legacy, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(protected_mapper, coarse, protected_state, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(thermodynamic_mapper, coarse, thermodynamic_state, children,
              ratio, coarse_geom, fine_geom, bcs);
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    for (int n = 0; n < ncomp; ++n) {
      require(protected_state(iterator(), n) == legacy(iterator(), n),
              "R-Z uniform interpolation is not exactly transparent");
      require(thermodynamic_state(iterator(), n) == legacy(iterator(), n),
              "R-Z thermodynamic-floor uniform interpolation is not exactly transparent");
    }
  }

  for (amrex::BoxIterator iterator(coarse_box); iterator.ok(); ++iterator) {
    const IntVect& iv = iterator();
    set_state(coarse, iv,
              {Real(0.2) + Real(0.03) * iv[0],
               Real(-0.1) + Real(0.02) * iv[1], Real(0.05),
               Real(50.0) + Real(0.4) * iv[0] - Real(0.2) * iv[1],
               Real(1.0) + Real(0.02) * iv[0] + Real(0.01) * iv[1]});
  }
  interpolate(amrex::lincc_interp, coarse, legacy, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(protected_mapper, coarse, protected_state, children, ratio,
              coarse_geom, fine_geom, bcs);
  interpolate(thermodynamic_mapper, coarse, thermodynamic_state, children,
              ratio, coarse_geom, fine_geom, bcs);
  for (amrex::BoxIterator iterator(children); iterator.ok(); ++iterator) {
    for (int n = 0; n < ncomp; ++n) {
      require(protected_state(iterator(), n) == legacy(iterator(), n),
              "admissible smooth R-Z interpolation is not bitwise transparent");
      require(thermodynamic_state(iterator(), n) == legacy(iterator(), n),
              "thermodynamic-floor smooth R-Z interpolation is not bitwise transparent");
    }
  }

  std::cout << std::setprecision(17)
            << "PASS rz_2d legacy_min_rhoe=" << minimum_legacy_rhoe
            << " protected_min_rhoe=" << minimum_protected_rhoe
            << " theta=" << common_theta
            << " annular_conservation_linf=" << conservation_residual
            << " thermo_min_ei=" << minimum_thermodynamic_ei
            << " thermo_annular_conservation_linf="
            << thermodynamic_conservation_residual << '\n';
#else
#error "Euler-admissible prolongation regression supports DIM=2 or DIM=3"
#endif
  amrex::Finalize();
  return EXIT_SUCCESS;
}
