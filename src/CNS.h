#ifndef CNS_H_
#define CNS_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Math.H>
#include <prob.h>
#include <CNSconstants.h>
#include <nscbc.h>
#include <RestartDtControl.h>

#include <Utilities.h>

#include <array>

using namespace amrex;

class CNS : public amrex::AmrLevel {
 public:
  // Init --------------------------------------------------------------------
  CNS();
  CNS(amrex::Amr& papa, int lev, const amrex::Geometry& level_geom,
      const amrex::BoxArray& bl, const amrex::DistributionMapping& dm,
      amrex::Real time);
  ~CNS();

  CNS(const CNS& rhs) = delete;
  CNS& operator=(const CNS& rhs) = delete;

  // Read parameters
  static void read_params();

  // Define data descriptors.
  static void variableSetUp();

  // Cleanup data descriptors at end of run.
  static void variableCleanUp();

  // Initialize data on this level from another CNS (during regrid).
  void init(amrex::AmrLevel& old) override;

  // Initialize data on this level after regridding if old level did not
  // previously exist
  void init() override;

  // Initialize grid data at problem start-up.
  virtual void initData() override;

  // Do work after init().
  virtual void post_init(amrex::Real stop_time) override;
  // -------------------------------------------------------------------------

  // Time-stepping -----------------------------------------------------------
  void compute_rhs(amrex::MultiFab& S, amrex::Real dt,
                   amrex::FluxRegister* fr_as_crse,
                   amrex::FluxRegister* fr_as_fine,
                   amrex::Real stage_time,
                   amrex::Real reflux_dt,
                   std::array<amrex::MultiFab*, AMREX_SPACEDIM>
                       captured_face_flux = {},
                   amrex::MultiFab* captured_rz_pressure_face_flux = nullptr,
                   std::array<amrex::iMultiFab*, AMREX_SPACEDIM>
                       captured_llf_fallback_mask = {});

#if NUM_SPECIES > 1
  void clip_species_state(amrex::MultiFab& S);                   
#endif  

  // void computeTemp(amrex::MultiFab& State, int ng);

  GpuArray<Real,AMREX_SPACEDIM> maxEigen();

  // Compute initial time step.
  amrex::Real initialTimeStep();

  void computeInitialDt(int finest_level, int sub_cycle,
                        amrex::Vector<int>& n_cycle,
                        const amrex::Vector<amrex::IntVect>& ref_ratio,
                        amrex::Vector<amrex::Real>& dt_level,
                        amrex::Real stop_time) override;

  void computeNewDt(int finest_level, int sub_cycle,
                    amrex::Vector<int>& n_cycle,
                    const amrex::Vector<amrex::IntVect>& ref_ratio,
                    amrex::Vector<amrex::Real>& dt_min,
                    amrex::Vector<amrex::Real>& dt_level, amrex::Real stop_time,
                    int post_regrid_flag) override;

  // Advance grids at this level in time.
  Real advance(amrex::Real time, amrex::Real dt, int iteration,
               int ncycle) override;

  // Do work after timestep().
  virtual void post_timestep(int iteration) override;

  virtual void postCoarseTimeStep(Real time) override;

  virtual void post_restart() override;

  void checkPointPost(const std::string& dir, std::ostream& os) override;

  void set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint) override;

  // -------------------------------------------------------------------------

  // Gridding ----------------------------------------------------------------
  virtual void post_regrid(int lbase, int new_finest) override;

#ifdef AMREX_USE_GPIBM
  void rebuildIBM(bool rebuild_surface = true);

  // Default-off, read-only context retained across the first advance after a
  // regrid.  It identifies fine valid cells initialized from the coarse level
  // and preserves the exact coarse source stencil used at the regrid time.
  void captureRegridProlongationContext(const CNS* old_level,
                                        amrex::Real time);
#endif

  // Error estimation for regridding.
  // virtual void errorEst (int lev, TagBoxArray& tags, Real time, int ngrow);
  virtual void errorEst(amrex::TagBoxArray& tb, int clearval, int tagval,
                        amrex::Real time, int n_error_buf = 0,
                        int ngrow = 0) override;

  // init
  CNS& getLevel(int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

  enum StateDataType { State_Type = 0, Stats_Type, Cost_Type };

  void buildMetrics();

  // Persistent ghost-cell NSCBC state.  The ghost state is integrated with
  // the same RK scheme as the interior conservative state.
  void initialize_nscbc_ghost_state(amrex::Real time);
  void copy_nscbc_ghost_to_state(amrex::MultiFab& state,
                                 amrex::MultiFab const& ghost) const;
  void compute_nscbc_ghost_rhs(amrex::MultiFab& state,
                               amrex::MultiFab& ghost_rhs) const;
  void fill_nscbc_gc_ghosts(amrex::MultiFab& state) const;

  static AMREX_FORCE_INLINE void rz_sanity_check(amrex::Geometry const& geom)
  {
    // RZ axis sanity check: for RZ the axis must be at r=0 and r-direction
    // cannot be periodic.
    if (geom.IsRZ()) {
#if (AMREX_SPACEDIM != 2)
      amrex::Abort("RZ requires AMREX_SPACEDIM=2 (axisymmetric r-z)");
#endif
      if (geom.isPeriodic(0)) {
        amrex::Abort(
            "RZ requires geometry.is_periodic[0]=0 (non-periodic r-direction)");
      }

      const amrex::Real rlo = geom.ProbLo(0);
      // Zero is exactly representable.  Requiring it here keeps every RZ
      // axis branch consistent and prevents division by a roundoff-sized
      // positive lower radius in face-flux assembly.
      if (rlo != amrex::Real(0.0)) {
        amrex::Abort("RZ requires geometry.prob_lo[0]=0 (axis at r=0)");
      }
    }
  }

  void avgDown();

#ifdef AMREX_USE_GPIBM
  // Re-publish the unique shared ghost-point state after AMR synchronization.
  // Reflux/average-down operate on the accepted conservative MultiFab after
  // the end-of-stage GP reconstruction, so solid-side GP cells must be
  // reconstructed again before they can serve as a coarse FillPatch source.
  // This routine writes GP cells only; active fluid cells are untouched.
  void refreshAcceptedSharedGPState(amrex::Real time);
#endif

  void printTotal() const;

  virtual void writePlotFile(const std::string& dir, std::ostream& os,
                             VisMF::How how = VisMF::NFiles) override;

  virtual void writePlotFilePost(const std::string& dir,
                                 std::ostream& os) override;

#if AMREX_USE_GPIBM
  virtual void writeSurfFile(bool force_compute = false,
                             bool force_output = false);
#endif

  // diagnostics
  static bool record_probe;
  void setupTimeProbe();
  void recordTimeProbe();
  static int time_probe_lev;
  static int time_probe_int;
  static amrex::Vector<std::string> time_probe_names;
  static amrex::Vector<amrex::Box> time_probe_boxes;


  // Parameters
  static int num_state_data_types;
  std::unique_ptr<amrex::FluxRegister> flux_reg;
  static int do_reflux;
  // Spatial interpolation used for state data at AMR coarse-fine interfaces.
  // 0: AMReX limited conservative linear prolongation.
  // 1: AMReX conservative quartic prolongation for smooth-flow verification.
  // This transfer choice does not change the point-sample finite-difference
  // semantics of the evolved CERISSE state. The default remains linear.
  static int amr_state_interp;
  static bool regrid_prolongation_diagnostics;

  static bool verbose;
  // static amrex::IntVect hydro_tile_size;
  static amrex::Real cfl;

  static int refine_max_dengrad_lev;
  static amrex::Real refine_dengrad;

  // Statistics
  static amrex::Real time_stats;
  static amrex::Real time_stat_level[10];
  static amrex::Real stats_start_time;
  static bool compute_stats, record_stats;
  static int INDEX_THERM;
  void setupStats();
  void computeStats();
  void debugStats(int is, int js, int ks);

  static amrex::Real gravity;

  static amrex::Real dt_constant;
  static bool dt_dynamic;
  static amrex::Real dt_max;        // optional absolute cap on dt (cns.dt_max)
  // Default-off exception for segmented checkpoint continuations.  Only the
  // first normal computeNewDt after a restart may bypass the legacy 1.1x
  // old-dt growth cap; CFL, dt_max, stop_time, and post-regrid caps remain.
  static bool restart_first_dt_from_cfl;
  cerisse::time_step::RestartFirstDtFromCflControl
      restart_first_dt_from_cfl_control;
  static int nstep_screen_output;
  static int dist_linear;
  static int order_rk;
  static int stages_rk;

  // cns.nscbc_{lo,hi}: 0=off, 1=relaxed inflow, 2=pure outflow,
  // 3=pressure-relaxed outflow.  The default implementation time-integrates
  // persistent ghost states; the opt-in GC implementation reconstructs them
  // from the current RK-stage interior state.
  static bool use_nscbc;
  static amrex::GpuArray<int, AMREX_SPACEDIM> nscbc_lo;
  static amrex::GpuArray<int, AMREX_SPACEDIM> nscbc_hi;
  static nscbc::Parm nscbc_parm;

  std::unique_ptr<amrex::MultiFab> nscbc_ghost_state;
  bool nscbc_ghost_initialized = false;

  // When true, the end-of-step IBM abort check uses "approaching the
  // smallr / ei_min clipping floors" as the failure criterion instead of
  // "already non-positive". Catches silent clipping in cons2prims that
  // would otherwise mask a numerical breakdown.
  static bool strict_positivity;

  // Opt-in conservative shared-face positivity limiter for the frozen
  // ideal-gas Euler-slip GP path. Cartesian 2-D/3-D stages blend one
  // face-flux family; paired R-Z stages blend the complete flux and pressure
  // companion with the same theta. Every SSPRK(4,3) Forward-Euler bracket is
  // checked independently. The default is false.
  static bool ibm_positivity_flux_limiter;
  static bool ibm_positivity_flux_limiter_verbose;
  // Reserved fail-closed switch for a future pure full-cell SSPRK retry path.
  // The current certified limiter performs local, global-theta, and all-low
  // fallback within each Forward-Euler bracket; whole-step retry is disabled.
  static bool ibm_positivity_retry;
  static int ibm_positivity_max_step_halvings;

#ifdef AMREX_USE_GPIBM
  // One means that this valid fine cell had no old-fine overlap and was
  // initialized through coarse-to-fine prolongation.  The coarse state uses
  // the coarsened fine BoxArray and the same DistributionMapping so every
  // bad child and its parent source are local to the same MPI rank.
  std::unique_ptr<amrex::iMultiFab> regrid_new_from_coarse_mask;
  std::unique_ptr<amrex::MultiFab> regrid_coarse_source_state;
#endif


  // When true, the end-of-step IBM check clips bad fluid cells (rho<=0
  // or E<=0) to (rho_floor, rho*ei_floor + KE) and continues, instead
  // of aborting. NaN/Inf cells still abort. Intended as an opt-in
  // robustness net for cases where WENO-Z5 briefly overshoots in
  // strong expansions at IBM sharp corners. Not conservative locally;
  // the number of clipped cells is reported every step it fires.
  static bool soft_positivity;

  // Pass 2 (flood-fill interior solid + zero momentum) always runs for FSI.
  // For static geometry, it is opt-in: some complex-geometry / multi-body
  // cases are actually *destabilised* by the flood-fill because the averaged
  // neighbour values spread post-shock states into the body, which the next
  // step's WENO stencil then reads back and oscillates on. Enable via
  // cns.pass2_static = 1 only when you explicitly want that behaviour.
  static bool pass2_static;

  // Utility-variables
  static bool use_utility;
  static Utility utilidades;

  // IBM-specific keywords
  static bool ib_move;
  static bool plot_surf;
  static int  surf_int;
  static std::string surf_filename;

  // EB-specific keywords
  static amrex::Real eb_weight;
  static bool eb_redistribution;
  static std::string eb_redistribution_type;

  // LES -variables
  static bool use_LES;

 public:
  PROB::ProbRHS prob_rhs{};   // per-level RHS object (Euler + diffusive + source functors)
  static PROB::ProbClosures* h_prob_closures;
  static PROB::ProbClosures* d_prob_closures;
  static PROB::ProbParm* h_prob_parm;       // host-resident objects used on CPU and as sources for copies
  static PROB::ProbParm* d_prob_parm;       // device-resident objects used on GPU (copied from host at initialization)
  static BCRec* h_phys_bc;
  static BCRec* d_phys_bc;
};

void cns_bcfill(amrex::Box const& bx, amrex::FArrayBox& data, const int dcomp,
                const int numcomp, amrex::Geometry const& geom,
                const amrex::Real time, const amrex::Vector<amrex::BCRec>& bcr,
                const int bcomp, const int scomp);

// declare main IB class instance
#ifdef AMREX_USE_GPIBM  
namespace IBM{
  inline PROB::ProbIB ib;
}
#endif

// declare main EB class instance
#ifdef CNS_USE_EB
namespace EBM{
  inline PROB::ProbEB eb;
}
#endif


#endif
