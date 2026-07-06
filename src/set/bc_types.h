#ifndef BCTYPES_H_
#define BCTYPES_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>

#include <cmath>


// shortcuts for different manual bc to be used within bcnormal in prob.h

// to activate add: 

template <typename cls_t>
class manual_bc_t
{
  public:

  struct nscbc_lodi_bc_parm_t {
    Real rho_inf = Real(1.0);
    Real u_inf = Real(0.0);
    Real v_inf = Real(0.0);
    Real w_inf = Real(0.0);
    Real p_inf = Real(1.0);

    // sigma in L_in = sigma (p - p_inf) - beta T_in.  The dimensional
    // coefficient is pressure_relax * c / l_ref.
    Real pressure_relax = Real(0.0);
    Real transverse_relax = Real(1.0);
    Real l_ref = Real(0.0);

    // Robustness controls for ghost-cell reconstruction.  These are local
    // limiters, not problem-specific wave-amplitude clamps.
    Real min_rho = Real(1.e-12);
    Real min_p = Real(1.e-12);
    Real max_rel_thermo_change = Real(0.50);
    Real max_vel_change_c = Real(2.0);
    Real acoustic_gate_threshold = Real(1.05);

    // Physical sensors for production mixed flows.  Strong compression shocks
    // and pressure-poor velocity/vorticity/wake packets should leave the domain
    // by low-constraint convection/extrapolation, not by acoustic LODI forcing.
    // Entropy/contact content is kept in the outgoing entropy characteristic.
    Real shock_pressure_jump = Real(1.e-2);
    Real shock_compression = Real(1.e-3);
    Real acoustic_coherence_min = Real(0.35);
    Real convective_velocity_jump = Real(1.e-5);
    Real convective_pressure_fraction = Real(0.25);
    Real pressure_relax_min_weight = Real(0.10);

    int use_transverse = 1;
    int use_acoustic_gate = 1;
    int use_ghost_projection = 1;
    int use_shock_sensor = 1;
    int use_convective_sensor = 1;
    int use_convective_extrapolation = 0;
    int use_pressure_relax_sensor = 1;
  };

    manual_bc_t() {}

    ~manual_bc_t() {}

  // .......................................................................//
  // \brief  fixed mass flow rate at given temperature anc composition
  //         pressure will adapt in ghost point 
  // \param  nx,ny,nz : normal face  pointing into domain
  // \param  rhoUfix  : fixed mass flow rate per unit area [kg / m2 s] 
  // \param  Tfix     : fixed Temperature  [K]
  // \param  Yfix     : fixed composition (mass fraction array)
  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_inlet_fixmassflow(
  const Real nx, const Real ny, Real nz, const cls_t* cls,
  const Real rhoUfix, const Real Tfix, const Real* Yfix,
  const Real Uinner[cls_t::NCONS], Real* Ughost )
  {
    // calculate primitive array
    Real Qinner[cls_t::NPRIM]={0.0};      
    cls->cons2prims_point(Uinner,Qinner);

    // pghost = pinner
    const Real P = Qinner[cls_t::QPRES];            

    // ensure mass flow rate is constant - no acoustic reflection at inlet      
    Ughost[cls_t::UMX] = nx*rhoUfix;
    Ughost[cls_t::UMY] = ny*rhoUfix;
    Ughost[cls_t::UMZ] = nz*rhoUfix;
        
    // adjust density ghost to match fix temperature with fix P (and Y)
#if NUM_SPECIES > 1       
    Real rho=0.0; cls->PYT2R(P,Yfix,Tfix,rho); // recalculate rho in case of multiple species
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Ughost[cls_t::UFS+n] = rho*Yfix[n];		
    } 
#else
    Real rho  = Qinner[cls_t::QRHO] * Qinner[cls_t::QT]/Tfix;  
    Ughost[cls_t::URHO] = rho;
#endif     
    // internal specific energy ghost point
    Real e_ext    = 0.0; cls->RYP2E(rho, Yfix, P, e_ext);      
    // kinetic energy ghost point
    Real rhoe_kin = 0.5*( Ughost[cls_t::UMX]*Ughost[cls_t::UMX] + 
                          Ughost[cls_t::UMY]*Ughost[cls_t::UMY] +
                          Ughost[cls_t::UMZ]*Ughost[cls_t::UMZ])/rho; 
    Ughost[cls_t::UET] = rho* e_ext + rhoe_kin;    
  }
  // .......................................................................//
  // \brief  fixed pressure
  // \param  nx,ny,nz: normal face  pointing into domain
  // \param  P0      : fix pressure [Pa]

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_fixP(
    const Real nx, const Real ny, Real nz, const cls_t* cls,
    const Real P0,
    const Real Uinner[cls_t::NCONS], Real* Ughost )
    {
    // calculate primitive array
    Real Qinner[cls_t::NPRIM]={0.0};      
    cls->cons2prims_point(Uinner,Qinner);

    // pghost = pfix
    const Real P = P0;      

    // T and species 0-gradient
    const Real T = Qinner[cls_t::QT];
    Real Y[NUM_SPECIES] = {1.0};
#if NUM_SPECIES > 1      
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = Qinner[cls_t::QFS+n];		
    } 
#endif
    // compute density and internal energy
    Real rho = 0.0; Real e_ext=0.0;
    cls->PYT2R(P,Y,T,rho); cls->RYP2E(rho,Y,P,e_ext);
  
    // assign rho
#if NUM_SPECIES > 1      
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Ughost[cls_t::UFS+n] = rho*Y[n];		
    } 
#else
    Ughost[cls_t::URHO] = rho;
#endif  
    //  assign momentum
    Ughost[cls_t::UMX] = rho*Qinner[cls_t::QU];
    Ughost[cls_t::UMY] = rho*Qinner[cls_t::QV];
    Ughost[cls_t::UMZ] = rho*Qinner[cls_t::QW];
  
    // kinetic energy ghost point
    Real rhoe_kin = 0.5*( Ughost[cls_t::UMX]*Ughost[cls_t::UMX] + 
                        Ughost[cls_t::UMY]*Ughost[cls_t::UMY] +
                        Ughost[cls_t::UMZ]*Ughost[cls_t::UMZ])/rho; 
    Ughost[cls_t::UET] = rho* e_ext + rhoe_kin;    
  }
  // .......................................................................//
  // \brief  subsonic outflow fixed pressure at outlet
  // \param  nx,ny,nz: normal face  pointing into domain
  // \param  P0      : outlet pressure
  // poor's man NSBC based on Whitfield et al. AIAA-84-1552(1984)
  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_subsonic_outflow_fixP(
    const Real nx, const Real ny, Real nz, const cls_t* cls,
    const Real P0,
    const Real Uinner[cls_t::NCONS], Real* Ughost )
    {
	
    // calculate primitive array
    Real Qinner[cls_t::NPRIM]={0.0};      
    cls->cons2prims_point(Uinner,Qinner);

    // pghost = poutlet
    const Real P    = P0;   
    // inner P and rho
    const Real Pi   = Qinner[cls_t::QPRES];  // pressure
    const Real one_over_rho0 = 1.0/Qinner[cls_t::QRHO];   // rho
    const Real one_over_c0   = 1.0/Qinner[cls_t::QC];     // sound speed

    const Real rho = Qinner[cls_t::QRHO]+ (P -Pi)*one_over_c0*one_over_c0;

    // Pout > Pinside : if nx> 0  du>0   , if nx < 0  du<0
    const Real u = Qinner[cls_t::QU] - (P -Pi)*one_over_c0*nx*one_over_rho0;
    const Real v = Qinner[cls_t::QV] - (P -Pi)*one_over_c0*ny*one_over_rho0;
    const Real w = Qinner[cls_t::QW] - (P -Pi)*one_over_c0*nz*one_over_rho0;

    // convert back to conservative vars

    //  assign momentum
    Ughost[cls_t::UMX] = rho*u;
    Ughost[cls_t::UMY] = rho*v;
    Ughost[cls_t::UMZ] = rho*w;
  
    // assign rho
    Real Y[NUM_SPECIES] = {1.0};
#if NUM_SPECIES > 1      
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = Qinner[cls_t::QFS+n];
      Ughost[cls_t::UFS+n] = rho*Y[n];		      
    } 
#else
    Ughost[cls_t::URHO] = rho;
#endif  

    // compute and assign energy
    const Real e_kin = 0.5*(u*u+ v*v + w*w);
    Real e_ext=0.0; cls->RYP2E(rho,Y,P,e_ext);
  
    Ughost[cls_t::UET] = rho* e_ext + rho*e_kin;   


  }

  // .......................................................................//
  // \brief  Local characteristic far-field boundary condition.
  //
  // This is an algebraic ghost-cell NSCBC-style closure for the inviscid
  // normal operator. The normal vector points outward from the domain.
  //
  // supersonic inflow:  all characteristics prescribed from freestream
  // supersonic outflow: all characteristics extrapolated from interior
  // subsonic inflow:    one acoustic characteristic extrapolated from interior,
  //                     remaining incoming data prescribed from freestream
  // subsonic outflow:   incoming acoustic characteristic prescribed from
  //                     freestream, outgoing acoustic/entropy/tangential data
  //                     extrapolated from interior
  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_nscbc_farfield(
    const Real nx, const Real ny, const Real nz, const cls_t* cls,
    const Real rho_inf, const Real u_inf, const Real v_inf,
    const Real w_inf, const Real p_inf,
    const Real Uinner[cls_t::NCONS], Real* Ughost,
    const Real* Yinf = nullptr)
  {
    Real Qinner[cls_t::NPRIM] = {0.0};
    cls->cons2prims_point(Uinner, Qinner);

    Real Yinner[NUM_SPECIES] = {1.0};
    Real Yfar[NUM_SPECIES] = {1.0};
#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yinner[n] = Qinner[cls_t::QFS + n];
      Yfar[n] = (Yinf == nullptr) ? Yinner[n] : Yinf[n];
    }
#endif

    const Real eps = Real(1.e-14);
    const Real nmag = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (nmag <= eps) {
      for (int n = 0; n < cls_t::NCONS; ++n) {
        Ughost[n] = Uinner[n];
      }
      return;
    }

    const Real nnx = nx / nmag;
    const Real nny = ny / nmag;
    const Real nnz = nz / nmag;

    const Real rhoi = amrex::max(Qinner[cls_t::QRHO], eps);
    const Real pi   = amrex::max(Qinner[cls_t::QPRES], eps);
    const Real ui   = Qinner[cls_t::QU];
    const Real vi   = Qinner[cls_t::QV];
    const Real wi   = Qinner[cls_t::QW];
    const Real ci   = amrex::max(Qinner[cls_t::QC], eps);
    const Real gam  = amrex::max(Qinner[cls_t::QG], Real(1.0) + eps);
    const Real gm1  = amrex::max(gam - Real(1.0), eps);

    const Real rhof = amrex::max(rho_inf, eps);
    const Real pf   = amrex::max(p_inf, eps);
    const Real cf   = std::sqrt(amrex::max(gam * pf / rhof, eps));

    const Real uni = ui * nnx + vi * nny + wi * nnz;
    const Real unf = u_inf * nnx + v_inf * nny + w_inf * nnz;

    Real rho = rhoi;
    Real p = pi;
    Real un = uni;
    Real utx = ui - uni * nnx;
    Real uty = vi - uni * nny;
    Real utz = wi - uni * nnz;
    const Real* Yout = Yinner;

    if (uni <= -ci) {
      rho = rhof;
      p = pf;
      un = unf;
      utx = u_inf - unf * nnx;
      uty = v_inf - unf * nny;
      utz = w_inf - unf * nnz;
      Yout = Yfar;
    } else if (uni >= ci) {
      for (int n = 0; n < cls_t::NCONS; ++n) {
        Ughost[n] = Uinner[n];
      }
      return;
    } else if (uni < Real(0.0)) {
      const Real jplus_i = uni + Real(2.0) * ci / gm1;
      const Real jminus_f = unf - Real(2.0) * cf / gm1;
      const Real cb = amrex::max(Real(0.25) * gm1 * (jplus_i - jminus_f), eps);
      const Real sfar = pf / std::pow(rhof, gam);

      un = Real(0.5) * (jplus_i + jminus_f);
      rho = std::pow(cb * cb / (gam * sfar), Real(1.0) / gm1);
      p = rho * cb * cb / gam;
      utx = u_inf - unf * nnx;
      uty = v_inf - unf * nny;
      utz = w_inf - unf * nnz;
      Yout = Yfar;
    } else {
      const Real jplus_i = uni + Real(2.0) * ci / gm1;
      const Real jminus_f = unf - Real(2.0) * cf / gm1;
      const Real cb = amrex::max(Real(0.25) * gm1 * (jplus_i - jminus_f), eps);
      const Real sint = pi / std::pow(rhoi, gam);

      un = Real(0.5) * (jplus_i + jminus_f);
      rho = std::pow(cb * cb / (gam * sint), Real(1.0) / gm1);
      p = rho * cb * cb / gam;
    }

    const Real u = utx + un * nnx;
    const Real v = uty + un * nny;
    const Real w = utz + un * nnz;

#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Ughost[cls_t::UFS + n] = rho * Yout[n];
    }
#else
    Ughost[cls_t::URHO] = rho;
#endif
    Ughost[cls_t::UMX] = rho * u;
    Ughost[cls_t::UMY] = rho * v;
    Ughost[cls_t::UMZ] = rho * w;

    Real e_ext = Real(0.0);
    cls->RYP2E(rho, Yout, p, e_ext);
    Ughost[cls_t::UET] =
      rho * e_ext + Real(0.5) * rho * (u * u + v * v + w * w);
  }

  template <typename Array4T>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void load_prims_point(
    const IntVect& iv, Array4T const& state, const cls_t* cls,
    Real q[cls_t::NPRIM])
  {
    Real u[cls_t::NCONS] = {Real(0.0)};
    for (int n = 0; n < cls_t::NCONS; ++n) {
      u[n] = state(iv, n);
    }
    cls->cons2prims_point(u, q);
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real clamp_value(
    const Real x, const Real lo, const Real hi)
  {
    return amrex::min(hi, amrex::max(lo, x));
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cons_from_prims(
    const cls_t* cls, const Real rho, const Real u, const Real v,
    const Real w, const Real p, const Real Y[NUM_SPECIES], Real* U)
  {
#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
      U[cls_t::UFS + n] = rho * Y[n];
    }
#else
    U[cls_t::URHO] = rho;
#endif
    U[cls_t::UMX] = rho * u;
    U[cls_t::UMY] = rho * v;
    U[cls_t::UMZ] = rho * w;

    Real e_ext = Real(0.0);
    cls->RYP2E(rho, Y, p, e_ext);
    U[cls_t::UET] =
      rho * e_ext + Real(0.5) * rho * (u * u + v * v + w * w);
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cons_rhs_from_prim_rhs(
    const Real rho, const Real u, const Real v, const Real w,
    const Real p, const Real gamma, const Real rho_t, const Real u_t,
    const Real v_t, const Real w_t, const Real p_t,
    const Real Y[NUM_SPECIES], const Real Y_t[NUM_SPECIES], Real* rhs)
  {
    const Real gm1 = amrex::max(gamma - Real(1.0), Real(1.e-14));

#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rhs[cls_t::UFS + n] = rho_t * Y[n] + rho * Y_t[n];
    }
#else
    rhs[cls_t::URHO] = rho_t;
#endif

    rhs[cls_t::UMX] = rho_t * u + rho * u_t;
    rhs[cls_t::UMY] = rho_t * v + rho * v_t;
    rhs[cls_t::UMZ] = rho_t * w + rho * w_t;
    rhs[cls_t::UET] =
      p_t / gm1 +
      Real(0.5) * rho_t * (u * u + v * v + w * w) +
      rho * (u * u_t + v * v_t + w * w_t);
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void prim_rhs_from_cons_rhs(
    const Real rho, const Real u, const Real v, const Real w,
    const Real gamma, const Real* rhs, const Real Y[NUM_SPECIES],
    Real& rho_t, Real& u_t, Real& v_t, Real& w_t, Real& p_t,
    Real Y_t[NUM_SPECIES])
  {
    const Real rho_safe = amrex::max(rho, Real(1.e-14));
    const Real gm1 = amrex::max(gamma - Real(1.0), Real(1.e-14));

#if NUM_SPECIES > 1
    rho_t = Real(0.0);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho_t += rhs[cls_t::UFS + n];
    }
#else
    rho_t = rhs[cls_t::URHO];
#endif

    u_t = (rhs[cls_t::UMX] - u * rho_t) / rho_safe;
    v_t = (rhs[cls_t::UMY] - v * rho_t) / rho_safe;
    w_t = (rhs[cls_t::UMZ] - w * rho_t) / rho_safe;

    p_t = gm1 *
      (rhs[cls_t::UET] -
       Real(0.5) * rho_t * (u * u + v * v + w * w) -
       rho * (u * u_t + v * v_t + w * w_t));

#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_t[n] = (rhs[cls_t::UFS + n] - Y[n] * rho_t) / rho_safe;
    }
#else
    (void)Y;
    Y_t[0] = Real(0.0);
#endif
  }

  // .......................................................................//
  // \brief  Multidimensional characteristic-amplitude LODI far-field ghost
  //         fill for Cartesian physical boundaries.
  //
  // This is the production-oriented NSCBC/LODI helper.  It uses a local
  // normal/tangential characteristic decomposition of the inviscid operator:
  //
  //   L+ = (u_n+c)(p_n + rho c u_{n,n})
  //   L- = (u_n-c)(p_n - rho c u_{n,n})
  //   Ls = u_n(rho_n - p_n/c^2)
  //   Lt = u_n u_{t,n}
  //
  // and includes transverse source terms in the incoming acoustic amplitude.
  // Supersonic and inflow cases are delegated to the fixed characteristic
  // far-field closure.  Subsonic outflow keeps outgoing entropy/tangential
  // gradients and prescribes the incoming acoustic amplitude by pressure
  // relaxation plus transverse relaxation.
  template <typename Array4T, typename GeomDataT>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool bc_nscbc_lodi_farfield(
    const IntVect& iv, Array4T const& state, Real* Ughost,
    const int idir, const int sgn, GeomDataT const& geomdata,
    const cls_t* cls, const nscbc_lodi_bc_parm_t& bp,
    const Real* Yinf = nullptr)
  {
    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();

    IntVect first(iv);
    first[idir] = (sgn == 1) ? domlo[idir] : domhi[idir];

    IntVect inward(first);
    inward[idir] += sgn;
    if (inward[idir] < domlo[idir] || inward[idir] > domhi[idir]) {
      Real Uinner[cls_t::NCONS] = {Real(0.0)};
      for (int n = 0; n < cls_t::NCONS; ++n) Uinner[n] = state(first, n);
      const Real normal[3] = {
        AMREX_D_DECL(idir == 0 ? -Real(sgn) : Real(0.0),
                     idir == 1 ? -Real(sgn) : Real(0.0),
                     idir == 2 ? -Real(sgn) : Real(0.0))};
      bc_nscbc_farfield(normal[0], normal[1], normal[2], cls,
                        bp.rho_inf, bp.u_inf, bp.v_inf, bp.w_inf, bp.p_inf,
                        Uinner, Ughost, Yinf);
      return true;
    }

    Real qi[cls_t::NPRIM] = {Real(0.0)};
    Real qin[cls_t::NPRIM] = {Real(0.0)};
    load_prims_point(first, state, cls, qi);
    load_prims_point(inward, state, cls, qin);

    Real Uinner[cls_t::NCONS] = {Real(0.0)};
    for (int n = 0; n < cls_t::NCONS; ++n) Uinner[n] = state(first, n);

    const Real nsgn = -Real(sgn);
    const Real eps = Real(1.e-14);
    const Real rho = amrex::max(qi[cls_t::QRHO], bp.min_rho);
    const Real p = amrex::max(qi[cls_t::QPRES], bp.min_p);
    const Real c = amrex::max(qi[cls_t::QC], eps);
    const Real vel[3] = {qi[cls_t::QU], qi[cls_t::QV], qi[cls_t::QW]};
    const Real un = nsgn * vel[idir];

    const Real normal[3] = {
      AMREX_D_DECL(idir == 0 ? nsgn : Real(0.0),
                   idir == 1 ? nsgn : Real(0.0),
                   idir == 2 ? nsgn : Real(0.0))};

    if (un <= Real(0.0) || un >= c) {
      bc_nscbc_farfield(normal[0], normal[1], normal[2], cls,
                        bp.rho_inf, bp.u_inf, bp.v_inf, bp.w_inf, bp.p_inf,
                        Uinner, Ughost, Yinf);
      return true;
    }

    const Real inv_dn = Real(1.0) / dx[idir];
    const Real rho_n_i = (qi[cls_t::QRHO] - qin[cls_t::QRHO]) * inv_dn;
    const Real p_n_i = (qi[cls_t::QPRES] - qin[cls_t::QPRES]) * inv_dn;
    Real vel_n_i[3] = {Real(0.0), Real(0.0), Real(0.0)};
    vel_n_i[0] = (qi[cls_t::QU] - qin[cls_t::QU]) * inv_dn;
    vel_n_i[1] = (qi[cls_t::QV] - qin[cls_t::QV]) * inv_dn;
    vel_n_i[2] = (qi[cls_t::QW] - qin[cls_t::QW]) * inv_dn;
    const Real un_n_i = nsgn * vel_n_i[idir];

    Real Y_i[NUM_SPECIES] = {Real(1.0)};
    Real Y_n_i[NUM_SPECIES] = {Real(0.0)};
    Real Y_far[NUM_SPECIES] = {Real(1.0)};
#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_i[n] = qi[cls_t::QFS + n];
      Y_far[n] = (Yinf == nullptr) ? Y_i[n] : Yinf[n];
      Y_n_i[n] = (qi[cls_t::QFS + n] - qin[cls_t::QFS + n]) * inv_dn;
    }
#else
    Y_i[0] = Real(1.0);
    Y_far[0] = Real(1.0);
#endif

    Real tan_adv_p = Real(0.0);
    Real tan_adv_un = Real(0.0);
    Real tan_div_u = Real(0.0);
    Real gradp_t2 = Real(0.0);
    Real grad_vel_t2 = Real(0.0);

    for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
      if (tdir == idir) continue;
      if (first[tdir] <= domlo[tdir] || first[tdir] >= domhi[tdir]) continue;

      IntVect jp(first);
      IntVect jm(first);
      jp[tdir] += 1;
      jm[tdir] -= 1;

      Real qp[cls_t::NPRIM] = {Real(0.0)};
      Real qm[cls_t::NPRIM] = {Real(0.0)};
      load_prims_point(jp, state, cls, qp);
      load_prims_point(jm, state, cls, qm);

      const Real inv_2dt = Real(0.5) / dx[tdir];
      const Real p_t = (qp[cls_t::QPRES] - qm[cls_t::QPRES]) * inv_2dt;
      const Real qp_unvel =
        (idir == 0 ? qp[cls_t::QU] : (idir == 1 ? qp[cls_t::QV] : qp[cls_t::QW]));
      const Real qm_unvel =
        (idir == 0 ? qm[cls_t::QU] : (idir == 1 ? qm[cls_t::QV] : qm[cls_t::QW]));
      const Real un_t = nsgn * (qp_unvel - qm_unvel) * inv_2dt;
      const Real ut_t =
        ((tdir == 0 ? qp[cls_t::QU] : (tdir == 1 ? qp[cls_t::QV] : qp[cls_t::QW])) -
         (tdir == 0 ? qm[cls_t::QU] : (tdir == 1 ? qm[cls_t::QV] : qm[cls_t::QW]))) *
        inv_2dt;

      const Real ut = vel[tdir];
      tan_adv_p += ut * p_t;
      tan_adv_un += ut * un_t;
      tan_div_u += ut_t;
      gradp_t2 += p_t * p_t;
      const Real u_t =
        (qp[cls_t::QU] - qm[cls_t::QU]) * inv_2dt;
      const Real v_t =
        (qp[cls_t::QV] - qm[cls_t::QV]) * inv_2dt;
      const Real w_t =
        (qp[cls_t::QW] - qm[cls_t::QW]) * inv_2dt;
      grad_vel_t2 += u_t * u_t + v_t * v_t + w_t * w_t;
    }

    const Real trans_minus =
      tan_adv_p + rho * c * c * tan_div_u - rho * c * tan_adv_un;

    const Real l_plus = (un + c) * (p_n_i + rho * c * un_n_i);
    const Real entropy_n_i = rho_n_i - p_n_i / (c * c);

    const Real lref =
      (bp.l_ref > Real(0.0)) ? bp.l_ref : (prob_hi[idir] - prob_lo[idir]);
    const Real sigma = bp.pressure_relax * c / amrex::max(lref, eps);
    const Real p_scale = amrex::max(p, bp.p_inf);
    const Real modal_floor =
      eps * p_scale / amrex::max(lref, eps);
    const Real dist = amrex::Math::abs(Real(iv[idir] - first[idir])) * dx[idir];

    const Real gradp_t = std::sqrt(gradp_t2);
    const Real grad_vel_t = std::sqrt(grad_vel_t2);
    const Real div_u = un_n_i + tan_div_u;
    const Real pressure_jump =
      amrex::max(amrex::Math::abs(p_n_i), gradp_t) * dx[idir] /
      amrex::max(p_scale, eps);
    const Real compression_jump =
      amrex::max(-div_u, Real(0.0)) * dx[idir] / c;
    if (bp.use_shock_sensor &&
        pressure_jump > bp.shock_pressure_jump &&
        compression_jump > bp.shock_compression) {
      for (int n = 0; n < cls_t::NCONS; ++n) Ughost[n] = Uinner[n];
      return true;
    }

    Real tangential_vel_n2 = Real(0.0);
    for (int vdir = 0; vdir < AMREX_SPACEDIM; ++vdir) {
      if (vdir == idir) continue;
      tangential_vel_n2 += vel_n_i[vdir] * vel_n_i[vdir];
    }
    const Real acoustic_p = amrex::Math::abs(p_n_i) * dx[idir];
    const Real acoustic_u = rho * c * amrex::Math::abs(un_n_i) * dx[idir];
    const Real acoustic_total = acoustic_p + acoustic_u;
    const Real acoustic_coherence =
      (acoustic_total > modal_floor * dx[idir])
      ? Real(2.0) * amrex::min(acoustic_p, acoustic_u) / acoustic_total
      : Real(0.0);
    const Real tangential_strength =
      rho * c * std::sqrt(tangential_vel_n2) * dx[idir];
    const Real transverse_vortical_strength =
      rho * c * grad_vel_t * dx[idir] * dx[idir] / amrex::max(lref, dx[idir]);
    const Real velocity_packet_strength =
      amrex::max(acoustic_u, tangential_strength + transverse_vortical_strength);
    const Real velocity_packet_floor =
      bp.convective_velocity_jump * p_scale;
    const bool pressure_poor_velocity_packet =
      acoustic_coherence < bp.acoustic_coherence_min &&
      velocity_packet_strength > velocity_packet_floor &&
      acoustic_p < bp.convective_pressure_fraction * velocity_packet_strength;
    if (bp.use_convective_sensor && pressure_poor_velocity_packet) {
      if (!bp.use_convective_extrapolation) {
        for (int n = 0; n < cls_t::NCONS; ++n) Ughost[n] = Uinner[n];
        return true;
      }

      const Real rho_scale = amrex::max(rho, bp.rho_inf);
      const Real vel_limit = bp.max_vel_change_c * c;
      const Real rho_delta =
        clamp_value(dist * rho_n_i,
                    -bp.max_rel_thermo_change * rho_scale,
                     bp.max_rel_thermo_change * rho_scale);
      const Real p_delta =
        clamp_value(dist * p_n_i,
                    -bp.max_rel_thermo_change * p_scale,
                     bp.max_rel_thermo_change * p_scale);

      const Real rho_g = amrex::max(rho + rho_delta, bp.min_rho);
      const Real p_g = amrex::max(p + p_delta, bp.min_p);
      Real vel_g[3] = {
        vel[0] + clamp_value(dist * vel_n_i[0], -vel_limit, vel_limit),
        vel[1] + clamp_value(dist * vel_n_i[1], -vel_limit, vel_limit),
        vel[2] + clamp_value(dist * vel_n_i[2], -vel_limit, vel_limit)};

      Real Y_g[NUM_SPECIES] = {Real(1.0)};
#if NUM_SPECIES > 1
      Real ysum = Real(0.0);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y_g[n] = clamp_value(Y_i[n] + dist * Y_n_i[n], Real(0.0), Real(1.0));
        ysum += Y_g[n];
      }
      if (ysum > eps) {
        for (int n = 0; n < NUM_SPECIES; ++n) Y_g[n] /= ysum;
      } else {
        Real yfar_sum = Real(0.0);
        for (int n = 0; n < NUM_SPECIES; ++n) yfar_sum += Y_far[n];
        if (yfar_sum > eps) {
          for (int n = 0; n < NUM_SPECIES; ++n) Y_g[n] = Y_far[n] / yfar_sum;
        } else {
          Y_g[0] = Real(1.0);
          for (int n = 1; n < NUM_SPECIES; ++n) Y_g[n] = Real(0.0);
        }
      }
#else
      Y_g[0] = Real(1.0);
#endif

      cons_from_prims(cls, rho_g, vel_g[0], vel_g[1], vel_g[2], p_g, Y_g, Ughost);
      return true;
    }

    Real acoustic_gate = Real(1.0);
    Real acoustic_projection = Real(1.0);
    if (bp.use_acoustic_gate || bp.use_ghost_projection) {
      const Real grad_floor =
        eps * amrex::max(p, bp.p_inf) / amrex::max(lref, eps);
      Real sin_theta = rho * c * amrex::Math::abs(tan_div_u) /
                       (gradp_t + grad_floor);
      if (bp.use_acoustic_gate &&
          sin_theta > amrex::max(bp.acoustic_gate_threshold, Real(1.0))) {
        acoustic_gate = Real(0.0);
      }
      sin_theta = clamp_value(sin_theta, Real(0.0), Real(0.98));
      const Real cos_theta =
        std::sqrt(amrex::max(Real(1.0) - sin_theta * sin_theta, Real(1.e-8)));
      if (bp.use_ghost_projection) {
        acoustic_projection =
          amrex::max(c - un, Real(0.0)) / c * cos_theta /
          (Real(1.0) + cos_theta);
      }
    }

    const Real transverse_coeff =
      (bp.use_transverse ? bp.transverse_relax : Real(0.0)) *
      acoustic_gate * acoustic_projection;
    const Real relax_weight =
      bp.use_pressure_relax_sensor
      ? amrex::max(bp.pressure_relax_min_weight, acoustic_coherence)
      : Real(1.0);
    const Real l_minus = sigma * relax_weight * (p - bp.p_inf) -
                         transverse_coeff * trans_minus;

    const Real a_minus = l_minus / amrex::min(un - c, -eps);
    const Real a_plus = l_plus / amrex::max(un + c, eps);
    const Real p_n = Real(0.5) * (a_minus + a_plus);
    const Real un_n = (a_plus - a_minus) / (Real(2.0) * rho * c);
    const Real rho_n = entropy_n_i + p_n / (c * c);

    Real vel_n[3] = {vel_n_i[0], vel_n_i[1], vel_n_i[2]};
    vel_n[idir] = nsgn * un_n;

    const Real rho_scale = amrex::max(rho, bp.rho_inf);
    const Real vel_limit = bp.max_vel_change_c * c;
    const Real rho_delta =
      clamp_value(dist * rho_n,
                  -bp.max_rel_thermo_change * rho_scale,
                   bp.max_rel_thermo_change * rho_scale);
    const Real p_delta =
      clamp_value(dist * p_n,
                  -bp.max_rel_thermo_change * p_scale,
                   bp.max_rel_thermo_change * p_scale);

    const Real rho_g = amrex::max(rho + rho_delta, bp.min_rho);
    const Real p_g = amrex::max(p + p_delta, bp.min_p);
    Real vel_g[3] = {
      vel[0] + clamp_value(dist * vel_n[0], -vel_limit, vel_limit),
      vel[1] + clamp_value(dist * vel_n[1], -vel_limit, vel_limit),
      vel[2] + clamp_value(dist * vel_n[2], -vel_limit, vel_limit)};

    Real Y_g[NUM_SPECIES] = {Real(1.0)};
#if NUM_SPECIES > 1
    Real ysum = Real(0.0);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_g[n] = clamp_value(Y_i[n] + dist * Y_n_i[n], Real(0.0), Real(1.0));
      ysum += Y_g[n];
    }
    if (ysum > eps) {
      for (int n = 0; n < NUM_SPECIES; ++n) Y_g[n] /= ysum;
    } else {
      Real yfar_sum = Real(0.0);
      for (int n = 0; n < NUM_SPECIES; ++n) yfar_sum += Y_far[n];
      if (yfar_sum > eps) {
        for (int n = 0; n < NUM_SPECIES; ++n) Y_g[n] = Y_far[n] / yfar_sum;
      } else {
        Y_g[0] = Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Y_g[n] = Real(0.0);
      }
    }
#else
    Y_g[0] = Real(1.0);
#endif

    cons_from_prims(cls, rho_g, vel_g[0], vel_g[1], vel_g[2], p_g, Y_g, Ughost);
    return true;
  }

  // .......................................................................//
  // \brief  RHS-level multidimensional LODI/NSCBC boundary update.
  //
  // This controls the incoming acoustic characteristic in the first interior
  // boundary cell for a subsonic outflow.  The pre-existing finite-volume RHS
  // is first converted to primitive time derivatives; only the L- acoustic
  // amplitude is corrected, preserving the FV outgoing acoustic, entropy,
  // tangential, and mixed-flow content as much as possible.  Transverse terms
  // enter the boundary characteristic relation explicitly:
  //
  //   q_t + A_n q_n = - sum_t A_t q_tan + S_relax .
  //
  // Current scope: inviscid characteristic operator for Cartesian faces.  If a
  // shock/wake sensor fires, the function leaves the finite-volume RHS already
  // computed by the solver in place.
  template <typename PrimArrayT, typename RhsArrayT>
  static void rhs_nscbc_lodi_farfield(
    const Box& bx, PrimArrayT const& prims, RhsArrayT const& rhs,
    const int idir, const int sgn, Geometry const& geom,
    const cls_t* cls, const nscbc_lodi_bc_parm_t& bp)
  {
    (void)cls;

    Box face = bx & geom.Domain();
    if (!face.ok()) return;

    const auto geomdata = geom.data();
    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();
    const int ib = (sgn == 1) ? domlo[idir] : domhi[idir];
    if (face.smallEnd(idir) > ib || face.bigEnd(idir) < ib) return;
    face.setSmall(idir, ib);
    face.setBig(idir, ib);

    amrex::ParallelFor(face, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int* dlo = geomdata.Domain().loVect();
      const int* dhi = geomdata.Domain().hiVect();
      const Real* prob_lo = geomdata.ProbLo();
      const Real* prob_hi = geomdata.ProbHi();
      const Real* dx = geomdata.CellSize();

      IntVect first(AMREX_D_DECL(i, j, k));
      IntVect inward(first);
      inward[idir] += sgn;
      if (inward[idir] < dlo[idir] || inward[idir] > dhi[idir]) return;

      Real qi[cls_t::NPRIM] = {Real(0.0)};
      Real qin[cls_t::NPRIM] = {Real(0.0)};
      for (int n = 0; n < cls_t::NPRIM; ++n) {
        qi[n] = prims(first, n);
        qin[n] = prims(inward, n);
      }

      const Real nsgn = -Real(sgn);
      const Real eps = Real(1.e-14);
      const Real rho = amrex::max(qi[cls_t::QRHO], bp.min_rho);
      const Real p = amrex::max(qi[cls_t::QPRES], bp.min_p);
      const Real c = amrex::max(qi[cls_t::QC], eps);
      const Real gamma = amrex::max(qi[cls_t::QG], Real(1.0) + eps);
      const Real vel[3] = {qi[cls_t::QU], qi[cls_t::QV], qi[cls_t::QW]};
      const Real un = nsgn * vel[idir];

      if (un <= Real(0.0) || un >= c) return;

      const Real inv_dn = Real(1.0) / dx[idir];
      const Real p_n_i = (qi[cls_t::QPRES] - qin[cls_t::QPRES]) * inv_dn;
      Real vel_n_i[3] = {Real(0.0), Real(0.0), Real(0.0)};
      vel_n_i[0] = (qi[cls_t::QU] - qin[cls_t::QU]) * inv_dn;
      vel_n_i[1] = (qi[cls_t::QV] - qin[cls_t::QV]) * inv_dn;
      vel_n_i[2] = (qi[cls_t::QW] - qin[cls_t::QW]) * inv_dn;
      const Real un_n_i = nsgn * vel_n_i[idir];

      Real Y_i[NUM_SPECIES] = {Real(1.0)};
#if NUM_SPECIES > 1
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y_i[n] = qi[cls_t::QFS + n];
      }
#else
      Y_i[0] = Real(1.0);
#endif

      Real tan_adv_p = Real(0.0);
      Real tan_adv_un = Real(0.0);
      Real tan_div_u = Real(0.0);
      Real gradp_t2 = Real(0.0);
      Real grad_vel_t2 = Real(0.0);

      for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
        if (tdir == idir) continue;
        if (first[tdir] <= dlo[tdir] || first[tdir] >= dhi[tdir]) continue;

        IntVect jp(first);
        IntVect jm(first);
        jp[tdir] += 1;
        jm[tdir] -= 1;

        Real qp[cls_t::NPRIM] = {Real(0.0)};
        Real qm[cls_t::NPRIM] = {Real(0.0)};
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          qp[n] = prims(jp, n);
          qm[n] = prims(jm, n);
        }

        const Real inv_2dt = Real(0.5) / dx[tdir];
        const Real p_tdir = (qp[cls_t::QPRES] - qm[cls_t::QPRES]) * inv_2dt;
        const Real qp_unvel =
          (idir == 0 ? qp[cls_t::QU] : (idir == 1 ? qp[cls_t::QV] : qp[cls_t::QW]));
        const Real qm_unvel =
          (idir == 0 ? qm[cls_t::QU] : (idir == 1 ? qm[cls_t::QV] : qm[cls_t::QW]));
        const Real un_tdir = nsgn * (qp_unvel - qm_unvel) * inv_2dt;
        const Real div_tdir =
          ((tdir == 0 ? qp[cls_t::QU] : (tdir == 1 ? qp[cls_t::QV] : qp[cls_t::QW])) -
           (tdir == 0 ? qm[cls_t::QU] : (tdir == 1 ? qm[cls_t::QV] : qm[cls_t::QW]))) *
          inv_2dt;

        const Real ut = vel[tdir];
        tan_adv_p += ut * p_tdir;
        tan_adv_un += ut * un_tdir;
        tan_div_u += div_tdir;
        gradp_t2 += p_tdir * p_tdir;

        const Real u_tdir = (qp[cls_t::QU] - qm[cls_t::QU]) * inv_2dt;
        const Real v_tdir = (qp[cls_t::QV] - qm[cls_t::QV]) * inv_2dt;
        const Real w_tdir = (qp[cls_t::QW] - qm[cls_t::QW]) * inv_2dt;
        grad_vel_t2 += u_tdir * u_tdir + v_tdir * v_tdir + w_tdir * w_tdir;
      }

      const Real lref =
        (bp.l_ref > Real(0.0)) ? bp.l_ref : (prob_hi[idir] - prob_lo[idir]);
      const Real sigma = bp.pressure_relax * c / amrex::max(lref, eps);
      const Real p_scale = amrex::max(p, bp.p_inf);
      const Real modal_floor =
        eps * p_scale / amrex::max(lref, eps);

      const Real gradp_t = std::sqrt(gradp_t2);
      const Real grad_vel_t = std::sqrt(grad_vel_t2);
      const Real div_u = un_n_i + tan_div_u;
      const Real pressure_jump =
        amrex::max(amrex::Math::abs(p_n_i), gradp_t) * dx[idir] /
        amrex::max(p_scale, eps);
      const Real compression_jump =
        amrex::max(-div_u, Real(0.0)) * dx[idir] / c;
      if (bp.use_shock_sensor &&
          pressure_jump > bp.shock_pressure_jump &&
          compression_jump > bp.shock_compression) {
        return;
      }

      Real tangential_vel_n2 = Real(0.0);
      for (int vdir = 0; vdir < AMREX_SPACEDIM; ++vdir) {
        if (vdir == idir) continue;
        tangential_vel_n2 += vel_n_i[vdir] * vel_n_i[vdir];
      }
      const Real acoustic_p = amrex::Math::abs(p_n_i) * dx[idir];
      const Real acoustic_u = rho * c * amrex::Math::abs(un_n_i) * dx[idir];
      const Real acoustic_total = acoustic_p + acoustic_u;
      const Real acoustic_coherence =
        (acoustic_total > modal_floor * dx[idir])
        ? Real(2.0) * amrex::min(acoustic_p, acoustic_u) / acoustic_total
        : Real(0.0);
      const Real tangential_strength =
        rho * c * std::sqrt(tangential_vel_n2) * dx[idir];
      const Real transverse_vortical_strength =
        rho * c * grad_vel_t * dx[idir] * dx[idir] / amrex::max(lref, dx[idir]);
      const Real velocity_packet_strength =
        amrex::max(acoustic_u, tangential_strength + transverse_vortical_strength);
      const bool pressure_poor_velocity_packet =
        acoustic_coherence < bp.acoustic_coherence_min &&
        velocity_packet_strength > bp.convective_velocity_jump * p_scale &&
        acoustic_p < bp.convective_pressure_fraction * velocity_packet_strength;
      if (bp.use_convective_sensor && pressure_poor_velocity_packet) {
        return;
      }

      Real rhs_fv[cls_t::NCONS] = {Real(0.0)};
      for (int n = 0; n < cls_t::NCONS; ++n) {
        rhs_fv[n] = rhs(i, j, k, n);
      }

      Real rho_time = Real(0.0);
      Real vel_time[3] = {Real(0.0), Real(0.0), Real(0.0)};
      Real p_time = Real(0.0);
      Real Y_time[NUM_SPECIES] = {Real(0.0)};
      prim_rhs_from_cons_rhs(rho, vel[0], vel[1], vel[2], gamma,
                             rhs_fv, Y_i, rho_time, vel_time[0],
                             vel_time[1], vel_time[2], p_time, Y_time);

      const Real trans_p =
        tan_adv_p + rho * c * c * tan_div_u;
      const Real trans_un = tan_adv_un;
      const Real trans_minus = trans_p - rho * c * trans_un;
      Real acoustic_gate = Real(1.0);
      Real acoustic_projection = Real(1.0);
      if (bp.use_acoustic_gate || bp.use_ghost_projection) {
        const Real grad_floor =
          eps * amrex::max(p, bp.p_inf) / amrex::max(lref, eps);
        Real sin_theta = rho * c * amrex::Math::abs(tan_div_u) /
                         (gradp_t + grad_floor);
        if (bp.use_acoustic_gate &&
            sin_theta > amrex::max(bp.acoustic_gate_threshold, Real(1.0))) {
          acoustic_gate = Real(0.0);
        }
        sin_theta = clamp_value(sin_theta, Real(0.0), Real(0.98));
        const Real cos_theta =
          std::sqrt(amrex::max(Real(1.0) - sin_theta * sin_theta, Real(1.e-8)));
        if (bp.use_ghost_projection) {
          acoustic_projection =
            amrex::max(c - un, Real(0.0)) / c * cos_theta /
            (Real(1.0) + cos_theta);
        }
      }
      const Real physical_transverse_coeff =
        bp.use_transverse ? Real(1.0) : Real(0.0);
      const Real target_transverse_coeff =
        (bp.use_transverse ? bp.transverse_relax : Real(0.0)) *
        acoustic_gate * acoustic_projection;
      const Real un_time_fv = nsgn * vel_time[idir];
      const Real l_minus_fv =
        -(p_time + physical_transverse_coeff * trans_p) +
        rho * c * (un_time_fv + physical_transverse_coeff * trans_un);

      const Real relax_weight =
        bp.use_pressure_relax_sensor
        ? amrex::max(bp.pressure_relax_min_weight, acoustic_coherence)
        : Real(1.0);
      const Real l_minus_target =
        sigma * relax_weight * (p - bp.p_inf) -
        target_transverse_coeff * trans_minus;
      const Real delta_lminus = l_minus_target - l_minus_fv;

      const Real delta_p_time = -Real(0.5) * delta_lminus;
      const Real delta_un_time = delta_lminus / (Real(2.0) * rho * c);
      p_time += delta_p_time;
      rho_time += delta_p_time / (c * c);
      vel_time[idir] += nsgn * delta_un_time;

      Real Urhs[cls_t::NCONS] = {Real(0.0)};
      cons_rhs_from_prim_rhs(rho, vel[0], vel[1], vel[2], p, gamma,
                             rho_time, vel_time[0], vel_time[1],
                             vel_time[2], p_time, Y_i, Y_time, Urhs);

      for (int n = 0; n < cls_t::NCONS; ++n) {
        rhs(i, j, k, n) = Urhs[n];
      }
    });
  }

};

#endif
