# Full axisymmetric R-Z Navier--Stokes manufactured solution.
TINY_PROFILE = FALSE
DEBUG        = FALSE
PRECISION    = DOUBLE
DIM          = 2
COMP         = gcc
CXXSTD       = c++20

USE_MPI         = TRUE
USE_CUDA        = FALSE
USE_OMP         = FALSE
USE_GPIBM       = FALSE
USE_EB          = FALSE
USE_PELEPHYSICS = FALSE

AMREX_NO_PROBINIT = TRUE

RZ_MMS_PHYSICS ?= navier-stokes
ifeq ($(RZ_MMS_PHYSICS),euler)
  MMS_VISCOUS := 0
  USERSuffix := $(USERSuffix).mms_rz_euler
else ifeq ($(RZ_MMS_PHYSICS),navier-stokes)
  MMS_VISCOUS := 1
  USERSuffix := $(USERSuffix).mms_rz_ns
else
  $(error RZ_MMS_PHYSICS must be one of: euler, navier-stokes)
endif

DEFINES += -DMMS_VISCOUS=$(MMS_VISCOUS)

# Central-difference order of the viscous face fluxes (2, 4, or 6).
RZ_MMS_VISCOUS_ORDER ?= 2
ifeq ($(filter $(RZ_MMS_VISCOUS_ORDER),2 4 6),)
  $(error RZ_MMS_VISCOUS_ORDER must be one of: 2, 4, 6)
endif
DEFINES += -DRZ_MMS_VISCOUS_ORDER=$(RZ_MMS_VISCOUS_ORDER)
ifneq ($(RZ_MMS_VISCOUS_ORDER),2)
  USERSuffix := $(USERSuffix).visc$(RZ_MMS_VISCOUS_ORDER)
endif

# Transport model: constant coefficients or nondimensional Sutherland laws.
RZ_MMS_TRANSPORT_MODEL ?= constant
ifeq ($(RZ_MMS_TRANSPORT_MODEL),constant)
  RZ_MMS_VARIABLE_TRANSPORT := 0
else ifeq ($(RZ_MMS_TRANSPORT_MODEL),sutherland)
  RZ_MMS_VARIABLE_TRANSPORT := 1
  USERSuffix := $(USERSuffix).sutherland
else
  $(error RZ_MMS_TRANSPORT_MODEL must be one of: constant, sutherland)
endif
DEFINES += -DRZ_MMS_VARIABLE_TRANSPORT=$(RZ_MMS_VARIABLE_TRANSPORT)

# WENO-Z nonlinear-weight exponent q (1 = production default).
CERISSE_WENOZ_POWER ?= 1
DEFINES += -DCERISSE_WENOZ_POWER=$(CERISSE_WENOZ_POWER)
ifneq ($(CERISSE_WENOZ_POWER),1)
  USERSuffix := $(USERSuffix).q$(CERISSE_WENOZ_POWER)
endif

RZ_MMS_EULER_SCHEME ?= llf-wenoz5
RZ_MMS_TENO5_STANDARD_TAU5 ?= 0
ifeq ($(RZ_MMS_EULER_SCHEME),llf-wenoz5)
  RZ_MMS_EULER_SCHEME_ID := 0
else ifeq ($(RZ_MMS_EULER_SCHEME),skew-central-o4)
  RZ_MMS_EULER_SCHEME_ID := 1
  USERSuffix := $(USERSuffix).skew_central_o4
else ifeq ($(RZ_MMS_EULER_SCHEME),skew-jst-o4)
  RZ_MMS_EULER_SCHEME_ID := 2
  USERSuffix := $(USERSuffix).skew_jst_o4
else ifeq ($(RZ_MMS_EULER_SCHEME),llf-teno5)
  RZ_MMS_EULER_SCHEME_ID := 3
  USERSuffix := $(USERSuffix).llf_teno5
  ifneq ($(strip $(RZ_MMS_TENO5_CUTOFF)),)
    RZ_MMS_TENO5_CUTOFF_TAG := $(subst .,_,$(subst -,_,$(RZ_MMS_TENO5_CUTOFF)))
    DEFINES += -DCERISSE_TENO5_SELECTION_CUTOFF=$(RZ_MMS_TENO5_CUTOFF)
    USERSuffix := $(USERSuffix).ct$(RZ_MMS_TENO5_CUTOFF_TAG)
  endif
  ifeq ($(RZ_MMS_TENO5_STANDARD_TAU5),1)
    DEFINES += -DCERISSE_TENO5_STANDARD_TAU5=1
    USERSuffix := $(USERSuffix).tau5
  else ifneq ($(RZ_MMS_TENO5_STANDARD_TAU5),0)
    $(error RZ_MMS_TENO5_STANDARD_TAU5 must be 0 or 1)
  endif
else ifeq ($(RZ_MMS_EULER_SCHEME),afd-hllc-wenoz5)
  RZ_MMS_EULER_SCHEME_ID := 4
  USERSuffix := $(USERSuffix).afd_hllc_wenoz5
else ifeq ($(RZ_MMS_EULER_SCHEME),afd-hllc-teno5)
  RZ_MMS_EULER_SCHEME_ID := 5
  USERSuffix := $(USERSuffix).afd_hllc_teno5
else ifeq ($(RZ_MMS_EULER_SCHEME),skew-central-o6)
  RZ_MMS_EULER_SCHEME_ID := 6
  USERSuffix := $(USERSuffix).skew_central_o6
else ifeq ($(RZ_MMS_EULER_SCHEME),llf-teno6)
  RZ_MMS_EULER_SCHEME_ID := 7
  USERSuffix := $(USERSuffix).llf_teno6
else
  $(error RZ_MMS_EULER_SCHEME must be one of: llf-wenoz5, llf-teno5, llf-teno6, afd-hllc-wenoz5, afd-hllc-teno5, skew-central-o4, skew-jst-o4, skew-central-o6)
endif

DEFINES += -DRZ_MMS_EULER_SCHEME_ID=$(RZ_MMS_EULER_SCHEME_ID)

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h

include $(AMR_SOLVER)/src/Make.CNS
