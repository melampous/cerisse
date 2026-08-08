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

MMS_PHYSICS ?= euler
ifeq ($(MMS_PHYSICS),euler)
  MMS_VISCOUS := 0
else ifeq ($(MMS_PHYSICS),navier-stokes)
  MMS_VISCOUS := 1
else
  $(error MMS_PHYSICS must be one of: euler, navier-stokes)
endif

MMS_VISCOUS_ORDER ?= 2
ifeq ($(filter $(MMS_VISCOUS_ORDER),2 4 6),)
  $(error MMS_VISCOUS_ORDER must be one of: 2, 4, 6)
endif

MMS_TRANSPORT_MODEL ?= constant
ifeq ($(MMS_TRANSPORT_MODEL),constant)
  MMS_VARIABLE_TRANSPORT := 0
  MMS_TRANSPORT_SUFFIX :=
else ifeq ($(MMS_TRANSPORT_MODEL),sutherland)
  MMS_VARIABLE_TRANSPORT := 1
  MMS_TRANSPORT_SUFFIX := .sutherland
else
  $(error MMS_TRANSPORT_MODEL must be one of: constant, sutherland)
endif

# Skew-symmetric central scheme options (used when MMS_EULER_SCHEME=skew).
MMS_SKEW_ORDER ?= 4
MMS_SKEW_DISSIPATION ?= 0
# Damping difference order (0 = historical default order-1; 3 or 5) and
# discontinuity-sensor exponent (1 = original JST scaling).
MMS_SKEW_DAMP_ORDER ?= 0
MMS_SKEW_SENSOR_POWER ?= 1

MMS_EULER_SCHEME ?= llf-wenoz5
MMS_TENO5_STANDARD_TAU5 ?= 0
ifeq ($(MMS_EULER_SCHEME),llf-wenoz5)
  MMS_EULER_SCHEME_ID := 0
  MMS_SCHEME_SUFFIX :=
else ifeq ($(MMS_EULER_SCHEME),llf-teno5)
  MMS_EULER_SCHEME_ID := 1
  MMS_SCHEME_SUFFIX := .llf_teno5
  ifneq ($(strip $(MMS_TENO5_CUTOFF)),)
    MMS_TENO5_CUTOFF_TAG := $(subst .,_,$(subst -,_,$(MMS_TENO5_CUTOFF)))
    DEFINES += -DCERISSE_TENO5_SELECTION_CUTOFF=$(MMS_TENO5_CUTOFF)
    MMS_SCHEME_SUFFIX := $(MMS_SCHEME_SUFFIX).ct$(MMS_TENO5_CUTOFF_TAG)
  endif
  ifeq ($(MMS_TENO5_STANDARD_TAU5),1)
    DEFINES += -DCERISSE_TENO5_STANDARD_TAU5=1
    MMS_SCHEME_SUFFIX := $(MMS_SCHEME_SUFFIX).tau5
  else ifneq ($(MMS_TENO5_STANDARD_TAU5),0)
    $(error MMS_TENO5_STANDARD_TAU5 must be 0 or 1)
  endif
else ifeq ($(MMS_EULER_SCHEME),afd-hllc-wenoz5)
  MMS_EULER_SCHEME_ID := 2
  MMS_SCHEME_SUFFIX := .afd_hllc_wenoz5
else ifeq ($(MMS_EULER_SCHEME),afd-hllc-teno5)
  MMS_EULER_SCHEME_ID := 3
  MMS_SCHEME_SUFFIX := .afd_hllc_teno5
else ifeq ($(MMS_EULER_SCHEME),llf-teno6)
  MMS_EULER_SCHEME_ID := 4
  MMS_SCHEME_SUFFIX := .llf_teno6
else ifeq ($(MMS_EULER_SCHEME),llf-wenoz5-old)
  MMS_EULER_SCHEME_ID := 5
  MMS_SCHEME_SUFFIX := .llf_wenoz5_old
else ifeq ($(MMS_EULER_SCHEME),llf-teno5-old)
  MMS_EULER_SCHEME_ID := 6
  MMS_SCHEME_SUFFIX := .llf_teno5_old
else ifeq ($(MMS_EULER_SCHEME),llf-teno6-old)
  MMS_EULER_SCHEME_ID := 7
  MMS_SCHEME_SUFFIX := .llf_teno6_old
else ifeq ($(MMS_EULER_SCHEME),skew)
  MMS_EULER_SCHEME_ID := 8
  ifeq ($(MMS_SKEW_DISSIPATION),1)
    MMS_SCHEME_SUFFIX := .skew$(MMS_SKEW_ORDER)_jst
  else
    MMS_SCHEME_SUFFIX := .skew$(MMS_SKEW_ORDER)_nodiss
  endif
  ifneq ($(MMS_SKEW_DAMP_ORDER),0)
    MMS_SCHEME_SUFFIX := $(MMS_SCHEME_SUFFIX)_d$(MMS_SKEW_DAMP_ORDER)
  endif
  ifneq ($(MMS_SKEW_SENSOR_POWER),1)
    MMS_SCHEME_SUFFIX := $(MMS_SCHEME_SUFFIX)_sp$(MMS_SKEW_SENSOR_POWER)
  endif
else
  $(error MMS_EULER_SCHEME must be one of: llf-wenoz5, llf-teno5, llf-teno6, afd-hllc-wenoz5, afd-hllc-teno5, skew, llf-wenoz5-old, llf-teno5-old, llf-teno6-old)
endif

ifeq ($(filter $(MMS_SKEW_ORDER),2 4 6),)
  $(error MMS_SKEW_ORDER must be one of: 2, 4, 6)
endif
DEFINES += -DMMS_SKEW_ORDER=$(MMS_SKEW_ORDER)
DEFINES += -DMMS_SKEW_DISSIPATION=$(MMS_SKEW_DISSIPATION)
DEFINES += -DMMS_SKEW_DAMP_ORDER=$(MMS_SKEW_DAMP_ORDER)
DEFINES += -DMMS_SKEW_SENSOR_POWER=$(MMS_SKEW_SENSOR_POWER)

DEFINES += -DMMS_VISCOUS=$(MMS_VISCOUS)
DEFINES += -DMMS_VISCOUS_ORDER=$(MMS_VISCOUS_ORDER)
DEFINES += -DMMS_EULER_SCHEME_ID=$(MMS_EULER_SCHEME_ID)
DEFINES += -DMMS_VARIABLE_TRANSPORT=$(MMS_VARIABLE_TRANSPORT)

# WENO-Z exponent q (1 = production default, 2 = squared-ratio variant).
CERISSE_WENOZ_POWER ?= 1
DEFINES += -DCERISSE_WENOZ_POWER=$(CERISSE_WENOZ_POWER)
MMS_VISCOUS_SUFFIX :=
ifneq ($(MMS_VISCOUS_ORDER),2)
  MMS_VISCOUS_SUFFIX := .visc$(MMS_VISCOUS_ORDER)
endif
USERSuffix := $(USERSuffix).mms_$(subst -,_,$(MMS_PHYSICS))$(MMS_SCHEME_SUFFIX)$(MMS_VISCOUS_SUFFIX)$(MMS_TRANSPORT_SUFFIX)
ifneq ($(CERISSE_WENOZ_POWER),1)
  USERSuffix := $(USERSuffix).q$(CERISSE_WENOZ_POWER)
endif

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h

include $(AMR_SOLVER)/src/Make.CNS
