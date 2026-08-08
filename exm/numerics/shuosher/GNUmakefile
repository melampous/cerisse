TINY_PROFILE = FALSE

DEBUG   = FALSE
PRECISION = DOUBLE
DIM     = 1
COMP     = gcc
CXXSTD   = c++20

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE
USE_GPIBM=FALSE
USE_EB = FALSE
USE_PELEPHYSICS = FALSE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

# Thesis baseline: local Lax--Friedrichs flux splitting with characteristic
# WENO-Z5 reconstruction.  The alternative schemes remain available through
# `make SHU_EULER_SCHEME=<name>` for method comparisons.
SHU_EULER_SCHEME ?= llf-wenoz5
SHU_SKEW_C2 ?= 1.5
SHU_SKEW_C4 ?= 0.016
SHU_TENO5_CUTOFF ?=
SHU_BUILD_TAG ?=
ifeq ($(SHU_EULER_SCHEME),skew)
  SHU_EULER_SCHEME_ID := 0
else ifeq ($(SHU_EULER_SCHEME),afd-hllc-teno5)
  SHU_EULER_SCHEME_ID := 1
else ifeq ($(SHU_EULER_SCHEME),afd-hllc-wenoz5)
  SHU_EULER_SCHEME_ID := 2
else ifeq ($(SHU_EULER_SCHEME),llf-teno5)
  SHU_EULER_SCHEME_ID := 3
else ifeq ($(SHU_EULER_SCHEME),llf-wenoz5)
  SHU_EULER_SCHEME_ID := 4
else ifeq ($(SHU_EULER_SCHEME),llf-teno6)
  SHU_EULER_SCHEME_ID := 5
else ifeq ($(SHU_EULER_SCHEME),llf-teno5-old)
  SHU_EULER_SCHEME_ID := 6
else ifeq ($(SHU_EULER_SCHEME),llf-wenoz5-old)
  SHU_EULER_SCHEME_ID := 7
else ifeq ($(SHU_EULER_SCHEME),llf-teno6-old)
  SHU_EULER_SCHEME_ID := 8
else
  $(error SHU_EULER_SCHEME must be one of: skew, afd-hllc-teno5, afd-hllc-wenoz5, llf-teno5, llf-wenoz5, llf-teno6, llf-teno5-old, llf-wenoz5-old, llf-teno6-old)
endif
DEFINES += -DSHU_EULER_SCHEME_ID=$(SHU_EULER_SCHEME_ID)
DEFINES += -DSHU_SKEW_C2_VALUE=$(SHU_SKEW_C2)
DEFINES += -DSHU_SKEW_C4_VALUE=$(SHU_SKEW_C4)

# Include the Skew--JST coefficients in that executable's name.  This prevents
# Make from reusing a binary built with a different coefficient pair.
ifeq ($(SHU_EULER_SCHEME),skew)
  SHU_SKEW_C2_TAG := $(subst .,_,$(SHU_SKEW_C2))
  SHU_SKEW_C4_TAG := $(subst .,_,$(SHU_SKEW_C4))
  USERSuffix := $(USERSuffix).shu_skew_c2_$(SHU_SKEW_C2_TAG)_c4_$(SHU_SKEW_C4_TAG)
else
  USERSuffix := $(USERSuffix).shu_$(SHU_EULER_SCHEME)
endif
ifeq ($(SHU_EULER_SCHEME),llf-teno5)
  ifneq ($(strip $(SHU_TENO5_CUTOFF)),)
    SHU_TENO5_CUTOFF_TAG := $(subst .,_,$(subst -,_,$(SHU_TENO5_CUTOFF)))
    DEFINES += -DCERISSE_TENO5_SELECTION_CUTOFF=$(SHU_TENO5_CUTOFF)
    USERSuffix := $(USERSuffix).ct$(SHU_TENO5_CUTOFF_TAG)
  endif
endif
ifneq ($(strip $(SHU_BUILD_TAG)),)
  USERSuffix := $(USERSuffix)_$(SHU_BUILD_TAG)
endif

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

CEXE_headers += prob.h

include ./Make.package
include $(AMR_SOLVER)/src/Make.CNS
