# Current all-fluid R-Z pressure-axis regression.  The scheme tag is part of
# both the object directory and the executable name, so an incremental build
# cannot silently reuse objects (or an executable) from a different spatial
# operator.
DEBUG     = FALSE
DIM       = 2
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE

USE_MPI   = FALSE
USE_CUDA  = FALSE
USE_OMP   = FALSE
USE_GPIBM = FALSE

AMREX_NO_PROBINIT = TRUE

RZ_AXIS_SCHEME ?= llf-wenoz5
ifeq ($(RZ_AXIS_SCHEME),llf-wenoz5)
  RZ_AXIS_SCHEME_ID := 0
  RZ_AXIS_SCHEME_TAG := llf_wenoz5
else ifeq ($(RZ_AXIS_SCHEME),llf-teno5)
  RZ_AXIS_SCHEME_ID := 1
  RZ_AXIS_SCHEME_TAG := llf_teno5
else ifeq ($(RZ_AXIS_SCHEME),llf-teno6)
  RZ_AXIS_SCHEME_ID := 2
  RZ_AXIS_SCHEME_TAG := llf_teno6
else ifeq ($(RZ_AXIS_SCHEME),afd-hllc-wenoz5)
  RZ_AXIS_SCHEME_ID := 3
  RZ_AXIS_SCHEME_TAG := afd_hllc_wenoz5
else ifeq ($(RZ_AXIS_SCHEME),afd-hllc-teno5)
  RZ_AXIS_SCHEME_ID := 4
  RZ_AXIS_SCHEME_TAG := afd_hllc_teno5
else ifeq ($(RZ_AXIS_SCHEME),skew4-central)
  RZ_AXIS_SCHEME_ID := 5
  RZ_AXIS_SCHEME_TAG := skew4_central
else ifeq ($(RZ_AXIS_SCHEME),skew4-jst)
  RZ_AXIS_SCHEME_ID := 6
  RZ_AXIS_SCHEME_TAG := skew4_jst
else
  $(error RZ_AXIS_SCHEME must be one of: llf-wenoz5, llf-teno5, llf-teno6, afd-hllc-wenoz5, afd-hllc-teno5, skew4-central, skew4-jst)
endif

DEFINES += -DRZ_AXIS_SCHEME_ID=$(RZ_AXIS_SCHEME_ID)
USERSuffix := $(USERSuffix).$(RZ_AXIS_SCHEME_TAG)

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
