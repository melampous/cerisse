# 2D Axisymmetric jet (CPU build for local math verification)
DEBUG   = FALSE
DIM     = 2
COMP    = gcc
CXXSTD  = c++23
PRECISION = DOUBLE
TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE

USE_GPIBM = TRUE
USE_EB = FALSE
IBM_EXPERIMENTAL_CUT_CONTROL = FALSE

ifneq ($(USE_EB),FALSE)
  $(error NPR-3 RZ is a pure GP-IBM case; USE_EB must be FALSE)
endif
ifneq ($(IBM_EXPERIMENTAL_CUT_CONTROL),FALSE)
  $(error NPR-3 RZ is a pure GP-IBM case; experimental cut-control must be FALSE)
endif

# Keep WENO-Z5 as the established default.  Give TENO5 a distinct target and
# object directory so the two template instantiations cannot share objects.
NPR3_EULER_SCHEME ?= wenoz5
ifeq ($(NPR3_EULER_SCHEME),wenoz5)
  NPR3_EULER_SCHEME_ID := 1
else ifeq ($(NPR3_EULER_SCHEME),teno5)
  NPR3_EULER_SCHEME_ID := 2
  USERSuffix := $(USERSuffix).npr3_rz_pure_gp_llf_teno5
else ifeq ($(NPR3_EULER_SCHEME),afd-hllc-wenoz5)
  NPR3_EULER_SCHEME_ID := 3
  USERSuffix := $(USERSuffix).npr3_rz_pure_gp_afd_hllc_wenoz5
else ifeq ($(NPR3_EULER_SCHEME),afd-hllc-teno5)
  NPR3_EULER_SCHEME_ID := 4
  USERSuffix := $(USERSuffix).npr3_rz_pure_gp_afd_hllc_teno5
else
  $(error NPR3_EULER_SCHEME must be one of: wenoz5, teno5, afd-hllc-wenoz5, afd-hllc-teno5)
endif
DEFINES += -DNPR3_EULER_SCHEME_ID=$(NPR3_EULER_SCHEME_ID)

AMREX_NO_PROBINIT = TRUE

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
