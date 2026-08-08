# Fast 2-D pure shared-GP/full-background-cell R-Z integration smoke.
DEBUG     = FALSE
DIM       = 2
COMP      = gcc
CXXSTD    = c++23
PRECISION = DOUBLE

TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE

USE_MPI   = FALSE
USE_CUDA  = FALSE
USE_OMP   = FALSE
USE_GPIBM = TRUE
USE_EB    = FALSE
IBM_EXPERIMENTAL_CUT_CONTROL = FALSE

ifneq ($(USE_EB),FALSE)
  $(error This test is pure GP-IBM; USE_EB must remain FALSE)
endif
ifneq ($(IBM_EXPERIMENTAL_CUT_CONTROL),FALSE)
  $(error This test forbids experimental cut-control flow updates)
endif

AMREX_NO_PROBINIT = TRUE
USERSuffix := $(USERSuffix).rz_pure_gp_positivity_smoke
ifeq ($(TEST_ACTIVE_RHS_HOOK),TRUE)
  DEFINES += -DCERISSE_TEST_ACTIVE_RHS_HOOK
  USERSuffix := $(USERSuffix).active_rhs_hook
endif
ifeq ($(TEST_RZ_GP_PARITY_PROBE),TRUE)
  ifneq ($(USE_CUDA),FALSE)
    $(error The direct pre/post GP parity probe is a CPU-only test fixture)
  endif
  DEFINES += -DCERISSE_TEST_RZ_GP_PARITY_PROBE
  USERSuffix := $(USERSuffix).generic_point_gp_axis_parity_probe
endif

AMR_SOLVER = $(abspath ../../../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
