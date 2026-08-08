TINY_PROFILE = FALSE
DEBUG        = FALSE
PRECISION    = DOUBLE
DIM          = 1
COMP         = gcc
CXXSTD       = c++20

USE_MPI         = TRUE
USE_CUDA        = FALSE
USE_OMP         = FALSE
USE_GPIBM       = FALSE
USE_EB          = FALSE
USE_PELEPHYSICS = FALSE

AMREX_NO_PROBINIT = TRUE

# Each method is compiled into a separate executable.  This prevents a
# runtime option from silently changing the spatial discretisation.
SHEAR_EULER_SCHEME ?= llf-wenoz5
ifeq ($(SHEAR_EULER_SCHEME),llf-wenoz5)
  SHEAR_EULER_SCHEME_ID := 0
else ifeq ($(SHEAR_EULER_SCHEME),llf-teno5)
  SHEAR_EULER_SCHEME_ID := 1
else ifeq ($(SHEAR_EULER_SCHEME),afd-hllc-wenoz5)
  SHEAR_EULER_SCHEME_ID := 2
else ifeq ($(SHEAR_EULER_SCHEME),skew-jst-o4)
  SHEAR_EULER_SCHEME_ID := 3
else
  $(error SHEAR_EULER_SCHEME must be one of: llf-wenoz5, llf-teno5, afd-hllc-wenoz5, skew-jst-o4)
endif

DEFINES += -DSHEAR_EULER_SCHEME_ID=$(SHEAR_EULER_SCHEME_ID)
USERSuffix := $(USERSuffix).shear_$(SHEAR_EULER_SCHEME)

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h

include $(AMR_SOLVER)/src/Make.CNS
