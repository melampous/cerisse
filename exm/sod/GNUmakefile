TINY_PROFILE = FALSE

DEBUG   = FALSE
DIM     = 3
COMP    = gcc
CXXSTD  = c++20

USE_MPI  = TRUE
USE_CUDA = TRUE
USE_OMP  = FALSE
USE_GPIBM= FALSE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

SOD_EULER_SCHEME ?= riemann
ifeq ($(SOD_EULER_SCHEME),riemann)
  SOD_EULER_SCHEME_ID := 0
else ifeq ($(SOD_EULER_SCHEME),afd-hllc-teno5)
  SOD_EULER_SCHEME_ID := 1
else ifeq ($(SOD_EULER_SCHEME),afd-hllc-wenoz5)
  SOD_EULER_SCHEME_ID := 2
else ifeq ($(SOD_EULER_SCHEME),llf-teno5)
  SOD_EULER_SCHEME_ID := 3
else ifeq ($(SOD_EULER_SCHEME),llf-wenoz5)
  SOD_EULER_SCHEME_ID := 4
else ifeq ($(SOD_EULER_SCHEME),llf-teno6)
  SOD_EULER_SCHEME_ID := 5
else ifeq ($(SOD_EULER_SCHEME),llf-teno5-old)
  SOD_EULER_SCHEME_ID := 6
else ifeq ($(SOD_EULER_SCHEME),llf-wenoz5-old)
  SOD_EULER_SCHEME_ID := 7
else ifeq ($(SOD_EULER_SCHEME),llf-teno6-old)
  SOD_EULER_SCHEME_ID := 8
else
  $(error SOD_EULER_SCHEME must be one of: riemann, afd-hllc-teno5, afd-hllc-wenoz5, llf-teno5, llf-wenoz5, llf-teno6, llf-teno5-old, llf-wenoz5-old, llf-teno6-old)
endif
DEFINES += -DSOD_EULER_SCHEME_ID=$(SOD_EULER_SCHEME_ID)
USERSuffix := $(USERSuffix).sod_$(SOD_EULER_SCHEME)

AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
