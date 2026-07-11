TINY_PROFILE = FALSE

DEBUG   = FALSE
DIM     = 1
COMP     = gcc
CXXSTD   = c++20

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE
USE_GPIBM=FALSE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

SHU_EULER_SCHEME ?= skew
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
else
  $(error SHU_EULER_SCHEME must be one of: skew, afd-hllc-teno5, afd-hllc-wenoz5, llf-teno5, llf-wenoz5)
endif
DEFINES += -DSHU_EULER_SCHEME_ID=$(SHU_EULER_SCHEME_ID)
USERSuffix := $(USERSuffix).shu_$(SHU_EULER_SCHEME)

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

include ./Make.package
include $(AMR_SOLVER)/src/Make.CNS
