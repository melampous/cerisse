TINY_PROFILE = FALSE
DEBUG        = FALSE
PRECISION    = DOUBLE
DIM          = 2
COMP         = gcc
CXXSTD       = c++20

USE_MPI       = TRUE
USE_CUDA      = FALSE
USE_OMP       = FALSE
USE_GPIBM     = FALSE
USE_EB        = FALSE
USE_PELEPHYSICS = FALSE

AMREX_NO_PROBINIT = TRUE

# Thesis baseline: local Lax--Friedrichs flux splitting with characteristic
# WENO-Z5 reconstruction.  The alternatives are separate executables so a
# method comparison cannot silently replace the baseline.
IV_EULER_SCHEME ?= llf-wenoz5
ifeq ($(IV_EULER_SCHEME),llf-wenoz5)
  IV_EULER_SCHEME_ID := 0
else ifeq ($(IV_EULER_SCHEME),llf-teno5)
  IV_EULER_SCHEME_ID := 1
else ifeq ($(IV_EULER_SCHEME),afd-hllc-wenoz5)
  IV_EULER_SCHEME_ID := 2
else ifeq ($(IV_EULER_SCHEME),afd-hllc-teno5)
  IV_EULER_SCHEME_ID := 3
else ifeq ($(IV_EULER_SCHEME),skew-central-o4)
  IV_EULER_SCHEME_ID := 4
else ifeq ($(IV_EULER_SCHEME),skew-jst-o4)
  IV_EULER_SCHEME_ID := 5
else
  $(error IV_EULER_SCHEME must be one of: llf-wenoz5, llf-teno5, afd-hllc-wenoz5, afd-hllc-teno5, skew-central-o4, skew-jst-o4)
endif
DEFINES += -DIV_EULER_SCHEME_ID=$(IV_EULER_SCHEME_ID)
USERSuffix := $(USERSuffix).iv_$(IV_EULER_SCHEME)

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h

include $(AMR_SOLVER)/src/Make.CNS
