# Reduced Cartesian underexpanded-jet A/B test.
#
# The problem definition is imported from the canonical NPR_0=3 viscous case.
# Only AFD-HLLC-WENO-Z5 is built because the study varies the runtime
# cns.afd_shock_llf switch with one executable.

DEBUG        = FALSE
DIM          = 3
COMP         = gcc
CXXSTD       = c++20
PRECISION    = DOUBLE
TINY_PROFILE = TRUE
CLIP_MINTEMP = FALSE

USE_MPI         = TRUE
USE_CUDA        = TRUE
USE_OMP         = FALSE
USE_GPIBM       = FALSE
USE_EB          = FALSE
USE_PELEPHYSICS = FALSE

AMREX_NO_PROBINIT = TRUE
CUDA_ARCH = 8.0 8.9

# JET_EULER_SCHEME_ID=2 selects afd_hllc_wenoz5_t in the canonical prob.h.
DEFINES += -DJET_EULER_SCHEME_ID=2
USERSuffix := $(USERSuffix).jet_afd-hllc-wenoz5_sensor-ab

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
