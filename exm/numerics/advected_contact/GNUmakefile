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

# Both values of cns.afd_shock_llf are run with this same executable.
USERSuffix := $(USERSuffix).contact_afd-hllc-wenoz5

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h

include $(AMR_SOLVER)/src/Make.CNS
