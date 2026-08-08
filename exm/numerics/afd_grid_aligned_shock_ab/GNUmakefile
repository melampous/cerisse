TINY_PROFILE = FALSE
DEBUG        = FALSE
PRECISION    = DOUBLE
DIM          = 2
COMP         = gcc
CXXSTD       = c++20

USE_MPI         = TRUE
USE_CUDA        = FALSE
USE_OMP         = FALSE
USE_GPIBM       = FALSE
USE_EB          = FALSE
USE_PELEPHYSICS = FALSE

AMREX_NO_PROBINIT = TRUE

# One executable is used for both members of the A/B comparison.  The only
# numerical-method change between the runs is the runtime option
# cns.afd_shock_llf.
USERSuffix := $(USERSuffix).afd_gridshock_wenoz5

AMR_SOLVER = $(abspath ../../..)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h

include $(AMR_SOLVER)/src/Make.CNS
