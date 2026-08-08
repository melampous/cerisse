DEBUG   = FALSE
DIM     = 2
COMP    = gcc
CXXSTD  = c++20
PRECISION = DOUBLE
TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE
USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE
USE_GPIBM = FALSE
AMREX_NO_PROBINIT = TRUE
# CPU case. Local test build: make AMR_SOLVER=/home/qiaoj/testcerisse/cerisse
AMR_SOLVER ?= /shared/cerisse_reflux
AMREX_HOME = $(AMR_SOLVER)/lib/amrex
CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
