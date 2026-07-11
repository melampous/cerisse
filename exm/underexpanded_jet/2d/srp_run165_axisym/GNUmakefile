# 2D Axisymmetric SRP Run-165 — CPU build
DEBUG   = FALSE
DIM     = 2
COMP    = gcc
CXXSTD  = c++20
PRECISION = DOUBLE
TINY_PROFILE = TRUE
CLIP_MINTEMP = FALSE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE

USE_GPIBM = FALSE

AMREX_NO_PROBINIT = TRUE

AMR_SOLVER = /home/qiaoj/testcerisse/cerisse
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
