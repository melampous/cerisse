# 2D Axisymmetric jet (CPU build for local math verification)
DEBUG   = FALSE
DIM     = 2
COMP    = gcc
CXXSTD  = c++23
PRECISION = DOUBLE
TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE

USE_GPIBM = TRUE

AMREX_NO_PROBINIT = TRUE

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
