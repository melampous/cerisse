DEBUG             = FALSE
DIM               = 3
COMP              = gnu
CXXSTD            = c++20
PRECISION         = DOUBLE
TINY_PROFILE      = TRUE
CLIP_MINTEMP      = TRUE

USE_MPI           = TRUE
USE_CUDA          = TRUE
CUDA_ARCH         = 90
USE_OMP           = FALSE
USE_GPIBM         = TRUE
USE_EB            = FALSE
AMREX_NO_PROBINIT = TRUE

IBM_EXPERIMENTAL_CUT_CONTROL = FALSE

ifneq ($(USE_EB),FALSE)
  $(error Test1853 repair validation requires pure GP-IBM)
endif
ifneq ($(IBM_EXPERIMENTAL_CUT_CONTROL),FALSE)
  $(error Test1853 repair validation forbids cut-control flow)
endif

AMR_SOLVER = /shared/cerisse_local_0808_thermodynamic_prolongation
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
