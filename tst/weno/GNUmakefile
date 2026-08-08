DIM        := 3
COMP       := gcc
CXXSTD     := c++20
DEBUG      := FALSE
USE_MPI    := FALSE
USE_OMP    := FALSE
USE_CUDA   := FALSE
USE_EB     := FALSE
BL_NO_FORT := TRUE

AMR_SOLVER := $(abspath ../..)
AMREX_HOME := $(AMR_SOLVER)/lib/amrex

CEXE_sources += main.cpp cuda_compile.cpp
CEXE_headers += Weno.h Weno_old.h

INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/rhs
INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/cls
INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/set

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
