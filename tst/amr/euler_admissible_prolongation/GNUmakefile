AMREX_HOME := $(abspath ../../../lib/amrex)
CERISSE_HOME := $(abspath ../../..)

DEBUG := FALSE
DIM ?= 3
COMP := gnu
USE_MPI := FALSE
USE_OMP := FALSE
BL_NO_FORT := TRUE

CEXE_sources += main.cpp
INCLUDE_LOCATIONS += $(CERISSE_HOME)/src/set

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules

warning_flags += -Werror=return-type
