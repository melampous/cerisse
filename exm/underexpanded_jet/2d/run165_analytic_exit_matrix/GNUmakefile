# Run165 RZ analytic-exit root-cause matrix (CPU, no IBM).
DEBUG          = FALSE
DIM            = 2
COMP           = gcc
CXXSTD         = c++20
PRECISION      = DOUBLE
TINY_PROFILE   = TRUE
CLIP_MINTEMP   = FALSE

USE_MPI   = TRUE
USE_CUDA  = FALSE
USE_OMP   = FALSE
USE_GPIBM = FALSE

AMREX_NO_PROBINIT = TRUE

JET_ENABLED   ?= 1
JET_SCALE     ?= 1.0
NS_ENABLED    ?= 0
LES_ENABLED   ?= 0
LES_CS        ?= 0.16
WENOZ_POWER   ?= 1
FACE_TRACE    ?= 0
EULER_SCHEME  ?= weno
SKEW_C2       ?= 1.5
SKEW_C4       ?= 0.016
ifeq ($(EULER_SCHEME),weno)
  EULER_SCHEME_ID := 0
else ifeq ($(EULER_SCHEME),skew)
  EULER_SCHEME_ID := 1
else ifeq ($(EULER_SCHEME),hllc-muscl)
  EULER_SCHEME_ID := 2
else ifeq ($(EULER_SCHEME),afd-hllc-weno)
  EULER_SCHEME_ID := 3
else ifeq ($(EULER_SCHEME),llf-teno)
  EULER_SCHEME_ID := 4
else ifeq ($(EULER_SCHEME),afd-hllc-teno)
  EULER_SCHEME_ID := 5
else
  $(error EULER_SCHEME must be weno, skew, hllc-muscl, afd-hllc-weno, llf-teno, or afd-hllc-teno)
endif
JET_SCALE_TAG  = $(subst .,p,$(JET_SCALE))
SKEW_C2_TAG    = $(subst .,_,$(SKEW_C2))
SKEW_C4_TAG    = $(subst .,_,$(SKEW_C4))
EBASE         ?= run165_analytic_exit_$(EULER_SCHEME)_jet$(JET_ENABLED)_js$(JET_SCALE_TAG)_cpu
TMP_BUILD_DIR ?= tmp_build_run165_analytic_exit_$(EULER_SCHEME)_jet$(JET_ENABLED)_js$(JET_SCALE_TAG)_c2$(SKEW_C2_TAG)_c4$(SKEW_C4_TAG)

DEFINES += -DRUN165_JET_ENABLED=$(JET_ENABLED)
DEFINES += -DRUN165_JET_SCALE=$(JET_SCALE)
DEFINES += -DRUN165_USE_NS=$(NS_ENABLED)
DEFINES += -DRUN165_USE_LES=$(LES_ENABLED)
DEFINES += -DRUN165_LES_CS=$(LES_CS)
DEFINES += -DCERISSE_WENOZ_POWER=$(WENOZ_POWER)
DEFINES += -DRUN165_EULER_SCHEME=$(EULER_SCHEME_ID)
DEFINES += -DRUN165_SKEW_C2=$(SKEW_C2)
DEFINES += -DRUN165_SKEW_C4=$(SKEW_C4)
ifeq ($(FACE_TRACE),1)
  DEFINES += -DCERISSE_ENABLE_FACE_TRACE=1
endif

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
