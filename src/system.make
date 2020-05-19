#

# Set compilers
FC=mpiifort
F77=mpiifort

# Linking flags
LINKFLAGS=
ARFLAGS=

# Compilation flags
COMPFLAGS= -O2 -g -traceback -fpp $(XC_COMPFLAGS)
COMPFLAGS_F77= $(COMPFLAGS)

# Set BLAS and LAPACK libraries
BLAS=

# Full library call; remove scalapack if using dummy diag module
LIBS= $(FFT_LIB) $(XC_LIB) $(BLAS) -mkl=cluster

# LibXC compatibility (LibXC below) or Conquest XC library

# Conquest XC library
XC_LIBRARY = CQ
XC_LIB =
XC_COMPFLAGS =

# LibXC compatibility
# Choose old LibXC (v2.x) or modern versions
#XC_LIBRARY = LibXC_v2
#XC_LIBRARY = LibXC
#XC_LIB = -lxcf90 -lxc
#XC_COMPFLAGS = -I/usr/local/include

# Set FFT library
FFT_LIB=
FFT_OBJ=fft_fftw3.o

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =

GRIDSOLVER = yes
ifeq ($(GRIDSOLVER), yes)
# DLmg gridsolver
  DLMG_DIR=$(HOME)/prog/dl_mg/2.0.2_intel_mpich
  LIBS+= -L$(DLMG_DIR)/lib -ldlmg
  COMPFLAGS+= -I$(DLMG_DIR)/lib -DGRIDSOLVER
endif
