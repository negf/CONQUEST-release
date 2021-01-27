#

# Set compilers
FC=mpiifort
F77=mpiifort

# Linking flags
LINKFLAGS=
ARFLAGS=

# Compilation flags
#CHECK_OPTS = -check all -check noarg_temp_created -O0
CHECK_OPTS = -xHost -O2
COMPFLAGS=$(CHECK_OPTS) -g -traceback -fpp $(XC_COMPFLAGS) -qopenmp
COMPFLAGS_F77= $(COMPFLAGS)

# Set BLAS and LAPACK libraries
BLAS=

# Full library call; remove scalapack if using dummy diag module
LIBS= $(FFT_LIB) $(XC_LIB) $(BLAS) -mkl=parallel  -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl

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
MULT_KERN = gemm
# Use dummy DiagModule or not
DIAG_DUMMY =

GRIDSOLVER = yes
ifeq ($(GRIDSOLVER), yes)
# DLmg gridsolver
  DLMG_DIR=$(HOME)/prog/dl_mg/git/dl_mg_code_public
  LIBS+= -L$(DLMG_DIR)/lib -ldlmg
  COMPFLAGS+= -I$(DLMG_DIR)/lib -DGRIDSOLVER
endif
