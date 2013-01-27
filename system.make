# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: system.make,v 1.2 2010/06/19 22:07:21 lt Exp $
#
# LCN Steinway system.make
#
#
FC = ftn
F77 = ftn
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS = -L/homes/lt/usr/steinway/lib -L/homes/lt/usr/steinway/acml5.2.0/gfortran64/lib -L/homes/lt/usr/steinway/gcc-4.6.3/lib64
# COMPFLAGS = -g
COMPFLAGS = -O3
# COMPFLAGS_F77 = -g
COMPFLAGS_F77 = -O3
ARFLAGS = 

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
BLACS = 
LIBS = $(FFT) -lacml -lgomp -lscalapack
# This line is for systems with NO ScaLAPACK and BLAS - replace DiagModule.f90 with DiagModule.f90_DUMMY
#LIBS = $(FFT) -llapack -lblas $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC77=$(F77)" "FC90=$(FC)" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

# Matrix multiplication kernel type
MULT_KERN = ompDojk
# Use dummy DiagModule or not
DIAG_DUMMY =
