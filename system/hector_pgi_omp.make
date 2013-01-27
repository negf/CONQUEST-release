# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: system.make,v 1.2 2010/06/19 22:07:21 lt Exp $
#
# Conquest example system.make
#
#
FC = ftn
F77 = ftn
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS =
# COMPFLAGS = -g
COMPFLAGS = -fastsse -mp
# COMPFLAGS_F77 = -g
COMPFLAGS_F77 = -fastsse
ARFLAGS =

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
BLACS =
LIBS = $(FFT) -lsci_pgi
# This line is for systems with NO ScaLAPACK and BLAS - replace DiagModule.f90 with DiagModule.f90_DUMMY
#LIBS = $(FFT) -llapack -lblas $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC77=$(F77)" "FC90=$(FC)" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

# Matrix multiplication kernel type
MULT_KERN = omp
# Use dummy DiagModule or not
DIAG_DUMMY =
