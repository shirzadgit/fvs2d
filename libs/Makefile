#
#  Top Level Makefile for LAPACK
#  Version 3.4.1
#  April 2012
#

include make.inc

lib: blaslib

clean: cleanblas


blaslib:
	( cd BLAS; $(MAKE) )

cleanblas:	
	( cd BLAS; $(MAKE) clean )

