SHELL = /bin/bash

#------------------------------------------------------
# Select Compiler
# FM:{INTEL, PGI, CRAY}
# FC:{mpiifort,mpif90,ftn}
#------------------------------------------------------
FM=INTEL
FC=mpiifort

#------------------------------------------------------
# Choose production or debug Code
# PRD: Production Code
# DBG: Debug Code
#------------------------------------------------------
DEBUG = PRD

#------------------------------------------------------
# Compiler and Options:
# FFLAGC: for code
# FFDEBG: for debugging
# FFLAGL: for library
#------------------------------------------------------
ifneq (,$(findstring INTEL,$(FM))) #-- INTEL COMPILERS
  FFLAGC = -r8 -Ofast -fpp -convert big_endian -traceback -qopenmp -mcmodel=medium
  FFLAGL = -r8 -Ofast -fPIC 
  FFDEBG = -g -shared-intel -debug inline-debug-info -check bounds
  FCOMPL = ${FC}
else ifneq (,$(findstring CRAY,$(FM))) #-- CRAY COMPILERS
  FFLAGC = -s real64 -Ofast -e P -h byteswapio
  FFLAGL = -s real64 -Ofast -fPIC 
  FFDEBG = -g 
  FCOMPL = ${FC}
else ifneq (,$(findstring PGI,$(FM))) #-- PGI GROUP COMPILERS
  FFLAGC = -r8 -Ofast -Mpreprocess -byteswapio -traceback -mp -mcmodel=medium
  FFLAGL = -r8 -Ofast -fPIC
  FFDEBG = -g
  FCOMPL = ${FC}
else ifneq (,$(findstring GNU,$(FM))) #-- PGI GROUP COMPILERS
  FFLAGC = -fdefault-real-8 -fdefault-double-8 -Ofast -cpp -fconvert=big-endian -fbacktrace -fopenmp -ffree-line-length-none --std=legacy
  FFLAGL = -fdefault-real-8 -fdefault-double-8 -Ofast -fPIC
  FFDEBG = -g
  FCOMPL = ${FC}
else
  FFLAGC = -r8 -Ofast -fpp -g -shared-intel -debug inline-debug-info -convert big_endian -traceback -mcmodel=large
  FFLAGL = -r8 -Ofast -fPIC -traceback
  FFDEBG = -g -shared-intel -debug inline-debug-info -check bounds
  FCOMPL = mpiifort
endif

#------------------------------------------------------
# Library
#------------------------------------------------------
LIBS    = -L ./../libs -lrefblas 

#------------------------------------------------------
# Executable directory
#------------------------------------------------------
RUN = ../run

#------------------------------------------------------
# Host
#------------------------------------------------------
host = $(shell hostname)


SOLVER_TYPE=NLIN
