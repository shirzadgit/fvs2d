#!/bin/bash
echo ' '
echo -n 'input compiler name (intel, pgi, cray): '; read F; FM="${F^^}"; 
echo -n 'input mpi-wrapper (mpiifort, mpif90, ftn): '; read F; FC="${F,,}"; 
echo ' '
echo ' '
sed -i "/FC=/c\FC=${FC}" "make.inc" 
sed -i "/FM=/c\FM=${FM}" "make.inc"

