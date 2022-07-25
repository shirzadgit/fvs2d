#!/bin/bash
echo ' '
echo -n ' input mpi-wrapper (e.g: mpiifort, mpif90, ftn): '; read F; FC="${F,,}";

echo ' '
echo " compiler: ${FC}"
echo " executable: fvs2d_nlin.exe"
echo " step 1 >>> create build directory"
echo " step 2 >>> cd build"
echo " step 3 >>> cmake -DCMAKE_Fortran_COMPILER=${FC} .."
echo ' step 4 >>> make'
echo ' executable file will be in build/run'
echo ' '
sleep 2

mkdir -p build
cd "build"
rm -rf *
touch cmake.sh
chmod +x cmake.sh
echo "cmake -DCMAKE_Fortran_COMPILER=${FC} .." >> cmake.sh
./cmake.sh
make
