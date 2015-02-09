# Build Fortran shared library
cd ./modules/fortran
make
# Build C++ shared library
cd ../../modules/cpp
make
cd ../..
# [base] directory
BDIR=$(pwd)
echo PYTHONPATH=${BDIR}/modules/python:${BDIR}/modules/fortran:${BDIR}/modules/cpp:\$PYTHONPATH >> ~/.bashrc
echo PATH=${BDIR}/exec:${BDIR}/util:\$PATH >> ~/.bashrc
