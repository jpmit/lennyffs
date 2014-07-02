#  Makefile for Fortran module.  This uses f2py, which is part of
#  numpy, to create a shared library (.so) that can be imported into
#  Python.

SRC = global/initsimf.f90 ops/bopsf.f90 len/len_energy.f90 \
      len/len_mccyclenvt.f90 len/len_mccyclenpt.f90 \
      gauss/gauss_energy.f90 gauss/gauss_mccyclenvt.f90 \
      gauss/gauss_mccyclenpt.f90 gauss/gauss_force.f90 \
      gauss/gauss_mdcyclenve.f90 clist/clist.f90 \
      util/util.f90     

all:
	f2py -c -m --f90flags="-O3" mcfuncslinux $(SRC)
