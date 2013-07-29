SRC = global/initsimf.f90 ops/bopsf.f90 len/len_energy.f90 len/len_mccyclenvt.f90 \
      len/len_mccyclenpt.f90

all:
	f2py -c -m --f90flags="-O3" mcfuncslinux $(SRC)
