# mcfuncs.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper for importing fortran library and C++ library.  By importing
this module we get access to all of the C++ and Fortran functions by
mcfuncs.[function name].
"""

# Fortran 90 extension module (built with f2py)
from mcfuncslinux import *
# C++ extension module (built with Boost.Python)
from op_ext import *
