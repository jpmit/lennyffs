# mcfuncs.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper for importing fortran library and C++ library.  By importing
this module we get access to all of the C++ and Fortran functions by
mcfuncs.[function name].  Note that these functions are wrapped for
convenience by the functions in modules/python/mccycle.py,
modules/python/energy.py and modules/python/mccycle.py.  They should
be called via these interfaces, and not by mcfuncs.[function name].
"""

# Fortran 90 extension module (built with f2py)
from mcfuncslinux import *

# C++ extension module (built with Boost.Python)
try:
    from op_ext import *
except ImportError:
    # print hostname: this is for debugging when the shared library is
    # not on the machine.
    import socket
    print socket.gethostname()
    # re-raise the previous error so that we still get the normal
    # exception information.
    raise
