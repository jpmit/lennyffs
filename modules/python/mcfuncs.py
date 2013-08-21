# mcfuncs.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper for importing correct fortran library
"""

import os

if os.uname()[0] == 'Linux':
    # Fortran 90 extension module (built with f2py)
    from mcfuncslinux import *
    # C++ extension module (built with Boost.Python)
    from op_ext import *
else:
    # this is now deprecated
    from mcfuncssun import *
