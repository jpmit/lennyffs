# mcfuncs.py
# 4th July 2012
# James Mithen

"""
Wrapper for importing correct fortran library
"""

import os

if os.uname()[0] == 'Linux':
    from mcfuncslinux import *
else:
    from mcfuncssun import *
