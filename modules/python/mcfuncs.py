# mcfuncs.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper for importing correct fortran library
"""

import os

if os.uname()[0] == 'Linux':
    from mcfuncslinux import *
else:
    from mcfuncssun import *
