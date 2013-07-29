# lenexceptions.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Exception classes used globally.
"""

# error with a utility (see util/ folder)
class UtilError(Exception): pass

# input parameter missing from file
class NoInputParamError(Exception): pass

# input parameter is not on an allowed list
class InputParamNotAllowedError(Exception): pass

# input parameter can't be converted properly
class InputParamUnconvertible(Exception): pass
