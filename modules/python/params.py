# params.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Param class for handling the parameters in simulation input file (in)
nicely.  The class will check if a value in the input file is (i)
allowed, (ii) of the correct type.  Examples of failing these checks
are given below for a couple of example lines from the input ('in')
file:
restrt new (parameter not allowed)
Tstar three (parameter cannot be converted to 'FLOAT', see paramdata.py).
"""

import sys
from paramdata import *

class Param(object):
    def __init__(self, name, svalue):
        """Create Param object from name and value, both strings"""
        
        assert isinstance(name, str) and isinstance(svalue, str),(""
                "name and value of parameter need to be passed as string")
        self.name = name
        self.svalue = svalue
        self.value = self.getvalue()

    def getvalue(self):
        """Return svalue, which is a string, as proper type e.g. int"""
        
        if self.name not in PDICT:
            print "Error: parameter '{0}' not allowed".format(self.name)
        else:
            try:
                value = FUNCS[PDICT[self.name]](self.svalue)
                return value
            except ValueError:
                print (("Error: parameter '{0}' could not be converted "
                        "to type {1}".format(self.name, PDICT[self.name])))
        return None

    def __str__(self):
        return ':'.join([self.name, self.svalue])

    def __repr__(self):
        return self.__str__()

def get_defaults():
    """Return dictionary of default parameters"""

    defaults = {}
    for (k, v) in DEFAULTS.values():
        p = Param(k, v)
        defaults[p.name] = p.value
    return defaults
