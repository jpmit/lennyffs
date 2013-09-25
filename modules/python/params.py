# params.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Param class for handling the parameters in simulation input file (in)
nicely.
"""

import sys
from paramdata import *

class Param(object):
    def __init__(self, name, svalue):
        """Create Param object from name and value, both strings"""
        
        assert isinstance(name,str) and isinstance(svalue,str),(""
                "name and value of parameter need to be passed as string")
        self.name = name
        self.svalue = svalue
        self.value = self.getvalue()

    def getvalue(self):
        """Return svalue, which is a string, as proper type e.g. int"""
        
        if self.name not in PDICT:
            print "Error: parameter '%s' not allowed" %self.name
        else:
            try:
                value = FUNCS[PDICT[self.name]](self.svalue)
                return value
            except ValueError:
                print (("Error: parameter '%s' could not be converted "
                        "to type %s" %(self.name,PDICT[self.name])))
        return None

    def __str__(self):
        return ':'.join([self.name,self.svalue])

    def __repr__(self):
        return self.__str__()
