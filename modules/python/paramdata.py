# paramdata.py
# 23 April 2013
# James Mithen

"""
Information on the allowed parameters,
and functions to convert them from a string to desired type
"""

INT = 'int'
INTLIST = 'intlist'
FLOAT = 'float'
STRING = 'str'
BOOL = 'bool'

# parameter dictionary contains type info for all allowed parameters
# key is name of parameter, value is type, either string, float, or int

PDICT = {# simulation params
         'simulation': STRING,
         'restartfile': STRING,
         # interaction potential params
         'potential' : STRING,
         'rcut': FLOAT,
         'r6mult': FLOAT,
         'r12mult': FLOAT,
         # fluid params
         'nstar': FLOAT,
         'Tstar': FLOAT,
         'nparfl': INT,
         'rcinit': FLOAT,
         'pressure': FLOAT,
         # surface params
         'surface': BOOL,
         'surftype': STRING,
         'lxsurf': INT,
         'lysurf': INT,
         'nlayersurf': INT,
         'nlatt': FLOAT,
         'flinit': STRING,
         'nlayerfl': INT,
         'fllayerspace': FLOAT,
         'plane': STRING,
         # box params
         'boxvol': FLOAT,
         'lboxx': FLOAT,
         'lboxy': FLOAT,
         'lboxz': FLOAT,
         # MC params
         'mctype': STRING,
         'ncycle': INT,
         'nsamp': INT,
         'maxdisp': FLOAT,
         'maxvol': FLOAT,
         # order params
         'stillsep': FLOAT,
         'q6link': FLOAT,
         'q6numlinks':INT,
         # FFS params
         'useffs': BOOL,
         'ffsname': STRING,
         'ffsrestart': INT,
         'lambdaA': INT,
         'numint': INT,
         'lambdas': INTLIST,
         'totalqhits': INT,
         'nshots': INT,
         'nbatch': INT,
         'minsuccess': INT,
         'lambdasamp': INT,
         'pruning': BOOL,
         'prunprob': FLOAT
         }

def __strtointlist(s):
    """Converts string e.g. '[0,2,3,4,6]' to a Python list of integers"""
    return [int(i) for i in s.split(',')]

def __strtobool(s):
    """Converts string to bool, note that we want to allow 'yes' and 'no'
    only, and throw an exception otherwise"""
    if s == 'yes':
        return True
    elif s == 'no':
        return False
    else:
        raise ValueError
        
# functions for converting from string to another type
# note that 'int', 'float', 'str' and 'bool' are built-in
    
FUNCS = {INT: int, FLOAT: float, STRING: str, BOOL: __strtobool,
         INTLIST: __strtointlist}
