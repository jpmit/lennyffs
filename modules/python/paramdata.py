# paramdata.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Information on the allowed input parameters, and functions to convert
them from a string to desired type.
"""

INT = 'int'
INTLIST = 'intlist'
FLOAT = 'float'
STRING = 'str'
BOOL = 'bool'

# parameter dictionary contains type info for all allowed parameters
# key is name of parameter, value is type, either string, float, int,
# or 'intlist'.

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
         'sameseed' : BOOL,

         # parameter for saving
         'nsave': INT,

         # order params
         'orderparam' : STRING,
         'opsamp' : INT,
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

def _strtointlist(s):
    """
    Converts string e.g. '[0,2,3,4,6]' to a Python list of integers.
    """
    try:
        ilist = [int(i) for i in s.split(',')]
    except:
        raise ValueError, 'could not convert "{0}" to intlist'.\
              format(s)
    
    return ilist

def _strtobool(s):
    """
    Converts string to bool, note that we want to allow 'yes' and 'no'
    only, and throw an exception otherwise.
    """
    
    if s == 'yes':
        return True
    elif s == 'no':
        return False
    else:
        raise ValueError, 'could not convert "{0}" to boolean'.\
              format(s)
        
# functions for converting from string to another type.  Note that
# 'int', 'float', 'str' and 'bool' are built-in.
    
FUNCS = {INT: int, FLOAT: float, STRING: str, BOOL: _strtobool,
         INTLIST: _strtointlist}
