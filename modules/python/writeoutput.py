# writeoutput.py
# 26th June 2012
# James Mithen

"""
Functions for writing output of FFS code
FUNCTIONS:
writepickparams - write parameters dictionary to pickle file
writeparams - write dictionary of parameters in human readable
            - form (i.e. not pickled)
writexyz - convenience function for writing xyz file
"""

import os
import bops
import readwrite
import pickle

def writepickparams(params, fname='params.pkl'):
    """Write out parameters dictionary as .pkl file"""
    fout = open(fname,'wb')
    pickle.dump(params,fout)
    fout.close()
    return

def writeparams(params, fname='params.out'):
    """Write out dictionary of parameters in human readable form"""
    fout = open(fname,'w')
    outstr = ('# Simulation parameters (Human readable version of params.pkl)\n'
             '# Note that only some of these parameters are used\n'
             '# Others are redundant and could prove downright misleading!\n')
    keys = params.keys()
    keys.sort()
    for k in keys:
        outstr = '%s%s %s\n' %(outstr,k,str(params[k]))
    fout.write(outstr)
    fout.close()
    return

def writexyz(fname,positions,params):
    """Write positions in xyz format with symbols as atom types"""
    # surface atoms are red (O) , xtal yellow (S), others blue (N)
    nparfl = params['npartot'] - params['nparsurf']
    symbols = ['O']*params['nparsurf'] + ['N']*nparfl

    # if params contains stillsep,etc,get xtal particles
    if 'stillsep' in params:
        xparnums = bops.getxpars(positions,params)
        for pnum in xparnums:
            symbols[pnum] = 'S'
        
    # write the file
    readwrite.wxyz(fname,positions,symbols)
    return
