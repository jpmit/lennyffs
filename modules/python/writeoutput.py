# writeoutput.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Functions for writing output of FFS simulations.

FUNCTIONS:
writepickparams - write parameters dictionary to pickle file.
writeparams     - write dictionary of parameters in human readable
                  form (i.e. not pickled).
writemdpick     - write pickle file for restarting MD simulations,
                  this includes the positions as well as the velocities.
readmdpick      - read pickle file for restarting MD simulations,
                  return positions and velocities.
writexyz_tf     - write xyz file with particle symbols as classified
                  by TF method.  Note this will compute order
                  parameters.
writexyz_ld     - write xyz file with particle symbols as classified
                  by LD method.  Note this will compute order
                  parameters. 
"""

import os
import pickle

import orderfuncs
import orderparam
import readwrite

# lookup table for Lechner-Dellago (LD) method symbols
_LD_SYMBOLS = {orderparam.LDFCC     : 'S', # FCC (yellow in jmol)
               orderparam.LDHCP     : 'P', # HCP (orange in jmol)
               orderparam.LDBCC     : 'F', # BCC (green in jmol)
               orderparam.LDLIQUID  : 'N', # LIQUID (blue in jmol)
               orderparam.LDICOS    : 'B', # ICOS (pink in jmol)
               orderparam.LDSURFACE : 'O'  # SURFACE (red in jmol)
               }
# lookup table for Ten-Wolde-Frenkel (TF) method symbols
_TF_SYMBOLS = {orderparam.TFLIQ     : 'N',
               orderparam.TFXTAL    : 'S',
               orderparam.TFSURF    : 'O'
               }

def writepickparams(params, fname='params.pkl'):
    """Write out parameters dictionary as .pkl file."""
    
    fout = open(fname, 'wb')
    pickle.dump(params, fout)
    fout.close()
    return


def writeparams(params, fname='params.out'):
    """Write out dictionary of parameters in human readable form."""
    
    fout = open(fname, 'w')
    outstr = ('# Simulation parameters (Human readable version of params.pkl)\n'
              '# Note that only some of these parameters are used\n'
              '# Others are redundant and could prove downright misleading!\n')
    keys = params.keys()
    keys.sort()
    for k in keys:
        outstr = '{0}{1} {2}\n'.format(outstr, k, str(params[k]))
    fout.write(outstr)
    fout.close()
    return


def writemdpick(fname, positions, velocities):
    """
    Write pickle file that includes both positions and velocities.
    """

    topick = [positions, velocities]
    
    fout = open(fname,'wb')
    pickle.dump(topick, fout)
    fout.close()


def readmdpick(fname):
    """
    Read file and return positions and velocities.
    """

    pfile = open(fname, 'rb')
    positions, velocities = pickle.load(pfile)
    pfile.close()
    return positions, velocities


def _get_box_kwargs(params):            
    """If npt simulation, return dictionary with box dims.

    This is used so that we can write box dimensions to the XYZ file.
    """

    if params['mctype'] == 'npt':
        return {'boxdims': [params['lboxx'], params['lboxy'],
                            params['lboxz']]}
    return {}


def writexyz_tf(fname, positions, params):
    """
    Write positions in xyz format with symbols as atom types, using
    TF approach for identifying crystalline particles.
    """

    # get TF classification
    tfclass = orderfuncs.tfclass(positions, params)
    # give each atom the correct symbol using lookup table
    symbols = [_TF_SYMBOLS[i] for i in tfclass]

    kwargs = _get_box_kwargs(params)
        
    # write the file
    readwrite.wxyz(fname, positions, symbols, **kwargs)
    return


def writexyz_ld(fname, positions, params):
    """
    Write positions in xyz format with symbols as atom types, using
    LD approach for identifying crystalline particles.
    """

    # get LD classification
    ldclass = orderfuncs.ldclass(positions, params)
    # give each atom the correct symbol using lookup table
    symbols = [_LD_SYMBOLS[i] for i in ldclass]

    kwargs = _get_box_kwargs(params)

    # write the file
    readwrite.wxyz(fname, positions, symbols, **kwargs)


def writexyz_noop(fname, positions, params):
    """
    Write positions in xyz format without computing any order
    parameter with which to classify the particles as xtal, liquid
    etc.  Here we simply make the symbol for every particle in the xyz
    file 'N' (Nitrogen), which is blue in jmol.
    """

    readwrite.wxyz(fname, positions, ['N']*len(positions))
