# orderparams.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Order parameter functions.  These use the extension module written in
C++.  Now deprecated functions that use the fortran extension module
can be found in orderfuncs.py .

Each of the functions below is designed to be used from the code as an
orderparameter in e.g. FFS or Umbrella Sampling simulations.  Each
returns a tuple, typically of length 1.

Functions that support the order parameters implemented below are in
orderfuncs.py.

FUNCTIONS:
allfracld_cpp   - Fractions of all polymorphs according to LD criterion.
allfracldtf_cpp - Same as above, but including TF xtal fraction.
fracld_cpp      - Fraction of xtal particles, according to LD criterion.
fractf_cpp      - Fraction of xtal particles, according to TF criterion.
nclusld_cpp     - Number of particles in largest cluster according to LD
                  criterion.
nclustf_cpp     - Number of particles in largest cluster according to TF
                  criterion.
q6global_cpp    - Global Q6 of the system.
"""

import numpy as np

import mcfuncs
import orderfuncs

# constants needed to interface with C++ extension module
# Mapping from C++ enum type LDCLASS to int
LDNUMPOLY = 6
LDFCC = 0
LDHCP = 1
LDBCC = 2
LDLIQUID = 3
LDICOS = 4
LDSURFACE = 5
# Mapping from C++ enum type TFCLASS to int
TFLIQ = 0
TFXTAL = 1
TFSURF = 2


def allfracld_cpp(positions, params):
    """Return fractions of all polymorphs."""
    
    # get LD classification of all particles, ignoring surface
    # particles
    pclass = orderfuncs.ldclass(positions, params)[params['nparsurf']:]

    fracs = np.zeros(LDNUMPOLY)
    for p in pclass:
        fracs[p] += 1

    return tuple(fracs / (params['npartot'] - params['nparsurf']))


def allfracldtf_cpp(positions, params):
    """
    Return fractions of all polymorphs AND TF crystal fraction, as
    a string.
    """
    
    return allfracld_cpp(positions, params) + ' {:.3f}'.\
           format(fractf_cpp(positions, params))


def fracld_cpp(positions, params):
    """
    Fraction of solid particles in system, according to
    Lechner-Dellago criterion.
    """

    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    frac = mcfuncs.fracsolidld(positions[:,0], positions[:,1],
                               positions[:,2], npar, nparsurf,
                               params['lboxx'], params['lboxy'],
                               params['lboxz'], zperiodic, nsep,
                               usenearest)
    
    return frac


def fractf_cpp(positions, params):
    """
    Fraction of particles in largest cluster, according to Ten-Wolde
    Frenkel criterion.
    """
    
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    nclus = mcfuncs.fracsolidtf(positions[:,0], positions[:,1],
                                positions[:,2], npar, nparsurf,
                                params['lboxx'], params['lboxy'],
                                params['lboxz'], zperiodic, nsep,
                                minlinks, thresh, usenearest)

    return nclus


def nclusld_cpp(positions, params):
    """
    Number of particles in largest cluster, according to
    Lechner-Dellago criterion.
    """
    
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    nclus = mcfuncs.nclusld(positions[:,0], positions[:,1],
                            positions[:,2], npar, nparsurf,
                            params['lboxx'], params['lboxy'],
                            params['lboxz'], zperiodic, nsep,
                            usenearest)

    return nclus


def nclustf_cpp(positions, params):
    """
    Number of particles in largest cluster, according to Ten-Wolde
    Frenkel criterion.
    """

    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    nclus = mcfuncs.nclustf(positions[:,0], positions[:,1],
                            positions[:,2], npar, nparsurf,
                            params['lboxx'], params['lboxy'],
                            params['lboxz'], zperiodic, nsep,
                            minlinks, thresh, usenearest)

    return nclus


def q6global_cpp(positions, params):
    """
    Global order parameter Q6.
    See ten Wolde et al. Faraday Discuss. 1996
    """

    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    q6 = mcfuncs.q6global(positions[:,0], positions[:,1],
                          positions[:,2], npar, nparsurf,
                          params['lboxx'], params['lboxy'],
                          params['lboxz'], zperiodic, nsep,
                          usenearest)

    return q6
