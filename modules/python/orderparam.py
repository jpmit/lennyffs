# orderparams.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Order parameter functions.  These use the extension module written in
C++.  Now deprecated functions that use the fortran extension module
can be found in orderfuncs.py .

Each of the functions below is designed to be used from the code as an
orderparameter in e.g. FFS or Umbrella Sampling simulations.  *Each
returns a tuple*, typically of length 1.

Functions that support the order parameters implemented below are in
orderfuncs.py.

FUNCTIONS:
stringify       - Convert OP tuple to a string for writing to file
allfracld_cpp   - Fractions of all polymorphs according to LD criterion.
allfracldtf_cpp - Same as above, but including TF xtal fraction.
fracld_cpp      - Fraction of xtal particles, according to LD criterion.
fractf_cpp      - Fraction of xtal particles, according to TF criterion.
nclusbcld_cpp   - Number of bcc-like and number of close-packed particles
                  in largest cluster (2d OP) according to LD criterion.
ncluscpld_cpp   - Number of close-packed particles in largest cluster
                  according to LD criterion.
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

def stringify(op):
    """Return OP tuple as a string for writing to file."""

    return ' '.join(['{:.3f}'.format(o) for o in op])

def allvx(positions, velocities, params):
    """All x co-ords of velocity."""

    return tuple(velocities[:, 0])

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
    """Return fractions of all polymorphs AND TF crystal fraction."""
    
    return allfracld_cpp(positions, params) + fractf_cpp(positions, params)


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
    
    return (frac,)


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

    return (nclus,)


def _ncluspolyld_cpp(positions, params):
    """
    Number of particles of each polymorph in largest cluster.
    """

    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    npoly = mcfuncs.ncluspolyld(positions[:,0], positions[:,1],
                                positions[:,2], npar, nparsurf,
                                params['lboxx'], params['lboxy'],
                                params['lboxz'], zperiodic, nsep,
                                usenearest)

    return npoly


def nclusbcld_cpp(positions, params):
    """
    Number of bcc-like particles in largest cluster and number of fcc
    OR hcp-like particles (i.e. number of 'close-packed' particles) in
    the largest cluster.
    """

    npoly = _ncluspolyld_cpp(positions, params)

    return (npoly[LDBCC], npoly[LDFCC] + npoly[LDHCP])


def ncluscpld_cpp(positions, params):
    """
    Number of close packed (i.e. identified as either fcc or hcp)
    particles in largest cluster.
    """

    npoly = _ncluspolyld_cpp(positions, params)

    return (npoly[LDFCC] + npoly[LDHCP],)


def nclusallcpld_cpp(positions, params):
    """
    Size of largest cluster where all particles in cluster are close
    packed (nb this is different to ncluscpld_cpp).
    """

    # get classification of all particles as np array so we can do
    # clever manipulations
    pclass = np.array(orderfuncs.ldclass(positions, params), dtype='int')

    cpclass = (pclass == LDFCC) + (pclass == LDHCP)
    npar = int(sum(cpclass)) # number of fcc or hcp particles
    zperiodic = params['zperiodic']
    nsep = params['stillsep']

    # get indices of particles in largest cluster
    clusnums = mcfuncs.largestcluster(positions[cpclass, 0], positions[cpclass, 1],
                                      positions[cpclass, 2], npar,
                                      params['lboxx'], params['lboxy'],
                                      params['lboxz'], zperiodic, nsep)
    return (len(clusnums),)


def nclusld_cpp(positions, params):
    """
    Number of particles in largest cluster, according to
    Lechner-Dellago criterion.
    """
    
    npar = params['npartot']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    nsep = params['stillsep']
    usenearest = params['usenearest']

    nclus = mcfuncs.nclusld(positions[:,0], positions[:,1],
                            positions[:,2], npar, nparsurf,
                            params['lboxx'], params['lboxy'],
                            params['lboxz'], zperiodic, nsep,
                            usenearest)

    return (nclus,)


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

    return (nclus,)


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

    return (q6,)
