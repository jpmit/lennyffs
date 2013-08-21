# orderparams.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Order parameters for FFS code.  These use the extension modules
written in Fortran and in C++.

FUNCTIONS:

"""

import numpy as np
import graph
import opfunctions
import mcfuncs

# Functions using Fortran extension module

def nclustf_fort(positions, params):
    """Number of particles in largest cluster, according to Ten-Wolde Frenkel
    criterion.  This uses the Fortran extension module"""
    # get all xtal particles
    xpars = opfunctions.getxpars(positions,params)
    nxtal = len(xpars)
    
    # To compute ntf, we first create a graph with the xtal particles
    # as nodes.  Edges of the graph are then drawn between
    # neighbouring particles i.e. particles whose separation is less
    # than the cutoff distance, which is params['stillsep'].  ntf is
    # the largest connected component of the graph.
    xgraph = opfunctions.getxgraph(positions,params,xpars)
    comps = graph.connected_comps(xgraph)
    if nxtal > 0:
        ntf = len(comps[0])
    else:
        ntf = 0
    return ntf

def fractf_fort(positions, params):
    """Fraction of solid particles in system, according to Ten-Wolde
    Frenkel criterion.  This uses the Fortran extension module."""
    # get all xtal particles
    xpars = opfunctions.getxpars(positions,params)
    nxtal = len(xpars)
    
    return float(nxtal) / (params['npartot'] - params['nparsurf'])

# Functions using C++ extension module
# Mapping from C++ enum type LDCLASS to int
LDNUMPOLY = 6
LDFCC = 0
LDHCP = 1
LDBCC = 2
LDLIQUID = 3
LDICOS = 4
LDSURFACE = 5

def nclustf_cpp(positions, params):
    """Number of particles in largest cluster, according to Ten-Wolde
    Frenkel criterion.  This uses the C++ extension module"""
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    nclus = mcfuncs.nclustf(positions[:,0], positions[:,1],
                            positions[:,2], npar, nparsurf,
                            params['lboxx'], params['lboxy'],
                            params['lboxz'],zperiodic, nsep,
                            minlinks,thresh)

    return nclus

def fractf_cpp(positions, params):
    """Fract of particles in largest cluster, according to Ten-Wolde
    Frenkel criterion.  This uses the C++ extension module"""
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    nclus = mcfuncs.fracsolidtf(positions[:,0], positions[:,1],
                                positions[:,2], npar, nparsurf,
                                params['lboxx'], params['lboxy'],
                                params['lboxz'],zperiodic, nsep,
                                minlinks,thresh)

    return nclus

def nclusld_cpp(positions, params):
    """Number of particles in largest cluster, according to
    Lechner-Dellago criterion.  This uses the C++ extension module"""
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    nclus = mcfuncs.nclusld(positions[:,0], positions[:,1],
                            positions[:,2], npar, nparsurf,
                            params['lboxx'], params['lboxy'],
                            params['lboxz'],zperiodic, nsep)

    return nclus

def fracld_cpp(positions, params):
    """Fraction of solid particles in system, according to
    Lechner-Dellago criterion"""
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    frac = mcfuncs.fracsolidld(positions[:,0], positions[:,1],
                               positions[:,2], npar, nparsurf,
                               params['lboxx'], params['lboxy'],
                               params['lboxz'],zperiodic, nsep)
    
    return frac

def allfracld_cpp(positions, params):
    """Return fractions of all polymorphs"""
    # get LD classification of all particles, ignore surface particles
    pclass = _ldclass(positions, params)[params['nparsurf']:]
    fracs = np.zeros(LDNUMPOLY)
    for p in pclass:
        fracs[p] += 1
    # return string representation for writing to file; may have to
    # change this at some point
    return ' '.join(['{:.3f}'.format(i) for i in
                     fracs / (params['npartot'] - params['nparsurf'])])

def allfracldtf_cpp(positions, params):
    """Fractions of all polymorphs and TF crystal fraction"""
    return allfracld_cpp(positions, params) + ' {:.3f}'.\
           format(fractf_cpp(positions, params))

def default(positions, params):
    """Default is no order parameter"""
    return 0.0

def _ldclass(positions, params):
    """List of particle classifications"""
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    parclass = mcfuncs.ldclass(positions[:,0], positions[:,1],
                               positions[:,2], npar, nparsurf,
                               params['lboxx'], params['lboxy'],
                               params['lboxz'],zperiodic, nsep)
    return parclass
