# orderparams.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Order parameters.  These use the extension modules written in Fortran
and in C++.  The suffix of each function indicates whether the
implementation is in fortran (fort) or C++ (cpp).

FUNCTIONS:
// FORTRAN FUNCTIONS (DEPRECATED)
_getxpars       - Return array of particle numbers that are xtal,
                  according to local bond order parameters.
_getxgraph      - Return Graph of xtal particles, with nodes that are xtal
                  pars, and edges between any two pars that are neighbours.
nclustf_fort    - Number of particles in largest cluster according to TF
                  criterion.
fractf_fort     - Fracition of xtal particles, according to TF criterion.
// CPP FUNCTIONS
nclustf_cpp     - Number of particles in largest cluster according to TF
                  criterion.
fractf_cpp      - Fraction of xtal particles, according to TF criterion.
nclusld_cpp     - Number of particles in largest cluster according to LD
                  criterion.
allfracld_cpp   - Fractions of all polymorphs according to LD criterion.
allfracldtf_cpp - Same as above, but including TF xtal fraction.
_ldclass        - List of particle classifications.
_ldclusnums     - Indices of particles in the largest cluster according
                  to LD criterion.
_clusnums       - Indices of particles in the largest cluster according
                  to TF criterion.
_q4w4q6w6       - Return q4bar, w4bar, q6bar, w6bar (the Lechner Dellago
                  versions) for all particles. 
"""

import numpy as np
import graph
import mcfuncs

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

##########################################
# Functions using Fortran extension module
##########################################

def _getxpars(positions,params):
    """Return array of xtal particle numbers."""
    
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    # call fortran routine to get xtal particles (see
    # modules/fortran/bopsf.f90)
    xpars, nxtal = mcfuncs.xpars(positions[:,0], positions[:,1],
                                 positions[:,2], nparsurf,
                                 params['lboxx'],
                                 params['lboxy'], params['lboxz'],
                                 zperiodic, nsep, minlinks, thresh)

    # Note that the fortran routine returns xpars(1:npar).  Most of
    # these entries will be zero, we only want the non-zero ones.
    # Also note we subtract 1 from ALL of the xtal particle numbers
    # returned from mcfuncs.xpars due to zero indexing in Python (and
    # one indexing in Fortran).
    return (xpars[:nxtal] - 1)

# since we don't have a Fortran function that computes the largest
# cluster, this is all done in Python (i.e. _getxgraph is not a
# wrapper to a Fortran function).
def _getxgraph(positions,params,xpars):
    """
    Create graph with xtal pars as nodes, edges between
    neighbouring pars.
    """

    xgraph = graph.Graph(xpars)
    stillsep = params['stillsep']
    stillsepsq = stillsep**2.0
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    zperiodic = params['zperiodic']

    # half times box width for computing periodic bcs
    p5lboxx = 0.5*lboxx
    p5lboxy = 0.5*lboxy
    p5lboxz = 0.5*lboxz
    
    for i in xpars:
        # get distances to all other solid particles
        for j in xpars:
            if i != j:
                sepx = positions[i][0] - positions[j][0]
                # periodic boundary conditions
                if (sepx > p5lboxx):
                    sepx = sepx - lboxx
                elif (sepx < -p5lboxx):
                    sepx = sepx + lboxx
                if (abs(sepx) < stillsep):
                    sepy = positions[i][1] - positions[j][1]
                    # periodic boundary conditions
                    if (sepy > p5lboxy):
                        sepy = sepy - lboxy
                    elif (sepy < -p5lboxy):
                        sepy = sepy + lboxy
                    if (abs(sepy) < stillsep):
                        sepz = positions[i][2] - positions[j][2]
                        # periodic boundary conditions
                        if (zperiodic):
                            if (sepz > p5lboxz):
                                sepz = sepz - lboxz
                            elif (sepz < -p5lboxz):
                                sepz = sepz + lboxz
                        if (abs(sepz) < stillsep):
                            # compute separation
                            rijsq = sepx**2 + sepy**2 + sepz**2
                            if rijsq < stillsepsq:
                                # particles are in same cluster
                                xgraph.add_edge(i,j)
    return xgraph

def nclustf_fort(positions, params):
    """
    Number of particles in largest cluster, according to Ten-Wolde
    Frenkel criterion.  This uses the Fortran extension module.
    Especially for large clusters, it is more efficient to use
    nclustf_cpp below; here the largest cluster calculation is done in
    Python which makes this slower.
    """
    
    # get all xtal particles
    xpars = _getxpars(positions,params)
    nxtal = len(xpars)
    
    # To compute ntf, we first create a graph with the xtal particles
    # as nodes.  Edges of the graph are then drawn between
    # neighbouring particles i.e. particles whose separation is less
    # than the cutoff distance, which is params['stillsep'].  ntf is
    # the largest connected component of the graph.
    xgraph = _getxgraph(positions, params, xpars)
    comps = graph.connected_comps(xgraph)
    if nxtal > 0:
        ntf = len(comps[0])
    else:
        ntf = 0
    return ntf

def fractf_fort(positions, params):
    """
    Fraction of solid (xtal) particles in system, according to
    Ten-Wolde Frenkel criterion.  This uses the Fortran extension
    module.
    """
    
    # get all xtal particles
    xpars = _getxpars(positions, params)
    nxtal = len(xpars)

    # we return the fraction of moving particles that are xtal
    return float(nxtal) / (params['npartot'] - params['nparsurf'])

##########################################
# Functions using C++ extension module
##########################################

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

def _tfclass(positions, params):
    """Return list of particle classifications according to TF method."""
    
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']
    usenearest = params['usenearest']

    parclass = mcfuncs.tfclass(positions[:,0], positions[:,1],
                               positions[:,2], npar, nparsurf,
                               params['lboxx'], params['lboxy'],
                               params['lboxz'], zperiodic, nsep,
                               minlinks, thresh, usenearest)
    return parclass


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

def allfracld_cpp(positions, params):
    """Return fractions of all polymorphs as a string."""
    
    # get LD classification of all particles, ignoring surface
    # particles
    pclass = _ldclass(positions, params)[params['nparsurf']:]
    fracs = np.zeros(LDNUMPOLY)
    for p in pclass:
        fracs[p] += 1
        
    # return string representation for writing to file; may have to
    # change this at some point.
    return ' '.join(['{:.3f}'.format(i) for i in
                     fracs / (params['npartot'] - params['nparsurf'])])

def allfracldtf_cpp(positions, params):
    """
    Return fractions of all polymorphs AND TF crystal fraction, as
    a string.
    """
    
    return allfracld_cpp(positions, params) + ' {:.3f}'.\
           format(fractf_cpp(positions, params))

def _ldclass(positions, params):
    """Return list of particle classifications according to LD method."""
    
    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    usenearest = params['usenearest']

    parclass = mcfuncs.ldclass(positions[:,0], positions[:,1],
                               positions[:,2], npar, nparsurf,
                               params['lboxx'], params['lboxy'],
                               params['lboxz'], zperiodic, nsep,
                               usenearest)
    return parclass

def _ldclusnums(positions, params):
    """
    Return indices of particles in largest cluster according to the LD
    criterion.
    """
    
    # get all xtal classifications
    ldclass = _ldclass(positions, params)

    # get indices of crystal particles
    xtalnums = []
    for (i, cl) in enumerate(ldclass):
        if cl in [LDFCC, LDHCP, LDBCC, LDICOS]:
            xtalnums.append(i)
    print len(xtalnums)
    
    # graph for computing largest cluster
    xgraph = _getxgraph(positions, params, xtalnums)
    comps = graph.connected_comps(xgraph)

    return comps[0]

def _clusnums(cpositions, params):
    """Return indices of particles in largest cluster."""

    # handle case of empty cpositions array
    if cpositions.shape == (0,):
        return []
    
    return mcfuncs.largestcluster(cpositions[:,0], cpositions[:,1],
                                  cpositions[:,2], len(cpositions),
                                  params['lboxx'], params['lboxy'],
                                  params['lboxz'], params['zperiodic'],
                                  params['stillsep'], params['usenearest'])

def _q4w4q6w6(positions, params):
    """Return LD order parameters."""

    ntot = params['npartot']
    q4w4q6w6 =  mcfuncs.q4w4q6w6(positions[:,0], positions[:,1],
                                 positions[:,2], ntot,
                                 params['nparsurf'],
                                 params['lboxx'], params['lboxy'],
                                 params['lboxz'], params['zperiodic'],
                                 params['stillsep'],
                                 params['usenearest'])

    # the C++ function returns a long list with q4 values followed by
    # w4 values followed by q6values followed by w6 values.  Split
    # this up here and return four lists.
    q4 = q4w4q6w6[:ntot]
    w4 = q4w4q6w6[ntot:2*ntot]
    q6 = q4w4q6w6[2*ntot:3*ntot]
    w6 = q4w4q6w6[3*ntot:4*ntot]

    return np.array(q4), np.array(w4), np.array(q6), np.array(w6)

def _numneighcut(positions, params):
    """Return list of number of neighbours of each particle."""

    nn = mcfuncs.numneighcut(positions[:,0], positions[:,1],
                             positions[:,2], params['npartot'],
                             params['lboxx'],
                             params['lboxy'],
                             params['lboxz'],
                             params['zperiodic'],
                             params['nsep'])
                             
    return nn

def nothing(positions, params):
    """No order parameter."""

    return 0.0
