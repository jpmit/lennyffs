# orderparams.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Local bond order parameters for FFS code
the xtal particles are identified in a Fortran routine
which is wrapped by the function getxpars() below.
The code for computing BopxBulk, the largest xtal cluster, is implemented here
This uses the Graph class (see graph.py)
FUNCTIONS:
bopxbulk - return size of largest xtal cluster using local bond order params
           see ten Wolde,Ruiz-Montero and Frenkel, Faraday Discuss. 104, 93
getxpars - return array of particle numbers that are xtal,
           according to local bond order parameters
getxgraph - return Graph of xtal particles, with nodes that are xtal pars, and
            edges between any two pars that are neighbours.
"""

import graph
import opfunctions
import mcfuncs

def ntf(positions, params):
    """Number of particles in largest cluster, according to Ten-Wolde Frenkel
    criterion"""
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

def nld(positions, params):
    """Number of particles in largest cluster, according to
    Lechner-Dellago criterion"""
    # to do
    return 0

def fractf(positions, params):
    """Fraction of solid particles in system, according to
    Ten-Wolde Frenkel criterion"""
    # get all xtal particles
    xpars = opfunctions.getxpars(positions,params)
    nxtal = len(xpars)
    
    return nxtal / (params['npartot'] - params['nparsurf'])

def fracld(positions, params):
    """Fraction of solid particles in system, according to
    Lechner-Dellago criterion"""
    # to do
    return 0.0

def default(positions, params):
    """Default is no order parameter"""
    return 0.0
