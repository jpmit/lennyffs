# bops.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Local bond order parameters for FFS code. The xtal particles are
identified in a Fortran routine which is wrapped by the function
getxpars() below.  The code for computing BopxBulk, the largest xtal
cluster, is implemented here. This uses the Graph class (see graph.py)

FUNCTIONS:
bopxbulk  - return size of largest xtal cluster using local bond order
            params see ten Wolde,Ruiz-Montero and Frenkel,
            Faraday Discuss. 104, 93
getxpars  - return array of particle numbers that are xtal,
            according to local bond order parameters
getxgraph - return Graph of xtal particles, with nodes that are xtal pars, and
            edges between any two pars that are neighbours.
"""

import graph
import mcfuncs

def bopxbulk(positions,params):
    """
    Return number of xtal particles (according to q6 criterion)
    and BOPxBulk, which is the largest cluster of xtal particles.
    """

    # get all xtal particles
    xpars = getxpars(positions,params)
    nxtal = len(xpars)
    
    # To compute BOPxBulk, we first create a graph with the xtal
    # particles as nodes.  Edges of the graph are then drawn between
    # neighbouring particles i.e. particles whose separation is less
    # than the cutoff distance, which is params['stillsep'].  Bopxbulk
    # is the largest connected component of the graph.
    xgraph = getxgraph(positions,params,xpars)
    comps = graph.connected_comps(xgraph)
    if nxtal > 0:
        bopxbulk = len(comps[0])
    else:
        bopxbulk = 0
    return nxtal,bopxbulk

def getxpars(positions,params):
    """Return array of xtal particle numbers."""

    npar = params['npartot']
    nsep = params['stillsep']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    thresh = params['q6link']
    minlinks = params['q6numlinks']

    # call fortran routine to get xtal particles (see bopsf.f90)
    xpars,nxtal = mcfuncs.xpars(positions[:,0],positions[:,1],
                                positions[:,2],nparsurf,params['lboxx'],
                                params['lboxy'],params['lboxz'],zperiodic,
                                nsep,minlinks,thresh)

    # Note that the fortran routine returns xpars(1:npar).  Most of
    # these entries will be zero, we only want the non-zero ones.
    # Also note we subtract 1 from ALL of the xtal particle numbers
    # returned from mcfuncs.xpars due to zero indexing in Python (and
    # one indexing in Fortran).
    return (xpars[:nxtal] - 1)

def getxgraph(positions,params,xpars):
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
