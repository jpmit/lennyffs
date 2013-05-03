# energy.py
# 26th June 2012
# James Mithen

"""
Wrapper to Fortran function for computing total energy (see energyf.f90)
FUNCTIONS:
totalenergy - computes total potential energy of system, including surface
"""

import mcfuncs

def totalenergy(positions,params):
    """Compute total energy of system, including surface"""
    npartot = params['npartot']
    rcut = params['rcut']
    rcsq = params['rcsq']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    r6mult = params['r6mult']
    r12mult = params['r12mult']

    etot = mcfuncs.totalenergy(positions[:,0],positions[:,1],positions[:,2],
                               rcut,rcsq,lboxx,lboxy,lboxz,vrc,vrc2,
                               nparsurf,zperiodic,r6mult,r12mult)
    return etot
