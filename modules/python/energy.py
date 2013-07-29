# energy.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper to Fortran function for computing total energy (see energyf.f90)
FUNCTIONS:
len_totalenergy - computes total potential energy of system, including
                  surface, for Lennard-Jones interaction potential.
gauss_totalenergy - computes total potential energy of system, including
                    surface, gor Gaussian interaction potential.
"""

import mcfuncs

def len_totalenergy(positions,params):
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

    etot = mcfuncs.len_totalenergy(positions[:,0],positions[:,1],
                                   positions[:,2],rcut,rcsq,lboxx,
                                   lboxy,lboxz,vrc,vrc2, nparsurf,
                                   zperiodic,r6mult,r12mult)
    return etot

def gauss_totalenergy(positions,params):
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

    etot = mcfuncs.len_totalenergy(positions[:,0],positions[:,1],
                                   positions[:,2],rcut,rcsq,lboxx,
                                   lboxy,lboxz,vrc,vrc2, nparsurf,
                                   zperiodic)
    return etot
