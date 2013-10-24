# energy.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper to Fortran function for computing total energy (see
modules/fortran/len/len_energyf.f90 and
modules/fortran/gauss/gauss_energyf.f90).

FUNCTIONS:
len_totalenergy   - computes total potential energy of system, including
                    surface, for Lennard-Jones potential.
len_energyipar    - computes potential energy between a single particle
                    and all others, for Lennard-Jones potential.
gauss_totalenergy - computes total potential energy of system, including
                    surface, for Gaussian potential.
gauss_energyipar  - computes potential energy between a single particle
                    and all others, for Gaussian potential.
"""

import mcfuncs

def len_totalenergy(positions,params):
    """Compute total energy of system, including surface."""

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

# this function is used for the utilites but isn't used by the main
# code.
def len_energyipar(ipar, positions, params, nparsurf = None):
    """Compute energy of particle ipar with all other particles."""

    rcut = params['rcut']
    rcsq = params['rcsq']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    zperiodic = params['zperiodic']
    r6mult = params['r6mult']
    r12mult = params['r12mult']

    # particles {1,...nsurf} use r6mult and r12 mult.
    if nparsurf is None:
        nparsurf = params['nparsurf']

    # ipar + 1 converts between zero indexing (Python) and one
    # indexing (Fortran).
    etot = mcfuncs.len_energyipar(ipar+1, positions[ipar][0],
                                  positions[ipar][1],
                                  positions[ipar][2],
                                  positions[:,0],
                                  positions[:,1],
                                  positions[:,2], rcut, rcsq, lboxx,
                                  lboxy, lboxz, vrc, vrc2, nparsurf,
                                  zperiodic, r6mult, r12mult)
    return etot

def gauss_totalenlist(positions, params):
    """
    Compute total energy of system, including surface, using cell
    list.
    """    
    
    rcut = params['rcut']
    rcsq = params['rcsq']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']

    etot = mcfuncs.gauss_totalencreatelist(positions[:,0],
                                           positions[:,1],
                                           positions[:,2], rcut,
                                           rcsq, lboxx,
                                           lboxy, lboxz, vrc,
                                           vrc2, nparsurf,
                                           zperiodic)
    return etot

# this is antiquated, use gauss_totalenlist above, which uses cell
# lists for efficiency.
def gauss_totalenergy(positions, params):
    """Compute total energy of system, including surface."""
    
    rcut = params['rcut']
    rcsq = params['rcsq']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']

    etot = mcfuncs.gauss_totalenergy(positions[:,0], positions[:,1],
                                     positions[:,2], rcut, rcsq, lboxx,
                                     lboxy,lboxz, vrc, vrc2, nparsurf,
                                     zperiodic)
    return etot

# this function is used for the utilites but shouldn't be needed
# otherwise.
def gauss_energyipar(ipar, positions, params, nparsurf = None):
    """Compute energy of particle ipar with all other particles."""
    
    rcut = params['rcut']
    rcsq = params['rcsq']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']


    if nparsurf is None:
        nparsurf = params['nparsurf']
    
    # ipar + 1 converts between zero indexing (Python) and one
    # indexing (Fortran).
    etot = mcfuncs.gauss_energyipar(ipar+1, positions[ipar][0],
                                    positions[ipar][1],
                                    positions[ipar][2],
                                    positions[:,0],
                                    positions[:,1],
                                    positions[:,2], rcut, rcsq, lboxx,
                                    lboxy, lboxz, vrc, vrc2, nparsurf,
                                    zperiodic)
    return etot
