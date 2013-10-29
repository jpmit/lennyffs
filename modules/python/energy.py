# energy.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper to Fortran functions for computing total potential energy (see
modules/fortran/len/len_energy.f90 and
modules/fortran/gauss/gauss_energy.f90 for the Fortran
implementations).

FUNCTIONS:
len_totalenergy   - computes total potential energy of system, including
                    surface, for Lennard-Jones potential.
len_energyipar    - computes potential energy between a single particle
                    and all others, for Lennard-Jones potential.
gauss_totalenergy - computes total potential energy of system, including
                    surface, for Gaussian potential.
gauss_energyipar  - computes potential energy between a single particle
                    and all others, for Gaussian potential.
gauss_totalenlist - computes total potential energy of system, including
                    surface, for Gaussian potential.  This uses cell lists
                    for efficiency, and therefore should be preferred to
                    gauss_totalenergy above.
"""

import mcfuncs

# this is deprecated, use len_totalenlist below, which uses cell lists
# for efficiency.
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

    etot = mcfuncs.len_totalenergy(positions[:,0], positions[:,1],
                                   positions[:,2], rcut, rcsq, lboxx,
                                   lboxy, lboxz, vrc, vrc2, nparsurf,
                                   zperiodic, r6mult, r12mult)
    return etot

# this function is used for the utilites but shouldn't be needed by
# the main code.
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

    # particles {1,...,nsurf} use r6mult and r12 mult.
    if nparsurf is None:
        nparsurf = params['nparsurf']

    # ipar + 1 converts between zero indexing (Python) and one
    # indexing (Fortran).
    etot = mcfuncs.len_energyipar(ipar + 1, positions[ipar][0],
                                  positions[ipar][1],
                                  positions[ipar][2],
                                  positions[:,0], positions[:,1],
                                  positions[:,2], rcut, rcsq, lboxx,
                                  lboxy, lboxz, vrc, vrc2, nparsurf,
                                  zperiodic, r6mult, r12mult)
    return etot

def len_totalenlist(positions, params):
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
    r6mult = params['r6mult']
    r12mult = params['r12mult']

    etot = mcfuncs.len_totalencreatelist(positions[:,0],
                                         positions[:,1],
                                         positions[:,2], rcut,
                                         rcsq, lboxx,
                                         lboxy, lboxz, vrc,
                                         vrc2, nparsurf,
                                         zperiodic, r6mult,
                                         r12mult)
    return etot

# this is deprecated, use gauss_totalenlist below, which uses cell
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
# by the main code.
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
