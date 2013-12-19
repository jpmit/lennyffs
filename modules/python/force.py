# force.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper to Fortran functions for computing force (see
modules/fortran/gauss/gauss_force.f90 for the Fortran implementation).

FUNCTIONS:
gauss_forceslist - computes force on every particle in the system
                   potential energy of system, including surface,
                   for Gaussian potential.  This uses cell lists
                   for efficiency.
"""

def gauss_forceslist(positions, params):
    """Compute force on every particle in system, including surface."""

    rcut = params['rcut']
    rcsq = params['rcsq']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    
    # note the function creates the cell list for us
    fx, fy, fz = mcfuncs.\
                 gauss_forcecreatelist(positions[:,0], positions[:,1],
                                       positions[:,2], rcut, rcsq,
                                       lboxx, lboxy, lboxz, nparsurf,
                                       zperiodic)

    # unsure that this is the most efficient way to create this array
    forces = np.empty([len(fx), 3])
    forces[:,0], forces[:,1], forces[:,2] = fx, fy, fz

    return forces
