# mccycle.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper to Fortran code for performing the Monte Carlo
simulation. This just calls the relevant fortran subroutine.

FUNCTIONS:
ipl_cyclenvt   - NVT MC for IPL potential.
ipl_cyclenpt   - NPT MC for IPL potential.
len_cyclenvt   - NVT MC for Lennard-Jones potential.
len_cyclenpt   - NPT MC for Lennard-Jones potential.
gauss_cyclenvt - NVT MC for Gaussian potential.
gauss_cyclenpt - NPT MC for Gaussian potential.
gauss_cyclemd  - NVE MD (not MC!) for Gaussian potential.
"""

import mcfuncs

def ipl_cyclenvt(positions, params, etot):
    """Performs the requested number of cycles of NVT MC."""

    # ncycle is the number of MC cycles we perform (one cycle consists
    # of on average a single displacement move per moving particle).
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    epsovert = 1.0/params['Tstar']    
    maxdisp = params['maxdisp']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    potexponent = params['potexponent']
    ss = params['sameseed']

    # setup and call the fortran subroutine
    xpos,ypos,zpos = positions[:,0],positions[:,1],positions[:,2]
    xpos, ypos, zpos, etot  = mcfuncs.\
                              ipl_executecyclesnvt(xpos, ypos, zpos,
                                                   ncycle, nsamp,
                                                   rc, rcsq, vrc,
                                                   vrc2, lboxx, lboxy,
                                                   lboxz, epsovert,
                                                   maxdisp, nparsurf,
                                                   zperiodic, ss,
                                                   potexponent,
                                                   etot)
    
    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos

    return positions, etot

def ipl_cyclenpt(positions, params, etot):
    """Performs the requested number of cycles of NPT MC."""

    # ncycle is the number of MC cycles we perform (one cycle consists
    # of on average a single displacement move per moving particle AND
    # a single volume move).
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    epsovert = 1.0/params['Tstar']    
    maxdisp = params['maxdisp']
    maxvol = params['maxvol']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    potexponent = params['potexponent']
    pressure = params['pressure']
    ss = params['sameseed']

    # setup and call the fortran subroutine
    xpos, ypos, zpos = positions[:,0], positions[:,1], positions[:,2]
    xpos,ypos,zpos,lx,ly,lz,etot  = mcfuncs.\
                                    ipl_executecyclesnpt(xpos, ypos,
                                                         zpos,
                                                         ncycle, nsamp,
                                                         rc, rcsq, vrc,
                                                         vrc2,
                                                         pressure,
                                                         lboxx, lboxy,
                                                         lboxz, epsovert,
                                                         maxdisp,
                                                         maxvol,
                                                         nparsurf,
                                                         zperiodic, ss,
                                                         potexponent, etot)
    # update box dimensions
    params['lboxx'] = lx
    params['lboxy'] = ly
    params['lboxz'] = lz
    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos

    return positions, etot

def len_cyclenvt(positions, params, etot):
    """Performs the requested number of cycles of NVT MC."""

    # ncycle is the number of MC cycles we perform (one cycle consists
    # of on average a single displacement move per moving particle).
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    eps4 = 4.0/params['Tstar']    
    maxdisp = params['maxdisp']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    r6mult = params['r6mult']
    r12mult = params['r12mult']
    ss = params['sameseed']

    # setup and call the fortran subroutine
    xpos,ypos,zpos = positions[:,0],positions[:,1],positions[:,2]
    xpos, ypos, zpos, etot  = mcfuncs.\
                              len_executecyclesnvt(xpos, ypos, zpos,
                                                   ncycle, nsamp,
                                                   rc, rcsq, vrc,
                                                   vrc2, lboxx, lboxy,
                                                   lboxz, eps4,
                                                   maxdisp, nparsurf,
                                                   zperiodic, ss,
                                                   r6mult, r12mult,
                                                   etot)
    
    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos

    return positions, etot

def len_cyclenpt(positions, params, etot):
    """Performs the requested number of cycles of NPT MC."""

    # ncycle is the number of MC cycles we perform (one cycle consists
    # of on average a single displacement move per moving particle AND
    # a single volume move).
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    eps4 = 4.0/params['Tstar']    
    maxdisp = params['maxdisp']
    maxvol = params['maxvol']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    r6mult = params['r6mult']
    r12mult = params['r12mult']
    pressure = params['pressure']
    ss = params['sameseed']

    # setup and call the fortran subroutine
    xpos, ypos, zpos = positions[:,0], positions[:,1], positions[:,2]
    xpos,ypos,zpos,lx,ly,lz,etot  = mcfuncs.\
                                    len_executecyclesnpt(xpos, ypos,
                                                         zpos,
                                                         ncycle, nsamp,
                                                         rc, rcsq, vrc,
                                                         vrc2,
                                                         pressure,
                                                         lboxx, lboxy,
                                                         lboxz, eps4,
                                                         maxdisp,
                                                         maxvol,
                                                         nparsurf,
                                                         zperiodic, ss,
                                                         r6mult,
                                                         r12mult, etot)
    # update box dimensions
    params['lboxx'] = lx
    params['lboxy'] = ly
    params['lboxz'] = lz
    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos

    return positions, etot

def gauss_cyclenvt(positions, params, etot):
    """Performs the requested number of cycles of NVT MC."""

    # ncycle is the number of MC cycles we perform (one cycle consists
    # of on average a single displacement move per moving particle).    
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    epsovert = 1.0/params['Tstar']    
    maxdisp = params['maxdisp']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    ss = params['sameseed']

    # setup and call the fortran subroutine
    xpos, ypos, zpos = positions[:,0], positions[:,1], positions[:,2]
    xpos, ypos, zpos, etot  = mcfuncs.\
                              gauss_executecyclesnvt(xpos, ypos, zpos,
                                                     ncycle, nsamp,
                                                     rc, rcsq, vrc,
                                                     vrc2, lboxx,
                                                     lboxy, lboxz,
                                                     epsovert,
                                                     maxdisp, nparsurf,
                                                     zperiodic, ss,
                                                     etot)
    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos

    return positions, etot

def gauss_cyclenpt(positions, params, etot):
    """Performs the requested number of cycles of NPT MC."""
    
    # ncycle is the number of MC cycles we perform (one cycle consists
    # of on average a single displacement move per moving particle AND
    # a single volume move).
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    epsovert = 1.0/params['Tstar']    
    maxdisp = params['maxdisp']
    maxvol = params['maxvol']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    pressure = params['pressure']
    ss = params['sameseed']

    # setup and call the fortran subroutine
    xpos, ypos, zpos = positions[:,0], positions[:,1], positions[:,2]
    xpos, ypos, zpos, lx, ly, lz, \
          etot  = mcfuncs.gauss_executecyclesnpt(xpos, ypos, zpos,
                                                 ncycle, nsamp, rc,
                                                 rcsq, vrc, vrc2,
                                                 pressure, lboxx,
                                                 lboxy, lboxz,
                                                 epsovert, maxdisp,
                                                 maxvol, nparsurf,
                                                 zperiodic, ss, etot)
    # update box dimensions
    params['lboxx'] = lx
    params['lboxy'] = ly
    params['lboxz'] = lz
    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos

    return positions, etot

def gauss_cyclemd(positions, params, velocities, forces):
    """Performs the requested number of cycles of NVE MD (not MC!)."""

    # ncycle is the number of MD cycles we perform (one cycle consists
    # of a single timestep).
    ncycle = params['cycle']
    nsamp = params['nsamp']
    rc = params['rcut']
    rcsq = params['rcsq']
    vrc = params['vrc']
    vrc2 = params['vrc2']    
    lboxx = params['lboxx']
    lboxy = params['lboxy']
    lboxz = params['lboxz']
    nparsurf = params['nparsurf']
    zperiodic = params['zperiodic']
    dt = params['dt']
    mass = params['mass']
    vscale = params['vscale']
    temp = params['Tstar']

    # setup and call the fortran subroutine
    xpos, ypos, zpos = positions[:,0], positions[:,1], positions[:,2]
    xvel, yvel, zvel = velocities[:,0], velocities[:,1], velocities[:,2]
    fx, fy, fz = forces[:,0], forces[:,1], forces[:,2]    
    xpos, ypos, zpos,\
    xvel, yvel, zvel,\
    fx, fy, fz  = mcfuncs.gauss_executecyclesnve(xpos, ypos, zpos,
                                                 xvel, yvel, zvel,
                                                 fx, fy, fz,
                                                 ncycle, nsamp,
                                                 dt, rc, rcsq,
                                                 vrc, vrc2,
                                                 lboxx,
                                                 lboxy, lboxz,
                                                 mass,
                                                 nparsurf,
                                                 zperiodic,
                                                 vscale, temp)

    positions[:,0], positions[:,1], positions[:,2] = xpos, ypos, zpos
    velocities[:,0], velocities[:,1], velocities[:,2] = xvel, yvel, zvel
    forces[:,0], forces[:,1], forces[:,2] = fx, fy, fz

    return positions, velocities, forces
