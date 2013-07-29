# mccycle.py
# 26th June 2012
# James Mithen

"""
Wrapper to Fortran code for performing the Monte Carlo simulation This
just calls the fortran subroutine 'executecycles' (see mccyclef.f90)
for nvt simulation, and 'executecyclesnpt' (see mccyclenptf.f90) for
npt simulation.
"""

import mcfuncs

def cycle(positions,params,etot):
    """Performs the requested number of cycles of MC"""
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

    # setup and call the fortran subroutine
    xpos,ypos,zpos = positions[:,0],positions[:,1],positions[:,2]
    xpos,ypos,zpos,etot  = mcfuncs.executecycles(xpos,ypos,zpos,ncycle,nsamp,
                                                 rc,rcsq,vrc,vrc2,lboxx,lboxy,
                                                 lboxz,eps4,maxdisp,nparsurf,
                                                 zperiodic,r6mult,r12mult,etot)
    positions[:,0],positions[:,1],positions[:,2] = xpos,ypos,zpos

    return positions, etot

def cyclenpt(positions,params,etot):
    """Performs the requested number of cycles of MC"""
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

    # setup and call the fortran subroutine
    xpos,ypos,zpos = positions[:,0],positions[:,1],positions[:,2]
    xpos,ypos,zpos,lx,ly,lz,etot  = mcfuncs.executecyclesnpt(xpos,ypos,zpos,
                                                             ncycle,nsamp,
                                                             rc,rcsq,vrc,vrc2,
                                                             pressure,
                                                             lboxx,lboxy,
                                                             lboxz,eps4,maxdisp,
                                                             maxvol,
                                                             nparsurf,
                                                             zperiodic,r6mult,
                                                             r12mult,etot)
    # update box dimensions
    params['lboxx'] = lx
    params['lboxy'] = ly
    params['lboxz'] = lz
    positions[:,0],positions[:,1],positions[:,2] = xpos,ypos,zpos

    return positions, etot
