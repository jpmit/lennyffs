# mccycle.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Wrapper to Fortran code for performing the Monte Carlo simulation This
just calls the fortran subroutine 'executecycles' (see mccyclef.f90)
for nvt simulation, and 'executecyclesnpt' (see mccyclenptf.f90) for
npt simulation.
"""

import mcfuncs

def len_cyclenvt(positions,params,etot):
    """Performs the requested number of cycles of NVT MC"""
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
    xpos,ypos,zpos,etot  = mcfuncs.len_executecyclesnvt(xpos,ypos,zpos,
                                                        ncycle,nsamp,
                                                        rc,rcsq,vrc,vrc2,
                                                        lboxx,lboxy,
                                                        lboxz,eps4,
                                                        maxdisp,nparsurf,
                                                        zperiodic,r6mult,
                                                        r12mult,etot)
    positions[:,0],positions[:,1],positions[:,2] = xpos,ypos,zpos

    return positions, etot

def len_cyclenpt(positions,params,etot):
    """Performs the requested number of cycles of NPT MC"""
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
    xpos,ypos,zpos,lx,ly,lz,etot  = mcfuncs.len_executecyclesnpt(xpos,ypos,zpos,
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

def gauss_cyclenvt(positions,params,etot):
    """Performs the requested number of cycles of NVT MC"""
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

    # setup and call the fortran subroutine
    xpos,ypos,zpos = positions[:,0],positions[:,1],positions[:,2]
    xpos,ypos,zpos,etot  = mcfuncs.gauss_executecyclesnvt(xpos,ypos,zpos,
                                                          ncycle,nsamp,
                                                          rc,rcsq,vrc,vrc2,
                                                          lboxx,lboxy,
                                                          lboxz,epsovert,
                                                          maxdisp,nparsurf,
                                                          zperiodic,
                                                          etot)
    positions[:,0],positions[:,1],positions[:,2] = xpos,ypos,zpos

    return positions, etot

def gauss_cyclenpt(positions,params,etot):
    """Performs the requested number of cycles of NPT MC"""
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

    # setup and call the fortran subroutine
    xpos,ypos,zpos = positions[:,0],positions[:,1],positions[:,2]
    xpos,ypos,zpos,lx,ly,lz,etot  = mcfuncs.gauss_executecyclesnpt(xpos,ypos,zpos,
                                                                   ncycle,nsamp,
                                                                   rc,rcsq,vrc,vrc2,
                                                                   pressure,
                                                                   lboxx,lboxy,
                                                                   lboxz,epsovert,maxdisp,
                                                                   maxvol,
                                                                   nparsurf,
                                                                   zperiodic,
                                                                   etot)
    # update box dimensions
    params['lboxx'] = lx
    params['lboxy'] = ly
    params['lboxz'] = lz
    positions[:,0],positions[:,1],positions[:,2] = xpos,ypos,zpos

    return positions, etot
