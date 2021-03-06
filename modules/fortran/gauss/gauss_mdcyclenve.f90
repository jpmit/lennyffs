! gauss_mccylenve.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutine for executing Molecular dynamics cycles.  We use
! the velocity Verlet algorithm.
!
! SUBROUTINES:
! gauss_executecyclesnve - execute ncycles molecular dynmics cycles
!                          note xpos,ypos,zpos and etot are returned 

subroutine gauss_executecyclesnve(xpos, ypos, zpos, xvel, yvel, zvel,&
                                  fx, fy, fz, ncycles, nsamp, dt, rc,&
                                  rcsq, vrc, vrc2, lboxx, lboxy, lboxz,&
                                  mass, npar, nsurf, zperiodic, vscale,&
                                  temp)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, intent(in) :: ncycles, nsamp, npar, nsurf
  real(kind=db), intent(in) :: dt, rc, rcsq, vrc, vrc2
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, mass, temp
  logical, intent(in) :: zperiodic, vscale

  ! outputs
  real(kind=db), dimension(npar), intent(inout) :: xpos, ypos, zpos,&
                                                   xvel, yvel, zvel,&
                                                   fx, fy, fz

  !f2py intent(in) :: ncycles, nsamp, dt, rc, rcsq, lboxx, lboxy, lboxz
  !f2py intent(in) :: vrc, vrc2, mass, npar, nparsuf, zperiodic
  !f2py intent(in,out) :: xpos, ypos, zpos, xvel, yvel, zvel, fx, fy, fz

  real(kind=db) :: rsc, xposi, yposi, zposi, xposinew, yposinew,&
                   zposinew, eold, enew, p5dt, p5dtsq
  integer :: cy, i
  real(kind=db), dimension(npar) :: newfx, newfy, newfz
  real(kind=db) :: sumv2, fs, epottot, epottot2, ekintot, ekintot2
  ! these are for cell lists
  integer :: ncelx, ncely, ncelz
  real(kind=db) :: rnx, rny, rnz
  integer, dimension(npar) :: ll
  integer, allocatable, dimension(:,:,:) :: hoc

  ! get the number of cells and build the cell list
  call getnumcells(lboxx, lboxy, lboxz, rc, ncelx, ncely, ncelz)
  write(*,*) 'num cells', ncelx, ncely, ncelz
  allocate( hoc(ncelx, ncely, ncelx) )
  ! construct the cell list
  call new_nlist(xpos, ypos, zpos, rc, lboxx, lboxy, lboxz, npar,&
                 ncelx, ncely, ncelz, ll, hoc, rnx, rny, rnz)
  
  ! precompute half of timestep and half of timestep squared
  p5dt = 0.5_db * dt
  p5dtsq = 0.5_db * (dt**2)

  ! output initial PE, KE and total energy per particle
  call gauss_totalenlist(ll, hoc, ncelx, ncely, ncelz, rnx, rny,&
                         rnz, xpos, ypos, zpos, rc, rcsq, lboxx,&
                         lboxy, lboxz, vrc, vrc2, npar, nsurf,&
                         zperiodic, epottot)
  call gauss_kineticen(xvel, yvel, zvel, npar, ekintot)

  ! we use ekintot2 and epottot2 to store total kinetic and total
  ! potential energy per particle.  Note that we output ekintot /
  ! (npar * T), which should fluctuate around the value 1.5.
  epottot2 = epottot / npar
  ekintot2 = ekintot / npar
  write(*,*) 0, epottot2, ekintot2 / temp, epottot2 + ekintot2

  do cy = 1, ncycles
     
     ! update positions using current values of velocity and force
     ! assume mass = 1 in our units for now
     xpos = xpos + xvel * dt + fx * p5dtsq
     ypos = ypos + yvel * dt + fy * p5dtsq
     zpos = zpos + zvel * dt + fz * p5dtsq

     ! apply periodic BCs
     do i = 1, npar

        ! x and y directions
        if (xpos(i) < 0.0_db) then
           xpos(i) = xpos(i) + lboxx
        else if (xpos(i) > lboxx) then
           xpos(i) = xpos(i) - lboxx
        end if
        if (ypos(i) < 0.0_db) then
           ypos(i) = ypos(i) + lboxy
        else if (ypos(i) > lboxy) then
           ypos(i) = ypos(i) - lboxy
        end if

        ! z direction
        if (zperiodic) then
           if (zpos(i) < 0.0_db) then
              zpos(i) = zpos(i) + lboxz
           else if (zpos(i) > lboxz) then
              zpos(i) = zpos(i) - lboxz
           end if
        else
           ! TODO: need to treat case of z not periodic here
        end if
        
     end do ! apply periodic BCs

     ! we need to rebuild the cell list with the new positions
     call new_nlist(xpos, ypos, zpos, rc, lboxx, lboxy, lboxz, npar,&
                    ncelx, ncely, ncelz, ll, hoc, rnx, rny, rnz)

     ! compute new forces using new positions (note we still need to
     ! keep old forces)
     call gauss_forcelist(ll, hoc, ncelx, ncely, ncelz, rnx, rny, rnz,&
                          xpos, ypos, zpos, rc, rcsq, lboxx, lboxy,&
                          lboxz, npar, nsurf, zperiodic, newfx, newfy,&
                          newfz)
  
     ! compute new velocities using new and old forces
     xvel = xvel + (newfx + fx) * p5dt
     yvel = yvel + (newfy + fy) * p5dt
     zvel = zvel + (newfz + fz) * p5dt     

     ! the new forces are now the current forces
     fx = newfx
     fy = newfy
     fz = newfz

     ! if set, we perform velocity scaling every timestep
     if (vscale) then
        sumv2 = sum(xvel**2 + yvel**2 + zvel**2) / npar
        fs = sqrt(3 * temp / sumv2)
        xvel = fs * xvel
        yvel = fs * yvel
        zvel = fs * zvel
     end if

     ! write out diagnostics after every nsamp cycles
     if (mod(cy, nsamp) == 0) then
        ! compute PE and KE
        call gauss_totalenlist(ll, hoc, ncelx, ncely, ncelz, rnx, rny,&
                               rnz, xpos, ypos, zpos, rc, rcsq, lboxx,&
                               lboxy, lboxz, vrc, vrc2, npar, nsurf,&
                               zperiodic, epottot)
        call gauss_kineticen(xvel, yvel, zvel, npar, ekintot)
        ! output PE/npar, KE/(npar*T) ~ 1.5 and total energy per particle
        epottot2 = epottot / npar
        ekintot2 = ekintot / npar
        write(*,*) cy, epottot2, ekintot2 / temp, epottot2 + ekintot2
     end if

  end do
  
end subroutine gauss_executecyclesnve
