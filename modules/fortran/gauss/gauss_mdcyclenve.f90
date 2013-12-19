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
                                  mass, npar, nsurf, zperiodic)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: ncycles, nsamp, npar, nsurf
  real(kind=db), intent(in) :: dt, rc, rcsq, vrc, vrc2
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, mass
  logical, intent(in) :: zperiodic
  ! outputs (note inout intent)
  real(kind=db), dimension(npar), intent(inout) :: xpos, ypos, zpos,&
                                                   xvel, yvel, zvel,&
                                                   fx, fy, fz

  !f2py intent(in) :: ncycles, nsamp, dt, rc, rcsq, lboxx, lboxy, lboxz
  !f2py intent(in) :: vrc, vrc2, mass, npar, nparsuf, zperiodic
  !f2py intent(in,out) :: xpos, ypos, zpos, xvel, yvel, zvel, fx, fy, fz

  real(kind=db) :: rsc, xposi, yposi, zposi, xposinew, yposinew,&
                   zposinew, eold, enew, p5dt, p5dtsq
  integer :: cy
  real(kind=db), dimension(npar) :: newfx, newfy, newfz
  real(kind=db) :: epottot, ekintot
  ! these are for cell lists
  integer :: ncelx, ncely, ncelz
  real(kind=db) :: rnx, rny, rnz
  integer, dimension(npar) :: ll
  integer, allocatable, dimension(:,:,:) :: hoc

  ! get the number of cells and build the cell list
  call getnumcells(lboxx, lboxy, lboxz, rc, ncelx, ncely, ncelz)
  write(*,*) 'num cells', ncelx, ncely, ncelz
  allocate( hoc(ncelx, ncely, ncelx) )
  
  ! precompute half of timestep and half of timestep squared
  p5dt = 0.5_db*dt
  p5dtsq = 0.5_db*(dt**2)

  do cy = 1, ncycles
     
     ! update positions using current values of velocity and force
     ! assume mass = 1 in our units for now
     xpos = xpos + xvel*dt + fx*p5dtsq
     ypos = ypos + yvel*dt + fy*p5dtsq
     zpos = zpos + zvel*dt + fz*p5dtsq

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
     xvel = xvel + (newfx + fx)*p5dt
     yvel = yvel + (newfy + fy)*p5dt
     zvel = zvel + (newfz + fz)*p5dt     

     ! the new forces are now the current forces
     fx = newfx
     fy = newfy
     fz = newfz

     ! write out diagnostics after every nsamp cycles
     if (mod(cy, nsamp) == 0) then
        ! todo: compute energy (say KE and PE) and output here (?)
        ! compute PE
        call gauss_totalenlist(ll, hoc, ncelx, ncely, ncelz, rnx, rny,&
                               rnz, xpos, ypos, zpos, rc, rcsq, lboxx,&
                               lboxy, lboxz, vrc, vrc2, npar, nsurf,&
                               zperiodic, epottot)
        write(*,*) cy, epottot
     end if

  end do
  
end subroutine gauss_executecyclesnve
