! util.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran utility functions.
!
! SUBROUTINES:
! get_nr - return number of particles with each separation 

subroutine get_nr(xpos, ypos, zpos, lboxx, lboxy, lboxz,&
                  npar, zperiodic, dr, nrvals, nprval)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, dr
  integer, intent(in) :: npar, nrvals
  logical, intent(in) :: zperiodic
  integer, dimension(nrvals), intent(out) :: nprval
  
  !f2py intent(in) :: xpos, ypos, zpos, lboxx, lboxy, lboxz, dr, npar, zperiodic, nrvals
  !f2py intent(out) :: nprval

  integer :: i, j, indx
  real(kind=db) :: sepx, sepy, sepz, p5lboxx, p5lboxy, p5lboxz, sepsq, sep

  p5lboxx = 0.5_db*lboxx
  p5lboxy = 0.5_db*lboxy
  p5lboxz = 0.5_db*lboxz

  do i=1,npar
     do j=i+1,npar
        sepx = xpos(i) - xpos(j)
        ! periodic boundary conditions
        if (sepx > p5lboxx) then
           sepx = sepx - lboxx
        else if (sepx < -p5lboxx) then
           sepx = sepx + lboxx
        end if
        sepy = ypos(i) - ypos(j)
        ! periodic boundary conditions
        if (sepy > p5lboxy) then
           sepy = sepy - lboxy
        else if (sepy < -p5lboxy) then
           sepy = sepy + lboxy
        end if
        sepz = zpos(i) - zpos(j)
        if (zperiodic) then
           ! periodic boundary conditions
           if (sepz > p5lboxz) then
              sepz = sepz - lboxz
           else if (sepz < -p5lboxz) then
              sepz = sepz + lboxz
           end if
        end if
        sepsq = sepx**2 + sepy**2 + sepz**2
        ! we could later on add a flag to this routine in case we want
        ! the squared separation rather than the separation itself
        sep = sepsq**0.5
        indx = int(sep / dr) + 1
        nprval(indx) = nprval(indx) + 1
     end do
  end do
end subroutine get_nr

!subroutine get_sk(xpos, ypos, zpos, lboxx, lboxy, lboxz,&
!                  npar, kmax)
  


!end subroutine get_sk
