! gauss_force.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Subroutines for computing force on every particle, using the cell
! list construction for efficiency.
!
! SUBROUTINES:
! gauss_fij             - compute force on particle i due to particle j
! gauss_forcecreatelist - create cell list then return force on each par
! gauss_forcelist       - compute force on every particle using cell list

subroutine gauss_fij(ipar, jpar, xposi, yposi, zposi, xposj, yposj, zposj,&
                     lboxx, lboxy, lboxz, rc, rcsq,  npar, nsurf,&
                     zperiodic, fx, fy, fz)
  !!! Compute force on particle i due to particle j, taking into
  !!! account periodic bcs.
  
  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)
  
  ! inputs
  integer, intent(in) :: ipar, jpar
  real(kind=db), intent(in) :: xposi, yposi, zposi, xposj, yposj, zposj
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, rc, rcsq
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic

  ! outputs
  real(kind=db), intent(out) :: fx, fy, fz

  !f2py intent(in) :: ipar, jpar, xposi, yposi, zposi, xposj, yposj, zposj
  !f2py intent(in) :: lboxx, lboxy, lboxz, rc, rcsq, npar, nsurf, zperiodic
  !f2py intent(out) :: fx, fy, fz

  real(kind=db) :: sepx, sepy, sepz, sepsq, prefac

  fx = 0.0_db
  fy = 0.0_db
  fz = 0.0_db
  ! note the ordering of the separation here; this is important for
  ! getting the direction of the force correct.  If sepx is positive,
  ! the force on particle i should be in the negative x direction.
  sepx = xposj - xposi
  ! periodic boundary conditions
  if (sepx > 0.5 * lboxx) then
     sepx = sepx - lboxx
  else if (sepx < -0.5 * lboxx) then
     sepx = sepx + lboxx
  end if

  if (abs(sepx) < rc) then
     sepy = yposj - yposi
     ! periodic boundary conditions
     if (sepy > 0.5*lboxy) then
        sepy = sepy - lboxy
     else if (sepy < -0.5*lboxy) then
        sepy = sepy + lboxy
     end if

     if (abs(sepy) < rc) then
        sepz = zposj - zposi
        if (zperiodic) then
           ! periodic boundary conditions
           if (sepz > 0.5 * lboxz) then
              sepz = sepz - lboxz
           else if (sepz < -0.5 * lboxz) then
              sepz = sepz + lboxz
           end if
        end if
        sepsq = sepx**2 + sepy**2 + sepz**2
        if (sepsq < rcsq) then
           ! add contribution to total potential energy
           prefac = -2*exp(-sepsq)
           if (ipar > nsurf .and. jpar > nsurf) then
              ! both particles are fluid particles
              fx = sepx * prefac
              fy = sepy * prefac
              fz = sepz * prefac
           else
              ! at least one particle is a surface particle
              ! currently we are doing the same thing, but may
              ! want to add a different potential later
              fx = sepx * prefac
              fy = sepy * prefac
              fz = sepz * prefac
           end if
        end if
     end if
  end if
end subroutine gauss_fij

subroutine gauss_forcecreatelist(xpos, ypos, zpos, rc, rcsq,&
                                 lboxx, lboxy, lboxz,&
                                 npar, nsurf, zperiodic, fx,&
                                 fy, fz)
  ! Create the cell list and then return force on each particle

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic

  ! outputs
  real(kind=db), dimension(npar), intent(out) :: fx, fy, fz

  !f2py intent(in) :: xpos, ypos, zpos, rc, rcsq, lboxx, lboxy, lboxz
  !f2py intent(in) :: npar, nsurf, zperiodic
  !f2py intent(out) :: fx, fy, fz
  !f2py depend(npar) :: xpos, ypox, zpos, fx, fy, fz

  integer, dimension(npar) :: ll
  integer, allocatable, dimension(:,:,:) :: hoc  
  integer :: ncelx, ncely, ncelz, status
  real(kind=db) :: rnx, rny, rnz
  
  ! get the number of cells and build the cell list
  call getnumcells(lboxx, lboxy, lboxz, rc, ncelx, ncely, ncelz)
  allocate(hoc(ncelx, ncely, ncelx))
  call new_nlist(xpos, ypos, zpos, rc, lboxx, lboxy, lboxz, npar,&
                 ncelx, ncely, ncelz, ll, hoc, rnx, rny, rnz)
  ! get total energy using cell list
  call gauss_forcelist(ll, hoc, ncelx, ncely, ncelz, rnx, rny, rnz,&
                       xpos, ypos, zpos, rc, rcsq, lboxx, lboxy,&
                       lboxz, npar, nsurf, zperiodic, fx, fy, fz)
  ! deallocate the head of cell array
  deallocate(hoc, STAT=status)

end subroutine gauss_forcecreatelist

subroutine gauss_forcelist(ll, hoc, ncelx, ncely, ncelz,&
                           rnx, rny, rnz, xpos, ypos, zpos,&
                           rc, rcsq, lboxx, lboxy, lboxz, npar,&
                           nsurf, zperiodic, fx, fy, fz)
  ! Compute force on every particle using cell list for efficiency

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, dimension(npar), intent(in) :: ll
  integer, intent(in) :: ncelx, ncely, ncelz
  real(kind=db), intent(in) :: rnx, rny, rnz
  integer, dimension(ncelx, ncely, ncelz), intent(in) :: hoc
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic
  
  ! outputs
  real(kind=db), dimension(npar), intent(out) :: fx, fy, fz

  real(kind=db) :: eij, xposi, yposi, zposi, fxij, fyij, fzij
  integer :: icelx, icely, icelz, cellnum, jpar, celx, cely, celz, &
             nceltot, ipar

  fx = 0.0_db
  fy = 0.0_db
  fz = 0.0_db

  do ipar = 1, npar

     xposi = xpos(ipar)
     yposi = ypos(ipar)
     zposi = zpos(ipar)
     
     ! determine cell that particle i is in
     icelx = int(xposi / rnx) + 1
     icely = int(yposi / rny) + 1
     icelz = int(zposi / rnz) + 1

     ! go through each cell in turn (27 in total in three dimensions),
     ! and add pot energy between particle i and all particles in the
     ! cell.  Note if
     if (ncelx == 1) then
        nceltot = 1
     else
        nceltot = 27
     end if

     do cellnum = 1, nceltot

        ! get the next cell indexes (ncelx, ncely, ncelz)
        call cellindx(cellnum, icelx, icely, icelz,& ! cell of particle i
                      ncelx, ncely, ncelz,&          ! total num cells in each dim
                      celx, cely, celz)              ! cell index we want

        jpar = hoc(celx, cely, celz)

        do while (jpar /= 0)
           if (ipar /= jpar) then

              ! get force on particle i due to particle j
              call gauss_fij(ipar, jpar, xposi, yposi, zposi, xpos(jpar),&
                             ypos(jpar), zpos(jpar),&
                             lboxx, lboxy, lboxz, rc, rcsq, npar, nsurf,&
                             zperiodic, fxij, fyij, fzij)
              fx(ipar) = fx(ipar) + fxij
              fy(ipar) = fy(ipar) + fyij
              fz(ipar) = fz(ipar) + fzij              
           end if
           ! next particle in this cell
           jpar = ll(jpar)
        end do
     end do
  end do

end subroutine gauss_forcelist
