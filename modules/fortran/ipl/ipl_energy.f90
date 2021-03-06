! ipl_energyf.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutines for computing potential energy (p.e.) for
! Inverse Power Law interaction.
!
! SUBROUTINES:
! ipl_totalenergy       - compute total p.e. of system by summing
!                         p.e. between all pairs of particles
! ipl_energyipar        - compute p.e. between a given particle
!                         and all other particles
! ipl_eij               - compute p.e. between particles i and j
! ipl_totalencreatelist - create the cell list and return total
!                         p.e.
! ipl_totalenlist       - compute total p.e. using cell lists
! ipl_enlist            - compute total p.e. of particle i using
!                         cell lists

subroutine ipl_totalenergy(xpos, ypos, zpos, rc, rcsq,&
                           lboxx, lboxy, lboxz, vrc, vrc2, npar,&
                           nsurf, zperiodic, potexponent, epot)
   !!! Computes total potential energy of system (in units of epsilon)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, intent(in) :: npar, nsurf
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz
  real(kind=db), intent(in) :: vrc, vrc2, potexponent
  logical, intent(in) :: zperiodic
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos

  ! outputs
  real(kind=db), intent(out) :: epot

  !f2py intent(in) :: xpos, ypos, zpos, rc, rcsq, lboxx, lboxy, lboxz
  !f2py intent(in) :: vrc, vrc2, npar, nsurf, zperiodic, potexponent
  !f2py intent(out) :: epot

  integer :: ipar, jpar
  real(kind=db) :: sepx, sepy, sepz, sepsq, rexponenti, halfexpminus2
  real(kind=db) :: p5lboxx, p5lboxy, p5lboxz

  halfexpminus2 = (potexponent - 2.0_db) / 2.0_db
  p5lboxx = 0.5_db * lboxx ! for periodic bcs
  p5lboxy = 0.5_db * lboxy
  p5lboxz = 0.5_db * lboxz
  epot = 0.0_db
  do ipar = 1, npar - 1
     do jpar = ipar + 1, npar
        sepx = xpos(ipar) - xpos(jpar)
        ! periodic boundary conditions
        if (sepx > p5lboxx) then
           sepx = sepx - lboxx
        else if (sepx < -p5lboxx) then
           sepx = sepx + lboxx
        end if
        !sepx = sepx - lboxx * nint(sepx/lboxx)
        if (abs(sepx) < rc) then
           sepy = ypos(ipar) - ypos(jpar)
           ! periodic boundary conditions
           if (sepy > p5lboxy) then
              sepy = sepy - lboxy
           else if (sepy < -p5lboxy) then
              sepy = sepy + lboxy
           end if
           !sepy = sepy - lboxy * nint(sepy/lboxy)         
           if (abs(sepy) < rc) then
              sepz = zpos(ipar) - zpos(jpar)
              if (zperiodic) then
                 ! periodic boundary conditions
                 if (sepz > p5lboxz) then
                    sepz = sepz - lboxz
                 else if (sepz < -p5lboxz) then
                    sepz = sepz + lboxz
                 end if
                 !sepz = sepz - lboxz*nint(sepz/lboxz)
              end if
              sepsq = sepx**2 + sepy**2 + sepz**2
              if (sepsq < rcsq) then
                 rexponenti = 1.0_db / (sepsq * (sepsq**halfexpminus2))
                 ! add contribution to total potential energy
                 if (ipar > nsurf .and. jpar > nsurf) then
                    ! both particles are fluid particles
                    epot = epot + rexponenti - vrc
                 else
                    ! at least one particle is a surface particle
                    ! currently we are doing the same thing, but may
                    ! want to add a different potential later
                    epot = epot + rexponenti - vrc2
                 end if
              end if
           end if
        end if

     end do
  end do

end subroutine ipl_totalenergy

subroutine ipl_energyipar(ipar, xposi, yposi, zposi, xpos, ypos, zpos, rc, rcsq,&
                          lboxx, lboxy, lboxz, vrc, vrc2, npar,&
                          nsurf, zperiodic, potexponent, epot)
   !!! Computes potential energy of particle ipar (in units of 4epsilon)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, intent(in) :: ipar, npar, nsurf
  real(kind=db), intent(in) :: xposi, yposi, zposi, rc, rcsq, lboxx, lboxy, lboxz
  real(kind=db), intent(in) :: vrc, vrc2, potexponent
  logical, intent(in) :: zperiodic
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos

  ! outputs
  real(kind=db), intent(out) :: epot

  !f2py intent(in) :: ipar, xposi, yposi, zposi, xpos, ypos, zpos, rc, rcsq, lboxx,
  !f2py intent(in) :: lboxy, lboxz, vrc, vrc2, npar, nsurf, zperiodic, potexponent
  !f2py intent(out) :: epot

  integer :: jpar
  real(kind=db) :: sepx, sepy, sepz, sepsq, rexponenti, halfexpminus2
  real(kind=db) :: p5lboxx, p5lboxy, p5lboxz

  halfexpminus2 = (potexponent - 2.0_db) / 2.0_db  
  p5lboxx = 0.5_db * lboxx ! for periodic bcs
  p5lboxy = 0.5_db * lboxy
  p5lboxz = 0.5_db * lboxz
  epot = 0.0_db
  do jpar = 1, npar
     if (ipar /= jpar) then
        sepx = xposi - xpos(jpar)
        ! periodic boundary conditions
        if (sepx > p5lboxx) then
           sepx = sepx - lboxx
        else if (sepx < -p5lboxx) then
           sepx = sepx + lboxx
        end if
        !sepx = sepx - lboxx * nint(sepx / lboxx)
        if (abs(sepx) < rc) then
           sepy = yposi - ypos(jpar)
           ! periodic boundary conditions
           if (sepy > p5lboxy) then
              sepy = sepy - lboxy
           else if (sepy < -p5lboxy) then
              sepy = sepy + lboxy
           end if
           !sepy = sepy - lboxy * nint(sepy / lboxy)         
           if (abs(sepy) < rc) then
              sepz = zposi - zpos(jpar)
              if (zperiodic) then
                 ! periodic boundary conditions
                 if (sepz > p5lboxz) then
                    sepz = sepz - lboxz
                 else if (sepz < -p5lboxz) then
                    sepz = sepz + lboxz
                 end if
                 !sepz = sepz - lboxz * nint(sepz / lboxz)
              end if
              sepsq = sepx**2 + sepy**2 + sepz**2
              if (sepsq < rcsq) then
                 rexponenti = 1.0_db / (sepsq * (sepsq**halfexpminus2))
                 ! add contribution to total potential energy
                 if (ipar > nsurf .and. jpar > nsurf) then
                    ! both particles are fluid particles
                    epot = epot + rexponenti - vrc
                 else
                    ! at least one particle is a surface particle
                    epot = epot + rexponenti - vrc2
                 end if
              end if
           end if
        end if
     end if
  end do

end subroutine ipl_energyipar

subroutine ipl_eij(ipar, jpar, xposi, yposi, zposi, xposj, yposj, zposj,&
                   lboxx, lboxy, lboxz, rc, rcsq, vrc, vrc2, npar, nsurf,&
                   zperiodic, potexponent, eij)
  !!! Compute potential energy between particles i and j, taking into
  !!! account periodic bcs.
  
  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)
  
  ! inputs
  integer, intent(in) :: ipar, jpar
  real(kind=db), intent(in) :: xposi, yposi, zposi, xposj, yposj, zposj
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, rc, rcsq, vrc, vrc2
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic
  real(kind=db), intent(in) :: potexponent

  ! outputs
  real(kind=db), intent(out) :: eij

  !f2py intent(in) :: ipar, jpar, xposi, yposi, zposi, xposj, yposj, zposj
  !f2py intent(in) :: lboxx, lboxy, lboxz, vrc, vrc2, npar, nsurf, zperiodic
  !f2py intent(in) :: potexponent
  !f2py intent(out) :: eij

  real(kind=db) :: sepx, sepy, sepz, sepsq, rexponenti, halfexpminus2
  halfexpminus2 = (potexponent - 2.0_db) / 2.0_db  

  eij = 0.0_db
  sepx = xposi - xposj
  ! periodic boundary conditions
  if (sepx > 0.5 * lboxx) then
     sepx = sepx - lboxx
  else if (sepx < -0.5 * lboxx) then
     sepx = sepx + lboxx
  end if

  if (abs(sepx) < rc) then
     sepy = yposi - yposj
     ! periodic boundary conditions
     if (sepy > 0.5 * lboxy) then
        sepy = sepy - lboxy
     else if (sepy < -0.5 * lboxy) then
        sepy = sepy + lboxy
     end if

     if (abs(sepy) < rc) then
        sepz = zposi - zposj
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
           rexponenti = 1.0_db / (sepsq * (sepsq**halfexpminus2))
           ! add contribution to total potential energy
           if (ipar > nsurf .and. jpar > nsurf) then
              ! both particles are fluid particles
              eij = rexponenti - vrc
           else
              ! at least one particle is a surface particle
              eij = rexponenti - vrc2
           end if
        end if
     end if
  end if
end subroutine ipl_eij

subroutine ipl_totalencreatelist(xpos, ypos, zpos, rc, rcsq,&
                                 lboxx, lboxy, lboxz, vrc, vrc2,&
                                 npar, nsurf, zperiodic, potexponent,&
                                 etot)
  ! Create the cell list and then return total energy

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz, vrc, vrc2
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic
  real(kind=db), intent(in) :: potexponent

  ! outputs
  real(kind=db), intent(out) :: etot

  !f2py intent(in) :: xpos, ypos, zpos, rc, rcsq, lboxx, lboxy, lboxz
  !f2py intent(in) :: vrc, vrc2, npar, nsurf, zperiodic, potexponent
  !f2py intent(out) :: etot

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
  call ipl_totalenlist(ll, hoc, ncelx, ncely, ncelz, rnx, rny, rnz,&
                       xpos, ypos, zpos, rc, rcsq, lboxx, lboxy,&
                       lboxz, vrc, vrc2, npar, nsurf, zperiodic,&
                       potexponent, etot)
  ! deallocate the head of cell array
  deallocate(hoc, STAT=status)

end subroutine ipl_totalencreatelist

subroutine ipl_totalenlist(ll, hoc, ncelx, ncely, ncelz,&
                           rnx, rny, rnz,&
                           xpos, ypos, zpos, rc, rcsq,&
                           lboxx, lboxy, lboxz, vrc, vrc2, npar, nsurf,&
                           zperiodic, potexponent, etot)
  ! Compute total potential energy of system using cell lists

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, dimension(npar), intent(in) :: ll
  integer, intent(in) :: ncelx, ncely, ncelz
  real(kind=db), intent(in) :: rnx, rny, rnz
  integer, dimension(ncelx, ncely, ncelz), intent(in) :: hoc
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz, vrc, vrc2
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic
  real(kind=db), intent(in) :: potexponent

  ! outputs
  real(kind=db), intent(out) :: etot

  real(kind=db) :: eij, xposi, yposi, zposi
  integer :: icelx, icely, icelz, cellnum, jpar, celx, cely, celz,&
             nceltot, ipar
  etot = 0.0_db

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
     ! cell.
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

              ! get p.e. between particles i and j
              call ipl_eij(ipar, jpar, xposi, yposi, zposi, xpos(jpar),&
                           ypos(jpar), zpos(jpar), lboxx, lboxy, lboxz,&
                           rc, rcsq, vrc, vrc2, npar, nsurf, zperiodic,&
                           potexponent, eij)           
              etot = etot + eij
           end if
           jpar = ll(jpar)

        end do
     end do
  end do

  ! double counting of potential energy
  etot = etot / 2.0_db

end subroutine ipl_totalenlist

subroutine ipl_enlist(ll, hoc, ncelx, ncely, ncelz, ipar, xposi,&
                      yposi, zposi, xpos, ypos, zpos, rc, rcsq,&
                      lboxx, lboxy, lboxz, vrc, vrc2, npar, nsurf,&
                      zperiodic, potexponent, epot)
  ! Compute potential energy of particle i using cell lists

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, dimension(npar), intent(in) :: ll
  integer, intent(in) :: ncelx, ncely, ncelz  
  integer, dimension(ncelx, ncely, ncelz), intent(in) :: hoc
  integer, intent(in) :: ipar
  real(kind=db), intent(in) :: xposi, yposi, zposi
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz, vrc, vrc2
  integer, intent(in) :: npar, nsurf
  logical, intent(in) :: zperiodic
  real(kind=db), intent(in) :: potexponent

  ! outputs
  real(kind=db), intent(out) :: epot

  real(kind=db) :: rnx, rny, rnz, eij
  integer :: icelx, icely, icelz, cellnum, jpar, celx, cely, celz, nceltot
  epot = 0.0_db

  ! cell dimension in x, y and z directions
  rnx = lboxx / ncelx
  rny = lboxy / ncely
  rnz = lboxz / ncelz

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
           
           ! get p.e. between particles i and j
           call ipl_eij(ipar, jpar, xposi, yposi, zposi, xpos(jpar),&
                        ypos(jpar), zpos(jpar), lboxx, lboxy, lboxz,&
                        rc, rcsq, vrc, vrc2, npar, nsurf, zperiodic,&
                        potexponent, eij)
           epot = epot + eij
        end if
        jpar = ll(jpar)

     end do
  enddo

end subroutine ipl_enlist
