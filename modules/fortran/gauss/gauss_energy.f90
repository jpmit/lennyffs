! gauss_energy.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutines for computing potential energy (p.e.) for
! Gaussian interaction.
!
! SUBROUTINES:
! gauss_totalenergy - compute total p.e. of system by summing
!                     p.e. between all pairs of particles
! gauss_energyipar  - compute p.e. between a given particle
!                     and all other particles

subroutine gauss_totalenergy(xpos, ypos, zpos, rc, rcsq,&
                             lboxx, lboxy, lboxz, vrc, vrc2, npar,&
                             nsurf, zperiodic, epot)
   !!! Computes total potential energy of system (in units of epsilon)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: npar, nsurf
  real(kind=db), intent(in) :: rc, rcsq, lboxx, lboxy, lboxz
  real(kind=db), intent(in) :: vrc,vrc2
  logical, intent(in) :: zperiodic
  real(kind=db), dimension(npar), intent(in) :: xpos,ypos,zpos
  ! outputs
  real(kind=db), intent(out) :: epot

  !f2py intent(in) :: xpos,ypos,zpos,rc,rcsq,lboxx,lboxy,lboxz
  !f2py intent(in) :: vrc,vrc2,npar,nsurf,zperiodic
  !f2py intent(out) :: epot

  integer :: ipar,jpar
  real(kind=db) :: sepx,sepy,sepz,sepsq,r2i,r6i,r12i
  real(kind=db) :: p5lboxx, p5lboxy, p5lboxz

  p5lboxx = 0.5_db*lboxx ! for periodic bcs
  p5lboxy = 0.5_db*lboxy
  p5lboxz = 0.5_db*lboxz
  epot = 0.0_db
  do ipar=1,npar-1
     do jpar=ipar+1,npar
        sepx = xpos(ipar) - xpos(jpar)
        ! periodic boundary conditions
        if (sepx > p5lboxx) then
           sepx = sepx - lboxx
        else if (sepx < -p5lboxx) then
           sepx = sepx + lboxx
        end if
        !sepx = sepx - lboxx*nint(sepx/lboxx)
        if (abs(sepx) < rc) then
           sepy = ypos(ipar) - ypos(jpar)
           ! periodic boundary conditions
           if (sepy > p5lboxy) then
              sepy = sepy - lboxy
           else if (sepy < -p5lboxy) then
              sepy = sepy + lboxy
           end if
           !sepy = sepy - lboxy*nint(sepy/lboxy)         
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
                 ! add contribution to total potential energy
                 if (ipar > nsurf .and. jpar > nsurf) then
                    ! both particles are fluid particles
                    epot = epot + exp(-sepsq) - vrc
                 else
                    ! at least one particle is a surface particle
                    ! currently we are doing the same thing, but may
                    ! want to add a different potential later
                    epot = epot + exp(-sepsq) - vrc2                    
                 end if
              end if
           end if
        end if

     end do
  end do

end subroutine gauss_totalenergy

subroutine gauss_energyipar(ipar,xposi,yposi,zposi,xpos,ypos,zpos,rc,rcsq,&
                            lboxx,lboxy,lboxz,vrc,vrc2,npar,&
                            nsurf,zperiodic,epot)
   !!! Computes potential energy of particle ipar (in units of epsilon)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: ipar,npar,nsurf
  real(kind=db), intent(in) :: xposi, yposi, zposi, rc, rcsq, lboxx, lboxy, lboxz
  real(kind=db), intent(in) :: vrc,vrc2
  logical, intent(in) :: zperiodic
  real(kind=db), dimension(npar), intent(in) :: xpos,ypos,zpos
  ! outputs
  real(kind=db), intent(out) :: epot

  !f2py intent(in) :: ipar,xposi,yposi,zposi,xpos,ypos,zpos,rc,rcsq,lboxx,
  !f2py intent(in) :: lboxy,lboxz,vrc,vrc2,npar,nsurf,zperiodic
  !f2py intent(out) :: epot

  integer :: jpar
  real(kind=db) :: sepx,sepy,sepz,sepsq,r2i,r6i,r12i
  real(kind=db) :: p5lboxx, p5lboxy, p5lboxz

  p5lboxx = 0.5_db*lboxx ! for periodic bcs
  p5lboxy = 0.5_db*lboxy
  p5lboxz = 0.5_db*lboxz
  epot = 0.0_db
  do jpar=1,npar
     if (ipar /= jpar) then
        sepx = xposi - xpos(jpar)
        ! periodic boundary conditions
        if (sepx > p5lboxx) then
           sepx = sepx - lboxx
        else if (sepx < -p5lboxx) then
           sepx = sepx + lboxx
        end if
        if (abs(sepx) < rc) then
           sepy = yposi - ypos(jpar)
           ! periodic boundary conditions
           if (sepy > p5lboxy) then
              sepy = sepy - lboxy
           else if (sepy < -p5lboxy) then
              sepy = sepy + lboxy
           end if
           if (abs(sepy) < rc) then
              sepz = zposi - zpos(jpar)
              if (zperiodic) then
                 ! periodic boundary conditions
                 if (sepz > p5lboxz) then
                    sepz = sepz - lboxz
                 else if (sepz < -p5lboxz) then
                    sepz = sepz + lboxz
                 end if
              end if
              sepsq = sepx**2 + sepy**2 + sepz**2
              if (sepsq < rcsq) then
                 ! add contribution to total potential energy
                 if (ipar > nsurf .and. jpar > nsurf) then
                    ! both particles are fluid particles
                    epot = epot + exp(-sepsq) - vrc
                 else
                    ! at least one particle is a surface particle
                    ! currently we are doing the same thing, but may
                    ! want to add a different potential later
                    epot = epot + exp(-sepsq) - vrc2
                 end if
              end if
           end if
        end if
     end if
  end do

end subroutine gauss_energyipar

subroutine gauss_eij(ipar, jpar, xposi, yposi, zposi, xposj, yposj, zposj,&
                     lboxx, lboxy, lboxz, rc, rcsq, vrc, vrc2, npar, nsurf,&
                     zperiodic, eij)
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
  ! outputs
  real(kind=db), intent(out) :: eij

  !f2py intent(in) :: ipar, jpar, xposi, yposi, zposi, xposj, yposj, zposj
  !f2py intent(in) :: lboxx, lboxy, lboxz, vrc, vrc2, npar, nsurf, zperiodic
  !f2py intent(out) :: eij

  real(kind=db) :: sepx, sepy, sepz, sepsq

  sepx = xposi - xposj
  ! periodic boundary conditions
  if (sepx > 0.5*lboxx) then
     sepx = sepx - lboxx
  else if (sepx < -0.5*lboxx) then
     sepx = sepx + lboxx
  end if

  if (abs(sepx) < rc) then
     sepy = yposi - yposj
     ! periodic boundary conditions
     if (sepy > 0.5*lboxy) then
        sepy = sepy - lboxy
     else if (sepy < -0.5*lboxy) then
        sepy = sepy + lboxy
     end if

     if (abs(sepy) < rc) then
        sepz = zposi - zposj
        if (zperiodic) then
           ! periodic boundary conditions
           if (sepz > 0.5*lboxz) then
              sepz = sepz - 0.5*lboxz
           else if (sepz < -0.5*lboxz) then
              sepz = sepz + lboxz
           end if
        end if
        sepsq = sepx**2 + sepy**2 + sepz**2
        if (sepsq < rcsq) then
           ! add contribution to total potential energy
           if (ipar > nsurf .and. jpar > nsurf) then
              ! both particles are fluid particles
              eij =  exp(-sepsq) - vrc
           else
              ! at least one particle is a surface particle
              ! currently we are doing the same thing, but may
              ! want to add a different potential later
              eij = exp(-sepsq) - vrc2
           end if
        end if
     end if
  end if
end subroutine gauss_eij

subroutine gauss_enlist(ll, hoc, ncelx, ncely, ncelz, ipar, xposi,&
                        yposi, zposi, xpos, ypos, zpos, rc, rcsq,&
                        lboxx, lboxy, lboxz, vrc, vrc2, npar, nsurf,&
                        zperiodic, epot)
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

  ! outputs
  real(kind=db), intent(out) :: epot

  real(kind=db) :: rnx, rny, rnz, eij
  integer :: icelx, icely, icelz, cellnum, i, j, celx, cely, celz
  epot = 0

  ! cell dimension in x, y and z directions
  rnx = lboxx / ncelx
  rny = lboxy / ncely
  rnz = lboxz / ncelz

  ! determine cell that particle i is in
  icelx = int(xposi / rnx)
  icely = int(yposi / rny)
  icelz = int(zposi / rnz)

  ! go through each cell in turn (27 in total in three dimensions),
  ! and add pot energy between particle i and all particles in the
  ! cell
  do cellnum = 1, 27
     ! get the next cell indexes (ncelx, ncely, ncelz)
     call cellindx(cellnum, icelx, icely, icelz,& ! cell of particle i
                   ncelx, ncely, ncelz,&          ! total num cells in each dim
                   celx, cely, celz)             ! cell index we want

     j = hoc(celx, cely, celz)

     do while (j /= 0)
        if (i /= j) then
           ! get p.e. between particles i and j
           call gauss_eij(i, j, xposi, yposi, zposi, xpos(j), ypos(j), zpos(j),&
                          lboxx, lboxy, lboxz, rc, rcsq, vrc, vrc2, npar, nsurf,&
                          zperiodic, eij)           
           epot = epot + eij
        end if
        j = ll(j)
     end do
  enddo

end subroutine gauss_enlist
