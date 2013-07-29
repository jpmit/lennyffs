! len_energyf.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutines for computing potential energy (p.e.) for
! Lennard Jones interaction.
!
! SUBROUTINES:
! len_totalenergy - compute total p.e. of system by summing
!                   p.e. between all pairs of particles
! len_energyipar  - compute p.e. between a given particle
!                   and all other particles

subroutine len_totalenergy(xpos,ypos,zpos,rc,rcsq,&
                           lboxx,lboxy,lboxz,vrc,vrc2,npar,&
                           nsurf,zperiodic,r6mult,r12mult,epot)
   !!! Computes total potential energy of system (in units of 4epsilon)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: npar, nsurf
  real(kind=db), intent(in) :: rc,rcsq,lboxx,lboxy,lboxz
  real(kind=db), intent(in) :: vrc,vrc2,r6mult,r12mult
  logical, intent(in) :: zperiodic
  real(kind=db), dimension(npar), intent(in) :: xpos,ypos,zpos
  ! outputs
  real(kind=db), intent(out) :: epot

  !f2py intent(in) :: xpos,ypos,zpos,rc,rcsq,lboxx,lboxy,lboxz
  !f2py intent(in) :: vrc,vrc2,npar,nsurf,zperiodic,r6mult,r12mult
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
                 r2i = 1.0_db/sepsq
                 r6i = r2i**3
                 r12i = r6i**2
                 ! add contribution to total potential energy
                 if (ipar > nsurf .and. jpar > nsurf) then
                    ! both particles are fluid particles
                    epot = epot + r12i - r6i - vrc
                 else
                    ! at least one particle is a surface particle
                    epot = epot + r12mult*r12i - r6mult*r6i - vrc2
                 end if
              end if
           end if
        end if

     end do
  end do

end subroutine len_totalenergy

subroutine len_energyipar(ipar,xposi,yposi,zposi,xpos,ypos,zpos,rc,rcsq,&
                          lboxx,lboxy,lboxz,vrc,vrc2,npar,&
                          nsurf,zperiodic,r6mult,r12mult,epot)
   !!! Computes potential energy of particle ipar (in units of 4epsilon)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: ipar,npar,nsurf
  real(kind=db), intent(in) :: xposi,yposi,zposi,rc,rcsq,lboxx,lboxy,lboxz
  real(kind=db), intent(in) :: vrc,vrc2,r6mult,r12mult
  logical, intent(in) :: zperiodic
  real(kind=db), dimension(npar), intent(in) :: xpos,ypos,zpos
  ! outputs
  real(kind=db), intent(out) :: epot

  !f2py intent(in) :: ipar,xposi,yposi,zposi,xpos,ypos,zpos,rc,rcsq,lboxx,
  !f2py intent(in) :: lboxy,lboxz,vrc,vrc2,npar,nsurf,zperiodic,r6mult,r12mult
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
        !sepx = sepx - lboxx*nint(sepx/lboxx)
        if (abs(sepx) < rc) then
           sepy = yposi - ypos(jpar)
           ! periodic boundary conditions
           if (sepy > p5lboxy) then
              sepy = sepy - lboxy
           else if (sepy < -p5lboxy) then
              sepy = sepy + lboxy
           end if
           !sepy = sepy - lboxy*nint(sepy/lboxy)         
           if (abs(sepy) < rc) then
              sepz = zposi - zpos(jpar)
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
                 r2i = 1.0_db/sepsq
                 r6i = r2i**3
                 r12i = r6i**2
                 ! add contribution to total potential energy
                 if (ipar > nsurf .and. jpar > nsurf) then
                    ! both particles are fluid particles
                    epot = epot + r12i - r6i - vrc
                 else
                    ! at least one particle is a surface particle
                    epot = epot + r12mult*r12i - r6mult*r6i - vrc2
                 end if
              end if
           end if
        end if
     end if
  end do

end subroutine len_energyipar
