! gauss_mccylefnpt.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutine for executing Monte carlo cycles.  A cycle
! consists of nparfl attempted positional moves and a single volume
! move (on average).  The Metropolis Monte Carlo algorithm is used.
!
! SUBROUTINES:
! gauss_executecyclesnpt - execute ncycles monte carlo cycles
!                          note xpos,ypos,zpos,lboxx,lboxy,lboxz and
!                          etot are returned. 

subroutine gauss_executecyclesnpt(xpos,ypos,zpos,ncycles,nsamp,rc,rcsq,&
                                  vrc,vrc2,press,lboxx,lboxy,lboxz,&
                                  epsovert,maxdisp,maxvol,npar,&
                                  nsurf,zperiodic,etot)
  ! execute ncycles MC cycles

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: ncycles,nsamp,npar,nsurf
  real(kind=db), intent(in) :: rc,rcsq,vrc,vrc2,press
  real(kind=db), intent(in) :: epsovert,maxdisp,maxvol
  logical, intent(in) :: zperiodic
  ! outputs (note inout intent)
  real(kind=db), dimension(npar), intent(inout) :: xpos,ypos,zpos
  real(kind=db), intent(inout) :: lboxx,lboxy,lboxz,etot

  !f2py intent(in) :: ncycles,nsamp,rc,rcsq,vrc,vrc2
  !f2py intent(in) :: epsovert, maxdisp,npar,nparsuf,zperiodic
  !f2py intent(in,out) :: xpos,ypos,zpos,lboxx,lboxy,lboxz,etot

  integer :: ipar,atmovdisp,acmovdisp,atmovvol,acmovvol,cy,it,nparfl,i,j
  real(kind=db) :: rsc,xposi,yposi,zposi,xposinew,yposinew,zposinew,eold,enew
  real(kind=db) :: lboxnew, lboxxold, lboxyold, lboxzold
  real(kind=db) :: vboxold, lnvold, lnvnew, vboxnew
  real(kind=db) :: scalefacx, scalefacy, scalefacz, arg, etotnew
  real(kind=db), dimension(3) :: rvec
  logical :: accept
  
  ! initialize random number generator
  call init_random_seed()

  atmovdisp = 0
  acmovdisp = 0
  atmovvol = 0
  acmovvol = 0
  nparfl = npar - nsurf

  write(*,'(F12.6, F12.6)') etot,lboxx
  do cy=1,ncycles
     ! each cycle is on average 1 move per fluid par + 1 vol move     
     do it=1,nparfl+1 

        ! pick a random number between [nsurf,ntot+1]
        call random_number(rsc)
        ipar = int(rsc*(nparfl + 1)) + nsurf + 1

        ! if ipar > npar, we attempt a volume move, else we attempt a
        ! positional move

        if (ipar > npar) then ! volume move
           atmovvol = atmovvol + 1

           ! old box volume
           lboxxold = lboxx
           lboxyold = lboxy
           lboxzold = lboxz
           vboxold = lboxx*lboxy*lboxz
           lnvold = log(vboxold)

           ! random number between 0 and 1 for attempted volume move
           call random_number(rsc)
              
           ! new box volume
           lnvnew = lnvold + maxvol*(rsc - 0.5_db)
           vboxnew = exp(lnvnew)

           ! scale factor for multiplying each dimension of simulation
           ! box. If z is periodic we multiply each of the three box
           ! dimensions by the same scale factor.  Otherwise, if z is
           ! not periodic, we make the box larger/smaller in the z
           ! direction only.
           
           if (zperiodic) then
              scalefacx = (vboxnew / vboxold) ** (1.0_db/3.0_db)
              scalefacy = scalefacx
              scalefacz = scalefacx
           else
              scalefacx = 1.0_db
              scalefacy = 1.0_db
              scalefacz = (vboxnew / vboxold)
           end if

           ! new box dimensions
           lboxx = lboxx*scalefacx
           lboxy = lboxy*scalefacy
           lboxz = lboxz*scalefacz

           ! rescale particle positions to new volume
           xpos = xpos*scalefacx
           ypos = ypos*scalefacy
           zpos = zpos*scalefacz

           call gauss_totalenergy(xpos,ypos,zpos,rc,rcsq,lboxx,lboxy,&
                                  lboxz,vrc,vrc2,npar,nsurf,zperiodic,&
                                  etotnew)

           ! See FS p122 (Algorithm 11) for this acceptance rule
           ! ignore - (Note that the factor in front of lnvnew - lnvold
           ! is 1 and not nparfl + 1.  This is since we are not using
           ! rescaled position coordinates and so the integral in
           ! FS Eq. (5.4.12) has a factor of V rather than V^{N+1})
           
           arg = epsovert*(etot - etotnew + press*(vboxold - vboxnew)) + &
                 (nparfl + 1)*(lnvnew - lnvold)
           accept = .True.
           if (arg < 0) then
              call random_number(rsc)
              if (rsc > exp(arg)) then
                 accept = .False.
                 ! back to old boxsize and positions
                 lboxx = lboxxold
                 lboxy = lboxyold
                 lboxz = lboxzold
                 xpos = xpos / scalefacx
                 ypos = ypos / scalefacy
                 zpos = zpos / scalefacz
              end if
           end if

           if (accept) then
              etot = etotnew
              acmovvol = acmovvol + 1
           end if

        else ! displacement move
           atmovdisp = atmovdisp + 1

           xposi = xpos(ipar)
           yposi = ypos(ipar)
           zposi = zpos(ipar)

           ! find old energy
           call gauss_energyipar(ipar,xposi,yposi,zposi,xpos,ypos,zpos,rc,&
                                 rcsq,lboxx,lboxy,lboxz,vrc,vrc2,npar,&
                                 nsurf,zperiodic,eold)

           ! displace particle
           call random_number(rvec)
           xposinew = xposi + maxdisp*(rvec(1) - 0.5_db)
           yposinew = yposi + maxdisp*(rvec(2) - 0.5_db)
           zposinew = zposi + maxdisp*(rvec(3) - 0.5_db)

           if (zperiodic) then
              if (zposinew < 0.0_db) then
                 zposinew = zposinew + lboxz
              else if (zposinew > lboxz) then
                 zposinew = zposinew - lboxz
              end if
           end if

           ! if z not periodic, we will reject
           ! the move if z > lboxz (hard wall boundary)
           if (zposinew < lboxz .and. zposinew > 0.0_db) then

              ! periodic boundary conditions in x and y
              if (xposinew < 0.0_db) then
                 xposinew = xposinew + lboxx
              else if (xposinew > lboxx) then
                 xposinew = xposinew - lboxx
              end if
              if (yposinew < 0.0_db) then
                 yposinew = yposinew + lboxy
              else if (yposinew > lboxy) then
                 yposinew = yposinew - lboxy
              end if

              ! find new energy
              call gauss_energyipar(ipar,xposinew,yposinew,zposinew,xpos,&
                                    ypos,zpos,rc,rcsq,lboxx,lboxy,lboxz,&
                                    vrc,vrc2,npar,nsurf,zperiodic,enew)

              ! choose whether to accept the move or not
              accept = .True.
              if (enew > eold) then
                 call random_number(rsc)
                 if (exp((eold - enew)*epsovert) < rsc) accept = .False.
              end if

              ! update positions if move accepted
              if (accept) then
                 xpos(ipar) = xposinew
                 ypos(ipar) = yposinew
                 zpos(ipar) = zposinew
                 etot = etot - eold + enew
                 acmovdisp = acmovdisp + 1
              end if
           
           end if
        end if
     end do
     
     ! write out energy after every nsamp cycles
     if (mod(cy,nsamp) == 0) write(*,'(I7,F12.6, F12.6, F12.6, F12.6)') cy,etot,lboxx,lboxy,lboxz
     
  end do

  ! write out acceptance ratio
  write(*,'("acceptance ratio", I7, I7, F7.3, I7, I7, F7.3)')&
       acmovdisp,atmovdisp,real(acmovdisp)/atmovdisp,&
       acmovvol,atmovvol,real(acmovvol)/atmovvol

end subroutine gauss_executecyclesnpt
