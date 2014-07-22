! initsimf.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutines for initialising particle positions and for
! initialising random number generator.
!
! SUBROUTINES:
! initpositionsnosurff - initialise particles positions in box
! initrandomseed       - initialise (seed) random number generator

subroutine initpositionsnosurff(npar, lboxx, lboxy, lboxz, rcinitsq,&
                                xpos, ypos, zpos, oregion, oxmin, oxmax,&
                                oymin, oymax, ozmin, ozmax, sameseed)

  !!! Initializes particles at random positions in cubic box of
  !!! dimensions lboxx*lboxy*lboxz

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, intent(in) :: npar
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, rcinitsq
  real(kind=db), intent(in) :: oxmin, oxmax, oymin, oymax, ozmin, ozmax  
  logical, intent(in) :: oregion, sameseed

  ! outputs
  real(kind=db), dimension(npar), intent(out) :: xpos, ypos, zpos

  !f2py intent(in) :: npar, lboxx, lboxy, lboxz, rcinitsq, sameseed
  !f2py intent(out) :: xpos, ypos, zpos

  integer :: i, j
  logical :: reject
  real(kind=db) :: sepx, sepy, sepz, sepsq
  real(kind=db), dimension(3) :: r

  ! initialize random number generator
  call init_random_seed(sameseed)

  do i = 1, npar
     reject = .TRUE. 
     do while (reject .eqv. .TRUE.)
        reject = .FALSE.
        call random_number(r)
        xpos(i) = lboxx * r(1)
        ypos(i) = lboxy * r(2)
        zpos(i) = lboxz * r(3)

        ! have we put the particle in the overlap region?
        if (oregion) then
           if ((xpos(i) > oxmin) .and. (xpos(i) < oxmax)&
                .and. (ypos(i) > oymin) .and. (ypos(i) < oymax)&
                .and. (zpos(i) > ozmin) .and. (zpos(i) < ozmax)) then
              reject = .TRUE.
              cycle ! next iteration of while loop
           end if
        end if

        ! have we put the particle too close to another?
        do j = 1, i - 1
           sepx = xpos(i) - xpos(j)
           sepy = ypos(i) - ypos(j)
           sepz = zpos(i) - zpos(j)
           ! periodic boundary conditions
           sepx = sepx - lboxx * nint(sepx / lboxx)
           sepy = sepy - lboxy * nint(sepy / lboxy)
           sepz = sepz - lboxz * nint(sepz / lboxz)
           sepsq = sepx**2 + sepy**2 + sepz**2
           if (sepsq < rcinitsq) then
              reject = .TRUE.
              exit
           end if
        end do
     end do
  end do
         
end subroutine initpositionsnosurff

subroutine init_random_seed(same)
  !!! Initialise random number generator This is a modified version
  !!! from the documentation for the GNU fortran compiler. See
  !!! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html. The
  !!! differences are:
  !!!
  !!! (1) I add the PID of the process to the
  !!! seed. This is so that if starting multiple simulations at the
  !!! same instant of time, the RNG is initialised differently.
  !!!
  !!! (2) The parameter same, when set to True, will initialise the
  !!! RNG with the same seed every time.  This is mainly for
  !!! debugging, so that we can repeat the same simulation.

  implicit none

  ! inputs
  logical, intent(in)  :: same
  
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  !f2py integer, optional, intent(in) :: s
     
  call random_seed(size = n)
  allocate(seed(n))

  if (same) then
     seed = (/ (i - 1, i = 1, n) /)
  else
     call system_clock(count=clock)
     seed = clock + getpid() + 37 * (/ (i - 1, i = 1, n) /)
  end if
     
  call random_seed(put = seed)

  deallocate(seed)
  
end subroutine
