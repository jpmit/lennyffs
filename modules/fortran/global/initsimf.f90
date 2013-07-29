! initsimf.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutines for initialising particle positions
! and for initialising random number generator
!
! SUBROUTINES:
! initpositionsnosurff - initialise particles positions in box
! initrandomseed       - initialise random number generator seed

subroutine initpositionsnosurff(npar,lboxx,lboxy,lboxz,rcinitsq,xpos,ypos,zpos)

  !!! Initializes particles at random positions in cubic box of length lbox

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! subroutine arguments
  ! inputs
  integer, intent(in) :: npar
  real(kind=db), intent(in) :: lboxx,lboxy,lboxz,rcinitsq
  ! outputs
  real(kind=db), dimension(npar), intent(out) :: xpos,ypos,zpos

  !f2py intent(in) :: npar,lbox,rcinitsq
  !f2py intent(out) :: xpos,ypos,zpos

  integer :: i,j
  logical :: overlap
  real(kind=db) :: sepx,sepy,sepz,sepsq
  real(kind=db), dimension(3) :: r

  ! initialize random number generator
  call init_random_seed()

  do i=1,npar
     overlap = .TRUE. 
     do while (overlap .eqv. .TRUE.)
        overlap = .FALSE.
        call random_number(r)
        xpos(i) = lboxx*r(1)
        ypos(i) = lboxy*r(2)
        zpos(i) = lboxz*r(3)
        do j=1,i-1
           sepx = xpos(i) - xpos(j)
           sepy = ypos(i) - ypos(j)
           sepz = zpos(i) - zpos(j)
           ! periodic boundary conditions
           sepx = sepx - lboxx*nint(sepx/lboxx)
           sepy = sepy - lboxy*nint(sepy/lboxy)
           sepz = sepz - lboxz*nint(sepz/lboxz)
           sepsq = sepx**2 + sepy**2 + sepz**2
           if (sepsq < rcinitsq) then
              overlap = .TRUE.
              exit
           end if
        end do
     end do
  end do
         
end subroutine initpositionsnosurff

subroutine init_random_seed()
  !!! Initialise random number generator
  !!! This comes from the documentation for the GNU fortran compiler
  !!! See http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
  !!! The only difference is that I add the PID of the process to the seed
  !!! This is so that if starting multiple simulations at
  !!!the same instant of time, the RNG is initialised differently
  
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
          
  call random_seed(size = n)
  allocate(seed(n))
          
  call system_clock(count=clock)

  seed = clock + getpid() + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
          
  deallocate(seed)
end subroutine
