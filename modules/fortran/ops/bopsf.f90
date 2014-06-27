! bopsf.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutines and functions for computing local bond order parameters
! See ten Wolde,Ruiz-Montero and Frenkel, Faraday Discuss. 104, 93  
!
! SUBROUTINES:
! xpars - compute xps, array of length npar, which contains indices of xtal
!         particles as identified by local bond order parameters.  Also return
!         nxtal, the number of xtal particles.  xps(1:nxtal) contains indexes
!         (into xpos,ypos,zpos) of xtal particles.
! ylm   - compute spherical harmonic y_lm
!
! FUNCTIONS:
! plm   - return legendre polynomial p_lm, note this works for positive m only
! ifac  - factorial function

subroutine xpars(xpos, ypos, zpos, npar, nsurf, lboxx, lboxy, lboxz,&
                 zperiodic, ncut, nlinks, linkval, xps, nxtal)
  !!! xpars returns array of indices of crystal particles xps,
  !!! as well as nxtal, the number of xtal particles.
  !!! this uses local bond order parameters.
  !!! see ten Wolde,Ruiz-Montero and Frenkel, Faraday Discuss. 104, 93  
  
  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)
  real(kind=db), parameter :: PI = 3.1415926535897932  

  ! inputs
  integer, intent(in) :: npar, nsurf
  real(kind=db), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=db), intent(in) :: linkval ! threshold for crystal link
  integer, intent(in) :: nlinks ! min number of links for crystal particle
  logical, intent(in) :: zperiodic
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, ncut

  ! outputs
  integer, intent(out), dimension(npar) :: xps ! array containing xtal par
                                               ! indices
  integer, intent(out) :: nxtal ! total number of crystal particles

  !f2py intent(in) :: zperiodic, lboxx, lboxy, lboxz, ncut
  !f2py intent(out) :: xps

  real(kind=db), dimension(npar,13) :: q6r, q6i
  integer, dimension(npar) :: numneigh ! number of neighbours for each par
  integer, dimension(npar) :: amxtal ! filled with 1 if crystal, 0 otherwise
  integer, dimension(npar,60) :: lneigh ! list of neighbours
  integer :: i, j, k, l, m, nlin
  integer :: maxneigh ! maximum number of neighbours
  real(kind=db) :: sepx, sepy, sepz, sepsq, ncutsq, r, costheta, rh, phi
  real(kind=db) :: ylmr, ylmi, qnorm, lval, p5lboxx, p5lboxy, p5lboxz

  q6r = 0.0_db
  q6i = 0.0_db
  numneigh = 0
  amxtal = 0
  ncutsq = ncut**2 ! square of neighbour cutoff distance
  l = 6 ! using q6 BOP
  maxneigh = 60
  p5lboxx = 0.5_db * lboxx ! for periodic bcs
  p5lboxy = 0.5_db * lboxy
  p5lboxz = 0.5_db * lboxz

  ! go through each particle in turn, and for each neighbour
  ! add contribution to q6r and q6i
  do i = 1, npar
     do j = 1, npar
        if (i /= j) then
           sepx = xpos(i) - xpos(j)
           ! periodic boundary conditions
           if (sepx > p5lboxx) then
              sepx = sepx - lboxx
           else if (sepx < -p5lboxx) then
              sepx = sepx + lboxx
           end if
           !sepx = sepx - lboxx * nint(sepx / lboxx)
           if (abs(sepx) < ncut) then
              sepy = ypos(i) - ypos(j)
              ! periodic boundary conditions
              if (sepy > p5lboxy) then
                 sepy = sepy - lboxy
              else if (sepy < -p5lboxy) then
                 sepy = sepy + lboxy
              end if
              !sepy = sepy - lboxy * nint(sepy / lboxy)         
              if (abs(sepy) < ncut) then
                 sepz = zpos(i) - zpos(j)
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
                 if (sepsq < ncutsq) then ! i and j are neighbours
                    if (numneigh(i) < maxneigh) then
                       numneigh(i) = numneigh(i) + 1
                       lneigh(i, numneigh(i)) = j ! j is next in list

                       ! compute angles cos(theta) and phi in spherical coords
                       r = sqrt(sepsq)
                       costheta = sepz / r
                       rh = sqrt(sepx**2 + sepy**2)
                       if ((sepy == 0.0_db) .and. (sepx == 0.0_db)) then
                          ! undefined value of phi in this coordinate system
                          phi = 0.0_db
                       else if (sepy > 0.0_db) then
                          phi = acos(sepx / rh)
                       else
                          phi = 2.0_db * PI - acos(sepx / rh)
                       end if
                    
                       ! compute contribution of particle j to q6 of particle i
                       do k = 1, 2 * l + 1
                          m = -l + k - 1
                          call ylm(l, m, costheta, phi, ylmr, ylmi)
                          q6r(i, k) = q6r(i, k) + ylmr
                          q6i(i, k) = q6i(i, k) + ylmi
                       end do
                    end if
                 end if
              end if
           end if
        end if
     end do

     ! Now we have an array of Nb(i)*\bar{qlm}(i) in Amanda notation
     ! Next, compute normalized values, i.e. \tidle{qlm} in Amanda notation
     if (numneigh(i) >= 1) then
        qnorm = 0.0_db
        do k = 1, 2 * l + 1
           qnorm = qnorm + q6r(i, k)**2 + q6i(i, k)**2
        end do
        qnorm = sqrt(qnorm)
        do k = 1, 2 * l + 1
           q6r(i, k) = q6r(i, k) / qnorm
           q6i(i, k) = q6i(i, k) / qnorm
        end do
     end if
  end do
  
  ! now we have \tilde{qlm}(i) for every particle i
  ! we need to compute the dot product \tilde{qlm}(i).\tilde{qlm}(j) for
  ! each neighbour pair, and whenever this is greater than the
  ! threshold (linkval, which usually is0.65), we call this a crystal link
  nxtal = 0
  do i=nsurf + 1, npar
     ! only interested in particle if it has at least nlinks neighbours
     if (numneigh(i) >= nlinks) then
        nlin = 0
        do j = 1, numneigh(i)
           k = lneigh(i, j)
           lval = 0.0_db
           do m = 1, 2 * l + 1
              lval = lval + q6r(i, m) * q6r(k, m) + q6i(i, m) * q6i(k, m)
           end do
           if (lval > linkval) then
              ! this is a crystal link
              nlin = nlin + 1
           end if
        end do
        ! if particle i has >=nlinks crystal links,
        ! it is in a crystal environment
        if (nlin >= nlinks) then
           amxtal(i) = 1
           nxtal = nxtal + 1
        end if
     end if
  end do

  ! finally, fill up xps array, which is an output, along with nxtal
  j = 1
  do i = nsurf + 1, npar
     if (amxtal(i) == 1) then ! we have a crystal particle
        xps(j) = i ! add this particle to the array
        j = j + 1
     end if
  end do

end subroutine xpars

subroutine ylm(l, m, costheta, phi, ylmr, ylmi)
  !!! Computes spherical harmonic Y(l,m,theta)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)
  real(kind=db), parameter :: PI = 3.1415926535897932

  ! inputs
  integer, intent(in) :: l, m
  real(kind=db), intent(in) :: costheta, phi

  ! outputs
  real(kind=db), intent(out) :: ylmr, ylmi

  !f2py intent(in) :: l, m, costheta, phi
  !f2py intent(out) :: ylmr, ylmi

  integer :: absm, ifac
  real(kind=db) :: coeff, plm, plmval

  absm = abs(m)
  coeff = sqrt((2.0_db * real(l, db) + 1.0_db) * real(ifac(l - absm), db)/ &
               (4.0_db * PI * real(ifac(l + absm), db)))

  if (m < 0) then
     coeff = coeff * (-1)**m
  end if

  plmval = plm(l,absm,costheta)
  ylmr = coeff * plmval * cos(m * phi)
  ylmi = coeff * plmval * sin(m * phi)

end subroutine ylm

function plm(l, m, x)
  !!! Computes the associated legendre polynomial P_l^m(x), for positive m only

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! function inputs
  integer, intent(in) :: l, m
  real(kind=db), intent(in) :: x

  ! function output
  real(kind=db) :: plm  

  !f2py intent(in) :: l, m

  integer :: i
  real(kind=db) :: pmm, fact, somx2, pm1m, pll

  ! First compute P_m^m
  pmm = 1.0_db
  fact = 1.0_db

  if (m > 0) then
     somx2 = sqrt(1.0_db - x**2)
     do i = 1, m
        pmm = -pmm * fact * somx2
        fact = fact + 2.0_db
     end do
  end if
  if (l == m) then
     plm = pmm ! return value
  else
     ! Compute P_m+1^m
     pm1m = (2 * m + 1) * x * pmm
     if (l == m + 1) then
        plm = pm1m ! return value
     else
        ! Compute P_l^m if not already returned
        do i = m + 2, l
           pll = (x * (2 * i - 1) * pm1m - (i + m - 1) * pmm) / (i - m)
           pmm = pm1m
           pm1m = pll
        end do
        plm = pll
     end if
  end if
end function plm

function ifac(n)
  !!! Computes n!, that is, n factorial
  
  implicit none
  integer :: ifac

  integer, intent(in) :: n

  integer :: i

  if (n < 2) then
     ifac = 1
  else
     ifac = 1
     do i = 2, n
        ifac = ifac * i
     end do
  end if
end function ifac
