! clist.f90
! James Mithen
! j.mithen@surrey.ac.uk
!
! Fortran subroutine for using cell lists. See Frenkel-Smit Appendix F.
!
! SUBROUTINES:

subroutine new_nlist(xpos, ypos, zpos, rc, lboxx, lboxy, lboxz, npar,&
                     ncelx, ncely, ncelz, ll, hoc, rnx, rny, rnz)

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)

  ! inputs
  integer, intent(in) :: npar
  real(kind=8), dimension(npar), intent(in) :: xpos, ypos, zpos
  real(kind=8), intent(in) :: rc, lboxx, lboxy, lboxz
  integer, intent(in) :: ncelx, ncely, ncelz

  ! outputs
  integer, dimension(npar), intent(out) :: ll
  integer, dimension(ncelx, ncely, ncelz), intent(out) :: hoc
  real(kind=8), intent(out) :: rnx, rny, rnz

  !f2py intent(in) :: npar, xpos, ypos, zpos, rc, lboxx, lboxy, lboxz
  !f2py intent(out) :: ll, hoc, rnx, rny, rnz
  
  integer :: i, j, k, icelx, icely, icelz

  do i = 1, ncelx
     do j = 1, ncely
        do k = 1, ncelz
           hoc(i, j, k) = 0
        end do
     end do
  end do

  ! cell dimension in x, y and z directions
  rnx = lboxx / ncelx
  rny = lboxy / ncely
  rnz = lboxz / ncelz

  ! go through each particle in turn and assign to cell
  do i = 1, npar
     icelx = int(xpos(i) / rnx) + 1
     icely = int(ypos(i) / rny) + 1
     icelz = int(zpos(i) / rnz) + 1

     ! 'backwards' way of building a linked list structure. At the end
     ! of this, hoc will point to a single particle, then we can get
     ! the other particles by traversing through ll until we hit the
     ! value of zero.
     ll(i) = hoc(icelx, icely, icelz)
     
     hoc(icelx, icely, icelz) = i
  end do

end subroutine new_nlist

subroutine getnumcells(lboxx, lboxy, lboxz, rc, ncelx, ncely, ncelz)
  !!! from box dimensions, get the number of cells in each dim.

  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)
  
  ! subroutine arguments
  ! inputs
  real(kind=db), intent(in) :: lboxx, lboxy, lboxz, rc

  ! outputs
  integer, intent(out) :: ncelx, ncely, ncelz

  !f2py intent(in) :: lboxx, lboxy, lboxz, rc
  !f2py intent(out) :: ncelx, ncely, ncelz

  ! build the cell list, this will fill ll, hoc, ncelx, ncely, ncelz
  ncelx = int(lboxx / rc)
  ncely = int(lboxy / rc)
  ncelz = int(lboxz / rc)
  
  ! if fewer than 3 cells in any direction, make box a single cell
  if (ncelx < 3 .or. ncely < 3 .or. ncelz < 3) then
     ncelx = 1
     ncely = 1
     ncelz = 1
  end if

end subroutine getnumcells

subroutine cellindx(celnum, icelx, icely, icelz, ncelx, ncely, ncelz,&
                    celx, cely, celz)
  !!! get cell indices celx, cely, celz.
  
  implicit none
  integer, parameter :: db = 8 !selected_real_kind(13)
  
  ! subroutine arguments
  ! inputs
  integer, intent(in) :: celnum, icelx, icely, icelz, ncelx, ncely, ncelz
  ! outputs
  integer, intent(out) :: celx, cely, celz

  !f2py intent(in) :: celnum, icelx, icely, icelz, ncelx, ncely, ncelz
  !f2py intent(out) :: celx, cely, celz
  
  integer :: bcelx, bcely, bcelz, zcelnum, dcelx, dcely, dcelz

  ! handle case of ncelx, ncely, ncelz = 1
  if (ncelx == 1) then
     celx = 1
     cely = 1
     celz = 1
  else
     
     ! base cell is the bottom left of the 27 cells
     bcelx = icelx - 1
     bcely = icely - 1
     bcelz = icelz - 1

     ! numbers to add to base cell index
     zcelnum = celnum - 1
     dcelx = mod(zcelnum, 3)
     dcely = mod(zcelnum / 3, 3)
     dcelz = mod(zcelnum / 9, 3)

     celx = bcelx + dcelx
     cely = bcely + dcely
     celz = bcelz + dcelz

     ! periodic boundary conditions
     if (celx <= 0) then
        celx = celx + ncelx
     else if (celx > ncelx) then
        celx = celx - ncelx
     end if

     if (cely <= 0) then
        cely = cely + ncely
     else if (cely > ncely) then
        cely = cely - ncely
     end if

     if (celz <= 0) then
        celz = celz + ncelz
     else if (celz > ncelz) then
        celz = celz - ncelz
     end if
  end if
end subroutine cellindx
