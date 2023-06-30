!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine sourceenergy

  use gridmod
  use sourcemod
  use totalsmod
  use physconstmod
  use inputparmod
  use groupmod
  implicit none
!##################################################
!This subroutine computes the distribution of source particles
!A fraction of the source particle number src_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################
  integer :: i
  ! real*8,dimension(grd_ncell) :: P_adia
  ! real*8,dimension(grp_ng) :: wlm,dwl

! Calculating fictitious emission energy per cell
!-- thermal source
  grd_emit = 0.d0
!-- surface source
  grd_emit = pc_pi*pc_c*grd_rcell**2*pc_acoef*grd_radtemp**4
!-- volume source (experimental)
  ! grd_emit = grd_vol*pc_acoef*grd_radtemp**4
!-- luminosity source (experimental)
  ! grd_emit = grd_lum !CAUTION:grd_emit can be >=0, will result in sampling problems
!-- emissivity source (experimental)
  ! do i=1,grp_ng
  !   dwl(i) = grp_wl(i+1) - grp_wl(i)
  !   wlm(i) = 0.5d0*(grp_wl(i+1) + grp_wl(i))
  ! enddo !i
  ! do i=1,grd_ncell-1
  !  grd_emit(i) = 4*pc_pi*grd_vol(i)*sum(grd_emiss(:,i)*dwl*pc_c/wlm**2)
!--  Fick's law (experimental)
  ! do i=1,grd_ncell-1
  !   grd_emit(i) = (-1.d0*pc_pi16*grd_rcell(i)**2/(3.d0*grd_capross(i)))* &
  !                 pc_acoef*pc_c*grd_radtemp(i)**3* &
  !                 (grd_radtemp(i+1)-grd_radtemp(i)) / &
  !                 (grd_rcell(i+1)-grd_rcell(i))
  ! enddo
!-- plus adiabatic compression/expansion source (experimental)
  ! P_adia = grd_natom * pc_kb / grd_tempinv
  ! grd_emit = grd_emit - P_adia*grd_divv

!-- renormalize source energy to the user input bolometric luminosity
  grd_emit = grd_emit / sum(grd_emit) * in_L_bol
  tot_sthermal = sum(grd_emit)

!-- inner boundary surface luminosity source (experimental)
  ! if(in_L_inner>0d0) tot_esurf = in_L_inner

end subroutine sourceenergy
! vim: fdm=marker
