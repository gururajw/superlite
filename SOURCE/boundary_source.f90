!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
subroutine boundary_source

  use randommod
  use particlemod
  use transportmod
  use sourcemod
  use physconstmod
  use groupmod
  use gridmod
  use totalsmod
  use inputparmod
  use miscmod
  implicit none

  integer :: ipart, ivac, ig, iig, i,j,k,l
  real*8 :: r1, r2, mu0, esurfpart, wl0, om0
  real*8 :: denom2, wl1, wl2
  real*8 :: srftemp = 1d4
  type(packet) :: ptcl
  real*8, dimension(grp_ng) :: emitsurfprobg  !surface emission probabilities

  esurfpart = tot_esurf/dble(src_nsurftot)

  select case(grd_igeom)
!-- 1D spherical
  case(11)
    i = 1 ! inner boundary
    j = 1
    k = 1
    l = grd_icell(i,j,k)
  case default
    stop 'boundary_source: invalid igeom'
  endselect

  if(in_srcmax>0d0) then
    srftemp = in_srcmax
  else
    srftemp = grd_radtemp(1)
  endif

  do ig = 1, grp_ng
    wl1 = pc_h*pc_c*grp_wlinv(ig+1)/(pc_kb*srftemp)
    wl2 = pc_h*pc_c*grp_wlinv(ig)/(pc_kb*srftemp)
    emitsurfprobg(ig) = 15d0*specint(wl1,wl2,3)/pc_pi**4
  enddo


!-- instantiating surface particles:
  do ipart=1,src_nsurf

!-- filling vacant spot in vacancy array
    ivac = src_ivacant(ipart)
    prt_isvacant(ivac) = .false.
!

!-- calculating wavelength
    denom2 = 0d0
    call rnd_r(r1,rnd_state)
    do ig = 1, grp_ng
      iig = ig
      if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
      denom2 = denom2+emitsurfprobg(ig)
    enddo

    call rnd_r(r1,rnd_state)
    wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))

!-- sampling surface projection
    if(in_surfsrcmu=='beam') then
      mu0 = 1d0
    elseif(in_surfsrcmu=='isot') then
      call rnd_r(r1,rnd_state)
      call rnd_r(r2,rnd_state)
      mu0 = max(r1,r2)
    endif
!
!-- selecting geometry and surface
    select case(grd_igeom)
!-- 1D (inner surface)
    case(11)
!-- calculating position
      ptcl%x = grd_xarr(i)
      ptcl%y = grd_yarr(1)
      ptcl%z = grd_zarr(1)
    endselect
!-- outward direction
    ptcl%mu = mu0

!-- particle properties are saved in comoving frame
    ptcl%e = esurfpart
    ptcl%e0 = esurfpart
    ptcl%wl = wl0

!-- save particle result
!-----------------------
     prt_particles(ivac) = ptcl


  enddo


end subroutine boundary_source
! vim: fdm=marker
