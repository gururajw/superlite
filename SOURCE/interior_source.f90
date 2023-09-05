!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine interior_source!(it)

  use randommod
  use sourcemod
  use miscmod
  use groupmod
  use gridmod
  use totalsmod
  use particlemod
  use physconstmod
  use inputparmod
  use countersmod

  implicit none
  ! integer, intent(in) :: it !experimental iterative approach
!##################################################
!This subroutine instantiates new volume (cell) particle properties.
!Composed of thermal source particle loop.
!##################################################
  integer :: i,j,k,ipart,ivac,ig,ii,iig
  integer :: nhere,iimpi,nemit
  integer*8,allocatable :: nvacantall(:)
  real*8 :: pwr
  real*8 :: r1, r2, r3, uul, uur, uumax
  real*8 :: om0, mu0, x0, ep0, wl0
  real*8 :: denom2,x1,x2
!-- neighbor emit values (for source tilting)
  integer :: icnb(6)
!-- jrad for sampling !experimental iterative approach
  ! real*8,dimension(grd_ncell) :: jradsum
  ! real*8,dimension(grp_ng) :: wlm,dwl
!
  real*8 :: emitprob(grp_ng)
  real*8 :: tradinv
  type(packet),pointer :: ptcl
!-- statement functions
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)

!-- shortcut
  pwr = in_srcepwr

  x1=grp_wlinv(grp_ng+1)
  x2=grp_wlinv(1)

  allocate(nvacantall(size(src_nvacantall)))

!-- Thermal volume particle instantiation: loop
  ! jradsum = 0.d0 !experimental iterative approach
  ! do i=1,grp_ng
  !   dwl(i) = grp_wl(i+1) - grp_wl(i)
  !   wlm(i) = 0.5d0*(grp_wl(i+1) + grp_wl(i))
  ! enddo !i
  ipart = src_nsurf
  iimpi = -1
  nvacantall = src_nvacantall
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     ! jradsum(l) = pc_c*sum(grd_jrad(:,l)*dwl(:)/wlm(:)**2) !experimental iterative approach
     call sourcenumbers_roundrobin_limit(iimpi,nvacantall, &
        grd_emit(l)**pwr,grd_nvol(l),nemit,nhere)
     if(nhere<1) cycle
!-- integrate planck function over each group
     tradinv = 1d0/grd_radtemp(l)
     if(tradinv/=tradinv) stop "initial_particles: tradinv NaN"
     call specintv(tradinv,grp_ng,emitprob)
!-- use this instead of blackbody distribution for emissivity based wavelength distribution (experimental)
     ! emitprob = emitprob + (grd_emiss(:,l)*dwl*pc_c/wlm**2) &
     !              /sum(grd_emiss(:,l)*dwl*pc_c/wlm**2)
     emitprob = emitprob/sum(emitprob)

!
!-- neighbors
     icnb(1) = grd_icell(max(i-1,1),j,k)      !left neighbor
     icnb(2) = grd_icell(min(i+1,grd_nx),j,k) !right neighbor
     icnb(3) = grd_icell(i,max(j-1,1),k)      !left neighbor
     icnb(4) = grd_icell(i,min(j+1,grd_ny),k) !right neighbor
     icnb(5) = grd_icell(i,j,max(k-1,1))      !left neighbor
     icnb(6) = grd_icell(i,j,min(k+1,grd_nz)) !right neighbor
!
!
  do ii=1,nhere
     ipart = ipart + 1
     ivac = src_ivacant(ipart)
     ptcl => prt_particles(ivac)

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.

!-- calculating wavelength
     denom2 = 0d0
     call rnd_r(r1,rnd_state)
     ! if(it>1) r1 = r1*jradsum(l) !experimental iterative approach
     do ig = 1, grp_ng
        iig = ig ! This is correct
        if (r1>=denom2.and.r1<denom2+emitprob(ig)) exit
        denom2 = denom2+emitprob(ig)
     enddo
!-- use this for experimental iterative approach
     ! do ig = 1, grp_ng
     !   iig = ig
     !   if(it==1) then
     !-- Blackbody
     !     if(r1>=denom2.and.r1<denom2 + specarr(ig)) exit
     !     denom2 = denom2 + specarr(ig)
     !   else
     !-- Jrad estimator for iteration
     !     if(r1>=denom2.and.r1<denom2 + &
     !      pc_c*grd_jrad(ig,l)*dwl(ig)/wlm(ig)**2) exit
     !      denom2 = denom2 + pc_c*grd_jrad(ig,l)*dwl(ig)/wlm(ig)**2
     !     endif
     ! enddo
     call rnd_r(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))
     ptcl%wl = wl0
!-- calculating particle energy
     ep0 = grd_emit(l)/dble(nemit)
     ptcl%e = ep0
     ptcl%e0 = ep0

!-- calculating direction cosine (comoving)
     call rnd_r(r1,rnd_state)
     mu0 = sqrt(r1) !for surface source
     ! mu0 = 1d0-2d0*r1 !for volume source
     ptcl%mu = mu0
!-- sampling azimuthal angle of direction
     call rnd_r(r1,rnd_state)
     om0 = pc_pi2*r1
     ptcl%om = om0

!
!-- x position:
     select case(grd_igeom)
!
!-- 1-D spherical
     case(11)
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_r(r1,rnd_state)
           x0 = (r1*grd_xarr(i+1)**3+(1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           call rnd_r(r2,rnd_state)
        enddo
        ptcl%x = x0
     endselect
!-- y,z position
        ptcl%y = grd_yarr(1)
        ptcl%z = grd_zarr(1)

  enddo !ipart
!
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=src_nnew) stop 'interior_source: n/=nnew'

  call counterreg(ct_npcreate, src_nnew)

end subroutine interior_source
! vim: fdm=marker
