!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
pure subroutine transport11(ptcl,ptcl2,vx,vy,vz,rndstate,&
  eraddens,eamp,jrad,totevelo,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(inout) :: vx,vy,vz
  real*8,intent(out) :: eraddens, eamp
  real*8,intent(out) :: jrad
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c

  logical :: lout
  integer :: ixnext
  real*8 :: r1, r2
  real*8 :: db, dcol, dthm, ddop
  real*8 :: darr(4)
  real*8 :: elabfact,elabfactold
  real*8 :: xold, muold
  real*8 :: help
!-- distance out of physical reach
  real*8 :: far
  real*8 :: emitlump
!-- interpolate velocity
  real*8 :: vhelp,vhelp1,vhelp2,help1,help2
  real*8 :: dummy

  integer,pointer :: ix,ic,ig
  integer,parameter :: iy=1, iz=1
  real*8,pointer :: x, mu, e, e0, wl, d
!-- anisotropic Thomson scattering
  real*8 :: mutemp,muhelp,munew1,munew2,B,C
!-- Doppler distance calculations
  real*8 :: wlhelp!,vxhelp,labfacthelp
  ! logical :: lsuplum
  real*8 :: shelp,xhelp,epstol,phi,dphi
  integer :: it

!-- statement function
  integer :: l
  real*8 :: dx
  real*8 :: dvx, dvdr
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dvx(l) = grd_vxarr(l+1) - grd_vxarr(l)
  dvdr(l) = dvx(l) / dx(l)

  ix => ptcl2%ix
  ic => ptcl2%ic
  ig => ptcl2%ig
  d => ptcl2%dist
  x => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

  dummy = vy
  dummy = vz

!-- no error by default
  ierr = 0
!-- init
  eraddens = 0d0
  eamp = 0d0
  jrad = 0d0
!
!-- calculating initial transformation factors
  elabfact = 1d0 - mu*vx*cinv
  elabfactold = elabfact

!-- distance longer than logest possible photon path
  far = 2d0*abs(sqrt(grd_xarr(ix+1)**2-grd_xarr(ix)**2))

!
!-- boundary distances
  if(ix==1 .or. mu>=-sqrt(1d0-(grd_xarr(ix)/x)**2)) then
!-- outer boundary
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
     ixnext = ix+1
  else
!-- inner boundary
     db = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
     ixnext = ix-1
  endif
!-- sanity check
  if(db/=db) then
!    stop 'transport11: db/=db'
     ierr = 1
     return
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = -log(r1)/(elabfact*grd_sig(ic))
  else
     dthm = far
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = far
  else
     call rnd_r(r1,rndstate)
     dcol = -log(r1)/(elabfact*grd_cap(ig,ic))
  endif
!
!-- Doppler shift distance
!-- \tilde{v} frame
  vhelp = grd_vxarr(ix) - dvdr(ix)*grd_xarr(ix)
  muhelp = (mu - vx*cinv) / (1 - mu*vx*cinv)
  wlhelp = (wl/elabfact) / (1 + (muhelp*(vx-vhelp)*cinv))
  ! lsuplum = .false.
  if(abs(vx - vhelp)/pc_c.ge.0.1d0) then  ! check for superluminal flow
     ! lsuplum = .true.
     ! ddop = far
     goto 67
  ! elseif(dvdr(ix)>0.and.ig<grp_ng) then !redshift
  !    ddop = pc_c*wlhelp/(wl*dvdr(ix))*(elabfact - wl*grp_wlinv(ig+1))
  ! elseif(dvdr(ix)<0.and.ig>1) then ! blueshift
  !    ddop = pc_c*wlhelp/(wl*dvdr(ix))*(elabfact - wl*grp_wlinv(ig))
  ! else
  !    ddop = far
  ! endif
!-- SuperNu method modified for non-homologous outflows
  elseif(dvdr(ix)>0.and.ig<grp_ng) then !redshift
     ddop = 3d0*pc_c/grd_divv(ix)*(elabfact - wl*grp_wlinv(ig+1))
  elseif(dvdr(ix)<0.and.ig>1) then ! blueshift
     ddop = 3d0*pc_c/grd_divv(ix)*(elabfact - wl*grp_wlinv(ig))
  else
     ddop = far
  endif

  if(ddop<0d0) then
     ddop = far
  endif
  goto 66

!-- Newton-Raphson Method (Wagle et al. 2023, ApJ, 953, 132)
67 continue
  epstol = 1d-3 ! tolerance for convergence
  it = 0
  shelp = db ! initial value (s_0)
  !-- new x & mu
  xhelp = sqrt(x**2 + shelp**2 + 2d0*shelp*x*mu)
  muhelp = (mu*x + shelp) / xhelp
  !-- phi
  phi = dvdr(ix)*mu*x - pc_c*(1-wl*grp_wlinv(ig+1))+&
        vhelp*muhelp + dvdr(ix)*shelp
  !-- derivative of phi
  dphi = dvdr(ix) + vhelp*(1 - muhelp**2)/xhelp
  !-- next value, s_k
  ddop = shelp - phi / dphi
  !-- iterate to find ddop
  do while((abs(ddop-shelp).gt.ddop*epstol).and.it.lt.100)
    it = it + 1
    shelp = ddop ! store previous value
    !-- new x & mu
    xhelp = sqrt(x**2 + shelp**2 + 2d0*shelp*x*mu)
    muhelp = (mu*x + shelp) / xhelp
    !-- phi
    if(dvdr(ix)>0.and.ig<grp_ng) then
      phi = dvdr(ix)*mu*x - pc_c*(1-wl*grp_wlinv(ig+1))+&
            vhelp*muhelp + dvdr(ix)*shelp
    elseif(dvdr(ix)<0.and.ig>1) then
      phi = dvdr(ix)*mu*x - pc_c*(1-wl*grp_wlinv(ig))+&
            vhelp*muhelp + dvdr(ix)*shelp
    else
      ddop = far
      exit
    endif
    !-- derivative of phi
    dphi = dvdr(ix) + vhelp*(1 - muhelp**2)/xhelp
    !-- update ddop
    ddop = shelp - phi / dphi
  enddo

  if(ddop.le.0d0.or.it.eq.100) then
     ddop = far
  endif

66 continue
!
!--finding minimum distance
  darr = [dcol,dthm,ddop,db]
  ptcl2%idist = minloc(darr,dim=1)
  d = minval(darr)
  if(any(darr/=darr)) then
     ierr = 3
     return
  endif
  if(d<0d0) then
     ierr = 4
     return
  endif

!-- updating position, angle
  xold = x
  x = sqrt(x**2 + d**2 + 2d0*d*x*mu)
!-- interpolate velocity
  vhelp1 = grd_vxarr(ix)
  vhelp2 = grd_vxarr(ix+1)
  help1 = grd_xarr(ix)
  help2 = grd_xarr(ix+1)
  vx = (vhelp1*(help2-x)+vhelp2*(x-help1))/dx(ix)
  muold = mu
  if(x==0d0) then
     mu = 1d0
  else
     mu = (xold*mu+d)/x
  endif

!-- tallying energy densities
  eraddens = e*elabfact**2*d*cinv

!-- updating raditation intensity
  jrad = pc_c*eraddens
!
!-- updating transformation factors
  elabfact = 1d0 - mu*vx*cinv

!
!-- instantiate new particle at inner boundary
  if(x < grd_xarr(1)) then
     x = grd_xarr(1)
     vx = grd_vxarr(1)
     ix = 1
     ic = grd_icell(ix,iy,iz)
     call rnd_r(r1,rndstate)
     !-- forward scattering 0 <= mu <=1
     mu = 1d0 - r1
     return
  endif
!

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
     call rnd_r(r1,rndstate)
!-- anisotropic scattering (experimental)
     if(d==dthm.and.trn_thm_aniso) then
       mutemp = ((4d0*r1-2d0+SQRT(1d0+16d0*(r1-0.5)**2d0))**(1d0/3d0) - &
             (SQRT(1d0+16d0*(r1-0.5)**2d0)-(4d0-2d0*r1))**(1d0/3d0))
       muhelp = mutemp ! need to verify
       B = 2*mutemp*muhelp
       C = mutemp**2 + muhelp**2 - 1
       munew1 = (-B + SQRT(B**2 - 4*C))/2
       munew2 = (-B - SQRT(B**2 - 4*C))/2
       !-- randomizing the angle sampling
       if(munew1.le.1.and.munew1.ge.-1.and.&
         munew2.le.1.and.munew2.ge.-1) then
         if(r1<0.5) then
           mu = munew1
         else
           mu = munew2
         endif
       elseif(munew1.le.1.and.munew1.ge.-1) then
         mu = munew1
       elseif(munew2.le.1.and.munew2.ge.-1) then
         mu = munew2
       else
         mu = 1d0-2d0*r1
       endif
!-- isotropic scattering
     else
       mu = 1d0-2d0*r1
     endif
!-- checking velocity dependence
     mu=(mu+vx*cinv)/(1d0+vx*mu*cinv)
!-- outer domain boundary check
  elseif(d==db) then
     lout = mu>=0d0.and.ix==grd_nx
     if(lout) then
!-- ending particle
        ptcl2%stat = 'flux'
        return
     endif
  endif

!
!-- Thomson scatter
  if(d == dthm) then
!-- lab wavelength
     wl = wl*(1d0-mu*vx*cinv)/elabfact
     help = elabfact/(1d0-mu*vx*cinv)
!-- velocity effects accounting
     totevelo=totevelo+e*(1d0-help)
!-- energy weight
     e = e*help
     e0 = e0*help

!
!-- effective scattering
  elseif(d == dcol) then
     call rnd_r(r1,rndstate)
!-- transforming to lab
     help = elabfact/(1d0-mu*vx*cinv)
!-- velocity effects accounting
     totevelo = totevelo+e*(1d0-help)
!-- energy weight
     e = e*help
     e0 = e0*help
!-- redistributing wavelength
     if(in_nlte) then !NLTE
        emitlump = grd_opaclump(9,ic)/grd_capemitgrey(ic)
     else !LTE
        emitlump = grd_opaclump(8,ic)/grd_capgrey(ic)
     endif
     if(grp_ng==1) then
        !-- do nothing
     else
        call rnd_r(r1,rndstate)
        ig = emitgroup(r1,ic)
!-- checking if DDMC in new group
        if(((grd_sig(ic)+grd_cap(ig,ic))*dx(ix) >= trn_tauddmc) &
            .and. .not.in_puretran) ptcl2%itype = 2
     endif
!-- checking for DDMC in new group
     if(ptcl2%itype==2) then
!-- transforming to cmf
!-- velocity effects accounting
        totevelo = totevelo+e*vx*mu*cinv
!-- energy weight
        e = e*(1d0-vx*mu*cinv)
        e0 = e0*(1d0-vx*mu*cinv)
        !wl = 0d0 !workaround ifort 13.1.3 bug
     else
!-- uniformly in new group
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- converting comoving wavelength to lab frame wavelength
        wl = wl*(1d0-vx*mu*cinv)
     endif

!
!-- outer radial bound
  elseif(d==db .and. ixnext>ix) then

     l = grd_icell(ix+1,iy,iz)
     if((grd_sig(l)+grd_cap(ig,l))*dx(ix+1)<trn_tauddmc &
          .or. in_puretran) then
!-- IMC in adjacent cell
        x = grd_xarr(ix+1)
        vx = grd_vxarr(ix+1)
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        mu = (mu-vx*cinv)/(1d0-vx*mu*cinv)
        help= (grd_cap(ig,l)+grd_sig(l))*dx(ix+1)
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5*abs(mu))) then
           ptcl2%itype = 2
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-elabfact)
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
!-- update
           ix = ix+1
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
           mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
           x = grd_xarr(ix+1)
           vx = grd_vxarr(ix+1)
        endif
     endif

!
!-- inner radial bound
  elseif(d==db) then

     l = grd_icell(ix-1,iy,iz)
     if((grd_sig(l)+grd_cap(ig,l))*dx(ix-1) &
           < trn_tauddmc .or. &
          in_puretran) then
!-- IMC in adjacent cell
        x = grd_xarr(ix)
        vx = grd_vxarr(ix)
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
!-- transforming x-cosine to cmf
        mu = (mu-vx*cinv)/(1d0-vx*mu*cinv)
!-- amplification factor
        if(.not.trn_noampfact .and. mu<0d0) then
            help = 1d0/abs(mu)
            help = min(100d0, help) !-- truncate singularity
!-- velocity effects accounting
            totevelo = totevelo-e*2d0*(0.55d0*help-1.25d0*abs(mu))*vx*cinv
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*vx*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0 = e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*vx*cinv)
              e = e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*vx*cinv)
        endif
        help = (grd_cap(ig,l)+grd_sig(l))*dx(ix-1)
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl2%itype = 2
!-- velocity effects accounting
           totevelo = totevelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
           ix = ix-1
           ic = grd_icell(ix,iy,iz)
        else
!-- resampling direction
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
!-- transforming mu to lab
           mu=(mu+vx*cinv)/(1d0+vx*mu*cinv)
           x = grd_xarr(ix)
           vx = grd_vxarr(ix)
        endif
     endif

!
!-- Doppler shift
  elseif(d == ddop) then
     if(dvdr(ix).gt.0) then ! redshit
        if(ig<grp_ng) then
!-- shifting group
          ig = ig+1
          wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfactold
        else
!-- resampling wavelength in highest group
          call rnd_r(r1,rndstate)
          wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
          wl = wl*elabfact
        endif
     elseif(dvdr(ix)<0) then ! blueshift
       if(ig>1) then
!-- shifting group
          wl = (grp_wl(ig)+1d-6*(grp_wl(ig-1)-grp_wl(ig)))*elabfactold
          ig = ig-1
        else
!-- resampling wavelength in lowest group
          call rnd_r(r1,rndstate)
          wl=1d0/(r1*grp_wlinv(1) + (1d0-r1)*grp_wlinv(2))
          wl = wl*elabfact
        endif
     endif
!-- check if ddmc region
     if (((grd_sig(ic)+grd_cap(ig,ic))*dx(ix) &
           >= trn_tauddmc) .and. .not.in_puretran) then
        ptcl2%itype = 2
!-- velocity effects accounting
        totevelo=totevelo+e*vx*mu*cinv
!
        e = e*elabfact
        e0 = e0*elabfact
        wl = wl/elabfact
     endif
  else
!    stop 'transport11: invalid distance'
     ierr = 17
     return
  endif

! !-- updating the co-moving group for superluminal flow for \tilde{v} frame
!   if(lsuplum.and.ptcl2%itype == 1) then
! !-- update velocity coordinate
!     vhelp1 = grd_vxarr(ix)
!     vhelp2 = grd_vxarr(ix+1)
!     help1 = grd_xarr(ix)
!     help2 = grd_xarr(ix+1)
!     vxhelp = (vhelp1*(help2-x)+vhelp2*(x-help1))/(help2-help1)
!     labfacthelp = 1d0-vxhelp*mu/pc_c
!     !-- comoving wavelength
!     wlhelp = wl/labfacthelp
!     !-- comoving group index
!     ig = binsrch(wlhelp,grp_wl,grp_ng+1,.false.)
! !-- check for ddmc again for new ig
!     if (((grd_sig(ic)+grd_cap(ig,ic))*dx(ix) &
!            >= trn_tauddmc) .and. .not.in_puretran) then
!       ptcl2%itype = 2
! !-- velocity effects accounting
!       e = e*labfacthelp
!       e0 = e0*labfacthelp
!       wl = wlhelp
!     endif
!   endif

end subroutine transport11
! vim: fdm=marker
