*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine physical_opacity_nlte(it)
c     ---------------------------
c$    use omp_lib
      use physconstmod
      use inputparmod
      use ffxsmod
      use bfxsmod, only:bfxs
      use bbxsmod, only:bb_xs,bb_nline
      use nltemod
      use ionsmod
      use gasmod
      use groupmod
      use miscmod
      use timingmod
      use elemdatamod, only:elem_data
      use mpimod
      implicit none
      integer,intent(in) :: it
************************************************************************
* compute NLTE bound-bound and LTE bound-free, free-free opacity.
************************************************************************
      integer :: i,j
      real*8 :: wlinv
c-- timing
      real*8 :: t0,t1,t2,t3,t4
c-- subgridded wavelength
      integer :: ngsub
      real*8 :: wlsub(grp_ng*grp_ngs+1)
c-- helper arrays
      real*8 :: grndlev(gas_ncell,ion_iionmax-1,gas_nelem)
      real*8 :: hckt(gas_ncell),tempinv(gas_ncell)
      real*8 :: hlparr(gas_ncell)
c-- nlte
      integer :: llw,lup,nlte_nlev,nlte_nlin,coll_nlin ! levels and lines
      integer :: g_llw,g_lup ! statistical weights
      real*8 :: Aul,Bul,Blu,Cul,Clu ! Einstein coefficients
      real*8 :: npop_llw,npop_lup
      real*8 :: u0,x_coll,f_coll,gamma,gamma_d ! Electron-Impact Excitation (EIE)
      real*8 :: x_PI,f_PI,x_RR,f_RR ! Photo-ionization and radiative recombination
      real*8,allocatable :: gamma_PI(:),gamma_RR(:)
      real*8,allocatable :: rcoeff(:,:),rmatrix(:,:),rinv(:,:) !(nlte_nlev,nlte_nlev) !rate matrices
      real*8,allocatable :: npophelp(:),nrate(:) !nlte_nlev !excitation level population column matrices n, and dn/dt
      real,parameter :: fconst = sngl(pc_pi*pc_e**2/(pc_me*pc_c))
      real*8,parameter :: kb_ev = pc_kb / pc_ev
c-- output for debugging
      logical :: do_output = .false.
c-- ffxs
      real*8,parameter :: c1 = 4d0*pc_e**6/(3d0*pc_h*pc_me*pc_c**4)*
     &  sqrt(pc_pi2/(3*pc_me*pc_h*pc_c))
      real*8 :: gg,u,gff,help
      real*8 :: rgg,ru,dgg,du
      integer :: iu,igg
c-- bfxs
      integer :: il,ig,igs,iigs,iz,ii,ie
      integer :: l,ll,iln,reclen
      logical :: dirty
      real*8 :: en,xs,wl
c-- bbxs
      real*8 :: phi,ocggrnd,wl0,dwl
      real*8 :: expfac(gas_ncell)
      real*8 :: caphelp,emithelp
      real*8 :: wlm,B_nu,J_nu
      real*8 :: specarr(grp_ng)
      real*8,parameter :: ftpi4 = 15d0/pc_pi**4
c-- temporary cap array in the right order
      real*8,allocatable :: cap(:,:) !(gas_ncell,grp_ng*grp_ngs)
      real*8,allocatable :: capemit(:,:) !(gas_ncell,grp_ng*grp_ngs)
      real*8,allocatable :: emiss(:,:) !(gas_ncell,grp_ng*grp_ngs)
      real*8 :: capp(gas_ncell),capr(gas_ncell)
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- warn once
      logical :: lwarn
      character(40) :: fname
c-- temporary population arrays
      type level_populations
        real*8,allocatable :: n(:) !nlev
      end type level_populations
      type(level_populations),
     &  dimension(gas_ncell,nlte_nspecies,nlte_nelem) :: npop
c
c-- subgroups
      ngsub = grp_ng*grp_ngs
c
c-- initialize
      allocate(cap(gas_ncell,ngsub))
      cap = 0d0
      allocate(capemit(gas_ncell,ngsub)) !"emission opacity"
      capemit = 0d0
      allocate(emiss(gas_ncell,ngsub)) !emissivity
      emiss = 0d0
c
c-- subgridded wavelength
      igs = 1
      wlsub(igs) = grp_wl(1)
      if(grp_ngs==1) then
       wlsub = grp_wl
      else
       help = 1d0/grp_ngs
       do ig=1,grp_ng
        do i=1,grp_ngs
         igs = igs + 1
         wlsub(igs) = (1d0 - i*help)*grp_wl(ig) + i*help*grp_wl(ig+1)
        enddo
       enddo
      endif
c
c-- helper arrays
      tempinv = 1.d0/gas_temp
      hckt = pc_h*pc_c/(pc_kb*gas_temp)
c
c-- thomson scattering
      if(in_nothmson) then
       gas_sig = 0d0
      else
       gas_sig = cthomson*gas_nelec*gas_natom/gas_vol
       where(gas_mass<=0d0) gas_sig = 0d0
      endif
c
      t0 = t_time()
c
c-- bound-bound
      if(.not. in_nobbopac) then
       do iz=1,gas_nelem
        do i=1,gas_ncell
         if(gas_mass(i)<=0d0) cycle
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &     grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)*
     &     ion_grndlev(iz,i)%ginv(ii)
        enddo !i
       enddo !iz
c
c-- NLTE opacity calculation
       do iz=1,nlte_nelem
        do ii=1,iz
         nlte_nlev = nlte_data(iz)%i(ii)%nlev
         nlte_nlin = nlte_data(iz)%i(ii)%nlin
         coll_nlin = nlte_data(iz)%i(ii)%coll_nlin
c-- allocate temporary arrays
         allocate(gamma_PI(nlte_nlev),gamma_RR(nlte_nlev))
         allocate(rcoeff(nlte_nlev,nlte_nlev))
         allocate(rmatrix(nlte_nlev,nlte_nlev))
         allocate(rinv(nlte_nlev,nlte_nlev))
         allocate(npophelp(nlte_nlev),nrate(nlte_nlev))
c-- rate and population calculations
         do i=1,gas_ncell
          if(gas_mass(i)<=0d0) cycle
          allocate(npop(i,ii,iz)%n(nlte_nlev))
          rcoeff = 0d0
          rmatrix = 0d0
          rinv = 0d0
          npophelp = 0d0
          nrate = 0d0
          nrate(1) = 1d0
          npop(i,ii,iz)%n = 0d0
c-- radiative terms
c
c$omp parallel do schedule(static)
c$omp& private(il,llw,lup,Aul,Bul,Blu,wl0,igs,J_nu)
c$omp& shared(nlte_data,gas_ncell,gas_mass,gas_temp,
c$omp&   gas_nelec,gas_natom,gas_vol,nlte_ndens)
          do il=1,nlte_nlin
           llw = nlte_data(iz)%i(ii)%llw(il)
           lup = nlte_data(iz)%i(ii)%lup(il)
           Aul = nlte_data(iz)%i(ii)%Aul(il)
           Bul = nlte_data(iz)%i(ii)%Bul(il)
           Blu = nlte_data(iz)%i(ii)%Blu(il)
           wl0 = nlte_data(iz)%i(ii)%wl0(il)*pc_ang  !in cm
c-- igs pointer
           do igs=0,ngsub
            if(wlsub(igs+1)>wl0) exit
           enddo !igs
c-- line in group
           if(igs<1) cycle
           if(igs>ngsub) cycle !can't exit in omp
c
           wlm = 0.5d0*(grp_wl(igs+1)+grp_wl(igs))
           J_nu = (wlm**2/pc_c)*planck(wlm,gas_temp(i))
c-- calculate rate coefficients r(u->l) and r(l->u) for each level pair l,u
           do l=1,nlte_nlev
            do ll=1,nlte_nlev
             if(llw==l.and.lup==ll) then
              rcoeff(l,ll) = Blu*J_nu
             elseif(llw==ll.and.lup==l) then
              rcoeff(l,ll) = Aul + Bul*J_nu
             endif
            enddo !ll
           enddo !l
          enddo !il
c$omp end parallel do
c
          if(any(rcoeff/=rcoeff)) stop 'nlte opacity:rcoeff NaN'
c
c-- collisional (EIE) terms
c$omp parallel do schedule(static)
c$omp& private(il,x_coll,f_coll,gamma,llw,lup,g_llw,g_lup,
c$omp&    gamma_d,Clu,Cul)
c$omp& shared(i,iz,ii,gas_ncell,gas_mass,gas_temp,gas_nelec,
c$omp&   gas_natom,gas_vol,tempinv,nlte_ndens,nlte_data,
c$omp&   nlte_nlev,rcoeff)
          do il=1,coll_nlin
           x_coll = nlte_data(iz)%i(ii)%delE(il)*tempinv(i)/kb_ev
           f_coll = nlte_data(iz)%i(ii)%C0(il) +
     &        nlte_data(iz)%i(ii)%C1(il)*log(x_coll) +
     &        nlte_data(iz)%i(ii)%C2(il)*log(x_coll)**2 +
     &        nlte_data(iz)%i(ii)%C3(il)*log(x_coll)**3
           f_coll = exp(f_coll)
           gamma = 1.58d-5*(f_coll/x_coll)*exp(-1d0*x_coll)*
     &        (tempinv(i)**1.5)/(kb_ev**1.5)
c-- de-excitation rates
           llw = nlte_data(iz)%i(ii)%llw_coll(il)
           lup = nlte_data(iz)%i(ii)%lup_coll(il)
           g_llw = nlte_data(iz)%i(ii)%glev(llw)
           g_lup = nlte_data(iz)%i(ii)%glev(lup)
           gamma_d = g_llw*gamma*exp(x_coll)/g_lup
c-- collisional rate coefficients
           Clu = gamma*gas_nelec(i)*gas_natom(i)/gas_vol(i)
           Cul = gamma_d*gas_nelec(i)*gas_natom(i)/gas_vol(i)
c-- update rate coefficients r(u->l) and r(l->u) for collision
           do l=1,nlte_nlev
            do ll=1,nlte_nlev
             if(llw==l.and.lup==ll) then
              rcoeff(l,ll) = rcoeff(l,ll) + Clu
             elseif(llw==ll.and.lup==l) then
              rcoeff(l,ll) = rcoeff(l,ll) + Cul
             endif
            enddo !ll
           enddo !l
          enddo !il
c$omp end parallel do
c
          if(any(rcoeff/=rcoeff)) stop 'nlte opacity:rcoeff NaN'
c-- PI/RR rates
          do ll=1,nlte_nlev
            l = nlte_data(iz)%i(ii)%lev_PI(ll)
            x_PI = nlte_data(iz)%i(ii)%delE_PI(ll)*tempinv(i)/kb_ev
            f_PI = nlte_data(iz)%i(ii)%C0_PI(ll) +
     &        nlte_data(iz)%i(ii)%C1_PI(ll)*log(x_PI) +
     &        nlte_data(iz)%i(ii)%C2_PI(ll)*log(x_PI)**2 +
     &        nlte_data(iz)%i(ii)%C3_PI(ll)*log(x_PI)**3
            f_PI = exp(f_PI)
            gamma_PI(l) = (f_PI/x_PI)*exp(-1d0*x_PI)*
     &        (tempinv(i)**1.5)/(kb_ev**1.5)
c
            x_RR = kb_ev*gas_radtemp(i)
            f_RR = nlte_data(iz)%i(ii)%C0_RR(ll) +
     &        nlte_data(iz)%i(ii)%C1_RR(ll)*log(x_RR) +
     &        nlte_data(iz)%i(ii)%C2_RR(ll)*log(x_RR)**2 +
     &        nlte_data(iz)%i(ii)%C3_RR(ll)*log(x_RR)**3
            f_RR = exp(f_RR)
            gamma_RR(l) = f_RR
            if(l>1) nrate(l) = -1d0*gamma_RR(l)*nlte_ndens(i,ii+1,iz)/
     &            nlte_ndens(i,ii,iz) ! RHS of rate matrix, gamma_RR*n_ionized/n_neutral
          enddo !ll
c-- populate rate matrix using rate coefficients
          do l=1,nlte_nlev
           do ll=1,nlte_nlev
            if(l==1) then
             rmatrix(l,ll) = 1d0
            elseif(l==ll) then
             rmatrix(l,ll) = -1d0*sum(rcoeff(l,:))
             rmatrix(l,ll) = rmatrix(l,ll) - gamma_PI(l)
            else
             rmatrix(l,ll) = rcoeff(ll,l)
            endif
           enddo !ll
          enddo !l
          if(any(rmatrix/=rmatrix)) stop 'nlte opacity:rmatrix NaN'
c-- calcualte inverse of the rate matrix
          rinv = inv(rmatrix,nlte_nlev,nlte_nlev)
          if(any(rinv/=rinv)) stop 'nlte opacity:rinv NaN'
c-- calculate level populations
          npophelp = matmul(rinv,nrate)
c-- save the populations for each gas cell
          npop(i,ii,iz)%n(:) = npophelp
c-- DEBUG
!           if(do_output) then
!             write(fname,'("npop_nlte.it",i1,".r",i1)') it,impi
!             open(2,file=fname,status='unknown',position='append')
!             write(2,'(*(es12.5))')
!      &       gas_rcell(i),gas_radtemp(i),npop(i,ii,iz)%n(:)
!             close(2)
!
!             write(fname,'("ndens_nlte.it",i1,".r",i1)') it,impi
!             open(3,file=fname,status='unknown',position='append')
!             write(3,'(*(es12.5))')
!      &       gas_rcell(i),gas_radtemp(i),nlte_ndens(i,ii,iz)
!             close(3)
!           endif
c
         enddo !i
         deallocate(gamma_PI,gamma_RR)
         deallocate(rcoeff,rmatrix,rinv,npophelp,nrate)
        enddo !ii
       enddo !iz
c
c
       dirty = .true.
       phi = 0d0
       do iz=1,nlte_nelem
        do ii=1,iz
         nlte_nlin = nlte_data(iz)%i(ii)%nlin
c
c$omp parallel do schedule (static)
c$omp& private(wl0,wl,wlinv,dwl,phi,caphelp,emithelp,igs,iigs,
c$omp&   dirty,B_nu,llw,lup,g_llw,g_lup,npop_llw,npop_lup)
c$omp& shared(ngsub,wlsub,gas_mass,gas_temp,gas_ncell,
c$omp&   npop,nlte_ndens,nlte_data,cap,capemit,emiss)
         do iln=1,nlte_nlin
          wl0 = nlte_data(iz)%i(ii)%wl0(iln)*pc_ang
c-- igs pointer
          do igs=0,ngsub
           if(wlsub(igs+1)>wl0) exit
           dirty = .true.
          enddo !igs
c-- line in group
          if(igs<1) cycle
          if(igs>ngsub) cycle !can't exit in omp
c
c-- update
          if(dirty) then
           dirty = .false.
           wl = .5*(wlsub(igs) + wlsub(igs+1))
           wlinv = 1d0/wl
           dwl = wlsub(igs+1) - wlsub(igs)  !in cm
c-- approximate profile function
           phi = wl**2/(dwl*pc_c)
          endif
          llw = nlte_data(iz)%i(ii)%llw(iln)
          lup = nlte_data(iz)%i(ii)%lup(iln)
          g_llw = nlte_data(iz)%i(ii)%gl(iln)
          g_lup = nlte_data(iz)%i(ii)%gu(iln)
          do i=1,gas_ncell
           if(gas_mass(i)<=0d0) cycle
           npop_llw = npop(i,ii,iz)%n(llw)*nlte_ndens(i,ii,iz) ! occupation number density of lower state
           npop_lup = npop(i,ii,iz)%n(lup)*nlte_ndens(i,ii,iz) ! occupation number density of upper state
c-- absorption opacity
           caphelp = npop_llw*fconst*nlte_data(iz)%i(ii)%f(iln)
     &      *phi*(1d0-(npop_lup*g_llw/(npop_llw*g_lup)))
           cap(i,igs) = cap(i,igs) + caphelp
c-- emission opacity
           B_nu = (wl**2/pc_c)*planck(wl,gas_temp(i))
c
           emithelp = npop_lup*nlte_data(iz)%i(ii)%Aul(iln)*
     &        phi*pc_h*pc_c/(4d0*pc_pi*wl0)
c
           capemit(i,igs) = capemit(i,igs) + emithelp/B_nu
c-- emissivity for output
           emiss(i,igs) = emiss(i,igs) + emithelp
          enddo !i
         enddo !iln
c$omp end parallel do
        enddo !ii
       enddo !iz
c-- LTE opacity calculations
       iigs = 0
       dirty = .true.
       phi = 0d0
c$omp parallel do schedule (static)
c$omp& private(wl0,wl,wlinv,dwl,phi,ocggrnd,expfac,igs,iigs,
c$omp&   caphelp,emithelp,dirty,B_nu)
c$omp& shared(ngsub,gas_mass,gas_temp,grndlev,hckt,wlsub,
c$omp&   gas_ncell,nlte_ndens,nlte_data,cap,capemit,emiss)
       do il=1,bb_nline
        wl0 = bb_xs(il)%wl0*pc_ang  !in cm
        iz = bb_xs(il)%iz
        if(iz.le.nlte_nelem) cycle ! precaution
        ii = bb_xs(il)%ii
c-- igs pointer
        do igs=iigs,ngsub
         if(wlsub(igs+1)>wl0) exit
         dirty = .true.
        enddo !igs
        iigs = igs
c-- line in group
        if(igs<1) cycle
        if(igs>ngsub) cycle !can't exit in omp
c
c-- update
        if(dirty) then
         dirty = .false.
         wl = .5*(wlsub(igs) + wlsub(igs+1))
         wlinv = 1d0/wl
         dwl = wlsub(igs+1) - wlsub(igs)  !in cm
c-- approximate expfac
         expfac = 1d0 - exp(-hckt*wlinv)
c-- approximate profile function
         phi = wl**2/(dwl*pc_c)
        endif
c
c-- profile function
*       phi = wl0**2/(dwl*pc_c)  !exact phi
c
c-- evaluate cap
        do i=1,gas_ncell
         if(gas_mass(i)<=0d0) cycle
         ocggrnd = grndlev(i,ii,iz)
         if(ocggrnd<=0d0) cycle
*        expfac = 1d0 - exp(-hckt(i)/wl0)  !exact expfac
         caphelp = phi*bb_xs(il)%gxs*ocggrnd*
     &    exp(-bb_xs(il)%chilw*hckt(i))*expfac(i)
         cap(i,igs) = cap(i,igs) + caphelp
c-- emission opacity
         capemit(i,igs) = capemit(i,igs) + caphelp
c-- emissivity for output
         B_nu = (wl**2/pc_c)*planck(wl,gas_temp(i)) !factor of lambda**2/c for B_nu to B_lambda conversion
         emithelp = caphelp*B_nu
         emiss(i,igs) = emiss(i,igs) + emithelp
        enddo !i
       enddo !il !bb_nline
c
c$omp end parallel do
c
      endif !in_nobbopac
c
      t1 = t_time()
c
c-- bound-free
      if(.not. in_nobfopac) then
c
       do iz=1,gas_nelem
        do i=1,gas_ncell
         if(gas_mass(i)<=0d0) cycle
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)
        enddo !i
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs,B_nu)
c$omp& shared(ngsub,gas_mass,grndlev,wlsub,cap,capemit,emiss)
       do igs=1,ngsub
        wl = wlsub(igs)  !in cm
        en = pc_h*pc_c/(pc_ev*wl) !photon energy in eV
        do iz=1,gas_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          do i=1,gas_ncell
            if(gas_mass(i)<=0d0) cycle
            cap(i,igs) = cap(i,igs) + xs*pc_mbarn*grndlev(i,ii,iz)
c-- emission opacity
            capemit(i,igs) = capemit(i,igs) +
     &          xs*pc_mbarn*grndlev(i,ii,iz)
c-- emissivity for output
            B_nu = (wl**2/pc_c)*planck(wl,gas_temp(i))
            emiss(i,igs) = emiss(i,igs) +
     &          xs*pc_mbarn*grndlev(i,ii,iz)*B_nu
          enddo !i
         enddo !ii
        enddo !iz
!       write(6,*) 'wl done:',igs !DEBUG
!       write(6,*) cap(:,igs) !DEBUG
       enddo !igs
c$omp end parallel do
c
      endif !in_nobfopac
c
      t2 = t_time()
c
c-- free-free
      if(.not. in_noffopac) then
c
c-- simple variant: nearest data grid point
       hlparr = (gas_natom/gas_vol)**2*gas_nelec
       lwarn = .true.
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,wlinv,u,ru,du,iu,help,gg,rgg,dgg,igg,gff,B_nu)
c$omp& shared(ngsub,lwarn,hckt,hlparr,gas_mass,wlsub,cap,capemit,emiss)
       do igs=1,ngsub
        wl = wlsub(igs)  !in cm
        wlinv = 1d0/wl  !in cm
c-- gcell loop
        do i=1,gas_ncell
         if(gas_mass(i)<=0d0) cycle
         u = hckt(i)*wlinv
         ru = 10d0*(log10(u) + 4d0) + 1d0
         iu = floor(ru)
         iu = max(iu,1)
         iu = min(iu,ff_nu-1)
         du = ru - iu
c
c-- element loop
         help = c1*sqrt(hckt(i))*(1d0 - exp(-u))*wl**3*hlparr(i)
         do iz=1,gas_nelem
          gg = iz**2*pc_rydberg*hckt(i)
          rgg = 5d0*(log10(gg) + 4d0) + 1d0
          igg = floor(rgg)
          igg = max(igg,1)
          igg = min(igg,ff_ngg-1)
          dgg = rgg - igg
c
c-- bilinear inter-extrapolation
          gff = (1d0-du)*(1d0-dgg)*ff_gff(iu,igg) +
     &          du*(1d0-dgg)*ff_gff(iu+1,igg) +
     &          (1d0-du)*dgg*ff_gff(iu,igg+1) +
     &          du*dgg*ff_gff(iu+1,igg+1)
          if(rgg>ff_ngg) gff = max(gff,1d0)
c-- cross section
          cap(i,igs) = cap(i,igs) +
     &      help*gff*iz**2*gas_natomfr(iz,i)
          capemit(i,igs) = capemit(i,igs) +
     &      help*gff*iz**2*gas_natomfr(iz,i)
c-- emissivity for output
          if(capemit(i,igs)/=capemit(i,igs))
     &       write(6,*) i,grp_wl(igs)/pc_ang,capemit(i,igs)
          B_nu = (wl**2/pc_c)*planck(wl,gas_temp(i))
          emiss(i,igs) = emiss(i,igs) +
     &      help*gff*iz**2*gas_natomfr(iz,i)*B_nu
         enddo !iz
        enddo !i
       enddo !igs
c$omp end parallel do
c
      endif !in_noffopac
c
      t3 = t_time()
c
c-- collapse subgridding
      if(grp_ngs/=1) then
       help = 1d0/grp_ngs  !assume evenly spaced subgroups
       do ig=1,grp_ng
        i = (ig-1)*grp_ngs + 1
        j = i + grp_ngs - 1
c-- absorption opacity
c-- first read
        capp = sum(cap(:,i:j), dim=2, mask=cap(:,i:j)>0d0)*help
        capr = sum(1d0/cap(:,i:j), dim=2, mask=cap(:,i:j)>0d0)*help
        where(capr>0d0) capr = 1d0/capr
c-- then overwrite
        cap(:,ig) = (1d0-in_opacmixrossel)*capp + in_opacmixrossel*capr
c-- "emission opacity"
c-- first read
        capp = 0d0
        capr = 0d0
        capp = sum(capemit(:,i:j), dim=2, mask=capemit(:,i:j)>0d0)*
     &    help
        capr = sum(1d0/capemit(:,i:j), dim=2,
     &    mask=capemit(:,i:j)>0d0)*help
        where(capr>0d0) capr = 1d0/capr
c-- overwrite
        capemit(:,ig) = (1d0-in_opacmixrossel)*capp +
     &    in_opacmixrossel*capr
c-- emissivity
        emiss(:,ig) = sum(emiss(:,i:j), dim=2, mask=emiss(:,i:j)>0d0)*
     &    help
       enddo !ig
      endif
c
c-- save
      gas_cap = transpose(sngl(cap(:,:grp_ng)))
      gas_capemit = transpose(sngl(capemit(:,:grp_ng)))
      gas_emiss = transpose(sngl(emiss(:,:grp_ng)))
c
c-- sanity check
      j = 0
      do i=1,gas_ncell
       if(gas_mass(i)<=0d0) cycle
       do ig=1,grp_ng
        if(gas_cap(ig,i)==0.) j = ior(j,1)
        if(gas_cap(ig,i)<0.) j = ior(j,2)
        if(gas_cap(ig,i)/=gas_cap(ig,i)) j = ior(j,4)
        if(gas_cap(ig,i)>huge(gas_cap)) j = ior(j,8)
       enddo !ig
      enddo !i
      if(iand(j,1)/=0) write(0,*) 'opacity_calc: some cap==0'
      if(iand(j,2)/=0) write(0,*) 'opacity_calc: some cap<0'
      if(iand(j,4)/=0) write(0,*) 'opacity_calc: some cap==NaN'
      if(iand(j,8)/=0) write(0,*) 'opacity_calc: some cap==inf'
c-- refine sanity check
      if(j>0) then
       j = 0
       do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
        do ig=1,grp_ng
         if(gas_cap(ig,i)/=0.) j = ior(j,1)
         if(gas_cap(ig,i)>=0.) j = ior(j,2)
         if(gas_cap(ig,i)==gas_cap(ig,i)) j = ior(j,4)
         if(gas_cap(ig,i)<=huge(gas_cap)) j = ior(j,8)
        enddo !ig
       enddo !i
       if(iand(j,1)==0) write(0,*) 'opacity_calc: all cap==0'
       if(iand(j,2)==0) write(0,*) 'opacity_calc: all cap<0'
       if(iand(j,4)==0) write(0,*) 'opacity_calc: all cap==NaN'
       if(iand(j,8)==0) write(0,*) 'opacity_calc: all cap==inf'
      endif
c
c-- sanity check ! capemit
      j = 0
      do i=1,gas_ncell
       if(gas_mass(i)<=0d0) cycle
       do ig=1,grp_ng
        if(gas_capemit(ig,i)==0.) j = ior(j,1)
        if(gas_capemit(ig,i)<0.) j = ior(j,2)
        if(gas_capemit(ig,i)/=gas_capemit(ig,i)) j = ior(j,4)
        if(gas_capemit(ig,i)>huge(gas_capemit)) j = ior(j,8)
       enddo !ig
      enddo !i
      if(iand(j,1)/=0) write(0,*) 'opacity_calc: some capemit==0'
      if(iand(j,2)/=0) write(0,*) 'opacity_calc: some capemit<0'
      if(iand(j,4)/=0) write(0,*) 'opacity_calc: some capemit==NaN'
      if(iand(j,8)/=0) write(0,*) 'opacity_calc: some capemit==inf'
c-- refine sanity check
      if(j>0) then
       j = 0
       do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
        do ig=1,grp_ng
         if(gas_capemit(ig,i)/=0.) j = ior(j,1)
         if(gas_capemit(ig,i)>=0.) j = ior(j,2)
         if(gas_capemit(ig,i)==gas_capemit(ig,i)) j = ior(j,4)
         if(gas_capemit(ig,i)<=huge(gas_capemit)) j = ior(j,8)
        enddo !ig
       enddo !i
       if(iand(j,1)==0) write(0,*) 'opacity_calc: all capemit==0'
       if(iand(j,2)==0) write(0,*) 'opacity_calc: all capemit<0'
       if(iand(j,4)==0) write(0,*) 'opacity_calc: all capemit==NaN'
       if(iand(j,8)==0) write(0,*) 'opacity_calc: all capemit==inf'
      endif
c
c-- sanity check ! emiss
      j = 0
      do i=1,gas_ncell
       if(gas_mass(i)<=0d0) cycle
       do ig=1,grp_ng
        if(gas_emiss(ig,i)==0.) j = ior(j,1)
        if(gas_emiss(ig,i)<0.) j = ior(j,2)
        if(gas_emiss(ig,i)/=gas_emiss(ig,i)) j = ior(j,4)
        if(gas_emiss(ig,i)>huge(gas_emiss)) j = ior(j,8)
       enddo !ig
      enddo !i
      ! if(iand(j,1)/=0) write(0,*) 'opacity_calc: some emiss==0'
      if(iand(j,2)/=0) write(0,*) 'opacity_calc: some emiss<0'
      if(iand(j,4)/=0) write(0,*) 'opacity_calc: some emiss==NaN'
      if(iand(j,8)/=0) write(0,*) 'opacity_calc: some emiss==inf'
c-- refine sanity check
      if(j>0) then
       j = 0
       do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
        do ig=1,grp_ng
         if(gas_emiss(ig,i)/=0.) j = ior(j,1)
         if(gas_emiss(ig,i)>=0.) j = ior(j,2)
         if(gas_emiss(ig,i)==gas_emiss(ig,i)) j = ior(j,4)
         if(gas_emiss(ig,i)<=huge(gas_emiss)) j = ior(j,8)
        enddo !ig
       enddo !i
       if(iand(j,1)==0) write(0,*) 'opacity_calc: all emiss==0'
       if(iand(j,2)==0) write(0,*) 'opacity_calc: all emiss<0'
       if(iand(j,4)==0) write(0,*) 'opacity_calc: all emiss==NaN'
       if(iand(j,8)==0) write(0,*) 'opacity_calc: all emiss==inf'
      endif
c
      deallocate(cap)
      deallocate(capemit)
      deallocate(emiss)
c
      t4 = t_time()
c-- register timing
      call timereg(t_opac,t4-t0)
      call timereg(t_bb,t1-t0)
      call timereg(t_bf,t2-t1)
      call timereg(t_ff,t3-t2)
c
      contains

      ! Returns the inverse of a matrix calculated by finding the LU
      ! decomposition.  Depends on LAPACK.
      function inv(A,size1,size2) result(Ainv)
      !     --------------------------------------------------
      implicit none
      integer,intent(in) :: size1,size2
      real*8,intent(in) :: A(size1,size2)
      real*8,dimension(size1,size2) :: Ainv

      real*8,dimension(size1) :: work  ! work array for LAPACK
      integer,dimension(size1) :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! SGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
        stop 'Matrix is numerically singular!'
      end if

      ! SGETRF computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
       stop 'Matrix inversion failed!'
      end if
      end function inv
c
      end subroutine physical_opacity_nlte
c vim: fdm=marker
