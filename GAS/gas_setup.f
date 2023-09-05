*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine gas_setup
c     --------------------
      use inputstrmod
      use physconstmod
      use inputparmod
      use gasmod
      use miscmod, only:warn
      use groupmod, only:grp_ng,grp_wl
      use gridmod
      implicit none
************************************************************************
* Initialize the gas grid
************************************************************************
      integer :: l,i
      integer :: ig
      ! real*8 :: W !dilution factor
      real*8 :: massfr(0:gas_nelem,gas_ncell) !0:container for unused elements
      real*8 :: wlm
c
c-- agnostic mass setup
      gas_mass = str_massdd
c
c-- volume from input.str if available
      if(str_lvol) gas_vol = str_voldd
c
c-- temperature
      gas_temp = str_tempdd
c-- radiation temperature
      if(str_lradtemp) then
c-- from input.str
        gas_radtemp = str_radtempdd
      else
c-- equal to gas_temp
        gas_radtemp = str_tempdd
      endif
c-- luminosity (experimental volume source)
      gas_lum = 0.d0
      if(str_llum) gas_lum = str_lumdd
c
c
c-- initial radiation
      gas_eraddens = pc_acoef*gas_radtemp**4
c-- disallow initial radiation in void cells
      where(gas_mass<=0d0) gas_eraddens = 0d0
c
c
c-- radiation energy
      gas_ur = pc_acoef*gas_radtemp**4
c--cell-centered radii
      i = 0
      gas_rcell = 0.d0
      gas_jrad = 0.d0
      do l=gas_icell1,gas_icell1+gas_ncell-1
        i = i+1
        gas_rcell(i) = grd_rcell(l)
      enddo !l
      if(any(gas_rcell<=0.d0)) stop "gas_setup: gas_rcell<=0"
      if(any(gas_rcell/=gas_rcell)) stop "gas_setup: NaN in gas_rcell"
c-- radiation intensity
c-- J_nu = B_nu for first iteration in NLTE
      if(in_nebular.or.in_nlte) then
        do i=1,gas_ncell
          !-- experimental dilution factor of J_nu in NLTE
          ! if(gas_rcell(i).gt.R_phot) then
          !   W = 0.5d0*(1d0-((1d0-(R_phot/gas_rcell(i))**2)**0.5))
          ! else
          !   W = 1d0
          ! endif
          !gas_jrad(:,i) = W*gas_jrad(:,i)
          do ig=1,grp_ng
            wlm = 0.5d0*(grp_wl(ig+1)+grp_wl(ig))
            gas_jrad(ig,i) = wlm**2*planck(wlm,gas_temp(i))/pc_c !first iteration
          enddo
        enddo
        if(any(gas_jrad<0.d0)) stop "gas_setup: gas_jrad<=0"
        if(any(gas_jrad/=gas_jrad)) stop "gas_setup: NaN in gas_jrad"
c
      endif
c-- adopt mass fractions from input file
      massfr = 0d0
      gas_massfr = 0d0
      if(.not.allocated(str_massfrdd))
     &        stop 'input.str mass fraction data not available'
      do l=1,str_nabund
        i = str_iabund(l)
        if(i>gas_nelem) i = 0 !divert to container
        massfr(i,:) = str_massfrdd(l,:)
      enddo
c
c-- convert mass fractions to # atoms
      call massfr2natomfr(massfr)
      gas_massfr = massfr(1:,:) !for output
c
c
      end subroutine gas_setup
c
c
c
      subroutine massfr2natomfr(massfr)
c     ----------------------------------
      use physconstmod
      use elemdatamod, only:elem_data
      use gasmod
      implicit none
      real*8,intent(inout) :: massfr(0:gas_nelem,gas_ncell)
************************************************************************
* convert mass fractions to natom fractions, and mass to natom.
************************************************************************
      integer :: i,l
      real*8 :: help
      logical :: lwarn
c
      lwarn = .true.
c
      do i=1,gas_ncell
       if(gas_mass(i)<=0) cycle
c-- sanity test
       if(all(massfr(1:,i)==0d0)) stop
     &    'massfr2natomfr: all mass fractions zero'
       if(any(massfr(1:,i)<0d0)) stop
     &    'massfr2natomfr: negative mass fractions'
c
c-- renormalize (the container fraction (unused elements) is taken out)
       help = sum(massfr(1:,i))
       massfr(:,i) = massfr(:,i)/help
c-- warn if renormalization factor is abnormal
       if(abs(help-1d0)>.01 .and. lwarn) then
        write(0,*) 'WARNING: large abundance normalization factor!',help
        write(6,*) 'WARNING: large abundance normalization factor!',help
        lwarn = .false. !warn once
       endif
c
c-- partial mass
       gas_natomfr(:,i) = massfr(:,i)*gas_mass(i)
c
c-- convert to natoms
       do l=1,gas_nelem
        gas_natomfr(l,i) = gas_natomfr(l,i)/(elem_data(l)%m*pc_amu)
       enddo !j
c
c-- calculate Ye
       gas_ye(i) = gas_natomfr(0,i)*.5d0 !container with unknown elements
c-- stable abundances
       do l=1,gas_nelem
        gas_ye(i) = gas_ye(i) + gas_natomfr(l,i)*l/elem_data(l)%m
       enddo
c-- normalize
       gas_ye(i) = gas_ye(i)/sum(gas_natomfr(:,i))
c
c-- total natom
       gas_natom(i) = sum(gas_natomfr(1:,i))
c
c-- convert natoms to natom fractions
       gas_natomfr(:,i) = gas_natomfr(:,i)/gas_natom(i)
c
      enddo !i
c
      end subroutine massfr2natomfr
c vim: fdm=marker
