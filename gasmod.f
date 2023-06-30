*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module gasmod
c
      implicit none
c***********************************************************************
c gas grid structure
c***********************************************************************
      integer,parameter :: gas_nelem=30
c
c-- input parameter (NLTE)
       real*8 :: R_phot !radius where Rossland mean optical depth = 1
c
c-- wavelength grid (gridmod has a copy as well)
      integer,private :: ng=0
c
c-- domain decomposed grid variables used to calculate the state of the material (gas)
      integer :: gas_ncell=0
      integer :: gas_icell1=0
      real*8,allocatable :: gas_temp(:)       !(ncell) !average temperature
      real*8,allocatable :: gas_radtemp(:)    !(ncell) !radiation temperature
      real*8,allocatable :: gas_lum(:)        !(ncell) !experimental cell luminosity as source
      real*8,allocatable :: gas_eraddens(:)   !radiation energy density [erg/cm^3]
      real*8,allocatable :: gas_ur(:)         !radiation energy [erg]
      real*8,allocatable :: gas_rho(:)        !density [g/cm^3]
      real*8,allocatable :: gas_vol(:)        !cell volume [cm^3]
      real*8,allocatable :: gas_rcell(:)      !cell center radius [cm]
      real*8,allocatable :: gas_divv(:)       !divergence of velocity
      real*8,allocatable :: gas_mass(:)       !cell mass [g]
      real*8,allocatable :: gas_massfr(:,:)   !(gas_nelem,gas_ncell) ! mass fractions for output
      real*8,allocatable :: gas_ye(:)         !electron fraction
      real*8,allocatable :: gas_natom(:)      !cell number of atoms
      real*8,allocatable :: gas_nelec(:)      !cell number of electrons per atom
      real*8,allocatable :: gas_natomfr(:,:)  !(0:gas_nelem,ncell)  !natom fractions (>0:stable, 0:container for unused elements)

c== DD copies
c-- Line+Cont. extinction coeff
      real*4,allocatable :: gas_cap(:,:) !(ng,ncell)
c-- Line+Cont. "emission opacity !NLTE
      real*4,allocatable :: gas_capemit(:,:) !(ng,ncell)
c-- Line+Cont. emissivity !OUTPUT
      real*4,allocatable :: gas_emiss(:,:) !(ng,ncell)
c-- scattering coefficient
      real*8,allocatable :: gas_sig(:) !(ncell)
c-- Planck opacity (gray)
      real*8,allocatable :: gas_capgrey(:) !(ncell)
c-- "Emission" opacity (gray)
      real*8,allocatable :: gas_capemitgrey(:) !(ncell)
c-- Rosseland mean opacity
      real*8,allocatable :: gas_capross(:) !(ncell)
c
      real*8,allocatable :: gas_emit(:) !(ncell) amount of fictitious thermal energy emitted per cell in a time step
c-- radiation intensity !NLTE
      real*8,allocatable :: gas_jrad(:,:) !(ng,ncell) radiation intensity per cell per group
c
      save
c
      contains
c
c
      subroutine gasmod_init(ltalk,icell1,ncell,ngin)
c----------------------------------------
      implicit none
      logical,intent(in) :: ltalk
      integer,intent(in) :: icell1,ncell,ngin
************************************************************************
* Allocate gas variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n
c
      ng = ngin
c
      gas_icell1 = icell1
      gas_ncell = ncell
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = gas_ncell*(8*(18 + (2*gas_nelem + 1)) +
     &   4*4*(2 + ng))/1024 !kB
       write(6,*) 'ALLOC gas      :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif !ltalk
c
c-- ndim=1 alloc
      allocate(gas_temp(gas_ncell))
      allocate(gas_radtemp(gas_ncell))
      allocate(gas_lum(gas_ncell))
      allocate(gas_ur(gas_ncell))
      allocate(gas_rho(gas_ncell))
      allocate(gas_vol(gas_ncell))
      allocate(gas_rcell(gas_ncell))
      allocate(gas_divv(gas_ncell))
      allocate(gas_mass(gas_ncell))
      allocate(gas_ye(gas_ncell))
      allocate(gas_natom(gas_ncell))
      gas_natom = 0d0
      allocate(gas_nelec(gas_ncell))
      gas_nelec = 1d0
      allocate(gas_sig(gas_ncell))
      allocate(gas_capgrey(gas_ncell))
      allocate(gas_capemitgrey(gas_ncell))
      allocate(gas_capross(gas_ncell))
c
      allocate(gas_eraddens(gas_ncell))
c
      allocate(gas_emit(gas_ncell))
c
c-- ndim=2 alloc small
      allocate(gas_natomfr(0:gas_nelem,gas_ncell))
      gas_natomfr = 0d0
      allocate(gas_massfr(gas_nelem,gas_ncell))
      gas_massfr = 0d0
c
c-- ndim=2 alloc big
      allocate(gas_cap(ng,gas_ncell))
      allocate(gas_capemit(ng,gas_ncell))
      allocate(gas_emiss(ng,gas_ncell))
      allocate(gas_jrad(ng,gas_ncell))

      end subroutine gasmod_init
c
c
      subroutine gas_dealloc
      deallocate(gas_temp)
      deallocate(gas_radtemp)
      deallocate(gas_lum)
      deallocate(gas_ur)
      deallocate(gas_rho)
      deallocate(gas_vol)
      deallocate(gas_rcell)
      deallocate(gas_divv)
      deallocate(gas_mass)
      deallocate(gas_massfr)
      deallocate(gas_ye)
      deallocate(gas_natom)
      deallocate(gas_nelec)
      deallocate(gas_sig)
      deallocate(gas_capgrey)
      deallocate(gas_capemitgrey)
      deallocate(gas_capross)
      deallocate(gas_eraddens)
      deallocate(gas_emit)
      deallocate(gas_jrad)
      deallocate(gas_natomfr)
      deallocate(gas_cap)
      deallocate(gas_capemit)
      deallocate(gas_emiss)
      end subroutine gas_dealloc
c
      end module gasmod
c vim: fdm=marker
