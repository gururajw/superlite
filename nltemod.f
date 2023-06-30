*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module nltemod
c     --------------
      implicit none
************************************************************************
* read bound-bound cross sections (oscillator strengths)
************************************************************************
      integer :: nlte_nlevel=0 ! number of nLTE levels
      integer :: nlte_nline=0 ! number of nLTE radiative lines
      integer :: coll_nline=0 ! number of nLTE electron-impact excitation lines
      integer :: nlte_nelem=0 ! number of species treated with nlte 1 - h, 2 - he
      integer :: nlte_nspecies=0 ! number of species treated with nlte 1 - h1, 2 - h1 & he1, 3 - h1, he1 & he2
      real*8,allocatable :: nlte_ndens(:,:,:) ! (gas_ncell,nlte_species,nlte_nel) occupation numbers for ionization states
c
c-- raw data read from file
c-- level data type
      type nlte_raw_level_data
       integer :: id,g !label,statistical weight
       real*8 :: chi     !in cm^-1
       integer :: n !princical quantum number
      end type nlte_raw_level_data
      type(nlte_raw_level_data),allocatable :: nlte_level(:) !nlte_nlevel
c-- line data type
      type nlte_raw_line_data
       integer :: lev1,lev2
       real*8 :: f,wl0 !in ang
      end type nlte_raw_line_data
      type(nlte_raw_line_data),allocatable :: nlte_line(:) !nlte_nline
c-- Electron-impact excitation (EIE) data type
      type coll_raw_line_data
       integer :: lev1,lev2
       real*8 :: C0,C1,C2,C3 ! 4 fit coefficients
       real*8 :: n_coll,delE !pre-factor coeff,transition energy (eV)
      end type coll_raw_line_data
      type(coll_raw_line_data),allocatable :: coll_line(:) !nlte_nline
c-- PI/RR transition data type
      type pi_rr_raw_line_data
       integer :: lev
       real*8 :: C0_RR,C1_RR,C2_RR,C3_RR ! 4 fit coefficients for RR rates
       real*8 :: n_RR !pre-factor coeff for RR rates
       real*8 :: C0_PI,C1_PI,C2_PI,C3_PI ! 4 fit coefficients for PI rates
       real*8 :: n_PI !pre-factor coeff for PI rates
       real*8 :: delE !transition energy (eV)
      end type pi_rr_raw_line_data
      type(pi_rr_raw_line_data),allocatable :: pi_rr_data(:) !nlte_nlevel
c-- permanent data type
c-- NLTE data type for ionizaton stages
      type nlte_species_data
       integer :: nlev,nlin,coll_nlin !number of nlte levels and radiative and collision lines per species
      !-- level data (nlev)
       real*8,allocatable :: chilev(:) !in cm^-1
       integer,allocatable :: levid(:),glev(:),nplev(:) !label,statistical weight,princical quantum number
      !-- radiative line data (nlin)
       integer,allocatable :: llw(:),lup(:) !lower and upper level id for radiative lines
       real*8,allocatable :: wl0(:) !in ang
       real*8,allocatable :: f(:) !oscillator strength
       integer,allocatable :: gl(:),gu(:) !g_lower_level,g_upper_level
       real*8,allocatable :: Aul(:),Bul(:),Blu(:) ! Einstein coefficients
       real*8,allocatable :: qlu(:) ! collisional excitation rates !van Regemorter
       !-- EIE line data (coll_nlin)
       integer,allocatable :: llw_coll(:),lup_coll(:) !lower and upper level id for collison lines
       real*8,allocatable :: C0(:),C1(:),C2(:),C3(:) !fit coefficients
       real*4,allocatable :: n_coll(:) !pre-factor coeff
       real*8,allocatable :: delE(:) ! EIE transition energy in eV
       !-- PI/RR transition data (nlev)
       integer,allocatable :: lev_PI(:) !level id for PI/RR rates
       real*8,allocatable :: C0_PI(:),C1_PI(:),C2_PI(:),C3_PI(:) !fit coefficients for PI rates
       real*4,allocatable :: n_PI(:) !pre-factor coeff for PI rates
       real*8,allocatable :: C0_RR(:),C1_RR(:),C2_RR(:),C3_RR(:) !fit coefficients for RR rates
       real*4,allocatable :: n_RR(:) !pre-factor coeff for RR rates
       real*8,allocatable :: delE_PI(:) ! EIE transition energy in eV (for both PI and RR transition)
      end type nlte_species_data
c-- NLTE data type for an element
      type nlte_elem_data
       type(nlte_species_data),allocatable :: i(:) !ionization states !excluding bare nucleus
      end type nlte_elem_data
      type(nlte_elem_data),allocatable :: nlte_data(:) !nlte_nelem

      save
c
      contains
c
c
c
      subroutine read_atom_nlte(iz,ii,istat,get_data)
c     ------------------------------------------
      use miscmod, only:lcase
      use elemdatamod, only:elem_data
      implicit none
      integer,intent(in) :: iz,ii
      integer,intent(out) :: istat
      logical,intent(in) :: get_data
************************************************************************
* Read a single .nlte file with atomic data for nlte calculations,
* or just poll the number of lines it contains
************************************************************************
      character(80) :: word
      character(13) :: fname
      integer :: l,istat2
      integer(1) :: byte
c
c-- filename
      write(fname,'("data.nlte.",a,i1)')
     &  lcase(trim(elem_data(iz)%sym)),ii
      open(4,file='Atoms/'//trim(fname),status='old',action='read',
     &  iostat=istat)
      if(istat/=0) goto 66
c
c-- read data size
      read(4,*)
      read(4,*)
      read(4,*)
      read(4,*,iostat=istat) nlte_nlevel,nlte_nline,coll_nline,word
      if(nlte_nlevel<=0 .or. nlte_nline<=0) istat = 1
      if(istat/=0) goto 66
c-- poll ready, return
      if(.not.get_data) goto 67
c
c-- allocate
      allocate(nlte_level(nlte_nlevel),nlte_line(nlte_nline))
      allocate(coll_line(coll_nline),pi_rr_data(nlte_nlevel))
c
c-- read data
c-- level data
      read(4,'(i5,9x,i5,f11.3,i5)',iostat=istat) nlte_level
      if(istat/=0) then
       write(6,*) fname
       stop 'read_atom_nlte: level data read error'
      endif
c-- line data
      read(4,'(i4,1x,i4,3x,es11.5,3x,es11.5)',iostat=istat) nlte_line
      if(istat/=0) then
       write(6,*) fname
       stop 'read_atom_nlte: line data read error'
      endif
c-- electron-impact excitation data
      do l=1,coll_nline
       read(4,'(i3,1x,i3,4(1x,es13.6),1x,f5.3,1x,es13.6)',iostat=istat)
     &     coll_line(l)
c       write(6,*) coll_line(l)
      enddo
c-- PI/RR transition data
      do l=1,nlte_nlevel
       read(4,'(i3,2(4(1x,es13.6),1x,f5.3),1x,es13.6)',iostat=istat)
     &     pi_rr_data(l)
c       write(6,*) pi_rr_data(l)
      enddo
c-- verify eof
      read(4,*,iostat=istat2) byte
      if(istat2>=0) then
       write(6,*) fname,byte
       stop 'read_atom_nlte: data remaining on input file'
      endif
c
67    continue
      close(4)
      return
c
66    continue
      if(.not. get_data) write(6,*) 'read_atom_nlte failed:',iz,ii
      close(4)
c
      end subroutine read_atom_nlte
c
c
c
      subroutine nlte_alloc_nlev(nelem,ncell)
c     -------------------------------------------
      implicit none
      integer,intent(in) :: nelem,ncell
************************************************************************
* nlte_ndens stores the occupation number density of the
* ionization states excluding bare nucleus for NLTE elements in all
* gas_cells.
************************************************************************
c
      allocate(nlte_ndens(ncell,nelem+1,nelem))
c
      end subroutine nlte_alloc_nlev
c
c
c
      subroutine get_rphot
      use gridmod
      use gasmod
      use miscmod
************************************************************************
* Find photospheric radius using Rosseland mean opacity for next iteration
************************************************************************
      integer :: i,iphot
      real*8 :: tau(grd_ncell)
      real*8 :: tau_phot ! to do: add a parameter in inputparmod
      do i=1,grd_ncell
        tau(i) = grd_capross(i)*(grd_xarr(i+1)-grd_xarr(i))
      enddo
      tau = cumsum(tau)
      tau_phot = 2.d0/3.d0
      iphot = binsrch_decr(tau_phot,tau,grd_ncell,.false.)
      R_phot = grd_rcell(iphot) ! update photospheric radius
      ! write(6,*) R_phot
      contains

      pure function cumsum(A) result(B)
      real*8,intent(in) :: A(:)
      real*8 :: B(size(A))
************************************************************************
* Calculate and return cumulative sum of A from outer bound inward
************************************************************************
      real*8 :: R(size(A))
      integer :: i,n
      n = size(A)
      do i=1,size(A)
        B(i) = A(n)
        n = n - 1
      enddo
      R(:) = [(sum(B(1:i)),i=1,size(A))]
      n = size(A)
      do i=1,size(A)
        B(i) = R(n)
        n = n - 1
      enddo
      end function cumsum
c
      end subroutine get_rphot
c
c
c
      subroutine nlte_rad_est
      use gridmod
      use groupmod
      use physconstmod
      implicit none
************************************************************************
* Use path-length estimator to update initial radiation energy for
* next iteration
************************************************************************
      integer :: i,ig
      real*8, allocatable :: wlm(:),dwl(:)
      real*8,allocatable :: Trad(:)
c
      allocate(Trad(grd_ncell))
      do i=1,grd_ncell
c-- temperature estimator
       Trad(i) = (grd_tally(i)/(pc_acoef*grd_vol(i)))**.25
c-- dampen the estimate for convergence
       grd_radtemp(i) = grd_radtemp(i) + 0.5*(Trad(i)-grd_radtemp(i))
c-- radiation energy density estimator
       grd_jrad(:,i) = wlm(:)**2*planckv(wlm(:),grd_radtemp(i))/pc_c
c
      enddo
c
      deallocate(Trad)
c
      end subroutine nlte_rad_est
c
c
c
      subroutine nlte_alloc_data(nelem,nspecies,nlevel,nline,coll_nline)
c     ------------------------------------------------
      implicit none
      integer,intent(in) :: nelem,nspecies
      integer,intent(in) :: nlevel(nspecies),nline(nspecies)
      integer,intent(in) :: coll_nline(nspecies)
************************************************************************
* allocate the nlte_data structure and keep track of the size in bytes.
* This is used in MPI mode before broadcasting nlte_data.
************************************************************************
      integer :: ii,iz,isp,nlev,nlin,coll_nlin
c
c-- copy
      nlte_nelem = nelem
      nlte_nspecies = nspecies
c
c-- allocate nlte data structure
      allocate(nlte_data(nelem))
      isp = 0
      do iz=1,nelem
       allocate(nlte_data(iz)%i(iz))
c-- allocate compressed data
       do ii=1,iz
        isp = isp + 1
        !-- level data
        nlev = nlevel(isp)
        nlte_data(iz)%i(ii)%nlev = nlev
        allocate(nlte_data(iz)%i(ii)%chilev(nlev))
        allocate(nlte_data(iz)%i(ii)%levid(nlev))
        allocate(nlte_data(iz)%i(ii)%glev(nlev))
        allocate(nlte_data(iz)%i(ii)%nplev(nlev))
        nlte_data(iz)%i(ii)%chilev = 0d0
        nlte_data(iz)%i(ii)%levid = 0
        nlte_data(iz)%i(ii)%glev = 0
        nlte_data(iz)%i(ii)%nplev = 0
        !-- radiative line data
        nlin = nline(isp)
        nlte_data(iz)%i(ii)%nlin = nlin
        allocate(nlte_data(iz)%i(ii)%wl0(nlin))
        allocate(nlte_data(iz)%i(ii)%f(nlin))
        allocate(nlte_data(iz)%i(ii)%gu(nlin))
        allocate(nlte_data(iz)%i(ii)%gl(nlin))
        allocate(nlte_data(iz)%i(ii)%llw(nlin))
        allocate(nlte_data(iz)%i(ii)%lup(nlin))
        allocate(nlte_data(iz)%i(ii)%Aul(nlin))
        allocate(nlte_data(iz)%i(ii)%Bul(nlin))
        allocate(nlte_data(iz)%i(ii)%Blu(nlin))
        allocate(nlte_data(iz)%i(ii)%qlu(nlin))
        nlte_data(iz)%i(ii)%wl0 = 0d0
        nlte_data(iz)%i(ii)%f = 0d0
        nlte_data(iz)%i(ii)%gl = 0
        nlte_data(iz)%i(ii)%gu = 0
        nlte_data(iz)%i(ii)%llw = 0
        nlte_data(iz)%i(ii)%lup = 0
        nlte_data(iz)%i(ii)%Aul = 0d0
        nlte_data(iz)%i(ii)%Bul = 0d0
        nlte_data(iz)%i(ii)%Blu = 0d0
        nlte_data(iz)%i(ii)%qlu = 0d0
        !-- collison line data
        coll_nlin = coll_nline(isp)
        nlte_data(iz)%i(ii)%coll_nlin = coll_nlin
        allocate(nlte_data(iz)%i(ii)%llw_coll(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%lup_coll(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%C0(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%C1(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%C2(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%C3(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%n_coll(coll_nlin))
        allocate(nlte_data(iz)%i(ii)%delE(coll_nlin))
        nlte_data(iz)%i(ii)%llw_coll = 0
        nlte_data(iz)%i(ii)%lup_coll = 0
        nlte_data(iz)%i(ii)%C0 = 0d0
        nlte_data(iz)%i(ii)%C1 = 0d0
        nlte_data(iz)%i(ii)%C2 = 0d0
        nlte_data(iz)%i(ii)%C3 = 0d0
        nlte_data(iz)%i(ii)%n_coll = 0d0
        nlte_data(iz)%i(ii)%delE = 0d0
        !-- PI/RR transition data
        allocate(nlte_data(iz)%i(ii)%lev_PI(nlev))
        allocate(nlte_data(iz)%i(ii)%C0_PI(nlev))
        allocate(nlte_data(iz)%i(ii)%C1_PI(nlev))
        allocate(nlte_data(iz)%i(ii)%C2_PI(nlev))
        allocate(nlte_data(iz)%i(ii)%C3_PI(nlev))
        allocate(nlte_data(iz)%i(ii)%n_PI(nlev))
        allocate(nlte_data(iz)%i(ii)%C0_RR(nlev))
        allocate(nlte_data(iz)%i(ii)%C1_RR(nlev))
        allocate(nlte_data(iz)%i(ii)%C2_RR(nlev))
        allocate(nlte_data(iz)%i(ii)%C3_RR(nlev))
        allocate(nlte_data(iz)%i(ii)%n_RR(nlev))
        allocate(nlte_data(iz)%i(ii)%delE_PI(nlev))
        nlte_data(iz)%i(ii)%lev_PI = 0
        nlte_data(iz)%i(ii)%C0_PI = 0d0
        nlte_data(iz)%i(ii)%C1_PI = 0d0
        nlte_data(iz)%i(ii)%C2_PI = 0d0
        nlte_data(iz)%i(ii)%C3_PI = 0d0
        nlte_data(iz)%i(ii)%n_PI = 0d0
        nlte_data(iz)%i(ii)%C0_RR = 0d0
        nlte_data(iz)%i(ii)%C1_RR = 0d0
        nlte_data(iz)%i(ii)%C2_RR = 0d0
        nlte_data(iz)%i(ii)%C3_RR = 0d0
        nlte_data(iz)%i(ii)%n_RR = 0d0
        nlte_data(iz)%i(ii)%delE_PI = 0d0
       enddo !ii
      enddo !iz
c
      end subroutine nlte_alloc_data
c
c
c
      subroutine nlte_dealloc
c     ----------------------
      implicit none
************************************************************************
* deallocate the complex nlte datastructure
************************************************************************
      integer :: iz,ii,i
c
      if(nlte_nelem.gt.0) then
       deallocate(nlte_ndens)
       do iz=1,nlte_nelem
        do ii=1,iz
          deallocate(nlte_data(iz)%i(ii)%chilev)
          deallocate(nlte_data(iz)%i(ii)%levid)
          deallocate(nlte_data(iz)%i(ii)%glev)
          deallocate(nlte_data(iz)%i(ii)%nplev)
          deallocate(nlte_data(iz)%i(ii)%wl0)
          deallocate(nlte_data(iz)%i(ii)%f)
          deallocate(nlte_data(iz)%i(ii)%gl)
          deallocate(nlte_data(iz)%i(ii)%gu)
          deallocate(nlte_data(iz)%i(ii)%llw)
          deallocate(nlte_data(iz)%i(ii)%lup)
          deallocate(nlte_data(iz)%i(ii)%Aul)
          deallocate(nlte_data(iz)%i(ii)%Bul)
          deallocate(nlte_data(iz)%i(ii)%Blu)
          deallocate(nlte_data(iz)%i(ii)%qlu)
          deallocate(nlte_data(iz)%i(ii)%llw_coll)
          deallocate(nlte_data(iz)%i(ii)%lup_coll)
          deallocate(nlte_data(iz)%i(ii)%C0)
          deallocate(nlte_data(iz)%i(ii)%C1)
          deallocate(nlte_data(iz)%i(ii)%C2)
          deallocate(nlte_data(iz)%i(ii)%C3)
          deallocate(nlte_data(iz)%i(ii)%n_coll)
          deallocate(nlte_data(iz)%i(ii)%delE)
          deallocate(nlte_data(iz)%i(ii)%lev_PI)
          deallocate(nlte_data(iz)%i(ii)%C0_PI)
          deallocate(nlte_data(iz)%i(ii)%C1_PI)
          deallocate(nlte_data(iz)%i(ii)%C2_PI)
          deallocate(nlte_data(iz)%i(ii)%C3_PI)
          deallocate(nlte_data(iz)%i(ii)%n_PI)
          deallocate(nlte_data(iz)%i(ii)%C0_RR)
          deallocate(nlte_data(iz)%i(ii)%C1_RR)
          deallocate(nlte_data(iz)%i(ii)%C2_RR)
          deallocate(nlte_data(iz)%i(ii)%C3_RR)
          deallocate(nlte_data(iz)%i(ii)%n_RR)
          deallocate(nlte_data(iz)%i(ii)%delE_PI)
        enddo !ii
        deallocate(nlte_data(iz)%i)
       enddo !iz
       deallocate(nlte_data)
      endif
      end subroutine nlte_dealloc
c
      end module nltemod
c vim: fdm=marker
