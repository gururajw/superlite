*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module mpimod
c     -------------
      implicit none
      INCLUDE 'mpif.h'
c
      integer,parameter :: impi0=0 !the master rank
      logical :: lmpi0 !true for the master rank
      integer :: impi !mpi rank
      integer :: nmpi !number of mpi tasks
      integer,private :: ierr
c
      integer,private,allocatable :: counts(:)
      integer,private,allocatable :: displs(:)
c
      save
c
      contains
c
c
c
      subroutine bcast_permanent
c     --------------------------
      use inputparmod
      use inputstrmod
      use sourcemod
      use ionsmod
      use ffxsmod
      use bfxsmod
      use bbxsmod
      use nltemod
      use gridmod
      use gasmod
      use groupmod
      use fluxmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
************************************************************************
      integer :: i,n
      integer :: il,ii,ir,ic
      integer :: nx,ny,nz
      integer :: ntemp
      logical,allocatable :: lsndvec(:)
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
      character(4),allocatable :: csndvec(:)
c
c-- inputparmod variables
c========================
c-- create pointer arrays
      call inputpar_create_pointers(il,ii,ir,ic)
c-- broadcast logicals
      n = il
      allocate(lsndvec(n))
      forall(i=1:n) lsndvec(i) = in_l(i)%p
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_l(i)%p = lsndvec(i)
      deallocate(lsndvec)
c-- broadcast integers
      n = ii
      allocate(isndvec(n))
      forall(i=1:n) isndvec(i) = in_i(i)%p
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_i(i)%p = isndvec(i)
      deallocate(isndvec)
c-- broadcast real*8
      n = ir
      allocate(sndvec(n))
      forall(i=1:n) sndvec(i) = in_r(i)%p
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_r(i)%p = sndvec(i)
      deallocate(sndvec)
c-- broadcast characters
      n = ic
      allocate(csndvec(n))
      forall(i=1:n) csndvec(i) = in_c(i)%p
      call mpi_bcast(csndvec,n*4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
!     forall(i=1:n) in_c(i)%p = csndvec(i) !avoid gfortran 4.6.3 compiler bug
      do i=1,n
       in_c(i)%p = csndvec(i)
      enddo
      deallocate(csndvec)
c
c-- set number of threads, i.e. non-automatic
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c
c-- everything else
c==================
c-- broadcast constants
c-- logical
      n = 5
      allocate(lsndvec(n))
      if(lmpi0) lsndvec = [str_lvoid,str_lvol,str_lye,str_lradtemp,
     &                     str_llum]
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      str_lvoid = lsndvec(1)
      str_lvol = lsndvec(2)
      str_lye = lsndvec(3)
      str_lradtemp = lsndvec(4)
      str_llum = lsndvec(5)
      deallocate(lsndvec)
c
c-- integer
      n = 8
      allocate(isndvec(n))
      if(lmpi0) ntemp = size(zeta_temp)
      if(lmpi0) isndvec = [ion_nion,ion_iionmax,bb_nline,
     &  nlte_nelem,nlte_nspecies,ntemp,str_nc,str_nabund]
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      ion_nion       = isndvec(1)
      ion_iionmax    = isndvec(2)
      bb_nline       = isndvec(3)
      nlte_nelem     = isndvec(4)
      nlte_nspecies  = isndvec(5)
      ntemp          = isndvec(6)
      str_nc         = isndvec(7)
      str_nabund     = isndvec(8)
      deallocate(isndvec)

c
c-- dimenstions
      nx = in_ndim(1)
      ny = in_ndim(2)
      nz = in_ndim(3)
c
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       if(bb_nline>0) allocate(bb_xs(bb_nline))
       allocate(str_xleft(nx+1))
       allocate(str_yleft(ny+1))
       allocate(str_zleft(nz+1))
       allocate(str_vxleft(nx+1))
       allocate(str_vyleft(ny+1))
       allocate(str_vzleft(nz+1))
       allocate(str_idcell(str_nc))
       if(str_nabund>0) allocate(str_iabund(str_nabund))
      endif
c
c-- inputstr
      call mpi_bcast(str_xleft,nx+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_yleft,ny+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_zleft,nz+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_vxleft,nx+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_vyleft,ny+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_vzleft,nz+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_idcell,str_nc,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      if(str_nabund>0) then
       call mpi_bcast(str_iabund,str_nabund,MPI_INTEGER,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c
c-- broadcast data
c-- bound-bound
      if(bb_nline>0) then
       call mpi_bcast(bb_xs,sizeof(bb_xs),MPI_BYTE,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c-- bound-free
      call mpi_bcast(bf_ph1,6*7*30*30,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(bf_ph2,7*30*30,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- free-free
      call mpi_bcast(ff_gff,ff_nu*ff_ngg,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call bcast_ions
      if(in_nlte) call bcast_nlte
c
c
      contains
c
      subroutine bcast_ions
c     ---------------------
      implicit none
************************************************************************
* broadcast the ions data structure
************************************************************************
      integer :: ii,iz,iion,n
      real*8 :: vec(11000)
      logical :: vecl(11000)
      integer :: nion(gas_nelem)
      integer :: nlev(ion_nion)
      real*8 :: e(ion_nion)
c
c-- evaluate shape info
      if(lmpi0) then
       iion = 0
       do iz=1,gas_nelem
        nion(iz) = ion_el(iz)%ni
        do ii=1,ion_el(iz)%ni
         iion = iion + 1
         nlev(iion) = ion_el(iz)%i(ii)%nlev
         e(iion) = ion_el(iz)%i(ii)%e
        enddo !ii
       enddo !iz
c-- sanity check
       if(iion/=ion_nion) stop "bcast_perm: ion_nion problem"
      endif
c
c-- bcast shape info and allocate
      call mpi_bcast(nion,gas_nelem,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nlev,ion_nion,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- allocate structure
      if(impi/=impi0) then
        call ions_alloc_el(gas_nelem,nion,ion_nion,nlev,ntemp)
      endif
c-- fill structure
      call mpi_bcast(e,ion_nion,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      iion = 0
      do iz=1,gas_nelem
       do ii=1,ion_el(iz)%ni
        iion = iion + 1
        n = nlev(iion)
c-- eion
        if(impi/=impi0) ion_el(iz)%i(ii)%e = e(iion)
c-- elev
        if(lmpi0) vec(:n) = ion_el(iz)%i(ii)%elev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%elev = vec(:n)
c-- glev
        if(lmpi0) vec(:n) = ion_el(iz)%i(ii)%glev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%glev = vec(:n)
c-- meta !nebular approximation
        if(lmpi0) vecl(:n) = ion_el(iz)%i(ii)%meta
        call mpi_bcast(vecl(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%meta = vecl(:n)
c-- zval !nebular approximation
        if(lmpi0) vec(:ntemp) = ion_el(iz)%i(ii)%zval
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%zval = vec(:n)
       enddo !ii
      enddo !iz
c-- temperature array for zeta values !nebular approximation
      call mpi_bcast(zeta_temp,ntemp,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c-- sanity check
      if(iion/=ion_nion) stop "bcast_perm: ion_nion problem"
c
      end subroutine bcast_ions

      subroutine bcast_nlte
      implicit none
************************************************************************
* broadcast the nlte data structure
************************************************************************
      integer :: ii,iz,isp,n
      real*8 :: vec(10000)
      integer :: vecn(10000)
      integer :: nlin(nlte_nspecies),nlev(nlte_nspecies) ! currently 1 species - h1
      integer :: coll_nlin(nlte_nspecies) ! currently 1 species - h1

c-- evaluate shape info
      if(lmpi0) then
c-- count nlte species
       isp = 0
       do iz=1,nlte_nelem
        do ii=1,iz
         isp = isp + 1
         nlev(isp) = nlte_data(iz)%i(ii)%nlev
         nlin(isp) = nlte_data(iz)%i(ii)%nlin
         coll_nlin(isp) = nlte_data(iz)%i(ii)%coll_nlin
        enddo !ii
       enddo !iz
c-- sanity check
       if(isp/=nlte_nspecies) stop "bcast_nlte: nlte_nspecies problem"
      endif !lmpi0
c-- broadcast shape info and allocate
      call mpi_bcast(nlev,nlte_nspecies,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nlin,nlte_nspecies,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(coll_nlin,nlte_nspecies,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- allocate structure
      if(impi/=impi0) then
       call nlte_alloc_data(nlte_nelem,nlte_nspecies,nlev,nlin,
     &   coll_nlin)
      endif
c-- fill structure
      isp = 0
      do iz=1,nlte_nelem
       do ii=1,iz
        isp = isp + 1
c-- level data
        n = nlev(isp)
c-- chilev
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%chilev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%chilev = vec(:n)
c-- glev
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%levid
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%levid = vecn(:n)
c-- glev
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%glev
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%glev = vecn(:n)
c-- nplev
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%nplev
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%nplev = vecn(:n)
c-- radiative line data
        n = nlin(isp)
c-- wl0
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%wl0
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%wl0 = vec(:n)
c-- f
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%f
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%f = vec(:n)
c-- gl
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%gl
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%gl = vecn(:n)
c-- gu
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%gu
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%gu = vecn(:n)
c-- llw
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%llw
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%llw = vecn(:n)
c-- lup
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%lup
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%lup = vecn(:n)
c-- Aul
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%Aul
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%Aul = vec(:n)
c-- Bul
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%Bul
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%Bul = vec(:n)
c-- Blu
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%Blu
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%Blu = vec(:n)
c-- qlu
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%qlu
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%qlu = vec(:n)
c-- EIE line data
        n = coll_nlin(isp)
c-- llw
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%llw_coll
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%llw_coll = vecn(:n)
c-- lup
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%lup_coll
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%lup_coll = vecn(:n)
c-- Fit coefficients
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C0
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C0 = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C1
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C1 = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C2
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C2 = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C3
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C3 = vec(:n)
c-- pre-factor coefficient
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%n_coll
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%n_coll = vec(:n)
c-- transition energy
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%delE
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%delE = vec(:n)
c-- PI & RR transition data
c-- level
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%lev_PI
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%lev_PI = vecn(:n)
c-- Fit coefficients for PI
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C0_PI
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C0_PI = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C1_PI
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C1_PI = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C2_PI
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C2_PI = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C3_PI
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C3_PI = vec(:n)
c-- pre-factor coefficient for PI
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%n_PI
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%n_PI = vecn(:n)
c-- Fit coefficients for RR
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C0_RR
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C0_RR = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C1_RR
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C1_RR = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C2_RR
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C2_RR = vec(:n)
c
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%C3_RR
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%C3_RR = vec(:n)
c-- pre-factor coefficient for RR
        if(lmpi0) vecn(:n) = nlte_data(iz)%i(ii)%n_RR
        call mpi_bcast(vecn(1),n,MPI_INTEGER,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%n_PI = vecn(:n)
c-- transition energy for RR/PI
        if(lmpi0) vec(:n) = nlte_data(iz)%i(ii)%delE_PI
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) nlte_data(iz)%i(ii)%delE_PI = vec(:n)
       enddo !ii
      enddo !iz
c-- sanity check
       if(isp/=nlte_nspecies) stop "bcast_nltem: nlte_nspecies problem"
c
      end subroutine bcast_nlte
c
c
      end subroutine bcast_permanent
c
c
c
      subroutine scatter_inputstruct(ndim,icell1,ncell)
c     -------------------------------------------------
      use inputstrmod
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3)
      integer,intent(out) :: icell1,ncell
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
************************************************************************
      integer :: i,n,nc,nx,ny,nz
c
      nx = ndim(1)
      ny = ndim(2)
      nz = ndim(3)
c
c-- calculate offsets
      allocate(counts(nmpi),displs(nmpi))
      displs(1) = 0
      nc = str_nc
      do i=1,nmpi
       n = ceiling(nc/(nmpi-i+1d0))
       counts(i) = n
       if(i-1==impi) ncell = n
       if(i-1==impi) icell1 = displs(i) + 1
       if(i<nmpi) displs(i+1) = displs(i) + n
       nc = nc - n
      enddo
      if(sum(counts)/=str_nc) stop 'scatter_inputstr: counts/=str_nc'
      if(nc/=0) stop 'scatter_inputstr: nc/=0'
c
c-- allocate domain decomposed and domain compressed
      if(impi/=impi0) allocate(str_massdc(str_nc))
      allocate(str_massdd(ncell))
      call mpi_scatterv(str_massdc,counts,displs,MPI_REAL8,
     &  str_massdd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c -- volume
      if(str_lvol) then
        if(impi/=impi0) allocate(str_voldc(str_nc))
        allocate(str_voldd(ncell))
        call mpi_scatterv(str_voldc,counts,displs,MPI_REAL8,
     &    str_voldd,ncell,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
      endif
c
c-- mass fractions if available
      if(str_nabund>0) then
       if(impi/=impi0) allocate(str_massfrdc(str_nabund,str_nc))
       allocate(str_massfrdd(str_nabund,ncell))
       n = str_nabund
       call mpi_scatterv(str_massfrdc,n*counts,n*displs,MPI_REAL8,
     &   str_massfrdd,n*ncell,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c
c-- gas temperature structure
!      if(str_ltemp) then
       if(impi/=impi0) allocate(str_tempdc(str_nc))
       allocate(str_tempdd(ncell))
       call mpi_scatterv(str_tempdc,counts,displs,MPI_REAL8,
     &  str_tempdd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
!      endif
c-- radiation temperature structure if available
      if(str_lradtemp) then
       if(impi/=impi0) allocate(str_radtempdc(str_nc))
       allocate(str_radtempdd(ncell))
       call mpi_scatterv(str_radtempdc,counts,displs,MPI_REAL8,
     &  str_radtempdd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      endif
c-- luminosity structure if available
      if(str_llum) then
       if(impi/=impi0) allocate(str_lumdc(str_nc))
       allocate(str_lumdd(ncell))
       call mpi_scatterv(str_lumdc,counts,displs,MPI_REAL8,
     &  str_lumdd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      endif
c-- ye structure if available
      if(str_lye) then
       if(impi/=impi0) allocate(str_yedc(str_nc))
       allocate(str_yedd(ncell))
       call mpi_scatterv(str_yedc,counts,displs,MPI_REAL8,
     &  str_yedd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      endif

      end subroutine scatter_inputstruct
c
c
c
      subroutine allgather_gridvolume
c     -------------------------------
      use gridmod
      use gasmod
      implicit none
************************************************************************
* gather gas_vol to grd_vol
************************************************************************
      call mpi_allgatherv(gas_vol,gas_ncell,MPI_REAL8,
     &  grd_vol,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      end subroutine allgather_gridvolume
c
c
      subroutine allgather_jrad
      use gridmod
      use gasmod
      use groupmod
      implicit none
************************************************************************
* gather initial gas_jrad to grd_jrad
************************************************************************
      call mpi_allgatherv(gas_jrad,grp_ng*gas_ncell,MPI_REAL8,
     &  grd_jrad,grp_ng*counts,grp_ng*displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      end subroutine allgather_jrad
c
c
      subroutine allgather_radtemp
c     -------------------------------
      use gridmod
      use gasmod
      implicit none
************************************************************************
* gather initial gas_radtemp to grd_radtemp
************************************************************************
      call mpi_allgatherv(gas_radtemp,gas_ncell,MPI_REAL8,
     &  grd_radtemp,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_lum,gas_ncell,MPI_REAL8,
     &  grd_lum,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      end subroutine allgather_radtemp
c
c
      subroutine bcast_nonpermanent
c     -----------------------------
      use gridmod
      use gasmod
      use groupmod
      use sourcemod
      use totalsmod
      use particlemod
      use timingmod
      use transportmod, only:trn_noampfact
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
************************************************************************
      real*8 :: t0,t1
      real*8 :: sndgas(gas_ncell)
      real*8 :: sndgrd(grd_ncell)
      integer :: nvacant
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
c-- gather
      sndgas = 1d0/gas_temp
      call mpi_allgatherv(sndgas,gas_ncell,MPI_REAL8,
     &  grd_tempinv,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_capgrey,gas_ncell,MPI_REAL8,
     &  grd_capgrey,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_capemitgrey,gas_ncell,MPI_REAL8,
     &  grd_capemitgrey,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_capross,gas_ncell,MPI_REAL8,
     &  grd_capross,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_natom,gas_ncell,MPI_REAL8,
     &  grd_natom,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_nelec,gas_ncell,MPI_REAL8,
     &  grd_nelec,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_emit,gas_ncell,MPI_REAL8,
     &  grd_emit,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_sig,gas_ncell,MPI_REAL8,
     &  grd_sig,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_cap,grp_ng*gas_ncell,MPI_REAL,
     &  grd_cap,grp_ng*counts,grp_ng*displs,MPI_REAL,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_capemit,grp_ng*gas_ncell,MPI_REAL,
     &  grd_capemit,grp_ng*counts,grp_ng*displs,MPI_REAL,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_emiss,grp_ng*gas_ncell,MPI_REAL,
     &  grd_emiss,grp_ng*counts,grp_ng*displs,MPI_REAL,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_jrad,grp_ng*gas_ncell,MPI_REAL8,
     &  grd_jrad,grp_ng*counts,grp_ng*displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_lum,gas_ncell,MPI_REAL8,
     &  grd_lum,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
c-- broadcast
      call mpi_bcast(tot_esurf,1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- allreduce
      if(.not.trn_noampfact) then
       sndgrd = grd_eamp
       call mpi_allreduce(sndgrd,grd_eamp,grd_ncell,MPI_REAL8,MPI_SUM,
     &   MPI_COMM_WORLD,ierr)
      endif
c
c-- allgather
      nvacant = count(prt_isvacant) !count returns the same type as prt_isvacant
      call mpi_allgather(nvacant,1,MPI_INTEGER8,
     &  src_nvacantall,1,MPI_INTEGER8,
     &  MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpibcast, t1-t0)
c
      end subroutine bcast_nonpermanent
c
c
      subroutine allgather_leakage
c     ------------------------------------------
      use gridmod
      use timingmod
      implicit none
************************************************************************
      integer :: n
      real*8,allocatable :: snd(:,:)
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
      allocate(snd(10,grd_ndd))
      snd = grd_opaclump(:,grd_idd1:grd_idd1+grd_ndd-1)
      call mpi_allgatherv(snd,11*grd_ndd,MPI_REAL8,
     &  grd_opaclump,11*counts,11*displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
      n = grd_nep
      allocate(snd(n,grd_ndd))
      snd = grd_emitprob(:,grd_idd1:grd_idd1+grd_ndd-1)
      call mpi_allgatherv(snd,n*grd_ndd,MPI_REAL8,
     &  grd_emitprob,n*counts,n*displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
      t1 = t_time()
      call timereg(t_mpimisc, t1-t0)
c
      end subroutine allgather_leakage
c
c
c
      subroutine reduce_gridtally
c     -----------------------
      use gridmod,nx=>grd_nx,ny=>grd_ny,nz=>grd_nz
      use gasmod
      use timingmod
      use fluxmod
      use sourcemod
      implicit none
************************************************************************
* Reduce the results from particle_advance that are needed for the
* temperature correction.
************************************************************************
      integer :: n
      integer :: isnd2f(flx_nmu,flx_nom)
      real*8 :: snd2f(flx_nmu,flx_nom)
      integer,allocatable :: isnd(:)
      real*8,allocatable :: snd(:)
      real*8,allocatable :: snd2(:,:)
      real*8 :: help
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
c-- dim==3
      allocate(isnd(grd_ncell))
      n = grd_ncell
c
      isnd = grd_methodswap
      call mpi_reduce(isnd,grd_methodswap,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(isnd)
c
      allocate(snd(grd_ncell))
      snd = grd_tally
      call mpi_reduce(snd,grd_tally,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
c
c-- scatter
c
      snd = grd_tally(:)
      call mpi_scatterv(snd,counts,displs,MPI_REAL8,
     &  gas_eraddens,gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
c-- timing statistics
      help = t_pckt_stat(1)
      call mpi_reduce(help,t_pckt_stat(1),1,MPI_REAL8,MPI_MIN,
     &  impi0,MPI_COMM_WORLD,ierr)
      help = t_pckt_stat(2)/nmpi
      call mpi_reduce(help,t_pckt_stat(2),1,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      help = t_pckt_stat(3)
      call mpi_reduce(help,t_pckt_stat(3),1,MPI_REAL8,MPI_MAX,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpireduc, t1-t0)
c
      end subroutine reduce_gridtally
c
c
c
      subroutine allreduce_jrad
c     -----------------------
      use gridmod
      use groupmod
      use timingmod
      implicit none
************************************************************************
* Reduce the results from particle_advance for radiation intensity
************************************************************************
      integer :: n
      real*8,allocatable :: snd2(:,:)
      real*8 :: help
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
      allocate(snd2(grp_ng,grd_ncell))
      n=grp_ng*grd_ncell
      snd2 = grd_jrad
      call mpi_allreduce(snd2,grd_jrad,n,MPI_REAL8,MPI_SUM,
     &  MPI_COMM_WORLD,ierr)
      deallocate(snd2)
c
      t1 = t_time()
      call timereg(t_mpireduc, t1-t0)
c
      end subroutine allreduce_jrad
c
c
c
      subroutine scatter_jrad
      use gridmod
      use groupmod
      use gasmod
      implicit none
************************************************************************
* scatter the results for radiation intensity to all ranks
************************************************************************
      integer :: n
      real*8,allocatable :: snd(:),snd2(:,:)
      real*8 :: help
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
      allocate(snd(grd_ncell),snd2(grp_ng,grd_ncell))
      n=grp_ng
c
      snd = grd_radtemp
      call mpi_scatterv(snd,counts,displs,MPI_REAL8,
     &  gas_radtemp,gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      snd2 = grd_jrad
      call mpi_scatterv(snd2,n*counts,n*displs,MPI_REAL8,
     &  gas_jrad,n*gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(snd,snd2)
      end subroutine scatter_jrad
c
c
c
      subroutine reduce_srctally
c     ---------------------------
      use gridmod
      use groupmod
      use sourcemod
      implicit none
************************************************************************
* Reduce flux arrays for output
************************************************************************
      integer :: n
      integer :: isnd2f(grd_ncell,grp_ng)
      real*8 :: snd2f(grd_ncell,grp_ng)
c
c-- flux dim==3
      n = grd_ncell*grp_ng
      isnd2f = src_number
      call mpi_reduce(isnd2f,src_number,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd2f = src_energy
      call mpi_reduce(snd2f,src_energy,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      end subroutine reduce_srctally
c
c
c
      subroutine reduce_fluxtally
c     ---------------------------
      use fluxmod
      implicit none
************************************************************************
* Reduce flux arrays for output
************************************************************************
      integer :: n
      integer :: isnd3f(flx_ng,flx_nmu,flx_nom)
      real*8 :: snd3f(flx_ng,flx_nmu,flx_nom)
c
c-- flux dim==3
      n = flx_ng*flx_nmu*flx_nom
      isnd3f = flx_lumnum
      call mpi_reduce(isnd3f,flx_lumnum,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_luminos
      call mpi_reduce(snd3f,flx_luminos,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_lumdev
      call mpi_reduce(snd3f,flx_lumdev,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      end subroutine reduce_fluxtally
c
c
c
      subroutine reduce_totals
c     -------------------------------------
      use gridmod
      use totalsmod
      use gasmod
      use timingmod
      implicit none
************************************************************************
* for output
************************************************************************
      integer :: n
      real*8,allocatable :: sndvec(:),rcvvec(:)
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
c-- dim==0
      n = 4
      allocate(sndvec(n))
      allocate(rcvvec(n))
      sndvec = [tot_eout,tot_evelo,tot_sflux]
c
      call mpi_reduce(sndvec,rcvvec,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      if(lmpi0) then
         tot_eout = rcvvec(1)
         tot_evelo = rcvvec(2)
         tot_sflux = rcvvec(3)
      else
c-- zero out cumulative values on all other ranks to avoid double counting.
         tot_eout = 0d0
         tot_evelo = 0d0
      endif
      deallocate(sndvec)
      deallocate(rcvvec)
c
      t1 = t_time()
      call timereg(t_mpimisc, t1-t0)
c
      end subroutine reduce_totals
c
c
      subroutine reduce_gasprop
      use gridmod
      use totalsmod
      use gasmod
      use timingmod
      use inputparmod
      use inputstrmod
      implicit none
************************************************************************
* for output profiles
************************************************************************
      integer :: n
      logical :: lstat = .false.
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
      call mpi_gatherv(gas_mass,gas_ncell,MPI_REAL8,
     &   grd_mass,counts,displs,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_gatherv(gas_rho,gas_ncell,MPI_REAL8,
     &   grd_rho,counts,displs,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_gatherv(gas_ye,gas_ncell,MPI_REAL8,
     &   grd_ye,counts,displs,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c
      lstat = allocated(grd_massfr)
      if(.not.lstat) then
        allocate(grd_massfr(gas_nelem,grd_ncell))
      endif
      n = gas_nelem
      call mpi_gatherv(gas_massfr,n*gas_ncell,MPI_REAL8,
     &   grd_massfr,n*counts,n*displs,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpimisc, t1-t0)

      end subroutine reduce_gasprop

      subroutine mpimod_dealloc
      deallocate(counts,displs)
      end subroutine mpimod_dealloc
c
      end module mpimod
c vim: fdm=marker
