*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
*See LANL_COPYING and LANL_README for details of LANL copyright assertion.
      module inputstrmod
c     ------------------
      implicit none
************************************************************************
* Supernova atmospheric stratification
************************************************************************
      character(9),private :: fname='input.str'
      integer :: str_nabund=0
      logical :: str_lradtemp=.false.,str_lvol=.false.
      logical :: str_lye=.false.
      logical :: str_llum=.false.
      integer,allocatable :: str_iabund(:) !(nabund)
c
      real*8,allocatable :: str_xleft(:) !(nx+1)
      real*8,allocatable :: str_yleft(:) !(ny+1)
      real*8,allocatable :: str_zleft(:) !(nz+1)
      real*8,allocatable :: str_vxleft(:) !(nx+1)
      real*8,allocatable :: str_vyleft(:) !(ny+1)
      real*8,allocatable :: str_vzleft(:) !(nz+1)
      real*8,allocatable :: str_mass(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_vol(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_temp(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_radtemp(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_lum(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_ye(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_massfr(:,:,:,:) !(nabund,nx,ny,nz)
c
c-- domain compression
      logical :: str_lvoid=.false.  !flag existence of void cells
      integer :: str_nc=0  !number of cells in compressed grid
      integer,allocatable :: str_idcell(:) !(nc)
      real*8,allocatable :: str_massdc(:) !(nc)
      real*8,allocatable :: str_voldc(:) !(nc)
      real*8,allocatable :: str_massfrdc(:,:) !(nabund,nc)
      real*8,allocatable :: str_tempdc(:) !(nc)
      real*8,allocatable :: str_radtempdc(:) !(nc)
      real*8,allocatable :: str_lumdc(:) !(nc)
      real*8,allocatable :: str_yedc(:) !(nc)
c
c-- domain decomposition
      real*8,allocatable :: str_massdd(:) !(gas_ncell)
      real*8,allocatable :: str_voldd(:) !(gas_ncell)
      real*8,allocatable :: str_massfrdd(:,:) !(nabund,gas_ncell)
      real*8,allocatable :: str_tempdd(:) !(gas_ncell)
      real*8,allocatable :: str_radtempdd(:) !(gas_ncell)
      real*8,allocatable :: str_lumdd(:) !(gas_ncell)
      real*8,allocatable :: str_yedd(:) !(gas_ncell)
c
      character(8),allocatable,private :: str_abundlabl(:) !(nabund)
c
      integer,private :: nx,ny,nz
      integer,private :: igeom
c
      save
      public
c
      contains
c
c
c
      subroutine inputstr_dealloc
c     ---------------------------
      implicit none
      if(allocated(str_iabund)) deallocate(str_iabund)
      deallocate(str_xleft,str_yleft,str_zleft)
      deallocate(str_vxleft,str_vyleft,str_vzleft)
      deallocate(str_idcell)
      deallocate(str_massdc,str_massdd)
      if(str_nabund>0) then
       deallocate(str_massfrdc,str_massfrdd)
       if(allocated(str_abundlabl)) deallocate(str_abundlabl) !only on impi0
      endif
      str_nabund=0
      deallocate(str_tempdc,str_tempdd)
      if(str_lvol) deallocate(str_voldc,str_voldd)
      if(str_lradtemp) deallocate(str_radtempdc,str_radtempdd)
      if(str_llum) deallocate(str_lumdc,str_lumdd)
      if(str_lye) deallocate(str_yedc,str_yedd)
      end subroutine inputstr_dealloc
c
c
c
      subroutine read_inputstr(igeomin,ndim,nmpi)
c     --------------------------------------------------------
      use physconstmod
      use miscmod
      implicit none
      integer,intent(in) :: igeomin,nmpi
      integer,intent(in) :: ndim(3)
************************************************************************
* Read the input structure file
************************************************************************
      integer :: i,j,k,l,ierr,nx_r,ny_r,nz_r,nvar,ncol
      integer :: jmass,jxleft,jye,jtemp
      integer :: jradtemp,jvol,jvxleft,jlum
      integer :: ncorner,nvoid,ncell,ncpr
      character(2) :: dmy
      character(8),allocatable :: labl(:)
      real*8,allocatable :: raw(:,:)
      real*8 :: mni56,help
      real*8 :: r,rs
c
c-- copy
      igeom = igeomin
      nx = ndim(1)
      ny = ndim(2)
      nz = ndim(3)
c
c-- open file
      open(4,file=fname,status='old',iostat=ierr)
      if(ierr/=0) stop 'read_inputstr: file missing: input.str'
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy, nx_r,ny_r,nz_r,ncol,str_nabund
      if(ierr/=0) stop 'read_inputstr: input.str fmt err: dimensions'
c-- verify dimension
      if(nx_r/=nx) stop 'read_inputstr: incompatible nx dimension'
      if(ny_r/=ny) stop 'read_inputstr: incompatible ny dimension'
      if(nz_r/=nz) stop 'read_inputstr: incompatible nz dimension'
c
c-- allocate label arrays
      nvar = ncol - str_nabund
      allocate(str_abundlabl(str_nabund))
      allocate(labl(nvar))
c
c-- read labels
      read(4,*,iostat=ierr) dmy, labl, str_abundlabl
      if(ierr/=0) stop 'read_inputstr: input.str fmt err: col labels'
c
c-- var pointers
      jxleft = 0
      jvxleft =0
      jmass = 0
      jvol = 0
      jye = 0
      jtemp = 0
      jradtemp = 0
      jlum = 0
      do i=1,nvar
       if(lcase(trim(labl(i)))=='x_left') jxleft = i
       if(lcase(trim(labl(i)))=='vx_left') jvxleft = i
       if(lcase(trim(labl(i)))=='mass') jmass = i
       if(lcase(trim(labl(i)))=='vol') jvol = i
       if(lcase(trim(labl(i)))=='ye') jye = i
       if(lcase(trim(labl(i)))=='temp') jtemp = i
       if(lcase(trim(labl(i)))=='rad_temp') jradtemp = i
       if(lcase(trim(labl(i)))=='lum') jlum = i
      enddo
      if(jmass==0) stop 'read_inputstr: mass label not found'
      if(jtemp==0) stop 'read_inputstr: temp label not found'
      if(jvol>0) str_lvol = .true.
      if(jradtemp>0) str_lradtemp = .true.
      if(jlum>0) str_llum = .true.
      if(jye>0) str_lye = .true.
c
c-- allocate data arrays
      allocate(str_xleft(nx+1))
      allocate(str_yleft(ny+1))
      allocate(str_zleft(nz+1))
      allocate(str_vxleft(nx+1))
      allocate(str_vyleft(ny+1))
      allocate(str_vzleft(nz+1))
      allocate(str_mass(nx,ny,nz))
      allocate(str_massfr(str_nabund,nx,ny,nz))
      allocate(str_temp(nx,ny,nz))
      if(str_lvol) allocate(str_vol(nx,ny,nz))
      if(str_lradtemp) allocate(str_radtemp(nx,ny,nz))
      if(str_llum) allocate(str_lum(nx,ny,nz))
      if(str_lye) allocate(str_ye(nx,ny,nz))
      allocate(raw(ncol,nx*ny*nz))
c
c-- read body
      read(4,*,iostat=ierr) raw
      if(ierr/=0) stop 'read_inputstr: input.str format err: body'
      read(4,*,iostat=ierr) dmy
      if(ierr/=-1) stop 'read_inputstr: input.str body too long'
      close(4)
c
c-- validity check
      if(any(raw/=raw)) stop 'read_inputstr: nan in input'
c
c-- transer data to final arrays
c-- dim 1
      if(igeom==11 .and. jxleft==0 .and. jvxleft==0
     &   .and. raw(1,1)/=0d0) then
c-- set inner boundary at 0 if only right bounds are provided
       str_xleft(1) = 0d0
       str_xleft(2:) = raw(1,:nx)
       str_vxleft(1) = 0d0
       str_vxleft(2:) = raw(2,:nx)
      else
c-- set inner boundary to first left bound
       str_xleft(1) = raw(1,1)
       str_xleft(2:) = raw(2,:nx)
       str_vxleft(1) = raw(jvxleft,1)
       str_vxleft(2:) = raw(jvxleft+1,:nx)
      endif
c-- dim 2
c-- set polar projection
      str_yleft = [-1d0,1d0]
c
c-- dim 3
c-- set azimuthal angle
      str_zleft = [0d0,pc_pi2]
c
c-- set non-radial velocity components to 0
      if(ny==1) str_vyleft = 0
      if(nz==1) str_vzleft = 0
c
c-- check grid monotonicity
      help = str_xleft(1)
      do i=2,nx+1
       if(str_xleft(i)<=help) stop
     &   'read_inputstr: x grid not increasing'
       help = str_xleft(i)
      enddo
c
c-- vars
      l = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
       l = l+1
       str_mass(i,j,k) = raw(jmass,l)
       str_temp(i,j,k)=raw(jtemp,l)
c-- if available in the input.str file
       if(str_lvol) str_vol(i,j,k) = raw(jvol,l)
       if(str_lradtemp) str_radtemp(i,j,k)=raw(jradtemp,l)
       if(str_llum) str_lum(i,j,k)=raw(jlum,l)
       if(str_lye) str_ye(i,j,k)=raw(jye,l)
      enddo
      enddo
      enddo
c-- sanity check
      if(any(str_mass<0d0)) stop 'read_inputstr: mass<0'
c
c-- abundances
      l = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
       l = l+1
       str_massfr(:,i,j,k) = raw(nvar+1:,l)
      enddo
      enddo
      enddo
c-- sanity check
      if(any(str_massfr<0d0)) stop 'read_inputstr: massfr<0'
c
      deallocate(raw,labl)
c
c-- convert abundlabl to element codes
      call elnam2elcode
c-- check codes are unique, i.e. no duplicate columns
      do i=1,str_nabund-1
       do j=i+1,str_nabund
        if(str_iabund(j)==str_iabund(i)) stop
     &    'read_inputstr: duplicate abund columns'
       enddo
      enddo
c
c
c-- count valid cells
      ncell = count(str_mass>0d0)
      if(ncell/=nx*ny*nz) ncell = ncell+1
      ncpr = ceiling(ncell/float(nmpi))
c
c
c-- output
      write(6,*)
      write(6,*) 'input structure:'
      write(6,*) '===================='
      write(6,*) 'igeom :',igeom
      write(6,*) 'ndim  :',nx,ny,nz
      write(6,*) 'ncell :',ncell
      write(6,*) 'nc/rnk:',ncpr,ncpr*nmpi-ncell,ncell/float(ncpr*nmpi)
      write(6,*) 'mass  :',sngl(sum(str_mass)/pc_msun), 'Msun'
      write(6,*)
c
      end subroutine read_inputstr
c
c
c
      subroutine inputstr_compress
c     ----------------------------
      implicit none
************************************************************************
* put valid (non-void) cells in sequence, link the other (void) cells
* to the dummy cell at the end of the sequence.
************************************************************************
      integer :: i,j,k,l
      integer :: idcell
c
      str_nc = count(str_mass>0d0)
c
c-- add void cell
      if(str_nc/=nx*ny*nz) then
       str_lvoid = .true.
       str_nc = str_nc+1
      endif
c
      allocate(str_idcell(str_nc))
      allocate(str_massdc(str_nc))

      if(str_nabund>0) then
       allocate(str_massfrdc(str_nabund,str_nc))
       str_massfrdc = 0d0
      endif
c
      allocate(str_tempdc(str_nc))
      if(str_lvol) allocate(str_voldc(str_nc))
      if(str_lradtemp) allocate(str_radtempdc(str_nc))
      if(str_llum) allocate(str_lumdc(str_nc))
      if(str_lye) allocate(str_yedc(str_nc))
c-- zero all, including the dummy cell
      str_idcell = 0
      str_massdc = 0d0
c-- void temp [K]
      str_tempdc = 1000d0
      if(str_lradtemp) str_radtempdc = 1000d0
      if(str_llum) str_lumdc = 0.d0
      if(str_lye) str_yedc = .5d0
c
      l = 0
      idcell = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
       idcell = idcell+1
c-- skip void cells
       if(str_mass(i,j,k)<=0d0) cycle
c-- insert
       l = l+1
       str_idcell(l) = idcell
       str_massdc(l) = str_mass(i,j,k)
       if(str_nabund>0) str_massfrdc(:,l) = str_massfr(:,i,j,k)
c
       str_tempdc(l) = str_temp(i,j,k)
       if(str_lvol) str_voldc(l) = str_vol(i,j,k)
       if(str_lradtemp) str_radtempdc(l) = str_radtemp(i,j,k)
       if(str_llum) str_lumdc(l) = str_lum(i,j,k)
       if(str_lye) str_yedc(l) = str_ye(i,j,k)
      enddo !i
      enddo !j
      enddo !k
c-- sanity check
      if(str_lvoid) l = l+1
      if(l/=str_nc) stop 'inputstr_compress: l/=str_nc' !one dummy cell
      if(idcell/=nx*ny*nz) stop 'inputstr_compress: idcell/=nx*ny*nz'
c
c-- deallocate full grid
      deallocate(str_mass)
      if(allocated(str_vol)) deallocate(str_vol)
      if(allocated(str_radtemp)) deallocate(str_radtemp)
      if(allocated(str_lum)) deallocate(str_lum)
      if(allocated(str_massfr)) deallocate(str_massfr)
      if(allocated(str_temp)) deallocate(str_temp)
c
      end subroutine inputstr_compress
c
c
c
      subroutine elnam2elcode
c     ------------------------------
      use miscmod, only:lcase
      use gasmod
      use elemdatamod
      implicit none
************************************************************************
* convert the abundlabl labels to element codes (atomic z number), which
* also serve as indexing pointers to the mass0fr array.
************************************************************************
      integer :: l,j,iabund
      character(4) :: elname
c
c-- allocate element code array (pointer to mass0fr)
      allocate(str_iabund(str_nabund))
c
c-- determine atomic z number
      do l=1,str_nabund
       iabund = 0
       elname = lcase(trim(str_abundlabl(l)))
c-- search elem_data for corresponding entry
       do j=1,elem_neldata
        if(lcase(trim(elem_data(j)%sym))==elname) exit
       enddo
c-- verify hit
       if(j>elem_neldata) then
        write(*,*) 'unknown chemical element name:',elname
        stop 'elnam2elcode: no such element found in elemdata'
       endif
       iabund = j
      !write(6,*) 'el found: ',elname,iabund
c
c-- store element code (pointer to mass0fr)
       str_iabund(l) = iabund
       enddo
      end subroutine elnam2elcode
c
      end module inputstrmod
c vim: fdm=marker
