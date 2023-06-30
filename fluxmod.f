*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module fluxmod
c     ------------------
      implicit none

      integer :: flx_ndim(3)=0  !number of flux groups (1D,2D,3D)
      real*8 :: flx_wlmin,flx_wlmax
      integer :: flx_ng=0  !number of flux groups (1D,2D,3D)
      integer :: flx_nmu=0 !number of polar bins (2D,3D)
      integer :: flx_nom=0 !number of azimuthal bins (3D)
c
c-- wavelength (wl), polar (mu), and azimuthal (om) bins
      real*8,allocatable :: flx_wl(:) !(flx_ng)
      real*8,allocatable :: flx_mu(:) !(flx_nmu)
      real*8,allocatable :: flx_om(:) !(flx_nom)
!
!-- radiative flux
!=================
!-- outbound grouped luminosity
      real*8,allocatable :: flx_luminos(:,:,:) !(flx_ng,flx_nmu,flx_nom)
!-- sampled devation of group luminosity
      real*8,allocatable :: flx_lumdev(:,:,:) !(flx_ng,flx_nmu,flx_nom)
!-- number of escaped particles per group
      integer,allocatable :: flx_lumnum(:,:,:) !(flx_ng,flx_nmu,flx_nom)
c
      save
c
      contains
c
c
      subroutine fluxgrid_setup
c     -------------------------
      implicit none
*******************************************************
* calculate flux grid and tally arrays:
* -- if any dimensions of flux grid are negative, a
* file could be read for the corresponding grid using
* the absolute value of the negative dimension as the
* index.
* -- if a dimension number is positive, it could be
* treated as the grid size.
* -- for wavelength: if nflx(1)=0, then grp_wl could
* be used as the tally grid, flx_wl.
*******************************************************
      integer :: i
      character(12) :: fname
c
      fname = 'input.fluxwl'
c
c-- check if bins are read or generated
      do i = 1,3
         if(i==1.and.flx_ndim(i)<0) then
            call read_fluxgrid(i,flx_ndim(i),fname)
         else
            call generate_fluxgrid(i,flx_ndim(i),flx_wlmin,flx_wlmax)
         endif
      enddo
c
c-- allocate flux tally arrays
      allocate(flx_luminos(flx_ng,flx_nmu,flx_nom))
      allocate(flx_lumdev(flx_ng,flx_nmu,flx_nom))
      allocate(flx_lumnum(flx_ng,flx_nmu,flx_nom))
c
      end subroutine fluxgrid_setup
c
c
      subroutine read_fluxgrid(iflx,idex,fname)
c     ------------------------------
      use physconstmod
      implicit none
      integer,intent(in) :: idex,iflx
      character(12),intent(in) :: fname
*************************************************************
* read lab wavelength grid for flux from file
*************************************************************
      logical :: lexists
      integer :: l,n
      real*8 :: help
c
c-- sanity check
      if(iflx>3.or.iflx<1) stop 'read_fluxgrid: invalid 1st arg'

      inquire(file=fname,exist=lexists)
      if(.not.lexists) stop 'read_fluxgrid: missing file'
c
c-- grid lookup (rtw: this routine is merely wlgrid_setup)
      open(4,file=fname,status='old')
c
c-- strip header
      do l=1,3
       read(4,*)
      enddo
c-- read dimensions
      read(4,*) !not used
c
      read(4,*) n
c
      allocate(flx_wl(n+1))
      read(4,*) flx_wl(:)
      close(4)
c-- check monotonicity
      help = 0d0
      do l=1,n+1
       if(flx_wl(l)<=help) stop 'read_wlgrid: wlgrid not increasing'
       help = flx_wl(l)
      enddo
c
      write(6,*)
      write(6,*) 'read_fluxgrid: n=',n
!     write(6,*) 'wavelength grid in [cm]'
!     write(6,*) flx_wl
c
      end subroutine read_fluxgrid
c
c
      subroutine generate_fluxgrid(iflx,idex,wlmin,wlmax)
c     ---------------------------------------
      use physconstmod
      use groupmod
      implicit none
      integer,intent(in) :: idex,iflx
      real*8,intent(in) :: wlmin,wlmax
*************************************************************
* generate lab wavelength grid for flux
*************************************************************
      integer :: i
      real*8 :: help
c
c-- wavelength
      if(iflx==1) then
         if(idex==0) then
c-- set wl to transport grid
            flx_ng = grp_ng
            allocate(flx_wl(flx_ng+1))
            flx_wl = grp_wl
c
         elseif(idex>0) then
c-- logarithmic wavelength
            flx_ng = idex
            allocate(flx_wl(flx_ng+1))
            help = wlmax/dble(wlmin)
            forall(i=1:flx_ng+1) flx_wl(i) =
     &        wlmin*help**((i-1d0)/flx_ng)
         else
            stop 'generate_fluxgrid: invalid nflx(1)'
         endif
c
c-- polar projection
      elseif(iflx==2) then
         if(idex>0) then
c-- uniform polar array
            flx_nmu = idex
            allocate(flx_mu(flx_nmu+1))
            help = 2d0/flx_nmu
            forall(i=1:flx_nmu+1) flx_mu(i) = -1d0+(i-1)*help
         else
            stop 'generate_fluxgrid: invalid nflx(2)'
         endif
c
c-- azimuthal angle
      else
         if(idex>0) then
c-- uniform azimuthal array
            flx_nom = idex
            allocate(flx_om(flx_nom+1))
            help = pc_pi2/flx_nom
            forall(i=1:flx_nom+1) flx_om(i) = (i-1)*help
         else
            stop 'generate_fluxgrid: invalid nflx(3)'
         endif
      endif
c
      end subroutine generate_fluxgrid
c
c
      subroutine flux_dealloc
      deallocate(flx_wl)
      deallocate(flx_mu)
      deallocate(flx_om)
      deallocate(flx_luminos)
      deallocate(flx_lumdev)
      deallocate(flx_lumnum)
      end subroutine flux_dealloc
c
c
      end module fluxmod
c vim: fdm=marker
