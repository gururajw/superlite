*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine output_grid
c     ----------------------
      use inputparmod
      use timingmod
      use gridmod
      use totalsmod
      implicit none
************************************************************************
* Write grid values to file
************************************************************************
      integer :: i,j,k
      integer,save :: nrow=0 !number of rows for printing grid variables
      integer,save :: ncpr=0 !number of cells per rank
      integer :: reclen
      logical,save :: lfirst=.true.
      character(16),save :: fstat='replace'
c
      integer,allocatable :: iarr(:)
      real*8,allocatable :: arr(:)
      real*8 :: t0,t1
c
      t0 = t_time()
c
      if(lfirst) then
c
       ncpr = grd_nx
       nrow = ceiling(float(grd_ncell)/ncpr)
c
c-- write once
c=============
       reclen = max(4,grd_nx+1,grd_ny+1,grd_nz+1,grd_ny*grd_nz)*12
       open(unit=4,file='output.grd_grid',status='replace',recl=reclen)
c-- header: dimension
       write(4,*) "#",grd_igeom
       write(4,*) "#",grd_nx,grd_ny,grd_nz
       write(4,*) "#",nrow*ncpr,nrow,ncpr
c-- body
       write(4,'(1p,10000e12.4)') grd_xarr(:)
       write(4,'(1p,10000e12.4)') grd_yarr(:)
       write(4,'(1p,10000e12.4)') grd_zarr(:)
c-- cell indices
       do k=1,grd_nz
       do j=1,grd_ny
          write(4,'(10000i12)') grd_icell(:,j,k)
       enddo
       enddo
       close(4)
      endif

c
c-- scalars
c==========
      open(unit=4,file='output.tot_energy',status=fstat,
     &  position='append',recl=5*12)
      if(lfirst) then
       write(4,'("#",i5)') 4 !number of columns
       write(4,'("#",16a12)')
     &   'eout','evelo','sflux','sthermal'
      endif
      write(4,'(1x,1p,16e12.4)')
     &  tot_eout,tot_evelo,tot_sflux,tot_sthermal
      close(4)

c
c-- grid arrays
c==============
      if(.not.in_io_nogriddump) then
       reclen = ncpr*12
c-- alloc
       allocate(iarr(ncpr*nrow),arr(ncpr*nrow))
       iarr(grd_ncell+1:) = 0
       arr(grd_ncell+1:) = 0d0
c
c-- grdinput values
       iarr(:grd_ncell) = grd_nvol
       open(unit=4,file='output.grd_nvol',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = 1d0/grd_tempinv
       open(unit=4,file='output.grd_temp',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = grd_radtemp
       open(unit=4,file='output.grd_radtemp',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       if(trim(in_io_opacdump)=='off') then
      ! c-- donothing
       else
         arr(:grd_ncell) = grd_capgrey
         open(unit=4,file='output.grd_capgrey',status=fstat,
     &      position='append',recl=reclen)
         do i=1,nrow
          write(4,'(1p,10000es12.4e3)') arr((i-1)*ncpr+1:i*ncpr)
         enddo
         close(4)
c
         arr(:grd_ncell) = grd_capemitgrey
         open(unit=4,file='output.grd_capemitgrey',status=fstat,
     &      position='append',recl=reclen)
         do i=1,nrow
          write(4,'(1p,10000es12.4e3)') arr((i-1)*ncpr+1:i*ncpr)
         enddo
         close(4)
c
         arr(:grd_ncell) = grd_capross
         open(unit=4,file='output.grd_capross',status=fstat,
     &      position='append',recl=reclen)
         do i=1,nrow
          write(4,'(1p,10000es12.4e3)') arr((i-1)*ncpr+1:i*ncpr) ! modified - gaw
         enddo
         close(4)
c
         arr(:grd_ncell) = grd_sig
         open(unit=4,file='output.grd_sig',status=fstat,
     &      position='append',recl=reclen)
         do i=1,nrow
          write(4,'(1p,10000e12.4e3)') arr((i-1)*ncpr+1:i*ncpr)
         enddo
         close(4)
       endif !opacdump
c
       arr(:grd_ncell) = grd_tally/grd_vol
       open(unit=4,file='output.grd_eraddens',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
c-- grdtally values
       if(in_io_dogrdtally) then
        iarr(:grd_ncell) = grd_methodswap
        open(unit=4,file='output.grd_methodswap',status=fstat,
     &   position='append',recl=reclen)
        do i=1,nrow
         write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
c
       endif !in_io_dogrdtally
c
       deallocate(iarr,arr)
      endif

c
c-- after the first iteration open files in append mode
      lfirst = .false.
      fstat = 'old'
c
c-- timing
      t1 = t_time()
      call timereg(t_output, t1-t0)
c
      end subroutine output_grid
c vim: fdm=marker
