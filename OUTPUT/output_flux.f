*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine output_flux
c     --------------------------
      use fluxmod
      use inputparmod, only:in_L_bol
      implicit none
************************************************************************
* Write flux arrays to file.
************************************************************************
      integer :: j,k
      logical,save :: lfirst=.true.
      character(16),save :: fstat='replace'
      integer :: reclen

      if(lfirst) then
       reclen = max(4,flx_ng+1,flx_nmu+1,flx_nom+1)*12
       open(unit=4,file='output.flx_grid',status='replace',recl=reclen)
c-- header: dimension
       write(4,*) "#",flx_ng,flx_nmu,flx_nom
c-- body
       write(4,'(1p,10000e12.4)') flx_wl(:)
       write(4,'(1p,100e12.4)') flx_mu(:)
       write(4,'(1p,100e12.4)') flx_om(:)
       close(4)
      endif
c
c-- flux arrays
c==============
      reclen = flx_ng*12
      open(unit=4,file='output.flx_luminos',status=fstat,
     &  position='append',recl=reclen)
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(1p,10000e12.4)') merge(flx_luminos(:,j,k),0d0,
     &   mask=flx_luminos(:,j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
      enddo
      enddo
      close(4)

      open(unit=4,file='output.flx_lumnum',status=fstat,
     &  position='append',recl=reclen)
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(10000i12)') flx_lumnum(:,j,k)
      enddo
      enddo
      close(4)

      open(unit=4,file='output.flx_lumdev',status=fstat,
     &  position='append',recl=reclen)
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(1p,10000e12.4)') merge(flx_lumdev(:,j,k),0d0,
     &   mask=flx_lumdev(:,j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
      enddo
      enddo
      close(4)

c
c-- after the first iteration open files in append mode
      lfirst = .false.
      fstat = 'old'
c
      end subroutine output_flux
c vim: fdm=marker
