*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine output_source
c     --------------------------
      use sourcemod
      use gridmod
      use groupmod
      implicit none
************************************************************************
* Write source arrays to a file.
************************************************************************
      integer :: i
      character(16),save :: fstat='replace'
      integer :: reclen
c
c-- source arrays
c==============
      reclen = grp_ng*14
c
      open(unit=4,file='output.src_number',status=fstat,
     &  position='append',recl=reclen)
      do i=1,grd_ncell
       write(4,'(1p,10000i12)') src_number(i,:)
      enddo
      close(4)
c
      open(unit=4,file='output.src_luminos',status=fstat,
     &  position='append',recl=reclen)
      do i=1,grd_ncell
       write(4,'(*(es14.4e3))') src_energy(i,:)
      enddo
      close(4)
c
c-- after the first iteration open files in append mode
      fstat = 'old'
c
      end subroutine output_source
c vim: fdm=marker
