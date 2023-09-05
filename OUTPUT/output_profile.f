*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine output_profile
c     ----------------------
      use inputparmod
      use timingmod
      use elemdatamod
      use miscmod, only:lcase
      use gridmod
      use gasmod, only:gas_nelem
      use physconstmod
      implicit none
************************************************************************
* Write profilles to files at given interval - added - gaw
************************************************************************
      integer :: j,l
c
      character(4) :: el_name
      character(16) :: fstat='replace'
      character(32) :: file_name
c
      real*8,allocatable :: arr(:)
      real*8 :: t0,t1
c
      t0 = t_time()
c
c-- write profile
c================
c
      file_name='output.profile'
      allocate(arr(grd_ncell))
      arr(:grd_ncell)=1d0/grd_tempinv
      open(unit=4,file=file_name,status=fstat,recl=10000)
      write(4,'(12(a12))',advance='no') '#     x_left','x_right',
     &        'vx_left','vx_right','mass','rho','vol','avg_temp',
     &        'rad_temp','ye','n_e','n_atom'
      do l=1,gas_nelem
        if(all(grd_massfr(l,:)==0d0)) cycle
        el_name = elem_data(l)%sym
        write(4,'(a12)',advance='no') lcase(trim(el_name))
      enddo !l
      write(4,*)
      do j=1,grd_ncell
        write(4,'(12(1p,e12.4))',advance='no')
     &          grd_xarr(j),grd_xarr(j+1),grd_vxarr(j),grd_vxarr(j+1),
     &          grd_mass(j),grd_rho(j),grd_vol(j),arr(j),grd_radtemp(j),
     &          grd_ye(j),grd_nelec(j)*grd_natom(j)/grd_vol(j),
     &          grd_natom(j)/grd_vol(j)
        do l=1,gas_nelem
          if(all(grd_massfr(l,:)==0d0)) cycle
          write(4,'(1p,e12.4)',advance='no') grd_massfr(l,j)
        enddo !l
        write(4,*)
      enddo !j
      close(4)
      deallocate(arr)
c
c-- timing
      t1 = t_time()
      call timereg(t_output, t1-t0)
c
      end subroutine output_profile
c vim: fdm=marker
