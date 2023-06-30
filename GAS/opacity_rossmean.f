*This file is part of SuperLite.  SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2020-2025 Gururaj A Wagle.  All rights reserved.
      subroutine opacity_rossmean
c     -----------------------------
      use inputparmod
      use groupmod
      use physconstmod
      use gasmod
      use miscmod
      implicit none
************************************************************************
* Convert opacity into Rosseland mean opacities
************************************************************************
      integer :: i
      real*8 :: specarr(grp_ng)
c
c-- Rossland opacity
      do i=1,gas_ncell
       call spectintv(1/gas_temp(i),grp_ng,specarr)
       gas_capross(i) = 1/sum((1/(gas_sig(i)+gas_cap(:,i)))*specarr)
      enddo !i
c
      end subroutine opacity_rossmean
c vim: fdm=marker
