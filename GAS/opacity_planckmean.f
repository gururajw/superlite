*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine opacity_planckmean
c     -----------------------------
      use inputparmod
      use groupmod
      use physconstmod
      use gasmod
      use miscmod
      implicit none
************************************************************************
* Convert opacity into grey opacities
************************************************************************
      integer :: i
      real*8 :: specarr(grp_ng)
c
c-- Planck opacity
      do i=1,gas_ncell
         call specintv(1/gas_temp(i),grp_ng,specarr)
         gas_capgrey(i) = sum(gas_cap(:,i)*specarr)
         if(in_nlte) then
          gas_capemitgrey(i) = sum(gas_capemit(:,i)*specarr)
         endif
      enddo !i
c
c
      end subroutine opacity_planckmean
c vim: fdm=marker
