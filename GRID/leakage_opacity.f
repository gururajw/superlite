*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine leakage_opacity
c     ----------------------------------------
      use gridmod
      use timingmod
      implicit none
************************************************************************
* wrapper
************************************************************************
      real*8 :: t0,t1
c
      t0 = t_time()
c
      select case(grd_igeom)
      case(11)
         call leakage_opacity11
      case default
         stop 'leakage_opacity: invalid igeom'
      endselect
c
      t1 = t_time()
      call timereg(t_opacleak,t1-t0)
c
      end subroutine leakage_opacity
c vim: fdm=marker
