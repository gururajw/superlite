*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module totalsmod
c     ----------------
      implicit none
************************************************************************
* Miscellaneous domain-integrated quantities.
************************************************************************
c-- energy conservation check quantities
      real*8 :: tot_eout = 0d0 !energy escaped
      real*8 :: tot_evelo = 0d0 !total energy change to rad field from fluid
      real*8 :: tot_esurf = 0d0
c
c-- sources (domain integrated, not time integrated)
c-- gas_ sources (need to be collected)
      real*8 :: tot_sthermal = 0d0
c-- flx_ sources (negative)
      real*8 :: tot_sflux = 0d0
c
      save
c
      end module totalsmod
c vim: fdm=marker
