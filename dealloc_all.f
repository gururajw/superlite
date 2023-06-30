*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine dealloc_all
c     ----------------------
      use mpimod
      use ionsmod
      use bbxsmod
      use inputparmod
      use nltemod
      use gridmod
      use groupmod
      use gasmod
      use particlemod
      use fluxmod
      use sourcemod
      use randommod
      use sourcemod
      implicit none
************************************************************************
* deallocate all used at the end of the program. Any allocatable arrays
* remaining allocated after this had to be dealt with earlier.
* This helps to catch memory leaks!
************************************************************************
c-- ionsmod
      call ions_dealloc
      call gas_dealloc
      call grid_dealloc
      call flux_dealloc
      deallocate(prt_particles,prt_isvacant)
      call mpimod_dealloc
      deallocate(grp_wl,grp_wlinv)
      deallocate(src_nvacantall)
      deallocate(rnd_states)
      if(allocated(bb_xs)) deallocate(bb_xs) !only impi==impi0, but only if nobbopac==f
      if(in_nlte) call nlte_dealloc
      call source_dealloc

      end subroutine dealloc_all
c vim: fdm=marker
