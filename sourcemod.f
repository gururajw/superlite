*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module sourcemod
c     ----------------
      implicit none
      integer :: src_ns
      integer :: src_nnew
      integer :: src_nsurf
      integer :: src_nsurftot
c
      integer,allocatable :: src_ivacant(:) !array of vacant particle array locations
      integer*8,allocatable :: src_nvacantall(:) !(nmpi) number of vacancies on each of the ranks
c
c-- source sampling for output
      real*8,allocatable :: src_energy(:,:) !(grd_ncell,grp_ng) spectral energy distribution of particles for each cell
      integer,allocatable :: src_number(:,:) !(grd_ncell,grp_ng) number distribution of source particles for each cell
c
      save
c
      contains

      subroutine sourcemod_init(nmpi)
c     ----------------------------------------
      use gridmod
      use groupmod
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* init particle module
************************************************************************
      allocate(src_nvacantall(nmpi))
      src_nvacantall = 0
      allocate(src_energy(grd_ncell,grp_ng))
      allocate(src_number(grd_ncell,grp_ng))
      end subroutine sourcemod_init
c
      subroutine source_dealloc
      deallocate(src_energy)
      deallocate(src_number)
      end subroutine source_dealloc
c
      end module sourcemod
c vim: fdm=marker
