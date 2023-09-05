*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine grid_volume(igeom)
c     ---------------------------------------------
      use gridmod
      use physconstmod
      implicit none
      integer,intent(in) :: igeom
************************************************************************
* calculate volumes for all grid cells
************************************************************************
      integer :: i,j,k,l
      real*8,allocatable :: vol(:,:,:) !too big for the stack
c
      allocate(vol(grd_nx,grd_ny,grd_nz))
c
      select case(igeom)
      case(11)
       forall(i=1:grd_nx,j=1:grd_ny,k=1:grd_nz)
        vol(i,j,k) =
     &    (grd_xarr(i+1)**3 - grd_xarr(i)**3) *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    (grd_zarr(k+1) - grd_zarr(k))/3d0
       endforall
      case default
       stop 'grid_volume: invalid igeom'
      endselect !igeom
c
c-- map to compressed grid
      grd_vol = 0d0
      do k=1,grd_nz
      do j=1,grd_ny
      do i=1,grd_nx
       l = grd_icell(i,j,k)
       grd_vol(l) = grd_vol(l) + vol(i,j,k) !multiple void cells are linked to the dummy cell
      enddo
      enddo
      enddo
c
      deallocate(vol)
c
      end subroutine grid_volume
c vim: fdm=marker
