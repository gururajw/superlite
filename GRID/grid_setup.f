*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine grid_setup
c     ---------------------
      use gridmod
      use inputparmod
      use inputstrmod
      use physconstmod
      implicit none
************************************************************************
* Setup the grid on the computational domain
************************************************************************
      integer :: i,j,k,l,idcell
c
c-- agnostic grid setup
      grd_xarr = str_xleft
      grd_yarr = str_yleft
      grd_zarr = str_zleft
      ! Velocity arrays
      grd_vxarr = str_vxleft
c
c-- cell pointers
c-- initialize void cells
      grd_icell = 0
c-- pointers into compressed grid
      l = 1
      idcell = 0
      loop_k: do k=1,grd_nz
       do j=1,grd_ny
       do i=1,grd_nx
        idcell = idcell + 1
        if(idcell == str_idcell(l)) then
         grd_icell(i,j,k) = l
         l = l + 1
        endif
        if(l>grd_ncell) exit loop_k
       enddo
       enddo
      enddo loop_k
      if(l/=grd_ncell+1) stop 'grid_setup: l/=grd_ncell+1'
c
c
c-- maximum grid radius and cell radii, velocity divergence
      select case(grd_igeom)
      case(1,11)
       grd_rcell = (grd_xarr(2:grd_nx+1)+grd_xarr(1:grd_nx))/2
       !-- \vec{\nabla} \cdot \vec{v} = 2*vx_center/x_center + dvx/dx
       grd_divv = 2*(grd_vxarr(2:grd_nx+1) + grd_vxarr(1:grd_nx))/
     &            (grd_xarr(2:grd_nx+1) + grd_xarr(1:grd_nx)) +
     &            (grd_vxarr(2:grd_nx+1) - grd_vxarr(1:grd_nx))/
     &            (grd_xarr(2:grd_nx+1) - grd_xarr(1:grd_nx))
      endselect
c
c-- sanity check
      if(minval(grd_vxarr)<0d0) stop 'grid_setup: grd_vxarr < 0'
      if(maxval(abs(grd_vxarr))>pc_c)
     &    stop 'grid_setup: grd_vxarr > pc_c'
c
c-- zero amplification-factor energy to begin with
      grd_eamp = 0d0
c
      end subroutine grid_setup
c vim: fdm=marker
