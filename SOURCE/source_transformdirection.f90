!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine source_transformdirection

  use particlemod
  use sourcemod
  use inputparmod
  use physconstmod
  use transportmod
  use gridmod
  use miscmod

  implicit none

  integer :: ipart,ivac
  integer :: ix,iy,iz
  real*8 :: vx,vy,vz,help1,help2,vhelp1,vhelp2

  type(packet),pointer :: p
  real*8,pointer :: x,y,z

!-- transform particle direction into lab frame
  do ipart=1,src_nnew
     ivac = src_ivacant(ipart)
     p => prt_particles(ivac)
     x => p%x
     y => p%y
     z => p%z

     ix = binsrch(x,grd_xarr,grd_nx+1,.false.)
     iy = binsrch(y,grd_yarr,grd_ny+1,.false.)
     iz = binsrch(z,grd_zarr,grd_nz+1,.false.)
!-- interpolate velocity coordinate
!-- x
     vhelp1 = grd_vxarr(ix)
     vhelp2 = grd_vxarr(ix+1)
     help1 = grd_xarr(ix)
     help2 = grd_xarr(ix+1)
     vx = (vhelp1*(help2-x)+vhelp2*(x-help1))/(help2-help1)
!-- y
     if(grd_igeom==11) vy=0 !1-D spherical
!-- z
     if(grd_igeom==11.or.grd_igeom==2) vz=0 !1-D spherical

     call direction2lab(vx,vy,vz,p%mu,p%om)
  enddo

end subroutine source_transformdirection
! vim: fdm=marker
