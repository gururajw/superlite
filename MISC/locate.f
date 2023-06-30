*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      pure integer function locate(xarray,n,x)
c     -------------------------------
      implicit none
      integer,intent(in) :: n
      real*8,intent(in) :: xarray(n),x
****************************************************************
* search in an ordered table.
* from numerical recipes
****************************************************************
      integer jl,jm,ju
c
      jl = 0
      ju = n + 1
      do
       jm = (ju + jl)/2
       if((xarray(n)>=xarray(1)) .eqv. (x>=xarray(jm))) then
        jl=jm
       else
        ju=jm
       endif
       if(ju-jl .le. 1) exit
      enddo
c--
      locate = jl
      end function locate
c vim: fdm=marker
