*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine gas_update(it)
c     -------------------------
      use gridmod
      use physconstmod
      use miscmod
      use ionsmod
      use totalsmod
      use gasmod
      use inputparmod
      use inputstrmod, only:str_lvol
      use timingmod
      use mpimod
      implicit none
      integer,intent(in) :: it
************************************************************************
* The work done here is:
* - temperature update (experimental iterative approach)
* - LTE/nebular EOS: ionization balance and electron density
* - LTE/NLTE opacities
************************************************************************
      logical,save :: lfirst=.true.
      logical :: do_output
      integer :: i,l
c-- timing
      real*8 :: t0,t1
c
c-- begin
      t0 = t_time()
c
      if(lfirst) then
c-- update gas volume
c========================================
        if(.not.str_lvol) then
          i = 0
          do l=gas_icell1,gas_icell1+gas_ncell-1
           i = i+1
           gas_vol(i) = grd_vol(l)
          enddo !l
        endif
c
c-- velocity divergence
c========================================
        if(in_nlte) then
          i = 0
          do l=gas_icell1,gas_icell1+gas_ncell-1
           i = i+1
           gas_divv(i) = grd_divv(l)
          enddo !l
        endif

c
c-- update gas density
c===============================================
        gas_rho = gas_mass/gas_vol
c
      endif !lfirst
c
c-- update gas temperature
c===============================================
      i = 0
      do l=gas_icell1,gas_icell1+gas_ncell-1
        i = i+1
        gas_radtemp(i) = grd_radtemp(l)
      enddo !l
c-- temperature
      gas_ur = pc_acoef*gas_radtemp**4
c
c-- sanity check temperatures
      if(any(gas_temp/=gas_temp)) stop 'gas_temp NaN'
      if(any(gas_temp<=0d0)) stop 'gas_temp<=0'
      if(any(gas_radtemp/=gas_radtemp)) stop 'gas_radtemp NaN'
      if(any(gas_radtemp<=0d0)) stop 'gas_radtemp<=0'
c
c
c
c-- solve LTE/nebular EOS
c===============================
      do_output = in_io_pdensdump=='on' .and. it==1
      if(.not.in_noeos) call eos_update(do_output)
c
c
c
c-- calculate opacities
c======================
c-- simple analytical group/grey opacities: Planck and Rosseland
      if(in_opacanaltype/='none') then
        call analytic_opacity
      else
c-- calculate physical opacities
        if(in_nlte) then !NLTE
          call physical_opacity_nlte(it)
        else
          call physical_opacity
        endif

      endif
      call opacity_planckmean
      call opacity_rossmean
c
c
      lfirst = .false.
c
c
      t1 = t_time()
      call timereg(t_gasupd,t1-t0)
c
      end subroutine gas_update
c vim: fdm=marker
