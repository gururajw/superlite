!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine analytic_opacity

  use inputparmod
  use groupmod
  use gridmod
  use gasmod
  use physconstmod
  use miscmod, only:specint
  implicit none

!#####################################
  !This subroutine computes Planck and Rosseland
  !opacities for several simple test group structures.
  !The power law coefficients are used in any structure
  !selection to either fully or partially determine
  !opacity dependence on temperature and density.
  !
  !calculates grey scattering opacity, gas_sig
!#####################################

  integer :: i, ig
  real*8 :: x1(grp_ng),x2(grp_ng)  !unitless energy group bounds
  real*8 :: capcoef(gas_ncell)

  capcoef =0d0
  gas_cap = 0.
  gas_sig = 0d0

  !Calculating grey scattering opacity
  gas_sig = in_gas_sigcoef*gas_temp**in_gas_sigtpwr* &
       gas_rho**in_gas_sigrpwr

  !Calculating grouped Planck and Rosseland opacities
  if(in_opacanaltype=='none') then
     return
  elseif(in_opacanaltype=='grey') then
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        gas_cap(:,i) = sngl(capcoef(i))
     enddo

  elseif(in_opacanaltype=='mono') then
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
!-- monotonic group dependence
     do i = 1, gas_ncell
        x1 = pc_h*pc_c*grp_wlinv(2:)/pc_kb
        x2 = pc_h*pc_c*grp_wlinv(:grp_ng)/pc_kb
        gas_cap(:,i) = sngl(0.5d0*capcoef(i)*(x1+x2)/(x1*x2)**2)
     enddo
  elseif(in_opacanaltype=='line') then
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr * &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        !
        !set odd group magnitudes (low)
        do ig = 1, grp_ng, 2
           gas_cap(ig,i) = sngl(capcoef(i)*in_ldisp1)
        enddo
        !set even group magnitudes (high)
        do ig = 2, grp_ng, 2
           gas_cap(ig,i) = sngl(capcoef(i)*in_ldisp2)
        enddo
        !
     enddo
  else
    stop 'analytic_opacity: in_opacanaltype invalid'
  endif

end subroutine analytic_opacity
! vim: fdm=marker
