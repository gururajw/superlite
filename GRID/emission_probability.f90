!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine emission_probability

  use miscmod
  use inputparmod
  use timingmod
  use gridmod
  use groupmod
  use physconstmod
  implicit none
!-----------------------
!multigroup volume emission probabilities
!-----------------------

  integer :: l,ig,igp1,iep,nepg
  real*8 :: t0,t1
  real*8 :: help
  real*8 :: specarr(grp_ng)

  t0 = t_time()

!-- grouped volume emission probabilities:

!-- one group
  if(grp_ng==1) then
     if(grd_nepg>1) stop 'emis_prob: grp_ng=1 & grd_nepg>1'
     if(in_nlte) then !NLTE
       grd_emitprob(1,grd_idd1:grd_idd1+grd_ndd-1) = &
            grd_capemitgrey(grd_idd1:grd_idd1+grd_ndd-1)
     else !LTE
       grd_emitprob(1,grd_idd1:grd_idd1+grd_ndd-1) = &
            grd_capgrey(grd_idd1:grd_idd1+grd_ndd-1)
     endif
     return
  endif

!-- multi-group
  do l=grd_idd1,grd_idd1+grd_ndd-1
!-- piecewise integration of planck function
     call specintv(grd_tempinv(l),grp_ng,specarr)
!-- cumulative sum of unnormalized emission probability
     ig = 1
     help = 0d0
     do iep=1,grd_nep-1
        nepg = min(iep*grd_nepg,grp_ng) - ig + 1
        igp1 = ig + nepg - 1
        if(in_nlte) then !NLTE
          help = help + sum(specarr(ig:igp1)* &
                            grd_capemit(ig:igp1,l))
        else !LTE
          help = help + sum(specarr(ig:igp1)*grd_cap(ig:igp1,l))
        endif
        grd_emitprob(iep,l) = help
        ig = igp1 + 1
     enddo !iep
     if(in_nlte) then !NLTE
       grd_emitprob(grd_nep,l)=grd_capemitgrey(l)
     else !LTE
       grd_emitprob(grd_nep,l)=grd_capgrey(l)
     endif
  enddo !l

  t1 = t_time()
  call timereg(t_emitp,t1-t0)

end subroutine emission_probability
! vim: fdm=marker
