!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
subroutine particle_advance

!$ use omp_lib
  use randommod
  use transportmod
  use miscmod
  use particlemod
  use totalsmod
  use groupmod
  use gridmod
  use physconstmod
  use inputparmod
  use timingmod
  use countersmod
  use fluxmod
  use sourcemod
  implicit none
!
!##################################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or escape as a flux particle.
  !DDMC and IMC particle events are being handled in separate subroutines.
!##################################################
  integer :: i
  integer :: npckt, nflux
  integer :: nstepddmc, nstepimc, nmethodswap
  real*8 :: r1, help, tau
  integer :: ipart
  integer, pointer :: ig, ic
  integer, pointer :: ix, iy, iz
  real*8, pointer :: x,y,z, mu, e, e0, wl, om
  real*8 :: vx,vy,vz,help1,help2,vhelp1,vhelp2
  real*8 :: t0,t1
  real*8 :: labfact, mu1, mu2
!
  integer :: icold, igold, ierr
  real*8 :: eraddens, eamp
  real*8 :: jrad
!-- specint cache
  integer :: nstepmax, ndist(-3:7)
!
  type(packet),target :: ptcl
  type(packet2),target :: ptcl2
  type(grp_t_cache),target :: cache
  real*8,target :: specarr(grp_ng)
  integer*2,target :: glumps(grp_ng)
  logical*2,target :: llumps(grp_ng)
!
  type(rnd_t) :: rndstate
  integer,save :: iomp=0
!
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
!
!-- init timers
  t_pckt_stat = (/1d30, 0d0, 0d0/) !min,mean,max
!
!-- energy tallies
  grd_tally = 0d0
  grd_jrad = 0.d0
  grd_eamp = 0d0

!-- optional grid tally
  if(in_io_dogrdtally) then
     grd_methodswap = 0
  endif

!-- Propagate all particles that are not considered vacant
  npckt = 0
  nflux = 0
  nstepddmc = 0
  nstepimc = 0
  nmethodswap = 0

  nstepmax = 0
  ndist = 0

!-- keep track of the number of particles in the particle array
  call counterreg(ct_nnonvacant,count(.not.prt_isvacant))

!$omp parallel &
!$omp private(ptcl,ptcl2,cache, &
!$omp    specarr,glumps,llumps, &
!$omp    x,y,z,mu,om,wl,e,e0,ix,iy,iz,ic,ig,icold,igold, &
!$omp    r1,vx,vy,vz,help1,help2,vhelp1,vhelp2, &
!$omp    i,t0,t1,help,tau,mu1,mu2,labfact, &
!$omp    rndstate,eraddens,eamp,ierr, iomp) &
!$omp reduction(+:grd_tally,grd_jrad, &
!$omp    tot_evelo,tot_eout,tot_sflux, &
!$omp    npckt,nflux, &
!$omp    nstepddmc,nstepimc,nmethodswap,ndist) &
!$omp reduction(max:nstepmax)

!-- thread id
!$ iomp = omp_get_thread_num()

!-- each thread uses its own rnd stream
  rndstate = rnd_states(iomp+1)

!$!-- thread timing
!$ if(.false.) then
      t0 = t_time() !serial version
!$ else
!$    t0 = omp_get_wtime()
!$ endif
!
!-- assigning pointers to corresponding particle properties
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  wl => ptcl%wl
  e => ptcl%e
  e0 => ptcl%e0
!-- secondary particle properties
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
!-- alloc
  cache%specarr => specarr
  cache%glumps => glumps
  cache%llumps => llumps
  cache%ic = 0

!$omp do schedule(static,1) !round-robin
  do ipart=1,prt_npartmax
     if(prt_isvacant(ipart)) cycle
!
!-- active particle
     ptcl = prt_particles(ipart) !copy properties out of array
     ptcl2%ipart = ipart
     npckt = npckt+1

!-- cell position
     ix = binsrch(x,grd_xarr,grd_nx+1,.false.)
     iy = binsrch(y,grd_yarr,grd_ny+1,.false.)
     iz = binsrch(z,grd_zarr,grd_nz+1,.false.)
!-- cell pointer
     ic = grd_icell(ix,iy,iz)

!-- group index
     ig = binsrch(wl,grp_wl,grp_ng+1,.false.) !co-moving frame
!
!-- determine particle type
     select case(grd_igeom)
     case(11)
        help = dx(ix)
     case default
        stop 'particle_advance: invalid igeom'
     endselect
!-- IMC or DDMC
     tau = (grd_sig(ic)+grd_cap(ig,ic))*help
     if(in_puretran .or. tau<trn_tauddmc) then
        ptcl2%itype = 1 !IMC
     else
        ptcl2%itype = 2 !DDMC
     endif
!-- interpolate velocity
!-- x
     vhelp1 = grd_vxarr(ix)
     vhelp2 = grd_vxarr(ix+1)
     help1 = grd_xarr(ix)
     help2 = grd_xarr(ix+1)
     vx = (vhelp1*(help2-x)+vhelp2*(x-help1))/(help2-help1)
!-- y
     vy = 0 ! 1-d sph
!-- z
     vz = 0 ! 1-d sph

!-- transform IMC particle into lab frame
     if(ptcl2%itype==1) then
        select case(grd_igeom)
        case(11) ! 1-d sph
           labfact = 1d0-vx*mu/pc_c
        endselect
!-- transform into lab frame
        wl = wl*labfact
        e = e/labfact
        e0 = e0/labfact
     endif

!-----------------------------------------------------------------------
!-- Advancing particle until interaction, or escape from domain
!Calling either diffusion or transport depending on particle type (ptcl2%itype)
     ptcl2%istep = 0
     ptcl2%idist = 0
     ptcl2%stat = 'live'

     do while (ptcl2%stat=='live')
        ptcl2%istep = ptcl2%istep + 1
        icold = ic
        igold = ig
        if(ptcl2%itype==1 .or. in_puretran) then
           nstepimc = nstepimc + 1
           call transport(ptcl,ptcl2,vx,vy,vz,rndstate, &
                eraddens,eamp,jrad,tot_evelo,ierr)
           if(ptcl2%itype/=1) then
              nmethodswap = nmethodswap + 1
              if(in_io_dogrdtally) grd_methodswap(icold) = grd_methodswap(icold) + 1
           endif
!-- tally eamp
           if(.not.trn_noampfact) grd_eamp(icold) = grd_eamp(icold) + eamp
        else
           nstepddmc = nstepddmc + 1
           call diffusion(ptcl,ptcl2,vx,vy,vz,cache,rndstate, &
                eraddens,jrad,tot_evelo,ierr)
           if(ptcl2%itype==1) then
              nmethodswap = nmethodswap + 1
              if(in_io_dogrdtally) grd_methodswap(icold) = grd_methodswap(icold) + 1
           endif
        endif
        i = max(-3,ptcl2%idist)
        ndist(i) = ndist(i) + 1
!-- tally rest
        grd_tally(icold) = grd_tally(icold) + eraddens
        grd_jrad(igold,icold) = grd_jrad(igold,icold) + jrad
!
!-- verify position
        if(ptcl2%itype==1 .and. ptcl2%stat=='live') then
           if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix) .or. x/=x) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv: x not in cell', &
                 ix,x,grd_xarr(ix),grd_xarr(ix+1)
           endif
           if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy) .or. y/=y) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv: y not in cell', &
                 iy,y,grd_yarr(iy),grd_yarr(iy+1)
           endif
           if(z>grd_zarr(iz+1) .or. z<grd_zarr(iz) .or. z/=z) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv: z not in cell', &
                 iz,z,grd_zarr(iz),grd_zarr(iz+1)
           endif
        endif

!-- check exit status
        if(ierr/=0) then !.or. ptcl2%istep>1000) then  !istep checker may cause issues in high-res simulations
           write(0,*) 'pa: ierr,ipart,istep,idist:',ierr,ptcl2%ipart,ptcl2%istep,ptcl2%idist
           write(0,*) 'dist:',ptcl2%dist
           write(0,*) 'taus,tauc:',grd_sig(ic)*help,grd_cap(ig,ic)*help
           write(0,*) 'ix,iy,iz,ic,ig:',ptcl2%ix,ptcl2%iy,ptcl2%iz,ptcl2%ic,ptcl2%ig
           write(0,*) 'x,y,z:',ptcl%x,ptcl%y,ptcl%z
           write(0,*) 'vx,vy,vz:',vx,vy,vz
           write(0,*) 'mu,om:',ptcl%mu,ptcl%om
           write(0,*)
           if(ierr>0) then
              if(trn_errorfatal) stop 'particle_advance: fatal transport error'
              ptcl2%stat = 'dead'
              exit
           endif
        endif
     enddo

!-- max step counter
     nstepmax = max(nstepmax,ptcl2%istep)


!-- mark particle slot occupied or vacant
     select case(ptcl2%stat)
     case('dead')
        prt_isvacant(ipart) = .true.
!
!-- outbound luminosity tally
     case('flux')
!-- register fresh flux particle (only once)
        nflux = nflux + 1
        tot_sflux = tot_sflux-e
        tot_eout = tot_eout+e
        ptcl%x = huge(help)  !-- mark particle as flux particle (stat is not saved in particle array)
!
!-- save particle properties
        prt_particles(ipart) = ptcl
!
    case default
       stop 'particle_advance: invalid particle status'
    end select

  enddo !ipart
!$omp end do nowait

! write(0,*) iomp,nstepmax, ndist

!-- save state
  rnd_states(iomp+1) = rndstate

!-- thread timing
!$ if(.false.) then
      t1 = t_time() !serial version
!$ else
!$    t1 = omp_get_wtime()
!$ endif

!$omp critical
  t1 = t1 - t0
  t_pckt_stat = (/min(t_pckt_stat(1),t1), t_pckt_stat(2)+t1/in_nomp, &
     max(t_pckt_stat(3),t1)/) !min,mean,max
!$omp end critical
!$omp end parallel

!-- print distance counters
! if(lmpi0) write(6,'(11(i6,"k"))',advance='no') ndist(-3:7)/1000

  call counterreg(ct_npactive, npckt)
  call counterreg(ct_npstepimc, nstepimc)
  call counterreg(ct_npstepddmc, nstepddmc)
  call counterreg(ct_npstepmax, nstepmax)
  call counterreg(ct_npmethswap, nmethodswap)
  call counterreg(ct_npflux, nflux)


end subroutine particle_advance
! vim: fdm=marker
