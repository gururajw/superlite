!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
program superlite

  use randommod
  use sourcemod
  use mpimod
  use transportmod
  use inputparmod
  use groupmod
  use gridmod
  use gasmod
  use particlemod
  use physconstmod
  use miscmod
  use totalsmod

  use inputstrmod
  use fluxmod

  use ionsmod, only:ions_read_data,ions_alloc_grndlev
  use nltemod
  use bfxsmod, only:bfxs_read_data
  use ffxsmod, only:ffxs_read_data
  use timingmod
  use countersmod

  implicit none
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer :: ierr, it,i
  integer :: icell1, ncell !number of cells per rank (gas_ncell)
  real*8 :: t0, t1 !timing
  character(25) :: msg,fname
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
!
!-- initialize timing module
  call timingmod_init
  call countersmod_init
!
!--
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
  if(lmpi0) then
     t0 = t_time()
!-- startup message
     call banner
!-- read runtime parameters
     call read_inputpars
!-- redirect stdout to file if selected
     call open_logfiles
!-- parse and verify runtime parameters
     call parse_inputpars(nmpi)
!
!-- read input structure
     call read_inputstr(in_grd_igeom,in_ndim,nmpi)
!-- compressed domain, serialize non-void cells
     call inputstr_compress

!-- READ DATA
!-- read ion and level data
     call ions_read_data(gas_nelem)  !ion and level data
!-- read bbxs data
     if(.not.in_nobbopac) then!bound-bound cross section data
       if(in_nlte) then !NLTE
         call read_nlte_data(gas_nelem,in_nlte_nelem)
       else !lte & nebular approximation
         call read_bbxs_data(gas_nelem)
       endif
     endif
!-- read bfxs data
     if(.not.in_nobfopac) call bfxs_read_data !bound-free cross section data
!-- read ffxs data
     if(.not.in_noffopac) call ffxs_read_data !free-free cross section data
!
     t1 = t_time()
     t_setup = t1-t0
  endif !impi
!-- broadcast init info from impi0 rank to all others
  call bcast_permanent !MPI
  call provide_inputpars(nmpi)
!-- domain-decompose input structure
  call scatter_inputstruct(in_ndim,icell1,ncell) !MPI

!--
!-- setup remaining modules
!==========================
!-- wlgrid (before grid setup)
  call groupmod_init
  call fluxgrid_setup

!-- setup spatial grid
  call gridmod_init(lmpi0,grp_ng,str_nc,str_lvoid,icell1,ncell)
  call grid_setup
!-- setup gas
  call gasmod_init(lmpi0,icell1,ncell,grp_ng)
  call gas_setup
!-- inputstr no longer needed
  call inputstr_dealloc

!-- create procedure pointers for the selected geometry
  call transportmod_init(grd_igeom)
!-- source particles
  call sourcemod_init(nmpi)

!-- allocate arrays of sizes retreived in bcast_permanent
  call ions_alloc_grndlev(gas_nelem,gas_ncell)  !ground state occupation numbers
  if(in_nlte) call nlte_alloc_nlev(nlte_nelem,gas_ncell) !NLTE
  call particle_alloc(lmpi0)

!-- initialize random number generator, use different seeds for each rank
  call rnd_init(in_nomp,impi)
!-- use extra stream for non-threaded code
  rnd_state = rnd_states(1)

!-- volume
  if(str_lvol) then ! gather grd_vol from gas_vol if vol in input.str
    call allgather_gridvolume
  else ! calculate grd_vol using grd_xarr
    call grid_volume(grd_igeom)
  endif

!-- gather initial radiation temperature
  call allgather_radtemp


!-- memory usage
  if(lmpi0) then
     msg = 'post setup:'
     write(6,*) 'memusg: ',msg,memusg()
  endif

!-- Monte Carlo loop
!=================
  if(lmpi0) then
     write(6,*)
     write(6,*) "MC calculations:"
     write(6,*) "===================="
  endif
!
  do it=1,in_niter ! iterations
     t_timelin(1) = t_time() !timeline

!-- write timestep
     if(lmpi0) then !NLTE
      write(6,*)
      write(6,'(1x,a4,1x,i2)') 'it: ',it
     endif
!-- update all non-permanent variables
     call gas_update(it)

     call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
     t_timelin(2) = t_time() !timeline

!-- gather from gas workers and broadcast to world ranks
     call bcast_nonpermanent !MPI

!-- source energy
     call sourceenergy

!-- update photospheric radius for NLTE calculations (experimental iterative method)
     if(in_nlte) call get_rphot
     if(lmpi0.and.in_nlte) then
       write(6,'(1x,a8,1x,es8.2)') "R_phot: ", R_phot
       write(6,*)
     endif
!
     if(lmpi0) write(6,*) "--> Calculating leakage opacity and Emission probability"
!-- source energy
     call leakage_opacity       !IMC-DDMC albedo coefficients and DDMC leakage opacities
     call emission_probability  !emission probabilities for ep-group in each cell
     call allgather_leakage !MPI

     if(lmpi0) write(6,*) "--> Assign source particles"
     call sourcenumbers         !number of source prt_particles

     if(src_nnew>0) then
        allocate(src_ivacant(src_nnew))
        call vacancies          !Storing vacant "prt_particles" indexes in ordered array "src_ivacant"
        call boundary_source    !properties of prt_particles on domain boundary
        call interior_source(it) !properties of prt_particles emitted in domain interior
        call source_transformdirection
        deallocate(src_ivacant)
     endif
!-- tally source particles for output
     call sourcetally
!-- write output for source particle distribution
     call reduce_srctally
     if(lmpi0) call output_source
!
!
!-- advance particles
     if(lmpi0) write(6,*) "--> Particle MC tranport"
     t_timelin(3) = t_time() !timeline
     call particle_advance
     call reduce_gridtally !MPI  !collect transport results from all workers
     call allreduce_jrad !MPI !NLTE !collect radiation intensity from all workers
     call scatter_jrad !MPI

!-- estimated radiation energy
     if(in_nlte) call nlte_rad_est

!-- print packet advance load-balancing info
     if(lmpi0) then
        call timereg(t_pcktmin,t_pckt_stat(1))
        call timereg(t_pcktmea,t_pckt_stat(2))
        call timereg(t_pcktmax,t_pckt_stat(3))
     endif
     t_timelin(4) = t_time() !timeline

!-- tally flux packets
     call fluxtally
!-- output
     if(lmpi0) write(6,*) "--> Tally flux packets"
     call reduce_fluxtally !MPI  !collect flux results from all workers
     if(lmpi0) call output_flux
!
     if(lmpi0) write(6,*) "--> Print output files"

     call reduce_gasprop
     call reduce_totals

!-- output
     if(lmpi0) then

!-- write output
        call output_grid
        call output_grp
        if(it>1.and.in_io_opacdump=='one') in_io_opacdump = 'off'
        if(it==in_niter.and.in_io_profdump) call output_profile

!-- write stdout
        write(6,*)
        write(6,*) "MC Results:"
        write(6,*) "===================="
        write(6,'(2(1x,a7))') 'nsrc','nflux'
        if(ct_nnonvacant(2)<1000000) then
           write(6,'(2(1x,i7))') ct_npcreate(2),ct_npflux(2)
        else
           write(6,'(2(i7,"k"))') ct_npcreate(2)/1000,ct_npflux(2)/1000
        endif
        write(6,*)
     endif !impi
!-- write timestep timing to file
     call timing_cycle(impi)
     call counters_cycle(impi)
     t_timelin(5) = t_time() !timeline
     t_timeline = t_timeline + (t_timelin(2:) - t_timelin(:4))

  enddo ! iterations
!
!
!--
!-- FINISH UP:
!=============
  call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output
  if(lmpi0) then
!
!-- print memory usage
     msg = 'post loop:'
     write(6,*)
     write(6,*) 'memusg: ',msg,memusg()

!-- print cpu timing usage
     t1 = t_time()
     t_all = t1 - t0
     call print_counters
     call print_timing
     write(6,*)
     write(6,*) 'SuperLite finished'
     if(in_io_grabstdout) write(0,'(a,f8.2,"s")')'SuperLite finished',t_all!repeat to stderr
  endif
!-- Clean up memory. (This helps to locate memory leaks)
  call dealloc_all
  call mpi_finalize(ierr) !MPI

end program superlite
! vim: fdm=marker
