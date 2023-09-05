*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module mpimod
c     -------------
      implicit none
      integer :: MPI_COMM_WORLD=0
      integer,parameter :: MPI_MAX_PROCESSOR_NAME=13
      integer,private :: ierr=0
      integer :: impi=0  !mpi rank
      integer,parameter :: impi0=0 !master mpi rank
      logical :: lmpi0   !true on the master mpi rank
      integer :: nmpi=1  !number of mpi tasks
c
      save
c
      contains
c
      subroutine bcast_permanent
      end subroutine bcast_permanent
c
c
      subroutine scatter_inputstruct(ndim,icell1,ncell)
      use inputstrmod
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3)
      integer,intent(out) :: icell1,ncell
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
************************************************************************
      integer :: dmy
c
      icell1 = 1
      ncell = str_nc
c
      dmy = ndim(1) !use the intent(in) variable
      allocate(str_massdd(ncell))
      str_massdd = reshape(str_massdc,[ncell])
      if(str_lvol) then
        allocate(str_voldd(ncell))
        str_voldd = reshape(str_voldc,[ncell])
      endif
      if(str_nabund>0) then
       allocate(str_massfrdd(str_nabund,ncell))
       str_massfrdd = reshape(str_massfrdc,[str_nabund,ncell])
      endif
c
      allocate(str_tempdd(ncell))
      str_tempdd = reshape(str_tempdc,[ncell])
c
      if(str_lradtemp) then
       allocate(str_radtempdd(ncell))
       str_radtempdd = reshape(str_radtempdc,[ncell])
      endif
      if(str_llum) then
       allocate(str_lumdd(ncell))
       str_lumdd = reshape(str_lumdc,[ncell])
      endif
      if(str_lye) then
       allocate(str_yedd(ncell))
       str_yedd = reshape(str_yedc,[ncell])
      endif
      end subroutine scatter_inputstruct
c
      subroutine allgather_gridvolume
c     -------------------------------
      use gridmod
      use gasmod
      grd_vol = reshape(gas_vol,[grd_ncell])
      end subroutine allgather_gridvolume
c
      subroutine allgather_radtemp
c     -------------------------------
      use gridmod
      use gasmod
      grd_radtemp = reshape(gas_radtemp,[grd_ncell])
      grd_lum = reshape(gas_lum,[grd_ncell])
      end subroutine allgather_radtemp
c
      subroutine allgather_jrad
c     -------------------------------
      use gridmod
      use gasmod
      use groupmod
      grd_jrad = reshape(gas_jrad,[grp_ng,grd_ncell])
      end subroutine allgather_jrad
c
c
      subroutine bcast_nonpermanent
      use gridmod
      use gasmod
      use groupmod
      use sourcemod
      use particlemod
************************************************************************
* Broadcast the data that changes with time.
* - stub
************************************************************************
c-- domain decomposition
      grd_tempinv = reshape(1d0/gas_temp,[grd_ncell])
      grd_emit = reshape(gas_emit,[grd_ncell])
      grd_jrad = reshape(gas_jrad,[grp_ng,grd_ncell])
      grd_lum = reshape(gas_lum,[grd_ncell])
c
      grd_cap = reshape(gas_cap,[grp_ng,grd_ncell])
      grd_capemit = reshape(gas_capemit,[grp_ng,grd_ncell])
      grd_emiss = reshape(gas_emiss,[grp_ng,grd_ncell])
      grd_sig = reshape(gas_sig,[grd_ncell])
      grd_capgrey = reshape(gas_capgrey,[grd_ncell])
      grd_capemitgrey = reshape(gas_capemitgrey,[grd_ncell])
      grd_capross = reshape(gas_capross,[grd_ncell])
      grd_nelec = reshape(gas_nelec,[grd_ncell])
      grd_natom = reshape(gas_natom,[grd_ncell])
c
      src_nvacantall(1) = count(prt_isvacant)
      end subroutine bcast_nonpermanent
c
c
      subroutine allgather_leakage
      end subroutine allgather_leakage
c
c
      subroutine reduce_gridtally
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      use gridmod
      use gasmod
      gas_eraddens = reshape(grd_tally,[grd_ncell])
      end subroutine reduce_gridtally
c
c
      subroutine allreduce_jrad !nLTE ! added - gaw
************************************************************************
* Reduce the results from particle_advance for radiation intensity
************************************************************************
      use gridmod
      use gasmod
      use groupmod
      gas_jrad = reshape(grd_jrad,[grp_ng,grd_ncell])
      end subroutine allreduce_jrad
c
c
      subroutine scatter_jrad
      end subroutine scatter_jrad
c
c
      subroutine reduce_srctally
      end subroutine reduce_srctally
c
c
      subroutine reduce_fluxtally
      end subroutine reduce_fluxtally
c
c
      subroutine reduce_totals
      end subroutine reduce_totals
c
      subroutine reduce_gasprop
      use gridmod
      use gasmod
      grd_mass = reshape(gas_mass,[grd_ncell])
      grd_rho = reshape(gas_rho,[grd_ncell])
      grd_ye = reshape(gas_ye,[grd_ncell])
      grd_massfr = reshape(gas_massfr,[gas_nelem,grd_ncell])
      end subroutine reduce_gasprop
c
      subroutine scatter_restart_data
      end subroutine scatter_restart_data
c
c
      subroutine collect_restart_data
      end subroutine collect_restart_data
c
c
      subroutine mpimod_dealloc
      end subroutine mpimod_dealloc
c
c-- MPI intrinsics
c-----------------
      subroutine mpi_init(ierr_)
      implicit none
      integer :: ierr_
      ierr_ = ierr
      end subroutine mpi_init
c
      subroutine mpi_comm_rank(mpi_comm,impi_,ierr_)
      implicit none
      integer :: mpi_comm
      integer :: impi_,ierr_
      ierr_ = ierr
      impi_ = impi
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_comm_rank
c
      subroutine mpi_comm_size(mpi_comm,nmpi_,ierr_)
      implicit none
      integer :: mpi_comm
      integer :: nmpi_,ierr_
      ierr_ = ierr
      nmpi_ = nmpi
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_comm_size
c
      subroutine mpi_get_processor_name(pname,ilen_,ierr_)
      implicit none
      character*(MPI_MAX_PROCESSOR_NAME) :: pname
      integer :: ilen_,ierr_
      pname = 'NOT AVAILABLE'
      ierr_ = ierr
      ilen_ = 1
      end subroutine mpi_get_processor_name
c
      subroutine mpi_barrier(mpi_comm,ierr_)
      implicit none
      integer :: mpi_comm,ierr_
      ierr_ = ierr
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_barrier
c
      subroutine mpi_finalize(ierr_)
      implicit none
      integer :: ierr_
      ierr_ = ierr
      end subroutine mpi_finalize
c
      end module mpimod
c vim: fdm=marker
