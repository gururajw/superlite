*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
*See LANL_COPYING and LANL_README for details of LANL copyright assertion.
      module inputparmod
c     ------------------
      implicit none
************************************************************************
* input parameters
* New parameters need to be added in this routine in four places:
*  1) variable declaration
*  2) namelist
*  3) pointer array
*  4) sanity conditions
************************************************************************
c-- model id and comments
      character(40) :: in_name = "spn"  !simulation name/title, for post-processing identification
      character(80) :: in_comment = ""  !comment if any
c-- parallelization
      integer :: in_nomp = 1  !number of openmp threads
c
c-- experimental iterative approach
      integer :: in_niter = 1 !number of iterations
c
c-- grid geometry and dimensions
      integer :: in_grd_igeom = 11 !only 11=1Dsph is accepted for now
      integer :: in_ndim(3) = [1, 1, 1]  ![nx,1,1] nx: number of x-direction cells
c
c
c
c-- group structure
      integer :: in_grp_ng = -1      !number of groups: 0 read wl-grid from file
      integer :: in_grp_ngs = 1      !number of subgroups per opacity group
      real*8 :: in_grp_wlmin =   100d-8  !lower wavelength boundary [cm]
      real*8 :: in_grp_wlmax = 32000d-8  !upper wavelength boundary [cm]
c
c
c-- outbound flux group and direction bins
      integer :: in_flx_ndim(3) = [0, 1, 1]
      real*8 :: in_flx_wlmin =  1000d-8  !lower wavelength flux boundary [cm]
      real*8 :: in_flx_wlmax = 32000d-8  !upper wavelength flux boundary [cm]
c
c-- interior source
      integer :: in_src_ns = 0     !number of source particles generated per time step (total over all ranks)
      integer :: in_src_n2s = -1   !2^n source particles generated per time step (total over all ranks)
      real*8 :: in_srcepwr = 1d0  !source particle number-energy slope, 1 is linear, equal number of packets per erg.
      real*8 :: in_L_bol = 1d0 ![erg/s] !bolometric luminosity
      real*8 :: in_L_inner = 0d0 ![erg/s] !experimental inner boundary source
c-- boundary source
      character(4) :: in_surfsrcmu = 'isot' !isot|beam: surface source direction distribution
      real*8 :: in_srcmax = 0d0 !peak source strength
c
c-- transport
      logical :: in_trn_errorfatal = .true. !stop on transport error, disable for production runs
      real*8 :: in_trn_tauddmc = 5d0  !number of mean free paths per cell required for DDMC
      real*8 :: in_taulump = 10d0 !number of of mean free paths needed to lump DDMC groups
      logical :: in_puretran = .false. !use IMC only instead of IMC+DDMC hybrid
      logical :: in_trn_noamp = .true.  !disable amplification factor
      logical :: in_thm_aniso = .false. !experimental anisotropic scattering
c
c-- debugging
      logical :: in_noeos = .false.     !don't use the EOS
c
c-- physical opacities
      real*8 :: in_opacmixrossel = 0d0 !mix rosseland with planck average, 1=pure rosseland
c-- physical opacities debuging
      logical :: in_nobbopac = .false.  !turn off bound-bound opacity
      logical :: in_nobfopac = .false.  !turn off bound-bound opacity
      logical :: in_noffopac = .false.  !turn off bound-bound opacity
      logical :: in_nothmson = .false.  !turn off thomson scattering
c
c-- experimental nebular approximation
      logical :: in_nebular = .false. ! turn on nebular approximation
c-- photospheric radius for dilution factor
      real*8 :: in_R_phot = -1 ![cm] photospheric radius where Rosseland mean optical depth = 1
c
c-- NLTE treatment for hydrogen b-b opacity calculations
      logical :: in_nlte = .false. ! turn on NLTE
      integer :: in_nlte_nelem = 1 ! number of elements treated with nLTE ! for now keep =1 for h
c
c-- analytic opacities
      character(4) :: in_opacanaltype = 'none'    !none|grey|mono|line: group opacity structure type
c-- line specific group structure
      real*8 :: in_ldisp1 = 1d0  !loosely speaking, the analytic odd group line strength
      real*8 :: in_ldisp2 = 1d0  !loosely speaking, the analytic even group line strength
c-- scattering terms:
      real*8 :: in_gas_sigcoef = 0d0  !power law absorption opacity coefficient
      real*8 :: in_gas_sigtpwr = 0d0  !power law absorption opacity temperature exponent
      real*8 :: in_gas_sigrpwr = 0d0  !power law absorption opacity density exponent
c-- absorption terms:
      real*8 :: in_gas_capcoef = 0d0  !power law absorption opacity coefficient
      real*8 :: in_gas_captpwr = 0d0  !power law absorption opacity temperature exponent
      real*8 :: in_gas_caprpwr = 0d0  !power law absorption opacity density exponent
c
c
c-- output
      logical :: in_io_grabstdout = .false.  !write stdout to file
      logical :: in_io_dogrdtally = .false. !write transport tallies per grid cell
      logical :: in_io_nogriddump = .false.  !don't write grid cell variables
      character(4) :: in_io_opacdump = 'off '    !off|one|all: write opacity data to file
      character(4) :: in_io_pdensdump = 'off '   !off|on: write partial densities to file
      logical :: in_io_profdump = .false.  !write profile to file at the end
      logical :: in_io_debugging = .false. !sets nogriddump to false, opacdump,pdensdump,dogrdtally,profdump to true
c
c-- runtime parameter namelist
      namelist /inputpars/
     & in_name,in_comment,
     & in_nomp,
     & in_niter,
!grd
     & in_grd_igeom,in_ndim,
!grp
     & in_grp_ng,in_grp_ngs,in_grp_wlmin,in_grp_wlmax,
!flx
     & in_flx_ndim,in_flx_wlmin,in_flx_wlmax,
!src
     & in_src_ns,in_src_n2s,in_srcepwr,in_srcmax,in_surfsrcmu,
     & in_L_bol,in_L_inner,
!trn
     & in_trn_errorfatal,in_puretran,
     & in_trn_noamp,
     & in_trn_tauddmc,in_taulump,
     & in_thm_aniso,
!gas
     & in_noeos,
     & in_nebular,in_R_phot,
     & in_nlte,in_nlte_nelem,
!phys opac
     & in_opacmixrossel,
     & in_nobbopac,in_nobfopac,in_noffopac,in_nothmson,
!anal opac
     & in_opacanaltype,
     & in_ldisp1,in_ldisp2,
     & in_gas_sigcoef,in_gas_sigtpwr,in_gas_sigrpwr,
     & in_gas_capcoef,in_gas_captpwr,in_gas_caprpwr,
!io
     & in_io_grabstdout,
     & in_io_nogriddump,in_io_dogrdtally,
     & in_io_opacdump,in_io_pdensdump,
     & in_io_profdump,in_io_debugging
c
c-- pointers
c
      integer,parameter,private :: npointers = 100
c
      type lptr
       logical,pointer :: p
      endtype lptr
      type(lptr) :: in_l(npointers)
c
      type iptr
       integer,pointer :: p
      endtype iptr
      type(iptr) :: in_i(npointers)
c
      type rptr
       real*8,pointer :: p
      endtype rptr
      type(rptr) :: in_r(npointers)
c
      type cptr
       character(4),pointer :: p
      endtype cptr
      type(cptr) :: in_c(npointers)
c
      public
      private inputpars
      save
c
c
      contains
c
      subroutine inputpar_create_pointers(il,ii,ir,ic)
c     ------------------------------------------------
      implicit none
************************************************************************
* create pointer arrays (not to be confused with a array pointers) to
* all input parameters.  These arrays are used in mpimod to broadcast
* all input parameters at once.
************************************************************************
      integer,intent(out) :: il,ii,ir,ic
c
c-- init
      il=0
      ii=0
      ir=0
      ic=0
c
      call insertl(in_io_grabstdout,in_l,il)
      call inserti(in_nomp,in_i,ii)
      call inserti(in_niter,in_i,ii)
      call inserti(in_grd_igeom,in_i,ii)
      call inserti(in_ndim(1),in_i,ii)
      call inserti(in_ndim(2),in_i,ii)
      call inserti(in_ndim(3),in_i,ii)
      call inserti(in_flx_ndim(1),in_i,ii)
      call inserti(in_flx_ndim(2),in_i,ii)
      call inserti(in_flx_ndim(3),in_i,ii)
      call insertr(in_flx_wlmin,in_r,ir)
      call insertr(in_flx_wlmax,in_r,ir)
      call insertl(in_io_nogriddump,in_l,il)
      call insertl(in_io_dogrdtally,in_l,il)
      call inserti(in_src_ns,in_i,ii)
      call inserti(in_src_n2s,in_i,ii)
      call insertr(in_L_bol,in_r,ir)
      call insertr(in_L_inner,in_r,ir)
      call insertl(in_puretran,in_l,il)
      call insertl(in_trn_errorfatal,in_l,il)
      call insertl(in_trn_noamp,in_l,il)
      call insertr(in_trn_tauddmc,in_r,ir)
      call insertr(in_taulump,in_r,ir)
      call insertl(in_thm_aniso,in_l,il)
      call inserti(in_grp_ng,in_i,ii)
      call inserti(in_grp_ngs,in_i,ii)
      call insertr(in_grp_wlmin,in_r,ir)
      call insertr(in_grp_wlmax,in_r,ir)
      call insertr(in_opacmixrossel,in_r,ir)
      call insertc(in_opacanaltype,in_c,ic)
      call insertr(in_ldisp1,in_r,ir)
      call insertr(in_ldisp2,in_r,ir)
      call insertr(in_gas_sigcoef,in_r,ir)
      call insertr(in_gas_sigtpwr,in_r,ir)
      call insertr(in_gas_sigrpwr,in_r,ir)
      call insertr(in_gas_capcoef,in_r,ir)
      call insertr(in_gas_captpwr,in_r,ir)
      call insertr(in_gas_caprpwr,in_r,ir)
      call insertc(in_surfsrcmu,in_c,ic)
      call insertr(in_srcmax,in_r,ir)
      call insertr(in_srcepwr,in_r,ir)
      call insertc(in_io_opacdump,in_c,ic)
      call insertc(in_io_pdensdump,in_c,ic)
      call insertl(in_noeos,in_l,il)
      call insertl(in_nebular,in_l,il)
      call insertr(in_R_phot,in_r,ir)
      call insertl(in_nlte,in_l,il)
      call inserti(in_nlte_nelem,in_i,ii)
      call insertl(in_nobbopac,in_l,il)
      call insertl(in_nobfopac,in_l,il)
      call insertl(in_noffopac,in_l,il)
      call insertl(in_nothmson,in_l,il)
      call insertl(in_io_profdump,in_l,il)
      call insertl(in_io_debugging,in_l,il)
c
      contains
c
      subroutine insertl(par,arr,i)
      implicit none
      logical,intent(in),target :: par
      type(lptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'insertl: i>npointers'
      arr(i)%p => par
      end subroutine insertl
c
      subroutine inserti(par,arr,i)
      implicit none
      integer,intent(in),target :: par
      type(iptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'inserti: i>npointers'
      arr(i)%p => par
      end subroutine inserti
c
      subroutine insertr(par,arr,i)
      implicit none
      real*8,intent(in),target :: par
      type(rptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'insertr: i>npointers'
      arr(i)%p => par
      end subroutine insertr
c
      subroutine insertc(par,arr,i)
      implicit none
      character(4),intent(in),target :: par
      type(cptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'inserti: i>npointers'
      arr(i)%p => par
      end subroutine insertc

      end subroutine inputpar_create_pointers
c
c
      subroutine read_inputpars
c     -------------------------
      implicit none
************************************************************************
* read the input parameter namelist
************************************************************************
      character(15),parameter :: fname='input.par'
c
c-- read namelist
      open(4,file=fname,status='old',err=66)
      read(4,nml=inputpars,end=67,err=68)
      close(4)
c
      return
66    stop 'read_inputpars: namelist input file missing: input.par'
67    stop 'read_inputpars: namelist missing or bad in input.par'
68    stop 'read_inputpars: ivalid parameters or values in namelist'

      end subroutine read_inputpars
c
c
c
      subroutine parse_inputpars(nmpi)
c     --------------------------------
      use miscmod, only:warn
c$    use omp_lib
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* parse the input parameter namelist
************************************************************************
      integer :: istat
c$    integer :: i
c
c-- dump namelist to stdout
      write(6,*) 'namelist read:'
      write(6,nml=inputpars)
      write(6,*)
c
c-- write simulation name to file
      open(4,file='output.name',iostat=istat)
      if(istat/=0) stop 'parse_inputpars: open output.name error'
      write(4,'(a)') trim(in_name)
      close(4)
c
c-- check input parameter validity
      if(in_nomp<0) stop 'in_nomp invalid'
      if(in_nomp==0 .and. nmpi>1) stop 'no in_nomp==0 in mpi mode'
c
      if(any(in_ndim<1)) stop 'in_ndim invalid'
c
      select case(in_grd_igeom)
      case(11)
       if(in_ndim(2)>1 .or. in_ndim(3)>1) stop 'in_ndim invalid'
       if(in_flx_ndim(2)/=1) stop 'in_flx_ndim(2) inval'
       if(in_flx_ndim(3)/=1) stop 'in_flx_ndim(3) inval'
      case default
       stop 'in_grd_igeom invalid'
      endselect
c
      if(in_io_debugging) then
        in_io_nogriddump = .false.
        in_io_profdump = .true.
        in_io_opacdump = 'one'
        in_io_pdensdump = 'on'
        in_io_dogrdtally = .true.
      endif
c
      if(in_io_nogriddump .and. in_io_dogrdtally) stop
     &   'dogridtally and !griddump'
      if(in_io_nogriddump .and. in_io_opacdump.ne.'off') stop
     &   'opacdump and !griddump'
c
      if(in_grp_ng<0) stop 'in_grp_ng invalid'
      if(in_grp_ngs<=0) stop 'in_grp_ngs invalid'
      if(in_grp_wlmin<0) stop 'in_grp_wlmin invalid'
      if(in_grp_wlmax<=in_grp_wlmin) stop 'in_grp_wlmax invalid'
c
      if(in_srcepwr<=0d0) stop 'in_srcpwr <= 0'
      if(in_src_ns<=0 .eqv. in_src_n2s<0) stop
     &  'use in_src_ns or in_src_n2s'
c
      if(in_L_bol.le.0) stop 'in_L_bol<=0'
c
      if(in_taulump<=.05d0*in_trn_tauddmc) stop 'in_taulump too small' !don't let scattering dominate
c
      select case(in_opacanaltype)
      case('none')
        if(in_nobbopac.and.in_nobfopac.and.in_noffopac)
     &   call warn('inputparmod','no phys opac + in_opacanaltype==none')
      case('mono')
      case('grey')
      case('line')
      case default
       stop 'in_opacanaltype unknown'
      end select
c-- disallow physical opacities when analytic opacities are selected
      if(in_opacanaltype/='none' .and. .not.(in_nobbopac .and.
     &  in_nobfopac .and. in_noffopac)) then
       stop 'in_no??opac: no physical opac allowed with analytic opac'
      endif
c
      if(in_nobbopac) call warn('read_inputpars','bb opacity disabled!')
      if(in_nobfopac) call warn('read_inputpars','bf opacity disabled!')
      if(in_noffopac) call warn('read_inputpars','ff opacity disabled!')
      if(in_nothmson) call warn('read_inputpars','Thomson disabled')
c
      select case(in_io_opacdump)
      case('off ')
      case('one ')
      case('all ') !-- for experimental iterative method
       call warn('read_inputpars',
     &   "in_io_opacdump=='all' will generate a big data file!")
      case default
       stop 'in_io_opacdump invalid'
      end select
c
      if(in_opacmixrossel<0d0 .or. in_opacmixrossel>1d0) then
       stop 'in_opacmixrossel invalid'
      endif
c
      if(in_nebular.and.in_R_phot.le.0) stop 'in_R_phot<=0'
c
      if(in_nlte.and.in_nlte_nelem.ne.1) stop 'in_nlte_nelem!=1'
c
      select case(in_io_pdensdump)
      case('off ')
      case('on ')
      case default
       stop 'in_io_pdensdump invalid'
      end select
c
c-- set the number of threads
c$    if(.false.) then
c-- serial run
       in_nomp = 1
c$    else
c-- openmp run
c$     if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c$omp parallel shared(in_nomp) private(i)
c$     i = omp_get_num_threads()
c$     if(in_nomp/=0 .and. i/=in_nomp)
c$   &   stop 'read_inputpars: in_nomp error'
c$     in_nomp = i
c$omp end parallel
c$    endif
      write(6,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &  nmpi,in_nomp,nmpi*in_nomp
      if(in_io_grabstdout) then
       write(0,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &   nmpi,in_nomp,nmpi*in_nomp
      endif
c
      end subroutine parse_inputpars
c
c
c
      subroutine provide_inputpars(nmpi)
c     -------------------------------------
      use physconstmod
      use particlemod
      use sourcemod
      use transportmod
      use fluxmod
      use groupmod
      use gasmod
      use gridmod
      use nltemod
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* Distribute the input parameter values to the respective modules.
* This needs to be called AFTER the values are bcast.
************************************************************************
      integer :: mpart
      integer :: ns
c
      mpart = int(int(2,8)**in_src_n2s/nmpi) !max number of particles
      mpart = max(mpart,in_src_ns/nmpi)
      prt_npartmax = mpart
c
      ns = int(int(2,8)**in_src_n2s/nmpi)
      ns = max(ns,in_src_ns/nmpi)
      src_ns = ns
c
      trn_tauddmc = in_trn_tauddmc
      trn_taulump = in_taulump
      trn_thm_aniso = in_thm_aniso
      trn_errorfatal = in_trn_errorfatal
      trn_noampfact = in_trn_noamp
c
      flx_ndim  = in_flx_ndim
      flx_wlmin = in_flx_wlmin
      flx_wlmax = in_flx_wlmax
c
      grp_ng    = in_grp_ng
      grp_ngs   = in_grp_ngs
      grp_wlmin = in_grp_wlmin
      grp_wlmax = in_grp_wlmax
c
      grd_igeom = in_grd_igeom
      grd_nx    = in_ndim(1)
      grd_ny    = in_ndim(2)
      grd_nz    = in_ndim(3)
c
      nlte_nelem = in_nlte_nelem
      R_phot = in_R_phot
c
      end subroutine provide_inputpars
c
      end module inputparmod
c vim: fdm=marker
