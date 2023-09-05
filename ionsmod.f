*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      module ionsmod
c     --------------
      implicit none
c
      integer,private :: nelem=0 !private copy (public in gasmod) value is obtained in ion_read_data
      integer :: ion_iionmax     !max number of ions an element has
      integer :: ion_nion=0      !total number of ions of all elements
      real*8,allocatable :: zeta_temp(:) !temperature array corresponding to zeta values for nebular approximation
c
      integer,private :: ncell=0 !number of cells
c
      type leveldata
       real*8 :: e     !ionization energy [erg]
       real*8 :: q     !partition function (sum)
       real*8 :: n     !number density
c
       integer :: nlev !number of levels
       real*8,allocatable :: elev(:) !exp(en) excitation energy [erg]
       real*8,allocatable :: glev(:) !degeneration count
       logical,allocatable :: meta(:) !metastability !nebular approximation
       real*8,allocatable :: zval(:) !zeta values !nebular approximation
      end type leveldata
c
c-- element configuration
      type elemconf
       integer :: ni   !number of ions
       real*8 :: ne    !number of electrons (relative to normalized total number of atoms (of this element))
       type(leveldata),allocatable :: i(:)
      end type elemconf
      type(elemconf),allocatable :: ion_el(:) !(nelem)
c
c-- ocupation number and g value of all ion's ground states
      type ocground
       integer :: ni   !number of ions
       real*8,allocatable :: oc(:)
       real*8,allocatable :: ginv(:)
      end type ocground
      type(ocground),allocatable :: ion_grndlev(:,:) !(nelem,gas_ncell)
c
*     private leveldata,elemconf,ocground
c
      save
c
      contains
c
c
c
      subroutine ions_alloc_grndlev(nelemin,ncellin)
c     -------------------------------------------
      implicit none
      integer,intent(in) :: nelemin,ncellin
************************************************************************
* ion_grndlev stores the occupation number density and g value of the
* ground states for all ions in all gas_vals cells.
************************************************************************
      integer :: iz,ni,i
c
      ncell = ncellin
c
      if(nelem==0) then
       nelem = nelemin
      elseif(nelem/=nelemin) then
       stop 'ions_solve_eos: nelem error'
      endif
c
      allocate(ion_grndlev(nelem,ncell))
      do i=1,ncell
       do iz=1,nelem
        ni = ion_el(iz)%ni
        ion_grndlev(iz,i)%ni = ni
        allocate(ion_grndlev(iz,i)%oc(ni))
        allocate(ion_grndlev(iz,i)%ginv(ni))
       enddo
      enddo !i
      end subroutine ions_alloc_grndlev
c
c
c
      subroutine ions_alloc_el(nelem,nion,nions,nlevel,ntemp)
c     ------------------------------------------------
      implicit none
      integer,intent(in) :: nelem,nions,ntemp
      integer,intent(in) :: nion(nelem),nlevel(nions)
************************************************************************
* allocate the ion_el structure and keep track of the size in bytes.
* This is used in MPI mode before broadcasting ion_el.
************************************************************************
      integer :: ii,iz,iion,nlev
c
c-- copy
      ion_nion = sum(nion)
c
c-- allocate temperature array
      allocate(zeta_temp(ntemp)) !nebular approximation
      zeta_temp = 0d0 !nebular approximation
c-- allocate data structure
      allocate(ion_el(nelem))
      iion = 0
      do iz=1,nelem
       allocate(ion_el(iz)%i(iz+1))
       ion_el(iz)%i%e = 0d0
       ion_el(iz)%i%q = 0d0
       ion_el(iz)%ni = nion(iz)
c-- allocate compressed data
       do ii=1,nion(iz)
        iion = iion + 1
        nlev = nlevel(iion)
        ion_el(iz)%i(ii)%nlev = nlev
        allocate(ion_el(iz)%i(ii)%elev(nlev))
        allocate(ion_el(iz)%i(ii)%glev(nlev))
        allocate(ion_el(iz)%i(ii)%meta(nlev)) !nebular approximation
        allocate(ion_el(iz)%i(ii)%zval(ntemp)) !nebular approximation
        ion_el(iz)%i(ii)%elev = 0d0
        ion_el(iz)%i(ii)%glev = 0d0
        ion_el(iz)%i(ii)%meta = .false. !nebular approximation
        ion_el(iz)%i(ii)%zval = 0d0 !nebular approximation
       enddo !ii
      enddo !iz
c
      end subroutine ions_alloc_el
c
c
c
      subroutine ions_dealloc
c     ----------------------
      implicit none
************************************************************************
* deallocate the complex ion_el datastructure, and ion_grndlev
************************************************************************
      integer :: iz,ii,i
c
      if(ion_nion>0) then
       do iz=1,nelem
        do ii=1,ion_el(iz)%ni
         deallocate(ion_el(iz)%i(ii)%elev)
         deallocate(ion_el(iz)%i(ii)%glev)
         deallocate(ion_el(iz)%i(ii)%meta) !nebular approximation
         deallocate(ion_el(iz)%i(ii)%zval) !nebular approximation
        enddo !ii
        deallocate(ion_el(iz)%i)
       enddo !iz
       deallocate(ion_el)
      endif
c
      if(ncell>0) then
       do i=1,ncell
        do iz=1,nelem
         deallocate(ion_grndlev(iz,i)%oc)
         deallocate(ion_grndlev(iz,i)%ginv)
        enddo
       enddo !i
       deallocate(ion_grndlev)
      endif
      end subroutine ions_dealloc
c
c
c
      subroutine ions_solve_eos(natomfr,temp,ndens,nelec,iconv)
c     --------------------------------------------------------
      use physconstmod
      use miscmod, only:warn
      implicit none
      real*8,intent(in) :: natomfr(nelem)
      real*8,intent(in) :: temp,ndens
      real*8,intent(inout) :: nelec
      integer,intent(out) :: iconv
************************************************************************
* solve LTE equation of state. Solve the Saha equation for all elements
* repeatedly to converge the electron density to a prescribed accuracy.
*
* input:
* - natomfr: normalized number fractions of each element
* - temp:    temperature [K]
* - ndens:   number density [cm^-3]
*
* output:
* - nelec: number of electrons relative to the total number of atoms,
*          where sum(natomfr)==1
************************************************************************
      integer,parameter :: nconv=40 !max number of convergence iterations
      real*8,parameter :: acc=1d-8 !accuracy requirment for convergence
      integer :: ii,iz,iprev,ihelp
      real*8 :: kti,sahac,sahac2 !constants
      real*8 :: help
c
      type nelec_conv !save two values for linear inter/extra-polation
       real*8 :: nel(2) !nelec
       real*8 :: err(2) !error
      end type nelec_conv
      type(nelec_conv) :: nec
      real*8 :: nelec_new,err,dxdy,ynew
c
c-- constant
      kti = 1d0/(pc_kb*temp)
      sahac = (pc_pi2*pc_me*pc_kb*temp/pc_h**2)**(1.5d0)
c
c-- default values
      nec%err(:) = (/1d10,-1d10/)
c
c-- sanity check
      if(nelem==0) stop 'ions_solve_eos: nelem error'
c
c
c-- evaluate partition functions: Q = Sum_i(g_i exp(e_i/kt)
c--
        do iz=1,nelem
         if(all(ion_el(iz)%i(:)%q==0)) then ! otherwise partition function evaluated in nebular eos subroutine
          do ii=1,ion_el(iz)%ni-1
            ion_el(iz)%i(ii)%q = sum(ion_el(iz)%i(ii)%glev(:)*
     &      exp(-kti*ion_el(iz)%i(ii)%elev(:))) !elev is h*c*chi, with [chi]=1/cm
          enddo !ii
          ion_el(iz)%i(ii)%q = 1d0
         endif
  !DEBUG>
  !     if(iz==2) then
  !      write(6,*) 'q:',ion_el(2)%i(:)%q !DEBUG
  !      do ii=1,ion_el(iz)%ni-1 !DEBUG
  !       write(6,*) 'glev:',ion_el(2)%i(ii)%glev(:)
  !      enddo
  !      do ii=1,ion_el(iz)%ni-1 !DEBUG
  !       write(6,*) 'elev:',ion_el(2)%i(ii)%elev(:)
  !      enddo
  !     endif
        enddo !iz
c
c
c-- repeat saha solver until convergence in nelec is achieved
c--
      iprev = 0
      do iconv=1,nconv !max number of convergence iterations
       sahac2 = 2d0/(nelec*ndens)*sahac !constant for saha_nelec
c-- solve saha equations for each element
       call saha_nelec(nelec_new)
       err = nelec/nelec_new - 1d0
!      write(6,*) iconv,nelec,nelec_new,err,iprev !useful debug output
c
c-- check for convergence
       if(abs(err) < acc) then
        if(iconv>1) nelec = nelec_new !leave untouched otherwise to save opacity recalculation
        exit !accuracy reached
       endif
c
c-- save results: replace the value on the same side of zero
       ihelp = maxloc(err*nec%err(:),dim=1) !equal signs yield positive number
       nec%nel(ihelp) = nelec
       nec%err(ihelp) = err
c
c-- check bracket
       if(iconv==2) then
        if(maxval(abs(nec%err),dim=1)==1d10) call warn('ions_solve_eos',
     &    'overshoot value not on other side')
       endif
c
c-- new guess
       if(iconv==1) then
c-- overshoot
        nelec = 1.5d0*nelec_new - .5d0*nelec !overshoot to obtain bracket values around 0 point
        nelec = min(3d0*nelec_new,nelec)   !limit
        nelec = max(.33d0*nelec_new,nelec) !limit
       else
c-- newton raphson, but not aways perfect, as it may repeatedly minimize one side, instead of close bracket
        ihelp = maxloc(abs(nec%err(:)),dim=1) !sides are flipping, newton-raphson works
        help = abs(nec%err(1)/nec%err(2))
        if((ihelp/=iprev) .and.
     &    help<1d5 .and. help>1d-5) then !bigger error is to too big
         ynew = 0d0 !perfect newton raphson
        else !sides are not flipping, depart from pure newton-raphson
         ynew = .001d0*nec%err(ihelp) !minimize the bigger gap irrespective of the smaller
        endif
        iprev = ihelp !remember which side was reduced
        dxdy = (nec%nel(2) - nec%nel(1))/(nec%err(2) - nec%err(1))
        nelec = nec%nel(1) - dxdy*(nec%err(1) - ynew)
c-- check new value
        if(nelec<=0d0) then
         write(6,*) 'iconv',iconv
         write(6,*) 'nelec,dxdy',nelec,dxdy
         write(6,*) 'nec'
         write(6,'(1p,4e12.4)') nec%nel,nec%err
         stop 'ions_solve_eos: nelec<=0'
        endif
c       nelec = min(1d5*nelec_new,nelec)  !limit
c       nelec = max(1d-5*nelec_new,nelec) !limit
       endif
      enddo !iconv
      if(iconv.gt.nconv) call warn('ions_solve_eos',
     &  'accuracy in nelec not reached')
c
c-- print ionization balance
c     do iz=1,nelem
c      ihelp = ion_el(iz)%ni
c      write(7,'(i3,1p,31e12.4)') iz,(ion_el(iz)%i(ii)%n,ii=1,ihelp)
c     enddo
c
      contains
c
      subroutine saha_nelec(nelec)
c     -----------------------------------
      implicit none
      real*8,intent(out) :: nelec
************************************************************************
* solve saha equations, updating the ionization balance, for a given
* nelec and temperature (both within sahac2) and return a new nelec.
************************************************************************
      integer :: ii,iz,nion,istart
      real*8 :: nsum
c
      do iz=1,nelem
       nion = ion_el(iz)%ni
c
c-- recursively compute n_(i+1)
       do istart=1,nion-1
        ion_el(iz)%i(istart)%n = 1d0 !use arbitrary start value for n_1
        do ii=istart+1,nion
         help = sahac2*exp(-kti*ion_el(iz)%i(ii-1)%e)*
     &     ion_el(iz)%i(ii-1)%n*ion_el(iz)%i(ii)%q/ion_el(iz)%i(ii-1)%q
         ion_el(iz)%i(ii)%n = help
        enddo !ii
        if(help<huge(nsum)) exit  !no overflow, we're done
       enddo !istart
c
c-- normalize n_i
       nsum = sum(ion_el(iz)%i(:nion)%n)
       if(nsum/=nsum .or. nsum>huge(nsum) .or. nsum<tiny(nsum)) then !verify nsum
        write(6,*) 'nsum,iz,nion,iconv=',nsum,iz,nion,iconv
        write(6,*) 'sahac,sahac2,kti',sahac,sahac2,kti
        write(6,*) 'nelec,ndens',nelec,ndens
        write(6,*) 'n'
        write(6,*) (ion_el(iz)%i(ii)%n,ii=1,nion)
        write(6,*) 'e'
        write(6,*) (ion_el(iz)%i(ii)%e,ii=1,nion)
        write(6,*) 'q'
        write(6,*) (ion_el(iz)%i(ii)%q,ii=1,nion)
        stop 'saha_nelec: nsum invalid'
       endif
       nsum = 1d0/nsum
       ion_el(iz)%i(:)%n = ion_el(iz)%i(:)%n*nsum
c
c-- compute number of electrons
       ion_el(iz)%ne = 0d0
       do ii=2,nion
        ion_el(iz)%ne = ion_el(iz)%ne + (ii-1)*ion_el(iz)%i(ii)%n
       enddo !ii
      enddo !iz
c
c-- compute new electron density
*     nelec = sum(natomfr(:)*ion_el(:)%ne)
      nelec = 0d0
      do iz=1,nelem
       nelec = nelec + natomfr(iz)*ion_el(iz)%ne
      enddo !iz
c
      end subroutine saha_nelec
c
      end subroutine ions_solve_eos
c
c
c-- experimental nebular approximation eos calculation
      subroutine ions_nebular(natomfr,temp,radtemp,rcell,ndens,nelec)
c     -----------------------------------
      use physconstmod
      use inputparmod
      use miscmod, only:warn
      implicit none
      real*8,intent(in) :: natomfr(nelem)
      real*8,intent(in) :: temp,radtemp,rcell,ndens
      real*8,intent(inout) :: nelec
************************************************************************
* correct the ionization balance using nebular approximation as
* described in Mazzali & Lucy 1993 A&A 279 447. Zeta values are read in
* and then interpolated for the cell (input) radiation temperature.
*
* input:
* - natomfr: normalized number fractions of each element
* - temp: gas kinetic temperature [K]
* - radtemp: radiation temperature [K]
* - rcell : cell center radius [cm]
* - ndens:   number density [cm^-3]
*
* output:
* - nelec: number of electrons relative to the total number of atoms,
*          where sum(natomfr)==1
************************************************************************
      integer,parameter :: nconv=40 !max number of convergence iterations
      real*8,parameter :: acc=1d-8 !accuracy requirment for convergence
      integer :: iconv,ihelp,iprev
      integer :: il,ii,iz,ize,istart,ntemp,niter,nion
      real*8 :: ktei,ktri,delta,chi_0,W,b,help,zeta
c-- ion population ratios
      type phi_vals ! values of phi for lte and nebular approximation
        real*8,allocatable :: phis(:) ! phi_{i,j+1/j}
      end type phi_vals
      type(phi_vals), allocatable :: phi_lte(:),phi_nlte(:),phi_cp(:)
c
      type nelec_conv !save two values for linear inter/extra-polation
       real*8 :: nel(2) !nelec
       real*8 :: err(2) !error
      end type nelec_conv
      type(nelec_conv) :: nec
      real*8 :: nelec_new,err,dxdy,ynew,ne_dens
c-- default values
      nec%err(:) = (/1d10,-1d10/)
c-- constants
      ktei = 1d0/(pc_kb*temp)
      ktri = 1d0/(pc_kb*radtemp)
      chi_0 = ion_el(1)%i(1)%e ! ionization edge for SNe II - 912 Angstrom
c-- evaluate dilution factor
      help = (1-(in_R_phot/rcell)**2)**0.5
      W = 0.5*(1-help)
c-- weighted mean departure coefficient
      b = 1d0/W
c-- evaluate partition functions:
      do iz=1,nelem
       help=0
       ion_el(iz)%i%q = 0d0
       do ii=1,ion_el(iz)%ni-1
c-- Q = Sum_i(g_i*exp(e_i/kt) [metastable] + W*Sum_i(g_i*exp(e_i/kt) [non-metastable]
         do il=1,ion_el(iz)%i(ii)%nlev-1
          if(ion_el(iz)%i(ii)%meta(il)) then
            help = help + ion_el(iz)%i(ii)%glev(il)*
     &                    exp(-ktri*ion_el(iz)%i(ii)%elev(il))
          else
            help = help + W*ion_el(iz)%i(ii)%glev(il)*
     &                    exp(-ktri*ion_el(iz)%i(ii)%elev(il))
          endif
         enddo !il
         ion_el(iz)%i(ii)%q = help
       enddo !ii
       ion_el(iz)%i(ii)%q = 1d0
      enddo !iz

c-- call LTE eos first to get LTE ion populations
      call ions_solve_eos(natomfr,temp,ndens,nelec,niter)

c-- determine zeta, delta, and then phi_nlte from phi_lte
      allocate(phi_lte(nelem)) ! lte ion population ratios
      allocate(phi_nlte(nelem)) ! nlte ion population ratios
      allocate(phi_cp(nelem)) ! cumulative product
      ntemp = size(zeta_temp)
      ize = itarget(radtemp,zeta_temp,ntemp) ! T_rad index to interpolate Zeta value
      do iz=1,nelem
        nion = ion_el(iz)%ni
        allocate(phi_lte(iz)%phis(nion))
        allocate(phi_nlte(iz)%phis(nion))
        allocate(phi_cp(iz)%phis(nion))
        phi_lte(iz)%phis = 1d0
        phi_nlte(iz)%phis = 1d0
        do istart=1,nion-1
         do ii=istart+1,nion
           phi_lte(iz)%phis(ii) = ion_el(iz)%i(ii)%n * nelec * ndens
     &                            / ion_el(iz)%i(ii-1)%n
           zeta = interp1d(ize,radtemp,zeta_temp,ion_el(iz)%i(ii)%zval)
           if(ion_el(iz)%i(ii)%e.ge.chi_0) then
              help = exp(ion_el(iz)%i(ii)%e*(ktri - ktei))
              delta = (temp*help)/(b*W*radtemp)
           else
              help = exp(ion_el(iz)%i(ii)%e*ktri - chi_0*ktei)
              help = (temp*help)/(b*W*radtemp)
              help = help - exp((ion_el(iz)%i(ii)%e - chi_0)*ktri)
              delta = 1 + help
           endif
           help = W*(delta*zeta+W*(1-zeta))*((temp/radtemp)**0.5)
           phi_nlte(iz)%phis(ii) = help * phi_lte(iz)%phis(ii)
         enddo !ii
        enddo !istart
      enddo !iz

c-- repeat until convergence in nelec is achieved
c--
      iprev = 0
      do iconv=1,nconv !max number of convergence iterations
c-- solve saha equations for each element
       ne_dens = nelec*ndens ! electron density
       call nebular_nelec(nelec_new)
       err = nelec/nelec_new - 1d0
!      write(6,*) iconv,nelec,nelec_new,err,iprev !useful debug output
c
c-- check for convergence
       if(abs(err) < acc) then
        if(iconv>1) nelec = nelec_new !leave untouched otherwise to save opacity recalculation
        exit !accuracy reached
       endif
c
c-- save results: replace the value on the same side of zero
       ihelp = maxloc(err*nec%err(:),dim=1) !equal signs yield positive number
       nec%nel(ihelp) = nelec
       nec%err(ihelp) = err
c
c-- check bracket
       if(iconv==2) then
        if(maxval(abs(nec%err),dim=1)==1d10) call warn('ions_solve_eos',
     &    'overshoot value not on other side')
       endif
c
c-- new guess
       if(iconv==1) then
c-- overshoot
        nelec = 1.5d0*nelec_new - .5d0*nelec !overshoot to obtain bracket values around 0 point
        nelec = min(3d0*nelec_new,nelec)   !limit
        nelec = max(.33d0*nelec_new,nelec) !limit
       else
c-- newton raphson, but not aways perfect, as it may repeatedly minimize one side, instead of close bracket
        ihelp = maxloc(abs(nec%err(:)),dim=1) !sides are flipping, newton-raphson works
        help = abs(nec%err(1)/nec%err(2))
        if((ihelp/=iprev) .and.
     &    help<1d5 .and. help>1d-5) then !bigger error is to too big
          ynew = 0d0 !perfect newton raphson
        else !sides are not flipping, depart from pure newton-raphson
          ynew = .001d0*nec%err(ihelp) !minimize the bigger gap irrespective of the smaller
        endif
        iprev = ihelp !remember which side was reduced
        dxdy = (nec%nel(2) - nec%nel(1))/(nec%err(2) - nec%err(1))
        nelec = nec%nel(1) - dxdy*(nec%err(1) - ynew)
c-- check new value
        if(nelec<=0d0) then
          write(6,*) 'iconv',iconv
          write(6,*) 'nelec,dxdy',nelec,dxdy
          write(6,*) 'nec'
          write(6,'(1p,4e12.4)') nec%nel,nec%err
          stop 'ions_solve_eos: nelec<=0'
        endif
c       nelec = min(1d5*nelec_new,nelec)  !limit
c       nelec = max(1d-5*nelec_new,nelec) !limit
       endif
      enddo !iconv
      if(iconv.gt.nconv) call warn('ions_solve_eos',
     &  'accuracy in nelec not reached')

      contains

      subroutine nebular_nelec(nelec)
      implicit none
      real*8,intent(out) :: nelec
************************************************************************
* solve for nelec, updating the ionization balance using nebular
* approximation and return a new nelec.
************************************************************************
      integer :: ii,iz,nion,istart
      real*8 :: help,nsum

      do iz=1,nelem
       nion = ion_el(iz)%ni
c-- compute cumulative product
       phi_cp(iz)%phis = 1d0
       do istart=1,nion-1
         do ii=istart+1,nion
           phi_cp(iz)%phis(ii) = phi_nlte(iz)%phis(ii)
     &                         * phi_cp(iz)%phis(ii-1) / ne_dens
           if(phi_cp(iz)%phis(ii)/=phi_cp(iz)%phis(ii).or.
     &        phi_cp(iz)%phis(ii)>huge(phi_cp(iz)%phis(ii)).or.
     &        phi_cp(iz)%phis(ii)<tiny(phi_cp(iz)%phis(ii))) then
c-- set cumulative product to 0 if NaN or inf
                phi_cp(iz)%phis(ii)=0d0
           endif
         enddo !ii
       enddo !istart
c-- compute n_i
       nsum = sum(ion_el(iz)%i(:nion)%n)
       do ii=1,nion
         help = nsum * phi_cp(iz)%phis(ii)
     &            / (sum(phi_cp(iz)%phis(:)))
         ion_el(iz)%i(ii)%n = help
       enddo !ii
c
c-- normalize n_i
       nsum = sum(ion_el(iz)%i(:nion)%n)
       if(nsum/=nsum .or. nsum>huge(nsum) .or. nsum<tiny(nsum)) then !verify nsum
         write(6,*) 'nsum,iz,nion,iconv=',nsum,iz,nion
         write(6,*) 'nelec,ndens',nelec,ndens
         write(6,*) 'n'
         write(6,*) (ion_el(iz)%i(ii)%n,ii=1,nion)
         write(6,*) 'e'
         write(6,*) (ion_el(iz)%i(ii)%e,ii=1,nion)
         write(6,*) 'q'
         write(6,*) (ion_el(iz)%i(ii)%q,ii=1,nion)
         stop 'nebular_nelec: nsum invalid'
       endif
       nsum = 1d0/nsum
       ion_el(iz)%i(:)%n = ion_el(iz)%i(:)%n*nsum
c
c-- compute number of electrons
       ion_el(iz)%ne = 0d0
       do ii=2,nion
         ion_el(iz)%ne = ion_el(iz)%ne + (ii-1)*ion_el(iz)%i(ii)%n
       enddo !ii
      enddo !iz
c
c-- compute new electron density
*     nelec = sum(natomfr(:)*ion_el(:)%ne)
      nelec = 0d0
      do iz=1,nelem
       nelec = nelec + natomfr(iz)*ion_el(iz)%ne
      enddo !iz
      end subroutine nebular_nelec

      pure integer function itarget(target,arr,size)
c     ----------------------------------
      integer, intent(in) :: size
      real*8,intent(in) :: target
      real*8,intent(in) :: arr(size)
************************************************************************
* search value closest to, but less than the target in an ordered array
************************************************************************
      integer :: il,ih,im
      il = 1
      ih = size
c
      do while(il<ih)
        im = (il+ih)/2
        if(target>=arr(im).and.target<arr(im+1)) then
          itarget = im
          return
        elseif(target>=arr(il).and.target<arr(im)) then
          ih = im
          im = (il+ih)/2
        elseif(target>=arr(im).and.target<arr(ih)) then
          il = im
          im = (il+ih)/2
        else
          if(target<arr(il)) then
            itarget = il
          else
            itarget = ih
          endif
        endif
      enddo

      end function itarget

      function interp1d(ix,x,x_arr,y_arr) result(y)
c     ----------------------------------
      integer, intent(in) :: ix
      real*8,intent(in) :: x,x_arr(:),y_arr(:)
      real*8 :: y
****************************************************************
* 1d interpolation, assumes monotonicity
****************************************************************
      real*8 :: help
c
      help = (x-x_arr(ix)) * (y_arr(ix+1)-y_arr(ix))
      help = help / (x_arr(ix+1)-x_arr(ix))

      y = y_arr(ix) + help

      end function interp1d

      end subroutine ions_nebular
c
c
c
      subroutine ions_read_data(nelem_in)
c     ----------------------------------
      use physconstmod
      use bbxsmod
      implicit none
      integer,intent(in) :: nelem_in
************************************************************************
* read ionization and level data from 'iondata.dat'
*
* exitation/ionization potentials are converted from [cm^-1] to [erg].
************************************************************************
      character(8),parameter :: fname='data.ion'
      integer :: l,ll,iz,ii,icod,nlevel
      integer :: ntemp,nzeta !nebular approximation
      integer :: istat
      real*8 :: help
      character(80) :: line
      character(80) :: line1
      real*8,allocatable :: rawzeta(:) !nebular approximation
c
c-- raw level data
      type raw_level_data
       real*8 :: chi
       character(12) :: label
       integer :: ilevel
       integer :: g
       logical :: meta !nebular approximation
      end type raw_level_data
      type(raw_level_data),allocatable :: rawlev(:)
c
c-- local private copy
      nelem = nelem_in
c
c-- allocate data structure
      allocate(ion_el(nelem))
      do iz=1,nelem
       ion_el(iz)%ni = iz+1 !include the bare nucleus
       allocate(ion_el(iz)%i(iz+1))
       forall(l=1:iz+1)
        ion_el(iz)%i(l)%nlev = 0
        ion_el(iz)%i(l)%e = 0d0
        ion_el(iz)%i(l)%q = 0d0
       endforall
      enddo
c
c-- open file
      open(4,file=fname,action='read',status='old',iostat=istat)
      if(istat/=0) stop 'ion_read_data: cannot read data.ion'
c
c
c-- read ions one by one
      do ll=1,1000
c
c-- read header
       read(4,*,iostat=istat) icod !ion code (HI = 0100)
       if(istat/=0) exit
       read(4,'(a80,a40)') line,line1 !contains ionization potential
       read(4,*)
       read(4,*) nlevel
c
c-- allocate temporary raw level data arrays
       allocate(rawlev(nlevel))
c
c-- read raw data
       read(4,'(f11.3,1x,a12,i5,i5,l6)') rawlev
c
c-- skip unused elements
       iz = icod/100  !element number
       ii = icod - 100*iz + 1 !ion number, starting from 1 for charge neutral atom, like in spectroscopic notation
       if(iz>nelem) then
        deallocate(rawlev)
        cycle
       endif
c
c-- extract ionization energy from 'line'
       if(index(line,'(Ky): ')>0) then
         l = index(line,'(Ky): ')
         read(line(l+6:),*) help
       elseif(index(line1,'(Ky): ')>0) then
         l = index(line1,'(Ky): ')
         read(line1(l+6:),*) help
       endif
       ion_el(iz)%i(ii)%e = pc_h*pc_c*help
c-- allocate data
       ion_el(iz)%i(ii)%nlev = nlevel
       allocate(ion_el(iz)%i(ii)%elev(nlevel))
       allocate(ion_el(iz)%i(ii)%glev(nlevel))
       allocate(ion_el(iz)%i(ii)%meta(nlevel))
c
c-- store data
       do l=1,nlevel
        ion_el(iz)%i(ii)%elev(l) = pc_h*pc_c*rawlev(l)%chi
        ion_el(iz)%i(ii)%glev(l) = rawlev(l)%g
        ion_el(iz)%i(ii)%meta(l) = rawlev(l)%meta
       enddo !l
c-- read ions done
       deallocate(rawlev)
      enddo !ll
c
c-- done
      close(4)
c
c-- open zeta data file !nebular approximation
      open(4,file='data.zeta',action='read',status='old',iostat=istat)
      if(istat/=0) stop 'ion_read_data: cannot read data.zeta'
c-- read header
      read(4,*)
      read(4,*) ntemp
c-- allocate temperature array for zeta values
      allocate(zeta_temp(ntemp))
c-- read tempurature data
      read(4,'(f9.1)') zeta_temp
c-- allocate zeta values
      do iz=1,nelem
        do ii=1,iz+1
          allocate(ion_el(iz)%i(ii)%zval(ntemp))
          ion_el(iz)%i(ii)%zval = 0d0
        enddo
      enddo
c-- read zeta values
      read(4,*) nzeta
      read(4,*)
      allocate(rawzeta(ntemp))
      do l=1,nzeta
        read(4,'(1x,i4,i4,*(1x,f7.4))') iz,ii,rawzeta
        ion_el(iz)%i(ii)%zval = rawzeta
      enddo
c-- done
      close(4)
c
c
c-- remove empty ions
      ion_iionmax = 0
      ion_nion = 0
      do iz=1,nelem
       do l=1,iz+1
        if(ion_el(iz)%i(l)%nlev==0) exit !this is considered the bare nucleus. todo: does this make sense?
       enddo
       ion_el(iz)%ni = l !cut off empty ions
       ion_iionmax = max(ion_iionmax,l)
       ion_nion = ion_nion + l
      enddo !iz
c
c-- asign a single state to the uppermost continuum
      do iz=1,nelem
       l = ion_el(iz)%ni !last ionization stage
       allocate(ion_el(iz)%i(l)%glev(1))
       allocate(ion_el(iz)%i(l)%elev(1))
       allocate(ion_el(iz)%i(l)%meta(1)) !nebular approximation
       ion_el(iz)%i(l)%nlev = 1
       ion_el(iz)%i(l)%glev(1) = 1d0
       ion_el(iz)%i(l)%elev(1) = 0d0
       ion_el(iz)%i(l)%meta(1) = .false. !nebular approximation
      enddo

c
c
c-- output
      write(8,*)
      write(8,*) 'ion_read_data:'
      write(8,*) '---------------------------'
      write(8,*) '# levels per ion'
      write(8,'(3x,"|",31i5)') (l,l=1,ion_iionmax)
      write(8,'(a3,"|",31a5)') '---',('-----',l=1,ion_iionmax)
      do iz=1,nelem
       write(8,'(i3,"|",91i5)') iz,
     &   (ion_el(iz)%i(l)%nlev,l=1,ion_el(iz)%ni)
      enddo
      write(8,*)
      write(8,*) 'ionization energy per ion [eV]'
      write(8,'(3x,"|",31i12)') (l,l=1,ion_iionmax)
      write(8,'(a3,"|",31a12)') '---',('------------',l=1,ion_iionmax)
      do iz=1,nelem
       write(8,'(i3,"|",1p,91e12.4)') iz,
     &   (ion_el(iz)%i(l)%e/pc_ev,l=1,ion_el(iz)%ni) !in eV units, because it's so common
      enddo
c
      end subroutine ions_read_data
c
      end module ionsmod
c vim: fdm=marker
