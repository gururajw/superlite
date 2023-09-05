*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine read_nlte_data(nelem,nlte_nel)
c     --------------------------------
      use physconstmod
      use ionsmod
      use bbxsmod
      use nltemod
      use timingmod
      use inputparmod
      use miscmod, only:warn
      implicit none
      integer,intent(in) :: nelem,nlte_nel
************************************************************************
* read all bound-bound cross sections
************************************************************************
      real,parameter :: fconst = sngl(pc_pi*pc_e**2/(pc_me*pc_c))
      real*8,parameter :: Econst = 4d0*pc_pi*fconst !nLTE
      ! real,parameter :: C0 = 5.465e-11 ! pi*a0^2*(8k/(m*pi))^0.5 !nlTE
      ! real*8,parameter :: IH = 2.178d-11 !Ionization energy of hydrogen !nLTE
      ! real*8 :: Elu
      integer :: nlinall,ilinall,nlinnlte,ilinnlte
      integer :: il,l,iz,ii,istat,llw,lhg,n!,n1
      real*8 :: t0,t1
c-- quick exit
      if(nlte_nel.gt.1) then
       stop 'read_nlte_data: no NLTE calculations for Z>1'
       return
      endif
c
c-- allocate permanent NLTE data structure
      allocate(nlte_data(nlte_nel))
      nlte_nelem = nlte_nel
      do iz=1,nlte_nel
        allocate(nlte_data(iz)%i(iz)) !exclude bare nucleus
      enddo

c-- determine total number of lines
c-- count NLTE lines
      nlinnlte = 0
      do iz=1,nlte_nel
       do ii=1,iz
        call read_atom_nlte(iz,ii,istat,get_data=.false.)
c-- test if succesfull
        if(istat/=0) stop "read_atom_nlte: read error"
        write(8,*) 'read_nlte successful:',iz,ii,
     &     nlte_nlevel,nlte_nline+coll_nline
        nlinnlte = nlinnlte + nlte_nline + coll_nline
        nlte_nspecies = nlte_nspecies + 1
c-- allocate permanent storage for NLTE data
        nlte_data(iz)%i(ii)%nlev=nlte_nlevel
        nlte_data(iz)%i(ii)%nlin=nlte_nline
        nlte_data(iz)%i(ii)%coll_nlin=coll_nline
        !-- level data
        allocate(nlte_data(iz)%i(ii)%chilev(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%levid(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%glev(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%nplev(nlte_nlevel))
        nlte_data(iz)%i(ii)%chilev = 0d0
        nlte_data(iz)%i(ii)%levid = 0
        nlte_data(iz)%i(ii)%glev = 0
        nlte_data(iz)%i(ii)%nplev = 0
        !-- radiative line data
        allocate(nlte_data(iz)%i(ii)%wl0(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%f(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%gl(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%gu(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%llw(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%lup(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%Aul(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%Bul(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%Blu(nlte_nline))
        allocate(nlte_data(iz)%i(ii)%qlu(nlte_nline))
        nlte_data(iz)%i(ii)%wl0 = 0d0
        nlte_data(iz)%i(ii)%f = 0d0
        nlte_data(iz)%i(ii)%gl = 0
        nlte_data(iz)%i(ii)%gu = 0
        nlte_data(iz)%i(ii)%llw = 0
        nlte_data(iz)%i(ii)%lup = 0
        nlte_data(iz)%i(ii)%Aul = 0d0
        nlte_data(iz)%i(ii)%Bul = 0d0
        nlte_data(iz)%i(ii)%Blu = 0d0
        nlte_data(iz)%i(ii)%qlu = 0d0
        !-- EIE line data
        allocate(nlte_data(iz)%i(ii)%llw_coll(coll_nline))
        allocate(nlte_data(iz)%i(ii)%lup_coll(coll_nline))
        allocate(nlte_data(iz)%i(ii)%C0(coll_nline))
        allocate(nlte_data(iz)%i(ii)%C1(coll_nline))
        allocate(nlte_data(iz)%i(ii)%C2(coll_nline))
        allocate(nlte_data(iz)%i(ii)%C3(coll_nline))
        allocate(nlte_data(iz)%i(ii)%n_coll(coll_nline))
        allocate(nlte_data(iz)%i(ii)%delE(coll_nline))
        nlte_data(iz)%i(ii)%llw_coll = 0
        nlte_data(iz)%i(ii)%lup_coll = 0
        nlte_data(iz)%i(ii)%C0 = 0d0
        nlte_data(iz)%i(ii)%C1 = 0d0
        nlte_data(iz)%i(ii)%C2 = 0d0
        nlte_data(iz)%i(ii)%C3 = 0d0
        nlte_data(iz)%i(ii)%n_coll = 0d0
        nlte_data(iz)%i(ii)%delE = 0d0
        !-- PI/RR transition data
        allocate(nlte_data(iz)%i(ii)%lev_PI(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C0_PI(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C1_PI(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C2_PI(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C3_PI(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%n_PI(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C0_RR(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C1_RR(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C2_RR(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%C3_RR(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%n_RR(nlte_nlevel))
        allocate(nlte_data(iz)%i(ii)%delE_PI(nlte_nlevel))
        nlte_data(iz)%i(ii)%lev_PI = 0
        nlte_data(iz)%i(ii)%C0_PI = 0d0
        nlte_data(iz)%i(ii)%C1_PI = 0d0
        nlte_data(iz)%i(ii)%C2_PI = 0d0
        nlte_data(iz)%i(ii)%C3_PI = 0d0
        nlte_data(iz)%i(ii)%n_PI = 0d0
        nlte_data(iz)%i(ii)%C0_RR = 0d0
        nlte_data(iz)%i(ii)%C1_RR = 0d0
        nlte_data(iz)%i(ii)%C2_RR = 0d0
        nlte_data(iz)%i(ii)%C3_RR = 0d0
        nlte_data(iz)%i(ii)%n_RR = 0d0
        nlte_data(iz)%i(ii)%delE_PI = 0d0
       enddo !ii
      enddo !iz
c-- count LTE lines
      nlinall = 0
      do iz=nlte_nel+1,nelem
       do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
        call read_atom(iz,ii,istat,get_data=.false.)
c-- test if succesfull
        if(istat/=0) cycle
        write(8,*) 'read_atom successful:',iz,ii,bb_nlevel,bb_nline
        nlinall = nlinall + bb_nline
       enddo !ii
      enddo !iz
      write(6,'(a,i8)') ' total number of lines read:',nlinall+nlinnlte
      write(8,'(a,i8)') ' total number of lines read:',nlinall+nlinnlte
c
c-- allocate permanent storage space for LTE line data
      if(nlinall+nlinnlte<=0) stop 'rd_bbxs_data: no sigle line read in'
      allocate(bb_xs(nlinall))
      n = int(sizeof(bb_xs)/1024) !kB
      !n1 = int(sizeof(nlte_data)/1024) !kB
      write(6,*) 'ALLOC bb_xs    :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      !write(6,*) 'ALLOC nlte_data:',n1,"kB",n1/1024,"MB",n1/1024**2,"GB"
c
      t0 = t_time()
c-- store data in permanent arrays for NLTE species
      ilinnlte = 0
      do iz=1,nlte_nel
       do ii=1,iz
        call read_atom_nlte(iz,ii,istat,get_data=.true.)
        do l=1,nlte_nlevel
c-- level data
         nlte_data(iz)%i(ii)%chilev(l) = nlte_level(l)%chi
         nlte_data(iz)%i(ii)%levid(l) = nlte_level(l)%id
         nlte_data(iz)%i(ii)%glev(l) = nlte_level(l)%g
         nlte_data(iz)%i(ii)%nplev(l) = nlte_level(l)%n
c-- PI/RR transition data
         nlte_data(iz)%i(ii)%lev_PI(l) = pi_rr_data(l)%lev
         nlte_data(iz)%i(ii)%C0_PI(l) = pi_rr_data(l)%C0_PI
         nlte_data(iz)%i(ii)%C1_PI(l) = pi_rr_data(l)%C1_PI
         nlte_data(iz)%i(ii)%C2_PI(l) = pi_rr_data(l)%C2_PI
         nlte_data(iz)%i(ii)%C3_PI(l) = pi_rr_data(l)%C3_PI
         nlte_data(iz)%i(ii)%n_PI(l) = pi_rr_data(l)%n_PI
         nlte_data(iz)%i(ii)%C0_RR(l) = pi_rr_data(l)%C0_RR
         nlte_data(iz)%i(ii)%C1_RR(l) = pi_rr_data(l)%C1_RR
         nlte_data(iz)%i(ii)%C2_RR(l) = pi_rr_data(l)%C2_RR
         nlte_data(iz)%i(ii)%C3_RR(l) = pi_rr_data(l)%C3_RR
         nlte_data(iz)%i(ii)%n_RR(l) = pi_rr_data(l)%n_RR
         nlte_data(iz)%i(ii)%delE_PI(l) = pi_rr_data(l)%delE
        enddo !l !nlte_level
c-- radiative line data
        do il=1,nlte_nline
         ilinnlte = ilinnlte + 1
         !-- level number for lower and upper levels of transition
         llw = nlte_line(il)%lev1
         lhg = nlte_line(il)%lev2
         !-- flip low<->high levels
         if(llw>lhg) then
          llw = lhg
          lhg = nlte_line(il)%lev1
         endif
         !-- wavelength for transition
         if((abs(nlte_level(lhg)%chi) - abs(nlte_level(llw)%chi)/=0d0))
     &    nlte_data(iz)%i(ii)%wl0(il) = 1e8/(abs(nlte_level(lhg)%chi) - !in ang
     &      abs(nlte_level(llw)%chi))
         if((nlte_data(iz)%i(ii)%wl0(il).le.0d0)) then
           nlte_data(iz)%i(ii)%wl0(il) = nlte_line(il)%wl0
         endif
c
c
         nlte_data(iz)%i(ii)%llw(il) = llw
         nlte_data(iz)%i(ii)%lup(il) = lhg
         !-- g of lower and upper level of transition
         nlte_data(iz)%i(ii)%gl(il) = nlte_level(llw)%g
         nlte_data(iz)%i(ii)%gu(il) = nlte_level(lhg)%g
         !-- oscillator strength l->u
         nlte_data(iz)%i(ii)%f(il) = nlte_line(il)%f
         !-- Einstein coefficients A & B
         nlte_data(iz)%i(ii)%Aul(il) =
     &      (2d0*Econst*nlte_level(llw)%g*nlte_line(il)%f)
     &       /(nlte_level(lhg)%g*(nlte_line(il)%wl0*pc_ang)**2)
         nlte_data(iz)%i(ii)%Bul(il) =
     &      (Econst*nlte_line(il)%wl0*pc_ang*nlte_level(llw)%g
     &       *nlte_line(il)%f)/(pc_h*pc_c*nlte_level(lhg)%g)
         nlte_data(iz)%i(ii)%Blu(il) =
     &      (Econst*nlte_line(il)%wl0*pc_ang*nlte_line(il)%f)
     &       /(pc_h*pc_c)
c-- Colisional excitaion, van Regemorter formula, eqn 9.58 of Hubeny and Mihalas (2014)
!         Elu = pc_h*pc_c/(nlte_line(il)%wl0*pc_ang)
!         nlte_data(iz)%i(ii)%qlu(il) = !remaining terms included in physical_opacity routine
!     &      C0*14.5*10**nlte_line(il)%f*(IH/Elu)**2
        enddo !il !nlte_line
        deallocate(nlte_level,nlte_line)
c-- EIE line data
        do il=1,coll_nline
         ilinnlte = ilinnlte + 1
         !-- level number for lower and upper levels of transition
         llw = coll_line(il)%lev1
         lhg = coll_line(il)%lev2
         !-- flip low<->high levels
         if(llw>lhg) then
          llw = lhg
          lhg = coll_line(il)%lev1
         endif
c
         nlte_data(iz)%i(ii)%llw_coll(il) = llw
         nlte_data(iz)%i(ii)%lup_coll(il) = lhg
         !-- energy of transition in eV
         nlte_data(iz)%i(ii)%delE(il) = coll_line(il)%delE
         !-- fit coefficients
         nlte_data(iz)%i(ii)%C0(il) = coll_line(il)%C0
         nlte_data(iz)%i(ii)%C1(il) = coll_line(il)%C1
         nlte_data(iz)%i(ii)%C2(il) = coll_line(il)%C2
         nlte_data(iz)%i(ii)%C3(il) = coll_line(il)%C3
         !-- pre-factor coefficient (currently not used)
         nlte_data(iz)%i(ii)%n_coll(il) = coll_line(il)%n_coll
        enddo !il
        deallocate(coll_line)
       enddo !ii
      enddo !iz
c-- store data in permanent array for LTE species
      ilinall = 0
      do iz=nlte_nel+1,nelem
       do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
        call read_atom(iz,ii,istat,get_data=.true.)
        if(istat/=0) cycle
c
        do l=1,bb_nline
         ilinall = ilinall + 1
         llw = bbxs_line(l)%lev1
         lhg = bbxs_line(l)%lev2
         !-- line center wavelength
         bb_xs(ilinall)%wl0 = 1e8/(abs(bbxs_level(lhg)%chi) - !in ang
     &     abs(bbxs_level(llw)%chi))
         !-- flip low<->high levels
         if(bb_xs(ilinall)%wl0 < 0.) then
          llw = lhg
          lhg = bbxs_line(l)%lev1
          bb_xs(ilinall)%wl0 = -bb_xs(ilinall)%wl0
         endif
         !-- g*xs
         bb_xs(ilinall)%gxs = fconst*bbxs_level(llw)%g*
     &     10.**bbxs_line(l)%f          !fconst = pi*e**2/(m_e*c)
         !-- exp(chi)
         bb_xs(ilinall)%chilw = abs(bbxs_level(llw)%chi)
         !-- ion code
         bb_xs(ilinall)%iz = int(iz,2)
         bb_xs(ilinall)%ii = int(ii,2)
        enddo !l
        !-- ready with raw data
        deallocate(bbxs_level,bbxs_line)
       enddo !ii
      enddo !iz
c-- verify counter
      if(ilinnlte/=nlinnlte)
     &  stop 'rd_nlte_data: NLTE line counter error'
      if(ilinall/=nlinall) stop 'rd_nlte_data: LTE line counter error'
c
c-- store counter globally
      bb_nline = nlinall
c
c-- sort lines - doesn't speed-up bb opacity.
      call sort_lines
c
      t1 = t_time()
!     write(6,'(a,f8.2,a)') ' time used for nlte reading:',t1-t0,'s'
c
      end subroutine read_nlte_data
c vim: fdm=marker
