*This file is part of SuperLite.
      subroutine output_grp
c     --------------------------
      use inputparmod
      use groupmod
      use gridmod
      use timingmod
      implicit none
************************************************************************
* Write flux arrays to file.
************************************************************************
      integer :: i
      character(16),save :: fstat='replace'
      integer :: reclen
      real*8 :: t0,t1
c
      t0 = t_time()
c--
      reclen = (grp_ng+1)*12
      open(unit=4,file='output.grp_grid',status=fstat,
     &    recl=reclen)
c-- header: dimension
      write(4,'("#",i5)') grp_ng
c-- body
      write(4,'(1p,10000e12.4)') grp_wl(:)
      close(4)
c
      if(trim(in_io_opacdump)=='off') then
      ! c-- donothing
      else
c-- write multigroup opacity to a file
        open(4,file='output.grp_cap',status=fstat,
     &      recl=reclen)
c-- write header
        write(4,'("#",a20)') 'Multigroup opacity'
        write(4,'("#",2a5)') 'ng','nr'
        write(4,'("#",2i5)') grp_ng,grd_ncell
        do i=1,grd_ncell
          write(4,'(1p,10000es12.4e3)') grd_cap(:grp_ng,i)
        enddo !i
        close(4)
c
c-- write multigroup emission opacity to a file
        open(4,file='output.grp_capemit',status=fstat,
     &      recl=reclen)
c-- write header
        write(4,'("#",a30)') 'Multigroup "emission" opacity'
        write(4,'("#",2a5)') 'ng','nr'
        write(4,'("#",2i5)') grp_ng,grd_ncell
        do i=1,grd_ncell
          write(4,'(1p,10000es12.4e3)') grd_capemit(:grp_ng,i)
        enddo !i
        close(4)
c
c-- write multigroup emission opacity to a file
        open(4,file='output.grp_emiss',status=fstat,
     &      recl=reclen)
c-- write header
        write(4,'("#",a30)') 'Multigroup emissivity'
        write(4,'("#",2a5)') 'ng','nr'
        write(4,'("#",2i5)') grp_ng,grd_ncell
        do i=1,grd_ncell
          write(4,'(1p,10000es12.4e3)') grd_emiss(:grp_ng,i)
        enddo !i
        close(4)
      endif
c-- timing
      t1 = t_time()
      call timereg(t_output, t1-t0)
c
      end subroutine output_grp
c vim: fdm=marker
