*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine banner
c     -----------------
      use mpimod
      implicit none
************************************************************************
* Print banner, start date/time, and code revision.
************************************************************************
      character*8 :: t_startdate
      character*10 :: t_starttime
c
      character*(MPI_MAX_PROCESSOR_NAME) :: pname
      integer :: ilen,ierr
c
      character*40 :: coderev_id='r1.0.0'
      character*28 :: build_date=''
c      include 'version.inc'
c
      call date_and_time(t_startdate,t_starttime)
      call mpi_get_processor_name(pname,ilen,ierr)
c
      write(6,'(13(1x,"===",a62,"===",/))')
     & "==============================================================",
     & "==============================================================",
     & "                                                              ",
     & "  #####                              #                        ",
     & " #     # #    # #####  ###### #####  #       # ###### ######  ",
     & " #       #    # #    # #      #    # #       #    #   #       ",
     & "  #####  #    # #    # #####  #    # #       #    #   #####   ",
     & "       # #    # #####  #      #####  #       #    #   #       ",
     & " #     # #    # #      #      #   #  #       #    #   #       ",
     & "  #####   ####  #      ###### #    # ######  #    #   ######  ",
     & "                                                              ",
     & "==============================================================",
     & "=============================================================="
c
      write(6,*) "Postprocessing Monte Carlo Spectral Synthesis Code"
      write(6,*) "by Gururaj A. Wagle"
      write(6,*) "Released (2023/Jun/30)"
      write(6,*) "Based on Radiation Transport Code SuperNu"
      write(6,*) "by Ryan T. Wollaeger and Daniel R. van Rossum"

c
      write(6,*) "code revision id: ", trim(coderev_id)
c
c     write(6,*) "build hostname  : ", trim(build_hostname)
c     write(6,*) "compiler        : ", trim(compiler_version)
c     write(6,*) "compiler flags  : ", trim(compiler_flags)
c
      write(6,*) 'simulation date : ', t_startdate//" / "//t_starttime
      write(6,*) 'processor name  : ', trim(pname)
      write(6,*)
c
      end subroutine banner
c vim: fdm=marker
