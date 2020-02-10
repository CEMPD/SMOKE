program convert_phour_cmv

  use m3utilio

  implicit none
  
  character*16 phour, infile
  integer phourout
  integer year, jdate, jtime, month, pmonth, day
  integer ios, srcidx, tstep, t
  
  integer, dimension(:), allocatable :: indxh
  real, dimension(:), allocatable :: tmpco
  real, dimension(:,:), allocatable :: co
  
  character*16 :: progname = 'convert_phour_cmv'
  
  phour = promptmfile('PHOUR01', fsread3, 'PHOUR01', progname)
  phourout = promptffile('PHOUR_OUT', .false., .true., 'PHOUR_OUT', progname)
  year = envint('YEAR', 'YEAR', 2017, ios)

  jdate = year * 1000 + 1
  jtime = 0
  
  if (.not. desc3(phour)) then
    call m3exit(progname, 0, 0, 'Could not get PHOUR description', -1)
  end if
  
  allocate(indxh(nrows3d))
  allocate(tmpco(nrows3d))
  allocate(co(nrows3d, 366*24))
  
  if (.not. read3(phour, 'INDXH', allays3, jdate, jtime, indxh)) then
    call m3exit(progname, 0, 0, 'Could not read INDXH from PHOUR', -1)
  end if

  if( .not. close3(phour) ) then
   call m3exit(progname, 0, 0, 'Could not close PHOUR', -1)
  end if

  tstep = 1
  do t = 1,365*24

    call daymon( jdate, month, day)

    if( month /= pmonth ) then
        write(infile,'(a,i2.2)') 'PHOUR',month
        phour = promptmfile(infile, fsread3, infile, progname)

        if (.not. desc3(phour)) then
           call m3exit(progname, 0, 0, 'Could not get PHOUR description', -1)
        end if
    end if

    if (.not. read3(phour, 'CO', allays3, jdate, jtime, tmpco)) then
      call m3exit(progname, jdate, jtime, 'Could not read CO from PHOUR', -1)
    end if

    do srcidx = 1, nrows3d
      co(srcidx, tstep) = max(tmpco(srcidx), 0.0)
    end do

    call nextime(jdate, jtime, 10000)
    tstep = tstep + 1
    pmonth = month

  end do
  
  do srcidx = 1, nrows3d
    write(phourout, '(I5, ",", I4, 8784(",", E14.7))') indxh(srcidx), year, co(srcidx, 1:tstep-1)
  end do
  
  call m3exit(progname, 0, 0, ' ', 0)

end program convert_phour_cmv
