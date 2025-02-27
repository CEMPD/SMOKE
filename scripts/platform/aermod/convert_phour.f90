program convert_phour

  use m3utilio

  implicit none
  
  character*16 phour, infile
  integer phourout
  integer year, jdate, jtime, month, pmonth, day
  integer ios, srcidx, tstep, t
  
  integer, dimension(:), allocatable :: indxh
  real, dimension(:), allocatable :: tmpannfac
  real, dimension(:,:), allocatable :: annfac
  
  character*16 :: progname = 'convert_phour'
  
  phour = promptmfile('PHOUR01', fsread3, 'PHOUR01', progname)
  phourout = promptffile('PHOUR_OUT', .false., .true., 'PHOUR_OUT', progname)
  year = envint('YEAR', 'YEAR', 2011, ios)

  jdate = year * 1000 + 1
  jtime = 0
  
  if (.not. desc3(phour)) then
    call m3exit(progname, 0, 0, 'Could not get PHOUR description', -1)
  end if
  allocate(indxh(nrows3d))
  allocate(tmpannfac(nrows3d))
  allocate(annfac(nrows3d, 366*24))
  
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

    if (.not. read3(phour, 'HOURACT', allays3, jdate, jtime, tmpannfac)) then
      call m3exit(progname, jdate, jtime, 'Could not read HOURACT from PHOUR', -1)
    end if

    do srcidx = 1, nrows3d
      annfac(srcidx, tstep) = tmpannfac(srcidx)
    end do

    call nextime(jdate, jtime, 10000)
    tstep = tstep + 1
    pmonth = month

  end do
 
  do srcidx = 1, nrows3d
    write(phourout, '(I5, ",", I4, 8784(",", E14.7))') indxh(srcidx), year, annfac(srcidx, 1:tstep-1)
  end do
  
  call m3exit(progname, 0, 0, ' ', 0)

end program convert_phour
