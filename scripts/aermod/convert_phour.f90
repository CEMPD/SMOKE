program convert_phour

  use m3utilio

  implicit none
  
  character*16 phour
  integer phourout
  integer year, jdate, jtime
  integer ios, srcidx, tstep
  
  integer, dimension(:), allocatable :: indxh
  real, dimension(:), allocatable :: tmpannfac
  real, dimension(:,:), allocatable :: annfac
  
  character*16 :: progname = 'convert_phour'
  
  phour = promptmfile('PHOUR', fsread3, 'PHOUR', progname)
  phourout = promptffile('PHOUR_OUT', .false., .true., 'PHOUR_OUT', progname)
  year = envint('YEAR', 'YEAR', 2011, ios)

  jdate = year * 1000 + 1
  jtime = 0
  
  if (.not. desc3(phour)) then
    call m3exit(progname, 0, 0, 'Could not get PHOUR description', -1)
  end if
  
  if (mxrec3d < 365*24) then
    call m3exit(progname, 0, 0, 'PHOUR must contain data for 365 days', -1)
  end if

  allocate(indxh(nrows3d))
  allocate(tmpannfac(nrows3d))
  allocate(annfac(nrows3d, 366*24))
  
  if (.not. read3(phour, 'INDXH', allays3, jdate, jtime, indxh)) then
    call m3exit(progname, 0, 0, 'Could not read INDXH from PHOUR', -1)
  end if

  tstep = 1
  do while (jdate / 1000 == year)
    if (.not. read3(phour, 'STKFL', allays3, jdate, jtime, tmpannfac)) then
      call m3exit(progname, jdate, jtime, 'Could not read STKFL from PHOUR', -1)
    end if

    do srcidx = 1, nrows3d
      annfac(srcidx, tstep) = tmpannfac(srcidx)
    end do

    call nextime(jdate, jtime, 10000)
    tstep = tstep + 1
  end do
  
  do srcidx = 1, nrows3d
    write(phourout, '(I5, ",", I4, 8784(",", F9.6))') indxh(srcidx), year, annfac(srcidx, 1:tstep-1)
  end do
  
  call m3exit(progname, 0, 0, ' ', 0)

end program convert_phour
