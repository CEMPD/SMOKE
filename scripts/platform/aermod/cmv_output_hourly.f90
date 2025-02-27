program cmv_output_hourly

  use m3utilio
  
  implicit none
  
  integer, parameter :: MAX_SMOKE_IDS_PER_LINE = 20
  
  integer :: srclist, ios, nlines, nsrc, i
  integer :: aermod_id, smoke_id
  integer, dimension(MAX_SMOKE_IDS_PER_LINE) :: smoke_ids
  real :: ann_emis
  character*2 :: state, prev_state
  character*20 :: file_prefix, run_group, fac_id, src_id, poll
  character*1000 :: output_dir, outfile
  
  character*16 :: phour, infile, hourly_poll
  integer year, jdate, jtime, month, pmonth, day
  integer srcidx, tstep, ntsteps
  character*4 :: year_str
  
  character*2, dimension(:), allocatable :: states
  character*20, dimension(:), allocatable :: file_prefixes, run_groups, fac_ids, src_ids
  real, dimension(:), allocatable :: annual_emis, hourly_emis
  integer, dimension(:), allocatable :: smk_to_aermod, indxh
  real, dimension(:,:), allocatable :: aermod_emis
  
  character*16 :: progname = 'cmv_output_hour'
  
  
  srclist = getefile('CMV_SRC_LIST', .true., .true., progname)
  if (srclist <= 0) then
    call m3exit(progname, 0, 0, 'Could not open CMV_SRC_LIST file', -1)
  end if
  call envstr('OUTPUT_DIR', 'OUTPUT_DIR', '', output_dir, ios)
  call envstr('HOURLY_POLL', 'HOURLY_POLL', 'VOC', hourly_poll, ios)
  
  year = envint('YEAR', 'YEAR', 2017, ios)
  write(year_str, '(i4)') year
  ntsteps = 365*24
  
  
  ! read list of AERMOD sources
  nlines = 0
  nsrc = 0
  read(srclist, *) ! skip header line
  do
    read(srclist, *, iostat=ios) state, file_prefix, run_group, fac_id, src_id, poll, ann_emis, smoke_ids
    if (ios /= 0) then
      exit
    end if
    nlines = nlines + 1
    do i = 1, MAX_SMOKE_IDS_PER_LINE
      if (smoke_ids(i) == 0) then
        exit
      end if
      nsrc = max(smoke_ids(i), nsrc)
    end do
  end do
  rewind(srclist)
  
  allocate(states(nlines), file_prefixes(nlines), run_groups(nlines))
  allocate(fac_ids(nlines), src_ids(nlines), annual_emis(nlines))
  
  allocate(aermod_emis(nlines, ntsteps))
  aermod_emis = 0.0
  
  allocate(smk_to_aermod(nsrc))
  smk_to_aermod = 0
  
  aermod_id = 0
  read(srclist, *) ! skip header line
  do
    aermod_id = aermod_id + 1
    read(srclist, *, iostat=ios) state, file_prefix, run_group, fac_id, src_id, poll, ann_emis, smoke_ids
    if (ios /= 0) then
      exit
    end if
    
    states(aermod_id) = state
    file_prefixes(aermod_id) = file_prefix
    run_groups(aermod_id) = run_group
    fac_ids(aermod_id) = fac_id
    src_ids(aermod_id) = src_id
    annual_emis(aermod_id) = ann_emis
    
    do i = 1, MAX_SMOKE_IDS_PER_LINE
      if (smoke_ids(i) == 0) then
        exit
      end if
      smk_to_aermod(smoke_ids(i)) = aermod_id
    end do
  end do
  
  
  allocate(indxh(1))
  allocate(hourly_emis(1))
  
  jdate = year * 1000 + 1
  jtime = 0
  do tstep = 1, ntsteps
  
    call daymon(jdate, month, day)
    
    if (month /= pmonth) then
      deallocate(indxh)
      deallocate(hourly_emis)
      
      write(infile, '(a, i2.2)') 'PHOUR', month
      phour = promptmfile(infile, fsread3, infile, progname)
      
      if (.not. desc3(phour)) then
        call m3exit(progname, 0, 0, 'Could not get PHOUR description', -1)
      end if
      
      allocate(indxh(nrows3d))
      allocate(hourly_emis(nrows3d))
    end if
    
    if (.not. read3(phour, 'INDXH', allays3, jdate, jtime, indxh)) then
      call m3exit(progname, jdate, jtime, 'Could not read INDXH from PHOUR', -1)
    end if
    
    if (.not. read3(phour, hourly_poll, allays3, jdate, jtime, hourly_emis)) then
      call m3exit(progname, jdate, jtime, 'Could not read ' // hourly_poll // ' from PHOUR', -1)
    end if
    
    do srcidx = 1, nrows3d
      smoke_id = indxh(srcidx)
      
      if (smoke_id == 0) then
        exit
      end if
      
      if (smoke_id > size(smk_to_aermod)) then
        cycle
      end if
      
      aermod_id = smk_to_aermod(smoke_id)
      if (aermod_id == 0) then
        cycle
      end if
      
      aermod_emis(aermod_id, tstep) = aermod_emis(aermod_id, tstep) + max(hourly_emis(srcidx), 0.0)
    end do
    
    call nextime(jdate, jtime, 10000)
    pmonth = month
  end do
  
  
  prev_state = '00';
  do aermod_id = 1, size(run_groups)
    if (states(aermod_id) /= prev_state) then
      if (prev_state /= '00') then
        close(50)
      end if
      
      outfile = trim(output_dir) // '/temporal/' // trim(file_prefixes(aermod_id)) // '_' // trim(states(aermod_id)) // '_hourly.csv'
      open(unit=50, file=outfile, access='append')
      
      prev_state = states(aermod_id)
    end if
  
    jdate = year * 1000 + 1
    jtime = 0
    do tstep = 1, ntsteps
      call daymon(jdate, month, day)
      
      write(50, '(a, ",", a, ",", a, ",", a, ",", i0, "," i0, ",", i0, ",", e18.12)') trim(run_groups(aermod_id)), trim(fac_ids(aermod_id)), trim(src_ids(aermod_id)), year_str(3:4), month, day, jtime/10000 + 1, aermod_emis(aermod_id, tstep) / annual_emis(aermod_id)
      
      call nextime(jdate, jtime, 10000)
    end do
  end do
  
  call m3exit(progname, 0, 0, '', 0)

end program cmv_output_hourly
