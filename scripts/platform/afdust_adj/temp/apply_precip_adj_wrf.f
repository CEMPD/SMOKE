        program apply_precip_adj

      use M3UTILIO

c kathy 17:
c ifort -fPIC -check -traceback -extend_source -zero -O3 -mp1 -o apply_precip_adj_wrf.int17 apply_precip_adj_wrf.f -L/home/local-rhel7/apps/netcdf-4.4.1/intel-17.0/lib -lnetcdf -lnetcdff -L/home/local-rhel7/apps/ioapi-3.2/intel-17.0/lib -lioapi -module /home/local-rhel7/apps/ioapi-3.2/intel-17.0/Linux2_x86_64ifort
c james 19: 
c ifort -Bstatic -static-intel -I/work/EMIS/smoke/test/repo/ioapi-3.2/ioapi/fixed_src -I/work/EMIS/smoke/test/ioapi -extend-source 132 -auto -zero -xHost -traceback -c ./apply_precip_adj_wrf.f
c ifort -Bstatic -static-intel -debug -o apply_precip_adj_wrf.x apply_precip_adj_wrf.o -extend-source 132 -auto -zero -xHost -traceback -L/work/EMIS/smoke/test/ioapi -L/work/EMIS/smoke/test/repo/netcdf/lib/ -lioapi -lnetcdff -lnetcdf

!
!   Purpose:           apply precip adjustment factor to afdust emissions using METCRO2D
!
!
!

         implicit none
      
c embedded within 'use M3UTILIO'      
C          include 'PARMS3.EXT'      ! I/O API constants
C          include 'FDESC3.EXT'      ! I/O API file description data structure
C          include 'IODECL3.EXT'     ! I/O API function declarations 

       integer logdev
       character(len=16) :: progname, infile, metcro2d, outfile

       INTEGER :: wrfvers = 0

       integer nsteps,nvars, nlays, ncols, nrows, tstep, sdate, stime
       integer istatus
       integer :: c_loop, r_loop, v_loop



       real, allocatable :: buffer(:,:,:,:)
       real, allocatable :: snow_cover(:,:)
       real, allocatable :: soim1(:,:),  sltyp(:,:) 
       
       real :: zero, ratio
       integer i_loop, t_loop, jdate, jtime, l_loop, edate, etime

       real, allocatable :: totals(:)

       character*80 mesg
       INTEGER, PARAMETER :: max_wrf_def = 11
       INTEGER, PARAMETER :: max_wrf_v4  = 12


C Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
C       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
C See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
        REAL sat_wrf_def( max_wrf_def )
        DATA sat_wrf_def/ 0.395, 0.410, 0.435, 0.485,
     &                    0.451, 0.420, 0.477, 0.476,
     &                    0.426, 0.482, 0.482        /

c wrf version 4 values
        REAL sat_wrf_v4( max_wrf_v4 )
        DATA sat_wrf_v4/ 0.395, 0.410, 0.435, 0.485,
     &                   0.480, 0.451, 0.420, 0.477,
     &                   0.476, 0.426, 0.482, 0.482  /

       INTEGER MAXSTYPES,IOS
       real, allocatable  :: SATURATION(:)

c now from 'use M3UTILIO' statement
C        integer trimlen, envint,getefile
C        external trimlen, envint, getefile
C        real envreal
C        external envreal

       integer trimlen

       logdev = init3()	!  initialization returns unit # for log

       wrfvers = ENVINT('WRF_VERSION', 'WRF version', 0, IOS )
       if (IOS .gt. 0) then
          print*,'must setenv WRF_VERSION'
          print*,'abort'
          goto 999
       endif

       if ( wrfvers .eq. 4 ) then
          MAXSTYPES = max_wrf_v4
          allocate(SATURATION(MAXSTYPES))
          SATURATION = sat_wrf_v4
       else
          MAXSTYPES = max_wrf_def
          allocate(SATURATION(MAXSTYPES))
          SATURATION = sat_wrf_def
       endif


       progname = 'APPLY_PRECIP_ADJ'
       infile        = 'INFILE'
       metcro2d      = 'METCRO2D'
       outfile       = 'OUTFILE'
       zero = 0.0D0
       
        if ( .not. open3( metcro2d, FSREAD3, progname ) ) THEN
           MESG = 'Could not open file "' //
     &     metcro2d( 1: TRIMLEN(metcro2d)) 
     &     // '" for reading'
           CALL M3EXIT( progname, 0, 0, MESG, 2 )
        end if 
	
        IF ( .NOT. DESC3(metcro2d))THEN
           MESG = 'Could not get description info for file "' //
     &              metcro2d( 1: TRIMLEN( metcro2d) ) //'"'
           CALL M3EXIT( progname, 0, 0, MESG, 2 )     
        ENDIF	
        jdate = SDATE3D
	jtime = STIME3D
            
	    

        if ( .not. open3( infile, FSREAD3, progname ) ) THEN
           MESG = 'Could not open file "' //
     &     infile( 1: TRIMLEN(infile)) 
     &     // '" for input'
           CALL M3EXIT( progname, 0, 0, MESG, 2 )
        end if        

        IF ( .NOT. DESC3(infile))THEN
           MESG = 'Could not get description info for file "' //
     &              infile( 1: TRIMLEN( infile) ) //'"'
           CALL M3EXIT( progname, 0, 0, MESG, 2 )     
        ENDIF

        nsteps  = MXREC3D
	nvars   = NVARS3D
        nlays   = NLAYS3D
	ncols   = NCOLS3D
	nrows   = NROWS3D

	tstep   = TSTEP3D	

        edate   = SDATE3D
	etime   = STIME3D


	
        allocate(buffer(ncols,nrows,nlays,nvars))
        allocate(rn(ncols, nrows))
	allocate(rc(ncols, nrows))
	allocate(snow_cover(ncols, nrows))
	allocate(totals(nvars))
	allocate(soim1(ncols, nrows))
	allocate(soit1(ncols, nrows))
	allocate(wspd10(ncols, nrows))
	allocate(sltyp(ncols, nrows))
	allocate(lai(ncols, nrows))
	
	




        SDATE3D = jdate
	STIME3D = jtime
	NVARS3D = nvars
	NLAYS3D = nlays
	NCOLS3D = ncols
	NROWS3D = nrows
	TSTEP3D = tstep

	    
	    
        if ( .not. open3( outfile, FSUNKN3, progname ) ) THEN
           MESG = 'Could not open file "' //
     &     outfile( 1: TRIMLEN(outfile)) 
     &     // '" for output'
           CALL M3EXIT( progname, 0, 0, MESG, 2 )
        end if 


		
	

        do t_loop = 1,nsteps        
	   if ( .not. 
     &        read3(infile,ALLVAR3,ALLAYS3,edate,etime,buffer)
     &        ) then
              mesg = 'Error reading from file '//
     &           infile( 1: TRIMLEN( infile ) ) 
              call m3exit( progname, 0, 0, MESG, 2 )
           end if 
           do v_loop = 1, nvars
	       totals(v_loop) = SUM(buffer(1:ncols,1:nrows,1,v_loop))
	   enddo


	   if ( .not. 
     &        read3(metcro2d,'SNOCOV',1,jdate,jtime,snow_cover)
     &        ) then
              mesg = 'Error reading from file '//
     &           metcro2d( 1: TRIMLEN( metcro2d ) ) 
              call m3exit( progname, 0, 0, MESG, 2 )
           end if 

	   	   	   
	   if ( .not. 
     &        read3(metcro2d,'SOIM1',1,jdate,jtime,soim1)
     &        ) then
              mesg = 'Error reading from file '//
     &           metcro2d( 1: TRIMLEN( metcro2d ) ) 
              call m3exit( progname, 0, 0, MESG, 2 )
           end if 


	   if ( .not. 
     &        read3(metcro2d,'SLTYP',1,jdate,jtime,sltyp)
     &        ) then
              mesg = 'Error reading from file '//
     &           metcro2d( 1: TRIMLEN( metcro2d ) ) 
              call m3exit( progname, 0, 0, MESG, 2 )
           end if 


           do v_loop = 1, 1
!	       write (*,'(A1,1X,I6,1X,F6.4)') 'R',jtime, SUM(buffer(1:ncols,1:nrows,1,v_loop))/totals(v_loop)
	   enddo


          do v_loop = 1, nvars
	     do l_loop = 1, nlays
	        do r_loop = 1, nrows
		   do c_loop = 1, ncols

                   if (snow_cover(c_loop,r_loop) .eq. 1.0) then
                       buffer(c_loop,r_loop,l_loop,v_loop) = 
     &            zero*buffer(c_loop,r_loop,l_loop,v_loop)
		   endif
		   
		   
		   
                   enddo
		enddo

             enddo
	  enddo 
	  
	  

           do v_loop = 1, 1
	       write (logdev,'(A1,1X,I6,1X,F6.4)') 'S',jtime, SUM(buffer(1:ncols,1:nrows,1,v_loop))/totals(v_loop)
	   enddo
	   

           do v_loop = 1, 1
!	       write (*,'(A1,1X,I6,1X,F6.4)') 'W',jtime, SUM(buffer(1:ncols,1:nrows,1,v_loop))/totals(v_loop)
	   enddo
	   

          do v_loop = 1, nvars
	     do l_loop = 1, nlays
	        do r_loop = 1, nrows
		   do c_loop = 1, ncols
		       if ((soim1(c_loop,r_loop) .lt. 1.0) .and. (soim1(c_loop,r_loop) .gt. 0.0)) then  ! check for lai > 0 to make sure we are over land
		          if ((sltyp(c_loop,r_loop) .ge. 1) .and. (sltyp(c_loop,r_loop)) .le. maxstypes) then
		             ratio = soim1(c_loop,r_loop)/saturation(INT(sltyp(c_loop,r_loop)))

                             if (ratio .ge. 0.5) then
                                buffer(c_loop,r_loop,l_loop,v_loop) = zero*buffer(c_loop,r_loop,l_loop,v_loop)
		             endif
			     
		   
		          endif
	               endif
                   enddo
		enddo

             enddo
	  enddo 
	  
	  

           do v_loop = 1, 1
	       write (logdev,'(A1,1X,I6,1X,F6.4)') 'M',jtime, SUM(buffer(1:ncols,1:nrows,1,v_loop))/totals(v_loop)
	   enddo
	   	   
	   
	   if ( .not. 
     &        write3(outfile,ALLVAR3,jdate,jtime,buffer)
     &        ) then
              mesg = 'Error writing to file '//
     &           outfile( 1: TRIMLEN( outfile ) ) 
              call m3exit( progname, 0, 0, MESG, 2 )
           end if 

           call nextime(jdate, jtime, tstep)
	   call nextime(edate, etime, tstep)
        enddo




 999    continue
	if (.not. shut3()) then
	   write (logdev,*) 'UNRECOVERABLE ERROR: could not do a shut3' 
	endif

	end
