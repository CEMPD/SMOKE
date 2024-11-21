        program mult
!
!   Purpose:           multiply two IOAPI files together using a one file as a scalar
!
!
!******** environmental variables used
!
!        
!         
!        FILENAME
!        
!

         implicit none
      
         include 'PARMS3.EXT'      ! I/O API constants
         include 'FDESC3.EXT'      ! I/O API file description data structure
         include 'IODECL3.EXT'     ! I/O API function declarations 

       integer logdev
       character(len=16) :: this_program, input, scalar
       character(len=16) :: output
       character(len=16) :: scalar_name
       

       integer nsteps,nvars, nlays, ncols, nrows, tstep, sdate, stime
       integer istatus
       integer :: c_loop, r_loop, v_loop



       real, allocatable :: buffer(:,:,:,:), scalar_var(:,:,:)
       

       integer i_loop, t_loop, jdate, jtime, l_loop



       character*80 mesg

       integer trimlen, envint,getefile
       external trimlen, envint, getefile
       real envreal
       external envreal


       logdev = init3()	!  initialization returns unit # for log

       this_program = 'MULT'
       scalar        = 'INFILE1'
       input         = 'INFILE2'
       output       = 'OUTFILE'
       
       scalar_name  = 'xportfrac'


        if ( .not. open3( input, FSREAD3, this_program ) ) THEN
           MESG = 'Could not open file "' //
     &     input( 1: TRIMLEN(input)) 
     &     // '" for input'
           CALL M3EXIT( this_program, 0, 0, MESG, 2 )
        end if        

        IF ( .NOT. DESC3(input))THEN
           MESG = 'Could not get description info for file "' //
     &              input( 1: TRIMLEN( input) ) //'"'
           CALL M3EXIT( this_program, 0, 0, MESG, 2 )     
        ENDIF

        nsteps  = MXREC3D
	nvars   = NVARS3D
        nlays   = NLAYS3D
	ncols   = NCOLS3D
	nrows   = NROWS3D
	sdate   = SDATE3D
	stime   = STIME3D
	tstep   = TSTEP3D	

	jdate    = sdate
	jtime    = stime


	
           allocate(buffer(ncols,nrows,nlays,nvars))



        if ( .not. open3( output, FSUNKN3, this_program ) ) THEN
           MESG = 'Could not open file "' //
     &     output( 1: TRIMLEN(output)) 
     &     // '" for output'
           CALL M3EXIT( this_program, 0, 0, MESG, 2 )
        end if 


        if ( .not. open3( scalar, FSREAD3, this_program ) ) THEN
           MESG = 'Could not open file "' //
     &     scalar( 1: TRIMLEN(scalar)) 
     &     // '" for input'
           CALL M3EXIT( this_program, 0, 0, MESG, 2 )
        end if  
            
        allocate (scalar_var(ncols,nrows,nlays))

       
	   if ( .not. 
     &        read3(scalar,scalar_name,ALLAYS3,0,0,scalar_var)
     &        ) then
              mesg = 'Error reading from file '//
     &           scalar( 1: TRIMLEN( scalar ) ) 
              call m3exit( this_program, 0, 0, MESG, 2 )
           end if 


		
	

        do t_loop = 1,nsteps        
	   if ( .not. 
     &        read3(input,ALLVAR3,ALLAYS3,jdate,jtime,buffer)
     &        ) then
              mesg = 'Error reading from file '//
     &           input( 1: TRIMLEN( input ) ) 
              call m3exit( this_program, 0, 0, MESG, 2 )
           end if 

          do v_loop = 1, nvars
	     do l_loop = 1, nlays
	        do r_loop = 1, nrows
		   do c_loop = 1, ncols

                buffer(c_loop,r_loop,l_loop,v_loop) = 
     &     scalar_var(c_loop,r_loop,l_loop)*buffer(c_loop,r_loop,l_loop,v_loop)
                   enddo
		enddo

             enddo
	  enddo 

	   if ( .not. 
     &        write3(output,ALLVAR3,jdate,jtime,buffer)
     &        ) then
              mesg = 'Error writing to file '//
     &           output( 1: TRIMLEN( output ) ) 
              call m3exit( this_program, 0, 0, MESG, 2 )
           end if 

           call nextime(jdate, jtime, tstep)
        enddo




	if (.not. shut3()) then
	   write (*,*) 'UNRECOVERABLE ERROR: could not do a shut3' 
	endif

	end
