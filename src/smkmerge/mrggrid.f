
        PROGRAM MRGGRID

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C    Program MRGGRID reads 2-D and 3-D I/O API files and merges them
C    into a single 2-D or 3-D file (depending on the inputs)
C    The time period merged is adjusted based on the latest
C    starting file and earliest ending file, unless MRG_DIFF_DAY is
C    set in which case the time period is based on the standard 
C    environment variables. All variables are merged, even if different 
C    variables are in each file.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    Original by M. Houyoux 4/98
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
 
C...........   MODULES for public variables
C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS, NLAYS, 
     &                     VGLVS, VGTYP, VGTOP

        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'
      
C...........   EXTERNAL FUNCTIONS
        CHARACTER(2)  CRLF
        LOGICAL       ENVYN
        INTEGER       GETFLINE
        LOGICAL       GETYN       
        INTEGER       INDEX1
        INTEGER       LBLANK
        INTEGER       PROMPTFFILE
        CHARACTER(16) PROMPTMFILE
        INTEGER       SEC2TIME
        INTEGER       SECSDIFF

        EXTERNAL CRLF, ENVYN, GETFLINE, GETYN, INDEX1, LBLANK,
     &           PROMPTFFILE, PROMPTMFILE, SEC2TIME, SECSDIFF

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Emissions arrays
        REAL, ALLOCATABLE :: E2D ( : )   ! 2-d emissions
        REAL, ALLOCATABLE :: EOUT( :,: ) ! output emissions

C...........   Input file descriptors
        INTEGER,       ALLOCATABLE :: DURATA( : ) ! no. time steps
        INTEGER,       ALLOCATABLE :: NCOLSA( : ) ! no. columns
        INTEGER,       ALLOCATABLE :: NROWSA( : ) ! no. rows
        INTEGER,       ALLOCATABLE :: NVARSA( : ) ! no. variables
        INTEGER,       ALLOCATABLE :: SDATEA( : ) ! start date
        INTEGER,       ALLOCATABLE :: STIMEA( : ) ! start time
        INTEGER,       ALLOCATABLE :: NLAYSA( : ) ! number of layers in the file
        CHARACTER(16), ALLOCATABLE :: FNAME ( : ) ! 2-d input file names
        LOGICAL,       ALLOCATABLE :: USEFIRST(:) ! true: use first time step of file

        LOGICAL,       ALLOCATABLE :: LVOUTA( :,: ) ! iff out var in input file
        CHARACTER(16), ALLOCATABLE :: VNAMEA( :,: ) ! variable names
        CHARACTER(16), ALLOCATABLE :: VUNITA( :,: ) ! variable units
        CHARACTER(80), ALLOCATABLE :: VDESCA( :,: ) ! var descrip
        CHARACTER(16)                 VNAMEP( MXVARS3 ) ! pt variable names
        CHARACTER(16)                 VUNITP( MXVARS3 ) ! pt variable units
        CHARACTER(80)                 VDESCP( MXVARS3 ) ! pt var descrip

C...........   Intermediate output variable arrays
        INTEGER       INDXN ( MXVARS3 ) ! sorting index for OUTIDX
        INTEGER       OUTIDX( MXVARS3 ) ! index to master model species list

        CHARACTER(16) OUTNAM( MXVARS3 ) ! unsorted output variable names
        CHARACTER(16) VUNITU( MXVARS3 ) ! unsorted output variable units
        CHARACTER(80) VDESCU( MXVARS3 ) ! unsorted output variable descriptions

        LOGICAL       LVOUTP( MXVARS3 ) ! iff output var exists in point input

C...........   Logical names and unit numbers

        INTEGER       IDEV            ! unit for logical names list for 2d files
        INTEGER       LDEV            ! unit for log file
        INTEGER       RDEV            ! unit for merge report file
        CHARACTER(16) ONAME           ! Merged output file name
        CHARACTER(16) PNAME           ! Point source input file name 

C...........   Other local variables 
        INTEGER       C, F, J, K, L, L1, L2, NL, V, T ! pointers and counters

        INTEGER       DUMMY                      ! dummy value for use with I/O API functions
        INTEGER       EDATE                      ! ending julian date
        INTEGER       ETIME                      ! ending time HHMMSS
        INTEGER    :: G_SDATE = 0                ! start date from environment
        INTEGER    :: G_STIME = 0                ! start time from environment
        INTEGER    :: G_NSTEPS = 1               ! number of time steps from environment
        INTEGER    :: G_TSTEP = 0                ! time step from environment
        INTEGER       IOS                        ! i/o status
        INTEGER       IREC                       ! line number count
        INTEGER       JDATE                      ! iterative julian date
        INTEGER       JTIME                      ! iterative time HHMMSS
        INTEGER       LB                         ! leading blanks counter
        INTEGER       LE                         ! location of end of string
        INTEGER       MXNFIL                     ! max no. of 2-d input files
        INTEGER       NFILE                      ! no. of 2-d input files
        INTEGER       NSTEPS                     ! no. of output time steps
        INTEGER       NVOUT                      ! no. of output variables
        INTEGER       RDATE                      ! reference date
        INTEGER       SAVLAYS                    ! number of layers
        INTEGER       SDATE                      ! starting julian date
        INTEGER       SECS                       ! tmp seconds
        INTEGER       SECSMAX                    ! seconds maximum
        INTEGER       SECSMIN                    ! seconds minimum
        INTEGER       STIME                      ! starting time HHMMSS
        INTEGER       STEPS                      ! tmp number of steps
        INTEGER       TIMET                      ! tmp time from seconds
        INTEGER       TSTEP                      ! time step
        INTEGER       VLB                        ! VGLVS3D lower bound 

        CHARACTER(16)  FDESC                     ! tmp file description
        CHARACTER(16)  NAM                       ! tmp file name
        CHARACTER(16)  VNM                       ! tmp variable name
        CHARACTER(256) LINE                      ! input buffer
        CHARACTER(256) MESG                      ! message field
        CHARACTER(15)  RPTCOL                    ! single column in report line
        CHARACTER(300) RPTLINE                   ! line of report file

        LOGICAL    :: EFLAG   = .FALSE.   ! error flag
        LOGICAL    :: FIRST3D = .TRUE.    ! true: first 3-d file not yet input
        LOGICAL    :: LFLAG   = .FALSE.   ! true iff 3-d file input
        LOGICAL    :: TFLAG   = .FALSE.   ! true: grid didn't match
        LOGICAL       MRGDIFF             ! true: merge files from different days

        CHARACTER(16) :: PROGNAME = 'MRGGRID' ! program name

C***********************************************************************
C   begin body of program MRGGRID
 
        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Read names of input files and open files
        MESG = 'Enter logical name for 2-D AND 3-D GRIDDED INPUTS list'

        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                      'FILELIST', PROGNAME   )

C.........  Get environment variables
        MESG = 'Merge files from different days into single file'
        MRGDIFF = ENVYN( 'MRG_DIFF_DAYS', MESG, .FALSE., IOS )

        IF( MRGDIFF ) THEN        
C.............  Get date and time settings from environment
            CALL GETM3EPI( -1, G_SDATE, G_STIME, G_TSTEP, G_NSTEPS )        
        END IF

C.........  Determine maximum number of input files in file
        MXNFIL = GETFLINE( IDEV, 'List of files to merge' )

C.........  Allocate memory for arrays that just depend on the maximum number
C           of files
        ALLOCATE( DURATA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DURATA', PROGNAME )
        ALLOCATE( NCOLSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCOLSA', PROGNAME )
        ALLOCATE( NROWSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NROWSA', PROGNAME )
        ALLOCATE( NLAYSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NLAYSA', PROGNAME )
        ALLOCATE( NVARSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NVARSA', PROGNAME )
        ALLOCATE( SDATEA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SDATEA', PROGNAME )
        ALLOCATE( STIMEA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STIMEA', PROGNAME )
        ALLOCATE( FNAME( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FNAME', PROGNAME )
        ALLOCATE( LVOUTA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LVOUTA', PROGNAME )
        ALLOCATE( VNAMEA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMEA', PROGNAME )
        ALLOCATE( VUNITA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITA', PROGNAME )
        ALLOCATE( VDESCA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCA', PROGNAME )
        
        IF( MRGDIFF ) THEN
            ALLOCATE( USEFIRST( MXNFIL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'USEFIRST', PROGNAME )
            USEFIRST = .TRUE.
        END IF

C.........  Allocate output layer structure
        ALLOCATE( VGLVS( 0:MXLAYS3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VGLVS', PROGNAME )
        VGLVS = 0.

C.........  Loop through input files and open them
        F = 0
        IREC = 0
        DO

C.............  Read file names - exit if read is at end of file
            READ( IDEV, 93000, END=27, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading file list at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( LINE .EQ. ' ' ) CYCLE

            F = F + 1

            IF( F .LE. MXNFIL ) THEN

                LB = LBLANK ( LINE )
                LE = LEN_TRIM( LINE )
                FNAME( F ) = LINE( LB+1:LE )

                IF ( .NOT. OPEN3( FNAME( F ), FSREAD3, PROGNAME )) THEN
 
                    MESG = 'Could not open file "' //
     &                     FNAME( F )( 1 : LEN_TRIM( FNAME(F) ) )// '".'
                    CALL M3MSG2( MESG )
                    MESG = 'Ending program "MRGGRID".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF      !  if open3() failed
            END IF

        END DO
27      CONTINUE

        NFILE = F

        IF( NFILE .GT. MXNFIL ) THEN
            WRITE( MESG,94010 )
     &        'INTERNAL ERROR: Dimension mismatch.  Input file count:',
     &        NFILE, 'program allows (MXNFIL):', MXNFIL
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSEIF( NFILE .EQ. 0 ) THEN
            MESG = 'No input files in list!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.........  Determine I/O API layer storage lower bound
        VLB = LBOUND( VGLVS3D,1 )

C.........  Get file descriptions and store for all input files
C.........  Loop through 2D input files
        NLAYS = 1
        DO F = 1, NFILE

            NAM = FNAME( F )
            IF ( .NOT. DESC3( NAM ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  NAM( 1:LEN_TRIM( NAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                NROWSA( F ) = NROWS3D
                NCOLSA( F ) = NCOLS3D
                NLAYSA( F ) = NLAYS3D
                NVARSA( F ) = NVARS3D
                SDATEA( F ) = SDATE3D
                STIMEA( F ) = STIME3D
                DURATA( F ) = MXREC3D
                
                IF( F == 1 ) TSTEP = TSTEP3D
                
                DO V = 1, NVARS3D
                    VNAMEA( V,F ) = VNAME3D( V )
                    VUNITA( V,F ) = UNITS3D( V )
                    VDESCA( V,F ) = VDESC3D( V )
                END DO
            END IF

C.............  Compare all other time steps back to first file.
C.............  They must match exactly.
            IF( TSTEP3D /= TSTEP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Time step', TSTEP3D,
     &                 'in file "' // TRIM( NAM ) //
     &                 '" is inconsistent with first file value of',
     &                 TSTEP
                CALL M3MSG2( MESG )
            END IF

C.............  Compare all other grids back to first grid.
C.............  They must match exactly.
            WRITE( FDESC, '(A,I3.3)' ) 'FILE', F
            TFLAG = .FALSE.
            CALL CHKGRID( FDESC, 'GRID', 0, TFLAG )

            IF( TFLAG ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( NAM )
                WRITE( MESG,94010 ) 'ERROR: File "' // NAM( 1:L ) //
     &            '" (NX,NY)  : (', NCOLSA( F ), ',', NROWSA( F ), ')'//
     &            CRLF() // BLANK10 // 'is inconsistent with first ' //
     &            'file (NX,NY) : (', NCOLS, ',', NROWS, ')'
                CALL M3MSG2( MESG )
            END IF

C.............  Compare layer structures for 3-d files. The number of layers do 
C               not need to match, but the layer structures do need to match.
            NLAYS = MAX( NLAYS, NLAYSA( F ) )
            IF( NLAYSA( F ) .GT. 1 ) THEN
                LFLAG = .TRUE.

C.................  For the first file that is 3-d, initialize output layer
C                   structure       
                IF ( FIRST3D ) THEN

                    NLAYS = NLAYSA( F )
                    VGTYP = VGTYP3D
                    VGTOP = VGTOP3D
                    VGLVS( 0:NLAYS ) = VGLVS3D( 0+VLB:NLAYS+VLB )   ! array
                    FIRST3D = .FALSE.

C.................  For additional 3-d files, compare the layer structures
                ELSE

C.....................  Check vertical type
                    IF( VGTYP3D .NE. VGTYP ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( NAM )
                        WRITE( MESG, 94010 ) 'ERROR: Vertical ' //
     &                         'coordinate type', VGTYP3D, 
     &                         'in file "'// NAM(1:L) //
     &                         '" is inconsistent with first 3-d'//
     &                         'file value of', VGTYP
                        CALL M3MSG2( MESG )
                    END IF

C.....................  Check vertical top
                    IF( VGTOP3D .NE. VGTOP ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( NAM )
                        WRITE( MESG, 94010 ) 'ERROR: Vertical ' //
     &                         'top value', VGTOP3D, 
     &                         'in file "'// NAM(1:L) //
     &                         '" is inconsistent with first 3-d'//
     &                         'file value of', VGTOP
                        CALL M3MSG2( MESG )
                    END IF

C.....................  Loop through layers of current file F
                    DO NL = 0, NLAYSA( F )

C.........................  For layers that are common to this file and previous
C                           files
                        IF( NL .LE. NLAYS ) THEN

                            IF( VGLVS3D( NL+VLB ) .NE. VGLVS( NL )) THEN
                                EFLAG = .TRUE.
                                L = LEN_TRIM( NAM )
                                WRITE( MESG, 94020 ) 'ERROR: Layer', NL,
     &                            'in file "'// NAM( 1:L ) // 
     &                            '" with level value', VGLVS3D(NL+VLB), 
     &                            CRLF()//BLANK10//'is inconsistent '//
     &                            'with first file value of', VGLVS(NL)
                                CALL M3MSG2( MESG )
                            END IF

C.........................  Add additional layers from current file to output 
C                           layer structure
                        ELSE
                            VGLVS( NL ) = VGLVS3D( NL+VLB )
                        END IF

                    END DO    ! End checking layers

C.....................  Reset the global number of layers as the maximum between
C                       the current file and all previous files
                    NLAYS = MAX( NLAYS, NLAYSA( F ) )

                END IF        ! End first 3-d file or not
            END IF            ! End 3-d files

        END DO                ! End loop through files

C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'Inconsistent time step, grid, or layers ' //
     &              'among the files!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( MRGDIFF ) THEN

C.............  Check that all files have same start time
            STIME = STIMEA( 1 )
            
            DO F = 2, NFILE
                NAM = FNAME( F )
            
                IF( STIMEA( F ) /= STIME ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Start time', 
     &                 STIMEA( F ), 'in file "' // TRIM( NAM ) //
     &                 '" is inconsistent with first file value of',
     &                 STIME
                    CALL M3MSG2( MESG )
                END IF
            END DO
            
            IF( EFLAG ) THEN
                MESG = 'Inconsistent start time among the files!'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check that environment settings are consistent with files
            IF( TSTEP /= G_TSTEP ) THEN
                WRITE( MESG,94010 ) 'ERROR: Value for G_TSTEP',
     &              G_TSTEP, 'is inconsistent' // CRLF() // BLANK10 //
     &              'with the time step of the files', TSTEP
                CALL M3EXIT( PROGNAME, 0 ,0, MESG, 2 )
            END IF

            IF( STIME /= G_STIME ) THEN
                WRITE( MESG,94010 ) 'ERROR: Value for G_STTIME',
     &              G_STIME, 'is inconsistent' // CRLF() // BLANK10 //
     &              'with the start time of the files', STIME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            DO F = 1, NFILE
                NAM = FNAME( F )
                
C.................  Check that all files are at least long enough to
C                   to cover requested output duration
                IF( DURATA( F ) < G_NSTEPS ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Number of time steps',
     &                  DURATA( F ), 'in file "' // TRIM( NAM ) // 
     &                  '" is insufficient' // CRLF() // BLANK10 //
     &                  'to cover the requsted number of output ' //
     &                  'time steps', G_NSTEPS
                    CALL M3MSG2( MESG )
                END IF

C.................  Check if file contains output start date and 
C                   enough data to cover output duration
                IF( .NOT. EFLAG ) THEN
                    SECS = SECSDIFF( SDATEA( F ), STIMEA( F ),
     &                               G_SDATE, G_STIME )
                    STEPS = SECS/3600

                    IF( STEPS > 0 .AND. STEPS < DURATA( F ) ) THEN
                        IF( DURATA(F) - STEPS >= G_NSTEPS ) THEN
                            USEFIRST( F ) = .FALSE.
                        ELSE
                            WRITE( MESG,94010 ) 'WARNING: File "' //
     &                          TRIM( NAM ) // '" contains the ' //
     &                          'requested output start date ', G_SDATE,
     &                          CRLF() // BLANK10 // 'and start time',
     &                          G_STIME, 'but does not contain ' //
     &                          'enough data ' 
     &                          // CRLF() // BLANK10 //
     &                          'to cover the requested number of ' //
     &                          'output time steps', G_NSTEPS, '.'
                            CALL M3MSG2( MESG )
                            WRITE( MESG,94010 ) BLANK5 //
     &                          'Data from the ' //
     &                          'start date of the file ', SDATEA( F ),
     &                          'will be used instead.'
                            CALL M3MSG2( MESG )
                        END IF
                    END IF
                END IF
            END DO
            
            IF( EFLAG ) THEN
                MESG = 'Problem with duration of files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            SDATE = G_SDATE
            NSTEPS = G_NSTEPS
        
        ELSE
        
C.............  Set start date/time as latest start date/time
            RDATE = ( SDATEA( 1 )/1000 ) * 1000 + 1  ! Set reference date

C.............  Find maximum seconds count b/w file and reference
            SECSMAX = 0
            DO F = 1, NFILE
                SECS = SECSDIFF( RDATE, STIMEA( 1 ), 
     &                           SDATEA( F ), STIMEA( F ) )
                SECSMAX = MAX( SECSMAX, SECS )

            END DO

C.............  Set latest start date/time of all files
            TIMET = SEC2TIME( SECSMAX )

            SDATE = RDATE
            STIME = STIMEA( 1 )
            CALL NEXTIME( SDATE, STIME, TIMET )

C.............  Set duration given shortest file, but initialize duration with
C             longest possible period across all durations.
            EDATE = SDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, MAXVAL( DURATA ) * 10000 )
            SECSMIN = SECSDIFF( SDATE, STIME, EDATE, ETIME ) ! for point

            DO F = 1, NFILE
                EDATE = SDATEA( F )
                ETIME = STIMEA( F )
                CALL NEXTIME( EDATE, ETIME, DURATA( F ) * 10000 )
                SECS = SECSDIFF( SDATE, STIME, EDATE, ETIME )
                SECSMIN = MIN( SECSMIN, SECS )

            END DO
           
            TIMET = SEC2TIME( SECSMIN )

            NSTEPS= TIMET / 10000 ! number of time steps 

        END IF

C.........  Build master output variables list
        NVOUT = 0

C.........  Loop through input files and build an output variable list
        DO F = 1, NFILE

C.............  Loop through variables in the files
            DO V = 1, NVARSA( F )

                VNM = VNAMEA( V,F )

C.................  Look for variable name in output list
                K = INDEX1( VNM, NVOUT, OUTNAM  )  ! look in output list

C.................  If its not in the output list, add it
                IF( K .LE. 0 ) THEN
                    NVOUT = NVOUT + 1
                    INDXN ( NVOUT ) = NVOUT
                    OUTNAM( NVOUT ) = VNM
                    VDESCU( NVOUT ) = VDESCA( V,F )
                    VUNITU( NVOUT ) = VUNITA( V,F )

C.................  If variable is in the output list, check the units
                ELSE
                    IF ( VUNITA( V,F ) .NE. VUNITU( K ) ) THEN
                        EFLAG = .TRUE.
                        L  = LEN_TRIM( VNM )
                        L1 = LEN_TRIM( VUNITA( V,F ) )
                        L2 = LEN_TRIM( VUNITU( K )   )
                        WRITE( MESG,94010 ) 'ERROR: Variable "' //
     &                         VNM( 1:L ) // '" in file', F,
     &                         'has units "'// VUNITA( V,F )( 1:L1 ) //
     &                         '"' // CRLF() // BLANK10 //
     &                         'that are inconsistent with a '//
     &                         'previous file that had units "' //
     &                         VUNITU(K)( 1:L2 )// '" for this variable'
                        CALL M3MSG2( MESG )

                    END IF  ! End check of units

                END IF      ! End variable in output list already or not

            END DO          ! End loop through variables in this file
        END DO              ! End loop through files.

C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'Inconsistent units for common variables among '//
     &             'the files!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort output variables into alphabetical order
        CALL SORTIC( NVOUT, INDXN, OUTNAM )

C.........  Set up for opening output file...
C.........  Get grid information
        IF( .NOT. DESC3( FNAME( 1 ) ) ) THEN
            MESG = 'Could not get description of file "'  //
     &              FNAME( 1 )( 1:LEN_TRIM( FNAME(1) ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        SDATE3D = SDATE
        STIME3D = STIME
        NVARS3D = NVOUT

C.........  Set up layer structure for output file
        NLAYS3D = NLAYS
        VGTOP3D = VGTOP
        VGTYP3D = VGTYP
        VGLVS3D = 0.     ! initialize array
        DO NL = 0, NLAYS
            VGLVS3D( NL+VLB ) = VGLVS( NL )
        END DO

C........  Update variable names in sorted order, and also 
C........  set up logical arrays for which files have which species
        DO V = 1, NVOUT
            VNM = OUTNAM( INDXN( V ) )

            VNAME3D( V ) = VNM                  ! store sorted output vars, etc.
            VDESC3D( V ) = VDESCU( INDXN( V ) )
            UNITS3D( V ) = VUNITU( INDXN( V ) )
            VTYPE3D( V ) = M3REAL

            DO F = 1, NFILE
                LVOUTA( V,F ) = .FALSE.

                J = INDEX1( VNM, NVARSA( F ), VNAMEA( 1,F ) )
                IF( J .GT. 0 ) LVOUTA( V,F ) = .TRUE.
                                
            END DO 
        END DO

C.........  Allocate memory for the number of grid cells and layers
        NGRID = NROWS * NCOLS
        ALLOCATE( E2D( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'E2D', PROGNAME )
        ALLOCATE( EOUT( NGRID, NLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOUT', PROGNAME )

C.........  Prompt for and open output file
        ONAME = PROMPTMFILE( 
     &          'Enter logical name for MERGED GRIDDED OUTPUT file',
     &          FSUNKN3, 'OUTFILE', PROGNAME )

C.........  Propmt for and open report file
        IF( MRGDIFF ) THEN
            RDEV = PROMPTFFILE(
     &             'Enter logical name for the MRGGRID REPORT file',
     &             .FALSE., .TRUE., 'REPMRGGRID', PROGNAME ) 

C.............  Write header line to report     
            WRITE( RPTLINE,93010 ) 'Output date'
            WRITE( RPTCOL,93010 ) 'Output time'
            RPTLINE = TRIM( RPTLINE ) // RPTCOL
            
            DO F = 1, NFILE
                NAM = FNAME( F )
                WRITE( RPTCOL,93010 ) TRIM( NAM ) // ' date'
                RPTLINE = TRIM( RPTLINE ) // RPTCOL
            END DO
            
            WRITE( RDEV,93000 ) TRIM( RPTLINE )
        END IF

C.........  Loop through hours
        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

C.............  Loop through species
            DO V = 1, NVOUT

                VNM = VNAME3D( V ) 

C.................  Output array
                EOUT = 0.   ! array

                DO F = 1, NFILE

C.....................  Set read date
                    IF( MRGDIFF .AND. USEFIRST( F ) ) THEN
                        DUMMY = 0
                        STEPS = SEC2TIME( 
     &                            SECSDIFF( 
     &                              SDATE, DUMMY, JDATE, DUMMY ) )
                        RDATE = SDATEA( F )
                        CALL NEXTIME( RDATE, DUMMY, STEPS )
                        
                    ELSE
                        RDATE = JDATE
                    END IF

C.....................  Set tmp variables
                    NAM = FNAME ( F )       ! input file name
                    NL  = NLAYSA( F )       ! number of layers

C.....................  If file has species, read (do this for all files)...
                    IF( LVOUTA( V,F ) ) THEN

C.........................  If 2-d input file, read, and add
                        IF( NL .EQ. 1 ) THEN
                            IF( .NOT. READ3( NAM, VNM, 1, RDATE,  
     &                                       JTIME, E2D          )) THEN

                                MESG = 'Could not read "' // VNM //
     &                                 '" from file "' //
     &                                 NAM( 1:LEN_TRIM( NAM ) )// '".'
                                CALL M3EXIT( PROGNAME, RDATE, JTIME, 
     &                                       MESG, 2 )
                            ENDIF

                            EOUT( 1:NGRID,1 ) = EOUT( 1:NGRID,1 ) + E2D

C.........................  If 3-d input file, allocate memory, read, and add
                        ELSE

                            DO K = 1, NL
                                IF( .NOT. READ3( NAM, VNM, K, RDATE,  
     &                                           JTIME, E2D      )) THEN

                                    MESG = 'Could not read "' // VNM //
     &                                     '" from file "' //
     &                                   NAM( 1:LEN_TRIM( NAM ) )// '".'
                                    CALL M3EXIT( PROGNAME, RDATE, JTIME,
     &                                           MESG, 2 )
                                END IF

                                EOUT( 1:NGRID,K )= EOUT( 1:NGRID,K ) + 
     &                                              E2D( 1:NGRID )
                            END DO

                        END IF  ! if 2-d or 3-d
                    END IF      ! if pollutant is in this file

C.....................  Build report line if needed
                    IF( MRGDIFF .AND. V == 1 ) THEN
                        IF( F == 1 ) THEN
                            WRITE( RPTLINE,93020 ) JDATE
                            WRITE( RPTCOL,93020 ) JTIME
                            RPTLINE = TRIM( RPTLINE ) // RPTCOL
                        END IF
                        
                        WRITE( RPTCOL,93020 ) RDATE
                        RPTLINE = TRIM( RPTLINE ) // RPTCOL
                    END IF

                END DO          ! loop over input files

C.................  Write species/hour to output file
                IF( .NOT. WRITE3( ONAME, VNM, JDATE, JTIME, EOUT )) THEN

                    MESG = 'Could not write "'// VNM// '" to file "'// 
     &                      ONAME( 1:LEN_TRIM( ONAME ) ) // '".'
     &                        
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF

            END DO   ! loop through variables

C.............  Write this time step to report
            IF( MRGDIFF ) THEN
                WRITE( RDEV,93000 ) TRIM( RPTLINE )
            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )
      
        END DO       ! loop through timesteps

C......... Normal Completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)
    
C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )
 
C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT(  A )

93010   FORMAT( A15 )

93020   FORMAT( I15 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, :, I3, :, 1X, 10 ( A, :, F8.5, :, 1X ) )

        END PROGRAM MRGGRID
