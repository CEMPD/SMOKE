
        PROGRAM MRGGRID

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C    Program MRGGRID reads 2-D area, biogenic, mobile, and 3-D
C    point source emissions and merges into a single 3-D file.
C    The time period merged is adjusted based on the latest
C    starting file and earliest ending file.  All variables are
C    merged, even if different variables are in each file.
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
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
 
        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'
      
C...........   EXTERNAL FUNCTIONS
        CHARACTER*2   CRLF
        LOGICAL       ENVYN
        INTEGER       GETFLINE
        LOGICAL       GETYN       
        INTEGER       INDEX1
        INTEGER       LBLANK
        INTEGER       PROMPTFFILE
        CHARACTER*16  PROMPTMFILE
        INTEGER       SEC2TIME
        INTEGER       SECSDIFF

        EXTERNAL CRLF, ENVYN, GETFLINE, GETYN, INDEX1, LBLANK,
     &           PROMPTFFILE, PROMPTMFILE, SEC2TIME, SECSDIFF

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C...........   LOCAL VARIABLES and their descriptions:

C...........   Emissions arrays
        REAL, ALLOCATABLE :: E2D( :,: ) ! all 2-d emissions
        REAL, ALLOCATABLE :: E3D( :,: ) ! 3-d emissions

C...........   Input file descriptors
        INTEGER     , ALLOCATABLE :: DURATA( : ) ! 2-d no. time steps
        INTEGER     , ALLOCATABLE :: NCOLSA( : ) ! 2-d no. columns
        INTEGER     , ALLOCATABLE :: NROWSA( : ) ! 2-d no. rows
        INTEGER     , ALLOCATABLE :: NVARSA( : ) ! 2-d no. variables
        INTEGER     , ALLOCATABLE :: SDATEA( : ) ! 2-d start date
        INTEGER     , ALLOCATABLE :: STIMEA( : ) ! 2-d start time
        CHARACTER*16, ALLOCATABLE :: FNAME ( : ) ! 2-d input file names
        INTEGER       DURATP           ! pt no. time steps
        INTEGER       NCOLSP           ! pt no. columns
        INTEGER    :: NLAYSP = 0       ! pt no. layers
        INTEGER       NROWSP           ! pt no. rows
        INTEGER       NVARSP           ! pt no. variables
        INTEGER       SDATEP           ! pt start date
        INTEGER       STIMEP           ! pt start time

        LOGICAL     , ALLOCATABLE :: LVOUTA( :,: ) ! iff out var in input file
        CHARACTER*16, ALLOCATABLE :: VNAMEA( :,: ) ! 2-d variable names
        CHARACTER*16, ALLOCATABLE :: VUNITA( :,: ) ! 2-d variable units
        CHARACTER*80, ALLOCATABLE :: VDESCA( :,: ) ! 2-d var descrip
        CHARACTER*16                 VNAMEP( MXVARS3 ) ! pt variable names
        CHARACTER*16                 VUNITP( MXVARS3 ) ! pt variable units
        CHARACTER*80                 VDESCP( MXVARS3 ) ! pt var descrip

C...........   Intermediate output variable arrays
        INTEGER       INDXN ( MXVARS3 ) ! sorting index for OUTIDX
        INTEGER       OUTIDX( MXVARS3 ) ! index to master model species list

        CHARACTER*16  OUTNAM( MXVARS3 ) ! unsorted output variable names
        CHARACTER*16  VUNITU( MXVARS3 ) ! unsorted output variable units
        CHARACTER*80  VDESCU( MXVARS3 ) ! unsorted output variable descriptions

        LOGICAL       LVOUTP( MXVARS3 ) ! iff output var exists in point input

C...........   Logical names and unit numbers

        INTEGER       IDEV            ! unit for logical names list for 2d files
        INTEGER       LDEV     
        CHARACTER*16  ONAME           ! Merged output file name
        CHARACTER*16  PNAME           ! Point source input file name 

C...........   Other local variables 
        INTEGER       C, F, J, K, L, V, T ! pointers and counters

        INTEGER       EDATE                      ! ending julian date
        INTEGER       ETIME                      ! ending time HHMMSS
        INTEGER       IOS                        ! i/o status
        INTEGER       IREC                       ! line number count
        INTEGER       JDATE                      ! iterative julian date
        INTEGER       JTIME                      ! iterative time HHMMSS
        INTEGER       LB                         ! leading blanks counter
        INTEGER       LE                         ! location of end of string
        INTEGER       MXNFIL                     ! max no. of 2-d input files
        INTEGER       NCOLS                      ! no. grid columns
        INTEGER       NF2D                       ! no. of 2-d input files
        INTEGER       NGRID                      ! no. grid cells
        INTEGER       NLAYS                      ! no. layers
        INTEGER       NROWS                      ! no. grid rows
        INTEGER       NSTEPS                     ! no. of output time steps
        INTEGER       NVOUT                      ! no. of output variables
        INTEGER       RDATE                      ! reference date
        INTEGER       SDATE                      ! starting julian date
        INTEGER       SECS                       ! tmp seconds
        INTEGER       SECSMAX                    ! seconds maximum
        INTEGER       SECSMIN                    ! seconds minimum
        INTEGER       STIME                      ! starting time HHMMSS
        INTEGER       TIMET                      ! tmp time from seconds
        INTEGER       TSTEP                      ! time step

        CHARACTER*16  NAM                        ! tmp file name
        CHARACTER*16  VNM                        ! tmp variable name
        CHARACTER*256 LINE                       ! input buffer
        CHARACTER*256 MESG                       ! message field

        LOGICAL    :: EFLAG = .FALSE.            ! error flag
        LOGICAL    :: LFLAG = .FALSE.            ! true iff 3-d file input

        CHARACTER*16  :: PROGNAME = 'MRGGRID' ! program name

C***********************************************************************
C   begin body of program MRGGRID
 
        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Retrieve values of environment variables
        LFLAG   = ENVYN( 'MRG_LAYERS_YN', 
     &                   'Use layer fractions or not', .FALSE., IOS )

C.........  Read names of input files and open files
        MESG = 'Enter logical name for 2-D GRIDDED INPUTS list'

        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                      'NLST2D', PROGNAME   )
        IF( LFLAG ) THEN
            PNAME = PROMPTMFILE(
     &         'Enter logical name for 3-D GRIDDED INPUT file',
     &         FSREAD3, 'LGTS', PROGNAME )
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

        NF2D = F

        IF( NF2D .GT. MXNFIL ) THEN
            WRITE( MESG,94010 )
     &        'INTERNAL ERROR: Dimension mismatch.  Input file count:',
     &        NF2D, 'program (MXNFIL):', MXNFIL
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSEIF( .NOT. LFLAG .AND. NF2D .EQ. 0 ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, 'No input files!', 2 )

        ENDIF

C.........  Get file descriptions and store for all input files
C.........  Loop through 2D input files
        DO F = 1, NF2D

            NAM = FNAME( F )
            IF ( .NOT. DESC3( NAM ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  NAM( 1:LEN_TRIM( NAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                NROWSA( F ) = NROWS3D
                NCOLSA( F ) = NCOLS3D
                NVARSA( F ) = NVARS3D
                SDATEA( F ) = SDATE3D
                STIMEA( F ) = STIME3D
                DURATA( F ) = MXREC3D
                DO V = 1, NVARS3D
                    VNAMEA( V,F ) = VNAME3D( V )
                    VUNITA( V,F ) = UNITS3D( V )
                    VDESCA( V,F ) = VDESC3D( V )
                END DO
            END IF

        END DO

        TSTEP  = TSTEP3D     ! Set for all program

C.........  Compare grids from all files with each other
        NROWS = NROWSA( 1 )
        NCOLS = NCOLSA( 1 )
        DO F = 1, NF2D

            IF( NROWSA( F ) .NE. NROWS .OR. 
     &          NCOLSA( F ) .NE. NCOLS      ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'File "' // FNAME( F ) //
     &            '" (NX,NY)  : (', NCOLSA( F ), ',', NROWSA( F ), ')'
                CALL M3MSG2( MESG )
            ENDIF

        END DO

C.........  3-D SOURCES - get description and compare grids
        IF( LFLAG ) THEN

            IF ( .NOT. DESC3( PNAME ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  PNAME( 1:LEN_TRIM( PNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSEIF( TSTEP3D .NE. TSTEP ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 'Bad point time step', 2 )

            ELSE
                NROWSP = NROWS3D
                NCOLSP = NCOLS3D
                NLAYSP = NLAYS3D
                NVARSP = NVARS3D
                SDATEP = SDATE3D
                STIMEP = STIME3D
                DURATP = MXREC3D
                DO V = 1, NVARSP
                    VNAMEP( V ) = VNAME3D( V )
                    VUNITP( V ) = UNITS3D( V )
                    VDESCP( V ) = VDESC3D( V )
                END DO
            ENDIF

            IF( NROWSP .NE. NROWS .OR. NCOLSP .NE. NCOLS ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &               'Point     (X,Y,Z): (', NCOLSP, ',', NROWSP, ',',
     &               NLAYSP, ')'
                CALL M3MSG2( MESG )
            ENDIF

        ELSE

            SDATEP = SDATEA( 1 )  ! For initializing references below
            STIMEP = STIMEA( 1 )
            DURATP = DURATA( 1 )

        END IF
        
        IF( EFLAG ) THEN
            WRITE( MESG,94010 )
     &           'First file (X,Y): (', NCOLS, ',', NROWS, ')' //
     &           CRLF() // BLANK5 //
     &           'The grids are inconsistent among the files!'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, 'Bad input', 2 )

        END IF

C.........  Set start date/time as latest start date/time

        RDATE = ( SDATEA( 1 )/1000 ) * 1000 + 1  ! Set reference date

C.........  Find maximum seconds count b/w file and reference, using point
C.........  sources to initialize maximum 

        SECSMAX = SECSDIFF( RDATE, STIMEA( 1 ), SDATEP, STIMEP )
        DO F = 1, NF2D
            SECS = SECSDIFF( RDATE, STIMEA( 1 ), 
     &                       SDATEA( F ), STIMEA( F ) )
            SECSMAX = MAX( SECSMAX, SECS )

        END DO

C.........  Set latest start date/time of all files
        TIMET = SEC2TIME( SECSMAX )

        SDATE = RDATE
        STIME = STIMEA( 1 )
        CALL NEXTIME( SDATE, STIME, TIMET )

C.........  Set duration given shortest file - initialize w/ point

        EDATE = SDATEP
        ETIME = STIMEP
        CALL NEXTIME( EDATE, ETIME, DURATP * 10000 )
        SECSMIN = SECSDIFF( SDATE, STIME, EDATE, ETIME ) ! for point

        DO F = 1, NF2D
            EDATE = SDATEA( F )
            ETIME = STIMEA( F )
            CALL NEXTIME( EDATE, ETIME, DURATA( F ) * 10000 )
            SECS = SECSDIFF( SDATE, STIME, EDATE, ETIME )
            SECSMIN = MIN( SECSMIN, SECS )

        END DO
           
        TIMET = SEC2TIME( SECSMIN )

        NSTEPS= TIMET / 10000 ! number of time steps 

C.........  Build master output variables list
        NVOUT = 0

C.........  Loop through 2-D input files
        DO F = 1, NF2D

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
                END IF 

            END DO
        END DO

C.........  POINT
        IF( LFLAG ) THEN
            DO V = 1, NVARSP

                VNM = VNAMEP( V )

C.................  Look for variable name in output list
                K = INDEX1( VNM, NVOUT, OUTNAM  )  ! look in output list

C.................  If its not in the output list, add it
                IF( K .LE. 0 ) THEN
                    NVOUT = NVOUT + 1
                    INDXN ( NVOUT ) = NVOUT
                    OUTNAM( NVOUT ) = VNM
                    VDESCU( NVOUT ) = VDESCP( V )
                    VUNITU( NVOUT ) = VUNITP( V )
                END IF 
            END DO          
        ENDIF 

C.........  Sort output variables into alphabetical order
        CALL SORTIC( NVOUT, INDXN, OUTNAM )

C.........  Set up for opening output file
        IF ( LFLAG ) THEN
            IF( .NOT. DESC3( PNAME ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  PNAME( 1:LEN_TRIM( PNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

        ELSE
            IF( .NOT. DESC3( FNAME( 1 ) ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  FNAME( 1 )( 1:LEN_TRIM( FNAME(1) ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

        ENDIF

        SDATE3D = SDATE
        STIME3D = STIME
        NVARS3D = NVOUT

C........  Update variable names in sorted order, and also 
C........  set up logical arrays for which files have which species
        DO V = 1, NVOUT
            VNM = OUTNAM( INDXN( V ) )

            VNAME3D( V ) = VNM                  ! store sorted output vars, etc.
            VDESC3D( V ) = VDESCU( INDXN( V ) )
            UNITS3D( V ) = VUNITU( INDXN( V ) )
            VTYPE3D( V ) = M3REAL

            IF( LFLAG ) THEN
                LVOUTP( V ) = .FALSE.
                J = INDEX1( VNM, NVARSP, VNAMEP )
                IF( J .GT. 0 ) LVOUTP( V ) = .TRUE.
            ENDIF
    
            DO F = 1, NF2D
                LVOUTA( V,F ) = .FALSE.

                J = INDEX1( VNM, NVARSA( F ), VNAMEA( 1,F ) )
                IF( J .GT. 0 ) LVOUTA( V,F ) = .TRUE.
                                
            END DO 
        END DO

C.........  Allocate memory for the number of grid cells and layers
        NGRID = NROWS * NCOLS
        NLAYS = MAX( NLAYSP, 1 )
        ALLOCATE( E2D( NGRID, MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'E2D', PROGNAME )
        ALLOCATE( E3D( NGRID, NLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'E3D', PROGNAME )

C.........  Prompt for and open output file

        ONAME = PROMPTMFILE( 
     &          'Enter logical name for MERGED GRIDDED OUTPUT file',
     &          FSUNKN3, 'EMIS3D', PROGNAME )

C.........  Loop through hours

        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

C.............  Loop through species
            DO V = 1, NVOUT

                VNM = VNAME3D( V ) 

C.................  Initialize 2-d variables to zero, in case any are not read
                E2D = 0.

                DO F = 1, NF2D

C.................  If file has species, read (do this for all files)...

                    NAM = FNAME( F )
                    IF( LVOUTA( V,F ) ) THEN

                        IF( .NOT. READ3( NAM  , VNM  , 1 , JDATE,  
     &                                  JTIME, E2D( 1,F ) ) ) THEN

                            MESG = 'Could not read "' // VNM //
     &                             '" from file "' //
     &                             NAM( 1 : LEN_TRIM( NAM ) ) // '".'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ENDIF
                    ENDIF

                END DO

C.................  Read point sources - if not there, then initialize
                IF( LFLAG .AND. LVOUTP( V ) ) THEN

                    IF( .NOT. READ3( PNAME, VNM, ALLAYS3, 
     &                               JDATE, JTIME, E3D    ) ) THEN

                        MESG = 'Could not read "'// VNM //
     &                         '" from file "' 
     &                         //PNAME( 1 : LEN_TRIM( PNAME ) ) // '".'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

C.................  Initialize 
                ELSE

                    E3D = 0. ! array

                ENDIF

C.................  Add up emissions in layer 1 for hour/species

                DO F = 1, NF2D
                    DO C = 1, NGRID

                        E3D( C,1 ) = E3D( C,1 ) + E2D( C,F )

                    END DO
                END DO

C.................  Write species/hour to output file
            IF( .NOT. WRITE3( ONAME, VNM, JDATE, JTIME, E3D ) ) THEN

                MESG = 'Could not write "' // VNM // '" to file "' // 
     &                 ONAME( 1:LEN_TRIM( ONAME ) ) // '".'
     &                        
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            ENDIF

            END DO   

            CALL NEXTIME( JDATE, JTIME, TSTEP )
      
        END DO        

C......... Normal Completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)
	
C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )
 
C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT(  A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

	END
