
        PROGRAM EMISFAC
        
C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: NEFS
                        
        USE MODMBSET

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         CVTRDTYPE
        INTEGER         CVTVEHTYPE
        CHARACTER*2     CRLF
        LOGICAL         OPNFULL3
        INTEGER         SECSDIFF
        LOGICAL         ISOPEN
        INTEGER         STR2INT
        LOGICAL         CHKINT
        INTEGER         JUNIT

        EXTERNAL        PROMPTFFILE, PROMPTMFILE, ENVYN, GETFLINE, 
     &                  CVTRDTYPE, CVTVEHTYPE, CRLF, OPNFULL3, 
     &                  SECSDIFF, ISOPEN, STR2INT, CHKINT, JUNIT

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: GRPLIST( :,: ) ! contents of GROUP file
        INTEGER, ALLOCATABLE :: TEMPCTY( : )   ! list of counties from temperature file

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )
        
C.........  Unit numbers and logical file names
        INTEGER         LDEV     ! unit number for log file
        INTEGER         SDEV     ! unit number for ASCII inventory file
        INTEGER         PDEV     ! unit number for speeds summary file (SPDSUM)
        INTEGER         GDEV     ! unit number for time period group file (GROUP)
        INTEGER         IDEV     ! unit number for county MOBILE6 scenarios file (M6LIST)
        INTEGER         MDEV     ! unit number for concatenated MOBILE6 input file (M6INPUT)
        INTEGER         TDEV     ! unit number for emission processes file
        
        CHARACTER*16    ANAME   !  logical name for ASCII inventory file
        CHARACTER*16    ENAME   !  logical name for I/O API inventory file   
        CHARACTER*16    TNAME   !  logical name for I/O API temperature file
        CHARACTER*16    FNAME   !  logical name for I/O API emission factors file

C.........   Other local variables
        INTEGER    I, L          ! counters and indices
        INTEGER    IOS           ! i/o status
        INTEGER    NINVARR       ! number inventory variables to input
        INTEGER    SDATE         ! current start date
        INTEGER    STIME         ! current start time
        INTEGER    EDATE         ! current end date
        INTEGER    TZONE         ! time zone (not used)
        INTEGER    TSTEP         ! time step of input temperature data (HHMMSS)
        INTEGER    TIMEDIFF      ! difference between file date and current sdate
        INTEGER    TEMPDATE      ! temporary date to pass to read routine
        INTEGER    TEMPTIME      ! temporary time to pass to read routine
        INTEGER    NSTEPS        ! no. time steps in temperature data
        INTEGER    NROWS         ! no. grid rows   
        INTEGER    NGRPLINES     ! no. lines in GROUP file     
        INTEGER    NUMSCEN       ! total number of M6 scenarios
        INTEGER    NUMSRC        ! total number of sources

        LOGICAL :: EFLAG    = .FALSE.    ! error flag
        LOGICAL :: TEMPFLAG = .TRUE.     ! true: replace temperatures in M6 scenarios 
        LOGICAL :: INITIAL  = .TRUE.     ! true: first time through loop
        LOGICAL :: FEXIST   = .FALSE.    ! true: file exists
        
        CHARACTER*20           MODELNAM  ! emission factor model name
        CHARACTER*20           GRP_NAME  ! temperature aggregation group
        CHARACTER(LEN=IOVLEN3) VOLNAM    ! volatile pollutant name
        CHARACTER(LEN=80)      M6INPUT   ! Mobile6 input file name
        CHARACTER(LEN=200)     TEMPDIR   ! location of hourly temperature files
        CHARACTER(LEN=256)     TEMPNAME  ! full temperature file name
        CHARACTER(LEN=50)      M6DIR     ! location of MOBILE6 files
        CHARACTER(LEN=200)     EMISDIR   ! directory for output EF files
        CHARACTER(LEN=20)      SEARCHSTR ! string used in search
        CHARACTER(LEN=MXDLEN3) TEMPLINE  ! line from file description
        CHARACTER*300          MESG      ! message buffer 

        CHARACTER*16  :: PROGNAME = 'EMISFAC' ! program name
        
C***********************************************************************
C   begin body of program EMISFAC

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  End program if source category is not mobile sources
        IF( CATEGORY /= 'MOBILE' ) THEN
            L = LEN_TRIM( PROGNAME )
            MESG = 'Program ' // PROGNAME( 1:L ) // ' does not ' //
     &             'support ' // CATEGORY( 1:CATLEN ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Obtain settings from the environment...

C.........  Get the name of the emission factor model to use
        MESG = 'Emission factor model'
        CALL ENVSTR( 'SMK_EF_MODEL', MESG, 'MOBILE6', MODELNAM, IOS )

C.........  Get environment variables that control program behavior
        TEMPFLAG = ENVYN( 'REPLACE_TEMPERATURES', 
     &                 'Replace temperatures in MOBILE6 scenarios',
     &                 .TRUE., IOS )
     
C.........  Get temperature aggregation length
        MESG = 'Temperature aggregation group'
        CALL ENVSTR( 'GROUP_TYPE', MESG, 'daily', GRP_NAME, IOS )
     
C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.......   Get file names and units; open input files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )
     
        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .TRUE., .TRUE., 'SPDSUM', PROGNAME )
        
        IDEV = PROMPTFFILE(
     &           'Enter logical name for M6LIST scenarios file',
     &           .TRUE., .TRUE., 'M6LIST', PROGNAME )
        
        TDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION PROCESSES file',
     &           .TRUE., .TRUE., 'MEPROC', PROGNAME )

        SELECT CASE( GRP_NAME )
        CASE( 'daily' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for DAILYGROUP file',
     &           .TRUE., .TRUE., 'DAILYGROUP', PROGNAME )
        CASE( 'weekly' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for WEEKLYGROUP file',
     &           .TRUE., .TRUE., 'WEEKLYGROUP', PROGNAME )
        CASE( 'monthly' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for MONTHLYGROUP file',
     &           .TRUE., .TRUE., 'MONTHLYGROUP', PROGNAME )
        CASE( 'episode' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for EPISODEGROUP file',
     &           .TRUE., .TRUE., 'EPISODEGROUP', PROGNAME )
        CASE DEFAULT
            MESG = 'ERROR: Unrecognized value for GROUP_TYPE.' //
     &             CRLF() // BLANK10 // 'Valid values are daily, ' //
     &             'weekly, monthly, or episode.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END SELECT
     
C.........  Get temperature file directory from the environment
        MESG = 'Location of hourly temperature files'
        CALL ENVSTR( 'SMK_TEMPATH', MESG, '.', TEMPDIR, IOS )
     
        TNAME = PROMPTMFILE(
     &          'Enter logical name for first hourly temperature file',
     &          FSREAD3, 'HOURLYT', PROGNAME )
        
C.........  Get MOBILE6 directory from the environment
        MESG = 'Location of MOBILE6 input and output files ' //
     &         '(50 characters or less)'
        CALL ENVSTR( 'SMK_M6PATH', MESG, '.', M6DIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: MOBILE6 files being placed in ' //
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_M6PATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF        

C.........  Get header description of inventory file 
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:LEN_TRIM( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

C.............  Ensure that there is at least one activity in the inventory 
C               file, or else this program does not need to be run
            IF( NIACT == 0 ) THEN
                MESG = 'ERROR: No activities are found in the ' //
     &                 'inventory file!  Program cannot be used.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

C.........  Set inventory variables to read
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'IVTYPE'
        NINVARR = 3

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Set up emission process variable names
        CALL EFSETUP( 'NONE', MODELNAM, MXVARS3, NEFS, VNAME3D, 
     &                 UNITS3D, VDESC3D, VOLNAM )

C.........  Read emission processes file.  Populate array in MODEMFAC.
        CALL RDEPROC( TDEV )

C.........  Get output directory information from the environment
        MESG = 'Path where emission factors files will be written'
        CALL ENVSTR( 'SMK_EMISPATH', MESG, '.', EMISDIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: Emission factors files being placed in ' //
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_EMISPATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF
        
C.........  Read the GROUP list file into an array
        MESG = 'Reading GROUP list file...'
        CALL M3MSG2( MESG )

C.........  Get the number of lines in the GROUP file     
        NGRPLINES = GETFLINE( GDEV, 'GROUP county list file' )

C.........  If no counties in group file, exit program        
        IF( NGRPLINES == 0 ) THEN
            MESG = 'ERROR: No counties listed in GROUP file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 
        
        ALLOCATE( GRPLIST( NGRPLINES,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPLIST', PROGNAME )

C.........  Read GROUP file
        CALL RDGRPLIST( GDEV, NGRPLINES, GRPLIST )
        
C.........  Read the M6LIST file into an array
        MESG = 'Reading M6LIST file...'
        CALL M3MSG2( MESG )

        CALL RDM6LIST( IDEV )

C.........  Allocate memory for the source/scenario number array
        ALLOCATE( SCENLIST( NSRC,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCENLIST', PROGNAME )
        
        SCENLIST = 0

C.........  Start loop over temperature files
        DO

            IF( TEMPFLAG ) THEN
C.................  Read header of temperature file
                IF ( .NOT. DESC3( TNAME ) ) THEN
                
                    MESG = 'Could not get description of file ' // TNAME
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                
C.................  Save header information that will be needed later
                ELSE
                    SDATE  = SDATE3D
                    STIME  = STIME3D
                    TSTEP  = TSTEP3D
                    NSTEPS = MXREC3D
                    NROWS  = NROWS3D
                    
C.....................  Find end date in file description
                    SEARCHSTR = '/END DATE/ '
                    L = LEN_TRIM( SEARCHSTR ) + 1
                    
                    DO I = 1, MXDESC3
                       IF( INDEX( FDESC3D( I ), 
     &                            SEARCHSTR( 1:L ) ) > 0 ) THEN
     	                       TEMPLINE = FDESC3D( I )
                       	   IF( CHKINT( TEMPLINE( L+1:L+8 ) ) ) THEN
                               EDATE = STR2INT( TEMPLINE( L+1:L+8 ) )
                               EXIT
                           ELSE
                               EFLAG = .TRUE.
                               EXIT
                           END IF
                       END IF
                       
                       IF( I == MXDESC3 ) THEN
                           EFLAG = .TRUE.
                       END IF
                    END DO
                    	
                    IF( EFLAG ) THEN
                        MESG = 'Could not get ending date of file ' //
     &                         TNAME
                    	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF

C.................  Allocate space for temperature info            
                IF( INITIAL ) THEN
                	ALLOCATE( TEMPCTY( NROWS ), STAT=IOS )
                	CALL CHECKMEM( IOS, 'TEMPCTY', PROGNAME )
                    ALLOCATE( TKHOUR( NROWS, 24 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'TKHOUR', PROGNAME )
                END IF
                
                TEMPCTY = 0
                TKHOUR = 0.

C.................  Read contents of hourly temperature file
                MESG = 'Reading HOURLY temperature file...'
                CALL M3MSG2( MESG )
            	
C.................  Read county list from file
                IF( .NOT. READ3( TNAME, 'COUNTIES', 1, SDATE, STIME, 
     &                           TEMPCTY ) ) THEN
     	            MESG = 'Could not read COUNTIES from ' // TNAME
     	            CALL M3EXIT( PROGNAME, SDATE, STIME, MESG, 2 )
     	        END IF
              
C.................  Read temperature data from file            
                TEMPDATE = SDATE
                TEMPTIME = STIME
                CALL RDHOURTEMP( TNAME, NROWS, TEMPDATE, TEMPTIME, 
     &                           TKHOUR )
            
            END IF
            
C.............  Create the concatenated MOBILE6 input file
            MESG = 'Writing MOBILE6 input file...'
            CALL M3MSG2( MESG )

C.............  Open new input file
            MDEV = JUNIT()
            WRITE( M6INPUT,94010 ) M6DIR( 1:LEN_TRIM( M6DIR ) ) // 
     &                             '/m6input.' // 
     &                             GRP_NAME( 1:LEN_TRIM( GRP_NAME ) ) //
     &                             '.', SDATE, '.txt'
            OPEN( UNIT=MDEV, FILE=M6INPUT, STATUS='REPLACE', 
     &            ACTION='WRITE', IOSTAT=IOS )
            
            IF( IOS /= 0 ) THEN
                MESG = 'Could not create file ' // M6INPUT
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF     
            
            NUMSCEN = 1
            NUMSRC  = 0
            CALL WRM6INPUT( GRPLIST, NGRPLINES, PDEV, MDEV, 
     &                      TEMPCTY, TKHOUR, NROWS, VOLNAM, 
     &                      NUMSCEN, NUMSRC, TEMPFLAG )
            
            CLOSE( MDEV )
            
            NUMSCEN = NUMSCEN - 1
            
C.............  Open file for storing emission factors (check this now rather than
C               waste time running Mobile6)
            FNAME = 'EMISFACS'
            
            IF( ISOPEN( FNAME ) ) THEN
                IF( .NOT. CLOSE3( FNAME ) ) THEN
                    MESG = 'Could not close file ' // 
     &                      FNAME( 1:LEN_TRIM( FNAME ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF
            
            CALL OPENSEF( NUMSRC, GRP_NAME, SDATE, EDATE, 
     &                    EMISDIR, FNAME )
            
C.............  Allocate space for storing emission factors
            IF( INITIAL ) THEN
                ALLOCATE( EMISSIONS( NUMSCEN*NM6PROFS*24 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EMISSIONS', PROGNAME )
                INITIAL = .FALSE.
            END IF
            
            EMISSIONS = 0.
            
C.............  Call Mobile6 with M6 input file name
C               Custom driver will use SMOKE log file for any screen output
            MESG = 'Running MOBILE6...' // CRLF()
            CALL M3MSG2( MESG )
            
            CALL SMKDRIVER( M6INPUT, LDEV, EFLAG )
            
            IF( EFLAG ) THEN
                MESG = 'Problem running MOBILE6'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
C.............  Match up emission factors with sources and write to file
            MESG = 'Writing emission factors to file...' // CRLF()
            CALL M3MSG2( MESG )
            
            CALL WREMFACS( FNAME, NUMSRC, SDATE )

            IF( TEMPFLAG ) THEN

C.................  Close current temperature file
                IF( .NOT. CLOSE3( TNAME ) ) THEN
                    MESG = 'Could not close file ' // 
     &                     TNAME( 1:LEN_TRIM( TNAME ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
C.................  Open next temperature file if available
             
C.................  Set start date of next file based on end date of current file
                SDATE = EDATE
                CALL NEXTIME( SDATE, STIME, 24*TSTEP )
             
C.................  Construct next file name
                WRITE( TEMPNAME,94010 ) 
     &                 TEMPDIR( 1:LEN_TRIM( TEMPDIR ) ) //
     &                 '/' // GRP_NAME ( 1:LEN_TRIM( GRP_NAME ) ) //
     &                 '.', SDATE, '.ncf'
             
C.................  Check if file exists
                INQUIRE( FILE=TEMPNAME, EXIST=FEXIST )
             
C.................  If file does not exist, we're done            
                IF( .NOT. FEXIST ) EXIT
                
C.................  Open temperature file
                IF( .NOT. OPNFULL3( TNAME, FSREAD3, 
     &                              TEMPNAME, PROGNAME ) ) THEN
                    MESG = 'Could not open temperature file ' // 
     &                     TEMPNAME( 1:LEN_TRIM( TEMPNAME ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

        END DO

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I7, A )

94020   FORMAT( A, I7, A1, I7, A )

        END PROGRAM EMISFAC
