
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

        EXTERNAL        PROMPTFFILE, PROMPTMFILE, ENVYN, GETFLINE, 
     &                  CVTRDTYPE, CVTVEHTYPE, CRLF

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: GRPLIST( :,: ) ! contents of GROUP file

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
        
        CHARACTER*16    ANAME   !  logical name for ASCII inventory file
        CHARACTER*16    ENAME   !  logical name for I/O API inventory file   
        CHARACTER*16    TNAME   !  logical name for I/O API temperature file
        CHARACTER*16    FNAME   !  logical name for I/O API emission factors file

C.........   Other local variables
        INTEGER    I, L              ! counters and indices
        INTEGER    IOS               ! i/o status
        INTEGER    NINVARR           ! number inventory variables to input
        INTEGER    SDATE_MET ! temperature file start date
        INTEGER    STIME_MET ! temperature file start time
        INTEGER    SDATE   !  episode start date
        INTEGER    STIME   !  episode start time
        INTEGER    TZONE   !  time zone (not used)
        INTEGER    TSTEP   !  time step of input temperature data (HHMMSS)
        INTEGER    NSTEPS  !  no. time steps in temperature data
        INTEGER    NROWS   !  no. grid rows   
        INTEGER    NGRPLINES  ! no. lines in GROUP file     
        INTEGER    NUMSCEN ! total number of M6 scenarios
        INTEGER    NUMSRC  ! total number of sources

        LOGICAL :: EFLAG    = .FALSE.      ! error flag
        LOGICAL :: TEMPFLAG = .TRUE.       ! true: replace temperatures in M6 scenarios 
        
        CHARACTER*20           MODELNAM  ! emission factor model name
        CHARACTER(LEN=80)      M6INPUT   ! Mobile6 input file name
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

C.........  Get the name of the emission factor model to use for one run
        MESG = 'Emission factor model'
        CALL ENVSTR( 'SMK_EF_MODEL', MESG, 'MOBILE6', MODELNAM, IOS )

C.........  End program if source category is not mobile sources
        IF( CATEGORY /= 'MOBILE' ) THEN
            L = LEN_TRIM( PROGNAME )
            MESG = 'Program ' // PROGNAME( 1:L ) // ' does not ' //
     &             'support ' // CATEGORY( 1:CATLEN ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Obtain settings from the environment...

C.........  Get environment variables that control program behavior
        TEMPFLAG = ENVYN( 'REPLACE_TEMPERATURES', 
     &                 'Replace temperatures in MOBILE6 scenarios',
     &                 .TRUE., IOS )
     
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

        GDEV = PROMPTFFILE(
     &           'Enter logical name for time period group file',
     &           .TRUE., .TRUE., 'GROUP', PROGNAME )
        
        IDEV = PROMPTFFILE(
     &           'Enter logical name for M6LIST scenarios file',
     &           .TRUE., .TRUE., 'M6LIST', PROGNAME )
        
        MDEV = PROMPTFFILE(
     &           'Enter logical name for MOBILE6 input file',
     &           .FALSE., .TRUE., 'M6INPUT', PROGNAME ) 
     
        TNAME = PROMPTMFILE(
     &          'Enter logical name for hourly temperature file',
     &          FSREAD3, 'HOURLYT', PROGNAME )
        
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
 
        END IF

C.........  Read header of temperature file
        IF ( .NOT. DESC3( TNAME ) ) THEN

            MESG = 'Could not get description of file ' // TNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Save header information that will be needed later
        ELSE
            SDATE_MET = SDATE3D
            STIME_MET = STIME3D
            TSTEP = TSTEP3D
            NSTEPS= MXREC3D
            NROWS = NROWS3D

        END IF

C.........  Set inventory variables to read
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'IVTYPE'
        NINVARR = 3

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Get default episode information for temperature data
        SDATE = SDATE_MET
        STIME = STIME_MET
        TZONE = 0
        CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

C.........  Set end date and time for run
C        EDATE = SDATE
C        ETIME = STIME
C        CALL NEXTIME( EDATE, ETIME, NSTEPS*10000 )

C.........  Read the hourly temperature file into an array
        MESG = 'Reading HOURLY temperature file...'
        CALL M3MSG2( MESG )

        ALLOCATE( TKCOUNTY( NROWS, 0:NSTEPS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKCOUNTY', PROGNAME )
        
        CALL RDHOURTEMP( TNAME, NROWS, NSTEPS, SDATE, STIME, TKCOUNTY )

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

C.........  Create the concatenated MOBILE6 input file
        MESG = 'Writing MOBILE6 input file...'
        CALL M3MSG2( MESG )
        
        NUMSCEN = 1
        NUMSRC  = 0
        
        CALL WRM6INPUT( GRPLIST, NGRPLINES, PDEV, MDEV, 
     &                  TKCOUNTY, NROWS, NSTEPS, NUMSCEN, NUMSRC )

        NUMSCEN = NUMSCEN - 1

C.........  Set up emission process variable names
        CALL EFSETUP( 'NONE', MODELNAM, MXVARS3, NEFS, VNAME3D, 
     &                 UNITS3D, VDESC3D )
     
C.........  Open file for storing emission factors (check this now rather than
C           waste time running Mobile6)
        FNAME = 'EMISFACS'
        CALL OPENSEF( NUMSRC, SDATE, STIME, FNAME )
        
C.........  Allocate space for storing emission factors
        ALLOCATE( EMISSIONS( NUMSCEN*161*24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISSIONS', PROGNAME )
        
        EMISSIONS = 0.

C.........  Get the full file name for the M6 input file and close the file
C           since MOBILE6 will reopen it
        INQUIRE( UNIT=MDEV, NAME=M6INPUT ) 
        CLOSE( MDEV )

C.........  Call Mobile6 with M6 input file name and unit number
C           Custom driver will use SMOKE log file for any screen output
        MESG = 'Running MOBILE6...' // CRLF()
        CALL M3MSG2( MESG )

        CALL SMKDRIVER( M6INPUT, MDEV, LDEV )

C.........  Match up emission factors with sources and write to file
        MESG = 'Writing emission factors to file...' // CRLF()
        CALL M3MSG2( MESG )
 
        CALL WREMFACS( FNAME, NUMSRC, SDATE )

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END PROGRAM EMISFAC