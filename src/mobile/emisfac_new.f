
        PROGRAM EMISFAC
        
C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO
        
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

        EXTERNAL        PROMPTFFILE, 
     &                  PROMPTMFILE, ENVYN

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C...........   LOCAL VARIABLES and their descriptions:

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

C.........   Other local variables
        INTEGER          I, L              ! counters and indices        
        INTEGER          IOS               ! i/o status
        INTEGER          NINVARR           ! number inventory variables to input
        
        LOGICAL :: EFLAG    = .FALSE.      ! error flag
        LOGICAL :: TEMPFLAG = .TRUE.       ! true: replace temperatures in M6 scenarios 
        
        CHARACTER*300          MESG      !  message buffer 
        
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

C.........  Set inventory variables to read
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'IVTYPE'
        NINVARR = 3

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Read the M6LIST file into an array
        MESG = 'Reading M6LIST file...'
        CALL M3MSG2( MESG )

        CALL RDM6LIST( IDEV )

C.........  Allocate memory for the source/scenario number array
        ALLOCATE( SCENLIST( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCENLIST', PROGNAME )

C.........  Create the concatenated MOBILE6 input file
        MESG = 'Writing MOBILE6 input file...'
        CALL M3MSG2( MESG )
        
        CALL WRM6INPUT( GDEV, PDEV, MDEV )

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END PROGRAM EMISFAC