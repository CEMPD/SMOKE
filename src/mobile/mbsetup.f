
        PROGRAM MBSETUP
        
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

        EXTERNAL        PROMPTFFILE, 
     &                  PROMPTMFILE

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C...........   LOCAL VARIABLES and their descriptions:

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )
        
C.........  Unit numbers and logical file names
        INTEGER         LDEV     ! unit number for log file
        INTEGER         XDEV     ! unit number for county cross-reference file
        INTEGER         VDEV     ! unit number for ref. county settings file
        INTEGER         SDEV     ! unit number for ASCII inventory file
        INTEGER         PDEV     ! unit number for speeds summary file
        
        INTEGER         DAYDEV   ! unit number for daily time period group file 
        INTEGER         WEEKDEV  ! unit number for weekly group file
        INTEGER         MONTHDEV ! unit number for monthly group file
        INTEGER         METDEV   ! unit number for met file length time period group
        
        CHARACTER*16    ANAME   !  logical name for ASCII inventory file
        CHARACTER*16    ENAME   !  logical name for I/O API inventory file

C.........   Other local variables
        INTEGER          I, L              ! counters and indices
        INTEGER          NINVARR           ! number inventory variables to input
        
        LOGICAL       :: EFLAG   = .FALSE. !  error flag
        
        CHARACTER*300          MESG      !  message buffer 
        
        CHARACTER*16  :: PROGNAME = 'MBSETUP' ! program name
        
C***********************************************************************
C   begin body of program MBSETUP

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
C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.......   Get file names and units; open input files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for MCREF cross-reference file',
     &           .TRUE., .TRUE., 'MCREF', PROGNAME )
     
        VDEV = PROMPTFFILE(
     &           'Enter logical name for MVREF settings file',
     &           .TRUE., .TRUE., 'MVREF', PROGNAME )
     
        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .FALSE., .TRUE., 'SPDSUM', PROGNAME )

C.........  Get file units for time period group files
        DAYDEV = PROMPTFFILE(
     &           'Enter logical name for DAILYGROUP file',
     &           .FALSE., .TRUE., 'DAILYGROUP', PROGNAME )
     
        WEEKDEV = PROMPTFFILE(
     &           'Enter logical name for WEEKLYGROUP file',
     &           .FALSE., .TRUE., 'WEEKLYGROUP', PROGNAME )
     
        MONTHDEV = PROMPTFFILE(
     &           'Enter logical name for MONTHLYGROUP file',
     &           .FALSE., .TRUE., 'MONTHLYGROUP', PROGNAME )
     
        METDEV = PROMPTFFILE(
     &           'Enter logical name for METGROUP file',
     &           .FALSE., .TRUE., 'METGROUP', PROGNAME )           

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
        IVARNAMS( 3 ) = 'SPEED'
        NINVARR = 3

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs and country/state/county codes
C           from the inventory arrays
        CALL GENUSLST
     
C.........  Read the county cross-reference file            
        MESG = 'Reading county cross-reference file...'
        CALL M3MSG2( MESG )
        
        CALL RDMCREF( XDEV )
        
C.........  Read the reference county settings file
        MESG = 'Reading reference county settings file...'
        CALL M3MSG2( MESG )
        
        CALL RDMVREF( VDEV )

C.........  Create speeds summary file
        MESG = 'Processing speed information...'
        CALL M3MSG2( MESG )

C.........  Loop through all the reference counties
        DO I = 1, NREFC
            CALL WRSPDSUM( PDEV, I )
        END DO

C THIS MAY BE MOVED TO PREMOBL LATER
C.........  Create time period group files
        MESG = 'Writing time period group files...'
        CALL M3MSG2( MESG )
        
        CALL WRTIMEGR( DAYDEV,   DAILY )
        CALL WRTIMEGR( WEEKDEV,  WEEKLY )
        CALL WRTIMEGR( MONTHDEV, MONTHLY )
        CALL WRTIMEGR( METDEV,   METLEN )

C END PREMOBL SECTION

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END PROGRAM MBSETUP