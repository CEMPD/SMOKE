
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

C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : ) ! Contiguous ungridding matrix

C.........  Arrays to hold counties inside grid
        INTEGER, ALLOCATABLE :: TMPCTY ( : )  ! temporary holding array
        INTEGER, ALLOCATABLE :: GRIDCTY( : )  ! counties inside grid 
        
C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )
        
C.........  Unit numbers and logical file names
        INTEGER         LDEV     ! unit number for log file
        INTEGER         XDEV     ! unit number for county cross-reference file
        INTEGER         VDEV     ! unit number for ref. county settings file
        INTEGER         SDEV     ! unit number for ASCII inventory file
        INTEGER         PDEV     ! unit number for speeds summary file
        
        CHARACTER*16    ANAME   !  logical name for ASCII inventory file
        CHARACTER*16    ENAME   !  logical name for I/O API inventory file
        CHARACTER*16    UNAME   ! logical name for ungridding-matrix input file

C.........   Other local variables
        INTEGER          I, K, L, S        ! counters and indices
        
        INTEGER          IOS               ! I/O status
        INTEGER          NINVARR           ! number inventory variables to input
        INTEGER          NMATX             ! size of ungridding matrix
        INTEGER          CURRCTY           ! current county FIPS code
        INTEGER          PREVCTY           ! previous county code
        INTEGER          NUMCTY            ! no. counties inside grid
        
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

        UNAME = PROMPTMFILE(
     &          'Enter logical name for UNGRIDDING MATRIX file',
     &          FSREAD3, CRL // 'UMAT', PROGNAME )
     
        XDEV = PROMPTFFILE( 
     &           'Enter logical name for MCREF cross-reference file',
     &           .TRUE., .TRUE., 'MCREF', PROGNAME )
     
        VDEV = PROMPTFFILE(
     &           'Enter logical name for MVREF settings file',
     &           .TRUE., .TRUE., 'MVREF', PROGNAME )
     
        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .FALSE., .TRUE., 'SPDSUM', PROGNAME )      

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

C.........  Read header of ungridding matrix...
        IF( .NOT. DESC3( UNAME ) ) THEN
            MESG = 'Could not get description for file ' // UNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store number of ungridding factors
        NMATX = NCOLS3D

C.........  Check the number of sources in ungridding matrix against inventory
        CALL CHKSRCNO( CATDESC, UNAME, NROWS3D, NSRC, EFLAG )

C......... If the dimensions were in error, abort
        IF( EFLAG ) THEN
            MESG = 'Ungridding matrix is inconsistent with inventory.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
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

C.........  Allocate memory for ungridding matrix
        ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UMAT', PROGNAME )
        ALLOCATE( TMPCTY( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPCTY', PROGNAME )

C.........  Read ungridding matrix 
        CALL RDUMAT( UNAME, NSRC, NMATX, NMATX, 
     &               UMAT( 1 ), UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ) )

C.........  Create list of counties inside grid
        TMPCTY  = 0   ! array
        PREVCTY = 0
        CURRCTY = 0
        NUMCTY  = 0

        DO S = 1, NSRC
            IF( UMAT( S ) == 0 ) CYCLE
            
            CURRCTY = IFIP( S )
            IF( CURRCTY /= PREVCTY ) THEN
                NUMCTY = NUMCTY + 1
                TMPCTY( NUMCTY ) = CURRCTY
                PREVCTY = CURRCTY
            END IF            
        END DO

        ALLOCATE( GRIDCTY( NUMCTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRIDCTY', PROGNAME )
        GRIDCTY = 0   ! array
        
        DO S = 1, NUMCTY
            GRIDCTY( S ) = TMPCTY( S )
        END DO
        	
        DEALLOCATE( TMPCTY )

C.........  Read the county cross-reference file            
        MESG = 'Reading county cross-reference file...'
        CALL M3MSG2( MESG )
        
        CALL RDMCREF( XDEV, GRIDCTY, NUMCTY )
        
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

        WRITE( PDEV,94010 ) ' '

C.........  Create time period group files
        MESG = 'Writing time period group files...'
        CALL M3MSG2( MESG )
        
        CALL WRTIMEGR

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END PROGRAM MBSETUP