
        PROGRAM MBSETUP

C***********************************************************************
C  program body starts at line 115
C
C  DESCRIPTION:
C       Sets up necessary data for MOBILE6 runs. Reads in MCREF and 
C       MVREF files. Removes any counties outside of the grid. Creates
C       the SPDSUM file to assocate each source with a road type and speed.
C       Groups counties by requested temporal averaging types.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
        
C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: IFIP, SPEED, VMT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CATDESC, CRL, MAPNAM, 
     &                     NMAP, NSRC, NIACT

C.........  This module is used for MOBILE6 setup information        
        USE MODMBSET, ONLY: NREFC

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        LOGICAL         ENVYN

        EXTERNAL        CRLF, INDEX1, PROMPTFFILE, PROMPTMFILE, ENVYN

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
        INTEGER         IDEV     ! tmp unit number if ENAME is map file
        INTEGER         LDEV     ! unit number for log file
        INTEGER         PDEV     ! unit number for speeds summary file
        INTEGER         SDEV     ! unit number for ASCII inventory file
        INTEGER         VDEV     ! unit number for ref. county settings file
        INTEGER         XDEV     ! unit number for county cross-reference file
        INTEGER         ZDEV     ! unit number for speed cross-reference file
        
        CHARACTER*16    ANAME    ! logical name for ASCII inventory file
        CHARACTER*16    INAME    ! tmp name for inven file of unknown fmt
        CHARACTER*16    ENAME    ! logical name for I/O API inventory file
        CHARACTER*16    UNAME    ! logical name for ungridding-matrix input file

C.........   Other local variables
        INTEGER          I, K, L, M, S     ! counters and indices
        
        INTEGER          IOS               ! I/O status
        INTEGER          NINVARR           ! number inventory variables to input
        INTEGER          NMATX             ! size of ungridding matrix
        INTEGER          CURRCTY           ! current county FIPS code
        INTEGER          PREVCTY           ! previous county code
        INTEGER          NGRDCTY           ! no. counties inside grid
        
        LOGICAL       :: EFLAG   = .FALSE. !  error flag
        LOGICAL       :: SPDFLAG = .FALSE. !  true: use speed profiles
        
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
     &             'support ' // TRIM( CATEGORY ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Obtain settings from the environment...
C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Check if speed profiles are to be used
        SPDFLAG = ENVYN( 'USE_SPEED_PROFILES', 
     &            'Use speed profiles instead of inventory speeds', 
     &            .FALSE., IOS )
     
C.......   Get file names and units; open input files
C.........  Prompt for and open inventory file 
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

        UNAME = PROMPTMFILE(
     &          'Enter logical name for UNGRIDDING MATRIX file',
     &          FSREAD3, CRL // 'UMAT', PROGNAME )
     
        XDEV = PROMPTFFILE( 
     &           'Enter logical name for MCREF cross-reference file',
     &           .TRUE., .TRUE., 'MCREF', PROGNAME )
     
        VDEV = PROMPTFFILE(
     &           'Enter logical name for MVREF settings file',
     &           .TRUE., .TRUE., 'MVREF', PROGNAME )

        IF( SPDFLAG ) THEN
            ZDEV = PROMPTFFILE(
     &           'Enter logical name for SPDREF speed profile ' //
     &           'cross-reference file',
     &           .TRUE., .TRUE., 'SPDREF', PROGNAME );
        END IF
             
        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .FALSE., .TRUE., 'SPDSUM', PROGNAME )      

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API header info is passed by include file and the
C           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

C.........  Ensure that there is at least one activity in the inventory 
C           file, or else this program does not need to be run
        IF( NIACT == 0 ) THEN
            MESG = 'No activities are found in the ' //
     &             'inventory file!  Program cannot be used.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
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
        IVARNAMS( 2 ) = 'CSOURC'
        IVARNAMS( 3 ) = 'CSCC'
        IVARNAMS( 4 ) = 'IRCLAS'
        IVARNAMS( 5 ) = 'IVTYPE'
        NINVARR = 5

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Read speed and VMT information from inventory
        IF( .NOT. SPDFLAG ) THEN

C.............  Make sure speed is in the inventory
            M = INDEX1( 'SPEED', NMAP, MAPNAM )
            IF( M <= 0 ) THEN
                MESG = 'Mobile inventory does not include speed ' //
     &                 'data' // CRLF() // BLANK5 // 'Set the ' //
     &                 'USE_SPEED_PROFILES environment variable ' //
     &                 'to Y and try again.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            ALLOCATE( SPEED( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPEED', PROGNAME )
            CALL RDMAPPOL( NSRC, 1, 1, 'SPEED', SPEED )
        END IF
        
C.........  Make sure VMT is in the inventory
        M = INDEX1( 'VMT', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include VMT data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
                   
        ALLOCATE( VMT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VMT', PROGNAME )
        CALL RDMAPPOL( NSRC, 1, 1, 'VMT', VMT )

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

C.........  Create list of counties inside grid that have VMT values
        TMPCTY  = 0   ! array
        PREVCTY = 0
        CURRCTY = 0
        NGRDCTY = 0

        DO S = 1, NSRC
            IF( UMAT( S ) == 0 .OR. VMT( S ) == 0 ) CYCLE
            
            CURRCTY = IFIP( S )
            IF( CURRCTY /= PREVCTY ) THEN
                NGRDCTY = NGRDCTY + 1
                TMPCTY( NGRDCTY ) = CURRCTY
                PREVCTY = CURRCTY
            END IF            
        END DO

        ALLOCATE( GRIDCTY( NGRDCTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRIDCTY', PROGNAME )
        GRIDCTY = 0   ! array
        
        DO S = 1, NGRDCTY
            GRIDCTY( S ) = TMPCTY( S )
        END DO
        	
        DEALLOCATE( TMPCTY )

C.........  Read the county cross-reference file            
        MESG = 'Reading county cross-reference file...'
        CALL M3MSG2( MESG )
        
        CALL RDMCREF( XDEV, GRIDCTY, NGRDCTY )
        
C.........  Read the reference county settings file
        MESG = 'Reading reference county settings file...'
        CALL M3MSG2( MESG )
        
        CALL RDMVREF( VDEV )

C.........  Create speeds summary file
        MESG = 'Processing speed information...'
        CALL M3MSG2( MESG )

C.........  Process speed profiles
        IF( SPDFLAG ) THEN
        
C.............  Read speed profile cross-reference file
            CALL RDSPDREF( ZDEV )

C.............  Assign speed profile to each source
            CALL ASGNSPDS
            
        END IF

C.........  Loop through all the reference counties
        DO I = 1, NREFC
            CALL WRSPDSUM( PDEV, I, SPDFLAG )
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
