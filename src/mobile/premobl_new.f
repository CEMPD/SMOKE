
        PROGRAM PREMOBL

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET
        
        USE MODMBSET

        IMPLICIT NONE
        
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:        
        CHARACTER*2     CRLF
        INTEGER         GETIFDSC
        REAL            ENVREAL
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         WKDAY
        
        EXTERNAL     CRLF, GETIFDSC, ENVREAL, MMDDYY, PROMPTFFILE, 
     &               PROMPTMFILE, WKDAY
        
C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Gridded meteorology data (dim: NGRID)
        REAL   , ALLOCATABLE :: TA( : )   !  one layer of temperature

C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : ) ! Contiguous ungridding matrix

C...........  Allocatable per-source arrays
        INTEGER, ALLOCATABLE :: DAYBEGT ( : ) ! daily start time HHMMSS
        INTEGER, ALLOCATABLE :: DAYENDT ( : ) ! daily end time HHMMSS
        LOGICAL, ALLOCATABLE :: LDAYSAV ( : ) ! true: src uses daylight time

C...........  Array that contains the names of the inventory variables needed 
C             for this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   File units and logical names:
        INTEGER      LDEV  ! unit number for log file
        INTEGER      RDEV  ! unit number for mobile codes conversions file
        INTEGER      SDEV  ! unit number for ASCII inventory file
        INTEGER      TDEV  ! unit number for emissions type by activity file
        INTEGER      PDEV  ! unit number for speeds summary file (SPDSUM)

        CHARACTER*16 ANAME ! logical name for mobile ASCII inventory file
        CHARACTER*16 ENAME ! logical name for mobile I/O API inventory file
        CHARACTER*16 HNAME ! logical name for output ungridded hourly temps
        CHARACTER*16 TNAME ! logical name for surface temp input file
        CHARACTER*16 UNAME ! logical name for ungridding-matrix input file
                
C...........   Other local variables:
        INTEGER    H, I, J, K, L, S, T, V  ! Counters and pointers

        INTEGER    DAY     !  tmp day of week number
        INTEGER    EDATE   !  ending input date counter (YYYYDDD) in GMT
        INTEGER    ENLEN   !  length of the emissions inven name
        INTEGER    ETIME   !  ending input time counter (HHMMSS)  in GMT
        INTEGER    IDATE   !  output date for min/max
        INTEGER    IOS     !  temporary I/O status
        INTEGER    ITIME   !  output time for min/max
        INTEGER    JDATE   !  input date counter (YYYYDDD) in GMT
        INTEGER    JTIME   !  input time counter (HHMMSS)  in GMT
        INTEGER    LDATE   !  date from previous loop iteration
        INTEGER    NCOLS   !  no. grid columns
        INTEGER    NCOLSU  !  no. grid columns in ungridding matrix
        INTEGER    NCOUNTY !  no. counties to process
        INTEGER    NGRID   !  no. grid cells
        INTEGER    NINVARR !  no. inventory variables to read
        INTEGER    NMATX   !  size of ungridding matrix
        INTEGER    NROWS   !  no. grid rows
        INTEGER    NROWSU  !  no. grid rows in ungridding matrix
        INTEGER    NSTEPS  !  number of time steps to process temperature data
        INTEGER    ODATE   !  output date
        INTEGER    OTIME   !  time in GMT for determining when to output
        INTEGER :: OSRC = 0!  number of sources outside grid
        INTEGER    SDATE   !  output start date
        INTEGER    SDATE_MET ! met file start date
        INTEGER    STIME   !  output start time
        INTEGER    STIME_MET ! met file start time
        INTEGER    TSTEP   !  time step of input temperature data (HHMMSS)
        INTEGER    TZONE   !  zone to determine output days
        
        LOGICAL :: EFLAG    = .FALSE.  !  true: error found
        LOGICAL :: LASTTIME = .FALSE.  !  true: final time step
        LOGICAL :: OFLAG    = .FALSE.  !  true: ungridding is 0 for some srcs
                
        CHARACTER(LEN=IOVLEN3) :: TVARNAME    !  temperature variable name
        CHARACTER(LEN=300) ::     MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'PREMOBL'   !  program name

C***********************************************************************
C   begin body of program PREMOBL

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
        
C.........  Get the name of the activity to use for one run
        MESG = 'Temperature variable name'
        CALL ENVSTR( 'TVARNAME', MESG, 'TEMP1P5', TVARNAME, IOS )

C.........  Set default name of meterology file, depending on the name of the
C           temperature variable
        TNAME = 'MET_CRO_2D'
        IF( TVARNAME == 'TA' ) TNAME = 'MET_CRO_3D'

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.......   Get file names and units; open input files

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENLEN = LEN_TRIM( ENAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        UNAME = PROMPTMFILE(
     &          'Enter logical name for UNGRIDDING MATRIX file',
     &          FSREAD3, CRL // 'UMAT', PROGNAME )

        TNAME = PROMPTMFILE(
     &         'Enter logical name for SURFACE TEMPERATURE file',
     &          FSREAD3, TNAME, PROGNAME )

        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .TRUE., .TRUE., 'SPDSUM', PROGNAME )
     
C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:ENLEN ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API header info is passed by include file and the
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

C.........  Create note about time zone expected in meteorology file
        WRITE( MESG, 94010 )
     &     'NOTE: Time stamps of input meteorology file are assumed ' //
     &     CRLF() // BLANK5 // '      to be in GMT'
        CALL M3MSG2( MESG )

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
            NCOLS = NCOLS3D
            NGRID = NROWS * NCOLS

        END IF

C.........  Read header of ungridding matrix...
        IF( .NOT. DESC3( UNAME ) ) THEN
            MESG = 'Could not get description for file ' // UNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store number of ungridding factors
        NMATX = NCOLS3D

C.........  Check dimensions of ungridding matrix...
C.........  Check the number of sources 
        CALL CHKSRCNO( CATDESC, UNAME, NROWS3D, NSRC, EFLAG )

C.........  Compare the gridded file settings with the ungridded file settings
        NCOLSU = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
        NROWSU = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
        CALL CHECK_GRID_DIMS( 'NCOLS', 'columns', NCOLSU, NCOLS )
        CALL CHECK_GRID_DIMS( 'NROWS', 'rows'   , NROWSU, NROWS )

C......... If the dimensions were in error, abort
        IF( EFLAG ) THEN
            MESG = 'Ungridding matrix is inconsistent with inventory '//
     &             'and/or gridded ' // CRLF() // BLANK5 //
     &             'meteorology data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set inventory variables to read
        NINVARR = 6
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'CSOURC'
        IVARNAMS( 4 ) = 'CSCC'
        IVARNAMS( 5 ) = 'CLINK'
        IVARNAMS( 6 ) = 'TZONES'

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs and country/state/county codes
C           from the inventory arrays
        CALL GENUSLST
        
C.........  Get sources and counties from SPDSUM file
        ALLOCATE( COUNTYSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'COUNTYSRC', PROGNAME )

        COUNTYSRC = 0
        CALL RDSPDSRC( PDEV, NSRC, COUNTYSRC )

C.........  Count number of unique counties in COUNTYSRC array
        CALL GETCTYNUM( NSRC, COUNTYSRC, NCOUNTY )

C.........  Retrieve environment variable settings for temperature range
        MESG = 'Minimum allowed daily temperature [deg F]'
        MINTEMP = ENVREAL( 'SMK_MINTEMP', MESG, 0., IOS )

        MESG = 'Maximum allowed daily temperature [deg F]'
        MAXTEMP = ENVREAL( 'SMK_MAXTEMP', MESG, 120., IOS )

C.........  Convert temperature parameters to proper units
C.........  Kelvin for now - future can be dependant on met input units

        MINTEMP = FTOC * ( MINTEMP - 32. ) + CTOK
        MAXTEMP = FTOC * ( MAXTEMP - 32. ) + CTOK

C.........  Allocate memory for other arrays in the program
        ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UMAT', PROGNAME )
        ALLOCATE( TA( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( DAYBEGT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        ALLOCATE( DAYENDT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        ALLOCATE( LDAYSAV( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDAYSAV', PROGNAME )
        ALLOCATE( TASRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASRC', PROGNAME )

C.........  Create array of which sources are affected by daylight savings
        CALL GETDYSAV( NSRC, IFIP, LDAYSAV )

C.........  Read ungridding matrix 
        CALL RDUMAT( UNAME, NSRC, NMATX, NMATX, 
     &               UMAT( 1 ), UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ) )

C.........  Get default episode information for processing met data
        SDATE = SDATE_MET
        STIME = STIME_MET
        TZONE = 0
        CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

C.........  Set end date and time for run
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS*10000 )

C.........  Allocate memory for storing temperature profiles
        ALLOCATE( TKHOUR( NSRC, 0:NSTEPS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKHOUR', PROGNAME )
        ALLOCATE( TKCOUNTY( NCOUNTY, 0:NSTEPS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKCOUNTY', PROGNAME )
        
C.........  Preprocess dates, times, and time zones.  Write report for
C           time zones here for zones that do not have a complete day of
C           data in beginning or end. 
        CALL CHKFULLDY( NSRC, SDATE, STIME, EDATE, ETIME, 
     &                  TZONES, LDAYSAV )

C.........  Open file for per-county hourly temperature profiles       
        HNAME = 'HOURLYT'
        CALL OPENSHOUR( ENAME, SDATE, STIME, TVARNAME, NCOUNTY, HNAME )

C.........  Process temperature information...

        L = LEN_TRIM( TVARNAME )
        MESG = 'Processing temperature data using variable "' //
     &         TVARNAME( 1:L ) // '" ...'
        CALL M3MSG2( MESG )

C.........  Loop through days/hours of temperature files
        IDATE = SDATE
        ITIME = 0
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        DO T = 1, NSTEPS

C.................  Read current temperature file
            IF ( .NOT. READ3( TNAME, TVARNAME, 1, 
     &                        JDATE, JTIME, TA ) ) THEN
                L = LEN_TRIM( TVARNAME )
                MESG = 'Could not read ' // TVARNAME( 1:L ) //
     &                 ' from ' // TNAME 
               CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

C.............  Apply ungridding matrix 
            CALL APPLUMAT( NSRC, NMATX, TA, UMAT(1), UMAT( NSRC+1 ),
     &                     UMAT( NSRC+NMATX+1 ), TASRC )

C.............  When new day...
            IF ( JDATE /= LDATE ) THEN

C.................  Write message for day of week and date
                DAY = WKDAY( JDATE )
                MESG = 'Processing ' // DAYS( DAY ) // ' ' // 
     &                 MMDDYY( JDATE )
                CALL M3MSG2( MESG )

C.................  Set start and end hours of day for all sources
                CALL SETSRCDY( NSRC, JDATE, TZONES, LDAYSAV, 
     &                         DAYBEGT, DAYENDT )

            END IF

C.............  First iteration in loop, set the output date/time
            IF( T == 1 ) THEN
                CALL SETOUTDT( NSRC, JDATE, TZONES, DAYBEGT, DAYENDT, 
     &                         ODATE, OTIME )
            END IF

C.............  Create hourly temperature array by source
            CALL HOURTEMP( NSRC, NSTEPS, T, JTIME, DAYBEGT, 
     &                     TASRC, TKHOUR )

            LASTTIME = ( T == NSTEPS )    

C.............  Adjust and output min/max data
            IF( LASTTIME .OR.
     &        ( JDATE >= ODATE .AND.
     &          JTIME == OTIME       ) ) THEN

C.................  Adjust hourly temperature data based on min and max
                CALL ADJSHOUR( NSRC, NSTEPS, MINTEMP, MAXTEMP, 
     &                         'temperature', TKHOUR )

C.................  Create county based 24-hour temperature profiles
                MESG = 'Averaging 24-hour temperature profiles...'
                CALL M3MSG2( MESG )

                CALL AVERTEMP( NSRC, NSTEPS, NCOUNTY, COUNTYSRC, 
     &                         TKHOUR, TKCOUNTY )

C.................  Count these sources (the first time)
                IF( .NOT. OFLAG ) THEN
                    DO S = 1, NSRC

                	IF( UMAT( S ) == 0 ) THEN
                            OFLAG = .TRUE.
                            OSRC = OSRC + 1
                	END IF

                    END DO
                END IF

C.................  Write hourly temperature data to file                
                MESG = 'Writing temperature profiles to output file...'
                CALL M3MSG2( MESG )
                
                CALL WRSHOUR( NCOUNTY, NSTEPS, IDATE, JTIME, 'HOURLYT',
     &                        TKCOUNTY )

C.................  Increment output time for per-day file
C                CALL NEXTIME( IDATE, ITIME, 240000 )

            END IF   ! End of section to output and store for emission factors

            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO   !  End loop on hours of temperature files
 
C.........  Write message when sources were excluded during ungridding
        IF( OFLAG ) THEN

            WRITE( MESG, 94010 )
     &             'NOTE: During ungridding, ', OSRC,
     &             'sources excluded from grid.' // CRLF() // BLANK10 //
     &             'These sources may be outside the grid.'

            CALL M3MESG( MESG )

         END IF
                    
C......... End program sucessfully

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
 
C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Formatted file I/O formats............ 93xxx
 
93010   FORMAT( I8, 1X, F13.5, 1X, F13.5, 1X, I8, 1X, A )
 
C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram checks the dimension of the grid

            SUBROUTINE CHECK_GRID_DIMS( VNAME, VDESC, UVAL, TVAL )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: VNAME   ! name of value
            CHARACTER(*), INTENT (IN) :: VDESC   ! description of value
            INTEGER     , INTENT (IN) :: UVAL    ! ungridding matrix value
            INTEGER     , INTENT (IN) :: TVAL    ! met data file value

C.............  Local variables
            INTEGER       L1, L2

C----------------------------------------------------------------------

            IF( UVAL /= TVAL ) THEN
                EFLAG = .TRUE.
                L1 = LEN_TRIM( VNAME )
                L2 = LEN_TRIM( VDESC )
                WRITE( MESG,94010 ) 'ERROR: Inconsistent number of ' //
     &                 VDESC( 1:L2 ) // '. In ungridding matrix, ' // 
     &                 VNAME( 1:L1 ) // '= ', UVAL,
     &                 CRLF() // BLANK5 // 
     &                 'but in meteorology file, ' // VNAME( 1:L1 ) // 
     &                 '= ', TVAL
                CALL M3MSG2( MESG )

            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_GRID_DIMS
                
        END PROGRAM PREMOBL