
        PROGRAM PREMOBL

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID
        
        USE MODMBSET, ONLY: COUNTYSRC, DAILY, WEEKLY, MONTHLY, EPISLEN

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
        INTEGER         GETFLINE
        INTEGER         ENVINT
        REAL            ENVREAL
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         WKDAY
        LOGICAL         OPNFULL3
        LOGICAL         ISOPEN
        
        EXTERNAL     CRLF, GETIFDSC, GETFLINE, ENVINT, 
     &               ENVREAL, MMDDYY, PROMPTFFILE, PROMPTMFILE, 
     &               SECSDIFF, WKDAY, OPNFULL3, ISOPEN
        
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

C...........  Allocatable arrays for met data
        INTEGER, ALLOCATABLE :: METCHECK( : ) ! dimension: nsteps in episode, 
                                              ! value indicates which met data file covers that hour
        CHARACTER(LEN=256), ALLOCATABLE :: METLIST( : ) ! listing of met file names
        CHARACTER(LEN=16) , ALLOCATABLE :: METLOGS( : ) ! listing of met logical names

        INTEGER, ALLOCATABLE :: NDAYSRC( :,: )  ! number of days to average by source

C...........  Array that contains the names of the inventory variables needed 
C             for this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   File units and logical names:
        INTEGER      LDEV  ! unit number for log file
        INTEGER      SDEV  ! unit number for ASCII inventory file
        INTEGER      TDEV  ! unit number for meteorology list file
        INTEGER      PDEV  ! unit number for speeds summary file (SPDSUM)
        INTEGER      DDEV  ! unit number for daily group file
        INTEGER      WDEV  ! unit number for weekly group file
        INTEGER      MDEV  ! unit number for monthly group file
        INTEGER      EDEV  ! unit number for episode group file

        CHARACTER*16 ANAME ! logical name for mobile ASCII inventory file
        CHARACTER*16 ENAME ! logical name for mobile I/O API inventory file
        CHARACTER*16 DNAME ! logical name for daily output ungridded hourly temps
        CHARACTER*16 WNAME ! logical name for weekly output hourly temps
        CHARACTER*16 MNAME ! logical name for monthly output hourly temps
        CHARACTER*16 PNAME ! logical name for episode output hourly temps
        CHARACTER*16 UNAME ! logical name for ungridding-matrix input file
                
C...........   Other local variables:
        INTEGER    J, K, L, N, S, T, T2, V  ! Counters and pointers

        INTEGER    EPI_SDATE      ! episode start date from E.V. (YYYYDDD)
        INTEGER    EPI_STIME      ! episode start time from E.V. (HHMMSS)
        INTEGER    EPI_RUNLEN     ! episode duration from E.V. (HHMMSS)
        INTEGER    EPI_NSTEPS     ! episode number of time steps
        INTEGER    EPI_EDATE      ! episode ending date based on ERUNLEN
        INTEGER    EPI_ETIME      ! episode ending time based on ERUNLEN
        
        INTEGER    ARRAYPOS    ! position in 24-hour arrays
        INTEGER    CURRMNTH    ! current month
        INTEGER    CURRDAY     ! current day
        INTEGER    DAY         ! tmp day of week number
        INTEGER    DDATE       ! output date for daily counties
        INTEGER    DTIME       ! output time in local time
        INTEGER    EDATE       ! ending input date counter (YYYYDDD) in GMT
        INTEGER    EDATE_NEW   ! ending date based on met data
        INTEGER    ENLEN       ! length of the emissions inven name
        INTEGER    ETIME       ! ending input time counter (HHMMSS)  in GMT
        INTEGER    ETIME_NEW   ! ending time based on met data
        INTEGER    HRPOS       ! starting position in METCHECK array
        INTEGER    IOS         ! temporary I/O status
        INTEGER    JDATE       ! input date counter (YYYYDDD) in GMT
        INTEGER    JTIME       ! input time counter (HHMMSS)  in GMT
        INTEGER    LDATE       ! date from previous loop iteration
        INTEGER    MDATE       ! output date for monthly counties
        INTEGER    METNGRID    ! no. grid cells in met data
        INTEGER    METSTART    ! starting point for good met data
        INTEGER    MISSTEPS    ! no. contiguous missing hours of met data
        INTEGER    MAX_MISS    ! maximum number of missed steps
        INTEGER :: NDYCNTY = 0 ! no. counties using day averaging
        INTEGER :: NWKCNTY = 0 ! no. counties using week averaging
        INTEGER :: NMNCNTY = 0 ! no. counties using month averaging
        INTEGER :: NEPCNTY = 0 ! no. counties using episode averaging
        INTEGER    NINVARR     ! no. inventory variables to read
        INTEGER    NLINES      ! no. lines in met list file
        INTEGER    NMATX       ! size of ungridding matrix
        INTEGER    NSTEPS      ! number of time steps to process temperature data
        INTEGER    NEWSTEPS    ! no. time steps in episode based on met data
        INTEGER    NSTEP_MET   ! no. time steps in met file
        INTEGER :: OSRC = 0    ! number of sources outside grid
        INTEGER    POS         ! position in time step loop
        INTEGER    RDATE       ! date to read met file
        INTEGER    SDATE       ! output start date
        INTEGER    SDATE_NEW   ! output start date based on met data
        INTEGER    SDATE_MET   ! met file start date
        INTEGER    STIME       ! output start time
        INTEGER    STIME_NEW   ! output start time based on met data
        INTEGER    STIME_MET   ! met file start time
        INTEGER    TMPDAY      ! temporary day
        INTEGER    TMPMNTH     ! temporary month
        INTEGER    TSPREAD     ! time spread: difference between TZMAX and TZMIN
        INTEGER    TZONE       ! zone to determine output days
        INTEGER    TZMIN       ! minimum time zone in inventory
        INTEGER    TZMAX       ! maximum time zone in inventory
        INTEGER    WDATE       ! output date for weekly counties
        
        LOGICAL :: EFLAG    = .FALSE.  !  true: error found
        LOGICAL :: FIRSTMET = .TRUE.   !  true: processing first met file
        LOGICAL :: DUPWARN  = .TRUE.   !  true: print warning about overlapping met data
        LOGICAL :: DAYAVER  = .FALSE.  !  true: daily averaging
        LOGICAL :: WEEKAVER = .FALSE.  !  true: weekly averaging 
        LOGICAL :: MONAVER  = .FALSE.  !  true: monthly averaging
        LOGICAL :: EPIAVER  = .FALSE.  !  true: episode averaging
        LOGICAL :: LASTTIME = .FALSE.  !  true: final time step
        LOGICAL :: OFLAG    = .FALSE.  !  true: ungridding is 0 for some srcs
                
        CHARACTER(LEN=IOVLEN3) :: TVARNAME    !  temperature variable name
        CHARACTER(LEN=256)     :: CURFNM      !  current met file name
        CHARACTER(LEN=16)      :: CURLNM      !  current met logical file name
        CHARACTER(LEN=256)     :: DUPNAME     !  name of overlapping met file
        CHARACTER(LEN=200)     :: TEMPDIR     !  directory for output files
        CHARACTER(LEN=14)      :: DTBUF       !  date buffer
        CHARACTER(LEN=300)     :: MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'PREMOBL'  !  program name

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
        
C.........  Get episode starting date and time and run length
        MESG = 'Episode start date (YYYYDDD)'
        EPI_SDATE = ENVINT( 'EPI_STDATE', MESG, 0, IOS )
        
        MESG = 'Episode start time (HHMMSS)'
        EPI_STIME = ENVINT( 'EPI_STTIME', MESG, 0, IOS )

        MESG = 'Episode duration (HHMMSS)'
        EPI_RUNLEN = ENVINT( 'EPI_RUNLEN', MESG, 0, IOS )
        
        EPI_NSTEPS = EPI_RUNLEN / 10000

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )
                
C.........  Get the name of the activity to use for one run
        MESG = 'Temperature variable name'
        CALL ENVSTR( 'TVARNAME', MESG, 'TEMP1P5', TVARNAME, IOS )

C.........  Set default name of meterology file, depending on the name of the
C           temperature variable
C       TNAME = 'MET_CRO_2D'
C       IF( TVARNAME == 'TA' ) TNAME = 'MET_CRO_3D'

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

        TDEV = PROMPTFFILE(
     &          'Enter logical name for METEOROLOGY LIST file',
     &          .TRUE., .TRUE., 'METLIST', PROGNAME )

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

C.........  Read header of ungridding matrix...
        IF( .NOT. DESC3( UNAME ) ) THEN
            MESG = 'Could not get description for file ' // UNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
 
C.........  Store number of ungridding factors
        NMATX = NCOLS3D
 
C.........  Check dimensions of ungridding matrix...
        CALL CHKSRCNO( CATDESC, UNAME, NROWS3D, NSRC, EFLAG )
        
C.........  If the dimensions were in error, abort
        IF( EFLAG ) THEN
            MESG = 'Ungridding matrix is inconsistent with ' //
     &             'inventory.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 

C.........  Initialize reference grid with ungridding matrix
        CALL CHKGRID( UNAME, 'GMAT', 0, EFLAG )
 
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

C.........  Define the minimum and maximum time zones in the inventory
        TZMIN = MINVAL( TZONES )
        TZMAX = MAXVAL( TZONES )

C.........  Adjust TZMIN for possibility of daylight savings
        TZMIN = MAX( TZMIN - 1, 0 )

C.........  Calculate time spread based on min and max time zone
        TSPREAD = TZMAX - TZMIN

C.........  Calculate required starting date and time based on episode settings
        SDATE = EPI_SDATE

C.........  Earliest time required will be 6 a.m. in time zone closest to GMT
        STIME = 6 + TZMIN - TZONE          ! starting time in output time zone
        
C.........  Make sure the starting time is between 0 and 23
        IF( STIME < 0 ) THEN
            STIME = STIME + 24
        ELSE
            STIME = MOD( STIME, 24 )
        END IF
        STIME = STIME*10000

C.........  If the episode start time is earlier than our calculated start time,
C           we need to set the starting date back one day
        IF( EPI_STIME < STIME ) THEN
            CALL NEXTIME( SDATE, STIME, -24*10000 ) 
        END IF
        
C.........  Calculate required ending date and time
        EPI_EDATE = EPI_SDATE
        EPI_ETIME = EPI_STIME
        CALL NEXTIME( EPI_EDATE, EPI_ETIME, EPI_NSTEPS*10000 )
        
        EDATE = EPI_EDATE
        
C.........  Latest time required will be 5 a.m. in time zone farthest from GMT
        ETIME = 5 + TZMAX - TZONE        ! ending time in output time zone
        
C.........  Make sure the ending time is between 0 and 23
        IF( ETIME < 0 ) THEN
            ETIME = ETIME + 24
        ELSE
            ETIME = MOD( ETIME, 24 )
        END IF
        ETIME = ETIME*10000

C.........  If the episode ending time is later than calculated end time,
C           set the ending date forward one day
        IF( EPI_ETIME > ETIME ) THEN
            CALL NEXTIME( EDATE, ETIME, 24*10000 )
        END IF

C.........  Convert start and end dates and times back to GMT
        CALL NEXTIME( SDATE, STIME, TZONE*10000 )
        CALL NEXTIME( EDATE, ETIME, TZONE*10000 )

C.........  Find the total number of time steps
        NSTEPS = 1 + SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600

C.........  Allocate met checking array
        ALLOCATE( METCHECK( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METCHECK', PROGNAME )
        
        METCHECK = 0    ! array

C.........  Get number of lines in met list file
        NLINES = GETFLINE( TDEV, 'METLIST file' )
        
C.........  Allocate memory for storing the file
        ALLOCATE( METLIST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METLIST', PROGNAME )
        ALLOCATE( METLOGS( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METLOGS', PROGNAME )
        
        METLIST = ' '   ! array
        METLOGS = ' '
        
C.........  Store lines of METLIST file
        CALL RDLINES( TDEV, 'METLIST file', NLINES, METLIST )

        MESG = 'Checking meteorology files...'
        CALL M3MSG2( MESG )

        FIRSTMET = .TRUE.

C.........  Loop through all meteorology files
        DO N = 1, NLINES
            CURFNM = METLIST( N )

C.............  Reset duplicate warning - prints error once per met file
            DUPWARN = .TRUE.    

C.............  Skip any blank lines
            IF( CURFNM == ' ' ) CYCLE

C.............  Assign and store logical file name
            WRITE( CURLNM,94030 ) 'MET_FILE_', N
            METLOGS( N ) = CURLNM
            
C.............  Try to open file   
            IF( .NOT. OPNFULL3( CURLNM, FSREAD3, CURFNM, 
     &                          PROGNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not open meteorology file ' //
     &                 CURFNM( 1:LEN_TRIM( CURFNM ) )
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Try to get description from file            
            IF( .NOT. DESC3( CURLNM ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not get description of ' //
     &                 'meteorology file ' // 
     &                 CURFNM( 1:LEN_TRIM( CURFNM ) )
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Save description parameters
            SDATE_MET = SDATE3D
            STIME_MET = STIME3D
            NSTEP_MET = MXREC3D

C.............  Compare the met grid with previous settings
C.............  The subgrid parameters will be set here, if there is a subgrid
            CALL CHKGRID( CURLNM, 'GRID', 1, EFLAG )
            METNGRID = NGRID

C............. If the dimensions were in error, abort
            IF( EFLAG ) THEN
                MESG = 'Grids in ungridding matrix and meteorology ' //
     &                 'data are inconsistent.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Use start date, start time, and no. time steps to fill in METCHECK array
            HRPOS = SECSDIFF( SDATE, STIME, SDATE_MET, STIME_MET )
            HRPOS = HRPOS / 3600

C.............  If met data starts before episode, skip ahead to start of episode
            IF( HRPOS <= 0 ) THEN
                NSTEP_MET = NSTEP_MET + HRPOS
                HRPOS = 1
            ELSE
                HRPOS = HRPOS + 1
            END IF

C.............  Loop through time steps in met data and fill in METCHECK array            
            DO J = HRPOS, HRPOS + NSTEP_MET
                IF( J > NSTEPS ) EXIT

C.................  Check if met data overlaps previous file and print warning                
                IF( METCHECK( J ) /= 0 .AND. DUPWARN ) THEN
                    DUPNAME = METLIST( METCHECK( J ) )
                    MESG = 'WARNING: Time period of meteorology file '//
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) ) //
     &                     'overlaps that of file ' //
     &                     DUPNAME( 1:LEN_TRIM( DUPNAME ) ) // '.' //
     &                     CRLF() // BLANK10 // 'Data from ' // 
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) ) //
     &                     'will be used.' 
                    CALL M3MESG( MESG )
                    DUPWARN = .FALSE.
                END IF
                
                METCHECK( J ) = N
            END DO 

        END DO

C.........  Exit if there was a problem with the meteorology files
        IF( EFLAG ) THEN
            MESG = 'Problem checking meteorology files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Open group files
        CALL OPENGROUP( DDEV, WDEV, MDEV, EDEV )

C.........  Check for missing met data
        MESG = 'Checking for missing meteorology data...'
        CALL M3MSG2( MESG )
        
C.........  Check returned unit number to find out which temporal 
C           averaging is required; set maximum number of missing
C           time steps and error/warning message
        IF( DDEV > 0 ) THEN
            DAYAVER  = .TRUE.
            MAX_MISS = 24
            MESG = 'WARNING: Missing more than a full day of ' //
     &             'meteorology data.' // CRLF() // BLANK5 // 
     &             'Episode ending date and time will be adjusted.'
        END IF
        
        IF( WDEV > 0 ) THEN
            WEEKAVER = .TRUE.
            MAX_MISS = 24*7
            MESG = 'WARNING: Missing more than 7 days of ' //
     &             'meteorology data.' // CRLF() // BLANK5 //
     &             'Episode ending date and time will be adjusted.'
        END IF
        
        IF( MDEV > 0 ) THEN
            MONAVER  = .TRUE.
            MAX_MISS = 24*30
            MESG = 'WARNING: Missing more than 30 days of ' //
     &             'meteorology data.' // CRLF() // BLANK5 //
     &             'Episode ending date and time will be adjusted.'
        END IF
        
        IF( EDEV > 0 ) THEN
            MAX_MISS = NSTEPS
            EPIAVER  = .TRUE.
        END IF

        NEWSTEPS = 0
        METSTART = 0
        MISSTEPS = 0
        
        DO T = 1, NSTEPS
        
            IF( METCHECK( T ) /= 0 ) THEN
                NEWSTEPS = NEWSTEPS + MISSTEPS
                MISSTEPS = 0
                IF( METSTART == 0 ) THEN
                    NEWSTEPS = 1
                    CALL NEXTIME( SDATE, STIME, (T-1)*10000 )
                    METSTART = T
                ELSE
                    NEWSTEPS = NEWSTEPS + 1
                END IF
            ELSE
                IF( METSTART /= 0 ) THEN
                    MISSTEPS = MISSTEPS + 1
                    IF( MISSTEPS == MAX_MISS ) THEN
                        CALL M3MSG2( MESG )
                        EXIT
                    END IF
                END IF
            END IF
            
        END DO

        IF( NEWSTEPS == 0 ) THEN
            MESG = 'Meteorology data does not cover ' //
     &             'requested time period'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Reset number of steps and print messages
        IF( NEWSTEPS /= NSTEPS ) THEN
            MESG = 'WARNING: Resetting output start date, start ' //
     &             'time, or duration to match available meteorology' //
     &             ' data.'
            CALL M3MSG2( MESG )

C.............  Calculate new end date and time
            EDATE = SDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, NEWSTEPS*10000 )

C.............  Reset start and end times to deal with 6 am - 5 am restrictions
            SDATE_NEW = SDATE
            STIME_NEW = 6 + TZMIN          ! starting time in GMT
            STIME_NEW = STIME_NEW*10000

            IF( STIME > STIME_NEW ) THEN
                CALL NEXTIME( SDATE_NEW, STIME_NEW, 24*10000 )
            END IF
            
C.............  Adjust metstart for difference between old and new start
            METSTART = METSTART + 
     &            SECSDIFF( SDATE, STIME, SDATE_NEW, STIME_NEW ) / 3600
            
            SDATE = SDATE_NEW
            STIME = STIME_NEW
            
            EDATE_NEW = EDATE
            ETIME_NEW = 5 + TZMAX          ! ending time in GMT
            ETIME_NEW = ETIME_NEW*10000
            
            IF( ETIME < ETIME_NEW ) THEN
                CALL NEXTIME( EDATE_NEW, ETIME_NEW, -24*10000 )
            END IF
            
            EDATE = EDATE_NEW
            ETIME = ETIME_NEW
            
C.............  Calculate number of time steps
            NSTEPS = SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600
                        
C.............  Get date buffer field
            DTBUF = MMDDYY( SDATE )
            L = LEN_TRIM( DTBUF )
                
            WRITE( MESG,94050 )
     &        'Output Time Zone:', TZONE,           CRLF() // BLANK5 //
     &        '  Met Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '  Met Start Time:', STIME + TZONE*10000,'HHMMSS'// 
     &                                              CRLF() // BLANK5 //
     &        '       Time Step:', 1,    'hour'  // CRLF() // BLANK5 //
     &        '    Met Duration:', NSTEPS, 'hours'   
     
            CALL M3MSG2( MESG )
        END IF
           
C.........  Get sources and counties from SPDSUM file
        ALLOCATE( COUNTYSRC( NSRC, 2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'COUNTYSRC', PROGNAME )

        COUNTYSRC = 0
        CALL RDSPDSRC( PDEV, NSRC, COUNTYSRC )

C.........  Assign temporal grouping values to each county
        IF( DAYAVER  ) CALL RDGROUPS( DDEV, COUNTYSRC, DAILY, NDYCNTY,
     &                                              .NOT. WEEKAVER .AND. 
     &                                              .NOT. MONAVER .AND. 
     &                                              .NOT. EPIAVER )
        IF( WEEKAVER ) CALL RDGROUPS( WDEV, COUNTYSRC, WEEKLY, NWKCNTY,
     &                                              .NOT. MONAVER .AND. 
     &                                              .NOT. EPIAVER )
        IF( MONAVER ) CALL RDGROUPS( MDEV, COUNTYSRC, MONTHLY, NMNCNTY,
     &                                              .NOT. EPIAVER )
        IF( EPIAVER ) CALL RDGROUPS( EDEV, COUNTYSRC, EPISLEN, NEPCNTY, 
     &                                              .TRUE. )

C.........  Allocate memory for storing county codes
        ALLOCATE( DYCODES( NDYCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DYCODES', PROGNAME )
        ALLOCATE( WKCODES( NWKCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WKCODES', PROGNAME )
        ALLOCATE( MNCODES( NMNCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MNCODES', PROGNAME )
        ALLOCATE( EPCODES( NEPCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EPCODES', PROGNAME )

C.........  Fill out individual FIPS arrays from COUNTYSRC array
        CALL ASGNGRPS( COUNTYSRC )

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
        ALLOCATE( TA( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( DAYBEGT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        ALLOCATE( DAYENDT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        ALLOCATE( LDAYSAV( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDAYSAV', PROGNAME )
        ALLOCATE( TASRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASRC', PROGNAME )
        ALLOCATE( NDAYSRC( NSRC,24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NDAYSRC', PROGNAME )

        NDAYSRC = 0   ! array

C.........  Create array of which sources are affected by daylight savings
        CALL GETDYSAV( NSRC, IFIP, LDAYSAV )

C.........  Read ungridding matrix 
        CALL RDUMAT( UNAME, NSRC, NMATX, NMATX, 
     &               UMAT( 1 ), UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ) )

C.........  Allocate memory for storing temperature profiles
        ALLOCATE( TKHOUR( NSRC, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKHOUR', PROGNAME )
        ALLOCATE( TDYCNTY( NDYCNTY, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TDYCNTY', PROGNAME )
        ALLOCATE( TWKCNTY( NWKCNTY, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TWKCNTY', PROGNAME )
        ALLOCATE( TMNCNTY( NMNCNTY, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMNCNTY', PROGNAME )
        ALLOCATE( TEPCNTY( NEPCNTY, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TEPCNTY', PROGNAME )
        
C.........  Preprocess dates, times, and time zones.  Write report for
C           time zones here for zones that do not have a complete day of
C           data in beginning or end. 
        CALL CHKFULLDY( NSRC, SDATE, STIME, EDATE, ETIME, 
     &                  TZONES, LDAYSAV )

C.........  Get output directory information from the environment
        MESG = 'Path where hourly temperature files will be written'
        CALL ENVSTR( 'SMK_TEMPATH', MESG, '.', TEMPDIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: Temperature files being placed in ' //
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_TEMPATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF
        
C.........  Open output files
        IF( DAYAVER ) THEN
            DNAME = 'DAYTEMP'
            CALL OPENSHOUR( ENAME, 'daily', SDATE, STIME, TVARNAME, 
     &                      NDYCNTY, TEMPDIR, DNAME )
        END IF
        IF( WEEKAVER ) THEN
            WNAME = 'WEEKTEMP'
            CALL OPENSHOUR( ENAME, 'weekly', SDATE, STIME, TVARNAME,
     &                      NWKCNTY, TEMPDIR, WNAME )
        END IF
        IF( MONAVER ) THEN
            MNAME = 'MONTEMP'
            CALL OPENSHOUR( ENAME, 'monthly', SDATE, STIME, TVARNAME,
     &                      NMNCNTY, TEMPDIR, MNAME )
        END IF
        IF( EPIAVER ) THEN
            PNAME = 'EPITEMP'
            CALL OPENSHOUR( ENAME, 'episode', SDATE, STIME, TVARNAME,
     &                      NEPCNTY, TEMPDIR, PNAME )
        END IF

C.........  Process temperature information...

        L = LEN_TRIM( TVARNAME )
        MESG = 'Processing temperature data using variable "' //
     &         TVARNAME( 1:L ) // '" ...'
        CALL M3MSG2( MESG )

C.........  Loop through days/hours of temperature files
        DDATE = SDATE
        DTIME = 0
        WDATE = SDATE
        MDATE = SDATE
        
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        
        DO T = METSTART, METSTART + NSTEPS

            POS = T - METSTART + 1
            RDATE = JDATE

C.............  Set correct met file for this time step

C.............  If no met file for this time step, use previous or next day
            IF( METCHECK( T ) == 0 ) THEN
                IF( T - 24 > 0 ) THEN
                    T2 = T - 24
                    RDATE = RDATE - 1
                ELSE
                    T2 = T + 24
                    RDATE = RDATE + 1
                END IF
            ELSE
                T2 = T
            END IF

C.............  If we still don't have met data for current day, quit
            IF( METCHECK( T2 ) == 0 ) THEN
                WRITE( MESG,94010 ) 'Could not get met file for ', 
     &                              JDATE, ' and ', JTIME
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            ELSE
                CURFNM = METLIST( METCHECK( T2 ) )
            END IF

C.............  Get logical file name 
            CURLNM = METLOGS( METCHECK( T2 ) )

C.............  Open the meteorology data file
            IF ( .NOT. OPNFULL3( CURLNM, FSREAD3, CURFNM, 
     &                           PROGNAME ) ) THEN
                MESG = 'Could not open meteorology file ' // CURFNM
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

C.............  Read current temperature file
            IF ( .NOT. READ3( CURLNM, TVARNAME, 1, 
     &                        RDATE, JTIME, TA ) ) THEN
                L = LEN_TRIM( TVARNAME )
                MESG = 'Could not read ' // TVARNAME( 1:L ) //
     &                 ' from ' // CURFNM 
                CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )

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

C.............  Create hourly temperature array by source
            CALL HOURTEMP( NSRC, JTIME, DAYBEGT, TASRC, COUNTYSRC, 
     &                     MINTEMP, MAXTEMP, NDAYSRC, TKHOUR )

C.............  Make sure we've waited long enough to catch all time zones
            IF( POS > TSPREAD ) THEN
            	
C.................  Adjust time step for 24-hour arrays
                ARRAYPOS = MOD( POS - TSPREAD, 24 )
                IF( ARRAYPOS == 0 ) ARRAYPOS = 24

C.................  Process daily averages
                IF( DAYAVER ) THEN

C.....................  Average temperatures across county group                
                    CALL AVERTEMP( NSRC, NDYCNTY, DYCODES, 
     &                             COUNTYSRC( :,1 ), ARRAYPOS, TKHOUR, 
     &                             TDYCNTY, NDAYSRC )

C.....................  Open new file if necessary (1 day per output file)
                    IF( POS > 24 .AND. 
     &                  MOD( POS - TSPREAD, 24 ) == 1 ) THEN
                        IF( .NOT. CLOSE3( DNAME ) ) THEN
                            MESG = 'Could not close file ' // 
     &                              DNAME( 1:LEN_TRIM( DNAME ) )
                            CALL M3EXIT( PROGNAME, DDATE, 0, MESG, 2 )
                        END IF
                        
                        CALL OPENSHOUR( ENAME, 'daily', DDATE, DTIME, 
     &                                  TVARNAME, NDYCNTY, TEMPDIR, 
     &                                  DNAME )
                    END IF                            

C.....................  Write time step to file
                    CALL WRSHOUR( DNAME, DDATE, DTIME, NDYCNTY, 
     &                            ARRAYPOS, DYCODES, TDYCNTY )
     
                END IF

C.............  If current output day is Saturday, process weekly averages
                IF( WEEKAVER .AND. WKDAY( DDATE ) == 6 ) THEN

C.....................  Average temperatures across county group 
                    CALL AVERTEMP( NSRC, NWKCNTY, WKCODES, 
     &                             COUNTYSRC( :,1 ), ARRAYPOS, TKHOUR,
     &                             TWKCNTY, NDAYSRC )

C.....................  Open new file if necessary (one week per output file)
                    IF( .NOT. ISOPEN( WNAME ) ) THEN
                        CALL OPENSHOUR( ENAME, 'weekly', WDATE, DTIME, 
     &                                  TVARNAME, NWKCNTY, TEMPDIR, 
     &                                  WNAME )
                    END IF       

C.....................  Write time step to file
                    CALL WRSHOUR( WNAME, WDATE, DTIME, NWKCNTY, 
     &                            ARRAYPOS, WKCODES, TWKCNTY )
                    
                    IF( DTIME == 230000 ) THEN
                    	
C.........................  Store current date to use for next week
                        WDATE = DDATE + 1
C                        CALL NEXTIME( WDATE, 0, 10000 )

C.........................  Close current output file                        
                        IF( .NOT. CLOSE3( WNAME ) ) THEN
                            MESG = 'Could not close file ' // 
     &                              WNAME( 1:LEN_TRIM( WNAME ) )
                            CALL M3EXIT( PROGNAME, WDATE, 0, MESG, 2 )
                        END IF
                    END IF                     
                END IF
                
C.................  If last day of month, process monthly averages
                CALL DAYMON( DDATE, CURRMNTH, CURRDAY )
                CALL DAYMON( DDATE + 1, TMPMNTH, TMPDAY )
            
                IF( MONAVER .AND. TMPMNTH /= CURRMNTH ) THEN

C.....................  Average temperatures across county group 
                    CALL AVERTEMP( NSRC, NMNCNTY, MNCODES, 
     &                             COUNTYSRC( :,1 ), ARRAYPOS, TKHOUR,
     &                             TMNCNTY, NDAYSRC )

C.....................  Open new file if necessary (one month per output file)
                    IF( .NOT. ISOPEN( MNAME ) ) THEN                        
                        CALL OPENSHOUR( ENAME, 'monthly', MDATE, DTIME, 
     &                                  TVARNAME, NMNCNTY, TEMPDIR, 
     &                                  MNAME )
                    END IF       

C.....................  Write time step to file                     
                    CALL WRSHOUR( MNAME, MDATE, DTIME, NMNCNTY, 
     &                            ARRAYPOS, MNCODES, TMNCNTY )               
                    
                    IF( DTIME == 230000 ) THEN
                    	
C.........................  Store current date to use for next month
                        MDATE = DDATE + 1
C                        CALL NEXTIME( MDATE, 0, 10000 )

C.........................  Close current output file   
                        IF( .NOT. CLOSE3( MNAME ) ) THEN
                            MESG = 'Could not close file ' // 
     &                              MNAME( 1:LEN_TRIM( MNAME ) )
                            CALL M3EXIT( PROGNAME, MDATE, 0, MESG, 2 )
                        END IF
                    END IF
                END IF

            END IF    ! time zone check

            LASTTIME = ( POS == NSTEPS )

C.............  Final steps on last time through
            IF( LASTTIME ) THEN

C.................  Process remaining week temperatures only if current date is 
C                   not Saturday; otherwise, this has already been handled above
                IF( WEEKAVER .AND. WKDAY( DDATE ) /= 6 ) THEN
                    IF( .NOT. ISOPEN( WNAME ) ) THEN                    
                        CALL OPENSHOUR( ENAME, 'weekly', WDATE, 0, 
     &                                  TVARNAME, NWKCNTY, TEMPDIR, 
     &                                  WNAME )
                    END IF
     
                    DO K = 1, 24                    
                        CALL AVERTEMP( NSRC, NWKCNTY, WKCODES, 
     &                                 COUNTYSRC( :,1 ), K, TKHOUR,
     &                                 TWKCNTY, NDAYSRC )
                       
                        CALL WRSHOUR( WNAME, WDATE, ( K-1 )*10000, 
     &                                NWKCNTY, K, WKCODES, TWKCNTY )
                    END DO
                END IF

C.................  Process remaining month temperatures only if current date is
C                   not the last day of the month
                IF( MONAVER .AND. TMPMNTH == CURRMNTH ) THEN
                    IF( .NOT. ISOPEN( MNAME ) ) THEN
                        CALL OPENSHOUR( ENAME, 'monthly', MDATE, 0, 
     &                                  TVARNAME, NMNCNTY, TEMPDIR, 
     &                                  MNAME )                    	
                    END IF
                	
                    DO K = 1, 24
                        CALL AVERTEMP( NSRC, NMNCNTY, MNCODES, 
     &                                 COUNTYSRC( :,1 ), K, TKHOUR,
     &                                 TMNCNTY, NDAYSRC )
                        
                        CALL WRSHOUR( MNAME, MDATE, ( K-1 )*10000, 
     &                                NMNCNTY, K, MNCODES, TMNCNTY )  
                    END DO
                END IF

C.................  Output episode averaged temperatures
                IF( EPIAVER ) THEN
                    DO K = 1, 24
                        CALL AVERTEMP( NSRC, NEPCNTY, EPCODES, 
     &                                 COUNTYSRC( :,1 ), K, TKHOUR,
     &                                 TEPCNTY, NDAYSRC )
     
                        CALL WRSHOUR( PNAME, SDATE, ( K-1 )*10000, 
     &                                NEPCNTY, K, EPCODES, TEPCNTY )
                    END DO    
                END IF

C.................  Count sources outside the grid
                IF( .NOT. OFLAG ) THEN
                    DO S = 1, NSRC

                	IF( UMAT( S ) == 0 ) THEN
                        OFLAG = .TRUE.
                        OSRC = OSRC + 1
                	END IF

                    END DO
                END IF

            END IF   ! last time through loop

C.............  Increment output time
            IF( POS > TSPREAD ) THEN
                CALL NEXTIME( DDATE, DTIME, 10000 )
            END IF

C.............  Increment loop time
            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, 10000 )

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

94030   FORMAT( A, I5 )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )
                
        END PROGRAM PREMOBL