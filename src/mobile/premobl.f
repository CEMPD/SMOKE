
        PROGRAM PREMOBL

C***********************************************************************
C  program body starts at line 205
C
C  DESCRIPTION:
C       Creates county-based 24-hour temperature and mixing ratio 
C       profiles based on gridded meteorology data. Temperatures and 
C       mixing ratios can be averaged across counties and different time periods. 
C       Also outputs hourly barometric pressure values. Tries to 
C       account for missing meteorology data as much as possible.
C
C  PRECONDITIONS REQUIRED:
C       Program Mbsetup has been run
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
C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: IFIP, TZONES

C...........   This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, CATDESC, CATLEN, NIACT, NSRC

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: MINTEMP, MAXTEMP, TASRC, QVSRC, PRESSRC,
     &                    TKHOUR, QVHOUR, BPHOUR, 
     &                    DYCODES, WKCODES, MNCODES, EPCODES,
     &                    TDYCNTY, TWKCNTY, TMNCNTY, TEPCNTY,
     &                    QVDYCNTY, QVWKCNTY, QVMNCNTY, QVEPCNTY,
     &                    BPDYCNTY, BPWKCNTY, BPMNCNTY, BPEPCNTY

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C.........  This module is used for MOBILE6 setup information        
        USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN

        IMPLICIT NONE
        
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
        INCLUDE 'M6CNST3.EXT'   !  MOBILE6 constants

C...........   EXTERNAL FUNCTIONS and their descriptions:        
        CHARACTER(2)    CRLF
        INTEGER         GETIFDSC
        INTEGER         GETFLINE
        INTEGER         ENVINT
        REAL            ENVREAL
        INTEGER         INDEX1
        CHARACTER(14)   MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         SECSDIFF
        LOGICAL         SETENVVAR
        INTEGER         WKDAY
        
        EXTERNAL     CRLF, GETIFDSC, GETFLINE, ENVINT, 
     &               ENVREAL, INDEX1, MMDDYY, PROMPTFFILE, PROMPTMFILE, 
     &               SECSDIFF, SETENVVAR, WKDAY
        
C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Gridded meteorology data (dim: NGRID)
        REAL   , ALLOCATABLE :: TA( : )   !  one layer of temperature
        REAL   , ALLOCATABLE :: QV( : )   !  water vapor mixing ratio
        REAL   , ALLOCATABLE :: PRES( : ) !  pressure

C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : ) ! Contiguous ungridding matrix

C...........  Allocatable per-source arrays
        INTEGER, ALLOCATABLE :: DAYBEGT ( : ) ! daily start time HHMMSS
        INTEGER, ALLOCATABLE :: DAYENDT ( : ) ! daily end time HHMMSS
        LOGICAL, ALLOCATABLE :: LDAYSAV ( : ) ! true: src uses daylight time
        INTEGER, ALLOCATABLE :: COUNTYSRC ( :,: ) ! FIPS code and averaging values

C...........  Allocatable arrays for met data
        INTEGER, ALLOCATABLE :: METDAYS( : ) ! dimension: nsteps in episode, 
                                             ! value indicates which met data file covers that hour
        CHARACTER(256), ALLOCATABLE :: METLIST( : ) ! listing of met file names

        INTEGER, ALLOCATABLE :: NDAYSRC( :,: )  ! number of days to average by source

C...........  Array that contains the names of the inventory variables needed 
C             for this program
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C...........   File units and logical names:
        INTEGER      DDEV  ! unit number for daily group file
        INTEGER      EDEV  ! unit number for episode group file
        INTEGER      IDEV  ! tmp unit number if ENAME is map file
        INTEGER      LDEV  ! unit number for log file
        INTEGER      MDEV  ! unit number for monthly group file
        INTEGER      PDEV  ! unit number for speeds summary file (SPDSUM)
        INTEGER      SDEV  ! unit number for ASCII inventory file
        INTEGER      TDEV  ! unit number for meteorology list file
        INTEGER      WDEV  ! unit number for weekly group file

        CHARACTER(16) ANAME ! logical name for mobile ASCII inventory file
        CHARACTER(16) ENAME ! logical name for mobile I/O API inventory file
        CHARACTER(16) DNAME ! logical name for daily output ungridded hourly temps
        CHARACTER(16) INAME ! tmp name for inven file of unknown fmt
        CHARACTER(16) METNAME ! logical name for meteorology files
        CHARACTER(16) MNAME ! logical name for monthly output hourly temps
        CHARACTER(16) PNAME ! logical name for episode output hourly temps
        CHARACTER(16) UNAME ! logical name for ungridding-matrix input file
        CHARACTER(16) WNAME ! logical name for weekly output hourly temps
                
C...........   Other local variables:
        INTEGER    I, J, K, L, N, S, T, T2, V  ! Counters and pointers

        INTEGER    EPI_SDATE      ! episode start date from E.V. (YYYYDDD)
        INTEGER    EPI_STIME      ! episode start time from E.V. (HHMMSS)
        INTEGER    EPI_RUNLEN     ! episode duration from E.V. (HHMMSS)
        INTEGER    EPI_NSTEPS     ! episode number of time steps
        INTEGER    EPI_EDATE      ! episode ending date based on ERUNLEN
        INTEGER    EPI_ETIME      ! episode ending time based on ERUNLEN
        
        INTEGER    ARRAYPOS    ! position in 24-hour arrays
        INTEGER    DAY         ! tmp day of week number
        INTEGER    DDATE       ! output date for daily counties
        INTEGER    DUMMYTIME   ! dummy time variable to use in calls to NEXTIME
        INTEGER    EDATE       ! ending input date counter (YYYYDDD) in GMT
        INTEGER    ETIME       ! ending input time counter (HHMMSS)  in GMT
        INTEGER    FILENUM     ! file number of current meteorology file
        INTEGER    IOS         ! temporary I/O status
        INTEGER    JDATE       ! input date counter (YYYYDDD) in GMT
        INTEGER    JTIME       ! input time counter (HHMMSS)  in GMT
        INTEGER    LDATE       ! date from previous loop iteration
        INTEGER    MDATE       ! output date for monthly counties
        INTEGER    MONTH       ! current month
        INTEGER    METNGRID    ! no. grid cells in met data
        INTEGER :: NDYCNTY = 0 ! no. counties using day averaging
        INTEGER :: NWKCNTY = 0 ! no. counties using week averaging
        INTEGER :: NMNCNTY = 0 ! no. counties using month averaging
        INTEGER :: NEPCNTY = 0 ! no. counties using episode averaging
        INTEGER    NINVARR     ! no. inventory variables to read
        INTEGER    NLINES      ! no. lines in met list file
        INTEGER    NMATX       ! size of ungridding matrix
        INTEGER    NSTEPS      ! number of time steps to process temperature data
        INTEGER :: OSRC = 0    ! number of sources outside grid
        INTEGER    OTIME       ! output time in local time
        INTEGER    POS         ! position in time step loop
        INTEGER    RDATE       ! date to read met file
        INTEGER    SDATE       ! output start date
        INTEGER    STIME       ! output start time
        INTEGER    TDATE       ! temporary julian date
        INTEGER    TTIME       ! temporary time
        INTEGER    TMPMNTH     ! temporary month
        INTEGER    TSPREAD     ! time spread: difference between TZMAX and TZMIN
        INTEGER    TZONE       ! zone to determine output days
        INTEGER    TZMIN       ! minimum time zone in inventory
        INTEGER    TZMAX       ! maximum time zone in inventory
        INTEGER    WDATE       ! output date for weekly counties
        
        LOGICAL :: EFLAG    = .FALSE.  !  true: error found
        LOGICAL :: GRID_ERR = .FALSE.  !  true: error found in grid settings
        LOGICAL :: DAYAVER  = .FALSE.  !  true: daily averaging
        LOGICAL :: WEEKAVER = .FALSE.  !  true: weekly averaging 
        LOGICAL :: MONAVER  = .FALSE.  !  true: monthly averaging
        LOGICAL :: EPIAVER  = .FALSE.  !  true: episode averaging
        LOGICAL :: OFLAG    = .FALSE.  !  true: ungridding is 0 for some srcs
        LOGICAL :: DAYOPEN  = .FALSE.  !  true: daily temp file is open
        LOGICAL :: WEEKOPEN = .FALSE.  !  true: weekly temp file is open
        LOGICAL :: MONOPEN  = .FALSE.  !  true: monthly temp file is open
        LOGICAL :: FILEOPEN = .FALSE.  !  true: met file is open
        LOGICAL :: FND_DATA = .FALSE.  !  true: found met data for this hour
        LOGICAL :: ALT_DATA = .FALSE.  !  true: using alternate data for this hour
        
        CHARACTER(512)     METFILE     !  tmp physical file name
        CHARACTER(512)     PREVFILE    !  previous physical file name
        CHARACTER(IOVLEN3) TVARNAME    !  temperature variable name
        CHARACTER(IOVLEN3) PRESNAME    !  pressure variable name
        CHARACTER(IOVLEN3) MIXNAME     !  mixing ratio name
        CHARACTER(200)     TEMPDIR     !  directory for output files
        CHARACTER(300)     MESG        !  message buffer

        CHARACTER(16) :: PROGNAME = 'PREMOBL'  !  program name

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
                
C.........  Get the name of the temperature variable
        MESG = 'Temperature variable name'
        CALL ENVSTR( 'TVARNAME', MESG, 'TEMPG', TVARNAME, IOS )

C.........  Set default names for additional variables
        PRESNAME = 'PRES'
        MIXNAME = 'QV'

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

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

        TDEV = PROMPTFFILE(
     &          'Enter logical name for METEOROLOGY LIST file',
     &          .TRUE., .TRUE., 'METLIST', PROGNAME )

        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .TRUE., .TRUE., 'SPDSUM', PROGNAME )
     
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

C.........  Create note about time zone expected in meteorology file
        WRITE( MESG, 94000 ) 'NOTE: Time stamps of input meteorology '//
     &                       'files are assumed to be in GMT.'
        CALL M3MSG2( MESG )
 
C.........  Set inventory variables to read
        NINVARR = 3
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'CSOURC'
        IVARNAMS( 3 ) = 'TZONES'

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

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

C.........  Get number of lines in met list file
        NLINES = GETFLINE( TDEV, 'METLIST file' )
        
C.........  Allocate memory for storing the met file information
        ALLOCATE( METLIST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METLIST', PROGNAME )
        ALLOCATE( METDAYS( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METDAYS', PROGNAME )
        
        METLIST = ' '    ! array
        METDAYS = 0
        
C.........  Store lines of METLIST file
        CALL RDLINES( TDEV, 'METLIST file', NLINES, METLIST )
        CLOSE( TDEV )

        MESG = 'Checking meteorology files...'
        CALL M3MSG2( MESG )

        METNAME = 'METFILE'

C.........  Loop through all meteorology files
        DO N = 1, NLINES

C.............  Close previous file if needed
            IF( FILEOPEN ) THEN
                IF( .NOT. CLOSE3( METNAME ) ) THEN
                    MESG = 'Could not close meteorology file ' //
     &                     TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    FILEOPEN = .FALSE.
                END IF
            END IF

C.............  Get physical file name for current iteration        
            METFILE = METLIST( N )  

C.............  Skip any blank lines
            IF( METFILE == ' ' ) CYCLE

C.............  Set logical file name
            IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not set logical file name ' //
     &                 'for file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Try to open file   
            IF( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not open meteorology file ' //
     &                 TRIM( METFILE )
                CALL M3MESG( MESG )
                CYCLE
            ELSE
                FILEOPEN = .TRUE.
            END IF

C.............  Try to get description from file            
            IF( .NOT. DESC3( METNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not get description of ' //
     &                 'meteorology file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Check that requested variables are in the file
            J = INDEX1( TVARNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TRIM( TVARNAME ) // 
     &                 ' in file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
            END IF
            
            J = INDEX1( PRESNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TRIM( PRESNAME ) //
     &                 ' in file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
            END IF
            
            J = INDEX1( MIXNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TRIM( MIXNAME ) //
     &                 ' in file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
            END IF

            IF( EFLAG ) CYCLE

C.............  Initialize (or check) reference grid with meteorology data
            CALL CHKGRID( METNAME, 'GRID', 0, GRID_ERR )
            METNGRID = NGRID

C............. If the dimensions were in error, print message and cycle
            IF( GRID_ERR ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Grid in meteorology file ' //
     &                 TRIM( METFILE ) // ' is inconsistent ' //
     &                 'with previous files.' 
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Build METDAYS array linking met files to a specific date and time

C............. Find starting position in METDAYS array
            POS = 1 + SECSDIFF( SDATE, STIME, SDATE3D, STIME3D ) / 3600
            
C.............  Make sure the met data is within the episode            
            IF( POS + MXREC3D - 1 <= 0 .OR. POS > NSTEPS ) CYCLE

C.............  Fill in array for number of steps in current met file            
            DO L = POS, POS + MXREC3D - 1
                IF( L <= 0 ) CYCLE
                IF( L > NSTEPS ) EXIT
                METDAYS( L ) = N
            END DO

        END DO

C.........  Exit if there was a problem with the meteorology files
        IF( EFLAG ) THEN
            MESG = 'Problem checking meteorology files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Open group files
        CALL OPENGROUP( DDEV, WDEV, MDEV, EDEV )

C.........  Check group file unit numbers to see which types of averaging are needed
        IF( DDEV > 0 ) DAYAVER  = .TRUE.
        IF( WDEV > 0 ) WEEKAVER = .TRUE.
        IF( MDEV > 0 ) MONAVER  = .TRUE.
        IF( EDEV > 0 ) EPIAVER  = .TRUE.

C.........  Check for missing meteorology data
        MESG = 'Checking for missing meteorology data...'
        CALL M3MSG2( MESG )

C.........  Check that all days are covered
        IF( DAYAVER ) THEN
            T = 1
            
C.............  Loop over all days in episode
            DO

C.................  Make sure we're within episode bounds
                IF( T > NSTEPS ) EXIT

C.................  Check all hours in current day, accounting for time zone spread            
                DO I = 0, 23 + TSPREAD
                    K = T + I

C.....................  Double check episode bounds
                    IF( K > NSTEPS ) EXIT

C.....................  If no met data for current step, try to find data
                    IF( METDAYS( K ) == 0 ) THEN
 
C.........................  Try 24 hours back
                        IF( K - 24 > 0 ) THEN
                            IF( METDAYS( K - 24 ) > 0 ) THEN
                                METDAYS( K ) = - METDAYS( K - 24 )
                                CYCLE
                            END IF

C.........................  Try 24 hours forward
                        ELSE IF( K + 24 < NSTEPS ) THEN
                            IF( METDAYS( K + 24 ) > 0 ) THEN
                                METDAYS( K ) = - METDAYS( K + 24 )
                                CYCLE
                            END IF
                        END IF

C.........................  No data available, exit with error                    
                        MESG = 'Meteorology data does not cover ' //
     &                         'requested episode.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
                END DO   ! end loop over hours in day

C.................  Skip to start of next day
                T = T + 24
            END DO   ! end loop over days in episode
        END IF   ! end day averaging check

C.........  Check that all weeks are covered
        IF( WEEKAVER ) THEN
            T = 1
            JDATE = SDATE
            JTIME = STIME
            
C.............  Loop over all weeks in episode
            DO

C.................  Check episode bounds
                IF( T > NSTEPS ) EXIT

C.................  Check hours in day                
                DO I = 0, 23 + TSPREAD
                    K = T + I
                    
C.....................  Double check episode bounds
                    IF( K > NSTEPS ) EXIT

C.....................  If no met data for current step, try to find data           
                    IF( METDAYS( K ) <= 0 ) THEN
                        FND_DATA = .FALSE.

C.........................  Loop through previous days in week
                        DO J = 1, 6
                            TDATE = JDATE
                            TTIME = JTIME

C.............................  Check episode bounds
                            IF( K - J*24 > 0 ) THEN

C.................................  Make sure it's still the same week
                                CALL NEXTIME( TDATE, TTIME, -J*240000 )
                                IF( WKDAY( TDATE ) /= 6 ) THEN
                                    IF( METDAYS( K - J*24 ) > 0 ) THEN
                                        FND_DATA = .TRUE.
                                        EXIT
                                    END IF

C.................................  Otherwise, it's the previous week, so exit
                                ELSE
                                    EXIT
                                END IF

C.............................  Otherwise, it's too far back, exit
                            ELSE
                                EXIT
                            END IF
                        END DO   ! end loop over previous days

C.........................  Skip rest of loop if we've found data                    
                        IF( FND_DATA ) CYCLE

C.........................  Loop through remaining days in week                    
                        DO J = 1, 6
                            TDATE = JDATE
                            TTIME = JTIME
          
C.............................  Check episode bounds              
                            IF( K + J*24 < NSTEPS ) THEN

C.................................  Make sure it's still the same week
                                CALL NEXTIME( TDATE, TTIME, J*240000 )
                                IF( WKDAY( TDATE ) /= 7 ) THEN
                                    IF( METDAYS( K + J*24 ) > 0 ) THEN
                                        FND_DATA = .TRUE.
                                        EXIT
                                    END IF
                                ELSE

C.................................  Otherwise, it's the next week, so exit
                                    EXIT
                                END IF  
                                
C.............................  Otherwise, it's too far forward, exit
                            ELSE
                                EXIT      
                            END IF
                        END DO   ! end loop over remaining days
                    
                        IF( FND_DATA ) CYCLE

C.....................  Still no data, exit with error                         
                        MESG = 'Meteorology data does not cover ' //
     &                         'requested episode.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
                END DO   ! end loop over hours in day
            
C.................  Advance to next Sunday at midnight            
                DO
                    T = T + 1
                    CALL NEXTIME( JDATE, JTIME, 10000 )
                    IF( WKDAY( JDATE ) == 7 .AND. JTIME == 0 ) EXIT
                END DO
            END DO   ! end loop over weeks in episode
        END IF   ! end week averaging check

C.........  Check that all months are covered
        IF( MONAVER ) THEN
            T = 1
            JDATE = SDATE
            JTIME = STIME
            CALL DAYMON( JDATE, MONTH, DAY )

C.............  Loop over all months in episode     
            DO

C.................  Check episode bounds
                IF( T > NSTEPS ) EXIT

C.................  Check hours in day            
                DO I = 0, 23 + TSPREAD
                    K = T + I
                    
C.....................  Check episode bounds
                    IF( K > NSTEPS ) EXIT

C.....................  Check for met data at this step                
                    IF( METDAYS( K ) <= 0 ) THEN
                        FND_DATA = .FALSE.

C.........................  Loop through previous days in month
                        DO J = 1, 30
                            TDATE = JDATE
                            TTIME = JTIME

C.............................  Check episode bounds                        
                            IF( K - J*24 > 0 ) THEN

C.................................  Make sure it's still the same month
                                CALL NEXTIME( TDATE, TTIME, -J*240000 )
                                CALL DAYMON( TDATE, TMPMNTH, DAY )
                                IF( TMPMNTH == MONTH ) THEN
                                    IF( METDAYS( K - J*24 ) > 0 ) THEN
                                        FND_DATA = .TRUE.
                                        EXIT
                                    END IF
                                    
C.................................  Otherwise, it's the previous month, so exit
                                ELSE
                                    EXIT
                                END IF
                                
C.............................  Otherwise, it's too far back, exit
                            ELSE
                                EXIT
                            END IF
                        END DO   ! end loop over previous days
 
C.........................  Skip rest of loop if we've found data                        
                        IF( FND_DATA ) CYCLE

C.........................  Loop through remaining days in month                    
                        DO J = 1, 30
                            TDATE = JDATE
                            TTIME = JTIME
                        
C.............................  Check episode bounds   
                            IF( K + J*24 < NSTEPS ) THEN

C.................................  Make sure it's still the same month
                                CALL NEXTIME( TDATE, TTIME, J*240000 )
                                CALL DAYMON( TDATE, TMPMNTH, DAY )
                                IF( TMPMNTH == MONTH ) THEN
                                    IF( METDAYS( K + J*24 ) > 0 ) THEN
                                        FND_DATA = .TRUE.
                                        EXIT
                                    END IF
                                    
C.................................  Otherwise, it's the next month, so exit
                                ELSE
                                    EXIT
                                END IF    
                                
C.............................  Otherwise, it's too far forward, exit
                            ELSE
                                EXIT        
                            END IF
                        END DO   ! end loop over remaining days
                    
                        IF( FND_DATA ) CYCLE

C.....................  Still no data, exit with error                         
                        MESG = 'Meteorology data does not cover ' //
     &                         'requested episode.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
                END DO   ! end loop over hours in day
            
C.................  Advance T to next first day of next month at midnight            
                DO
                    T = T + 1
                    CALL NEXTIME( JDATE, JTIME, 10000 )
                    CALL DAYMON( JDATE, MONTH, DAY )
                    IF( DAY == 1 .AND. JTIME == 0 ) EXIT
                END DO
            END DO   ! end loop over months in episode
        END IF   ! end month averaging check
        
C.........  Check that the episode is covered
        IF( EPIAVER ) THEN

C.............  Loop over hours in day accounting for time zones
            DO T = 1, 24 + TSPREAD
                K = T

C.................  Check episode bounds
                IF( K > NSTEPS ) EXIT
                
C.................  Check for met data at this step
                IF( METDAYS( K ) <= 0 ) THEN
                    FND_DATA = .FALSE.

C.....................  Loop through additional days in episode
                    DO
                        K = K + 24

C.........................  Double check episode bounds
                        IF( K > NSTEPS ) EXIT

C.........................  If found data, go on to next hour                        
                        IF( METDAYS( K ) > 0 ) THEN
                            FND_DATA = .TRUE.
                            EXIT
                        END IF
                    END DO   ! end loop over days in episode
                
                    IF( FND_DATA ) CYCLE

C.....................  Still no data, exit with error                
                    MESG = 'Meteorology data does not cover ' //
     &                     'requested episode.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END DO   ! end loop over hours in day
        END IF   ! end episode averaging check

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

C.........  Make sure min and max fall within range allowed by MOBILE6
        IF( MINTEMP < M6MINTEMP ) THEN
            MINTEMP = M6MINTEMP
            WRITE( MESG,94070 ) 
     &          'WARNING: Adjusting minimum temperature to ', 
     &          M6MINTEMP, ' degrees F.'
            CALL M3MSG2( MESG )
        END IF
            
        IF( MAXTEMP > M6MAXTEMP ) THEN
            MAXTEMP = M6MAXTEMP
            WRITE( MESG,94070 ) 
     &          'WARNING: Adjusting maximum temperature to ',
     &          M6MAXTEMP, ' degrees F.'
            CALL M3MSG2( MESG )
        END IF

C.........  Convert temperature parameters to proper units
C.........  Kelvin for now - future can be dependant on met input units
        MINTEMP = FTOC * ( MINTEMP - 32. ) + CTOK
        MAXTEMP = FTOC * ( MAXTEMP - 32. ) + CTOK

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

C.........  Compare the ungridding matrix with previous settings
C.........  The subgrid parameters will be set here, if there is a subgrid
        CALL CHKGRID( UNAME, 'GMAT', 1, EFLAG )
        
        IF( EFLAG ) THEN
            MESG = 'Ungridding matrix is inconsistent with ' //
     &             'meteorology data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF       

C.........  Allocate memory for other arrays in the program
        ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UMAT', PROGNAME )
        
        ALLOCATE( TA( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( QV( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QV', PROGNAME )
        ALLOCATE( PRES( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRES', PROGNAME )
        
        ALLOCATE( DAYBEGT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        ALLOCATE( DAYENDT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        ALLOCATE( LDAYSAV( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDAYSAV', PROGNAME )
        
        ALLOCATE( TASRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASRC', PROGNAME )
        ALLOCATE( QVSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVSRC', PROGNAME )
        ALLOCATE( PRESSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRESSRC', PROGNAME )
        
        ALLOCATE( NDAYSRC( NSRC,24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NDAYSRC', PROGNAME )

        NDAYSRC = 0   ! array

C.........  Create array of which sources are affected by daylight savings
        CALL GETDYSAV( NSRC, IFIP, LDAYSAV )

C.........  Read ungridding matrix 
        CALL RDUMAT( UNAME, NSRC, NMATX, NMATX, 
     &               UMAT( 1 ), UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ) )

C.........  Allocate memory for storing meteorology profiles
        ALLOCATE( TKHOUR( NSRC, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKHOUR', PROGNAME )
        ALLOCATE( QVHOUR( NSRC, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVHOUR', PROGNAME )
        ALLOCATE( BPHOUR( NSRC, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BPHOUR', PROGNAME )
        
        ALLOCATE( TDYCNTY( NDYCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TDYCNTY', PROGNAME )
        ALLOCATE( QVDYCNTY( NDYCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVDYCNTY', PROGNAME )
        ALLOCATE( BPDYCNTY( NDYCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BPDYCNTY', PROGNAME )
        
        ALLOCATE( TWKCNTY( NWKCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TWKCNTY', PROGNAME )
        ALLOCATE( QVWKCNTY( NWKCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVWKCNTY', PROGNAME )
        ALLOCATE( BPWKCNTY( NWKCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BPWKCNTY', PROGNAME )
        
        ALLOCATE( TMNCNTY( NMNCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMNCNTY', PROGNAME )
        ALLOCATE( QVMNCNTY( NMNCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVMNCNTY', PROGNAME )
        ALLOCATE( BPMNCNTY( NMNCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BPMNCNTY', PROGNAME )
        
        ALLOCATE( TEPCNTY( NEPCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TEPCNTY', PROGNAME )
        ALLOCATE( QVEPCNTY( NEPCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVEPCNTY', PROGNAME )
        ALLOCATE( BPEPCNTY( NEPCNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BPEPCNTY', PROGNAME )

C.........  Get output directory information from the environment
        MESG = 'Path where hourly meteorology files will be written'
        CALL ENVSTR( 'SMK_METPATH', MESG, '.', TEMPDIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: Meteorology files being placed in ' //
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_METPATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF
        
C.........  Open output files
        IF( DAYAVER ) THEN
            DNAME = 'DAYTEMP'
            CALL OPENSHOUR( ENAME, 'daily', SDATE, EDATE, TVARNAME, 
     &                      NDYCNTY, TEMPDIR, DNAME )
            DAYOPEN = .TRUE.
        END IF
        IF( WEEKAVER ) THEN
            WNAME = 'WEEKTEMP'
            CALL OPENSHOUR( ENAME, 'weekly', SDATE, EDATE, TVARNAME,
     &                      NWKCNTY, TEMPDIR, WNAME )
            WEEKOPEN = .TRUE.
        END IF
        IF( MONAVER ) THEN
            MNAME = 'MONTEMP'
            CALL OPENSHOUR( ENAME, 'monthly', SDATE, EDATE, TVARNAME,
     &                      NMNCNTY, TEMPDIR, MNAME )
            MONOPEN = .TRUE.
        END IF
        IF( EPIAVER ) THEN
            PNAME = 'EPITEMP'
            CALL OPENSHOUR( ENAME, 'episode', SDATE, EDATE, TVARNAME,
     &                      NEPCNTY, TEMPDIR, PNAME )
        END IF

C.........  Process meteorology data...
        MESG = 'Processing meteorology data using variables ' //
     &         TRIM( TVARNAME ) // ', ' // TRIM( MIXNAME ) // 
     &         ', ' // TRIM( PRESNAME ) // '...'
        CALL M3MSG2( MESG )

C.........  Loop through days/hours of meteorology files
        DDATE = SDATE
        OTIME = 0
        WDATE = SDATE
        MDATE = SDATE
        
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        
        DO T = 1, NSTEPS

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

C.............  Determine input file for this hour
            POS = T
            
C.............  Get file number for current iteration
            FILENUM = METDAYS( POS )
            IF( FILENUM <= 0 ) THEN
                ALT_DATA = .TRUE.
            ELSE
                ALT_DATA = .FALSE.
            END IF

C.............  Skip file opening when not doing day averaging and using alternate data
            IF( .NOT. ALT_DATA .OR. DAYAVER ) THEN
            
C.................  Get file name
                METFILE = METLIST( ABS( FILENUM ) )

C.................  Close previous file if needed
                IF( METFILE .NE. PREVFILE ) THEN
                    IF( FILEOPEN ) THEN
                        IF( .NOT. CLOSE3( METNAME ) ) THEN
                            MESG = 'Could not close meteorology ' //
     &                             'file ' // TRIM( PREVFILE )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ELSE
                            FILEOPEN = .FALSE.
                        END IF
                    END IF

                    PREVFILE = METFILE

                END IF

C.................  Set logical file name
                IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                    MESG = 'Could not set logical file name for ' //
     &                     'file ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Open the meteorology data file
                IF ( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                    MESG = 'Could not open meteorology file ' // 
     &                     TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                ELSE
                    FILEOPEN = .TRUE.
                END IF

C.................  Reset read date when using alternate data
                IF( ALT_DATA ) THEN
                    IF( .NOT. DESC3( METNAME ) ) THEN
                        MESG = 'Could not get description of ' //
     &                         'meteorology file ' // TRIM( METFILE )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    ELSE
                        IF( SDATE3D < JDATE ) THEN
                            RDATE = JDATE - 1
                        ELSE
                            RDATE = JDATE + 1
                        END IF
                    END IF
                ELSE
                    RDATE = JDATE
                END IF
            
C.................  Read current meteorology file
                IF ( .NOT. READ3( METNAME, TVARNAME, 1, 
     &                            RDATE, JTIME, TA ) ) THEN
                    MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                     ' from ' // TRIM( METFILE ) 
                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                END IF
                
                IF ( .NOT. READ3( METNAME, MIXNAME, 1,
     &                            RDATE, JTIME, QV ) ) THEN
                    MESG = 'Could not read ' // TRIM( MIXNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                END IF
                
                IF ( .NOT. READ3( METNAME, PRESNAME, 1,
     &                            RDATE, JTIME, PRES ) ) THEN
                    MESG = 'Could not read ' // TRIM( PRESNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                END IF

C.................  Apply ungridding matrix 
                CALL APPLUMAT( NSRC, NMATX, TA, UMAT(1), UMAT( NSRC+1 ),
     &                         UMAT( NSRC+NMATX+1 ), TASRC )
                CALL APPLUMAT( NSRC, NMATX, QV, UMAT(1), UMAT( NSRC+1 ),
     &                         UMAT( NSRC+NMATX+1 ), QVSRC )
                CALL APPLUMAT( NSRC, NMATX, PRES, UMAT(1), 
     &                         UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ), 
     &                         PRESSRC )

C.................  Create hourly meteorology arrays by source
                CALL HOURTEMP( NSRC, JTIME, DAYBEGT, COUNTYSRC, 
     &                         MINTEMP, MAXTEMP, ALT_DATA, NDAYSRC )
     
            END IF  ! check for using alternate data or day averaging

C.............  Make sure we've waited long enough to catch all time zones
            IF( POS > TSPREAD ) THEN

C.................  Adjust time step for 24-hour arrays
                ARRAYPOS = MOD( POS - TSPREAD, 24 )
                IF( ARRAYPOS == 0 ) ARRAYPOS = 24

C.................  Process daily averages
                IF( DAYAVER ) THEN

C.....................  Average temperatures across county group                
                    CALL AVERTEMP( NSRC, NDYCNTY, DYCODES, 
     &                             COUNTYSRC( :,1 ), ARRAYPOS, 
     &                             TDYCNTY, QVDYCNTY, BPDYCNTY, 
     &                             NDAYSRC )

C.....................  Open new file if necessary (1 day per output file)
                    IF( .NOT. DAYOPEN ) THEN
                        CALL OPENSHOUR( ENAME, 'daily', DDATE, EDATE, 
     &                                  TVARNAME, NDYCNTY, TEMPDIR, 
     &                                  DNAME )
                        DAYOPEN = .TRUE.
                    END IF  

C.....................  Write time step to file
                    CALL WRSHOUR( DNAME, DDATE, OTIME, NDYCNTY, 
     &                            DYCODES, TDYCNTY, QVDYCNTY, BPDYCNTY )
     
                    IF( OTIME == 230000 ) THEN

C.........................  Close current output file                        
                        IF( .NOT. CLOSE3( DNAME ) ) THEN
                            MESG = 'Could not close file ' // 
     &                              DNAME( 1:LEN_TRIM( DNAME ) )
                            CALL M3EXIT( PROGNAME, DDATE, 0, MESG, 2 )
                        ELSE
                            DAYOPEN = .FALSE.
                        END IF
                    END IF     
                END IF

C.............  If current output day is Saturday, process weekly averages
                IF( WEEKAVER .AND. WKDAY( DDATE ) == 6 ) THEN

C.....................  Average temperatures across county group 
                    CALL AVERTEMP( NSRC, NWKCNTY, WKCODES, 
     &                             COUNTYSRC( :,1 ), ARRAYPOS,
     &                             TWKCNTY, QVWKCNTY, BPWKCNTY, 
     &                             NDAYSRC )

C.....................  Open new file if necessary (one week per output file)
                    IF( .NOT. WEEKOPEN ) THEN
                        CALL OPENSHOUR( ENAME, 'weekly', WDATE, EDATE, 
     &                                  TVARNAME, NWKCNTY, TEMPDIR, 
     &                                  WNAME )
                        WEEKOPEN = .TRUE.
                    END IF       

C.....................  Write time step to file
                    CALL WRSHOUR( WNAME, WDATE, OTIME, NWKCNTY, 
     &                            WKCODES, TWKCNTY, QVWKCNTY, BPWKCNTY )
                    
                    IF( OTIME == 230000 ) THEN

C.........................  Store current date to use for next week
                        WDATE = DDATE + 1
                        DUMMYTIME = 0
                        CALL NEXTIME( WDATE, DUMMYTIME, 0 )

C.........................  Close current output file                        
                        IF( .NOT. CLOSE3( WNAME ) ) THEN
                            MESG = 'Could not close file ' // 
     &                              WNAME( 1:LEN_TRIM( WNAME ) )
                            CALL M3EXIT( PROGNAME, WDATE, 0, MESG, 2 )
                        ELSE
                            WEEKOPEN = .FALSE.
                        END IF
                    END IF                     
                END IF
                
C.................  If last day of month, process monthly averages
                CALL DAYMON( DDATE, MONTH, DAY )
                CALL DAYMON( DDATE + 1, TMPMNTH, DAY )
            
                IF( MONAVER .AND. TMPMNTH /= MONTH ) THEN

C.....................  Average temperatures across county group 
                    CALL AVERTEMP( NSRC, NMNCNTY, MNCODES, 
     &                             COUNTYSRC( :,1 ), ARRAYPOS,
     &                             TMNCNTY, QVMNCNTY, BPMNCNTY,
     &                             NDAYSRC )

C.....................  Open new file if necessary (one month per output file)
                    IF( .NOT. MONOPEN ) THEN                        
                        CALL OPENSHOUR( ENAME, 'monthly', MDATE, EDATE,
     &                                  TVARNAME, NMNCNTY, TEMPDIR, 
     &                                  MNAME )
                        MONOPEN = .TRUE.
                    END IF       

C.....................  Write time step to file                     
                    CALL WRSHOUR( MNAME, MDATE, OTIME, NMNCNTY, 
     &                            MNCODES, TMNCNTY, QVMNCNTY, BPMNCNTY )               
                    
                    IF( OTIME == 230000 ) THEN

C.........................  Store current date to use for next month
                        MDATE = DDATE + 1
                        DUMMYTIME = 0
                        CALL NEXTIME( MDATE, DUMMYTIME, 0 )

C.........................  Close current output file   
                        IF( .NOT. CLOSE3( MNAME ) ) THEN
                            MESG = 'Could not close file ' // 
     &                              MNAME( 1:LEN_TRIM( MNAME ) )
                            CALL M3EXIT( PROGNAME, MDATE, 0, MESG, 2 )
                        ELSE
                            MONOPEN = .FALSE.
                        END IF
                    END IF
                END IF

            END IF    ! time zone check

C.............  Final steps on last time through
            IF( T == NSTEPS ) THEN

C.................  Process remaining week temperatures only if current date is 
C                   not Saturday; otherwise, this has already been handled above
                IF( WEEKAVER .AND. WKDAY( DDATE ) /= 6 ) THEN
                    IF( .NOT. WEEKOPEN ) THEN                    
                        CALL OPENSHOUR( ENAME, 'weekly', WDATE, EDATE, 
     &                                  TVARNAME, NWKCNTY, TEMPDIR, 
     &                                  WNAME )
                        WEEKOPEN = .TRUE.
                    END IF
     
                    DO K = 1, 24                    
                        CALL AVERTEMP( NSRC, NWKCNTY, WKCODES, 
     &                                 COUNTYSRC( :,1 ), K,
     &                                 TWKCNTY, QVWKCNTY, BPWKCNTY, 
     &                                 NDAYSRC )
                       
                        CALL WRSHOUR( WNAME, WDATE, ( K-1 )*10000, 
     &                                NWKCNTY, WKCODES, TWKCNTY, 
     &                                QVWKCNTY, BPWKCNTY )
                    END DO
                END IF

C.................  Process remaining month temperatures only if current date is
C                   not the last day of the month
                IF( MONAVER .AND. TMPMNTH == MONTH ) THEN
                    IF( .NOT. MONOPEN ) THEN
                        CALL OPENSHOUR( ENAME, 'monthly', MDATE, EDATE,
     &                                  TVARNAME, NMNCNTY, TEMPDIR, 
     &                                  MNAME )
                        MONOPEN = .TRUE.
                    END IF

                    DO K = 1, 24
                        CALL AVERTEMP( NSRC, NMNCNTY, MNCODES, 
     &                                 COUNTYSRC( :,1 ), K,
     &                                 TMNCNTY, QVMNCNTY, BPMNCNTY, 
     &                                 NDAYSRC )
                        
                        CALL WRSHOUR( MNAME, MDATE, ( K-1 )*10000, 
     &                                NMNCNTY, MNCODES, TMNCNTY, 
     &                                QVMNCNTY, BPMNCNTY )  
                    END DO
                END IF

C.................  Output episode averaged temperatures
                IF( EPIAVER ) THEN
                    DO K = 1, 24
                        CALL AVERTEMP( NSRC, NEPCNTY, EPCODES, 
     &                                 COUNTYSRC( :,1 ), K,
     &                                 TEPCNTY, QVEPCNTY, BPEPCNTY, 
     &                                 NDAYSRC )
     
                        CALL WRSHOUR( PNAME, SDATE, ( K-1 )*10000, 
     &                                NEPCNTY, EPCODES, TEPCNTY,
     &                                QVEPCNTY, BPEPCNTY )
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
                CALL NEXTIME( DDATE, OTIME, 10000 )
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

94000   FORMAT( A )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( A, I5 )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )
     
94070   FORMAT( A, F5.1, A )
                
        END PROGRAM PREMOBL
