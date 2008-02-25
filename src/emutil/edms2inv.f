      PROGRAM EDMS2INV

C***********************************************************************
C  program body starts at line  
C
C  DESCRIPTION:
C       This program converts the output from the EDMS model developed to
C       assess the air quality of airport emission sources into hourly 
C       and annual inventories that can be read by SMOKE.
C
C  PRECONDITIONS REQUIRED:
C       EDMS model has been run and provided airport sources emissions
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 11/2006 by B.H. Baek
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
C.......  MODULES
      USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, CNTYTZNM

      USE MODLISTS, ONLY: ITMSPC, ITCASA, ITVTSA, ITEXPL, ITNAMA,
     &                    NINVTBL, MXIDAT, INVDNAM


      IMPLICIT NONE

C.......  INCLUDES
      INCLUDE 'IODECL3.EXT'               ! I/O API function declarations
      INCLUDE 'EMCNST3.EXT'               ! emissions constant parameters
      
C.......  EXTERNAL FUNCTIONS and their descriptions:
      CHARACTER(2)  CRLF
      INTEGER       ENVINT
      INTEGER       FIND1
      INTEGER       FINDC
      INTEGER       GETFLINE
      INTEGER       GETNLIST
      INTEGER       INDEX1
      INTEGER       JULIAN
      INTEGER       PROMPTFFILE
      INTEGER       STR2INT
      REAL          STR2REAL
      LOGICAL       STRLIST
      LOGICAL       BLKORCMT
      LOGICAL       ENVYN
      
      EXTERNAL      CRLF, ENVINT, FIND1, INDEX1, GETFLINE, PROMPTFFILE, 
     &              STR2INT, STR2REAL, STRLIST, BLKORCMT, ENVYN, JULIAN,
     &              GETNLIST

C.......  LOCAL PARAMETERS
      CHARACTER(50), PARAMETER :: 
     & CVSW = '$Name$'  ! CVS release tag
      INTEGER, PARAMETER :: MXSEG = 15    ! number of segments in line
      INTEGER, PARAMETER :: MXOUT = 11    ! number of output values

C.......  LOCAL VARIABLES

C.......  Static arrays
      CHARACTER(150)   SEGMENT( MXSEG )    ! parsed input line

C.......  Allocatable arrays
      CHARACTER(15), ALLOCATABLE :: APRTID( : )
      CHARACTER(15), ALLOCATABLE :: APTLOC( : )
      CHARACTER(10), ALLOCATABLE :: APTSCC( : )
      CHARACTER(10), ALLOCATABLE :: TSCC  ( : )
      INTEGER,       ALLOCATABLE :: LOCID ( : )
      INTEGER,       ALLOCATABLE :: STATE ( : )
      INTEGER,       ALLOCATABLE :: COUNTY( : )
      INTEGER,       ALLOCATABLE :: POINTS( : )
      REAL,          ALLOCATABLE :: LAT   ( : )
      REAL,          ALLOCATABLE :: LON   ( : )
      REAL,          ALLOCATABLE :: AREA  ( : )
      REAL,          ALLOCATABLE :: X     ( : )
      REAL,          ALLOCATABLE :: Y     ( : )
      INTEGER,       ALLOCATABLE :: HEIGHT( : )
      INTEGER,       ALLOCATABLE :: ACRFT_NOX( : )
      INTEGER,       ALLOCATABLE :: ACRFT_PMT( : )
      CHARACTER(15), ALLOCATABLE :: OUTVAR( : )
      CHARACTER(8),  ALLOCATABLE :: TYPE  ( : )
      CHARACTER(8),  ALLOCATABLE :: CDATE ( : )   ! current DATE

      REAL,          ALLOCATABLE :: EFTOTAL( : )
      CHARACTER(15), ALLOCATABLE :: NPTOTAL( : )

      REAL,          ALLOCATABLE :: CVFACT( : )
      CHARACTER(15), ALLOCATABLE :: NPFACT( : )

      REAL,          ALLOCATABLE :: ALLVAL( :,: )    ! output variable values
      REAL,          ALLOCATABLE :: NONHAP( :,:,: )    ! output NONHAPTOG variable values

C.......  File units and logical names
      INTEGER      :: ADEV = 0            ! output annual inventory
      INTEGER      :: CDEV = 0            ! country, state, and county file
      INTEGER      :: DDEV = 0            ! output daily inventory
      INTEGER      :: EDEV = 0            ! list of main EDMS input files
      INTEGER      :: FDEV = 0            ! input file
      INTEGER      :: KDEV = 0            ! EDMS source id and SCC
      INTEGER      :: HDEV = 0            ! list of hourly emission input files
      INTEGER      :: LDEV = 0            ! unit number for log file
      INTEGER      :: RDEV = 0            ! report file
      INTEGER      :: SDEV = 0            ! ICAO airport code, FIPS, LAT and LONG
      INTEGER      :: PDEV = 0            ! unit number for inventory data table
      INTEGER      :: TDEV = 0            ! EDMS species conversion factors

C.......  Other local variables
      INTEGER         I, II, J, K, L, L1, M, NN, NF, T ! counters and indices
      INTEGER         IOS                 ! i/o status
      INTEGER         IREC                ! line counter
      INTEGER         IHOUR               ! current time step
      INTEGER         IDATE               ! current date
      INTEGER         NDATE               ! total processing days
      INTEGER         NDY                 ! number of processing days
      INTEGER         IDAY                ! current day
      INTEGER         IMON                ! current month
      INTEGER         IYR                 ! current year
      INTEGER         SDATE               ! start date
      INTEGER         STDAY               ! start day
      INTEGER         STMON               ! start month
      INTEGER         STYR                ! start year
      INTEGER         EDATE               ! start date
      INTEGER         ENDAY               ! end day
      INTEGER         ENMON               ! end month
      INTEGER         ENYR                ! end year
      INTEGER         STID                ! state code
      INTEGER         CYID                ! county code
      INTEGER         LHAP                ! HAPs identifier
      INTEGER         LOCVAL              ! location identifier
      INTEGER         FIP                 ! state and county code
      INTEGER         NOUTVAR             ! number of output variables
      INTEGER         MXSRC               ! Max number of SRCPARAM in input file
      INTEGER         NAPT                ! A number of EDMS_SCCs list
      INTEGER         NFAC                ! A number of conversion factors
      INTEGER         MXLOC               ! Max number of LOCATION in input file
      INTEGER         MXFILES             ! Max number of input files listed in FILELIST
      INTEGER         UTMZONE             ! airport UTM zone
      INTEGER         SEGS                ! number of x, y segments
      INTEGER         NSEGS               ! total number of x, y segments
      INTEGER         MXSEGS              ! max number of x, y segments
      INTEGER         NRPT                ! number of repeat
      INTEGER         POS                 ! position in INVDNAM array

      REAL            DAYTOT              ! Daily total emission
      REAL            LATICAO             ! ICAO latitude
      REAL            LONICAO             ! ICAO longitude
      REAL         :: LATORG  = 9999.0    ! EDMS origin latitude
      REAL         :: LONORG  = 9999.0    ! EDMS origin longitude
      REAL            XORG                ! EDMS origin UTM x-coor
      REAL            YORG                ! EDMS origin UTM y-coor
      REAL            XLOC                ! relative UTM
      REAL            YLOC                ! relative UTM
      REAL            XVAL                ! absolute UTM
      REAL            YVAL                ! absolute UTM
      REAL            LATVAL              ! absolute latitude
      REAL            LONVAL              ! absolute longitude
      REAL            ZLOC                ! elevated height
      REAL            WIDTH               ! source width
      REAL            LENGTH              ! source length
      REAL            SUM                 ! sum to determine area
      REAL            CONVF               ! HAPs conversion factors
      REAL            TOTAL               ! emission total of period

      LOGICAL         ADD                 ! true: add current airport source to master list
      LOGICAL      :: NO_FLAG = .FALSE.   ! true: flag for processing NO
      LOGICAL      :: NO2_FLAG= .FALSE.   ! true: flag for processing NO2
      LOGICAL      :: PM_FLAG = .FALSE.   ! true: flag for processing PMTOTAL
      LOGICAL         NPFLAG              ! true: flag for new pollutants
      LOGICAL      :: AFLAG = .FALSE.     ! true: read AIRPORT ICAO aiprot code from MAIN_EDMS file
      LOGICAL      :: OFLAG = .FALSE.     ! true: read ORIGIN cord from MAIN_EDMS file
      LOGICAL      :: COUNT = .TRUE.      ! true: count total nubmer of sources
      LOGICAL      :: SIZES = .FALSE.     ! true: define total nubmer of sources
      LOGICAL      :: FIRSTIME = .TRUE.   ! true: write YEAR into ptday output file

      CHARACTER(15)   APRT                ! airport identifier
      CHARACTER(10)   SCC                 ! SCC
      CHARACTER(8)    TDATE               ! Target DATE
      CHARACTER(8)    STDATE              ! Start DATE
      CHARACTER(8)    ENDATE              ! Ending DATE
      CHARACTER(8)    PDATE               ! previous DATE
      CHARACTER(4)    YEAR                ! inventory year
      CHARACTER(2)    CMONTH              ! data start month
      CHARACTER(8)    DEFDATE             ! default date for fake daily entries
      CHARACTER(4)    ICAOCODE            ! Airport ICAO code
      CHARACTER(5)    FIPSID              ! Airport FIPS code
      CHARACTER(3)    TZONE               ! time zone
      CHARACTER(8)    SRCTYPE             ! tmp source type

      CHARACTER(256)  LINE                ! input line buffer
      CHARACTER(256)  MESG                ! message buffer

      CHARACTER(IOVLEN3)  POLNAM          ! tmp pollutant name
      CHARACTER(IOVLEN3)  LPOLNAM         ! tmp previous pollutant name
      
      CHARACTER(16) :: PROGNAME = 'EDMS2INV' ! program name

C***********************************************************************
C   begin body of program EDMS2INV

      LDEV = INIT3()

C.......  Write copyright, version, web address, header info, and prompt
C         to continue running the program.
      CALL INITEM( LDEV, CVSW, PROGNAME )

C.......  Get environment variable settings
      CALL ENVSTR( 'START_DATE', MESG, ' ', STDATE, IOS )
      MESG ='Processing start date ' // STDATE
      CALL M3MSG2( MESG )
      STYR  = STR2INT( STDATE( 7:8 ) )
      STMON = STR2INT( STDATE( 1:2 ) )
      STDAY = STR2INT( STDATE( 4:5 ) )
      SDATE = JULIAN( STYR, STMON, STDAY )   ! convert to julian date

      CALL ENVSTR( 'END_DATE', MESG, ' ', ENDATE, IOS )
      MESG ='Processing ending date ' // ENDATE
      CALL M3MSG2( MESG )
      ENYR  = STR2INT( ENDATE( 7:8 ) )
      ENMON = STR2INT( ENDATE( 1:2 ) )
      ENDAY = STR2INT( ENDATE( 4:5 ) )
      EDATE = JULIAN( ENYR, ENMON, ENDAY )
      
      NDATE = EDATE - SDATE + 1    ! total processing days

C......  Open unit numbers of input files
      MESG = 'Enter logical name for INVENTORY DATA TABLE file'
      PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE', PROGNAME )

      MESG = 'Enter logical name for a main EDMS input file'
      EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MAIN_EDMS', PROGNAME )

      MESG = 'Enter logical name for hourly EDMS inputs list'
      HDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FILELIST', PROGNAME )
      
      MESG = 'Enter logical name for ICAO airport name, country, ' //
     &       'state, county, latitude, and longitude'
      SDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'ICAO_FIPS', PROGNAME )

      MESG = 'Enter logical name for a file of EDMS source and SCC'
      KDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'EDMS_SCC', PROGNAME )

      MESG = 'Enter logical name for a file of conversion factors'
      TDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'EDMS_FACT', PROGNAME )

      MESG = 'Enter logical name for country, state, and county ' //
     &       'file'
      CDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'COSTCY', PROGNAME )
      
      MESG = 'Enter logical name for output annual inventory'
      ADEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTINV', PROGNAME )
      
      MESG = 'Enter logical name for output hourly inventory'
      DDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTHOUR', PROGNAME )

      MESG = 'Enter logical name for report file'
      RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'REPORT', PROGNAME )

C.......  initialize array for storing pollutant names
      ALLOCATE( OUTVAR( 50 ), STAT=IOS )
      CALL CHECKMEM( IOS, 'OUTVAR', PROGNAME )

      OUTVAR = ' '
      FIPSID = ' '

C.......  Read the country, state, and county file
      CALL RDSTCY( CDEV, 1, 0 )

C........ Read inventory table
      CALL RDCODNAM( PDEV )

C.......  Read main EDMS input file
      DO I = 1, 2

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
          IF( SIZES ) THEN
      
              ALLOCATE( TSCC( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'TSCC', PROGNAME )
              ALLOCATE( APRTID( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'APRTID', PROGNAME )
              ALLOCATE( ACRFT_NOX( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'ACRFT_NOX', PROGNAME )
              ALLOCATE( ACRFT_PMT( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'ACRFT_PMT', PROGNAME )
              ALLOCATE( LOCID( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'LOCID', PROGNAME )
              ALLOCATE( STATE( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'STATE', PROGNAME )
              ALLOCATE( COUNTY( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'COUNTY', PROGNAME )
              ALLOCATE( LAT( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'LAT', PROGNAME )
              ALLOCATE( LON( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'LON', PROGNAME )
              ALLOCATE( HEIGHT( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'HEIGHT', PROGNAME )
              ALLOCATE( TYPE( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'TYPE', PROGNAME )
              ALLOCATE( AREA( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'AREA', PROGNAME )
              ALLOCATE( POINTS( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'POINTS', PROGNAME )
      
              TSCC   = ''   ! arrays
              APRTID = ''   ! arrays
              LOCID  = 0
              STATE  = 0
              COUNTY = 0
              LAT    = 0.
              LON    = 0.
              HEIGHT = 0
              TYPE   = ''
              AREA   = 0.
              POINTS = 0
              ACRFT_NOX = 0
              ACRFT_PMT = 0

          END IF

          IF( I == 2 ) REWIND( EDEV )

C...........  Initialize values before reading file
          IREC   = 0
          MXLOC  = 0
          LOCVAL = 1

C...........  Read through input file          
          DO
          
              READ( EDEV, 93000, IOSTAT=IOS ) LINE
              IREC = IREC + 1
              
              IF( IOS > 0 ) THEN
                  WRITE( MESG,94010 ) 'I/O error', IOS,
     &                'reading input file at line', IREC
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

C...............  Check for end of file              
              IF( IOS < 0 ) THEN
                  IREC = IREC - 1
                  EXIT
              END IF

C...............  Read data from LINE
              CALL PARSLINE( LINE, MXSEG, SEGMENT )

C...............  Look for airport ICAO code and location coordinate
              IF( INDEX( LINE, 'AIRPORT' ) > 0 .AND. COUNT ) THEN
                  ICAOCODE = ADJUSTL( SEGMENT( 3 ) )
                  L = LEN_TRIM( ICAOCODE )
                  IF( L .NE. 4 ) THEN
                      MESG = 'ERROR : Incorrect ICAO airport code ' //
     &                       'format or not listed in EDMS_MAIN file.'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

C...................  Look for airport FIPS, XLOC and YLOC using ICAO code
                  CALL READ_ICAO_FIPS( ICAOCODE )

C...................  stores EDMS source IDs and SCCs
                  CALL READ_EDMS_FACTOR

C...................  stores EDMS source IDs and SCCs
                  CALL READ_EDMS_SCC

C...................  Write airport ICAO code and assigned fips code
                  MESG = 'Current AIRPORT ICAO CODE : ' //  ICAOCODE // 
     &                   '  ::  ' // 'co/st/ct (FIPS) : ' // FIPSID
                  IF( .NOT. AFLAG ) CALL M3MSG2( MESG )

                  AFLAG = .TRUE.

                  CYCLE

              END IF

C...............  Look for airport ICAO code and location coordinate
              IF( INDEX( LINE, 'ORIGIN' ) > 0 .AND. COUNT ) THEN
                  LATORG = STR2REAL( SEGMENT( 3 ) )
                  LONORG = STR2REAL( SEGMENT( 4 ) )

                  WRITE( MESG, 94011) 'Current AIRPORT LAT & LON ' //
     &                   ' Coordinate : ', LATORG , ' X ', LONORG

                  IF( .NOT. OFLAG ) CALL M3MSG2( MESG )
                  
                  OFLAG = .TRUE.
                  CYCLE

              END IF

C...............  To ensure that it reads airport info from EDMS file first
              IF( FIPSID == ' ' .AND. SIZES ) THEN 
                  MESG = 'ERROR: Airport ICAO code is not found : ' //
     &               'Please make sure that you list a main EDMS file'//
     &               ' prior to hourly EDMS input file in FILELIST'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

              IF( ( LATORG == 9999. .OR. LONORG == 9999. )
     &              .AND. SIZES                         ) THEN 
                  MESG = 'CRITICAL WARNING: Airport **ORIGIN '//
     &               'coordinate was not found in main EDMS file :' //
     &               CRLF() // BLANK5 // 'Using default Lat and Lon ' //
     &               'from a ICAO_FIPS file.'
                  CALL M3MSG2( MESG )

                  LATORG = LATICAO  ! replace with ICAO latitude
                  LONORG = LONICAO  ! replace with ICAO longitude

              END IF
              
C...............  Get state and county code from airport source ID
              STID = STR2INT( FIPSID( 1:2 ) )
              CYID = STR2INT( FIPSID( 3:5 ) )

              IF( INDEX( LINE, 'LOCATION' ) > 0 ) THEN

C...................  To ensure that AIRPORT and ORIGIN info were existed in EDMS_MAIN file
                  IF( .NOT. AFLAG ) THEN
                      MESG = 'ERROR: Could not find the ICAO airport '//
     &                       'code (AIRPORT) from MAIN_EMDS file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  IF( .NOT. OFLAG ) THEN
                      MESG = 'ERROR: Could not find the ORIGIN ' //
     &                       'airport location from MAIN_EMDS file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  APRT = ICAOCODE // '-' // ADJUSTL( SEGMENT( 2 ) )

C...................  Define the location of airport id in main EDMS source list
                  K = INDEX1( APRT, NAPT, APTLOC )
                  IF( K < 1 ) THEN
                      MESG = 'ERROR: Could not find the matched EDMS '//
     &                   'source ID (' // TRIM( SEGMENT( 2 ) ) // 
     &                   ') from EDMS_SCC file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  SCC = APTSCC( K )

C...............  Save the source type
                  SRCTYPE = SEGMENT( 3 )

C...............  Convert rel. x,y to abs. lat and long
                  XLOC = STR2REAL( SEGMENT( 4 ) )
                  YLOC = STR2REAL( SEGMENT( 5 ) )
             
C...................  Convert rel x,y to absolute x,y
                  CALL LL2UTM( LONORG, LATORG, UTMZONE, XORG, YORG)
                  XVAL = XLOC + XORG
                  YVAL = YLOC + YORG

C...................  Convert abs. x,y to lat and lon
                  CALL UTM2LL( XVAL, YVAL, UTMZONE, LONVAL, LATVAL )

C...................  Check if we already have this airport in master list
                  ADD = .TRUE.
                  IF( SIZES ) THEN
                      DO J = 1, MXLOC
                          IF( APRT   == APRTID( J ) .AND.
     &                        LATVAL == LAT( J )    .AND.
     &                        LONVAL == LON( J )          ) THEN
                              ADD = .FALSE.
                              CYCLE
                          ELSE
                              LOCVAL = LOCID( J ) + 1
                          END IF
                      END DO
                  END IF
              
                  IF( ADD ) THEN
                      MXLOC = MXLOC + 1
                      IF( SIZES ) THEN
                          TSCC  ( MXLOC ) = SCC
                          APRTID( MXLOC ) = APRT
                          LOCID ( MXLOC ) = LOCVAL
                          STATE ( MXLOC ) = STID
                          COUNTY( MXLOC ) = CYID
                          LAT   ( MXLOC ) = LATVAL
                          LON   ( MXLOC ) = LONVAL
                          TYPE  ( MXLOC ) = SRCTYPE
                      END IF
                  END IF

              END IF

            IF( INDEX( LINE, 'SRCPARAM' ) > 0 .AND. .NOT. COUNT ) THEN
                APRT = ICAOCODE // '-' // ADJUSTL( SEGMENT( 2 ) )
                ZLOC = STR2REAL( SEGMENT ( 4 ) )

                K = INDEX1( APRT, MXSRC, APRTID )
                IF( K .LE. 0 ) THEN
                    MESG = 'ERROR: Could not find matching source '
     &                       // TRIM( SEGMENT( 2 ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                HEIGHT( K ) = ZLOC * 3.28084  ! convert m to ft

                IF( TYPE( K ) .EQ. 'AREA' ) THEN
                    WIDTH = STR2REAL( SEGMENT( 5 ) )
                    LENGTH = STR2REAL( SEGMENT( 6 ) )
                    AREA( K ) = WIDTH * LENGTH

                ELSE IF( TYPE( K ) .EQ. 'VOLUME' ) THEN
                    AREA( K ) = 1.0

                ELSE IF( TYPE( K ) .EQ. 'AREAPOLY' ) THEN
                    POINTS( K ) = STR2INT( SEGMENT( 5 ) )

                ENDIF

            END IF

            IF( INDEX( LINE, 'AREAVERT' ) > 0 .AND. .NOT. COUNT ) THEN
                APRT = ICAOCODE // '-' // ADJUSTL( SEGMENT( 2 ) )

                SEGS = GETNLIST( 256, LINE )
                SEGS = SEGS - 2
                IF( SEGS .EQ. 0 ) CYCLE

                K = INDEX1( APRT, MXSRC, APRTID )
                IF( K .LE. 0 ) THEN
                    MESG = 'ERROR: Could not find matching source '
     &                       // TRIM( SEGMENT( 2 ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Allocate memory
                ALLOCATE( X( POINTS( K ) ) )
                ALLOCATE( Y( POINTS( K ) ) )
                X = 0.
                Y = 0.

                IF( SEGS .EQ. ( POINTS( K ) * 2 ) ) THEN
                    M = 2
                    DO J = 1, POINTS( K )
                        M = M + 1
                        X( J ) = STR2REAL( SEGMENT( M ) )
                        M = M + 1
                        Y( J ) = STR2REAL( SEGMENT( M ) )
                    END DO

                ELSE IF( SEGS .LT. ( POINTS( K ) * 2 ) ) THEN
                    M = 2
                    DO J = 1, SEGS / 2
                        M = M + 1
                        X( J ) = STR2REAL( SEGMENT( M ) )
                        M = M + 1
                        Y( J ) = STR2REAL( SEGMENT( M ) )
                    END DO

                    NSEGS = SEGS
                    NRPT  = INT( POINTS( K ) / 5 )

                    DO II = 1,  NRPT 

                        READ( EDEV, 93000, IOSTAT=IOS ) LINE
                        IREC = IREC + 1

                        IF( IOS > 0 ) THEN
                            WRITE( MESG,94010 ) 'I/O error', IOS,
     &                          'reading input file at line', IREC
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF

C.........................  Check for end of file              
                        IF( IOS < 0 ) THEN
                            IREC = IREC - 1
                            EXIT
                        END IF

C.........................  Skip if no data available
                        SEGS = GETNLIST( 256, LINE )
                        SEGS = SEGS - 2
                        IF( SEGS .NE. 0 ) THEN

C.............................  Read data from LINE
                            CALL PARSLINE( LINE, MXSEG, SEGMENT )

                            M = 2
                            MXSEGS = ( NSEGS/2 ) + ( SEGS/2 )

                            DO J = ( NSEGS/2 ) + 1, MXSEGS
                                M = M + 1
                                X( J ) = STR2REAL( SEGMENT( M ) )
                                M = M + 1
                                Y( J ) = STR2REAL( SEGMENT( M ) )
                            END DO

                            NSEGS = NSEGS + SEGS

                        END IF

                    END DO   ! end of looping polygon source vertice coordinates

                END IF

                SUM = 0.
                DO J = 1, POINTS( K )
                    IF( J .LT. POINTS( K ) ) THEN
                        SUM = SUM + ( ( X( J ) * Y( J + 1 ) ) -
     &                              ( X( J + 1 ) * Y( J ) ) )
                    ELSE
                        SUM = SUM + ( ( X( J ) * Y( 1 ) ) -
     &                              ( X( 1 ) * Y( J ) ) )
                    END IF
                END DO

                AREA( K ) = ABS( 0.5 * SUM )

                DEALLOCATE( X, Y )

            END IF

          END DO

          COUNT = .FALSE.
          SIZES = .TRUE.
          MXSRC = MXLOC

      END DO

C.......  Look up time zone based on state and county
      FIP = STID * 1000 + CYID
      K = FIND1( FIP, NCOUNTY, CNTYCOD )
             
      IF( K < 1 ) THEN
          WRITE( MESG,94010 ) 'Could not find FIPS code',
     &        FIP, 'in COSTCY file'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END IF
              
      TZONE = CNTYTZNM( K )

C.......  Write hourly inventory header
      WRITE( DDEV,93000 ) '#TYPE    Hourly Airport EDMS Inventory'
      WRITE( DDEV,93000 ) '#COUNTRY US'

C.......  Reading a input list file
      MXFILES = GETFLINE( HDEV, 'FILELIST input file' )

C.......  Allocate arrays to store output values
      ALLOCATE( ALLVAL( MXSRC, 24 ), STAT=IOS )
      CALL CHECKMEM( IOS, 'ALLVAL', PROGNAME )
      ALLOCATE( NONHAP( NDATE, MXSRC, 24 ), STAT=IOS )
      CALL CHECKMEM( IOS, 'NONHAP', PROGNAME )
      ALLOCATE( CDATE( NDATE ), STAT=IOS )
      CALL CHECKMEM( IOS, 'CDATE', PROGNAME )
      ALLOCATE( NPTOTAL( MXFILES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'NPTOTAL', PROGNAME )
      ALLOCATE( EFTOTAL( MXFILES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'EFTOTAL', PROGNAME )

      ALLVAL = 0.  ! array
      NONHAP = 0.  ! array
      CDATE  = ' ' ! array
      EFTOTAL = 0. ! array
      NPTOTAL = ' ' ! array

C.......  Process files in input list

      I = 0
      LHAP = 0
      LPOLNAM = ' '
      POLNAM = ' '
      DO NN = 1, MXFILES

          READ( HDEV, 93000, IOSTAT=IOS ) LINE

C...........  Skip blank and comment lines
          IF( BLKORCMT( LINE ) ) CYCLE

C...........  Number of files
          I = I + 1

C...........  Define current processing pollutant
          CALL PARSLINE( LINE, 2, SEGMENT )
          POLNAM = TRIM( SEGMENT( 1 ) )
          LINE = ADJUSTL( SEGMENT( 2 ) )

C...........  One cond:1st pol has to be THC to reduce usage of memory
          IF( I == 1 ) THEN
              IF( POLNAM /= 'THC' ) THEN
                  MESG = 'ERROR: You MUST process THC first '
     &                    // 'before processing other pollutants.'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF
          END IF

          IF( POLNAM == 'VOC' ) THEN
              MESG = 'ERROR: Can not process VOC BUT can ' //
     &               'convert THC to VOC.'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C...........  One cond:1st pol has to be THC to reduce usage of memory
          ELSE IF( POLNAM == 'THC' ) THEN
              IF( I == 1 ) THEN
                  MESG = 'NOTE: Pollutant THC is internally ' //
     &                   'converted to TOG.'
                  CALL M3MSG2( MESG )
                  POLNAM = 'TOG'
              END IF
          END IF

C...........  Rename total aircraft PM species
          IF( POLNAM == 'PM25' ) POLNAM = 'PM2_5'

C...........  look for pol names in a list of pols
          POS  = INDEX1( POLNAM, NINVTBL, ITNAMA )

          IF( POS < 1 ) THEN
              MESG = 'ERROR: Pollutant ID ' // TRIM( POLNAM )//
     &               ' is not found in the INVTABLE file.'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

C...........  Check for i/o errors
          IF( IOS /= 0 ) THEN
              WRITE( MESG,94010 ) 'I/O error', IOS,
     &            'reading input file at line', I
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

C...........  Open input pollutant file      
          OPEN( FDEV, FILE = LINE, STATUS='OLD', IOSTAT=IOS )

          IF( IOS /= 0 ) THEN
              MESG = 'Could not open file:' // 
     &               CRLF() // BLANK5 // TRIM( LINE )
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

          MESG = 'Successfully opened EDMS hourly input file:' //
     &           CRLF() // BLANK5 // TRIM( LINE )
          CALL M3MSG2( MESG )

C...........  Initialize values before reading file
          NDY    = 1
          IREC   = 0
          NPFLAG = .TRUE.
          PDATE  = ' '
          TDATE  = ' '
          ALLVAL = 0.0
          TOTAL  = 0.0          

C...........  Read through input file          
          DO

              READ( FDEV, 93000, IOSTAT=IOS ) LINE
              IREC = IREC + 1

              IF( IOS > 0 ) THEN
                  WRITE( MESG,94010 ) 'I/O error', IOS,
     &                'reading input file at line', IREC
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

C...............  Check for end of file              
              IF( IOS < 0 ) THEN
                  IREC = IREC - 1
                  EXIT
              END IF

C...............  Skip blank/comment/non-HOUREMIS lines
              IF( INDEX( LINE, 'HOUREMIS' ) < 1 ) CYCLE

C...............  Read data from LINE
              CALL PARSLINE( LINE, MXSEG, SEGMENT )

C...............  Skip when it is out of range of episode dates
              IYR   = STR2INT( SEGMENT( 3 ) )    ! integer current year
              IMON  = STR2INT( SEGMENT( 4 ) )    ! integer current month
              IDAY  = STR2INT( SEGMENT( 5 ) )    ! integer current hour
              IDATE = JULIAN( IYR, IMON, IDAY )

C...............  Skip non-target episode date
              IF( IDATE < SDATE ) CYCLE
              IF( IDATE > EDATE ) GOTO 555

C...............  Ensure that it read pollutant name first
              IF( POLNAM == ' ' ) THEN
                  MESG = 'ERROR: Pollutant name (POLLUTANT) was NOT ' //
     &                   'defined in EDMS hourly files.'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

C...............  temporarily store source characteristics
              IHOUR = STR2INT( SEGMENT( 6 ) )    ! integer current hour
              APRT = ICAOCODE // '-' // ADJUSTL( SEGMENT( 7 ) )

C...............  Define the location of airport id in main EDMS source list
              K = INDEX1( APRT, MXSRC, APRTID )
              IF( K < 1 ) THEN
                  MESG='ERROR: Hourly EDMS source ID ( ' // 
     &               TRIM(SEGMENT( 7 )) // ' ) is not listed in the ' //
     &               'MAIN_EDMS file.'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

C...............  Store the date of hourly emissions
              TDATE = TRIM(SEGMENT( 4 )) // '/' // TRIM(SEGMENT( 5 )) //
     &                '/' // TRIM(SEGMENT( 3 ))

              IF( TDATE == ' ' ) THEN
                  MESG = 'INTERNAL ERROR: Incorrect format of HRE file'
              END IF

              IF( FIRSTIME ) THEN
                  YEAR = '20' // TRIM( SEGMENT( 3 ) )
                  WRITE( DDEV,93000 ) '#YEAR    ' // YEAR
                  FIRSTIME = .FALSE.
              END IF

C...............  Reinitialize date for new hourly pollutant file
              IF( NPFLAG ) THEN
                  PDATE = TDATE
                  NPFLAG = .FALSE.
              END IF

C...............  Write out all stored output values to PTHOUR hourly inv file
              IF( PDATE /= TDATE ) THEN

                  DO J = 1, MXSRC

C.......................  Writing current processing hour message
                      MESG = 'Processing hourly ' // TRIM( POLNAM ) // 
     &                       ' emission on date : '// PDATE
                      IF( J == 1 ) CALL M3MSG2( MESG )

                      DAYTOT = 0.0

C.......................  compute daily total by summing 24hrs hourly emis
                      DO T = 1, 24
                          DAYTOT = DAYTOT + ALLVAL( J,T )
                      END DO

C.......................  Skip any zero daily total
                      IF( DAYTOT == 0.0 ) CYCLE 

                      IF( POLNAM /= 'TOG' ) THEN 
                      IF( POLNAM /= 'PMTOTAL' ) THEN 
                         WRITE( DDEV,93020 ) STATE( J ), COUNTY( J ), 
     &                      APRTID( J ),LOCID( J ), HEIGHT( J ), 
     &                      POLNAM(1:5), PDATE, TZONE, 
     &                     ( ALLVAL(J,T), T = 1,24 ), DAYTOT, TSCC( J ),
     &                     POLNAM
                      END IF
                      END IF
                  END DO

                  CDATE( NDY ) = PDATE    ! Store current date 
                  NDY = NDY + 1   ! count number of proc days

                  PDATE = TDATE           ! Store previous date
                  ALLVAL = 0.0

              END IF

C...............  Store output values : convert metric g/sec/m2 to short tons/hr
              ALLVAL( K, IHOUR ) =  STR2REAL( SEGMENT( 8 ) ) 
     &                              * 0.003968254 * AREA( K )

C...............  Store output values for report: convert metric g/sec/m2 to kg/hr
              TOTAL = TOTAL + STR2REAL( SEGMENT( 8 ) ) * 3.6 * AREA(K)

C...............  Convert THC(g in CH4) to TOG(g) 
C                 Note: This conversion factor 1.148106 is based on
C                 GSPRO #1098 for SCC 227502000 only.
C...............  Store converted TOG (THC-->TOG) to NONHAP(:,:,:) initially
              IF( POLNAM == 'TOG' ) THEN

C................... 0.947 = mass conv factor (0.865) * 1.0947 (THC to VOC)
                  ALLVAL( K, IHOUR ) = ALLVAL( K, IHOUR ) * 0.947  !CH4 to C

C...............  Convert TOG(g in C) = 1.148106 * VOC(g in C) 
                  ALLVAL( K, IHOUR ) = ALLVAL( K, IHOUR ) * 1.148106

                  NONHAP( NDY, K, IHOUR ) = ALLVAL( K, IHOUR )

C............... Convert NO emission in NO2 equivalcy to NO equivalency (30/44)
              ELSE IF( POLNAM == 'NO' ) THEN
                  NO_FLAG = .TRUE.
                  IF( ALLVAL( K, IHOUR ) > 0 ) ACRFT_NOX( K ) = K
                  ALLVAL( K, IHOUR ) = ALLVAL( K, IHOUR ) * 0.652

C............... Convert NO2 emission in NO2 equivalcy to NO equivalency (30/44)
              ELSE IF( POLNAM == 'NO2' ) THEN
                  NO2_FLAG = .TRUE.
                  IF( ALLVAL( K, IHOUR ) > 0 ) ACRFT_NOX( K ) = K

C............... Check NO/NO2 emission has been processed for NOX renoralization
C............... which means excluding aircraft emissions from NOX HRE.
              ELSE IF( POLNAM == 'NOX' ) THEN
                  IF( .NOT. NO_FLAG .OR. .NOT. NO2_FLAG ) THEN
                      MESG = 'ERROR: Please process both NO and NO2 '//
     &                       'HRE files before processing NOX HRE file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  IF( K == ACRFT_NOX( K ) ) THEN
                      ALLVAL( K, IHOUR ) = 0.0
                      CYCLE   ! Skip if it is aircraft source
                  END IF

C............... Check PMTOTAL has been processed for PM10/PM25 renormalization
              ELSE IF( POLNAM == 'PMTOTAL' ) THEN
                  PM_FLAG = .TRUE.
                  IF( ALLVAL( K, IHOUR ) > 0 )  ACRFT_PMT( K ) = K
                  CYCLE         ! Skip processing PMTOTAL

C............... Check PM25 has been processed for PM10/PM25 renormalization
C............... which means excluding aircraft emissions from PM10/PM2_5 HRE.
              ELSE IF( POLNAM == 'PM10' .OR. POLNAM == 'PM2_5' ) THEN
                  IF( .NOT. PM_FLAG ) THEN
                      MESG = 'ERROR: Please process PMTOTAL HRE files'//
     &                  ' before processing NOX, PM10 or PM2_5 HRE file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  IF( K == ACRFT_PMT( K ) ) THEN
                      ALLVAL( K, IHOUR ) = 0.0
                      CYCLE   ! Skip if it is aircraft source
                  END IF


C............... Temporarly adjustment of PEC emission
              ELSE IF( POLNAM == 'PEC' ) THEN
                  ALLVAL( K, IHOUR ) = ALLVAL( K, IHOUR ) * 0.2623

C..............  Conversion all HAPs (g in CH4) to (g)
              ELSE IF( ITVTSA(POS)=='V' .OR. ITVTSA(POS)=='T' )THEN

C...................  Apply individual HAP conversion factor to convert
C                     CH4 equivalent to species mass equivalent.
                  NF = INDEX1( POLNAM, NFAC, NPFACT )
                  IF( NF < 1 ) THEN
                      MESG = 'ERROR: Could not find the matched EDMS '//
     &                   'conversion factor of ( ' // TRIM( POLNAM ) // 
     &                   ' ) from EDMS_FACTOR file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF
                  CONVF = CVFACT( NF )

                  ALLVAL( K,IHOUR ) = ALLVAL( K,IHOUR ) * CONVF

C...................  Count a number of HAPs
                  IF( POLNAM /= LPOLNAM ) THEN
                      LHAP = LHAP + 1   ! counts HAPs 
                      LPOLNAM = POLNAM
                  END IF

C...................  Substract all HAPs from CAP TOG (stored in NONHAP array)
C                     to compute NONHAPTOG
                  NONHAP( NDY,K,IHOUR ) = NONHAP( NDY,K,IHOUR )
     &                                    - ALLVAL( K, IHOUR )

C..................  Error msg when total HAPs is greater than TOG
                 IF( NONHAP( NDY,K,IHOUR ) < 0.0 ) THEN
                    WRITE( MESG, 94010 )'ERROR: Toal sum of HAPs is '//
     &                 'greater then CAP TOG emission on ' // PDATE //
     &                 ' at timestep ::', IHOUR 
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

              END IF

          END DO

C...........  Write out last date stored output values to PTHOUR hourly inv file
C             before processing next pollutant file.
555       DO J = 1, MXSRC

C...............  Writing current processing hour message
              MESG = 'Processing hourly ' // TRIM( POLNAM ) // 
     &               ' emission on date : '// TDATE
              IF( J == 1 ) CALL M3MSG2( MESG )

              DAYTOT = 0.0
C...............  compute daily total by summing 24hrs hourly emis
              DO K = 1, 24
                  DAYTOT = DAYTOT + ALLVAL( J,K )
              END DO

C...............  Skip any zero daily total
              IF( DAYTOT == 0.0 ) CYCLE
              IF( POLNAM /= 'TOG' ) THEN 
              IF( POLNAM /= 'PMTOTAL' ) THEN 

                WRITE( DDEV,93020 ) STATE( J ), COUNTY(J), APRTID(J),
     &            LOCID( J ),HEIGHT( J ), POLNAM, TDATE, TZONE,
     &            ( ALLVAL( J, K ), K = 1,24 ), DAYTOT, TSCC( J ),
     &            POLNAM
              END IF
              END IF

              CDATE( NDY ) = TDATE    ! Store current date 

          END DO

          CLOSE( FDEV )

C...........  Store total emission factor and EDMS pollutant names
          NPTOTAL( I ) = ADJUSTR( POLNAM )
          IF( POLNAM == 'TOG' ) NPTOTAL( I ) = '          ' //'  THC'
          EFTOTAL( I ) = TOTAL

C...........  Adding new pollutants
          ADD = .TRUE.

          DO J = 1, NOUTVAR
              IF( POLNAM == OUTVAR( J ) ) THEN
                  ADD = .FALSE.
                  CYCLE
              END IF
          END DO

          IF( ADD ) THEN
              IF( POLNAM /= 'TOG' ) THEN
              IF( POLNAM /= 'PMTOTAL' ) THEN
                  NOUTVAR = NOUTVAR + 1
                  OUTVAR( NOUTVAR ) = POLNAM
              END IF
              END IF
          END IF

      END DO

C.......  Total number of processing species      
      NN = I

C.......  Write out computed NONHAPTOG = TOG - total HAPs
      IF( LHAP > 0 ) THEN
          POLNAM = 'NONHAPTOG'
      ELSE
          POLNAM  = 'TOG'
      END IF

      NOUTVAR = NOUTVAR + 1
      OUTVAR( NOUTVAR ) = POLNAM

      DO I = 1, NDATE
          DO J = 1, MXSRC

C...............  Writing current processing hour message
              MESG = 'Processing hourly ' // TRIM( POLNAM ) // 
     &               ' emission on date : '// CDATE( I )
              IF( J == 1 .AND. LHAP > 0 ) CALL M3MSG2( MESG )

              DAYTOT = 0.0
C...............  compute daily total by summing 24hrs hourly emis
              DO K = 1, 24
                  DAYTOT = DAYTOT + ALLVAL( J,K )
              END DO

C...............  Skip any zero daily total
              IF( DAYTOT == 0.0 ) CYCLE

              WRITE( DDEV,93020 ) STATE(J), COUNTY(J), APRTID( J ),
     &            LOCID( J ),HEIGHT( J ), POLNAM, CDATE( I ), TZONE,
     &            ( NONHAP( I, J, K ), K = 1,24 ), DAYTOT, TSCC( J ),
     &            POLNAM

          END DO
      END DO

C...........  Writing current processing hour message
      WRITE( MESG,94010 ) 'NOTE : Total ', LHAP, ' HAPs were ' //
     &    'processed to compute ' // TRIM( POLNAM ) 
      IF( LHAP > 0 ) CALL M3MSG2( MESG )

C.......  Write ORL format annual inventory header for CAP/HAPs
      WRITE( ADEV,93000 ) '#ORL'
      WRITE( ADEV,93000 ) '#TYPE    Airport EDMS Inventory'
      WRITE( ADEV,93000 ) '#COUNTRY US'
      WRITE( ADEV,93000 ) '#YEAR    ' // YEAR

C.......  Write annual inventory values
      DO I = 1, MXSRC
           DO J = 1, NOUTVAR
               WRITE( ADEV,93010 ) STATE( I ), COUNTY( I ), APRTID( I ), 
     &             LOCID( I ), HEIGHT( I ), TSCC( I ), HEIGHT( I ),
     &             LON( I ), LAT( I ), UTMZONE, OUTVAR( J )
          END DO

      END DO

C...........  Writing reports for total emission EDMS pollutant
      WRITE( MESG,94010 ) 'Writing total emissions of EDMS pollutants'
      CALL M3MSG2( MESG )

C.......  Write ORL format annual inventory header for CAP/HAPs
      WRITE( RDEV,93000 ) '# Report of total emission factors for '//
     &                    'the FAA EDMS pollutants'
      WRITE( RDEV,93000 ) '# Period : ' // STDATE // ' - ' // ENDATE
      
C.......  Write total emission factors for EDMS pollutants
      WRITE( RDEV,93001 ) ( NPTOTAL( J ), J = 1, NN )      
      WRITE( RDEV,94012 ) ( EFTOTAL( J ), J = 1, NN )

C.......  End program successfully
      CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C.......  Formatted file I/O formats...... 93xxx
93000 FORMAT( A )
93001 FORMAT( 50(A,) )
93010 FORMAT( I2.2,I3.3,',"',A,'","',I5,'","',I10,'",,"","',A10,
     &        '",,,"',I10,'",,,,,,,,"L",',F9.4,',',F9.4,',',I2,
     &        ',"',A,'",,,,,,,,,,,,,,,,,,' )
93020 FORMAT(I2.2, I3.3, A15, 2I12,12X, A5, A8, A3, 24E7.1,E8.2,1X,A10,
     &       1X, A16)

C.......  Internal buffering formats...... 94xxx
94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )
94011 FORMAT( 10 ( A, :, F10.3, :, 2X  ) )
94012 FORMAT( 50F15.3 )

C******************  INTERNAL SUBPROGRAMS  *****************************

       CONTAINS

C      This subroutine looks for a co/st/ct (FIPS) code for current airport
C      from ICAO_FIPS file

          SUBROUTINE READ_ICAO_FIPS( ICAOCODE )

C...........  Subroutine arguments
          CHARACTER(*), INTENT ( IN ) :: ICAOCODE

C...........   Local variables
          INTEGER         I, N                  ! indices and counters

          INTEGER      :: NLINES = 0            ! number of lines in input file

          CHARACTER( 4 )  ICAOTMP               ! Read tmp ICAO code
          CHARACTER(256)  LINE                  ! Read buffer for a line
          CHARACTER(300)  MESG                  ! Message buffer
          CHARACTER(60)   SEGMENT( MXSEG )      ! line parsing array

C......................................................................

C.........  Get the number of lines
        NLINES = GETFLINE( SDEV, 'ICAO_FIPS input file' )
        
C.........  Read the ICAO_FIPS file and find FIPS for ICAOCODE

        N = 0
        IREC  = 0
        
        DO I = 1, NLINES

            READ ( SDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading ICAO_FIPS '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = TRIM( LINE )

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

            ICAOTMP = ADJUSTL( SEGMENT( 1 ) )    ! tmp ICAO code

            IF( ICAOTMP == ICAOCODE ) THEN
                N = N + 1
                FIPSID  =  ADJUSTL( SEGMENT( 2 ) )
                LATICAO = STR2REAL( SEGMENT( 3 ) )
                LONICAO = STR2REAL( SEGMENT( 4 ) )
                UTMZONE = STR2INT ( SEGMENT( 5 ) )
            END IF

        END DO    ! end of loop

        IF( N == 0 ) THEN 
            MESG = 'ERROR: Not available FIPS code for airport ' //
     &             'ICAO : ' // ICAOCODE
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        
        ELSE IF( N > 1 ) THEN
            MESG = 'WARNING: Duplicate ICAO airport entries in ' //
     &             'ICAO_FIPS file'
            CALL M3MSG2( MESG)

        END IF

        RETURN
        
C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000   FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

        END SUBROUTINE READ_ICAO_FIPS

C******************  INTERNAL SUBPROGRAMS  *****************************

C      This subroutine stores a list of EDMS CH4 equilvalent pollutant 
C      with conversion factors to convert to species mass emission factor unit

       SUBROUTINE READ_EDMS_FACTOR

C...........   Local variables
          INTEGER         I, N                  ! indices and counters

          INTEGER      :: NLINES = 0            ! number of lines in input file

          CHARACTER(256)  LINE                  ! Read buffer for a line
          CHARACTER(300)  MESG                  ! Message buffer
          CHARACTER(15 )  SEGMENT( 2 )          ! line parsing array

C......................................................................

C.........  Get the number of lines
        NLINES = GETFLINE( TDEV, 'EDMS_FACTOR input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
        ALLOCATE( CVFACT( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CVFACT', PROGNAME )
        ALLOCATE( NPFACT( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NPFACT', PROGNAME )

        NPFACT = ' '
        CVFACT = 0.0
        
C.........  Read the ICAO_FIPS file and find FIPS for ICAOCODE
        N = 0
        IREC  = 0

        DO I = 1, NLINES

            READ ( TDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading EDMS_FACTOR '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = TRIM( LINE )

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
            CALL PARSLINE( LINE, 2, SEGMENT )

            N = N + 1
            NPFACT( N ) = SEGMENT( 1 )
            CVFACT( N ) = STR2REAL( SEGMENT( 2 ) )
        END DO    ! end of loop

        IF( N == 0 ) THEN 
            MESG = 'ERROR: No entries of EDMS source with SCC'
            CALL M3MSG2( MESG )
        END IF

        NFAC = N

        RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000   FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

        END SUBROUTINE READ_EDMS_FACTOR

C******************  INTERNAL SUBPROGRAMS  *****************************

C      This subroutine stores a list of EDMS with SCC 

       SUBROUTINE READ_EDMS_SCC

C...........   Local variables
          INTEGER         I, N                  ! indices and counters

          INTEGER      :: NLINES = 0            ! number of lines in input file

          CHARACTER(15 )  PAPTLOC               ! previous airport src ID
          CHARACTER(15 )  CAPTLOC               ! current airport src ID
          CHARACTER(256)  LINE                  ! Read buffer for a line
          CHARACTER(300)  MESG                  ! Message buffer
          CHARACTER(60)   SEGMENT( MXSEG )      ! line parsing array

C......................................................................

C.........  Get the number of lines
        NLINES = GETFLINE( KDEV, 'EDMS_SCC input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
        ALLOCATE( APTLOC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'APTLOC', PROGNAME )
        ALLOCATE( APTSCC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'APTSCC', PROGNAME )

        APTLOC = ' '
        APTSCC = ' '
        
C.........  Read the ICAO_FIPS file and find FIPS for ICAOCODE

        N = 0
        IREC  = 0
        PAPTLOC = ' '

        DO I = 1, NLINES

            READ ( KDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading EDMS_SCC '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = TRIM( LINE )

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

            CAPTLOC = ICAOCODE // '-' // ADJUSTL( SEGMENT( 1 ) )

            IF( CAPTLOC .NE. PAPTLOC ) THEN
                N = N + 1
                APTLOC( N ) = CAPTLOC
                APTSCC( N ) = ADJUSTL( SEGMENT( 2 ) )
                PAPTLOC = CAPTLOC

            ELSE
                MESG = 'WARNING: Duplicate entries of airport source:'//
     &                 TRIM( SEGMENT( 1 ) )
                CALL M3MSG2( MESG )

            END IF

        END DO    ! end of loop

        IF( N == 0 ) THEN 
            MESG = 'ERROR: No entries of EDMS source with SCC'
            CALL M3MSG2( MESG )
        END IF

        NAPT = N

        RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000   FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

        END SUBROUTINE READ_EDMS_SCC

C......................................................................

      END PROGRAM EDMS2INV
