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
      INTEGER       INDEX1
      INTEGER       JULIAN
      INTEGER       PROMPTFFILE
      INTEGER       STR2INT
      REAL          STR2REAL
      LOGICAL       STRLIST
      LOGICAL       BLKORCMT
      LOGICAL       ENVYN
      
      EXTERNAL      CRLF, ENVINT, FIND1, INDEX1, GETFLINE, PROMPTFFILE, 
     &              STR2INT, STR2REAL, STRLIST, BLKORCMT, ENVYN, JULIAN

C.......  LOCAL PARAMETERS
      CHARACTER(50), PARAMETER :: 
     & CVSW = '$Name$'  ! CVS release tag
      INTEGER, PARAMETER :: MXSEG = 15    ! number of segments in line
      INTEGER, PARAMETER :: MXOUT = 11    ! number of output values

C.......  LOCAL VARIABLES

C.......  Static arrays
      CHARACTER(20)   SEGMENT( MXSEG )    ! parsed input line

C.......  Allocatable arrays
      CHARACTER(15), ALLOCATABLE :: APRTID( : )
      CHARACTER(15), ALLOCATABLE :: APTLOC( : )
      CHARACTER(10), ALLOCATABLE :: APTSCC( : )
      CHARACTER(10), ALLOCATABLE :: TSCC  ( : )
      INTEGER,       ALLOCATABLE :: LOCID ( : )
      INTEGER,       ALLOCATABLE :: STATE ( : )
      INTEGER,       ALLOCATABLE :: COUNTY( : )
      REAL,          ALLOCATABLE :: LAT   ( : )
      REAL,          ALLOCATABLE :: LON   ( : )
      INTEGER,       ALLOCATABLE :: HEIGHT( : )
      CHARACTER(5),  ALLOCATABLE :: OUTVAR( : )

      REAL,          ALLOCATABLE :: ALLVAL( :,: )    ! output variable values

C.......  File units and logical names
      INTEGER         ADEV                ! output annual inventory
      INTEGER         CDEV                ! country, state, and county file
      INTEGER         DDEV                ! output daily inventory
      INTEGER         EDEV                ! list of main EDMS input files
      INTEGER         FDEV                ! input file
      INTEGER         KDEV                ! EDMS source id and SCC
      INTEGER         HDEV                ! list of hourly emission input files
      INTEGER         LDEV                ! unit number for log file
      INTEGER         SDEV                ! ICAO airport code, FIPS, LAT and LONG

C.......  Other local variables
      INTEGER         I, J, K, L, L1, L2  ! counters and indices
      INTEGER         IOS                 ! i/o status
      INTEGER         IREC                ! line counter
      INTEGER         IHOUR               ! current time step
      INTEGER         IDATE               ! current date
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
      INTEGER         LOCVAL              ! location identifier
      INTEGER         FIP                 ! state and county code
      INTEGER         NOUTVAR             ! number of output variables
      INTEGER         MXSRC               ! Max number of SRCPARAM in input file
      INTEGER         NAPT                ! A number of EDMS_SCCs list
      INTEGER         MXLOC               ! Max number of LOCATION in input file
      INTEGER         MXFILES             ! Max number of input files listed in FILELIST
      INTEGER         UTMZONE             ! airport UTM zone

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

      LOGICAL         ADD                 ! true: add current airport source to master list
      LOGICAL         NPFLAG              ! true: flag for new pollutants
      LOGICAL      :: AFLAG = .FALSE.     ! true: read AIRPORT ICAO aiprot code from MAIN_EDMS file
      LOGICAL      :: OFLAG = .FALSE.     ! true: read ORIGIN cord from MAIN_EDMS file
      LOGICAL      :: COUNT = .TRUE.      ! true: count total nubmer of sources
      LOGICAL      :: SIZES = .FALSE.     ! true: define total nubmer of sources
      LOGICAL      :: FIRSTIME = .TRUE.   ! true: write YEAR into ptday output file

      CHARACTER(15)   APRT                ! airport identifier
      CHARACTER(10)   SCC                 ! SCC
      CHARACTER(8)    STDATE              ! Start DATE
      CHARACTER(8)    ENDATE              ! Ending DATE
      CHARACTER(8)    CDATE               ! current DATE
      CHARACTER(8)    PDATE               ! previous DATE
      CHARACTER(11)   CHOUR               ! current DATE//HOUR
      CHARACTER(11)   PCHOUR              ! previous DATE//HOUR
      CHARACTER(4)    YEAR                ! inventory year
      CHARACTER(2)    CMONTH              ! data start month
      CHARACTER(8)    DEFDATE             ! default date for fake daily entries
      CHARACTER(4)    ICAOCODE            ! Airport ICAO code
      CHARACTER(5)    FIPSID              ! Airport FIPS code
      CHARACTER(5)    POLNAM              ! tmp pollutant name
      CHARACTER(3)    TZONE               ! time zone

      CHARACTER(256)  LINE                ! input line buffer
      CHARACTER(256)  MESG                ! message buffer
      
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

      MESG = 'Enter logical name for a main EDMS input file'
      EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MAIN_EDMS', PROGNAME )

      MESG = 'Enter logical name for hourly EDMS inputs list'
      HDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FILELIST', PROGNAME )
      
      MESG = 'Enter logical name for ICAO airport name, country, ' //
     &       'state, county, latitude, and longitude'
      SDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'ICAO_FIPS', PROGNAME )

      MESG = 'Enter logical name for a file of EDMS source and SCC'
      KDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'EDMS_SCC', PROGNAME )

      MESG = 'Enter logical name for country, state, and county ' //
     &       'file'
      CDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'COSTCY', PROGNAME )
      
      MESG = 'Enter logical name for output annual inventory'
      ADEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTINV', PROGNAME )
      
      MESG = 'Enter logical name for output daily inventory'
      DDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTHOUR', PROGNAME )


C.......  initialize array for storing pollutant names
      ALLOCATE( OUTVAR( 50 ), STAT=IOS )
      CALL CHECKMEM( IOS, 'OUTVAR', PROGNAME )

      OUTVAR = ' '
      FIPSID = ' '

C.......  Read the country, state, and county file
      CALL RDSTCY( CDEV, 1, 0 )

C.......  Read main EDMS input file
      DO I = 1, 2

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
          IF( SIZES ) THEN
      
              ALLOCATE( TSCC( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'TSCC', PROGNAME )
              ALLOCATE( APRTID( MXSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'APRTID', PROGNAME )
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
      
              TSCC   = ''   ! arrays
              APRTID = ''   ! arrays
              LOCID  = 0
              STATE  = 0
              COUNTY = 0
              LAT    = 0.
              LON    = 0.
              HEIGHT = 0

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
                      MESG = 'ERROR: Could not find the CIAO airport '//
     &                       'code (AIRPORT) from MAIN_EMDS file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  IF( .NOT. OFLAG ) THEN
                      MESG = 'ERROR: Could not find the ORIGIN ' //
     &                       'airport location from MAIN_EMDS file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  APRT = ICAOCODE // '-' // ADJUSTL( SEGMENT( 2 ) )
                  ZLOC = STR2REAL( SEGMENT( 6 ) )

C...................  Define the location of airport id in main EDMS source list
                  K = INDEX1( APRT, NAPT, APTLOC )
                  IF( K < 1 ) THEN
                      MESG = 'ERROR: Could not find the matched EDMS '//
     &                   'source ID (' // TRIM( SEGMENT( 2 ) ) // 
     &                   ') from EDMS_SCC file'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  SCC = APTSCC( K )

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
                          HEIGHT( MXLOC ) = ZLOC * 3.28084  ! conver meter to ft
                      END IF
                  END IF

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

C.......  Allocate arrays to store output values
      ALLOCATE( ALLVAL( MXSRC, 24 ), STAT=IOS )
      CALL CHECKMEM( IOS, 'ALLVAL', PROGNAME )

      ALLVAL = 0.  ! array

C.......  Process files in input list
      MXFILES = GETFLINE( HDEV, 'FILELIST input file' )

      DO I = 1, MXFILES

          READ( HDEV, 93000, IOSTAT=IOS ) LINE

C...........  Skip blank and comment lines
          IF( BLKORCMT( LINE ) ) CYCLE

C...........  Check for i/o errors
          IF( IOS /= 0 ) THEN
              WRITE( MESG,94010 ) 'I/O error', IOS,
     &            'reading input file at line', I
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

C...........  Open input file      
          OPEN( FDEV, FILE=LINE, STATUS='OLD', IOSTAT=IOS )

          IF( IOS /= 0 ) THEN
              MESG = 'Could not open file:' // 
     &               CRLF() // BLANK5 // TRIM( LINE )
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

          MESG = 'Successfully opened EDMS hourly input file:' //
     &           CRLF() // BLANK5 // TRIM( LINE )
          CALL M3MSG2( MESG )

C...........  Initialize values before reading file
          IREC   = 0
          NPFLAG = .TRUE.

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

C...............  Read name for current pollutant
              L = INDEX( LINE, 'POLLUTANT' )
              IF ( L > 0 ) THEN
                  POLNAM = ADJUSTL( LINE( L+10:L+50 ) )
                  IF( POLNAM == 'PM25' ) POLNAM = 'PM2_5'
              END IF

C...............  Skip blank/comment/non-HOUREMIS lines
              IF( INDEX( LINE, 'HOUREMIS' ) < 1 ) CYCLE

C...............  Skip when it is out of range of episode dates
              IYR   = STR2INT( LINE( 13:14 ) )    ! integer current hour
              IMON  = STR2INT( LINE( 16:17 ) )    ! integer current hour
              IDAY  = STR2INT( LINE( 19:20 ) )    ! integer current hour
              IDATE = JULIAN( IYR, IMON, IDAY )
              IF( IDATE < SDATE .OR. IDATE > EDATE ) CYCLE

C...............  Read data from LINE
              CALL PARSLINE( LINE, MXSEG, SEGMENT )

C...............  Ensure that it read pollutant name first
              IF( POLNAM == ' ' ) THEN
                  MESG = 'ERROR: Pollutant name (POLLUTANT) was NOT ' //
     &                   'defined in EDMS hourly files.'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

C...............  Store date info  'yy mm dd'
              CHOUR = LINE( 13:23 )

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
              CDATE = CHOUR( 4:5 ) // '/' // CHOUR( 7:8 ) //
     &                    '/' // CHOUR( 1:2 )

              IF( FIRSTIME ) THEN
                  YEAR = '20' // CHOUR( 1:2 )
                  WRITE( DDEV,93000 ) '#YEAR    ' // YEAR
                  FIRSTIME = .FALSE.
              END IF

C...............  Reinitialize date for new hourly pollutant file
              IF( NPFLAG ) THEN
                  PDATE = CDATE
                  NPFLAG = .FALSE.
              END IF

C...............  Write out all stored output values to PTHOUR hourly inv file
              IF( PDATE /= CDATE ) THEN

                  DO J = 1, MXSRC

C.......................  Writing current processing hour message
                      MESG = 'Processing hourly ' // TRIM( POLNAM ) // 
     &                       ' emission on date : '// PDATE
                      IF( J == 1 ) CALL M3MSG2( MESG )
                      DAYTOT = 0.0

C.......................  compute daily total by summing 24hrs hourly emis
                      DO K = 1, 24
                          DAYTOT = DAYTOT + ALLVAL( J,K )
                      END DO
                      
C.......................  Skip any zero daily total
                      IF( DAYTOT == 0.0 )CYCLE

                      WRITE( DDEV,93020 ) STATE( J ), COUNTY( J ), 
     &                   APRTID( J ),LOCID( J ), HEIGHT( J ), POLNAM,
     &                   PDATE, TZONE, 
     &                  ( ALLVAL( J, K ), K = 1,24 ), DAYTOT, TSCC( J )

                  END DO

                  PDATE = CDATE
                  ALLVAL = 0.0

              END IF

C..............  Store output values : convert metric g/sec to short tons/hr
              ALLVAL( K, IHOUR ) =  STR2REAL( SEGMENT( 8 ) ) * 0.003968254

          END DO

C...........  Write out last date stored output values to PTHOUR hourly inv file\
C             before processing next pollutant file.
          DO J = 1, MXSRC

C...............  Writing current processing hour message
              MESG = 'Processing hourly ' // TRIM( POLNAM ) // 
     &               ' emission on date : '// CDATE
              IF( J == 1 ) CALL M3MSG2( MESG )

              DAYTOT = 0.0
C...............  compute daily total by summing 24hrs hourly emis
              DO K = 1, 24
                  DAYTOT = DAYTOT + ALLVAL( J,K )
              END DO

C...............  Skip any zero daily total
              IF( DAYTOT == 0.0 ) CYCLE

              WRITE( DDEV,93020 ) STATE( J ), COUNTY( J ), 
     &            APRTID( J ),LOCID( J ),HEIGHT( J ), POLNAM, CDATE,
     &            TZONE, ( ALLVAL( J, K ), K = 1,24 ), DAYTOT, TSCC( J )

          END DO

          CLOSE( FDEV )

C...........  Adding new pollutants
          ADD = .TRUE.

          DO J = 1, NOUTVAR
              IF( POLNAM == OUTVAR( J ) ) THEN
                  ADD = .FALSE.
                  CYCLE
              END IF
          END DO

          IF( ADD ) THEN
              NOUTVAR = NOUTVAR + 1
              OUTVAR( NOUTVAR ) = POLNAM
          END IF

      END DO

C.......  Write annual inventory header
      WRITE( ADEV,93000 ) '#IDA'
      WRITE( ADEV,93000 ) '#TYPE    Airport EDMS Inventory'
      WRITE( ADEV,93000 ) '#COUNTRY US'
      WRITE( ADEV,93000 ) '#YEAR    ' // YEAR
      
      MESG = '#DATA'
      DO I = 1, NOUTVAR
          MESG = TRIM( MESG ) // ' ' // TRIM( OUTVAR( I ) )
      END DO
      WRITE( ADEV,93000 ) TRIM( MESG )
      
C.......  Write annual inventory values
      DO I = 1, MXSRC
          WRITE( ADEV,93010 ) STATE( I ), COUNTY( I ), APRTID( I ), 
     &        LOCID( I ), HEIGHT( I ), TSCC( I ), HEIGHT( I ),
     &        LAT( I ), LON( I )
      END DO

C.......  End program successfully
      CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C.......  Formatted file I/O formats...... 93xxx
93000 FORMAT( A )
93010 FORMAT( I2.2, I3.3, A15, I15, I12, 54X, A10, 8x, I4, 107X, 2F9.4 )
93020 FORMAT( I2.2, I3.3, A15, 2I12,12X, A5, A8, A3, 24E7.1,E8.2,1X,A10)

C.......  Internal buffering formats...... 94xxx
94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )
94011 FORMAT( 10 ( A, :, F10.3, :, 2X  ) )

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
                LATICAO = STR2REAL( SEGMENT( 3) )
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
