      PROGRAM BLUESKY2INV

C***********************************************************************
C  program body starts at line  
C
C  DESCRIPTION:
C       This program converts the output from the BlueSky CONSUME model
C       into daily and annual inventories that can be read by SMOKE.
C
C  PRECONDITIONS REQUIRED:
C       CONSUME model has been run
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 1/2005 by C. Seppanen
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
      INTEGER       GETFLINE
      INTEGER       INDEX1
      INTEGER       PROMPTFFILE
      INTEGER       STR2INT
      REAL          STR2REAL
      LOGICAL       STRLIST
      
      EXTERNAL      CRLF, ENVINT, FIND1, INDEX1, GETFLINE, PROMPTFFILE, 
     &              STR2INT, STR2REAL, STRLIST

C.......  LOCAL PARAMETERS
      CHARACTER(50), PARAMETER :: 
     & CVSW = '$Name SMOKEv4.9_Jun2022$'  ! CVS release tag
      INTEGER, PARAMETER :: MXSEG = 15    ! number of segments in line
      INTEGER, PARAMETER :: MXOUT = 11    ! number of output values
      REAL,    PARAMETER :: MTON2STON = 1.10231131  ! conversion factor for metric to short tons

C.......  LOCAL VARIABLES

C.......  Static arrays
      CHARACTER(80)   SEGMENT( MXSEG )    ! parsed input line
      CHARACTER(5) :: ALLNAM ( MXOUT ) =
     &    ( / 'AREA ', 'HFLUX', 'PM2_5', 'PM10 ', 'PM   ', 'PMC  ',
     &        'CO   ', 'CO2  ', 'CH4  ', 'NMHC ', 'TOG  ' / )
      CHARACTER(5)    VARLIST( MXOUT )    ! output variable list
      REAL            ALLVAL ( MXOUT )    ! output variable values

C.......  Allocatable arrays
      CHARACTER(15), ALLOCATABLE :: FIREID( : )
      INTEGER,       ALLOCATABLE :: LOCID ( : )
      INTEGER,       ALLOCATABLE :: STATE ( : )
      INTEGER,       ALLOCATABLE :: COUNTY( : )
      REAL,          ALLOCATABLE :: LAT   ( : )
      REAL,          ALLOCATABLE :: LON   ( : )
      INTEGER,       ALLOCATABLE :: OUTVAR( : )

C.......  File units and logical names
      INTEGER         LDEV                ! unit number for log file
      INTEGER         IDEV                ! list of input files
      INTEGER         FDEV                ! input file
      INTEGER         CDEV                ! country, state, and county file
      INTEGER         ADEV                ! output annual inventory
      INTEGER         DDEV                ! output daily inventory

C.......  Other local variables
      INTEGER         I, J, K             ! counters and indices
      INTEGER         IOS                 ! i/o status
      INTEGER         IREC                ! line counter
      INTEGER         MXFIRES             ! max. number of fires
      INTEGER         STID                ! state code
      INTEGER         CYID                ! county code
      INTEGER         LOCVAL              ! location identifier
      INTEGER         STORED              ! number of stored fires in master list
      INTEGER         FIP                 ! state and county code
      INTEGER         IDATE               ! start date of data
      INTEGER         IMONTH              ! start month of data
      INTEGER         IDAY                ! start day of data
      INTEGER         NOUTVAR             ! number of output variables

      REAL            LATVAL              ! fire latitude
      REAL            LONVAL              ! fire longitude

      LOGICAL         ADD                 ! true: add current fire to master list

      CHARACTER(15)   FIRE                ! fire identifier
      CHARACTER(10)   SCC                 ! SCC
      CHARACTER(3)    TZONE               ! time zone
      CHARACTER(8)    DATE                ! fire date
      CHARACTER(7)    CDATE               ! data start date
      CHARACTER(4)    YEAR                ! inventory year
      CHARACTER(2)    CMONTH              ! data start month
      CHARACTER(2)    CDAY                ! data start day
      CHARACTER(8)    DEFDATE             ! default date for fake daily entries
      
      CHARACTER(256)  LINE                ! input line buffer
      CHARACTER(256)  MESG                ! message buffer
      
      CHARACTER(16) :: PROGNAME = 'BLUESKY2INV' ! program name

C***********************************************************************
C   begin body of program BLUESKY2INV

      LDEV = INIT3()

C.......  Write copyright, version, web address, header info, and prompt
C         to continue running the program.
      CALL INITEM( LDEV, CVSW, PROGNAME )

C.......  Get environment variable settings
      MESG = 'Enter logical name for inputs list'
      IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FILELIST', PROGNAME )
      
      MESG = 'Enter logical name for country, state, and county ' //
     &       'file'
      CDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'COSTCY', PROGNAME )
      
      MESG = 'Enter logical name for output annual inventory'
      ADEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTINV', PROGNAME )
      
      MESG = 'Enter logical name for output daily inventory'
      DDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTDAY', PROGNAME )

C.......  Get list of output variables
      MESG = 'Output variable names'
      IF( .NOT. STRLIST( 'VARLIST', MESG, MXOUT, 
     &                   NOUTVAR, VARLIST ) ) THEN
          MESG = 'Could not read list of output variable names'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END IF

      ALLOCATE( OUTVAR( NOUTVAR ), STAT=IOS )
      CALL CHECKMEM( IOS, 'OUTVAR', PROGNAME )

C.......  Check that requested variables are valid
      J = 0
      DO I = 1, NOUTVAR
          K = INDEX1( VARLIST( I ), MXOUT, ALLNAM )
          IF( K > 0 ) THEN
              J = J + 1
              OUTVAR( J ) = K
          ELSE
              MESG = 'WARNING: Requested variable "' // 
     &               TRIM( VARLIST( I ) ) //'" is not valid'
              CALL M3MSG2( MESG )
          END IF
      END DO
      NOUTVAR = J

      IF( NOUTVAR == 0 ) THEN
          MESG = 'No valid output variables'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END IF

C.......  Get start date of data      
      IDATE = ENVINT( 'G_STDATE', 'Start date (YYYYDDD)', 0, IOS )

C.......  Set inventory year based on start date      
      WRITE( CDATE, '(I7)' ) IDATE
      YEAR = CDATE( 1:4 )

C.......  Default date will be day before start date; convert date to string
      IDATE = IDATE - 1
      CALL DAYMON( IDATE, IMONTH, IDAY )
      
      WRITE( CDATE, '(I7)' ) IDATE
      WRITE( CMONTH,'(I2)' ) IMONTH
      IF( CMONTH( 1:1 ) == ' ' ) CMONTH( 1:1 ) = '0'
      WRITE( CDAY,  '(I2)' ) IDAY
      IF( CDAY( 1:1 ) == ' ' ) CDAY( 1:1 ) = '0'
      
      DEFDATE = CMONTH // '/' // CDAY // '/' // CDATE( 3:4 )

C.......  Read the country, state, and county file
      CALL RDSTCY( CDEV, 1, 0 )

C.......  Write daily inventory header
      WRITE( DDEV,93000 ) '#COUNTRY US'

C.......  Determine number of lines in filelist; this will be the maximum
C         number of fires
      MXFIRES = GETFLINE( IDEV, 'List of input files' )
      
      ALLOCATE( FIREID( MXFIRES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'FIREID', PROGNAME )
      ALLOCATE( LOCID( MXFIRES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'LOCID', PROGNAME )
      ALLOCATE( STATE( MXFIRES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'STATE', PROGNAME )
      ALLOCATE( COUNTY( MXFIRES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'COUNTY', PROGNAME )
      ALLOCATE( LAT( MXFIRES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'LAT', PROGNAME )
      ALLOCATE( LON( MXFIRES ), STAT=IOS )
      CALL CHECKMEM( IOS, 'LON', PROGNAME )
      
      FIREID = ''   ! arrays
      LOCID  = 0
      STATE  = 0
      COUNTY = 0
      LAT    = 0.
      LON    = 0.
      
      STORED = 0

C.......  Process files in input list
      DO I = 1, MXFIRES

          READ( IDEV, 93000, IOSTAT=IOS ) LINE

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

          MESG = 'Successfully opened input file:' //
     &           CRLF() // BLANK5 // TRIM( LINE )
          CALL M3MESG( MESG )

C...........  Initialize values before reading file
          IREC   = 0
          LOCVAL = 1
          ALLVAL = 0.  ! array

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
              
              CALL PARSLINE( LINE, MXSEG, SEGMENT )
              
C...............  Get state and county code from fire ID
              FIRE = SEGMENT( 2 )(  1:12 )
              STID = STR2INT( SEGMENT( 2 )( 14:15 ) )
              CYID = STR2INT( SEGMENT( 2 )( 16:18 ) )
              
              LATVAL = STR2REAL( SEGMENT( 3 ) )
              LONVAL = STR2REAL( SEGMENT( 4 ) )

C...............  Check if we already have this fire in master list
              IF( IREC == 1 ) THEN
                  ADD = .TRUE.
              
                  DO J = 1, STORED
                      IF( FIRE == FIREID( J ) .AND.
     &                    STID == STATE ( J ) .AND.
     &                    CYID == COUNTY( J ) ) THEN
                          
                          IF( LATVAL == LAT( J ) .AND.
     &                        LONVAL == LON( J ) ) THEN
                              ADD = .FALSE.
                              EXIT
                          ELSE
                              LOCVAL = LOCID( J ) + 1
                          END IF
                      END IF
                  END DO
              
                  IF( ADD ) THEN
                      STORED = STORED + 1
                      FIREID( STORED ) = FIRE
                      LOCID ( STORED ) = LOCVAL
                      STATE ( STORED ) = STID
                      COUNTY( STORED ) = CYID
                      LAT   ( STORED ) = LATVAL
                      LON   ( STORED ) = LONVAL
                  END IF
              END IF
              
              SCC = '2810001000'
              
              DATE = SEGMENT( 1 )( 5:6 ) // '/' //
     &               SEGMENT( 1 )( 7:8 ) // '/' //
     &               SEGMENT( 1 )( 3:4 )

C...............  Look up time zone based on state and county
              FIP = STID * 1000 + CYID
              K = FIND1( FIP, NCOUNTY, CNTYCOD )
              
              IF( K <= 0 ) THEN
                  WRITE( MESG,94010 ) 'Could not find FIPS code',
     &                FIP, 'in COSTCY file'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF
              
              TZONE = CNTYTZNM( K )

C...............  Store output values; only store area and heat flux from
C                 first record
              IF( IREC == 1 ) THEN
                  ALLVAL( 1 )  = STR2REAL( SEGMENT( 6 ) )              ! area
                  ALLVAL( 2 )  = STR2REAL( SEGMENT( 8 ) )              ! hflux
              END IF
              ALLVAL( 3 )  = ALLVAL( 3 )  + STR2REAL( SEGMENT( 9 ) )   ! pm2.5
              ALLVAL( 4 )  = ALLVAL( 4 )  + STR2REAL( SEGMENT( 10 ) )  ! pm10
              ALLVAL( 5 )  = ALLVAL( 5 )  + STR2REAL( SEGMENT( 11 ) )  ! pm
             !ALLVAL( 6 )                                              ! pmc
              ALLVAL( 7 )  = ALLVAL( 7 )  + STR2REAL( SEGMENT( 12 ) )  ! co
              ALLVAL( 8 )  = ALLVAL( 8 )  + STR2REAL( SEGMENT( 13 ) )  ! co2
              ALLVAL( 9 )  = ALLVAL( 9 )  + STR2REAL( SEGMENT( 14 ) )  ! ch4
              ALLVAL( 10 ) = ALLVAL( 10 ) + STR2REAL( SEGMENT( 15 ) )  ! nmhc
             !ALLVAL( 11 )                                             ! tog
              
          END DO

C...........  Don't write output if no lines in file
          IF( IREC == 0 ) CYCLE

C...........  Calculate PMC and TOG
          ALLVAL( 6 )  = ALLVAL( 4 ) - ALLVAL( 3 )
          ALLVAL( 11 ) = ALLVAL( 9 ) + ALLVAL( 10 )

C...........  Convert values from metric tons to short tons
          DO J = 3, MXOUT 
              ALLVAL( J ) = ALLVAL( J ) * MTON2STON
          END DO

C...........  Write default daily entry; all sources in the daily
C             inventory must have an entry for the same day - use
C             the day before the start date
          IF( ADD ) THEN
              WRITE( DDEV,93020 ) STATE( STORED ), COUNTY( STORED ),
     &            FIREID( STORED ), LOCID( STORED ), 
     &            ALLNAM( OUTVAR( 1 ) ), DEFDATE, TZONE, 0.,
     &            '2810001000'
          END IF

C...........  Write daily inventory values; use a different format for
C             heat flux since the values are so large
          DO J = 1, NOUTVAR
              IF( ALLNAM( OUTVAR( J ) ) == 'HFLUX' ) THEN
                  WRITE( DDEV,93030 ) STATE( STORED ), COUNTY( STORED ),
     &                FIREID( STORED ), LOCID( STORED ), 
     &                ALLNAM( OUTVAR( J ) ), DATE, TZONE,
     &                ALLVAL( OUTVAR( J ) ), '2810001000'
              ELSE
                  WRITE( DDEV,93020 ) STATE( STORED ), COUNTY( STORED ),
     &                FIREID( STORED ), LOCID( STORED ), 
     &                ALLNAM( OUTVAR( J ) ), DATE, TZONE,
     &                ALLVAL( OUTVAR( J ) ), '2810001000'
              END IF
          END DO
      
      END DO

C.......  Write annual inventory header
      WRITE( ADEV,93000 ) '#IDA'
      WRITE( ADEV,93000 ) '#TYPE    Fire Event Inventory'
      WRITE( ADEV,93000 ) '#COUNTRY US'
      WRITE( ADEV,93000 ) '#YEAR    ' // YEAR
      
      MESG = '#DATA'
      DO I = 1, NOUTVAR
          MESG = TRIM( MESG ) // ' ' // TRIM( ALLNAM( OUTVAR( I ) ) )
      END DO
      WRITE( ADEV,93000 ) TRIM( MESG )
      
C.......  Write annual inventory values
      DO I = 1, STORED
          WRITE( ADEV,93010 ) STATE( I ), COUNTY( I ), FIREID( I ), 
     &        LOCID( I ), '2810001000', LAT( I ), LON( I )
      END DO

C.......  End program successfully
      CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C.......  Formatted file I/O formats...... 93xxx
93000 FORMAT( A )
93010 FORMAT( I2, I3, A15, I15, 66X, A10, 119X, F9.4, F9.4 )
93020 FORMAT( I2, I3, A15, I12, 24X, A5, A8, A3, F18.5, 1X, A10 )
93030 FORMAT( I2, I3, A15, I12, 24X, A5, A8, A3, F18.1, 1X, A10 )

C.......  Internal buffering formats...... 94xxx
94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )

      END PROGRAM BLUESKY2INV
