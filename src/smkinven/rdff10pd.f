
        SUBROUTINE RDFF10PD( FDEV, TZONE, TSTEP, MXPDSRC, GETSIZES, 
     &                      GETCOUNT, FIRSTCALL, DAYFLAG, SDATE, STIME, 
     &                      EDATE, ETIME, EASTAT, SPSTAT )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads the day- or hour-specific emissions in
C      FF10_HOURLY and FF10_DAILY format. It appends the records to 
C      the global storage from the MODDAYHR
C
C  PRECONDITIONS REQUIRED:
C      Must complete processing annual/avg inventory (SMK_AVEINV_YN=Y)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created by B.H. Baek on 8/2011
C
C***************************************************************************
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is the inventory arrays
        USE MODSOURC, ONLY: CIFIP, CSOURC, INTGRFLAG, CINTGR

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVTBL, ITFACA, ITNAMA,
     &                      ITKEEPA, SORTCAS, SCASIDX, NUNIQCAS, INVDVTS,
     &                      UCASNPOL, UNIQCAS, UCASIDX, UCASNKEP, MXIDAT,
     &                      INVDNAM

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NIPPA, NSRC, EANAM, NCHARS, INV_MON

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, IDXSRC, SPDIDA,
     &                      CODEA, EMISVA, DYTOTA, CIDXA

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, USEDAYLT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        INTEGER      ENVINT
        LOGICAL      ENVYN, CHKINT
        INTEGER      FIND1
        INTEGER      FINDC
        INTEGER      INDEX1
        INTEGER      JULIAN
        INTEGER      SECSDIFF
        INTEGER      STR2INT
        REAL         STR2REAL
        REAL         YR2DAY
        INTEGER      YEAR4
        INTEGER      GETTZONE
        LOGICAL      ISDSTIME
        LOGICAL      USEEXPGEO

        EXTERNAL     CRLF, ENVINT, ENVYN, FIND1, FINDC, INDEX1, JULIAN, 
     &               SECSDIFF, STR2INT, STR2REAL, YEAR4, YR2DAY, CHKINT,
     &               GETTZONE, ISDSTIME, USEEXPGEO

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: FDEV           ! file unit no.
        INTEGER, INTENT (IN)  :: TZONE          ! output time zone
        INTEGER, INTENT (IN)  :: TSTEP          ! time step HHMMSS
        INTEGER, INTENT (IN)  :: MXPDSRC        ! max. day- or hr-specific source
        LOGICAL, INTENT (IN)  :: GETSIZES       ! true: get no. time steps & pols
        LOGICAL, INTENT (IN)  :: GETCOUNT       ! true: get max no. srcs per time
        LOGICAL, INTENT (IN)  :: FIRSTCALL      ! true: first call of a loop
        LOGICAL, INTENT (IN)  :: DAYFLAG        ! true: day-, false: hour-spec
        INTEGER, INTENT(INOUT):: SDATE          ! Julian starting date in TZONE
        INTEGER, INTENT(INOUT):: STIME          ! start time of data in TZONE
        INTEGER, INTENT(OUT)  :: EDATE          ! Julian ending date in TZONE
        INTEGER, INTENT(OUT)  :: ETIME          ! ending time of data in TZONE
        INTEGER, INTENT(OUT)  :: EASTAT( NIPPA ) ! true: pol/act appears in data
        INTEGER, INTENT(OUT)  :: SPSTAT( MXSPDAT ) ! true: special in data

C...........   Local list of bad sources to prevent duplicate writing of error
C              messages
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )

C...........   Local parameters
        INTEGER, PARAMETER :: MXSEG = 60     ! max no of segments

C...........   Local segment arrays
        CHARACTER( 30 )    :: SEGMENT( MXSEG ) = ' '  ! temporary line segments

C...........   Local list of FIPS start/end positions to facilitate
C              faster lookups
        INTEGER, ALLOCATABLE, SAVE :: STARTSRC( : )
        INTEGER, ALLOCATABLE, SAVE :: ENDSRC( : )

C...........   Local list of arrays for warning handling
        LOGICAL, ALLOCATABLE, SAVE :: WARNKEEP( : ) ! true: write warning for Keep = N
        LOGICAL, ALLOCATABLE, SAVE :: WARNMULT( : ) ! true: write warning for Multiple pollutants from a single pollutant in Inventory Table

C...........   Temporary read arrays
        REAL            TDAT( 31,24 )       ! temporary data values

C...........   Other local variables
        INTEGER          D, H, HS, I, J, N, NV, L, LL, L1, L2, S, T    ! counters and indices
        INTEGER          ES, NS, SS    ! end src, tmp no. src, start sourc

        INTEGER          CIDX             ! tmp data index
        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER, SAVE :: ICC = 0          ! tmp country code from header
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          JTIME            ! tmp HHMMSS time
        INTEGER          LYEAR            ! leap year
        INTEGER, SAVE :: LOOPNO = 0       ! no. of loops
        INTEGER, SAVE :: MAXPTR           ! maximum time step reference pointer
        INTEGER, SAVE :: MINPTR           ! minimum time step reference pointer
        INTEGER          MONTH            ! tmp month number
        INTEGER, SAVE :: MXWARN       	  ! max no. warnings
        INTEGER, SAVE :: NBADSRC = 0      ! no. bad sources
        INTEGER, SAVE :: NFIELD = 1       ! number of data fields
        INTEGER       :: NPOA   = 0       ! unused header number of pol/act
        INTEGER, SAVE :: NSTEPS = 0       ! number of time steps
        INTEGER, SAVE :: NWARN( 5 )       ! warnings counter
        INTEGER          PTR              ! tmp time step pointer
        INTEGER       :: RDATE = 1980001  ! reference date: Jan 1, 1980
        INTEGER       :: RTIME = 0        ! reference time
        INTEGER, SAVE :: S1 = 0           ! saved 1st position of extra field
        INTEGER, SAVE :: S2 = 0           ! saved 2nd position of extra field
        INTEGER, SAVE :: NDAYS    = 0     ! saved processing days
        INTEGER, SAVE :: B_YEAR   = 0     ! saved processing base year
        INTEGER, SAVE :: SDATESAV = 0     ! saved start date
        INTEGER, SAVE :: STIMESAV = 0     ! saved start time
        INTEGER, SAVE :: TDIVIDE  = 1     ! time step divisor
        INTEGER          WD               ! tmp field width
        INTEGER          YEAR             ! 4-digit year
        INTEGER       :: YR4 = 0          ! unused header year
        INTEGER          ZONE             ! source time zones
        INTEGER          DZONE            ! time shift (ZONE-TZONE)
        INTEGER          FSTPRV, FSTDATE, LSTDATE, FSTLOC, LSTLOC   ! processing start/end dates

        REAL             CONVFAC          ! tmp conversion factor from Inventory Table
        REAL             TOTAL            ! tmp daily total of hourly file

        LOGICAL, SAVE :: DFLAG  = .FALSE. ! true: dates set by data
        LOGICAL       :: EFLAG  = .FALSE. ! TRUE iff ERROR
        LOGICAL       :: WARNOUT = .FALSE.! true: then output warnings
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called
        LOGICAL, SAVE :: LFLAG = .FALSE.  ! true: output daily/hourly inv in local time
        LOGICAL, SAVE :: SFLAG = .FALSE.  ! true: use daily total from hourly
        LOGICAL, SAVE :: TFLAG = .FALSE.  ! true: use SCCs for matching with inv

        CHARACTER(256) :: BUFFER = ' '    ! src description buffer 
        CHARACTER(1920):: LINE   = ' '    ! line buffer 
        CHARACTER(512) :: MESG   = ' '    ! message buffer

        CHARACTER(FIPLEN3) CFIP      ! tmp co/st/cy code
        CHARACTER(FIPLEN3) LFIP      ! previous st/co FIPS code
        CHARACTER(CASLEN3) CDAT      ! tmp Inventory data (input) name
        CHARACTER(IOVLEN3) CNAM,PNAM ! tmp SMOKE name
        CHARACTER(PLTLEN3) FCID      ! tmp facility ID
        CHARACTER(CHRLEN3) SKID      ! tmp stack ID
        CHARACTER(CHRLEN3) DVID      ! tmp device ID
        CHARACTER(CHRLEN3) PRID      ! tmp process ID
        CHARACTER(SCCLEN3) TSCC      ! tmp source category code
        CHARACTER(ALLLEN3) CSRC      ! tmp source string

        CHARACTER(16) :: PROGNAME = 'RDFF10PD' !  program name

C***********************************************************************
C   begin body of program RDFF10PD
C.........  First time routine called
        IF( FIRSTIME ) THEN

C.............  Get environment variable using an hourly file as a daily file
C.............  NOTE - the hourly file will have been assigned as a daily
C               file when it was opened.
            MESG = 'Use daily totals only from hourly data file'
            SFLAG = ENVYN( 'HOURLY_TO_DAILY', MESG, .FALSE., IOS )

C.............  No time zone shift for AERMOD support
            MESG = 'Outputs local time daily and/or hourly inventories (No time shift)'
            LFLAG = ENVYN( 'OUTPUT_LOCAL_TIME', MESG, .FALSE., IOS )

C.............  Get processing base year info
            MESG = 'Define Processing Base Year for daily/hourly-specific inventory'
            B_YEAR = ENVINT( 'BASE_YEAR', MESG, 0, IOS )
            IF( B_YEAR == 0 ) THEN
                MESG = 'ERROR: MUST define the processing base year for '//
     &                 'daily/hourly-specific inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET , ' ', 100, I )

C.............  Give note if file is being read as a daily file
            IF( DAYFLAG .AND. SFLAG ) THEN
                MESG = 'NOTE: Daily data only being used from an ' //
     &                 'hourly emissions file'
                CALL M3MSG2( MESG )

C.............  Otherwise, ignore setting because it is an hourly file
            ELSE IF( SFLAG ) THEN
                SFLAG = .FALSE.
                MESG = 'NOTE: Ignoring HOURLY_TO_DAILY setting for ' //
     &                 'reading hourly emissions data'
                CALL M3MSG2( MESG )
            END IF

C.............  Allocate memory for bad source storage
            ALLOCATE( BADSRC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

C.............  Create unique list of FIPS codes and other things
            CALL GENUSLST

C.............  Build helper arrays for making searching faster
            ALLOCATE( STARTSRC( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STARTSRC', PROGNAME )
            ALLOCATE( ENDSRC( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ENDSRC', PROGNAME )
            STARTSRC = 0
            ENDSRC = 0
            S = 0
            DO I = 1, NINVIFIP
                DO
                    S = S + 1
                    IF ( S .GT. NSRC ) EXIT
                    IF( CIFIP( S ) .EQ. INVCFIP( I ) ) THEN
                        IF( STARTSRC( I ) .EQ. 0 ) STARTSRC( I ) = S
                        ENDSRC( I ) = S
                    ELSE
                        S = S - 1
                        EXIT   
                    END IF
                END DO
            END DO

C.............  Initialize warnings counter
            NWARN = 0  ! array

            FIRSTIME = .FALSE.

        END IF

C.........  For the first call in a loop of files, initialize variables
        IF( FIRSTCALL ) THEN
            MINPTR  = 99999999
            MAXPTR  = 0

C.............  Set time step divisor
            TDIVIDE = 3600 * TSTEP / 10000

C.............  If dates have been set by the data, set the number of steps
C               steps
            IF( DFLAG ) THEN
                NSTEPS = 1+ SECSDIFF( SDATE,STIME,EDATE,ETIME ) / TDIVIDE
                SDATESAV = SDATE
                STIMESAV = STIME
            END IF

C.............  Set switch for printing errors only the first loop through all
C               of the input files.  The second time through is indicated
C               for the second time that FIRSTCALL is true.  
C.............  Reset loop counter if call is to get dimensions only (because
C               this means it is the first call or daily or hourly)
            IF( GETSIZES ) LOOPNO = 0
            LOOPNO = LOOPNO + 1
            WARNOUT = ( LOOPNO .EQ. 1 )

C.............  Deallocate warning arrays
            IF( ALLOCATED( WARNKEEP ) ) DEALLOCATE( WARNKEEP, WARNMULT )
            ALLOCATE( WARNKEEP( NUNIQCAS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WARNKEEP', PROGNAME )
            ALLOCATE( WARNMULT( NUNIQCAS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WARNMULT', PROGNAME )
            WARNKEEP = .TRUE.
            WARNMULT = .TRUE.

        END IF

C.........  Loop through file and read it. In the first section, determine
C           the minimum and maximum date. Use a reference date to do this. In
C           the second section, determine the number of records per time 
C           step. In the third section, read and store the data.  When storing
C           data, time step index is computed from the start date/time instead
C           of the reference date/time so that the indexing will work properly.
        IREC = 0
        TDAT = 0.0   !  array
        DO         !  Head of period-specific file read loop

C.............  Read first line of file
            READ( FDEV, 93000, END=299 ) LINE
            IREC = IREC + 1

            L = LEN_TRIM( LINE )

C.............  Skip blank lines 
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set
C               properly
            CALL GETHDR( 1, .FALSE., .FALSE., .FALSE.,
     &                   LINE, ICC, YR4, NPOA, IOS )

C.............  Interpret error status
            IF( IOS .EQ. 4 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: DATA header entry should not be used'//
     &                 'for EMS-95 day- or hour-specific files.'
                CALL M3MSG2( MESG )

            ELSE IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Parse line into segments
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Determine if file is day- or hour-specific by the length of the
C               lines. Make sure day- and hour-specific data are not in the
C               same file.
C.............  If the file is hourly but the only the daily is to be read, then
C               behave as if it is a daily file.

C.............  Skip column header line
            IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) CYCLE 

C.............  Set Julian day from MMDDYY8 SAS format
            IF( DAYFLAG ) THEN
                YEAR  = YR4
                MONTH = STR2INT( SEGMENT( 13 ) )
                DAY   = 1
            ELSE
                YEAR  = STR2INT( TRIM( SEGMENT(13)( 1:4 ) ) )
                MONTH = STR2INT( TRIM( SEGMENT(13)( 5:6 ) ) )
                DAY   = STR2INT( TRIM( SEGMENT(13)( 7:8 ) ) )
            END IF

            JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )
            JTIME = 0

C.............  Set the number of fields, depending on day- or hour-specific
            IF( DAYFLAG ) THEN
                NFIELD  = MON_DAYS( MONTH )
                LYEAR =  INT( 1 / YR2DAY( YEAR ) )   ! convert year to days
                IF( LYEAR > 365 .AND. MONTH == 2 ) NFIELD = 29
                FSTLOC = 1
                LSTLOC = NFIELD
            ELSE
                NFIELD = 1
                FSTLOC = 1
                LSTLOC = 24
            END IF

C.............  Skip non-processing month/day 
            IF( INV_MON > 0 ) THEN
                NDAYS  = MON_DAYS( INV_MON )
                LYEAR =  INT( 1 / YR2DAY( B_YEAR ) )   ! convert year to days
                IF( LYEAR > 365 .AND. MONTH == 2 ) NDAYS = 29
                FSTDATE = 1000 * B_YEAR + JULIAN( B_YEAR, INV_MON, 1 )
                LSTDATE = 1000 * B_YEAR + JULIAN( B_YEAR, INV_MON, NDAYS )
                CALL NEXTIME( FSTDATE, JTIME, -240000 )     ! include the last day of previous month
                CALL NEXTIME( LSTDATE, JTIME,  240000 )     ! include the first day of next month 

                IF( DAYFLAG ) THEN

                    FSTPRV = FSTDATE  ! reset the last of previous month to first day of prv month for a proper daily inv 
                    CALL DAYMON ( FSTPRV, MONTH, N )
                    N = MON_DAYS( MONTH )
                    CALL NEXTIME( FSTPRV, JTIME, -240000 * (N-1) )

                    IF( LSTDATE == JDATE ) THEN
                        NFIELD = 1
                        FSTLOC = 1
                        LSTLOC = 1
                        JDATE  = LSTDATE
                    ELSE IF( FSTPRV == JDATE ) THEN
                        NFIELD = 1
                        FSTLOC = N
                        LSTLOC = N
                        JDATE  = FSTDATE
                    ELSE
C.........................  Check start/end dates with emissions within the month
                        S1 = 15   ! pollutant field start position
                        FSTLOC = 0
                        LSTLOC = 0

                        DO J = 1, NFIELD
                            IF( STR2REAL( SEGMENT( S1-1+J ) ) > 0.0  ) THEN
                                FSTLOC = J
                                EXIT
                            END IF
                        END DO

                        DO J = NFIELD, 1, -1
                            IF( STR2REAL( SEGMENT( S1-1+J ) ) > 0.0 ) THEN
                                LSTLOC = J
                                EXIT
                            END IF
                        END DO

                        IF( FSTLOC == 0 .OR. LSTLOC == 0 ) THEN
                            NFIELD = 0
                        ELSE
                            NFIELD = ( LSTLOC - FSTLOC ) + 1
                            JDATE = JDATE + FSTLOC - 1           ! Update actual processing date
                        END IF

                    END IF

                END IF

                IF( .NOT. ( FSTDATE <= JDATE .AND. JDATE <= LSTDATE ) ) CYCLE

            END IF

            IF( NFIELD < 1 ) CYCLE

C.............  Read FIPS code
            CFIP = REPEAT( '0', FIPLEN3 )
            IF( USEEXPGEO() ) THEN
                CFIP(  1: 3 ) = ADJUSTR( SEGMENT( 1 )( 1:3 ) )
                CFIP(  4: 9 ) = ADJUSTR( SEGMENT( 2 )( 1:6 ) )
                CFIP( 10:12 ) = ADJUSTR( SEGMENT( 3 )( 1:3 ) )
            ELSE
                WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC  ! country code of FIPS        
                CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 2 )( 1:5 ) )  ! state/county code
            END IF

C.............  Search for time zone for current county
            I = FINDC( CFIP, NCOUNTY, CNTYCOD )

C.............  If time zone name is not found, thenoutput error
            IF( I .LE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find time zone for county: '//
     &                 CFIP // ' from COSTCY file'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Local time shift flag (county-specific)
            IF( .NOT. LFLAG ) THEN

C.................  Set time zone number
                ZONE = GETTZONE( CFIP )
 
C.................  If daily emissions are not in the output time zone 
                IF( WARNOUT .AND. .NOT. DAYFLAG .AND. ZONE .NE. TZONE .AND .
     &              NWARN( 1 ) .LE. MXWARN ) THEN
                    WRITE( MESG,94010 ) 
     &                  'WARNING: Time zone ', ZONE, 'in hourly-specific ' //
     &                  'file at line', IREC, CRLF() // BLANK10 //  
     &                  'does not match output time zone', TZONE
                    CALL M3MESG( MESG )
                    NWARN( 1 ) = NWARN( 1 ) + 1
                END IF

C.................  Check if date is in daylight time, if local zone has
C                   already been converted, and if this FIPS code is
C                   exempt from daylight time or not.
                IF( ISDSTIME( JDATE ) .AND. USEDAYLT( I ) ) THEN
                    ZONE = ZONE - 1
                END IF
               
                DZONE = ZONE - TZONE

C.............  Reset time shift to 0 to correctly compute local time zone
            ELSE

                DZONE = 0
                
            END IF

C.............  Convert date and time to output time zone.
            CALL NEXTIME( JDATE, JTIME, DZONE * 10000 )

C.............  Determine time step pointer based on reference time
            PTR = SECSDIFF( RDATE, RTIME, JDATE, JTIME ) / TDIVIDE + 1
            
C.............  Store minimum time step number as compared to reference
            IF( PTR .LT. MINPTR ) MINPTR = PTR

C.............  Store maximum time step number as compared to reference
            IF( DAYFLAG ) THEN
                PTR = SECSDIFF( RDATE, RTIME, JDATE+NFIELD-1, JTIME ) / TDIVIDE + 1
            END IF

            IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

C.............  Find source ID
C.............  Set key for searching sources
            IF( CATEGORY == 'POINT' ) THEN
                FCID = ADJUSTL( SEGMENT( 4 ) )   ! EIS_FACILITY_ID in FF10&IDA (PlantID in ORL)
                SKID = ADJUSTL( SEGMENT( 5 ) )   ! EIS_UNIT_ID in FF10&IDA (PointID in ORL)
                DVID = ADJUSTL( SEGMENT( 6 ) )   ! EIS_REL_POINT_ID in FF10&IDA (StackID in ORL)
                PRID = ADJUSTL( SEGMENT( 7 ) )   ! EIS_PROCESS_ID in FF10&IDA (SegmentID in ORL)
            END IF

            TSCC = ' '

C.............  If FIPS code is not the same as last time, then
C               look it up and get indidies
            IF( CFIP .NE. LFIP ) THEN
                J = FINDC( CFIP, NINVIFIP, INVCFIP )
                IF( J .LE. 0 ) THEN
                    WRITE( MESG,93000 ) 'INTERNAL ERROR: Could not '//
     &                     'find FIPS ' // CFIP // ' in internal list.'
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                END IF

                SS = STARTSRC( J )
                ES = ENDSRC( J )
                NS = ES - SS + 1
                LFIP = CFIP

            END IF

C.............  If SCCs are needed for matching...
            IF ( TFLAG ) THEN

                TSCC = ADJUSTL( SEGMENT( 8 ) )     ! SCC from FF10 HRDAY format
                IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )

C.................  Build source characteristics field for searching inventory
                IF( CATEGORY == 'POINT' ) THEN
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                ELSE
                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CSRC )
                END IF
                
C.................  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )

C.............  If SCCs are not being used for matching (at least not yet)...
            ELSE

C.................  Build source characteristics field for searching inventory
                IF( CATEGORY == 'POINT' ) THEN
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                ELSE
                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CSRC )
                END IF

C.................  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )

C.................  If source is not found for day-specific processing, see 
C                   if reading the SCC in helps (needed for IDA format)
                IF( J .LE. 0 ) THEN

                    TSCC = ADJUSTL( SEGMENT( 8 ) )     ! SCC from FF10 HRDAY format
                    IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )

C.....................  Build source characteristics field for searching inventory
                    IF( CATEGORY == 'POINT' ) THEN
                        CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                    ELSE
                        CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, 
     &                               CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                               CHRBLNK3, CSRC )
                    END IF

C.....................  Search for this record in sources
                    J = FINDC( CSRC, NS, CSOURC( SS ) )
                    IF ( J .GT. 0 ) TFLAG = .TRUE.

                END IF

            END IF

C.............  Store source in list of bad sources
C.............  Print warning about sources not found in the inventory
            IF( J .LE. 0 ) THEN

C.................  Search for source in list of bad sources
                J = INDEX1( CSRC, NBADSRC, BADSRC )

C.................  If source is not found, give a message.  Don't need the
C                   WARNOUT controller because this section only gets
C                   invoked once.
                IF( J .LE. 0 ) THEN

                    NBADSRC = NBADSRC + 1
                    BADSRC( NBADSRC ) = CSRC

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                    IF( NWARN( 3 ) .LE. MXWARN ) THEN
                        MESG = 'WARNING: Period-specific record does '//
     &                         'not match inventory sources: '//
     &                         CRLF() // BLANK10 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN( 3 ) = NWARN( 3 ) + 1
                    END IF

                END IF

                CYCLE               !  to head of read loop

C.............  Otherwise, update master list of sources in the inventory
            ELSE
                S = SS - 1 + J         ! calculate source number

            END IF

C.............  Check pollutant code and set index I
            CDAT = TRIM( SEGMENT( 9 ) )     ! pollutant name

C.............  Left justify and convert pollutant name to upper case
            CDAT = ADJUSTL( CDAT ) 
            CALL UPCASE( CDAT ) 

C.............  Look up pollutant name in unique sorted array of
C               Inventory pollutant names
            CIDX  = FINDC( CDAT, NUNIQCAS, UNIQCAS )

C.............  Check to see if data name is in inventory list
            COD  = INDEX1( CDAT, NIPPA, EANAM )

C.............  If pollutant name is not in Inventory Table list
            IF ( CIDX .LE. 0 ) THEN

C.................  Check to see if data name is in list of special names
                CIDX= INDEX1( CDAT, MXSPDAT, SPDATNAM )

C.................  Store status of special data and flag code with
C                   special integer so can ID these records later.
                IF( CIDX .GT. 0 ) THEN
                    SPSTAT( CIDX ) = CIDX 
                    COD = CODFLAG3 + CIDX

                ELSE IF ( CIDX .LE. 0 ) THEN
                    IF( WARNOUT .AND. NWARN( 2 ) .LE. MXWARN ) THEN
                        WRITE( MESG,94010 )
     &                   'WARNING: Skipping pollutant "'// TRIM(CDAT)//
     &                   '" at line', IREC, '- not in Inventory Table'
                        CALL M3MESG( MESG )
                        NWARN( 2 ) = NWARN( 2 ) + 1
                    END IF
                    CYCLE      !  to head of loop

                END IF

C.............  Otherwise, pollutant is in list of Inventory Data Names
            ELSE

C.................  Write warning if pollutant is not kept.  Write only
C                   one time.
               IF( UCASNKEP(CIDX) .LE. 0 ) THEN
                   IF( GETSIZES .AND. WARNKEEP(CIDX) ) THEN 
                       WRITE( MESG,94010 )
     &                   'WARNING: Skipping all lines for pollutant "'//
     &                   TRIM( CDAT )// '" because pollutant is not '//
     &                   'kept by Inventory Table.'
                       CALL M3MESG( MESG )
                   END IF 
                   WARNKEEP( CIDX ) = .FALSE.
                   CYCLE
               END IF

C................  Get Inventory Data SMOKE name from Inventory Table arrays/indices
               CNAM = ITNAMA( SCASIDX( UCASIDX( CIDX ) ) )

C................  Look up SMOKE name in list of annual EI pollutants
               COD = INDEX1( CNAM, NIPPA, EANAM )

C................  Check to ensure that it handles NOI and NONHAP pollutants
C                  while combining VOC + HAPs
               IF( INTGRFLAG ) THEN

C....................  Preventing processing precomputed NONHAP[VOC|TOG]
                   IF( INDEX( CNAM,'NONHAP' ) > 0 ) THEN
                       MESG = 'ERROR: Can NOT process precomputed '// TRIM(CNAM)//
     &                     ' when SMK_PROCESS_HAPS was set to process anuual inventory'
                       CALL M3EXIT( PROGNAME, 0, 0, MESG , 2 )
                   END IF

                   NV = INDEX1( CNAM, MXIDAT, INVDNAM )

                   IF( CINTGR( S ) == 'N' .AND. INVDVTS( NV ) /= 'N' ) THEN
                       PNAM = TRIM( CNAM ) // '_NOI'
                       COD = INDEX1( PNAM, NIPPA, EANAM )

                   ELSE IF( CINTGR( S ) == 'Y' ) THEN
                       L = INDEX( CNAM, ETJOIN )
                       LL= LEN_TRIM( CNAM )
                       PNAM = CNAM
                       IF( L > 0 ) PNAM = CNAM( L+2:LL )
                       IF( PNAM == 'VOC' .OR. PNAM == 'TOG' ) THEN
                           IF( L > 0 ) THEN
                               PNAM = CNAM(1:L+1) // 'NONHAP' //
     &                                CNAM(L+2:LL)
                           ELSE
                               PNAM = 'NONHAP' // TRIM( CNAM )
                           END IF
                           COD = INDEX1( PNAM, NIPPA, EANAM )
                       END IF

                   END IF

               END IF

C................  Check to ensure that the SMOKE intermediate name
C                  set by the Inventory Table is actually in the annual
C                  inventory.  If not, write warning message and cycle.
               IF( COD .LE. 0 ) THEN
                   IF( WARNOUT .AND. NWARN( 5 ) .LE. MXWARN ) THEN
                       WRITE( MESG,94010 )
     &                   'WARNING: Skipping pollutant "'// TRIM(CNAM)//
     &                   '" at line', IREC, '- not in annual inventory.'
                       CALL M3MESG( MESG )
                       NWARN( 5 ) = NWARN( 5 ) + 1
                   END IF
                   CYCLE

C................  If it's found, then record that this pollutant was found
               ELSE
                   EASTAT( COD ) = CIDX
               END IF

            END IF  ! if cidx le 0 or not

C.............  If only getting dates and pollutant information, go 
C               to next loop iteration
            IF( GETSIZES ) CYCLE

C.............  Determine time step pointer based on actual start time
            PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME ) / TDIVIDE + 1

C.............  Skip record if it is out of range of output file
C.............  NOTE - this is only useful if reading only part of data
            IF( PTR. LT. 1 .OR. PTR .GT. NSTEPS ) CYCLE

C.............  Count estimated record count per time step
            DO T = PTR, MIN( PTR + 23, NSTEPS )
                MXPDPT( T ) = MXPDPT( T ) + 1
            END DO

C.............  If only counting records per time step, go to next loop
C               iteration
            IF( GETCOUNT ) CYCLE

C.............  Store source ID
            LPDSRC( S ) = .TRUE.

C.............  Check and set emissions values
            S1 = 15   ! pollutant field start position

            N = 0
            DO J = FSTLOC, LSTLOC
                IF( DAYFLAG ) THEN
                    N = N + 1                   ! N is equal to NFIELD for day-specific
                    TDAT( N,: )  = STR2REAL( SEGMENT( S1-1+J ) )
                ELSE
                    TDAT( :,J )  = STR2REAL( SEGMENT( S1-1+J ) )
                ENDIF 
            END DO

C.............  If available, set total value from hourly file
            TOTAL = 0.
            IF( SFLAG .OR. .NOT. DAYFLAG ) THEN

                IF( SEGMENT( S1-1 ) .NE. ' ' ) THEN
                    TOTAL = STR2REAL( SEGMENT( S1-1 ) )
                    IF( TOTAL < 0.0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Bad line', IREC,
     &                    ': total value "'//TRIM(SEGMENT( S1-1 ))//'"' 
                        CALL M3MESG( MESG )
                        CYCLE  ! to head of read loop
                    END IF
                END IF
            END IF

C.............  Set conversion factor from Inventory Table. Default is
C               1., which is also what is used in all but a handful of
C               special toxics cases.
            CONVFAC = ITFACA( SCASIDX( UCASIDX( CIDX ) ) )

C.............  Record needed data for this source and time step
            DO D = 1, NFIELD

                H = 0
                DO T = PTR, MIN( PTR + 23, NSTEPS )

                    H = H + 1
                    NPDPT( T ) = NPDPT( T ) + 1

                    HS = NPDPT( T )

                    IF( HS .LE. MXPDSRC ) THEN

                        IDXSRC( HS,T ) = HS
                        SPDIDA( HS,T ) = S
                        CIDXA ( HS,T ) = CIDX
                        CODEA ( HS,T ) = COD
                        EMISVA( HS,T ) = CONVFAC * TDAT( D,H )  ! Store data in emissions
                        DYTOTA( HS,T ) = CONVFAC * TOTAL
                    END IF

                END DO
                
                PTR = PTR + H
            
            END DO

        END DO

299     CONTINUE   ! Exit from read loop

C.........  Abort if error found while reading file
        IF( EFLAG ) THEN
            MESG = 'Problem processing day- or hour-specific data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Update output starting date/time and ending date/time
        DFLAG = .TRUE.
        SDATE = RDATE
        STIME = RTIME
        DO I = 1, MINPTR - 1
            CALL NEXTIME( SDATE, STIME, TSTEP )
        END DO

        EDATE = RDATE
        ETIME = RTIME
        
        DO I = 1, MAXPTR - 1
            CALL NEXTIME( EDATE, ETIME, TSTEP )
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( I6.6 )

        END SUBROUTINE RDFF10PD
