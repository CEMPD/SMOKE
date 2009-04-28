
        SUBROUTINE RDEPSPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, 
     &                      GETSIZES, GETCOUNT, FIRSTCALL, DAYFLAG, 
     &                      SDATE, STIME, EDATE, ETIME, EASTAT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDCSRC, CHECKMEM 
C      Functions: I/O API functions, GETFLINE, YR2DAY
C
C  REVISION  HISTORY:
C
C****************************************************************************
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

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: MXIDAT, INVDCOD, INVDNAM

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, USEDAYLT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NSRC, EANAM, NCHARS

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, IDXSRC, SPDIDA,
     &                      CODEA, EMISVA

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL                CHKINT
        LOGICAL                CHKREAL
        CHARACTER(2)           CRLF
        LOGICAL                ENVYN
        INTEGER                FIND1
        INTEGER                FINDC   !  returns -1 for failure
        INTEGER                GETTZONE
        INTEGER                INDEX1
        LOGICAL                ISDSTIME
        INTEGER                JULIAN 
        INTEGER                SECSDIFF
        INTEGER                STR2INT
        REAL                   STR2REAL
        INTEGER                YEAR4  

        EXTERNAL CHKINT, CHKREAL, CRLF, ENVYN, FIND1, FINDC, GETTZONE,
     &           INDEX1, ISDSTIME, JULIAN, SECSDIFF, STR2INT, 
     &           STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV           ! input file unit no.
        INTEGER, INTENT (IN) :: TZONE          ! output time zone
        INTEGER, INTENT (IN) :: INSTEP         ! expected data time step HHMMSS
        INTEGER, INTENT (IN) :: OUTSTEP        ! output time step HHMMSS
        INTEGER, INTENT (IN) :: MXPDSRC        ! max. day- or hr-specific source
        LOGICAL, INTENT (IN) :: GETSIZES       ! true: get no. time steps & pols
        LOGICAL, INTENT (IN) :: GETCOUNT       ! true: get max no. srcs per time
        LOGICAL, INTENT (IN) :: FIRSTCALL      ! true: first call of a loop
        LOGICAL, INTENT (IN) :: DAYFLAG        ! true: day-, false: hour-spec
        INTEGER,INTENT(INOUT):: SDATE          ! Julian starting date in TZONE
        INTEGER,INTENT(INOUT):: STIME          ! start time of data in TZONE
        INTEGER, INTENT(OUT) :: EDATE          ! Julian ending date in TZONE
        INTEGER, INTENT(OUT) :: ETIME          ! ending time of data in TZONE
        LOGICAL, INTENT(OUT) :: EASTAT( NIPPA ) ! true: pol/act appears in data

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: SRTDCOD( : )
        INTEGER, ALLOCATABLE, SAVE :: SRTDIDX( : )

C.........  Local list of bad sources to prevent duplicate writing of error
C           messages
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )
        
C.........  Temporary variables for storing source characteristics.  These
C           variables must be the width of the fields for global source
C           characteristics definition for use in BLDCSRC.
        CHARACTER(PLTLEN3) FCID  ! tmp plant ID
        CHARACTER(CHRLEN3) PTID  ! tmp point ID
        CHARACTER(CHRLEN3) SKID  ! tmp stack ID
        CHARACTER(CHRLEN3) SGID  ! tmp segment ID

C...........   Other local variables

        REAL             EDIV    !  emissions divisor
        REAL             EMIS    !  tmp emission value
        REAL             RSTEPS  !  number of time steps
        REAL             TFAC    !  premultiplied time step factor

        INTEGER          HS, I, J, K, L, L2, S, T  ! counters and indices

        INTEGER          COD     !  Temporary pollutant code number
        INTEGER          DD      !  tmp day 
        INTEGER          EDT     !  tmp end date
        INTEGER          ETM     !  tmp end time
        INTEGER          FIP     !  tmp state/county code
        INTEGER          HH      !  tmp hour
        INTEGER, SAVE :: ICC     !  position of CNTRY in CTRYNAM
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER, SAVE :: INY     !  tmp emissions year
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          JTIME            ! tmp HHMMSS time
        INTEGER          LDATE            ! previous Julian date
        INTEGER, SAVE :: LOOPNO = 0       ! no. of loops
        INTEGER, SAVE :: MAXPTR  !  maximum time step reference pointer
        INTEGER, SAVE :: MINPTR  !  minimum time step reference pointer
        INTEGER          MM      !  tmp month
        INTEGER, SAVE :: MXWARN  !  maximum number of warnings
        INTEGER, SAVE :: NBADSRC = 0      ! no. bad sources
        INTEGER          NS      !  tmp no. record's time steps
        INTEGER, SAVE :: NWARN =0!  number of warnings in this routine
        INTEGER, SAVE :: NSRCPOL = 0 ! cumulative source x data variables
        INTEGER, SAVE :: NSTEPS = 0       ! number of time steps
        INTEGER          PTR, PTR1, PTR2  ! tmp time step pointers
        INTEGER       :: RDATE = 1980001  ! reference date: Jan 1, 1980
        INTEGER       :: RTIME = 0        ! reference time
        INTEGER, SAVE :: SDATESAV = 0     ! saved start date
        INTEGER, SAVE :: STIMESAV = 0     ! saved start time
        INTEGER          SDT     !  tmp start date
        INTEGER          SEGMENT !  tmp integer segment number
        INTEGER          STACK   !  tmp integer stack ID
        INTEGER          STM     !  tmp start time
        INTEGER, SAVE :: TDIVIDE  = 1     ! time step divisor
        INTEGER          YY               ! 2-digit year
        INTEGER          ZONE             ! tmp source time zone

        LOGICAL       :: DAYLIT   = .FALSE.  ! true: date in daylight time
        LOGICAL, SAVE :: DFLAG    = .FALSE.  ! true: dates set by data
        LOGICAL       :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL       :: WARNOUT  = .FALSE.  ! true: then output warnings

        CHARACTER(2)       TMPAA !  tmp time period code
        CHARACTER(100)  :: BUFFER = ' '     ! src description buffer 
        CHARACTER(300)     LINE  !  Input line from POINT file
        CHARACTER(300)     MESG  !  Text for M3EXIT()
        CHARACTER(IOVLEN3) CPOL  !  Temporary pollutant code
        CHARACTER(FIPLEN3) CFIP  !  Character FIP code
        CHARACTER(SCCLEN3) TSCC  !  Temporary character SCC
        CHARACTER(CHRLEN3) CHAR4     ! tmp plant characteristic 4
        CHARACTER(ALLLEN3) CSRC      ! tmp source string

        CHARACTER(16) :: PROGNAME = 'RDEPSPD' ! Program name

C***********************************************************************
C   begin body of subroutine RDEPSPD

C.........  Set up settings the first time the subroutine is called
        IF( FIRSTIME ) THEN

C.............  Create sorted data codes
            ALLOCATE( SRTDIDX( MXIDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRTDIDX', PROGNAME )
            ALLOCATE( SRTDCOD( MXIDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRTDCOD', PROGNAME )

            DO I = 1, MXIDAT
                SRTDIDX( I ) = I
            END DO

            CALL SORTI1( MXIDAT, SRTDIDX, INVDCOD )

            DO I = 1, MXIDAT
                J = SRTDIDX( I )
                SRTDCOD( I ) = INVDCOD( J )
            END DO

C.............  Allocate memory for bad source storage
            ALLOCATE( BADSRC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

            FIRSTIME = .FALSE.

        END IF

C.........  For the first call in a loop of files, initialize variables
        IF( FIRSTCALL ) THEN
            MINPTR  = 99999999
            MAXPTR  = 0

C.............  Set time step divisor
            TDIVIDE = 3600 * OUTSTEP / 10000
            TFAC    = REAL( OUTSTEP ) / REAL( INSTEP * TDIVIDE )

C.............  If dates have been set by the data, set number of steps
            IF( DFLAG ) THEN
                NSTEPS = 1+ SECSDIFF( SDATE,STIME,EDATE,ETIME )/ TDIVIDE
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

        END IF

C........................................................................
C.............  Head of the FDEV-read loop  .............................
C........................................................................

C.........  Loop through file and read it. In the first section, determine
C           the minimum and maximum date. Use a reference date to do this. In
C           the second section, determine the number of records per time 
C           step. In the third section, read and store the data.  When storing
C           data, time step index is computed from the start date/time instead
C           of the reference date/time so that the indexing will work properly.
        IREC = 0
        ICC  = -9
        DO         !  Head of period-specific file read loop

C.............  Read a line of input file and check input status

            READ( FDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS, 
     &                 'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            L = LEN_TRIM( LINE )  ! store width of line and check

C.............  Skip blank lines
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set 
C               properly.  Note that data value names are not read.
            CALL GETHDR( MXIDAT, .TRUE., .TRUE., .FALSE., 
     &                   LINE, ICC, INY, I, IOS )

C.............  Interpret error status
            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
            END IF                

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Check and set time period type (Year/day/hourly)
            TMPAA = ADJUSTL( LINE( 57:58 ) )
            CALL UPCASE( TMPAA )

C.............  Emissions are not over a special interval, skip line
            IF( TMPAA( 1:1 ) .NE. 'S' ) CYCLE

C.............  Set pollutant name
            CPOL = ADJUSTL( LINE( 156:160 ) )
            L    = LEN_TRIM( CPOL )

C.............  Ensure that pollutant field is not blank
C.............  If it is blank, move on to next record
            IF( CPOL .EQ. ' ' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Pollutant field ' //
     &                 'is missing at line', IREC
                CALL M3MESG( MESG )
                CYCLE

C.............  Determine if pollutant ID is supplied as integer or string
C.............  If supplied as integer, look up name
            ELSE IF( CHKINT( CPOL ) ) THEN
                COD = STR2INT( CPOL )
                J = FIND1( COD, MXIDAT, SRTDCOD )

C.................  Error if pollutant ID is not in master list
                IF( J .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'ERROR: Source dropped: ' //
     &                     'pollutant number "' // CPOL( 1:L ) //
     &                     '" at line', IREC, CRLF() // BLANK10 //
     &                     'of emission file is not in master ' //
     &                     'pollutants list'
                    CALL M3MESG( MESG )
                    CYCLE

C.................  Store pollutant position in master list and name
                ELSE
                    COD  = SRTDIDX( J )            ! position in original list
                    CPOL = INVDNAM( COD )
                END IF

C.............  Otherwise, make sure name is in list of valid pollutants
            ELSE
                COD = INDEX1( CPOL, MXIDAT, INVDNAM )

C..................  Error if pollutant name is not in master list
                IF( COD .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'ERROR: Source dropped: ' //
     &                     'pollutant "' // CPOL( 1:L ) //
     &                     '" at line', IREC, CRLF() // BLANK10 //
     &                     'of emission file is not in master ' //
     &                     'pollutants list'
                    CALL M3MESG( MESG )
                    CYCLE

                END IF 

            END IF

C.............  Check pollutant code for being in inventory files
            COD  = INDEX1( CPOL, NIPPA, EANAM )

            IF ( COD .LE. 0 ) THEN

                IF( WARNOUT ) THEN
                    WRITE( MESG,94010 ) 
     &                 'WARNING: Skipping pollutant "'// CPOL( 1:L )//
     &                 '" at line', IREC, '- not in inventory'
                    CALL M3MESG( MESG )
                END IF

                CYCLE      !  to head of loop
            END IF

            EASTAT( COD ) = .TRUE.

C.............  Check state/county codes, error for missing
            IF( .NOT. CHKINT( LINE( 12:16 ) ) .OR.
     &          LINE( 12:16 ) .EQ. ' '             ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: State and county ' //
     &                 'code is non-integer or missing at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Retrieve country, state, county code
            FIP  = ICC * 100000 + STR2INT( LINE( 12:16 ) )

C.............  Determine local time zone (non-daylight)
            ZONE = GETTZONE( FIP )

C.............  Set Julian day from file, and set time 
            SDT  = STR2INT( LINE( 60:67 ) )
            STIME = 10000 * MOD( SDT,100 ) 
            SDT   = SDT / 100
            DD    = MOD( SDT , 100 )           ! Set 2-digit day
            MM    = MOD( SDT / 100, 100 )      ! Set 2-digit month
            YY    = SDT/10000                  ! Set 2-digit year 

            YY    = YEAR4( YY )                ! Convert to 4-digit year
            SDT   = 1000 * YY + JULIAN( YY, MM, DD )
            SDATE = SDT

            EDT  = STR2INT( LINE( 69:76 ) )
            ETIME = 10000 * MOD( EDT , 100 )
            EDT   = EDT / 100
            DD    = MOD( EDT , 100 )
            MM    = MOD( EDT / 100, 100 )
            YY    = EDT/10000

            YY    = YEAR4( YY )                ! Convert to 4-digit year
            EDATE = 1000 * YY + JULIAN( YY, MM, DD )

C.............  Determine the float number of time steps (time step set as 
C               24 or 1 by calling routine)
C.............  NOTE - because the ending time as defined by EPS input files
C               is set as the time the emissions apply until, there is no
C               need for a +1 in the following formula.
            RSTEPS = SECSDIFF( SDATE,STIME,EDATE,ETIME ) * TFAC

C.............  Check validity of dates/times for current call of routine...
C.............  Error if partial time step data found
            NS     = INT( RSTEPS )
            IF( MOD( RSTEPS, REAL( NS ) ) .GT. 0. ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Data record at line', IREC,
     &                 'contains less than 1 time step of data.'
                CALL M3MESG( MESG )
                CYCLE

C.............  Warning if hourly data used for more than one day in a single
C               record
            ELSE IF ( INSTEP .EQ. 10000 .AND. RSTEPS .GT. 24. ) THEN
                WRITE( MESG,94010 ) 'WARNING: Data record at line', 
     &                 IREC, 'has hourly data for more than one day'
                CALL M3MESG( MESG )

            END IF

C.............  Set divisor for special time period to determine if hourly,
C               daily, or part-day data
            EDIV = 1. / RSTEPS

C.............  Determine if this FIPS code is exempt from daylight time
            K = FIND1( FIP, NCOUNTY, CNTYCOD )
            
C.............  Loop through time steps for this record
            JDATE  = SDATE
            JTIME  = STIME
            LDATE  = -9
            DAYLIT = .FALSE.
            DO T = 1, NS

C.................  Convert the time zone for this record into local time
                IF( JDATE .NE. LDATE ) THEN

C.....................  Check if date is in daylight time, if local zone has
C                       already been converted, and if this FIPS code is
C                       exempt from daylight time or not.
                    IF(       ISDSTIME( JDATE ) .AND. 
     &                  .NOT. DAYLIT            .AND.
     &                  USEDAYLT( K )                   ) THEN
                        DAYLIT = .TRUE.
                        ZONE = ZONE - 1

                    ELSE IF( .NOT. ISDSTIME( JDATE ) .AND. 
     &                             DAYLIT            .AND.
     &                       USEDAYLT( K )                  ) THEN
                        DAYLIT = .FALSE.
                        ZONE = ZONE + 1
                    END IF

C.....................  Convert date and time to output time zone
                    CALL NEXTIME( JDATE, JTIME, ( ZONE-TZONE ) * 10000 )

C.....................  Update previous iteration date for next iteration
                    LDATE = JDATE

                END IF

C.................  Determine time step pointer based on reference time
                PTR = SECSDIFF( RDATE,RTIME,JDATE,JTIME ) / TDIVIDE + 1

C.................  Store minimum time step number as compared to reference
                IF( PTR .LT. MINPTR ) MINPTR = PTR

C.................  Store maximum time step number as compared to reference
                IF( PTR .GT. MAXPTR ) MAXPTR = PTR

C.................  If not only getting dates and pollutant information
                IF( .NOT. GETSIZES ) THEN

C.....................  Determine time step index based on actual time
                    PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME ) / 
     &                    TDIVIDE + 1

C.....................  Count estimated record count per time step
                    MXPDPT( PTR ) = MXPDPT( PTR ) + 1

C.....................  Store first and last pointers in loop
                    IF( T .EQ. 1  ) PTR1 = PTR
                    IF( T .EQ. NS ) PTR2 = PTR

                END IF

C.................  Increment time step
                CALL NEXTIME( JDATE, JTIME, OUTSTEP )

            END DO

C.............  If only getting dates and pollutant information, go 
C               to next loop iteration
C.............  If only counting records per time step, go to next loop
C               iteration
            IF( GETSIZES .OR. GETCOUNT ) CYCLE

C.............  Check emission value
            IF( .NOT. CHKREAL( LINE( 162:171 ) ) ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( CPOL )
                WRITE( MESG,94010 ) 'ERROR: Emissions data value for '//
     &                 CPOL( 1:L ) // ' is not a number ' // CRLF()//
     &                 BLANK10 // 'or has bad formatting at line', IREC
                CALL M3MESG( MESG )

            ELSE IF( LINE( 162:171 ) .EQ. ' ' ) THEN
                L = LEN_TRIM( CPOL )
                WRITE( MESG,94010 ) 'WARNING: Emissions data for ' //
     &                 CPOL( 1:L ) // ' are missing at line', IREC
                CALL M3MESG( MESG )
                LINE( 162:171 ) = '0.'  ! to prevent STR2REAL warning

            END IF

C.............  If there has been an error, do not try to store any of the
C               records.  Instead  go to next line of file.
            IF( EFLAG ) CYCLE
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields...

            FCID = ADJUSTL ( LINE(  40: 44 ) )  ! plant ID
            SKID = ADJUSTL ( LINE(  45: 48 ) )  ! stack ID
            PTID = ADJUSTL ( LINE(  50: 52 ) )  ! point ID
            SGID = ADJUSTL ( LINE(  54: 55 ) )  ! segment ID
            TSCC = ADJUSTL ( LINE(  29: 38 ) )  ! SCC code
            EMIS = STR2REAL( LINE( 162:171 ) )

C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Convert stack and segment to integer and then character to
C               remove leading zeros
            IF( CHKINT( SKID ) .AND. SKID .NE. ' ' ) THEN
                STACK   = STR2INT( SKID )
                WRITE( SKID, '(I4)' ) STACK
            END IF

            IF( CHKINT( SGID ) .AND. SGID .NE. ' ' ) THEN
                SEGMENT = STR2INT( SGID )
                WRITE( SGID, '(I2)' ) SEGMENT
            END IF

            CHAR4 = TSCC
            CALL BLDCSRC( CFIP, FCID, SKID, PTID, SGID, 
     &                    CHAR4, CHRBLNK3, POLBLNK3, CSRC )

C.............  Search for this record in sources
            S = FINDC( CSRC, NSRC, CSOURC )

C.............  Store source in list of bad sources
C.............  Print warning about sources not found in the inventory
            IF( S .LE. 0 ) THEN

C.................  Search for source in list of bad sources
                S = INDEX1( CSRC, NBADSRC, BADSRC )

                IF( WARNOUT .AND. S .LE. 0 ) THEN                

                    NBADSRC = NBADSRC + 1
                    BADSRC( NBADSRC ) = CSRC

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                    MESG = 'WARNING: Source will be dropped since it '//
     &                     'is not in the inventory:' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                END IF

                CYCLE               !  to head of read loop

C.............  Otherwise, update master list of sources in the inventory
            ELSE
                LPDSRC( S ) = .TRUE.

            END IF

            DO T = PTR1, PTR2

                NPDPT( T ) = NPDPT( T ) + 1

                HS = NPDPT( T )

                IF( HS .LE. MXPDSRC ) THEN

                    IDXSRC( HS,T ) = HS
                    SPDIDA( HS,T ) = S
                    CODEA ( HS,T ) = COD
                    EMISVA( HS,T ) = EMIS * EDIV

                END IF

            END DO      !  end of loop on time steps for current record

        END DO          !  to head of FDEV-read loop

199     CONTINUE        !  end of the FDEV-read loop

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
            CALL NEXTIME( SDATE, STIME, OUTSTEP )
        END DO

        EDATE = RDATE
        ETIME = RTIME
        DO I = 1, MAXPTR - 1
            CALL NEXTIME( EDATE, ETIME, OUTSTEP )
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDEPSPD
