
        SUBROUTINE RDORLFR( FDEV, TZONE, TSTEP, MXPDSRC, GETSIZES, 
     &                      GETCOUNT, FIRSTCALL, DAYFLAG, SDATE, STIME, 
     &                      EDATE, ETIME, EASTAT, SPSTAT )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads the day- emissions in ORL FIREEMIS format.
C      It appends the records to the global storage from the MODDAYHR.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 03/06 by B.H. Baek (based on RDEMSPD.F)
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CIFIP, CSOURC, HEATCONTENT, INTGRFLAG, CINTGR

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, UCASNKEP, NUNIQCAS,
     &                      UNIQCAS, NINVTBL, ITNAMA, ITCASA, FIREFLAG,
     &                      UCASIDX, SCASIDX, INVDVTS, MXIDAT, INVDNAM

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NSRC, EANAM, NCHARS, NMAP, MAPNAM,
     &                     MAPFIL

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, NPDPTP, IDXSRC, 
     &                      SPDIDA, CODEA, CIDXA, EMISVA, DYTOTA

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  EXTERNAL FUNCTIONS
C       CHARACTER(2) CRLF
C       LOGICAL      ENVYN
C       INTEGER      ENVINT
C       INTEGER      FIND1
C       INTEGER      FINDC
C       INTEGER      INDEX1
C       INTEGER      JULIAN
C       INTEGER      SECSDIFF
C       INTEGER      STR2INT
C       REAL         STR2REAL
C       INTEGER      YEAR4
        LOGICAL      BLKORCMT
        INTEGER      GETTZONE
C       LOGICAL      SETENVVAR

C        EXTERNAL     BLKORCMT, CRLF, ENVYN, FIND1, FINDC, INDEX1, 
C     &               JULIAN, SECSDIFF, STR2INT, STR2REAL, YEAR4,
C     &               GETTZONE, SETENVVAR, ENVINT
        EXTERNAL     BLKORCMT, GETTZONE

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT( IN ):: FDEV           ! file unit no.
        INTEGER, INTENT( IN ):: TZONE          ! output time zone
        INTEGER, INTENT( IN ):: TSTEP          ! time step HHMMSS
        INTEGER, INTENT( IN ):: MXPDSRC        ! max. day- or hr-specific source
        LOGICAL, INTENT( IN ):: GETSIZES       ! true: get no. time steps & pols
        LOGICAL, INTENT( IN ):: GETCOUNT       ! true: get max no. srcs per time
        LOGICAL, INTENT( IN ):: FIRSTCALL      ! true: first call of a loop
        LOGICAL, INTENT( IN ):: DAYFLAG        ! true: day-specific wildfire data
        INTEGER,INTENT(INOUT):: SDATE          ! Julian starting date in TZONE
        INTEGER,INTENT(INOUT):: STIME          ! start time of data in TZONE
        INTEGER, INTENT(OUT) :: EDATE          ! Julian ending date in TZONE
        INTEGER, INTENT(OUT) :: ETIME          ! ending time of data in TZONE
        INTEGER, INTENT(OUT) :: EASTAT( NIPPA ) ! true: pol/act appears in data
        INTEGER, INTENT(OUT) :: SPSTAT( MXSPDAT ) ! true: special in data

C...........   SUBROUTINE PARAMETERS
        INTEGER      , PARAMETER :: NSEG = 9        ! number of fields for ORL FIREDATA input format
        REAL         , PARAMETER :: TON2LB = 2000.  ! pounds per short ton
        CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

C...........   Local list of bad sources to prevent duplicate writing of error
C              messages
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )

C...........   Local list of whether or not a warning was written for a pollutant or not
        LOGICAL, ALLOCATABLE, SAVE :: LCODWARN ( : )

C...........  Local list for erroneous pollutant names in the file
        INTEGER, SAVE :: NBADPOLS
        CHARACTER(IOVLEN3), ALLOCATABLE, SAVE :: BADPOLS( : )

C...........   Local list of FIPS start/end positions to facilitate
C              faster lookups
        INTEGER, ALLOCATABLE, SAVE :: STARTSRC( : )
        INTEGER, ALLOCATABLE, SAVE :: ENDSRC( : )

C...........   Local list of arrays for warning handling
        LOGICAL, ALLOCATABLE, SAVE :: WARNKEEP( : ) ! true: write warning for Keep = N
        LOGICAL, ALLOCATABLE, SAVE :: WARNMULT( : ) ! true: write warning for Multiple pollutants from a single pollutant

C...........   Temporary read arrays
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line

C...........   Local arrays
        REAL              , ALLOCATABLE, SAVE :: DTACBRN( : )    ! storing acre burned value (acre/day) for computing HFLUX
        REAL              , ALLOCATABLE, SAVE :: DTFUELD( : )    ! storing fuel loading value (tons/acre) for computing HFLUX

        INTEGER           , ALLOCATABLE, SAVE :: NSRCPDDAT( :,: )    ! counting number of sources per day/pollutant
        INTEGER           , ALLOCATABLE, SAVE :: IDXSD    ( : )      ! sorting index for CSRCDAYA
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: CSRCDAYA ( : )      ! unsorted source/day array
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: CSRCDAY  ( : )      ! sorted source/day array

C...........   Other local variables
        INTEGER          H, HS, I, II, J, K, L, LL, N, NV, S, T, V1, V2    ! counters and indices
        INTEGER          L0, L1, L2, L3, L4, L5
        INTEGER          ES, NS, SS       ! end src, tmp no. src, start sourc

        INTEGER          D, SD
        INTEGER       :: N1 = 0
        INTEGER       :: N2 = 0

        INTEGER          CIDX             ! CAS data index
        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER          FIP              ! tmp co/st/cy code
        INTEGER, SAVE :: ICC = 0          ! tmp country code from header
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          JD               ! Julian day number 1...365,366
        INTEGER          JTIME            ! tmp HHMMSS time
        INTEGER          ESTIME           ! tmp HHMMSS episode start time
        INTEGER          EETIME           ! tmp HHMMSS episode end time
        INTEGER       :: LFIP = 0         ! previous st/co FIPS code
        INTEGER, SAVE :: LOOPNO = 0       ! no. of loops
        INTEGER, SAVE :: MAXPTR           ! maximum time step reference pointer
        INTEGER, SAVE :: MINPTR           ! minimum time step reference pointer
        INTEGER          MONTH            ! tmp month number
        INTEGER, SAVE :: MXWARN           !  maximum number of warnings
        INTEGER, SAVE :: NWARN( 5 )       ! warnings counter
        INTEGER, SAVE :: NBADSRC = 0      ! no. bad sources
        INTEGER, SAVE :: NACRBND = 0      ! no. of acres burned var
        INTEGER, SAVE :: NFUELD  = 0      ! no. of fuel loading var
        INTEGER       :: NPOA    = 0      ! unused header number of pol/act
        INTEGER, SAVE :: NSRCDAY = 0      ! no. of src/day combos for computed vars
        INTEGER, SAVE :: NSTEPS  = 0      ! number of time steps
        INTEGER          PTR              ! tmp time step pointer
        INTEGER       :: RDATE = 1980001  ! reference date: Jan 1, 1980
        INTEGER       :: RTIME = 0        ! reference time
        INTEGER, SAVE :: SDATESAV = 0     ! saved start date
        INTEGER, SAVE :: STIMESAV = 0     ! saved start time
        INTEGER, SAVE :: TDIVIDE  = 1     ! time step divisor
        INTEGER          WD               ! tmp field width
        INTEGER          YEAR             ! 4-digit year
        INTEGER       :: YR4 = 0          ! unused header year
        INTEGER          DZONE            ! time shift (ZONE-TZONE)
        INTEGER          ZONE             ! source time zones

        REAL             TDAT             ! temporary data values

        LOGICAL, SAVE :: TFLAG  = .FALSE. ! true: use SCCs for matching with inv
        LOGICAL, SAVE :: DFLAG  = .FALSE. ! true: dates set by data
        LOGICAL       :: EFLAG  = .FALSE. ! TRUE iff ERROR
        LOGICAL       :: WARNOUT= .FALSE. ! true: then output warnings
        LOGICAL, SAVE :: LFLAG = .FALSE.  ! true: output daily/hourly inv in local time
        LOGICAL, SAVE :: PRCHFX = .FALSE. ! true: skip adding HFLUX due to precomputed heat flux
        LOGICAL       :: HFXFLAG= .FALSE. ! true: adding HFLUX into a list
        LOGICAL       :: BNHRFLAG=.FALSE. ! true: adding BEGHOUR into a list
        LOGICAL       :: ENHRFLAG=.FALSE. ! true: adding ENDHOUR into a list
        LOGICAL, SAVE :: FIRSTCOUNT = .TRUE.! true: until after first time routine is called with GETCOUNT=TRUE
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called

        CHARACTER(256) :: BUFFER = ' '    ! src description buffer 
        CHARACTER(300) :: LINE   = ' '    ! line buffer 
        CHARACTER(300) :: MESG   = ' '    ! message buffer

C.........  Temporary local character variables
        CHARACTER(FIPLEN3) CFIP      ! tmp co/st/cy code
        CHARACTER(CASLEN3) CDAT      ! tmp data name (*16)
        CHARACTER(IOVLEN3) CNAM,PNAM ! tmp SMOKE name
        CHARACTER(IOVLEN3) CTMP      ! tmp data name (*16)
        CHARACTER(PLTLEN3) FCID      ! tmp facility ID (*15)
        CHARACTER(CHRLEN3) SKID      ! tmp stack ID (*15) = LocID
        CHARACTER(CHRLEN3) DVID      ! dummy device ID
        CHARACTER(CHRLEN3) PRID      ! dummy process ID
        CHARACTER(SCCLEN3) TSCC      ! tmp source category code (*10)
        CHARACTER(ALLLEN3) CSRC      ! tmp source string
        CHARACTER(ALLLEN3) CSRCD     ! tmp source/date string
        CHARACTER(ALLLEN3) TSRC      ! tmp source string
        CHARACTER( 8 )     DATE      ! tmp date string
        CHARACTER(IOVLEN3) PNAME     ! logical file name for data files

        CHARACTER(16) :: PROGNAME = 'RDORLFR' !  program name

C***********************************************************************
C   begin body of program RDORLFR

C.........  First time routine called
        IF( FIRSTIME ) THEN

C.............  No time zone shift for AERMOD support
            MESG = 'Outputs local time daily and/or hourly inventories (No time shift)'
            LFLAG = ENVYN( 'OUTPUT_LOCAL_TIME', MESG, .FALSE., IOS )

C.............  Allocate memory for storing counting number of sources per day/pollutant
            ALLOCATE( NSRCPDDAT( 366,NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NSRCPDDAT', PROGNAME )
            NSRCPDDAT = 0  ! array

C.............  Allocate memory for bad source storage
            ALLOCATE( BADSRC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

C.............  Allocate memory for bad pollutant issues
            ALLOCATE( LCODWARN( NINVTBL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LCODWARN', PROGNAME )
            ALLOCATE( BADPOLS( NINVTBL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BADPOLS', PROGNAME )
            LCODWARN = .FALSE.
            BADPOLS = ' '

C.............  Create unique list of FIPS codes and other things
            CALL GENUSLST

C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

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

C............  Open I/O API inventory HEATCONTENT file and store to use for 
C              computing HFLUX in PDAY intermediate output file.
            MESG = 'Reading HEATCONTENT data from inventory file...'
            CALL M3MSG2( MESG )

C.............  Open I/O API inventory HEATCONTENT file and store
            ALLOCATE( HEATCONTENT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HEATCONENT', PROGNAME )

            CALL RDMAPPOL( NSRC, 1, 1, 'HEATCONTENT', HEATCONTENT )

            FIRSTIME = .FALSE.
            FIREFLAG = .TRUE.


C.............  Initialize warnings counter
            NWARN = 0  ! array

        END IF  ! End first time subroutine is called

C.........  For the first call in a loop of files, initialize variables
        IF( FIRSTCALL ) THEN
            MINPTR  = 99999999
            MAXPTR  = 0

C.............  Set time step divisor
            TDIVIDE = 3600 * TSTEP / 10000

C.............  If dates have been set by the data, set the number of steps
C               steps
            IF( DFLAG ) THEN
                NSTEPS = 1+ SECSDIFF( SDATE,STIME,EDATE,ETIME )/TDIVIDE
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

C.........  For the second pass of this routine, allocate memory needed for calculating
C           values from day-specific data on the fly.
        IF ( GETCOUNT .AND. FIRSTCOUNT ) THEN

C.............  Determine how much memory is needed for allocating arrays
C               This should be the maximum number of 
C               source/days of fuel load, and acres burned.
            DO D = 1, 366
                DO I = 1, NIPPA
                    IF( EANAM(I)=='FUEL_LOAD' )   N1= N1+ NSRCPDDAT(D,I)
                    IF( EANAM(I)=='ACRESBURNED' ) N2= N2+ NSRCPDDAT(D,I)
                END DO
            END DO
            N = MAX( N1, N2 )

            ALLOCATE( IDXSD( N ), STAT=IOS )          ! Sorting index
            CALL CHECKMEM( IOS, 'IDXSD', PROGNAME )
            ALLOCATE( CSRCDAYA( N ), STAT=IOS )       ! Unsorted SOURCE/DAY combos
            CALL CHECKMEM( IOS, 'CSRCDAYA', PROGNAME )
            ALLOCATE( CSRCDAY( N ), STAT=IOS )        ! Sorted SOURCE/DAY combos
            CALL CHECKMEM( IOS, 'CSRCDAY', PROGNAME )
            ALLOCATE( DTACBRN( N ), STAT=IOS )  ! To store acres burned
            CALL CHECKMEM( IOS, 'DTACBRN', PROGNAME )
            ALLOCATE( DTFUELD( N ), STAT=IOS )  ! To store fuel load
            CALL CHECKMEM( IOS, 'DTFUELD', PROGNAME )

            IDXSD    = 0
            CSRCDAYA = ' '
            CSRCDAY  = ' '
            DTACBRN  = BADVAL3
            DTFUELD  = BADVAL3

            FIRSTCOUNT = .FALSE.

        END IF 

C.........  For the third pass of this routine, create sorted CSRCDAY
C           routine.  This are used for caculating new values (e.g., PMC) from day-specific
C           data on the fly.
        IF ( .NOT. GETSIZES .AND. .NOT. GETCOUNT ) THEN

            CALL SORTIC( NSRCDAY, IDXSD, CSRCDAYA ) 

            DO SD = 1, NSRCDAY
                K = IDXSD( SD )
                CSRCDAY( SD ) = CSRCDAYA( K )
            END DO

        END IF

C.........  Loop through file and read it. In the first section, determine
C           the minimum and maximum date. Use a reference date to do this. In
C           the second section, determine the number of records per time 
C           step. In the third section, read and store the data.  When storing
C           data, time step index is computed from the start date/time instead
C           of the reference date/time so that the indexing will work properly.
        LFIP = 0 
        IREC = 0
        TDAT = 0

        DO         !  Head of period-specific file read loop

C.............  Read first line of file
            READ( FDEV, 93000, END=299 ) LINE
            IREC = IREC + 1

            L = LEN_TRIM( LINE )

C.............  Skip blank lines 
            IF( L == 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set
C               properly
            CALL GETHDR( 1, .TRUE., .FALSE., .FALSE., LINE, ICC, YR4,
     &                   NPOA, IOS )

C.............  Interpret error status
            IF( IOS .GT. 0 ) EFLAG = .TRUE.

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Get lines
            CALL PARSLINE( LINE, NSEG, SEGMENT )

C.............  Use the file format definition to parse the line into
C               the various data fields
            CFIP = REPEAT( '0', FIPLEN3 )
            WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC  ! country code of FIPS
            CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )  ! state/county code
            FIP    = STR2INT( CFIP )          ! FIP codes
            FCID   = SEGMENT( 2 )   ! fire ID
            SKID   = SEGMENT( 3 )   ! location ID
            TSCC   = SEGMENT( 4 )   ! SCC
            DVID   = ' '            ! dummy device id
            PRID   = ' '            ! dummy process id
            CDAT   = SEGMENT( 5 )   ! Pollutants(FUEL_LOAD, ACRESBURNED,,,)
            DATE   = SEGMENT( 6 )   ! Date of episode
            ESTIME = STR2INT( SEGMENT( 8 ) ) * 10000 ! episode start time
            EETIME = STR2INT( SEGMENT( 9 ) ) * 10000 ! episode end time

            IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )

C.............  Conver CAS number to pollutant names if available
            I = INDEX1( CDAT, NINVTBL, ITCASA )
            IF( I > 0 ) CDAT = ITNAMA( I )

C............. Check fire beginning and ending time format and print warning if necessary
            IF( EETIME > 230000 .OR. EETIME < 0 ) THEN
                MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //
     &                 ' Date: ' // DATE  // ' :: Fire ending' //
     &                 ' time not in the acceptable range of 0 to 23'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

            IF( ESTIME > 230000 .OR. ESTIME < 0  ) THEN
                MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //
     &                 ' Date: ' // DATE  // ' :: Fire starting' //
     &                 ' time not in the acceptable range of 0 to 23'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

            IF( ESTIME > EETIME ) THEN
                MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //
     &                 ' Date: ' // DATE  // ' :: Fire ending time' //
     &                 ' can not be earlier then a begining time'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

C.............  Check date format
            IF( DATE( 3:3 ) /= '/' .OR.
     &          DATE( 6:6 ) /= '/'      ) THEN
                MESG = 'ERROR: Incorrect date format ( MM/DD/YY ) :' //
     &                 ' Region: '// CFIP // ' SCC: ' // TSCC // 
     &                 ' Date: ' // DATE
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

C.............  Check and Set emissions values
            TDAT = STR2REAL( SEGMENT( 7 ) )     ! Day-specific total emission

            IF ( TDAT .LT. 0.0 )  THEN
                EFLAG = .TRUE.
                WRITE( MESG,94030 ) 'ERROR: Bad data value "',
     &              TDAT, '" of ' // TRIM( CDAT ) //
     &              ' for region ' // CFIP //
     &              ' and SCC ' // TSCC // ' on date '// DATE
                CALL M3MESG( MESG )
                CYCLE  ! to head of read loop
            END IF

C.............  Counting the number of times precomputed HFLUX
C               values appear in the input file
            IF( GETSIZES .AND. CDAT == 'HFLUX' ) PRCHFX = .TRUE.

C.............  Counting and adding HFLUX, BEGHOUR, and ENDHOUR
C               building a list of source characteristics and store
            IF( GETCOUNT .AND. .NOT. PRCHFX ) THEN

                IF( HFXFLAG  ) CDAT = 'HFLUX'
                IF( BNHRFLAG ) CDAT = 'BEGHOUR'
                IF( ENHRFLAG ) CDAT = 'ENDHOUR'

                HFXFLAG = .FALSE.
                IF( CDAT == 'ACRESBURNED' ) THEN
                    NACRBND = NACRBND + 1
                    HFXFLAG = .TRUE.    ! indicating adding HFLUX
                    BACKSPACE( FDEV )
                END IF

                BNHRFLAG = .FALSE.
                IF( CDAT == 'HFLUX' ) THEN
                    BNHRFLAG = .TRUE.    ! indicating adding BEGHOUR
                    BACKSPACE( FDEV )
                END IF

                ENHRFLAG = .FALSE.
                IF( CDAT == 'BEGHOUR' ) THEN
                    ENHRFLAG = .TRUE.    ! indicating adding ENDHOUR
                    BACKSPACE( FDEV )
                END IF

                IF( CDAT == 'FUEL_LOAD' ) NFUELD = NFUELD + 1

            END IF

C.............  Set Julian day from MMDDYY8 SAS format
            MONTH = STR2INT( DATE( 1:2 ) )
            DAY   = STR2INT( DATE( 4:5 ) )
            YEAR  = YEAR4( STR2INT( DATE( 7:8 ) ) )

            JD = JULIAN( YEAR, MONTH, DAY )
            JDATE = 1000 * YEAR + JD
            JTIME = 0


C.............  Local time shift flag (county-specific)
            IF( .NOT. LFLAG ) THEN
            
C.................  Set time zone number
                ZONE = GETTZONE( CFIP )
            
C.................  If daily emissions are not in the output time zone, print 
C                   warning
                IF( GETCOUNT ) THEN
                   IF( WARNOUT .AND. ( ZONE .NE. TZONE ) .AND. 
     &                ( NWARN( 1 ) .LE. MXWARN )               ) THEN
                       WRITE( MESG,94010 ) 'WARNING: Time zone ', ZONE, 
     &                   'in day-specific file at line of pollutant ' //
     &                   TRIM( CDAT ) // ' on ' // TRIM( DATE ) // 
     &                   ' does not match output time zone', TZONE
                       CALL M3MESG( MESG )
                       NWARN( 1 ) = NWARN( 1 ) + 1
                   END IF
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

C.............  Store maximum time step number as compared to rference
            IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

C.............  If FIPS code is not the same as last time, then
C               look it up and get indidies
            IF( FIP .NE. LFIP ) THEN
                J = FINDC( CFIP, NINVIFIP, INVCFIP )
                IF( J .LE. 0 ) THEN
                    WRITE( MESG,93000 ) 'INTERNAL ERROR: Could not '//
     &                     'find FIPS code ' // CFIP // 'in internal list.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                END IF

                SS = STARTSRC( J )
                ES = ENDSRC( J )
                NS = ES - SS + 1
                LFIP = FIP

            END IF

C.............  If SCCs are needed for matching...
            IF ( TFLAG ) THEN

C.................  Build source characteristics field for searching inventory
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                        TSCC, CHRBLNK3, POLBLNK3, CSRC )

C.................  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )

C.............  If SCCs are not being used for matching (at least not yet)...
            ELSE

C.................  Build source characteristics field for searching inventory
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                        TSCC, CHRBLNK3, POLBLNK3, CSRC )

C.................  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )

C.................  If source is not found for day-specific processing, see 
C                   if reading the SCC in helps (needed for IDA format)
                IF( J .LE. 0 ) THEN

C.....................  Build source characteristics field for searching inventory
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )

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

C.............  Look up pollutant name in unique sorted array of
C               Inventory pollutant names
            CDAT  = SEGMENT( 5 )
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

C................  If not in list of special names, check to see
C                  if it's a SMOKE pollutant name (intermediate name)
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
               IF( UCASNKEP(CIDX) .LE. 0 .AND. WARNKEEP(CIDX) ) THEN
                   WARNKEEP( CIDX ) = .FALSE.
                   IF( GETSIZES ) THEN
                       WRITE( MESG,94010 )
     &                   'WARNING: Skipping all lines for pollutant "'//
     &                   TRIM( CDAT )// '" because pollutant is not '//
     &                   'kept by Inventory Table.'
                       CALL M3MESG( MESG )
                   END IF
                   CYCLE
               ELSE IF ( UCASNKEP(CIDX) .GT. 1 .AND.
     &                   WARNMULT(CIDX)              ) THEN
                   WARNMULT( CIDX ) = .FALSE.
                   IF( GETSIZES ) THEN
                       WRITE( MESG,94010 )
     &                   'WARNING: Skipping all lines for pollutant "'//
     &                   TRIM( CDAT )// '" because Inventory Table '//
     &                   'splits it into',UCASNKEP(CIDX),'pollutants.'//
     &                   CRLF()//BLANK10//'The SMOKE code needs to '//
     &                   'be enhanced to support this approach for '//
     &                   'day- and hour-specific data.'
                       CALL M3MESG( MESG )
                   END IF
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

C.............  Count the number of sources per day & pollutant/variable
C.............  This will give us how many source/date combos there are for 
C               any variables, including HFLUX
            NSRCPDDAT( JD, COD ) = NSRCPDDAT( JD, COD ) + 1
            
C.............  If only getting dates and pollutant information, go 
C               to next loop iteration
            IF( GETSIZES ) CYCLE

C.............  Determine time step pointer based on actual start time
            PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME )/TDIVIDE + 1

C.............  Skip record if it is out of range of output file
C.............  NOTE - this is only useful if reading only part of data
            IF( PTR. LT. 1 .OR. PTR .GT. NSTEPS ) CYCLE
            
C.............  Count estimated record count per time step
            DO T = PTR, MIN( PTR + 23, NSTEPS )
                MXPDPT( T ) = MXPDPT( T ) + 1
            END DO

C.............  Store variable values.  Only need to do this on the the second
C               pass.  Need to do this before the third pass through the data because
C               that is when the calculation is made.           
            IF( GETCOUNT .AND. .NOT. PRCHFX ) THEN          ! No precomputed formula/heat flux
                IF( ( HFXFLAG .OR. CDAT == 'FUEL_LOAD' ) ) THEN  ! Acres burned value or fuel load value

C.....................  Figure out which source/day this is for storing in correct source/day
C.....................  This code does *not* assume that the data have been sorted first.
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                            TSCC, DATE, POLBLNK3, CSRCD )

C.....................  Build unsorted arrays of source/days and emissions for calculating formula
                    SD = 0
                    IF( NSRCDAY > 0 ) THEN
                        SD= INDEX1( CSRCD, NSRCDAY, CSRCDAYA )
                    END IF

                    IF( SD <= 0 ) THEN
                        NSRCDAY = NSRCDAY + 1
                        SD = NSRCDAY
                        CSRCDAYA( SD ) = CSRCD
                        IDXSD   ( SD ) = SD
                    END IF

                    IF( HFXFLAG ) DTACBRN( SD ) = TDAT        ! storing acres burned
                    IF( CDAT == 'FUEL_LOAD' ) DTFUELD( SD ) = TDAT ! storing fuel load

                END IF
            END IF    ! Second pass only
            
C.............  If only counting records per time step, go to next loop
C               iteration
            IF( GETCOUNT ) CYCLE

C.............  Store source ID
            LPDSRC( S ) = .TRUE.

C.............  Computing HFLUX, BEGHOUR, ENDHOUR (as a default)
            IF( CDAT == 'HFLUX' .AND. .NOT. PRCHFX ) THEN

                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, TSCC, DATE,
     &                        POLBLNK3, CSRCD )

C.................  Lookup source/date string in master list to get position
                SD = FINDC( CSRCD, NSRCDAY, CSRCDAY )
                K  = IDXSD( SD )

                IF( DTFUELD( K ) < AMISS3 ) THEN
                    LL = LEN_TRIM( CSRCD )
                    CALL FMTCSRC( CSRCD, 6, BUFFER, L2 )

                    MESG = 'ERROR: Missing value of '//
     &                     'fuel load for source:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) // 
     &                     ' on date ' // CSRCD( LL-7: LL )
                    CALL M3MSG2( MESG )
                    EFLAG = .TRUE.

                ELSE 
C.....................  Compute Heat Flux value
                    TDAT = DTACBRN( K ) * DTFUELD( K ) * ! computing HFLUX (BTU/day)
     &                     HEATCONTENT( S ) * TON2LB     ! HEATCONTENT(BTU/lb)=8000

                END IF

            END IF

            IF( CDAT == 'BEGHOUR' ) TDAT = REAL( ESTIME )  ! storing BEGHOUR
            IF( CDAT == 'ENDHOUR' ) TDAT = REAL( EETIME )  ! storing ENDHOUR

C.............  Record needed data for this source and time step
            H = 0
            DO T = PTR, MIN( PTR + 23, NSTEPS )
                H = H + 1
                NPDPT( T ) = NPDPT( T ) + 1
                NPDPTP( T,COD ) = NPDPTP( T,COD ) + 1

                HS = NPDPT( T )

                IF( HS .LE. MXPDSRC ) THEN

                    IDXSRC( HS,T ) = HS
                    SPDIDA( HS,T ) = S
                    CODEA ( HS,T ) = COD
                    CIDXA ( HS,T ) = CIDX
                    EMISVA( HS,T ) = TDAT  ! Store data in emissions
                    DYTOTA( HS,T ) = TDAT

                END IF
            END DO

        END DO     ! Main read loop of day-specific data

299     CONTINUE   ! Exit from read loop

C.........  Warning messages for HFLUX 
        IF( GETCOUNT ) THEN

            IF( PRCHFX ) THEN
                MESG = 'WARNING: Skipping internal heat flux '//
     &                 'computation due to the existence of '//
     &                 'precomputed HFLUX in PTDAY file'
                CALL M3MSG2( MESG )
            END IF

            IF( NACRBND .NE. NFUELD .AND. .NOT. PRCHFX ) THEN
                MESG = 'ERROR: No of ACRESBURNED and FUEL_LOAD are' //
     &                 ' not matched for heat flux computation.'
                CALL M3MSG2( MESG )
            END IF

            IF( NACRBND < 1 .AND. .NOT. PRCHFX ) THEN
                MESG = 'FATAL ERROR: No ACRESBURNED data are available'
     &                 // ' for internal heat flux computation.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( NFUELD < 1 .AND. .NOT. PRCHFX ) THEN
                MESG = 'FATAL ERROR: No FUEL_LOAD data are available'
     &                 // ' for internal heat flux computation.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

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
94030   FORMAT( 10( A, :, F15.0, :, 1X ) )

94020   FORMAT( I6.6 )

        END SUBROUTINE RDORLFR
