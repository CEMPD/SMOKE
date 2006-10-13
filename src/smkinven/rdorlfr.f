
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is the inventory arrays
        USE MODSOURC, ONLY: IFIP, CSOURC, HEATCONTENT

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP, UCASNKEP, NUNIQCAS,
     &                      UNIQCAS, NINVTBL, ITNAMA, ITCASA, FIREFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NSRC, EANAM, NCHARS, NMAP, MAPNAM,
     &                     MAPFIL

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, NPDPTP, IDXSRC, 
     &                      SPDIDA, CODEA, EMISVA, DYTOTA

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        LOGICAL      ENVYN
        INTEGER      ENVINT
        INTEGER      FIND1
        INTEGER      FINDC
        INTEGER      INDEX1
        INTEGER      JULIAN
        INTEGER      SECSDIFF
        INTEGER      STR2INT
        REAL         STR2REAL
        INTEGER      YEAR4
        LOGICAL      BLKORCMT
        INTEGER      GETTZONE
        LOGICAL      SETENVVAR


        EXTERNAL     BLKORCMT, CRLF, ENVYN, FIND1, FINDC, INDEX1, 
     &               JULIAN, SECSDIFF, STR2INT, STR2REAL, YEAR4,
     &               GETTZONE, SETENVVAR, ENVINT

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT( IN ):: FDEV           ! file unit no.
        INTEGER, INTENT( IN ):: TZONE          ! output time zone
        INTEGER, INTENT( IN ):: TSTEP          ! time step HHMMSS
        INTEGER, INTENT( IN ):: MXPDSRC        ! max. day- or hr-specific source
        LOGICAL, INTENT( IN ):: GETSIZES       ! true: get no. time steps & pols
        LOGICAL, INTENT( IN ):: GETCOUNT       ! true: get max no. srcs per time
        LOGICAL, INTENT( IN ):: FIRSTCALL      ! true: first call of a loop
        LOGICAL, INTENT( IN ):: DAYFLAG        ! true: day-, false: hour-spec
        INTEGER,INTENT(INOUT):: SDATE          ! Julian starting date in TZONE
        INTEGER,INTENT(INOUT):: STIME          ! start time of data in TZONE
        INTEGER, INTENT(OUT) :: EDATE          ! Julian ending date in TZONE
        INTEGER, INTENT(OUT) :: ETIME          ! ending time of data in TZONE
        LOGICAL, INTENT(OUT) :: EASTAT( NIPPA ) ! true: pol/act appears in data
        LOGICAL, INTENT(OUT) :: SPSTAT( MXSPDAT ) ! true: special in data

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
        
C...........   Temporary read arrays
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line

C...........   Local arrays
        REAL              , ALLOCATABLE, SAVE :: DTACBRN( : )    ! storing acre burned value (acre/day) for computing HFLUX
        REAL              , ALLOCATABLE, SAVE :: DTFUELD( : )    ! storing fuel loading value (tons/acre) for computing HFLUX
        REAL              , ALLOCATABLE, SAVE :: DTVAR1 ( : )    ! storing variable 1 value for formula
        REAL              , ALLOCATABLE, SAVE :: DTVAR2 ( : )    ! storing variable 2 value for formula

        INTEGER           , ALLOCATABLE, SAVE :: NSRCPDDAT( :,: )    ! counting number of sources per day/pollutant
        INTEGER           , ALLOCATABLE, SAVE :: IDXSD    ( : )      ! sorting index for CSRCDAYA
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: CSRCDAYA ( : )      ! unsorted source/day array
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: CSRCDAY  ( : )      ! sorted source/day array

C...........   Other local variables
        INTEGER          H, HS, I, II, J, K, L, LL, N, S, T    ! counters and indices
        INTEGER          L0, L1, L2, L3, L4, L5
        INTEGER          ES, NS, SS       ! end src, tmp no. src, start sourc

        INTEGER          D, SD, N1, N2, N3, N4

        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER          FIP              ! tmp co/st/cy code
        INTEGER       :: ICC = 0          ! tmp country code from header
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          JD               ! Julian day number 1...365,366
        INTEGER          JTIME            ! tmp HHMMSS time
        INTEGER          ESTIME           ! tmp HHMMSS episode start time
        INTEGER          EETIME           ! tmp HHMMSS episode end time
        INTEGER          LFIP             ! previous st/co FIPS code
        INTEGER, SAVE :: LOOPNO = 0       ! no. of loops
        INTEGER, SAVE :: MAXPTR           ! maximum time step reference pointer
        INTEGER, SAVE :: MINPTR           ! minimum time step reference pointer
        INTEGER          MONTH            ! tmp month number
        INTEGER, SAVE :: MXWARN           !  maximum number of warnings
        INTEGER, SAVE :: NWARN = 0        !  current number of warnings
        INTEGER, SAVE :: NBADSRC = 0      ! no. bad sources
        INTEGER, SAVE :: NCOMP = 0        ! no. of formulas
        INTEGER, SAVE :: NFIELD = 0       ! number of data fields
        INTEGER, SAVE :: NHFLX   = 0      ! no. of heat flux var 2
        INTEGER, SAVE :: NACRBND = 0      ! no. of acres burned var
        INTEGER, SAVE :: NFUELD  = 0      ! no. of fuel loading var
        INTEGER, SAVE :: NPRCHFX = 0      ! no. of precomputed HFLUX
        INTEGER, SAVE :: NPRCFRM = 0      ! no. of precomputed formula results
        INTEGER, SAVE :: NFUEL  = 0       ! no. of FUEL_LOAD
        INTEGER, SAVE :: NFRM   = 0       ! no. of records with precomputed formula value
        INTEGER, SAVE :: NVAR1  = 0       ! no. of first variables in formula
        INTEGER, SAVE :: NVAR2  = 0       ! no. of second variables in formula
        INTEGER       :: NPOA   = 0       ! unused header number of pol/act
        INTEGER, SAVE :: NSRCDAY = 0      ! no. of src/day combos for computed vars
        INTEGER, SAVE :: NSTEPS = 0       ! number of time steps
        INTEGER          PTR              ! tmp time step pointer
        INTEGER       :: RDATE = 1980001  ! reference date: Jan 1, 1980
        INTEGER       :: RTIME = 0        ! reference time
        INTEGER, SAVE :: SDATESAV = 0     ! saved start date
        INTEGER, SAVE :: STIMESAV = 0     ! saved start time
        INTEGER, SAVE :: TDIVIDE  = 1     ! time step divisor
        INTEGER          WD               ! tmp field width
        INTEGER          YEAR             ! 4-digit year
        INTEGER       :: YR4 = 0          ! unused header year
        INTEGER          ZONE             ! source time zones

        REAL             TDAT             ! temporary data values

        LOGICAL, SAVE :: CHKMINUS 
        LOGICAL, SAVE :: IFLAG  = .FALSE. ! true: Open annual/average inventory
        LOGICAL, SAVE :: DFLAG  = .FALSE. ! true: dates set by data
        LOGICAL       :: EFLAG  = .FALSE. ! TRUE iff ERROR
        LOGICAL       :: WARNOUT= .FALSE. ! true: then output warnings
        LOGICAL, SAVE :: PRCHFX = .FALSE. ! true: skip adding HFLUX due to precomputed heat flux
        LOGICAL, SAVE :: PRCFRM = .FALSE. ! true: skip computing formula due to precomputed values
        LOGICAL       :: HFXFLAG= .FALSE. ! true: adding HFLUX into a list
        LOGICAL       :: FRMFLAG= .FALSE. ! true: adding formula values into a list
        LOGICAL       :: BNHRFLAG=.FALSE. ! true: adding BEGHOUR into a list
        LOGICAL       :: ENHRFLAG=.FALSE. ! true: adding ENDHOUR into a list
        LOGICAL, SAVE :: FIRSTCOUNT = .TRUE.! true: until after first time routine is called with GETCOUNT=TRUE
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called

        CHARACTER(256) :: BUFFER = ' '    ! src description buffer 
        CHARACTER(300) :: LINE   = ' '    ! line buffer 
        CHARACTER(300) :: MESG   = ' '    ! message buffer

C.........  Saved local character variables
        CHARACTER(IOVLEN3), SAVE :: FVAR      ! name of formula resulting variable
        CHARACTER(IOVLEN3), SAVE :: INAMVAR1  ! tmp variable1 cas number (*16)
        CHARACTER(IOVLEN3), SAVE :: INAMVAR2  ! tmp variable2 cas number (*16)
        CHARACTER(IOVLEN3), SAVE :: VAR1      ! formula input variable 1
        CHARACTER(IOVLEN3), SAVE :: VAR2      ! formula input variable 2

C.........  Temporary local character variables
        CHARACTER(FIPLEN3) CFIP      ! tmp co/st/cy code
        CHARACTER(IOVLEN3) CDAT      ! tmp data name (*16)
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
        CHARACTER(512)     VAR_FORMULA !  formula string

        CHARACTER(16) :: PROGNAME = 'RDORLFR' !  program name

C***********************************************************************
C   begin body of program RDORLFR

C.........  First time routine called
        IF( FIRSTIME ) THEN

C.............  Get value of these controls from the environment
            IFLAG = ENVYN ( 'IMPORT_AVEINV_YN', ' ', .TRUE., IOS )

            MESG = 'Inventory formula'
            CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )

C.............  Figure out how many variables there are based on the
C               number of commas found in the string.
            IF( LEN( VAR_FORMULA ) > 0 ) NCOMP = 1
            DO I = 1, L
                IF( VAR_FORMULA( I:I ) == ',' ) NCOMP = NCOMP + 1
            ENDDO

C.............  Only one formula is allowed
            IF ( NCOMP > 1 ) THEN
                WRITE( MESG,94010 ) 'Only one formula (not the',NCOMP,
     &                 'given) are allowed when running day-'//
     &                 'specific fires.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Obtain names for formula variables and whether formula is minus or plus
            L1 = INDEX( VAR_FORMULA, '=' )
            L2 = INDEX( VAR_FORMULA, '+' )
            L3 = INDEX( VAR_FORMULA, '-' )

            CHKMINUS = .FALSE.
            CHKMINUS = ( L3 .GT. 0 )

            IF( L1 .LE. 0 .OR. ( L2 .LE. 0 .AND. L3 .LE. 0 ) ) THEN
                MESG = 'Could not interpret formula for extra ' //
     &                 'pollutant from environment variable ' //
     &                 CRLF() // BLANK10 // '"' // TRIM( FORMEVNM ) //
     &                 '": ' // TRIM( VAR_FORMULA )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.................  Extract formula variable names
            L4 = L2
            IF( CHKMINUS ) L4 = L3

            L      = LEN_TRIM( VAR_FORMULA )
            FVAR = ADJUSTL( VAR_FORMULA(    1:L1-1 ) )
            VAR1 = ADJUSTL( VAR_FORMULA( L1+1:L4-1 ) )
            VAR2 = ADJUSTL( VAR_FORMULA( L4+1:L    ) )

            ALLOCATE( NSRCPDDAT( 366,NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NSRCPDDAT', PROGNAME )
            NSRCPDDAT = 0  ! array

C.............  Check Pollutant Code Numbers for formula pollutants
            I = INDEX1( VAR1, NINVTBL, ITNAMA )
            INAMVAR1 = ITCASA( I )
            IF( INAMVAR1 == ' ' ) THEN
                MESG = 'ERROR: Enter Inventory Polluant Code for '//
     &                 TRIM( VAR1 ) // ' pollutant ' //
     &                 'in the inventory table file ($INVTABLE).'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = INDEX1( VAR2, NINVTBL, ITNAMA )
            INAMVAR2 = ITCASA( I )
            IF( INAMVAR2 == ' ' ) THEN
                MESG = 'ERROR: Enter Inventory Polluant Code for '//
     &                 TRIM( VAR2 ) // ' pollutant ' //
     &                 'in the inventory table file ($INVTABLE).'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

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
                    IF( IFIP( S ) .EQ. INVIFIP( I ) ) THEN
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

        END IF  ! End first time subroutine is called

C.........  For the first call in a loop of files, initialize variables
        IF( FIRSTCALL ) THEN
            MINPTR  = 99999999
            MAXPTR  = 0

C.............  Set time step divisor
            TDIVIDE = 3600 * TSTEP / 10000

C.............  Set the number of day-specific fields
            IF( DAYFLAG ) NFIELD = 1

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

        END IF

C.........  For the second pass of this routine, allocate memory needed for calculating
C           values from day-specific data on the fly.
        IF ( GETCOUNT .AND. FIRSTCOUNT ) THEN

C.............  Determine how much memory is needed for allocating arrays
C               for computing the formula.  This should be the maximum number of 
C               source/days of formula variables, fuel load, and acres burned.
            DO D = 1, 366
                DO I = 1, NIPPA
                    IF( EANAM(I)== VAR1 )         N1= N1+ NSRCPDDAT(D,I)
                    IF( EANAM(I)== VAR2 )         N2= N2+ NSRCPDDAT(D,I)
                    IF( EANAM(I)=='FUEL_LOAD' )   N3= N3+ NSRCPDDAT(D,I)
                    IF( EANAM(I)=='ACRESBURNED' ) N4= N4+ NSRCPDDAT(D,I)
                END DO
            END DO
            N = MAX( N1, N2, N3, N4 )

            ALLOCATE( IDXSD( N ), STAT=IOS )          ! Sorting index
            CALL CHECKMEM( IOS, 'IDXSD', PROGNAME )
            ALLOCATE( CSRCDAYA( N ), STAT=IOS )       ! Unsorted SOURCE/DAY combos
            CALL CHECKMEM( IOS, 'CSRCDAYA', PROGNAME )
            ALLOCATE( CSRCDAY( N ), STAT=IOS )        ! Sorted SOURCE/DAY combos
            CALL CHECKMEM( IOS, 'CSRCDAY', PROGNAME )
            ALLOCATE( DTVAR1( N ), STAT=IOS )         ! To store 1st formula variable values
            CALL CHECKMEM( IOS, 'DTVAR1', PROGNAME )
            ALLOCATE( DTVAR2( N ), STAT=IOS )         ! To store 2nd formula variable values
            CALL CHECKMEM( IOS, 'DTVAR2', PROGNAME )
            ALLOCATE( DTACBRN( N ), STAT=IOS )  ! To store acres burned
            CALL CHECKMEM( IOS, 'DTACBRN', PROGNAME )
            ALLOCATE( DTFUELD( N ), STAT=IOS )  ! To store fuel load
            CALL CHECKMEM( IOS, 'DTFUELD', PROGNAME )

            IDXSD    = 0
            CSRCDAYA = ' '
            CSRCDAY  = ' '
            DTVAR1   = BADVAL3
            DTVAR2   = BADVAL3
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
            CALL PADZERO( SEGMENT( 1 )( 1:5 ) )
            WRITE( CFIP, '(I1,A)' ) ICC, SEGMENT( 1 )( 1:5 )  ! country code of FIPS     
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

C.............  Counting the number of times precomputed HFLUX and precomputed formula
C               values appear in the input file
            IF( CDAT == 'HFLUX' .AND. .NOT. HFXFLAG ) THEN
                NPRCHFX = NPRCHFX + 1
            END IF

            IF( CDAT == FVAR .AND. .NOT. FRMFLAG ) THEN
                NPRCFRM = NPRCFRM + 1
            END IF

C.............  Counting and adding HFLUX, BEGHOUR, ENDHOUR and/or formula variables
C               building a list of source characteristics and store
            IF( HFXFLAG .AND. .NOT. PRCHFX ) CDAT = 'HFLUX'
            IF( BNHRFLAG ) CDAT = 'BEGHOUR'
            IF( ENHRFLAG ) CDAT = 'ENDHOUR'
            IF( FRMFLAG .AND. .NOT. PRCFRM ) CDAT = FVAR

C.............  Adding additional variables and lines if necessary
            HFXFLAG = .FALSE.
            IF( CDAT == 'ACRESBURNED' .AND. .NOT. PRCHFX ) THEN
                HFXFLAG = .TRUE.    ! indicating adding HFLUX
                BACKSPACE( FDEV )
                NACRBND = NACRBND + 1
            END IF

            IF( CDAT == 'FUEL_LOAD' ) THEN
                NFUELD = NFUELD + 1
            END IF

            FRMFLAG = .FALSE.
            IF( CDAT == INAMVAR1 .AND. .NOT. PRCFRM ) THEN
                FRMFLAG = .TRUE.   ! indicating adding formula var
                BACKSPACE( FDEV )
                NVAR1 = NVAR1 + 1
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

C.............  Set Julian day from MMDDYY8 SAS format
            MONTH = STR2INT( DATE( 1:2 ) )
            DAY   = STR2INT( DATE( 4:5 ) )
            YEAR  = YEAR4( STR2INT( DATE( 7:8 ) ) )

            JD = JULIAN( YEAR, MONTH, DAY )
            JDATE = 1000 * YEAR + JD
            JTIME = 0
            
C.............  Set time zone number
            ZONE = GETTZONE( FIP )
            
C.............  If daily emissions are not in the output time zone, print 
C               warning
            IF( WARNOUT .AND. DAYFLAG .AND. ZONE .NE. TZONE 
     &          .AND. GETCOUNT ) THEN
                WRITE( MESG,94010 ) 
     &                'WARNING: Time zone ', ZONE, 'in day-specific ' //
     &                'file at line of pollutant ' // TRIM( CDAT ) //
     &                ' on ' // TRIM( DATE ) // 
     &                ' does not match output time zone', TZONE
                CALL M3MESG( MESG )

            END IF

C.............  Convert date and time to output time zone.
            CALL NEXTIME( JDATE, JTIME, ( ZONE - TZONE ) * 10000 )

C.............  Determine time step pointer based on reference time
            PTR = SECSDIFF( RDATE, RTIME, JDATE, JTIME ) / TDIVIDE + 1

C.............  Store minimum time step number as compared to reference
            IF( PTR .LT. MINPTR ) MINPTR = PTR

C.............  Store maximum time step number as compared to rference
            IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

C.............  Check pollutant code and set index I
            COD  = INDEX1( CDAT, NIPPA, EANAM )
            IF( COD > 0  ) II = 1

            IF( COD .LE. 0 ) THEN
                II = INDEX1( CDAT, NINVTBL, ITCASA )

C.................  This code resets the pollutant name CDAT to the SMOKE name
C                   instead of the inventory name (which could be the CAS number)
                IF( II > 1 ) THEN
                    CTMP = CDAT
                    CDAT = ADJUSTL( ITNAMA( II ) )
                    COD  = INDEX1( CDAT, NIPPA, EANAM )
                END IF
            END IF

C.............  Check to see if data name is in Inventory Table list
            IF( II <= 0 ) THEN

                J = 0
                IF( NBADPOLS > 0 ) J = INDEX1( CDAT, NBADPOLS, BADPOLS )

                IF ( J <= 0 ) THEN

                    WRITE( MESG,93000 ) 
     &                 'WARNING: Skipping pollutant "'// TRIM( CDAT ) //
     &                 '" in PTDAY file because it is not in the '//
     &                 'Inventory Table (INVTABLE) file.'
                    CALL M3MSG2( MESG )

                    NBADPOLS = NBADPOLS + 1
                    BADPOLS( NBADPOLS ) = CDAT

                END IF

                CYCLE      !  to head of read loop

C.............  Check to see if data name is in inventory data list
            ELSE IF( COD .LE. 0 ) THEN

                IF ( .NOT. LCODWARN( II ) ) CYCLE  ! If have already warned, then cycle

                LCODWARN( II ) = .TRUE.

                WRITE( MESG,93000 ) 
     &                'WARNING: Skipping pollutant "'// TRIM( CTMP ) //
     &                '" in PTDAY file because it is not in PTINV file.'
                CALL M3MSG2( MESG )

                CYCLE      !  to head of read loop

C.............  If it is, store status of inventory data
            ELSE 
                EASTAT( COD ) = .TRUE.

            END IF

C.............  Count the number of sources per day & pollutant/variable
C.............  This will give us how many source/date combos there are for 
C               any variables, including HFLUX, and formula variables
            NSRCPDDAT( JD, COD ) = NSRCPDDAT( JD, COD ) + 1
            
C.............  If only getting dates and pollutant information, go 
C               to next loop iteration
            IF( GETSIZES ) CYCLE
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

C.............  Determine time step pointer based on actual start time
            PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME )/TDIVIDE + 1

C.............  Skip record if it is out of range of output file
C.............  NOTE - this is only useful if reading only part of data
            IF( PTR. LT. 1 .OR. PTR .GT. NSTEPS ) CYCLE
            
C.............  Count estimated record count per time step
            DO T = PTR, MIN( PTR + 23, NSTEPS )
                MXPDPT( T ) = MXPDPT( T ) + 1
            END DO

C.............  Store formula variable values.  Only need to do this on the the second
C               pass.  Need to do this before the third pass through the data, because
C               that is when the formula calculation is made.           
            IF( GETCOUNT .AND.                                  ! Second pass through code ONLY
     &        ( .NOT. PRCFRM .AND.                              ! No precomputed formula AND      
     &           ( FRMFLAG .OR. CDAT == VAR2 ) ) .OR.           ! either formula values
     &        ( .NOT. PRCHFX .AND.                              ! No precomputed Heat flux
     &           ( HFXFLAG .OR. CDAT == 'FUEL_LOAD' ) )  ) THEN ! Acres burned value or fuel load value

C.................  Figure out which source/day this is for storing in correct source/day
C.................  This code does *not* assume that the data have been sorted first.
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                        TSCC, DATE, POLBLNK3, CSRCD )

C.................  Build unsorted arrays of source/days and emissions for calculating formula
                SD = 0
                IF( NSRCDAY > 0 ) SD= INDEX1( CSRCD, NSRCDAY, CSRCDAYA )
                IF( SD <= 0 ) THEN
                    NSRCDAY = NSRCDAY + 1
                    SD = NSRCDAY
                    CSRCDAYA( SD ) = CSRCD
                    IDXSD   ( SD ) = SD
                END IF

                IF( FRMFLAG ) DTVAR1( SD ) = TDAT         ! storing variable 1 for formula calc
                IF( CDAT == VAR2 ) DTVAR2( SD ) = TDAT    ! storing variable 2 for formula calc
                IF( HFXFLAG ) DTACBRN( SD ) = TDAT        ! storing acres burned
                IF( CDAT == 'FUEL_LOAD' ) DTFUELD( SD ) = TDAT ! storing fuel load

            END IF
            
C.............  If only counting records per time step, go to next loop
C               iteration
            IF( GETCOUNT ) CYCLE
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

C.............  If FIPS code is not the same as last time, then
C               look it up and get indidies
            IF( FIP .NE. LFIP ) THEN
                J = FIND1( FIP, NINVIFIP, INVIFIP )
                IF( J .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'INTERNAL ERROR: Could not ',
     &                     'find FIPS code', FIP, 'in internal list.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                END IF

                SS = STARTSRC( J )
                ES = ENDSRC( J )
                NS = ES - SS + 1
                LFIP = FIP

            END IF

C.............  Include hack to change how CSRC is constructed
C               when CSOURCE has been read from inventory PNTS file
C               instead of created in the same SMOKE run.
            IF( .NOT. IFLAG ) THEN
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID//'     ', 
     &                        TSCC, CHRBLNK3, POLBLNK3, CSRC )
            ELSE
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                        TSCC, CHRBLNK3, POLBLNK3, CSRC )
            END IF

C.............  Search for this record in sources
            J = FINDC( CSRC, NS, CSOURC( SS ) )

C.............  If source not found, store source in list of bad sources
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
                    MESG = 'WARNING: Period-specific record does ' //
     &                     'not match inventory sources: '//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                END IF

                CYCLE               !  to head of read loop

C.............  Otherwise, update master list of sources in the inventory
            ELSE
                S = SS - 1 + J         ! calculate source number
                LPDSRC( S ) = .TRUE.

            END IF

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

C.............  Compute formula.  This calculation uses stored values of formula inputs that
C               were stored on the second pass through this subroutine.
            IF( CDAT == FVAR .AND. NVAR1 > 0 .AND. .NOT. PRCFRM ) THEN

C.................  Build source/date string to lookup position for doing calculation
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, TSCC, DATE,
     &                        POLBLNK3, CSRCD )

C.................  Lookup source/date string in master list to get position
                SD = FINDC( CSRCD, NSRCDAY, CSRCDAY )
                K  = IDXSD( SD )

C.....................  If PM2.5 value is missing, then assume zero
                IF( DTVAR2( K ) < AMISS3 ) THEN
                    LL = LEN_TRIM( CSRCD )
                    CALL FMTCSRC( CSRCD, 6, BUFFER, L2 )

                    MESG = 'WARNING: Resetting missing value of '//
     &                     'PM2.5 to 0. for source:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) // 
     &                     ' on date ' // CSRCD( LL-7: LL )
                    CALL M3MSG2( MESG )

                    DTVAR2( K ) = 0.0

                END IF

C.....................  Compute formula value
               IF( CHKMINUS ) THEN
                   TDAT = DTVAR1( K ) - DTVAR2( K )  ! computing formula result
               ELSE
                   TDAT = DTVAR1( K ) + DTVAR2( K )  ! computing formula result
               END IF

               IF( TDAT < 0 ) THEN
                    LL = LEN_TRIM( CSRCD )
                    CALL FMTCSRC( CSRCD, 6, BUFFER, L2 )

                    MESG = 'WARNING: Resetting negative value of '//
     &                     'computed variable ' // TRIM( FVAR )// 
     &                     ' to 0. for source:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) // 
     &                     ' on date ' // CSRCD( LL-7: LL )
                    CALL M3MSG2( MESG )

                    TDAT = 0.0

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
                    EMISVA( HS,T ) = TDAT  ! Store data in emissions
                    DYTOTA( HS,T ) = TDAT

                END IF
            END DO

        END DO     ! Main read loop of day-specific data

299     CONTINUE   ! Exit from read loop

C.........  Warning messages for HFLUX and formula result (QA checks)
        IF( GETCOUNT ) THEN
            IF( NPRCHFX > 0 .AND. .NOT. PRCHFX ) THEN
                MESG = 'WARNING: Skipping internal heat flux '//
     &                 'computation due to the existence of '//
     &                 'precomputed HFLUX in PTDAY'
                CALL M3MSG2( MESG )
                PRCHFX = .TRUE.
            END IF

            IF( NPRCFRM > 0 .AND. .NOT. PRCFRM ) THEN
                MESG = 'WARNING: Skipping formula computation for '//
     &                 TRIM(FVAR)// ' due to the existence of '//
     &                 'precomputed '//TRIM(FVAR)// ' in PTDAY'
                CALL M3MSG2( MESG )
                PRCFRM = .TRUE.
            END IF

C.............  Give warning if variables needed for formula are not present  
            IF( NVAR1 < 1 .AND. .NOT. PRCFRM ) THEN
                MESG = 'WARNING: No ' // TRIM( VAR1 ) //
     &                 ' data are available to '//
     &                 'compute '//TRIM(FVAR)// '- values will be 0.0'
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
