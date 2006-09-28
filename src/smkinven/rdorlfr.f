
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
        USE MODSOURC, ONLY: IFIP, CSOURC, HEATCONTENT, FIREPOL, NFRPOL

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP, UCASNKEP, NUNIQCAS,
     &                      UNIQCAS, NINVTBL, ITNAMA, ITCASA, FIREFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NSRC, EANAM, NCHARS, NMAP, MAPNAM,
     &                     MAPFIL

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, IDXSRC, SPDIDA,
     &                      CODEA, EMISVA, DYTOTA

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

C...........   Local list of bad sources to prevent duplicate writing of error
C              messages
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )

C...........   Local list of FIPS start/end positions to facilitate
C              faster lookups
        INTEGER, ALLOCATABLE, SAVE :: STARTSRC( : )
        INTEGER, ALLOCATABLE, SAVE :: ENDSRC( : )
        
C...........   Temporary read arrays
        REAL            TDAT( 24 )       ! temporary data values

C...........   Local arrays
        INTEGER                               :: IDXFIRE( NSEG ) ! index for wildfire pollutants
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: HFXBSRC( : )    ! Build source characteristics field for HFLUX
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: PMCBSRC( : )    ! Build source characteristics field for PMC
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: FULBSRC( : )    ! Build source characteristics field for FUEL_LOAD
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: P25BSRC( : )    ! Build source characteristics field for PM2_5
        REAL              , ALLOCATABLE, SAVE :: DTACBRN( : )    ! storing acre burned value (acre/day) for computing HFLUX
        REAL              , ALLOCATABLE, SAVE :: DTFUELD( : )    ! storing fuel loading value (tons/acre) for computing HFLUX
        REAL              , ALLOCATABLE, SAVE :: DTPM10 ( : )    ! storing PM10 value for computing PMC
        REAL              , ALLOCATABLE, SAVE :: DTPM25 ( : )    ! storing PM2_5 value for computing PMC

C...........   Other local variables
        INTEGER          H, HS, I, II, ID, IC, J, K, L, LL, N, S, T    ! counters and indices
        INTEGER          L0, L1, L2, L3, L4, L5
        INTEGER          ES, NS, SS       ! end src, tmp no. src, start sourc

        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER          FIP              ! tmp co/st/cy code
        INTEGER       :: ICC = 0          ! tmp country code from header
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          EPSDATE           ! tmp episode julian date
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
        INTEGER, SAVE :: NFIELD = 0       ! number of data fields
        INTEGER, SAVE :: NHFLX   = 0      ! no. of heat flux var 2
        INTEGER, SAVE :: NACRBND = 0      ! no. of acres burned var
        INTEGER, SAVE :: NFUELD  = 0      ! no. of fuel loading var
        INTEGER, SAVE :: NPRCHFX = 0      ! no. of precomputed HFLUX
        INTEGER, SAVE :: NPRCPMC = 0      ! no. of precomputed PMC
        INTEGER, SAVE :: NFUEL   = 0      ! no. of FUEL_LOAD
        INTEGER, SAVE :: NPMC    = 0      ! no. of coarse particulate matter var
        INTEGER, SAVE :: NPM25   = 0      ! no. of PM2_5
        INTEGER, SAVE :: NPM10  = 0       ! no. of PM10 var
        INTEGER       :: NPOA   = 0       ! unused header number of pol/act
        INTEGER       :: NSEG   = 32      ! maximum no of segments
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

        REAL             TOTAL            ! tmp daily total of hourly file
        REAL          :: CONVT = 2000     ! conversion factor (ton to pounds)

        LOGICAL, SAVE :: IFLAG  = .FALSE. ! true: Open annual/average inventory
        LOGICAL, SAVE :: DFLAG  = .FALSE. ! true: dates set by data
        LOGICAL       :: EFLAG  = .FALSE. ! TRUE iff ERROR
        LOGICAL       :: WARNOUT= .FALSE. ! true: then output warnings
        LOGICAL, SAVE :: PRCHFX = .FALSE. ! true: skip adding HFLUX due to precomputed heat flux
        LOGICAL, SAVE :: PRCPMC = .FALSE. ! true: skip adding PMC due to precomputed PMC
        LOGICAL       :: HFXFLAG= .FALSE. ! true: adding HFLUX into a list
        LOGICAL       :: PMCFLAG= .FALSE. ! true: adding PMC into a list
        LOGICAL       :: BNHRFLAG=.FALSE. ! true: adding BEGHOUR into a list
        LOGICAL       :: ENHRFLAG=.FALSE. ! true: adding ENDHOUR into a list
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called
        LOGICAL, SAVE :: TFLAG  = .FALSE. ! true: use SCCs for matching with inv

        CHARACTER(256) :: BUFFER = ' '    ! src description buffer 
        CHARACTER(300) :: LINE   = ' '    ! line buffer 
        CHARACTER(300) :: MESG   = ' '    ! message buffer

        CHARACTER(FIPLEN3) CFIP      ! tmp co/st/cy code
        CHARACTER(IOVLEN3) CDAT      ! tmp data name (*16)
        CHARACTER(IOVLEN3) PCAS10    ! tmp PM10 cas number (*16)
        CHARACTER(IOVLEN3) PCAS25    ! tmp PM25 cas number (*16)
        CHARACTER(IOVLEN3) CTMP      ! tmp data name (*16)
        CHARACTER(CHRLEN3) CHAR4     ! tmp characteristic 4  (*15)
        CHARACTER(PLTLEN3) FCID      ! tmp facility ID (*15)
        CHARACTER(CHRLEN3) SKID      ! tmp stack ID (*15) = LocID
        CHARACTER(CHRLEN3) DVID      ! dummy device ID
        CHARACTER(CHRLEN3) PRID      ! dummy process ID
        CHARACTER(SCCLEN3) TSCC      ! tmp source category code (*10)
        CHARACTER(ALLLEN3) CSRC      ! tmp source string
        CHARACTER(ALLLEN3) TSRC      ! tmp source string
        CHARACTER( 8 )     DATE      ! tmp date string
        CHARACTER(IOVLEN3) PNAME     ! logical file name for data files
        CHARACTER(IOVLEN3) INVAR     ! tmp inventory pollutant name
        CHARACTER(40)      SEGMENT( 10 ) ! segments of line

        CHARACTER(16) :: PROGNAME = 'RDORLFR' !  program name

C***********************************************************************
C   begin body of program RDORLFR

C.........  First time routine called
        IF( FIRSTIME ) THEN

C.........  Get value of these controls from the environment
            IFLAG = ENVYN ( 'IMPORT_AVEINV_YN', ' ', .TRUE., IOS )

C.........  allocate memories for building source characteristics for 
C           heat flux and coarse PMC (if applicable) and for storing 
C           acre burned, and fuel loading for computing HFLUX and PM10 
C           PM2_5 for computing PMC

            ALLOCATE( HFXBSRC( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HFXBSRC', PROGNAME )
            ALLOCATE( FULBSRC( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FULBSRC', PROGNAME )
            ALLOCATE( DTACBRN( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DTACBRN', PROGNAME )
            ALLOCATE( DTFUELD( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DTFUELD', PROGNAME )

            HFXBSRC = ' '
            FULBSRC = ' '
            DTACBRN = 0.
            DTFUELD = 0.

            ALLOCATE( PMCBSRC( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PMCBSRC', PROGNAME )
            ALLOCATE( P25BSRC( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'P25BSRC', PROGNAME )
            ALLOCATE( DTPM10( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DTPM10', PROGNAME )
            ALLOCATE( DTPM25( NSRC*10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DTPM25', PROGNAME )

            PMCBSRC = ' '
            P25BSRC = ' '
            DTPM10 = 0.
            DTPM25 = 0.

C.............  Check CAS numbers for pollutants, fuel_load, acresburned
            IC = INDEX1( 'PM2_5', NINVTBL, ITNAMA )
            PCAS25 = ITCASA( IC )
            IF( PCAS25 == ' ' ) THEN
                MESG = 'ERROR: Enter CAS nubmer for PM2_5 pollutant ' //
     &                 'in an inventory table ($INVTABLE).'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IC = INDEX1( 'PM10', NINVTBL, ITNAMA )
            PCAS10 = ITCASA( IC )
            IF( PCAS10 == ' ' ) THEN
                MESG = 'ERROR: Enter CAS nubmer for PM10 pollutant ' //
     &                 'in an inventory table ($INVTABLE).'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Allocate memory for bad source storage
            ALLOCATE( BADSRC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

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

C............  Loop through all pollutant files...
            DO N = 1, NMAP
            
                INVAR = MAPNAM( N )
                IF( INVAR .NE. 'HEATCONTENT' ) CYCLE

C............  Set environment variable for input file
                PNAME = 'TMP_POL_FILE'
                IF( .NOT. SETENVVAR( PNAME, MAPFIL( N ) ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not set logical file name ' //
     &                     CRLF() // BLANK10 // 'for file ' //
     &                     TRIM( MAPFIL( N ) )
                    CALL M3MSG2( MESG )

C................  Open pollutant file with physical file name
                ELSE IF( .NOT. OPENSET( PNAME, FSREAD3, PROGNAME ) )THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open "'// TRIM( MAPNAM(N) )
     &                   // '" file listed in ' // CRLF() // BLANK10 // 
     &                   'map-formatted inventory for file name: ' //
     &                   CRLF() // BLANK10 // TRIM( MAPFIL( N ) )
                    CALL M3MSG2( MESG )

                END IF

                SELECT CASE( INVAR )

C............  Open I/O API inventory HEATCONTENT file and store
                CASE( 'HEATCONTENT' )

                    ALLOCATE( HEATCONTENT( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'HEATCONENT', PROGNAME )

                    IF( .NOT. READSET( PNAME, 'HEATCONTENT', ALLAYS3, 
     &                  1, 0, 0, HEATCONTENT ) ) THEN
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

                END SELECT

                IF( .NOT. CLOSESET( PNAME ) ) THEN
                    MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( MAPFIL( N ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END DO

            FIRSTIME = .FALSE.
            FIREFLAG = .TRUE.

        END IF

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

C.........  Loop through file and read it. In the first section, determine
C           the minimum and maximum date. Use a reference date to do this. In
C           the second section, determine the number of records per time 
C           step. In the third section, read and store the data.  When storing
C           data, time step index is computed from the start date/time instead
C           of the reference date/time so that the indexing will work properly.
        IREC = 0
        TDAT = 0   !  array

        DO         !  Head of period-specific file read loop

C.............  Read first line of file
            READ( FDEV, 93000, END=299 ) LINE
            IREC = IREC + 1

            L = LEN_TRIM( LINE )
            L0 = INDEX( LINE, '#DATA' )

C.............  Skip blank lines 
            IF( L == 0 ) CYCLE
            IF( L0 > 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set
C               properly
            CALL GETHDR( 1, .FALSE., .FALSE., .FALSE., LINE, ICC, YR4,
     &                   NPOA, IOS )

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

C.............  Get lines
            CALL PARSLINE( LINE, 9, SEGMENT )

C.............  Use the file format definition to parse the line into
C           the various data fields
            WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS     
            CFIP( 2:6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )  ! state/county code

C.............  Replace blanks with zeros        
            DO I = 1,FIPLEN3
                IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
            END DO

            FIP    = STR2INT( CFIP )          ! FIP codes
            FCID   = ADJUSTL( SEGMENT( 2 ) )  ! fire ID
            SKID   = ADJUSTL( SEGMENT( 3 ) )  ! location ID
            TSCC   = ADJUSTL( SEGMENT( 4 ) )  ! SCC
            DVID   = ADJUSTL( ' ' )           ! dummy device id
            PRID   = ADJUSTL( ' ' )           ! dummy process id
            CDAT   = ADJUSTL( SEGMENT( 5 ) )  ! Pollutants(FUEL_LOAD, ACRESBURNED,,,)
            DATE   = ADJUSTL( SEGMENT( 6 ) )  ! Date of episode
            ESTIME = STR2INT( SEGMENT( 8 ) ) * 10000 ! episode start time
            EETIME = STR2INT( SEGMENT( 9 ) ) * 10000 ! episode end time

C............. Check fire beginning and ending time format and print warning if necessary
            IF( EETIME > 230000 ) THEN
                MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //
     &                 ' Date: ' // DATE  // ' :: Fire beginning' //
     &                 ' and ending time not in the acceptable range '//
     &                 'of 0 to 23'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE

            ELSE IF( ESTIME > EETIME ) THEN
                MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //
     &                 ' Date: ' // DATE  // ' :: Fire ending time' //
     &                 ' can not be earlier then a begining time'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE

            END IF

C.............  Check date format
            ID = INDEX( DATE, '/' )
            IF( ID <= 0 ) THEN
                MESG = 'ERROR: Incorrect date format ( MM/DD/YY ) :' //
     &                 ' Region: '// CFIP // ' SCC: ' // TSCC // 
     &                 ' Date: ' // DATE
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

C.............  Check and Set emissions values
            TDAT( 1 )  = STR2REAL( SEGMENT( 7 ) )     ! Day-specific total emission

            IF ( TDAT( 1 ) .LT. 0.0 )  THEN
                EFLAG = .TRUE.
                WRITE( MESG,94030 ) 'ERROR: Bad data value "',
     &              TDAT( 1 ), '" of ' // TRIM( CDAT ) //
     &              ' for region ' // CFIP //
     &              ' and SCC ' // TSCC // ' on date '// DATE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                CYCLE  ! to head of read loop
            END IF

C.............  Check a list of pol names
            IF( CDAT == 'HFLUX' .AND. .NOT. HFXFLAG ) THEN
                NPRCHFX = NPRCHFX + 1
            END IF

            IF( CDAT == 'PMC' .AND. .NOT. PMCFLAG ) THEN
                NPRCPMC = NPRCPMC + 1
            END IF

C.............  Counting and adding HFLUX, BEGHOUR, ENDHOUR and/or PMC
C               building a list of source characteristics and store
            IF( HFXFLAG .AND. .NOT. PRCHFX ) THEN
                CDAT = ADJUSTL( 'HFLUX' )
                IF( GETSIZES ) THEN
                    NHFLX = NHFLX + 1
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                            TSCC, DATE, POLBLNK3, CSRC )
                    HFXBSRC( NHFLX ) = CSRC
                    DTACBRN( NHFLX ) = TDAT( 1 )    ! storing ACRESBURNED for HFLUX	calculation
                END IF
            END IF

            IF( BNHRFLAG ) CDAT = ADJUSTL( 'BEGHOUR' )

            IF( ENHRFLAG ) CDAT = ADJUSTL( 'ENDHOUR' )

            IF( PMCFLAG .AND. .NOT. PRCPMC ) THEN
                CDAT = ADJUSTL( 'PMC' )
                IF( GETSIZES ) THEN
                    NPMC = NPMC + 1
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                            TSCC, DATE, POLBLNK3, CSRC )
                    PMCBSRC( NPMC ) = CSRC
                    DTPM10 ( NPMC ) = TDAT( 1 )  ! storing PM10 for PMC calculation
                END IF
            END IF


            IF( GETCOUNT ) THEN
C..................  Build source characteristics field for searching a matching inventory
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, TSCC, DATE,
     &                        POLBLNK3, CSRC )

C.................  Looking for FUEL_LOAD and PM2_5 for HFLUX and PMC computation, respectively
                IF( CDAT == 'FUEL_LOAD' .AND. .NOT. PRCHFX ) THEN
                    DO K = 1, NHFLX
                        TSRC = HFXBSRC( K )
                        IF( TSRC == CSRC ) THEN
                            NFUEL = NFUEL + 1
                            FULBSRC( NFUEL ) = CSRC   ! storing FUEL_LOAD source characteristics
                            DTFUELD( K ) = TDAT( 1 )  ! storing FUEL_LOAD for HFLUX calculation
                        END IF
                    END DO
                END IF

                IF( CDAT == PCAS25 .AND. .NOT. PRCPMC ) THEN
                    DO K = 1, NPMC
                        TSRC = PMCBSRC( K )
                        IF( TSRC == CSRC ) THEN
                            NPM25 = NPM25 + 1
                            P25BSRC( NPM25 ) = CSRC   ! storing PM2_5 source characteristics
                            DTPM25( K ) = TDAT( 1 )   ! storing PM2_5 for PMC calculation
                        END IF
                    END DO
                END IF

            END IF

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

            PMCFLAG = .FALSE.
            IF( CDAT == PCAS10 .AND. .NOT. PRCPMC ) THEN
                PMCFLAG = .TRUE.   ! indicating adding PMC var
                BACKSPACE( FDEV )
                NPM10 = NPM10 + 1
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

            JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )
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

C.............  Convert start and end time from BEGHOUR and ENDHOUR to output time zone
            EPSDATE = JDATE
            CALL NEXTIME( EPSDATE, ESTIME,( ZONE - TZONE ) * 10000 )
            CALL NEXTIME( EPSDATE, EETIME,( ZONE - TZONE ) * 10000 )

C.............  Determine time step pointer based on reference time
            PTR = SECSDIFF( RDATE, RTIME, JDATE, JTIME ) / TDIVIDE + 1

C.............  Store minimum time step number as compared to reference
            IF( PTR .LT. MINPTR ) MINPTR = PTR

C.............  Store maximum time step number as compared to rference
            IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

C.............  Check pollutant code and set index I
            COD  = INDEX1( CDAT, NIPPA, EANAM )

            IF( COD .LE. 0 ) THEN
                II = INDEX1( CDAT, NINVTBL, ITCASA )
                
                
                IF( II > 1 ) THEN
                    CTMP = CDAT
                    CDAT = ADJUSTL( ITNAMA( II ) )
                    COD  = INDEX1( CDAT, NIPPA, EANAM )
                END IF
            END IF

C.............  Check to see if data name is in inventory list
            IF( COD .LE. 0 ) THEN
                IF( WARNOUT .AND. NWARN < MXWARN ) THEN
                    L = LEN_TRIM( CTMP )
                    WRITE( MESG,93000 ) 
     &                 'WARNING: Skipping pollutant "'// CTMP //
     &                 '" : Region: '// CFIP // ' SCC: ' // TSCC //
     &                 ' Date: ' //TRIM( DATE )//' - not in inventory'
                    CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF
                CYCLE      !  to head of loop

C.............  If it is, store status of inventory data
            ELSE 
                EASTAT( COD ) = .TRUE.

            END IF

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

C.............  If only counting records per time step, go to next loop
C               iteration
            IF( GETCOUNT ) CYCLE

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

C.............  If SCCs are needed for matching...
            IF ( TFLAG ) THEN

                IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )
                CHAR4 = TSCC

                IF( .NOT. IFLAG ) THEN
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID//'     ', 
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                ELSE
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                END IF

C.................  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )

C.............  If SCCs are not being used for matching (at least not yet)...
            ELSE

C.................  Build source characteristics field for searching inventory
                IF( .NOT. IFLAG ) THEN
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID//'     ', 
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                ELSE
                    CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                            TSCC, CHRBLNK3, POLBLNK3, CSRC )
                END IF

C.................  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )

C.................  If source is not found for day-specific processing, see 
C                   if reading the SCC in helps
                IF( J .LE. 0 ) THEN

                    IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )
                    CHAR4 = TSCC

                    IF( .NOT. IFLAG ) THEN
                        CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID//
     &                       '     ', TSCC, CHRBLNK3, POLBLNK3, CSRC )
                    ELSE
                        CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, 
     &                                TSCC, CHRBLNK3, POLBLNK3, CSRC )
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

C.............  Computing HFLUX, BEGHOUR, ENDHOUR (as a default) and PMC (if applicable)
            IF( CDAT == 'HFLUX' .AND. .NOT. PRCHFX ) THEN
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, TSCC, DATE,
     &                        POLBLNK3, CSRC )

                DO K = 1, NHFLX
                    TSRC = HFXBSRC( K )
                    IF( TSRC == CSRC ) THEN
                        TDAT( 1 ) = DTACBRN( K ) * DTFUELD( K )  ! computing HFLUX (BTU/day)
     &                              * HEATCONTENT( S ) * CONVT   ! HEATCONTENT(BTU/lb)=8000
                    END IF
                END DO
            END IF

            IF( CDAT == 'PMC' .AND. NPM10 > 0 .AND. .NOT. PRCPMC ) THEN
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, TSCC, DATE,
     &                        POLBLNK3, CSRC )
                DO K = 1, NPMC
                    TSRC = PMCBSRC( K )
                    IF( TSRC == CSRC ) THEN
                        TDAT( 1 ) = DTPM10( K ) - DTPM25( K )  ! computing PMC
                    END IF
                    IF( TDAT( 1 ) < 0 ) THEN
                        MESG = 'ERROR: Computed PMC = PM10 - PM2_5 ' //
     &                     'is negative '//TSRC( LL-SCCLEN3-8 : LL-8 )//
     &                     ' on date ' // TSRC( LL-7: LL )
                        CALL M3MSG2( MESG )
                        EFLAG = .TRUE.
                    END IF
                END DO
            END IF

            IF( CDAT == 'BEGHOUR' ) TDAT( 1 ) = REAL( ESTIME )  ! storing BEGHOUR
            IF( CDAT == 'ENDHOUR' ) TDAT( 1 ) = REAL( EETIME )  ! storing ENDHOUR

C.............  If daily data, set all TDATs with daily value
            IF( DAYFLAG ) TDAT( 2:24 ) = TDAT( 1 )  ! array

C.............  Record needed data for this source and time step
            H = 0
            TOTAL = 0.0
            DO T = PTR, MIN( PTR + 23, NSTEPS )
                H = H + 1
                NPDPT( T ) = NPDPT( T ) + 1

                HS = NPDPT( T )

                IF( HS .LE. MXPDSRC ) THEN

                    IDXSRC( HS,T ) = HS
                    SPDIDA( HS,T ) = S
                    CODEA ( HS,T ) = COD
                    EMISVA( HS,T ) = TDAT( H )  ! Store data in emissions
                    DYTOTA( HS,T ) = TOTAL

                END IF
            END DO

        END DO

299     CONTINUE   ! Exit from read loop

C.........  Warning messages for HFLUX and PMC (QA checks)
        IF( GETCOUNT ) THEN
            IF( NPRCHFX > 0 .AND. .NOT. PRCHFX ) THEN
                MESG = 'WARNING: Skipping internal heat flux '//
     &                 'computation due to the existence of '//
     &                 'precomputed HFLUX in PTDAY'
                CALL M3MSG2( MESG )
                PRCHFX = .TRUE.
            END IF

            IF( NPRCPMC > 0 .AND. .NOT. PRCPMC ) THEN
                MESG = 'WARNING: Skipping internal coare PM (PMC)'//
     &                 'computation due to the existence of '//
     &                 'precomputed PMC in PTDAY'
                CALL M3MSG2( MESG )
                PRCPMC = .TRUE.
            END IF
            
            IF( NPM10 < 1 .AND. .NOT. PRCPMC ) THEN
                MESG = 'WARNING: No PM10 data are available to '//
     &                 'compute coarse PM (PMC).'
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

C.........  Warning msg of no matching inventory available in day-specific
C           inventory file for HFLUX and/or PMC (if applicable).
        
        IF( .NOT. GETSIZES .AND. .NOT. GETCOUNT ) THEN
            IF( NFUEL .NE. NHFLX .AND. .NOT. PRCHFX ) THEN
                DO I = 1, NHFLX
                    TSRC = HFXBSRC( I )
                    LL = LEN_TRIM( TSRC )
                    K = INDEX1( TSRC, NFUEL, FULBSRC )

                    IF( K < 1 ) THEN
                        MESG = 'ERROR: Missing fuel loading data ' //
     &                    'for Heat Flux computation of SCC ' //
     &                    TSRC( LL-SCCLEN3-8 : LL-8 ) // ' on date ' //
     &                    TSRC( LL-7: LL )
                        CALL M3MSG2( MESG )
                        EFLAG = .TRUE.
                    END IF
                END DO
            END IF

            IF( NPM25 .NE. NPMC .AND. .NOT. PRCPMC ) THEN
                DO I = 1, NPMC
                    TSRC = PMCBSRC( I )
                    LL = LEN_TRIM( TSRC )
                    K = INDEX1( TSRC, NPM25, P25BSRC )

                    IF( K < 1 ) THEN
                        MESG = 'ERROR: Missing PM2_5 data' //
     &                    ' for PMC computation of SCC ' //
     &                    TSRC( LL-SCCLEN3-8 : LL-8 ) // ' on date ' //
     &                    TSRC( LL-7: LL )
                        CALL M3MSG2( MESG )
                        EFLAG = .TRUE.
                    END IF
                END DO
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
