
        SUBROUTINE RDCEMPD( FDEV, TZONE, TSTEP, MXPDSRC, GETSIZES, 
     &                      GETCOUNT, FIRSTCALL, DAYFLAG, SDATE, STIME, 
     &                      EDATE, ETIME, EASTAT, SPSTAT )

C***************************************************************************
C  subroutine body starts at line 189
C
C  DESCRIPTION:
C      This subroutine reads the hourly CEM data, matches it with 
C      the inventory data, and gets set up to output PHOUR and ASCII formats
C      with hourly data.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 10/01 by M. Houyoux
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
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: ORISFLAG, NINVORIS, INVORIS,
     &                      NORISBLR, ORISBLR, IORSMTCH,
     &                      INVORFP, OBSRCNT, OBSRCBG, OBSRCNM

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, EANAM, CATEGORY, BYEAR, NSRC, JSTACK

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: NUNFDORS, UNFDORS, MXPDPT, NPDPT, LPDSRC,
     &                      IDXSRC, SPDIDA, CODEA, EMISVA,
     &                      NOBRLIST, OBRLIST, ANNNOX, ANNSO2, ANNGLOAD,
     &                      ANNSLOAD, ANNHEAT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C.........  EXTERNAL FUNCTIONS
        LOGICAL      BLKORCMT
        CHARACTER(2) CRLF
        REAL         ENVREAL
        LOGICAL      ENVYN
        INTEGER      FINDC
        INTEGER      GETTZONE
        INTEGER      INDEX1
        INTEGER      JULIAN
        INTEGER      SECSDIFF
        INTEGER      YEAR4
        INTEGER      STR2INT

        EXTERNAL     BLKORCMT, CRLF, ENVREAL, ENVYN, FINDC, GETTZONE, 
     &               INDEX1, JULIAN, SECSDIFF, YEAR4, STR2INT

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV           ! file unit no.
        INTEGER, INTENT (IN) :: TZONE          ! output time zone
        INTEGER, INTENT (IN) :: TSTEP          ! time step HHMMSS
        INTEGER, INTENT (IN) :: MXPDSRC        ! max. day- or hr-specific source
        LOGICAL, INTENT (IN) :: GETSIZES       ! true: get no. time steps & pols
        LOGICAL, INTENT (IN) :: GETCOUNT       ! true: get max no. srcs per time
        LOGICAL, INTENT (IN) :: FIRSTCALL      ! true: first call of a loop
        LOGICAL, INTENT (IN) :: DAYFLAG        ! true: day-, false: hour-spec
        INTEGER,INTENT(INOUT):: SDATE          ! Julian starting date in TZONE
        INTEGER,INTENT(INOUT):: STIME          ! start time of data in TZONE
        INTEGER,INTENT(INOUT):: EDATE          ! Julian ending date in TZONE
        INTEGER,INTENT(INOUT):: ETIME          ! ending time of data in TZONE
        INTEGER, INTENT(OUT) :: EASTAT( NIPPA ) ! true: pol/act appears in data
        INTEGER, INTENT(OUT) :: SPSTAT( MXSPDAT ) ! true: special in daa

C...........   Local parameters
        INTEGER, PARAMETER :: NCEMPOL  = 2   ! number of pollutants in CEM data
        INTEGER, PARAMETER :: NOXIDX   = 1   ! index to CEMPOLS for NOX
        INTEGER, PARAMETER :: SO2IDX   = 2   ! index to CEMPOLS for SO2

        CHARACTER(IOVLEN3), PARAMETER :: CEMPOLS( NCEMPOL ) = 
     &                                      ( / 'NOX             ',
     &                                          'SO2             ' / )

C...........   Local allocatable
        REAL, ALLOCATABLE, SAVE :: EMIS( :,: ) ! annual inventory emissions
        REAL, ALLOCATABLE       :: STKVAL( : ) ! summed stack flow rates
        
        CHARACTER(PLTLEN3+CHRLEN3), ALLOCATABLE :: STKNAM( : ) ! stack names by source
        
        LOGICAL, ALLOCATABLE, SAVE :: MASLIST( :,: ) ! flag if have data for ORIS/boiler
                                                     ! and time step

C...........   Temporary read arrays
        INTEGER, SAVE :: CEMPIDX( NCEMPOL )   ! Index from CEMPOLS to EASTAT
        REAL             CEMEMIS( NCEMPOL )   ! CEM emissions for line

C.........  File names and unit numbers
        CHARACTER(IOVLEN3) :: ENAME  ! emis i/o api inven logical name
        CHARACTER(IOVLEN3) :: ANAME  ! emis ASCII inven logical name

C...........   Other local variables
        INTEGER          I, J, L, N, S, S1, V, T  ! counters and indices

        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER, SAVE :: FLOWPOS          ! position of hourly flow rate variable
        INTEGER          HH               ! 2-digit hour (0-23) from CEM data
        INTEGER          INVOBPOS         ! position in inventory ORIS/boiler list
        INTEGER          INVORPOS         ! position in invenotry ORIS list
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          JTIME            ! tmp HHMMSS time
        INTEGER, SAVE :: LOOPNO = 0       ! no. of loops
        INTEGER          MASOBPOS         ! position in master ORIS/boiler list
        INTEGER, SAVE :: MAXPTR           ! maximum time step reference pointer
        INTEGER, SAVE :: MINPTR           ! minimum time step reference pointer
        INTEGER          MONTH            ! tmp month number
        INTEGER, SAVE :: MXCALL           ! max times the routine is called
        INTEGER, SAVE :: MXUNFDORS = 0    ! max no. of bad ORIS IDs
        INTEGER          MXSTKS           ! max stacks per source
        INTEGER, SAVE :: NCALL            ! number of times the routine has been called
        INTEGER          NS               ! tmp no. sources per ORIS/boiler
        INTEGER, SAVE :: NSTEPS = 0       ! number of time steps
        INTEGER          NSTK             ! number of stacks per source
        INTEGER          PTR              ! tmp time step pointer
        INTEGER       :: RDATE = 1980001  ! reference date: Jan 1, 1980
        INTEGER       :: RTIME = 0        ! reference time
        INTEGER, SAVE :: SDATESAV = 0     ! saved start date
        INTEGER, SAVE :: STIMESAV = 0     ! saved start time
        INTEGER, SAVE :: TDIVIDE  = 1     ! time step divisor
        INTEGER          YEAR             ! 4-digit year
        INTEGER          YY               ! 2-digit year
        INTEGER          YYMMDD           ! year, month, and day
        INTEGER          DZONE            ! time shift (ZONE-TZONE)
        INTEGER          ZONE             ! source time zones

        REAL             DENOM            ! denominator for assignment weighting
        REAL             GLOAD            ! tmp gross load
        REAL             HTINPUT          ! tmp heat input
        REAL             NOXRATE          ! tmp NOx emission rate
        REAL             OPTIME           ! tmp operating time (unused)
        REAL             SLOAD            ! tmp steam load
        REAL             ANNFAC           ! factor for calculating hourly emissions
        REAL             EMISVAL          ! calculated hourly emissions
        REAL, SAVE    :: FLOWFAC          ! factor for calculating hourly flow rate
        REAL             FLOWVAL          ! calculate hourly flow rate
        REAL             SUMFLOW          ! summed hourly flow rate

        LOGICAL       :: ALLZERO            ! true: set all values to zero
        LOGICAL       :: CEMPOL   = .FALSE. ! true: processing CEM pollutant
        LOGICAL, SAVE :: DFLAG    = .FALSE. ! true: dates have been set by the data
        LOGICAL       :: EFLAG    = .FALSE. ! true: an error has occurred
        LOGICAL       :: MATCHFLG = .FALSE. ! true: a CEM/inventory match found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL       :: FIRSTLINE          ! true: just read first line of data
        LOGICAL, SAVE :: LFLAG    = .FALSE. ! true: output daily/hourly inv in local time
        LOGICAL, SAVE :: SFLAG              ! true: use daily total from hourly
        LOGICAL       :: WARNOUT  = .FALSE. ! true: output warnings
        LOGICAL, SAVE :: YFLAG    = .FALSE. ! true: year mismatch found

        CHARACTER( 6 )   CYYMMDD            ! tmp YYMMDD
        CHARACTER(256) :: MESG    = ' '     ! message buffer

        CHARACTER(BLRLEN3) BLID       ! tmp boiler ID
        CHARACTER(BLRLEN3) PBLID      ! previous boiler ID
        CHARACTER(IOVLEN3) CBUF       ! tmp variable name
        CHARACTER(FIPLEN3) CFIP       ! tmp co/st/cy code
        CHARACTER(ORSLEN3) CORS       ! tmp DOE plant ID
        CHARACTER(ORSLEN3) PCORS      ! previous DOE plant ID
        CHARACTER(SCCLEN3) TSCC       ! tmp source category code
        CHARACTER(ALLLEN3) CSRC       ! tmp source string
        CHARACTER(PLTLEN3) PLT        ! tmp plant code
        CHARACTER(CHRLEN3) STK        ! tmp stack code

        CHARACTER(16) :: PROGNAME = 'RDCEMPD' !  program name

C***********************************************************************
C   begin body of program RDCEMPD

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

C.............  Give note if file is being read as a daily file
            IF( SFLAG ) THEN
                SFLAG = .FALSE.
                MESG = 'NOTE: Ignoring HOURLY_TO_DAILY setting for ' //
     &                 'reading CEM emissions data'
                CALL M3MSG2( MESG )
            END IF

C.............  Get environment variable for calculating flow rate
            MESG = 'Heat input factor for calculating flow rate ' //
     &             '(ft^3/MMBTU)'
            FLOWFAC = ENVREAL( 'FLOW_RATE_FACTOR', MESG, 0.0, IOS )
            FLOWPOS = MXSPDAT + CODFLAG3
            
            IF( FLOWFAC > 0. ) THEN
                MESG = 'NOTE: Hourly flow rates will be calculated ' //
     &                 'from CEM data'
                CALL M3MSG2( MESG )
            ELSE
                MESG = 'NOTE: Hourly flow rates will not be ' //
     &                 'calculated from CEM data'
                CALL M3MSG2( MESG )
            END IF

C.............  Store location of CEM pollutants in master pollutants array
            CEMPIDX = 0        ! array
            DO V = 1, NCEMPOL
                I = INDEX1( CEMPOLS( V ), NIPPA, EANAM )
                CEMPIDX( V ) = I
            END DO

C.............  Create unique lists with ORIS
            ORISFLAG = .TRUE.
            CALL GENUSLST

C.............  Read CEM summary information
            CALL RDCEMSUM
            
C.............  Get output inventory file names given source category
            CALL GETINAME( CATEGORY, ENAME, ANAME )

C.............  Allocate memory for inventory emissions
            ALLOCATE( EMIS( NSRC, NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMIS', PROGNAME )
            EMIS = BADVAL3  ! array

C.............  Read emissions from inventory file
            CALL RDMAPPOL( NSRC, NIPPA, 1, EANAM, EMIS )

            FIRSTIME = .FALSE.

        END IF

C.........  Set variable status
        EASTAT = 1 
        SPSTAT( MXSPDAT ) = 1

C.............  Allocate memory for bad ORIS IDs
        IF ( .NOT. GETSIZES             .AND. 
     &       .NOT. ALLOCATED( UNFDORS )       ) THEN
            ALLOCATE( UNFDORS( MXUNFDORS ), STAT=IOS )   ! count of sources
            CALL CHECKMEM( IOS, 'UNFDORS', PROGNAME )
        END IF

C.........  For the first call in a loop of files, initialize variables
        IF( FIRSTCALL ) THEN
            MINPTR  = 99999999
            MAXPTR  = 0

C.............  Set time step divisor
            TDIVIDE = 3600 * TSTEP / 10000

C.............  If SDATE and EDATE have been set by the data, then compute 
C               total time steps in the data
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

C.............  Set up counters for number of times this routine has been called
            IF( GETSIZES ) MXCALL = 0
            NCALL = 0
            
        END IF

        IF( GETSIZES ) MXCALL = MXCALL + 1
        NCALL = NCALL + 1

C.........  Get count of records per time step if we're ready
C.........  Every sources that matches the master list of ORIS/boilers
C           will have emissions for every time step
        IF( GETCOUNT ) THEN
            IF( MXPDPT( 1 ) == 0 ) THEN
                DO MASOBPOS = 1, NOBRLIST
                    INVOBPOS = FINDC( OBRLIST( MASOBPOS ), 
     &                                NORISBLR, ORISBLR )
                    IF( INVOBPOS <= 0 ) CYCLE
                    
                    NS = OBSRCNT( INVOBPOS )
                    
                    IF( INVOBPOS > 0 ) THEN
                        DO PTR = 1, NSTEPS
                            MXPDPT( PTR ) = MXPDPT( PTR ) + 
     &                                      NS * ( NIPPA + 1 )
                        END DO
                    END IF
                END DO
            END IF

C.............  Allocate memory to indicate if the CEM data contains entries
C               for a given ORIS/boiler and time step
            IF( .NOT. ALLOCATED( MASLIST ) ) THEN
                ALLOCATE( MASLIST( NOBRLIST,NSTEPS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MASLIST', PROGNAME )
                MASLIST = .FALSE.
            END IF
            
            RETURN
        END IF
            
C.........  Loop through file and read it. 
C           1. In the first section, determine the minimum and maximum date. 
C              Use a reference date to do this.
C           2. In the second section, read and store the data.  When storing
C              data, time step index is computed from the start date/time instead
C              of the reference date/time so that the indexing will work properly.

        IREC = 0
        PCORS = ' '
        FIRSTLINE = .TRUE.
        
        DO         !  Head of period-specific file read loop

C.............  Skip comments or blank lines (also includes #CEM header)
            READ( FDEV, *, END=299 ) MESG
            IF( BLKORCMT( MESG ) ) THEN
                IREC = IREC + 1
                CYCLE
            ELSE
                BACKSPACE( FDEV )
            END IF

C.............  There is no error checking to help speed things up
            READ( FDEV, *, ERR=1001, END=299 ) 
     &            CORS, BLID, CYYMMDD, HH, CEMEMIS( NOXIDX ),
     &            CEMEMIS( SO2IDX ), NOXRATE, OPTIME,
     &            GLOAD, SLOAD, HTINPUT

            IREC = IREC + 1

C.............  Write message about which record number we're on
            IF ( FIRSTLINE .OR. MOD( IREC,500000 ) .EQ. 0 ) THEN

                IF( GETSIZES ) THEN
                    WRITE( MESG, 94010 ) 
     &                 BLANK5// 'Determining dates for record numbers',
     &                 IREC, 'to', IREC + 500000 - 1, '...'
                    CALL M3MSG2( MESG )

                ELSE
                    WRITE( MESG, 94010 ) 
     &                 BLANK5// 'Processing data for record numbers',
     &                 IREC, 'to', IREC + 500000 - 1, '...'
                    CALL M3MSG2( MESG )

                END IF

                FIRSTLINE = .FALSE.

            END IF

C.............  Format strings properly
            CORS = ADJUSTR( CORS )
            BLID = ADJUSTR( BLID )

C.............  Skip entry if not in master list of ORIS/boilers
C.............  If CEM summary and input CEM data are consistent, the only
C               non-matches will be entries skipped by CEMScan
            MASOBPOS = FINDC( CORS // BLID, NOBRLIST, OBRLIST )
            IF( MASOBPOS <= 0 ) CYCLE

C.............  Search for ORIS ID in ORIS/FIPs code list and then get time
C               zone from FIPS/Time zone list.
            INVORPOS = FINDC( CORS, NINVORIS, INVORIS )

C.............  If ORIS ID is not found in list...
            IF( INVORPOS <= 0 ) THEN

                IF( GETSIZES ) THEN

                    IF ( CORS /= PCORS ) MXUNFDORS = MXUNFDORS + 1

                ELSE IF( NUNFDORS > 0 ) THEN
                
C.....................  Check to see if ORIS in in list
                    I = INDEX1( CORS, NUNFDORS, UNFDORS )

C.....................  Add to list if it's not there
                    IF ( I <= 0 ) THEN 
                        NUNFDORS = NUNFDORS + 1
                        UNFDORS( NUNFDORS ) = CORS
                    END IF

C.................  If ORIS is first UNFDORS one
                ELSE
                    NUNFDORS = 1
                    UNFDORS( 1 ) = CORS

                END IF
                PCORS = CORS

C.................  Skip it
                CYCLE
            END IF
            
            IORSMTCH( INVORPOS ) = .TRUE.

C.............  Skip entry if not in inventory list of ORIS/boilers
            INVOBPOS = FINDC( CORS // BLID, NORISBLR, ORISBLR )
            IF( INVOBPOS <= 0 ) CYCLE

            MATCHFLG = .TRUE.

C.............  Set Julian day
            YYMMDD = STR2INT( CYYMMDD )
            YY    = YYMMDD/10000
            MONTH = ( YYMMDD - YY * 10000 ) / 100
            DAY   = YYMMDD - YY * 10000 - MONTH * 100
            YEAR  = YEAR4( YY )

C.............  Ensure that year of data is consistent with base year
            IF ( YEAR .NE. BYEAR ) YFLAG = .TRUE. 

            JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )

C.............  Convert from times 0 through 23 to HHMMSS
            JTIME = HH * 10000

C.............  Convert date and time to output time zone
            CFIP = INVORFP ( INVORPOS )

C.............  Local time shift flag (county-specific)
            IF( .NOT. LFLAG ) THEN

                ZONE = GETTZONE( CFIP )
                DZONE = ZONE - TZONE

C.............  Reset time shift to 0 to correctly compute local time zone
            ELSE
                DZONE  = 0

            END IF

C.............  Convert date and time to output time zone.
            CALL NEXTIME( JDATE, JTIME, DZONE * 10000 )

C.............  Determine time step pointer based on reference time
            PTR = SECSDIFF( RDATE, RTIME, JDATE, JTIME ) / TDIVIDE + 1

C.............  Store minimum time step number as compared to reference
            IF( PTR < MINPTR ) MINPTR = PTR

C.............  Store maximum time step number as compared to reference
            IF( PTR > MAXPTR ) MAXPTR = PTR

C.............  If only getting dates, go to next loop iteration
            IF( GETSIZES ) CYCLE

C.............  Determine time step pointer based on actual start time
            PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME )/ TDIVIDE + 1

C.............  Skip record if it is out of range of output file
C.............  NOTE - this is only useful if reading only part of data
            IF( PTR < 1 .OR. PTR > NSTEPS ) CYCLE

C.............  Store indicator for this ORIS/boiler and time step
            MASLIST( MASOBPOS, PTR ) = .TRUE.

C.............  Get number of inventory sources for this ORIS/boiler
            NS = OBSRCNT( INVOBPOS )
            S1 = OBSRCBG( INVOBPOS )

            ALLZERO = .FALSE.

C.............  Compute factor for calculate hourly emissions; check that
C               hourly value is valid
            IF( ANNHEAT( MASOBPOS ) > 0. .AND. HTINPUT > 0. ) THEN
                ANNFAC = HTINPUT / ANNHEAT( MASOBPOS )
            ELSE IF( ANNSLOAD( MASOBPOS ) > 0. .AND. SLOAD > 0. ) THEN
                ANNFAC = SLOAD / ANNSLOAD( MASOBPOS )
            ELSE IF( ANNGLOAD( MASOBPOS ) > 0. .AND. GLOAD > 0. ) THEN
                ANNFAC = GLOAD / ANNGLOAD( MASOBPOS )
            ELSE
                IF( CEMEMIS( NOXIDX ) <= 0. .AND. CEMEMIS( SO2IDX ) <= 0. ) THEN
                    ALLZERO = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'WARNING: All emissions set to 0. for ORIS: ' //
     &              CORS // ' Boiler: ' // BLID // CRLF() // BLANK10 //
     &              'for date: ', YYMMDD, 'and hour: ', HH,
     &              'since hourly heat input,' // CRLF() // BLANK10 //
     &              'gross load, and steam load are missing'
                    CALL M3MESG( MESG )
                END IF
            END IF

            IF( ALLZERO ) THEN
                DO I = 0, NS - 1
                    S = OBSRCNM( S1 + I )
                    
                    DO V = 1, NIPPA
                        EMISVAL = 0.
                        CALL STORE_SOURCE_DATA( S, PTR, V, EMISVAL )
                    END DO
                    
                    FLOWVAL = 0.0
                    CALL STORE_SOURCE_DATA( S, PTR, FLOWPOS, FLOWVAL )

                END DO
                
                CYCLE
            END IF

C.............  Compute emissions from rate and heat input and adjust
            CEMEMIS( NOXIDX ) = CEMEMIS( NOXIDX ) * LB2TON
            CEMEMIS( SO2IDX ) = CEMEMIS( SO2IDX ) * LB2TON

C.............  Record needed data for this source and time step
            DO V = 1, NIPPA

C.................  Check annual emissions for each source; annual emissions
C                   are always used to calculate hourly emissions, so check
C                   for any missing values
                DO I = 0, NS - 1
                    S = OBSRCNM( S1 + I )
                    
                    IF( EMIS( S,V ) < 0. ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Missing "' //
     &                      TRIM( EANAM( V ) ) // '" emissions ' //
     &                      'for source number ', S
                        CALL M3MESG( MESG )
                    END IF
                END DO
                
                IF( EFLAG ) CYCLE

C.................  Check if pollutant is CEM pollutant
                CEMPOL = .FALSE.
                DO N = 1, NCEMPOL
                    IF( CEMPIDX( N ) == V ) THEN

C.........................  Only use CEM data for NOX and SO2 if annual emissions are non-zero
                        IF( N == SO2IDX ) THEN
                          IF( ANNSO2( MASOBPOS ) > 0. ) THEN
                            CEMPOL = .TRUE.
                            IF( CEMEMIS( N ) < 0. ) THEN
                              CEMEMIS( N ) = 0.
                              WRITE( MESG,94010 ) 
     &                          'WARNING: Missing SO2 emissions for ' //
     &                          'ORIS: ' // CORS // ' Boiler: ' // 
     &                          BLID // CRLF() // BLANK10 // 
     &                          'for date: ', YYMMDD, 'and hour: ', HH,
     &                          'reset to 0.'
                              CALL M3MESG( MESG )
                            END IF
                          END IF
                        ELSE IF( N == NOXIDX ) THEN
                          IF( ANNNOX( MASOBPOS ) > 0. ) THEN
                            CEMPOL = .TRUE.
                            IF( CEMEMIS( N ) < 0. ) THEN
                              CEMEMIS( N ) = 0.
                              WRITE( MESG,94010 ) 
     &                          'WARNING: Missing NOX emissions for ' //
     &                          'ORIS: ' // CORS // ' Boiler: ' // 
     &                          BLID // CRLF() // BLANK10 // 
     &                          'for date: ', YYMMDD, 'and hour: ', HH,
     &                          'reset to 0.'
                              CALL M3MESG( MESG )
                            END IF
                          END IF
                        ELSE
                          CEMPOL = .TRUE.
                        END IF
                        
                        EXIT
                    END IF
                END DO

                IF( CEMPOL ) THEN

C.....................  Create weights for allocating CEM emissions to the
C                       sources based on the emissions in the annual inventory
                    DENOM = 0.

                    DO I = 0, NS - 1
                        S = OBSRCNM( S1 + I )
                        DENOM = DENOM + EMIS( S,V )
                    END DO

C.....................  If annual emissions are all zero, use uniform weights
                    IF( DENOM == 0. ) THEN
                        DO I = 0, NS - 1
                            S = OBSRCNM( S1 + I )
                            EMIS( S,V ) = 1.
                            DENOM = DENOM + EMIS( S,V )
                        END DO
                    END IF

                    DENOM = 1. / DENOM
                END IF

C.................  Loop through sources that match this ORIS/boiler
C.................  Store hourly emissions using weights for CEM pollutants;
C                   otherwise calculate hourly emissions using annual CEM factor
                DO I = 0, NS - 1
                    S = OBSRCNM( S1 + I )

                    IF( CEMPOL ) THEN
                        EMISVAL = CEMEMIS( N ) * EMIS( S,V ) * DENOM
                    ELSE
                        EMISVAL = ANNFAC * EMIS( S,V )
                    END IF

                    CALL STORE_SOURCE_DATA( S, PTR, V, EMISVAL )

                END DO  ! end loop over sources that match ORIS/boiler

            END DO      ! end loop over pollutants

C.............  Calculate and store flow rate
            IF( FLOWFAC > 0. .AND. HTINPUT > 0. ) THEN
                FLOWVAL = FLOWFAC * HTINPUT * HR2SEC * FT2M3

C.................  Calculate rate for each stack corresponding to ORIS/boiler
                MXSTKS = OBSRCNT( INVOBPOS )

                IF( ALLOCATED( STKNAM ) ) DEALLOCATE( STKNAM, STKVAL )
                ALLOCATE( STKNAM( MXSTKS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STKNAM', PROGNAME )
                ALLOCATE( STKVAL( MXSTKS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STKVAL', PROGNAME )
                STKNAM = ' '
                STKVAL = 0.

                NSTK = 0
                DO I = 0, NS - 1
                    S = OBSRCNM( S1 + I )
                    CSRC = CSOURC( S )
                    PLT = CSRC( PTBEGL3( 2 ):PTENDL3( 2 ) )
                    STK = CSRC( PTBEGL3( JSTACK ):PTENDL3( JSTACK ) )

                    J = INDEX1( PLT // STK, NSTK, STKNAM )

                    IF( J > 0 ) THEN
                        STKVAL( J ) = STKVAL( J ) + FLOWVAL
                    ELSE
                        NSTK = NSTK + 1
                        STKNAM( NSTK ) = PLT // STK
                        STKVAL( NSTK ) = FLOWVAL
                    END IF
                END DO

                DO I = 0, NS - 1
                    S = OBSRCNM( S1 + I )
                    CSRC = CSOURC( S )
                    PLT = CSRC( PTBEGL3( 2 ):PTENDL3( 2 ) )
                    STK = CSRC( PTBEGL3( JSTACK ):PTENDL3( JSTACK ) )

                    J = INDEX1( PLT // STK, NSTK, STKNAM )
                    
                    CALL STORE_SOURCE_DATA( S, PTR, FLOWPOS, 
     &                                      STKVAL( J ) )
                END DO

            ELSE                     ! Store annfac to LAY1F special variable for AERMOD support
                DO I = 0, NS - 1

                    S = OBSRCNM( S1 + I )
                    CALL STORE_SOURCE_DATA( S, PTR, FLOWPOS, ANNFAC )

                END DO
            END IF

        END DO          ! end loop over lines in file

299     CONTINUE   ! Exit from read loop

        IF ( .NOT. MATCHFLG ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No matches found between CEM data and ' //
     &             'inventory'
            CALL M3MSG2( MESG )
        END IF

        IF ( YFLAG ) THEN
            WRITE( MESG,94010 ) 'WARNING: Some excluded CEM records '//
     &             'contained years other than the base year', BYEAR
            CALL M3MSG2( MESG )
        END IF

C.........  Abort if error found while reading file
        IF( EFLAG ) THEN
            MESG = 'Problem processing CEM data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If this is the last time through this routine, make sure that
C           all ORIS/boilers are accounted for
        IF( .NOT. GETSIZES .AND. 
     &      .NOT. GETCOUNT .AND. 
     &      NCALL == MXCALL      ) THEN
            
            DO MASOBPOS = 1, NOBRLIST
                DO PTR = 1, NSTEPS
                    IF( .NOT. MASLIST( MASOBPOS,PTR ) ) THEN
                        INVOBPOS = FINDC( OBRLIST( MASOBPOS ), 
     &                                    NORISBLR, ORISBLR )
                        IF( INVOBPOS <= 0 ) CYCLE
                        
                        NS = OBSRCNT( INVOBPOS )
                        S1 = OBSRCBG( INVOBPOS )
                        
                        DO I = 0, NS - 1
                            S = OBSRCNM( S1 + I )
                            
                            DO V = 1, NIPPA
                                EMISVAL = 0.
                                CALL STORE_SOURCE_DATA( S, PTR, 
     &                                                  V, EMISVAL )
                            END DO
                            
                            FLOWVAL = 0.
                            CALL STORE_SOURCE_DATA( S, PTR, FLOWPOS, 
     &                                              FLOWVAL )

                        END DO
                    END IF
                END DO
            END DO
            
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
        T = 1
        IF( LFLAG ) T = 0      ! add additional hour when outputing in local time to overwrap next day
        DO I = 1, MAXPTR - T
            CALL NEXTIME( EDATE, ETIME, TSTEP )
        END DO

        RETURN

1001    WRITE( MESG,94010 ) 'Problem reading CEM data at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
  
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( I6.6 )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS
        
            SUBROUTINE STORE_SOURCE_DATA( SRC, PTR, VAR, VAL )
            
C.............  Subroutine arguments
            INTEGER, INTENT(IN) :: SRC
            INTEGER, INTENT(IN) :: PTR
            INTEGER, INTENT(IN) :: VAR
            REAL,    INTENT(IN) :: VAL
            
            INTEGER  HRSRC
            
C----------------------------------------------------------------------

            LPDSRC( SRC ) = .TRUE.
            NPDPT ( PTR ) = NPDPT( PTR ) + 1
            
            HRSRC = NPDPT( PTR )

            IF( HRSRC <= MXPDSRC ) THEN
            
                IDXSRC( HRSRC,PTR ) = HRSRC
                SPDIDA( HRSRC,PTR ) = SRC
                CODEA ( HRSRC,PTR ) = VAR
                EMISVA( HRSRC,PTR ) = VAL
                
            END IF

            RETURN
            
            END SUBROUTINE STORE_SOURCE_DATA

        END SUBROUTINE RDCEMPD
