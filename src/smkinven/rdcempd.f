
        SUBROUTINE RDCEMPD( FDEV, TZONE, TSTEP, MXPDSRC, GETSIZES, 
     &                      GETCOUNT, FIRSTCALL, DAYFLAG, SDATE, STIME, 
     &                      EDATE, ETIME, EASTAT )

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
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: ORISFLAG, NINVORIS, INVORIS, NORISPNT,
     &                      ORISPNT, NORISBLR, ORISBLR, IORSMTCH,
     &                      INVORFP, OBSRCNT, OPSRCNT, OBSRCBG,
     &                      OPSRCBG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, EANAM, CATEGORY, BYEAR, NSRC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: NUNFDORS, UNFDORS, MXPDPT, NPDPT, LPDSRC,
     &                      IDXSRC, SPDIDA, CODEA, EMISVA

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  EXTERNAL FUNCTIONS
        LOGICAL      CHKREAL
        CHARACTER(2) CRLF
        LOGICAL      ENVYN
        INTEGER      FINDC
        INTEGER      GETTZONE
        INTEGER      INDEX1
        INTEGER      JULIAN
        INTEGER      SECSDIFF
        INTEGER      STR2INT
        REAL         STR2REAL
        INTEGER      YEAR4

        EXTERNAL     CHKREAL, CRLF, ENVYN, FINDC, GETTZONE, INDEX1, 
     &               JULIAN, SECSDIFF, STR2INT, STR2REAL, YEAR4

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
        LOGICAL, INTENT(OUT) :: EASTAT( NIPPA ) ! true: pol/act appears in data

C...........   Local parameters
        INTEGER, PARAMETER :: NCEMPOL  = 3   ! number of pollutants in CEM data
        INTEGER, PARAMETER :: CO2IDX   = 1   ! index to CEMPOLS for CO2
        INTEGER, PARAMETER :: SO2IDX   = 2   ! index to CEMPOLS for SO2
        INTEGER, PARAMETER :: NOXIDX   = 3   ! index to CEMPOLS for NOX

        CHARACTER(IOVLEN3), PARAMETER :: CEMPOLS( NCEMPOL ) = 
     &                                      ( / 'CO2             '
     &                                        , 'SO2             '
     &                                        , 'NOX             ' / )

C...........   Local allocatable
        REAL, ALLOCATABLE, SAVE :: EMIS( :,: ) ! inven emissions for CEM pols

C...........   Temporary read arrays
        INTEGER, SAVE :: CEMPIDX( NCEMPOL )   ! Index from CEMPOLS to EASTAT
        REAL             CEMEMIS( NCEMPOL )   ! CEM emissions for line

C.........  File names and unit numbers
        CHARACTER(IOVLEN3) :: ENAME  ! emis i/o api inven logical name
        CHARACTER(IOVLEN3) :: ANAME  ! emis ASCII inven logical name

C...........   Other local variables
        INTEGER          HS, I, J, L, L1, L2, S, S1, S2, T, V  ! counters and indices

        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER          FIP              ! tmp co/st/cy code
        INTEGER          HH               ! 2-digit hour (0-23) from CEM data
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER          JDATE            ! tmp Julian date
        INTEGER          JTIME            ! tmp HHMMSS time
        INTEGER, SAVE :: LOOPNO = 0       ! no. of loops
        INTEGER, SAVE :: MAXPTR           ! maximum time step reference pointer
        INTEGER, SAVE :: MINPTR           ! minimum time step reference pointer
        INTEGER          MONTH            ! tmp month number
        INTEGER, SAVE :: MXUNFDORS = 0    ! max no. of bad ORIS IDs
        INTEGER, SAVE :: NFIELD = 0       ! number of data fields
        INTEGER, SAVE :: NFM1   = 0       ! number of data fields minus 1
        INTEGER          NS               ! tmp no. sources per ORIS
        INTEGER, SAVE :: NSTEPS = 0       ! number of time steps
        INTEGER, SAVE :: NRECS  = 0       ! number of CEM records
        INTEGER          PFIP             ! previous iteration FIPS code
        INTEGER          PTR              ! tmp time step pointer
        INTEGER       :: RDATE = 1980001  ! reference date: Jan 1, 1980
        INTEGER       :: RTIME = 0        ! reference time
        INTEGER, SAVE :: SDATESAV = 0     ! saved start date
        INTEGER, SAVE :: STIMESAV = 0     ! saved start time
        INTEGER, SAVE :: TDIVIDE  = 1     ! time step divisor
        INTEGER          WD               ! tmp field width
        INTEGER          YEAR             ! 4-digit year
        INTEGER          YY               ! 2-digit year
        INTEGER          YYMMDD           ! year, month, and day
        INTEGER          ZONE             ! source time zones

        REAL             DENOM            ! denominator for assignment weighting
        REAL             GLOAD            ! unused CEM data field
        REAL             HTINPUT          ! tmp heat input
        REAL             NSDV1            ! no. src. per ORIS / 1.
        REAL             OPTIME           ! tmp operating time (unused)
        REAL             SLOAD            ! unused CEM data field

        LOGICAL       :: BLRFLAG          ! true: matched rec to inv oris/blr
        LOGICAL, SAVE :: DFLAG = .FALSE.  ! true: dates have been set by the data
        LOGICAL       :: EFLAG = .FALSE.  ! TRUE iff ERROR
        LOGICAL       :: MATCHFLG = .FALSE.! true: a CEM/inventory match found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called
        LOGICAL, SAVE :: SFLAG            ! true: use daily total from hourly
        LOGICAL       :: WARNOUT = .FALSE.! true: then output warnings
        LOGICAL, SAVE :: YFLAG  = .FALSE. ! true: year mismatch found

        CHARACTER(5)      BBUF             ! tmp boiler ID from CEM data
        CHARACTER(256) :: MESG   = ' '     ! message buffer

        CHARACTER(BLRLEN3) BLID       ! tmp boiler ID with inventory length
        CHARACTER(BLRLEN3) PBLID      ! previous boiler ID
        CHARACTER(IOVLEN3) CBUF       ! tmp variable name
        CHARACTER(FIPLEN3) CFIP       ! tmp co/st/cy code
        CHARACTER(CHRLEN3) PNT        ! tmp characteristic 1
        CHARACTER(ORSLEN3) CORS       ! tmp DOE plant ID
        CHARACTER(ORSLEN3) PCORS      ! previous DOE plant ID
        CHARACTER(SCCLEN3) TSCC       ! tmp source category code
        CHARACTER(ALLLEN3) CSRC       ! tmp source string

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

C.............  Give note if file is being read as a daily file
            IF( SFLAG ) THEN
                SFLAG = .FALSE.
                MESG = 'NOTE: Ignoring HOURLY_TO_DAILY setting for ' //
     &                 'reading CEM emissions data'
                CALL M3MSG2( MESG )
            END IF

C.............  Allocate memory for bad source storage
c            ALLOCATE( BADSRC( NSRC ), STAT=IOS )
c            CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

C.............  Set status of input variables based purely on what we know is 
C               contained in the input CEM data (NOx, CO2, and SO2)
            CEMPIDX = 0        ! array
            DO V = 1, NCEMPOL
                I = INDEX1( CEMPOLS( V ), NIPPA, EANAM )

                IF ( I .GT. 0 ) EASTAT( I ) = .TRUE.
                CEMPIDX( V ) = I

            END DO

C.............  Create unique lists with ORIS
            ORISFLAG = .TRUE.
            CALL GENUSLST

C.............  Get output inventory file names given source category
            CALL GETINAME( CATEGORY, ENAME, ANAME )

C.............  Allocate memory for inventory emissions for creating weights
            ALLOCATE( EMIS( NSRC, NCEMPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMIS', PROGNAME )
            EMIS = BADVAL3  ! array

C.............  Read emissions from inventory file
            DO V = 1, NCEMPOL

                IF ( CEMPIDX( V ) .LE. 0 ) CYCLE

                CBUF = EANAM ( CEMPIDX( V ) )

                CALL RDMAPPOL( NSRC, 1, 1, CBUF, EMIS( 1,V ) )

            END DO

            FIRSTIME = .FALSE.

        END IF

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

        END IF

C.........  Loop through file and read it. In the first section, determine
C           the minimum and maximum date. Use a reference date to do this. In
C           the second section, determine the number of records per time 
C           step. In the third section, read and store the data.  When storing
C           data, time step index is computed from the start date/time instead
C           of the reference date/time so that the indexing will work properly.

        IREC = 0
        PCORS = ' '
        
C.........  Skip #CEM header line if it is present
        READ( FDEV, * ) MESG
        IF( INDEX( MESG, '#CEM' ) < 0 ) THEN
            BACKSPACE( FDEV )
        ELSE
            IREC = IREC + 1
        END IF
        
c        TDAT = 0   !  array
        DO         !  Head of period-specific file read loop

C.............  There is no error checking to help speed things up
            READ( FDEV, *, ERR=1001, END=299 ) 
     &            CORS, BBUF, YYMMDD, HH, CEMEMIS( CO2IDX ),
     &            CEMEMIS( SO2IDX ), CEMEMIS( NOXIDX ), OPTIME,
     &            GLOAD, SLOAD, HTINPUT
     
            IREC = IREC + 1

C.............  Write message about which record number we're on
            IF ( IREC <= 2 .OR. MOD( IREC,500000 ) .EQ. 0 ) THEN

                IF( NRECS .GT. 0 .AND. GETCOUNT ) THEN
                    WRITE( MESG, 94010 ) 
     &                 BLANK5// 'Counting record numbers', IREC, 
     &                 'to', MIN( IREC + 500000 - 1, NRECS ), '...'
                    CALL M3MSG2( MESG )

                ELSE IF ( NRECS .GT. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                 BLANK5// 'Processing record numbers', IREC, 
     &                 'to', MIN( IREC + 500000 - 1, NRECS ), '...'
                    CALL M3MSG2( MESG )

                END IF

            END IF

C.............  Initialize flag to indicate oris/boiler match or oris/point
            BLRFLAG = .FALSE.

C.............  Format strings properly. Assumed BLRLEN3 < 5.
            CORS = ADJUSTR( CORS )
            BBUF = ADJUSTL( BBUF )
            BLID = BBUF( 1:BLRLEN3 )
            BLID = ADJUSTR( BLID )
            PNT  = BBUF
            PNT  = ADJUSTR( PNT )

C.............  Compute emissions from rate and heat input and adjust
C               units from pounds to short tons
            CEMEMIS( CO2IDX ) = CEMEMIS( CO2IDX ) * LB2TON
            CEMEMIS( NOXIDX ) = CEMEMIS( NOXIDX ) * HTINPUT * LB2TON
            CEMEMIS( SO2IDX ) = CEMEMIS( SO2IDX ) * LB2TON

C.............  Set Julian day from MMDDYY6 SAS format
            YY    = YYMMDD/10000
            MONTH = ( YYMMDD - YY * 10000 ) / 100
            DAY   = YYMMDD - YY * 10000 - MONTH * 100
            YEAR  = YEAR4( YY )

C.............  Ensure that year of data is consistent with base year
            IF ( YEAR .NE. BYEAR ) YFLAG = .TRUE. 

            JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )

C.............  Skip dates that are outside the range allowed
            IF ( JDATE .LT. SDATE .AND. JDATE .GT. EDATE ) CYCLE

C.............  Convert from times 0 through 23 to HHMMSS
            JTIME = HH * 10000

C.............  Search for ORIS ID in ORIS/FIPs code list and then get time
C               zone from FIPS/Time zone list.
            I = FINDC( CORS, NINVORIS, INVORIS )

C.............  If ORIS ID is not found in list...
            IF( I .LE. 0 ) THEN

                IF( GETSIZES ) THEN

                    IF ( CORS .NE. PCORS ) MXUNFDORS = MXUNFDORS + 1

                ELSE IF( NUNFDORS .GT. 0 ) THEN
C.....................  Check to see if ORIS in in list
                    J = INDEX1( CORS, NUNFDORS, UNFDORS )

C.....................  Add to list if it's not there
                    IF ( J .LE. 0 ) THEN 
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

C.............  Otherwise, get FIPs code and time zone
            ELSE
                MATCHFLG = .TRUE.
                FIP   = INVORFP( I )
                ZONE  = GETTZONE( FIP )
            END IF
 
C.............  Convert date and time to output time zone.
            CALL NEXTIME( JDATE, JTIME, ( ZONE - TZONE ) * 10000 )

C.............  Determine time step pointer based on reference time
            PTR = SECSDIFF( RDATE, RTIME, JDATE, JTIME ) / TDIVIDE + 1

C.............  Store minimum time step number as compared to reference
            IF( PTR .LT. MINPTR ) MINPTR = PTR

C.............  Store maximum time step number as compared to reference
            IF( PTR .GT. MAXPTR ) MAXPTR = PTR

C.............  If only getting dates, go to next loop iteration
            IF( GETSIZES ) CYCLE

C.............  Determine time step pointer based on actual start time
            PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME )/ TDIVIDE + 1

C.............  Skip record if it is out of range of output file
C.............  NOTE - this is only useful if reading only part of data
            IF( PTR. LT. 1 .OR. PTR .GT. NSTEPS ) CYCLE

C.............  Try to match using ORIS and points (from IDA only!)
            J = FINDC( CORS // PNT, NORISPNT, ORISPNT )

C.................  If no match, search for this record in list of 
C                   ORIS // boilers
            IF ( J .LE. 0 ) THEN
                J = FINDC( CORS // BLID, NORISBLR, ORISBLR )

C.................  Skip sources that don't match
                IF ( J .LE. 0 ) CYCLE
                
                BLRFLAG = .TRUE.
                NS = OBSRCNT( J )  ! number of sources from oris/boiler
                IORSMTCH( I ) = .TRUE.

C.............  Get number of sources for this ORIS//point
            ELSE
                NS = OPSRCNT( J )  ! number of sources from oris/point
                IORSMTCH( I ) = .TRUE.

            END IF

C.............  Count estimated source count per time step. Actual may not
C               include all pollutants, but this should be
            MXPDPT( PTR ) = MXPDPT( PTR ) + NS * NCEMPOL

C.............  If only counting records per time step, go to next loop
C               iteration
            IF( GETCOUNT ) CYCLE

C.............  Get source list from those records that do match
            IF ( BLRFLAG ) THEN
                S1 = OBSRCBG( J )  ! first source
            ELSE
                S1 = OPSRCBG( J )  ! first source
            END IF

            S2 = S1 + NS - 1   ! last source

C.............  Record needed data for this source and time step
            T = PTR
            DO V = 1, NCEMPOL

C.................  Skip CEM pollutants that are not in the inventory
                IF ( CEMPIDX( V ) .LE. 0 ) CYCLE

C.............  Create weights for allocating CEM emissions to
C               these sources based on the emissions in the annual inventory.
C.............  Make sure that if the emissions are zero, weights are uniform.

                DENOM = SUM( EMIS( S1:S2,V ) )
                IF ( DENOM .EQ. 0. ) THEN
                    EMIS( S1:S2,V ) = 1.  ! array
                    DENOM = SUM( EMIS( S1:S2,V ) )

                ELSE IF ( DENOM .LT. 0. ) THEN
                    EFLAG = .TRUE.
                    CBUF = EANAM( CEMPIDX( V ) )
                    L = LEN_TRIM( CBUF ) 
                    WRITE( MESG,94010 ) 'ERROR: Missing "' // 
     &                     CBUF( 1:L ) // '" emissions at or ' //
     &                     'near source', S
                    CALL M3MESG( MESG )

                END IF

                DENOM = 1. / DENOM


C.................  Loop through sources that match this ORIS/boiler
C.................  Store hourly emissions using weights if multiple sources
C                   match the ORIS/boiler. Weights are based on annual
C                   inventory emissions.
                DO S = S1, S2

                    LPDSRC( S ) = .TRUE.
                    NPDPT( T ) = NPDPT( T ) + 1

                    HS = NPDPT( T )

                    IF( HS .LE. MXPDSRC ) THEN

                        IDXSRC( HS,T ) = HS
                        SPDIDA( HS,T ) = S
                        CODEA ( HS,T ) = CEMPIDX(V)
                        IF( CEMEMIS( V ) .LT. 0 ) THEN
                            EMISVA( HS,T ) = BADVAL3
                        ELSE
                            EMISVA( HS,T )= CEMEMIS(V)* EMIS(S,V)* DENOM
                        END IF

                    END IF

                END DO  ! end loop over sources that match ORID/Boiler ID

            END DO      ! end loop over CEM pollutants

        END DO          ! end loop over lines in file

299     CONTINUE   ! Exit from read loop

        NRECS = IREC

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

1001    WRITE( MESG,94010 ) 'Problem reading CEM data at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
  

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( I6.6 )

        END SUBROUTINE RDCEMPD
