
      SUBROUTINE PROCTPRO( NFLAG, METFLAG, METNAME )

C*******************************************************************************
C    DESCRIPTION
C      Processes (read/sort/filter) GENTPRO-style (ASCII CSV)
C      cross-references and profiles and MET-based profiles.
C
C    PRECONDITIONS REQUIRED:
C      setenv  TREF            <path for cross-reference file>
C      setenv  TPRO_MONTHLY    <path for month-of-year profiles file>
C      setenv  TPRO_WEEKLY     <path for day-of-week   profiles file>
C      setenv  TPRO_DAILY      <path for day-of-month  pfrofiles file>
C      setenv  TPRO_HOUR       <path for M3IO met based        profiles file>
C
C    INTERNAL SUBROUTINES AND FUNCTIONS:
C      TIMEFLAGS() sets active month/day-of-week flags and SDATE:EDATE
C      for the set of episodes being run
C
C      CSVPROF() opens and counts CSV-profile files, allocates arguments,
C      and filters and sorts the input table
C
C      CSVDOMP() opens and counts TPRO_DAILY day-of-month CSV-profile files,
C      allocates arguments, and filters and sorts the input table
C
C      CSVOPEN() opens and counts lines and data-lines in the
C      CSV-file FNAME, rewinding after completion
C
C      SORTPRO() sorts profiles for CSVPRROF() and CSVDOMP()
C
C      SORTREF() allocates sorted XREF arrays and sorts XREF tables
C      into them.
C
C      ISLEAP()  is TRUE iff its argument date is a leap-year
C
C    REVISION HISTORY:
C        Adapted 7/2014 by Carlie J. Coats, Jr. from SMOKE 3.5.1
C        RDTREF(), RDTPROF(), and ASGNTPRO() for new GENTPRO CSV profiles
C        and cross-references
C*****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004-2014, Environmental Modeling for Policy Development
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
C.........  MODINFO  contains the information about the source category
C.........  MODSOURC contains the inventory arrays
C.........  MODTMPRL contains the temporal profile tables

        USE M3UTILIO

        USE MODINFO,  ONLY: CATEGORY, NIPPA, EANAM, NSRC

        USE MODSOURC, ONLY: CIFIP, CSCC, TZONES

        USE MODTMPRL, ONLY: NMON,   NWEK,   NHRL,   NDOM,
     &                      MONFAC, WEKFAC, HRLFAC, DOMFAC,
     &                      MONREF, WEKREF, HRLREF, DAYREF,
     &                      METPROTYPE, METREFFLAG, HOUR_TPROF,
     &                      MONFLAG, DAYFLAG, WKDFLAG, WKEFLAG,
     &                      ITDATE, RUNLEN, POLREFFLAG,
     &                      MTHIDP, WEKIDP, DOMIDP, HRLIDP,
     &                      MTHCOUNT, WEKCOUNT, DOMCOUNT, MONCOUNT,
     &                      TUECOUNT, WEDCOUNT, THUCOUNT, FRICOUNT,
     &                      SATCOUNT, SUNCOUNT, METCOUNT,
     &                      MTHKEYS, WEKKEYS, DOMKEYS, MONKEYS,
     &                      TUEKEYS, WEDKEYS, THUKEYS, FRIKEYS,
     &                      SATKEYS, SUNKEYS, METKEYS,
     &                      MTHPDEX, WEKPDEX, DOMPDEX, MONPDEX,
     &                      TUEPDEX, WEDPDEX, THUPDEX, FRIPDEX,
     &                      SATPDEX, SUNPDEX, NMETPROF, METPROF,
     &                      MTHPROF, WEKPROF, DOMPROF, HRLPROF

        IMPLICIT NONE

C.........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters
C        INCLUDE 'FDESC3.EXT'    !  i/o api file-header structures
C        INCLUDE 'IODECL3.EXT'   !  i/o api subroutine declarations

C.........   Arguments:

        LOGICAL,       INTENT(IN   ) :: NFLAG       !  no xref:  constant profiles
        LOGICAL,       INTENT(  OUT) :: METFLAG     !  use hour-specific met-based profiles
        CHARACTER(16), INTENT(  OUT) :: METNAME     !  logical name for hour-specific file

C.........   EXTERNAL FUNCTIONS and their descriptions:

C       CHARACTER(2), EXTERNAL :: CRLF

        LOGICAL, EXTERNAL :: CHKINT
        LOGICAL, EXTERNAL :: BLKORCMT
C       INTEGER, EXTERNAL :: ENVINT
C       LOGICAL, EXTERNAL :: ENVYN
C       INTEGER, EXTERNAL :: FIND1
C       INTEGER, EXTERNAL :: FINDC
        INTEGER, EXTERNAL :: GETNLIST
        INTEGER, EXTERNAL :: GETFLINE
C       INTEGER, EXTERNAL :: SECSDIFF
C       INTEGER, EXTERNAL :: INDEX1
        LOGICAL, EXTERNAL :: SETSCCTYPE
C       INTEGER, EXTERNAL :: STR2INT

C.........   Local parameters

        INTEGER, PARAMETER :: AREATYP  =  1
        INTEGER, PARAMETER :: MOBILTYP =  2
        INTEGER, PARAMETER :: POINTTYP =  3
        INTEGER, PARAMETER :: MXTCOL   = 15

        CHARACTER( 1),      PARAMETER :: BLANK   = ' '
        CHARACTER(16),      PARAMETER :: PNAME   = 'PROCTPRO'                   !  subroutine name
        CHARACTER(24),      PARAMETER :: ZEROS   = '000000000000000000000000'   !  "all-zeros"
        CHARACTER(24),      PARAMETER :: CMISS   = '????????????????????????'   !  "not a legal string-entry"

        CHARACTER( 6), PARAMETER :: LOCCATS( 3 ) = ( / 'AREA  ', 'MOBILE', 'POINT ' / )

C.........   Local variables

        INTEGER     F, I, J, J1, J2, J3, J4, J5, JJ, K, L, M, N, NN, S, V, W    !  counters and indices
        INTEGER     ISTAT
        INTEGER     SDATE, STIME, EDATE !  starting, ending date for this set of episodes
        INTEGER     NLINES, NDATA       !  from CSVOPEN()

        INTEGER     ICAT    !  category subscript in LOCCATS

        INTEGER     XDEV    !  unit number for TREF file

        INTEGER     COD     !  temporary pollutant/emission type code
        INTEGER     FIP     !  temporary FIPS code
        INTEGER     IDIU    !  temporary diurnal profile code
        INTEGER     IDUM    !  tmp dummy integer
        INTEGER     IMON    !  temporary monthly profile code
        INTEGER     IOS     !  i/o status
        INTEGER     IWEK    !  temporary weekly profile code
        INTEGER     IREC    !  record counter
        INTEGER     JS      !  position of SCC in source chars in x-ref file
        INTEGER     JSPC    !  tmp index to master pollutant/etype list
        INTEGER     LINTYPE !  temporary source category code
        INTEGER     LPCK    !  length of point definition packet
        INTEGER     NCP     !  input point source header parm
        INTEGER     NFIELD  !  tmp number of fields in LINE
        INTEGER     NREF    !  number of x-ref entries before filtering
        INTEGER     NXREF   !  number of valid x-ref entries
        INTEGER     RDT     !  temporary road class code
        INTEGER     TMON    !  temporary met-based monthly profile code
        INTEGER     TDAY    !  temporary met-based daily profile code
        INTEGER     THRS    !  temporary met-based hourly profile code
        INTEGER     VTYPE   !  temporary vehicle type number

        LOGICAL     EFLAG   !  true: error occurred
        LOGICAL     AFLAG   !  true: error occurred
        LOGICAL     PFLAG   !  true: pol/act-spec entries skipped
        LOGICAL     SKIPREC !  true: skip this x-ref entry

        INTEGER, ALLOCATABLE :: COUNTIES( : )
        INTEGER, ALLOCATABLE ::  METFIPS( : )

C.........  Unsorted, sorted cross-reference data structures filtered
C.........  from the input file
C.........  Note that week-day, week-end, and all-days references are
C.........  mapped onto the day-specific MON, ... SUN diurnal-profile
C.........  references:
C.........  Doc claims IDs are CHAR*15.  Here, we embed them into CHAR*16
C.........  for better memory alignment

        INTEGER     DAYCOUNT

        CHARACTER(16), ALLOCATABLE :: MTHIDU( : )       !  unsorted
        CHARACTER(16), ALLOCATABLE :: WEKIDU( : )
        CHARACTER(16), ALLOCATABLE :: DOMIDU( : )
        CHARACTER(16), ALLOCATABLE :: MONIDU( : )
        CHARACTER(16), ALLOCATABLE :: TUEIDU( : )
        CHARACTER(16), ALLOCATABLE :: WEDIDU( : )
        CHARACTER(16), ALLOCATABLE :: THUIDU( : )
        CHARACTER(16), ALLOCATABLE :: FRIIDU( : )
        CHARACTER(16), ALLOCATABLE :: SATIDU( : )
        CHARACTER(16), ALLOCATABLE :: SUNIDU( : )
        CHARACTER(16), ALLOCATABLE :: METIDU( : )

        CHARACTER(16), ALLOCATABLE :: MTHIDS( : )       !  sorted
        CHARACTER(16), ALLOCATABLE :: WEKIDS( : )
        CHARACTER(16), ALLOCATABLE :: DOMIDS( : )
        CHARACTER(16), ALLOCATABLE :: MONIDS( : )
        CHARACTER(16), ALLOCATABLE :: TUEIDS( : )
        CHARACTER(16), ALLOCATABLE :: WEDIDS( : )
        CHARACTER(16), ALLOCATABLE :: THUIDS( : )
        CHARACTER(16), ALLOCATABLE :: FRIIDS( : )
        CHARACTER(16), ALLOCATABLE :: SATIDS( : )
        CHARACTER(16), ALLOCATABLE :: SUNIDS( : )
        CHARACTER(16), ALLOCATABLE :: METIDS( : )

        CHARACTER(ALLLEN3), ALLOCATABLE :: MTHKEYU( : )     !  unsorted cross-reference keys
        CHARACTER(ALLLEN3), ALLOCATABLE :: WEKKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: DOMKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: MONKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: TUEKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: WEDKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: THUKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: FRIKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: SATKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: SUNKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: WKEKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: WKDKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: ALDKEYU( : )
        CHARACTER(ALLLEN3), ALLOCATABLE :: METKEYU( : )

        LOGICAL             LEAPYEAR

        CHARACTER(2)        MBUF        !  2-digit month "MM"
        CHARACTER(1)        SCC1        !  1st character of SCC
        CHARACTER(5)        CPOS        !  tmp sorted position of pol/act
        CHARACTER(8)        LASTMET, MISSMET

        CHARACTER(SICLEN3)  CDUM        !  dummy character field for SIC
        CHARACTER(MACLEN3)  CDUM2       !  dummy charecter field for MACT
        CHARACTER(LNKLEN3)  CLNK        !  temporary link code
        CHARACTER(ALLLEN3)  CSRCALL     !  buffer for source char, incl pol/act
        CHARACTER(FIPLEN3)  CFIP        !  buffer for CFIPS code
        CHARACTER(SCCLEN3)  TSCC        !  temporary SCC
        CHARACTER(PLTLEN3)  CPLT        !  tmp plant ID
        CHARACTER(CHRLEN3)  CPNT        !  tmp point ID
        CHARACTER(CHRLEN3)  CSTK        !  tmp stack ID
        CHARACTER(CHRLEN3)  CSEG        !  tmp segment ID
        CHARACTER(CHRLEN3)  CPL5        !  tmp plt char 5
        CHARACTER(IOVLEN3)  CPOA        !  temporary pollutant/emission type
        CHARACTER(RWTLEN3)  CRWT        !  roadway type no.
        CHARACTER(VIDLEN3)  CVID        !  vehicle type ID no.
        CHARACTER(CHRLEN3)  CHARS( 5 )  !  temporary plant characteristics
        CHARACTER(FIPLEN3)  PREVFIP     !  previous region code

        CHARACTER(16)       MTHNAME, WEKNAME, DOMNAME, DIUNAME

        CHARACTER(16)       ANAME
        CHARACTER(16)       THISID, LASTID
        CHARACTER(20)       FIELD( MXTCOL )   !  Array for reading temporal x-ref fields

C......... Data structures for filtering XREF file by actually-occurring
C......... and default FIPS and pollutants.
C......... NOTE:  FLTRXREF() does *not* filter by FIP!

        INTEGER             NFIPKEY
        INTEGER             INDXKEY( NIPPA+2*NSRC+1 )   !  sorting index for pol
        CHARACTER(FIPLEN3)  FIPKEYU( 2*NSRC+1 )         !  unsorted, with duplicates
        CHARACTER(FIPLEN3)  FIPKEYS( 2*NSRC+1 )         !  sorted, duplicate-free

        CHARACTER(256)      MESG
        CHARACTER(512)      LINE

        INTEGER, SAVE :: MXWARN = -9999     !  from env vble SMK_MAXWARNING

C...........   body   ......................................................

        MESG = 'Processing temporal XREFS and PROFs for '//CATEGORY
        CALL M3MESG( MESG )

        IF ( NFLAG ) THEN       !  no time-dependence in emissions profiles

            ALLOCATE( METPROF( NSRC,  NIPPA ),
     &                MTHPROF( NSRC,  NIPPA ),
     &                WEKPROF( NSRC,  NIPPA ),
     &                DOMPROF( NSRC,  NIPPA ),
     &                HRLPROF( NSRC,7,NIPPA ),
     &                 MONFAC( 12,0:1 ),
     &                 WEKFAC(  7,0:1 ),
     &                 HRLFAC( 24,0:1 ), STAT = IOS )
            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' ) 'ERROR:  allocation failure.  STAT=', IOS
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            !  array assignments for uniform profiles:

            METPROF = IMISS3
            MTHPROF = 1
            WEKPROF = 1
            DOMPROF = 1
            HRLPROF = 1
            MONFAC  = 1.0
            WEKFAC  = 1.0
            HRLFAC  = 1.0

            CALL M3MSG2( 'PROCTPRO:  initialized uniform-profile references and profiles' )

            RETURN

        END IF      !  if nflag:  no time-dependence in emissions profiles


        IF ( MXWARN .LT. 0 ) THEN
            MXWARN = ENVINT( WARNSET , BLANK, 100, ISTAT )
            IF ( ISTAT .GT. 0 ) THEN
                MESG = 'ERROR:  bad env vble "SMK_MAXWARNING"'
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF
        END IF

C.........  Ensure that the CATEGORY is valid

        ICAT = INDEX1( CATEGORY, 3, LOCCATS )

        IF ( ICAT .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: category "' // TRIM( CATEGORY ) //
     &             '" is not valid in routine ' // PNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PNAME, 0, 0, ' ', 2 )
        ENDIF

C.........  Sort the FIPS codes and defaults, and store them for use in
C.........  filtering the XREF file:

        N = 1
        INDXKEY( N ) = N        !  ultimate default
        FIPKEYU( N ) = REPEAT( '0' , FIPLEN3)
        DO I = 1, NSRC
            N = N + 1
            INDXKEY( N ) = N
            FIPKEYU( N ) = CIFIP( I )( 1:STALEN3 ) // '000' !  state only (default value)
            N = N + 1
            INDXKEY( N ) = N
            FIPKEYU( N ) = CIFIP( I )
        END DO

        CALL SORTIC( N, INDXKEY, FIPKEYU )

        PREVFIP = ''       !  now construct duplicate-free sorted list:
        M = 0
        DO I = 1, N
            K = INDXKEY( I )
            IF ( FIPKEYU( K ) .NE. PREVFIP ) THEN
                M = M + 1
                FIPKEYS( M ) = FIPKEYU( K )
            END IF
        END DO
        NFIPKEY = M


C.........  Set active day-of-week, month-of-year flags for SDATE:EDATE

        CALL TIMEFLAGS( )


C.........  Open and count TREF file

        CALL M3MSG2( 'Reading temporal cross-reference file...' )
        ANAME = CATEGORY(1:1) // 'TREF'
        XDEV  = CSVOPEN( ANAME, NLINES, NDATA )

C.........  Allocate/initialize scratch data structures

        MTHCOUNT = 0
        WEKCOUNT = 0
        DOMCOUNT = 0
        MONCOUNT = 0
        TUECOUNT = 0
        WEDCOUNT = 0
        THUCOUNT = 0
        FRICOUNT = 0
        SATCOUNT = 0
        SUNCOUNT = 0
        METCOUNT = 0

        ALLOCATE( MTHIDU( NDATA ),
     &            WEKIDU( NDATA ),
     &            DOMIDU( NDATA ),
     &            MONIDU( NDATA ),
     &            TUEIDU( NDATA ),
     &            WEDIDU( NDATA ),
     &            THUIDU( NDATA ),
     &            FRIIDU( NDATA ),
     &            SATIDU( NDATA ),
     &            SUNIDU( NDATA ),
     &            METIDU( NDATA ),
     &           MTHKEYU( NDATA ),
     &           WEKKEYU( NDATA ),
     &           DOMKEYU( NDATA ),
     &           MONKEYU( NDATA ),
     &           TUEKEYU( NDATA ),
     &           WEDKEYU( NDATA ),
     &           THUKEYU( NDATA ),
     &           FRIKEYU( NDATA ),
     &           SATKEYU( NDATA ),
     &           SUNKEYU( NDATA ),
     &           METKEYU( NDATA ),
     &        METREFFLAG( NIPPA ),
     &        POLREFFLAG( NIPPA ),
     &      METPROF( NSRC,NIPPA ),STAT = IOS )
        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG, '( A, I10 )' ) 'ERROR:  xref allocation failure.  STAT=', IOS
            CALL M3EXIT( PNAME, 0,0, MESG, 2 )
        END IF

        EFLAG      = .FALSE.
        PFLAG      = .FALSE.
        METREFFLAG = .FALSE.
        POLREFFLAG = .FALSE.
        METPROF    = IMISS3

C.........  Read and process the TREF file

        LASTMET = CMISS
        MISSMET = LASTMET
        CDUM    = BLANK
        CDUM2   = BLANK
        IDUM    = 0
        M       = 0
        W       = 0

        DO L = 1, NLINES

            READ( XDEV, '( A )', IOSTAT=ISTAT ) LINE
            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10, 1X, A, I10 )' )
     &              'ERROR: reading "TREF" at line', L,
     &              'IOSTAT=', ISTAT
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                W     = W + 1
                IF ( W .GT. MXWARN ) EXIT
                CYCLE
            ELSE IF ( BLKORCMT( LINE ) ) THEN
                CYCLE
            END IF

            CALL PARSLINE( LINE, MXTCOL, FIELD )    !  does adjustl() on all fields

            IF ( FIELD(9) .EQ. BLANK ) THEN
                WRITE( MESG, '( A, I10, 2X, A )' )
     &              'ERROR:  Bad xref in "TREF" at line', L,
     &              'for field 9 (profile-ID)'
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                W     = W + 1
                IF ( W .GT. MXWARN ) EXIT
                CYCLE
            END IF

            TSCC = FIELD(1)
            CFIP = FIELD(2)
            CPOA = FIELD(7)

            CALL FLTRNEG( TSCC )     ! Filter 0 and -9 to blank
            CALL PADZERO( TSCC )     ! Pad with zeros
            CALL FLTRNEG( CFIP )     ! Filter 0 and -9 to blank
            CALL PADZERO( CFIP )     ! Pad with zeros

C.............  Skip lines that are not valid for this FIP

            JJ = FINDC( CFIP, NFIPKEY, FIPKEYS )
            IF ( JJ .LE. 0 )  CYCLE

C.............  Post-process x-ref information to scan for '-9',
C.............  pad with zeros, compare SCC version master list,
C.............  compare SIC version to master list, and compare
C.............  pol/act name with master list.
C.............  NOTE:  FLTRXREF() does *not* filter by FIP,
C.............  and does LINEAR searches for pollutant!

            CALL FLTRXREF( CFIP, CDUM, TSCC, CPOA, CDUM2,
     &                     IDUM, IDUM, JSPC, PFLAG, SKIPREC )

C.............  Skip lines that are not valid for this inven and src cat

            IF ( SKIPREC )  CYCLE

C.............  Write pol/act position to a character string

            IF ( JSPC .EQ. 0 ) THEN
                CPOS = ZEROS
            ELSE
                WRITE( CPOS, '(I5.5)' ) JSPC      ! species index into EANAM, from FLTRXREF
                POLREFFLAG( JSPC ) = .TRUE.
            END IF

            IF ( ICAT .EQ. AREATYP ) THEN

                IF ( FIELD(3) .NE. BLANK )  CYCLE
                CALL BLDCSRC( CFIP, TSCC, BLANK,
     &                        BLANK, BLANK, BLANK,
     &                        BLANK, CPOS, CSRCALL )

            ELSE IF ( ICAT .EQ. POINTTYP ) THEN

                CPLT = FIELD(3)
                CPNT = FIELD(4)
                CSTK = FIELD(5)
                CSEG = FIELD(6)
                CPL5 = TSCC                     !  padded from length=10 to length=15
                CALL BLDCSRC( CFIP, CPLT, CPNT, CSTK, CSEG, CPL5, TSCC, CPOS, CSRCALL )

            ELSE IF ( ICAT .EQ. MOBILTYP ) THEN

C.................   M Houyoux note: TSCC has been put in here instead of road type
C.................  and link has been removed.  These were breaking the county-SCC specific
C.................  assignments by setting CNFIP in xreftbl.f to be non-blank and not the SCC.
C.................  However, this change breaks link-specific profile assignments, which
C.................  are not likely to be used anyway.  I suggest that we just remove
C.................  link-specific assignments from the documentation for Spcmat.

                CALL BLDCSRC( CFIP, TSCC, BLANK,
     &                        BLANK, BLANK, BLANK,
     &                        BLANK, CPOS, CSRCALL )

            END IF


            SELECT CASE( FIELD(8) )

                CASE( BLANK )

                    CYCLE

                CASE( 'MONTHLY' )

                    MTHCOUNT = MTHCOUNT + 1
                    MTHIDU ( MTHCOUNT ) = FIELD(9)
                    MTHKEYU( MTHCOUNT ) = CSRCALL

                CASE( 'WEEKLY' )

                    WEKCOUNT = WEKCOUNT + 1
                    WEKIDU ( WEKCOUNT ) = FIELD(9)
                    WEKKEYU( WEKCOUNT ) = CSRCALL

                CASE( 'DAILY' )

                    DOMCOUNT = DOMCOUNT + 1
                    DOMIDU ( DOMCOUNT ) = FIELD(9)
                    DOMKEYU( DOMCOUNT ) = CSRCALL

                CASE( 'MONDAY' )
                    IF ( .NOT.DAYFLAG(1) ) CYCLE

                    MONCOUNT = MONCOUNT + 1
                    MONIDU ( MONCOUNT ) = FIELD(9)
                    MONKEYU( MONCOUNT ) = CSRCALL

                CASE( 'TUESDAY' )
                    IF ( .NOT.DAYFLAG(2) ) CYCLE

                    TUECOUNT = TUECOUNT + 1
                    TUEIDU ( TUECOUNT ) = FIELD(9)
                    TUEKEYU( TUECOUNT ) = CSRCALL

                CASE( 'WEDNESDAY' )
                    IF ( .NOT.DAYFLAG(3) ) CYCLE

                    WEDCOUNT = WEDCOUNT + 1
                    WEDIDU ( WEDCOUNT ) = FIELD(9)
                    WEDKEYU( WEDCOUNT ) = CSRCALL

                CASE( 'THURSDAY' )
                    IF ( .NOT.DAYFLAG(4) ) CYCLE

                    THUCOUNT = THUCOUNT + 1
                    THUIDU ( THUCOUNT ) = FIELD(9)
                    THUKEYU( THUCOUNT ) = CSRCALL

                CASE( 'FRIDAY' )
                    IF ( .NOT.DAYFLAG(5) ) CYCLE

                    FRICOUNT = FRICOUNT + 1
                    FRIIDU ( FRICOUNT ) = FIELD(9)
                    FRIKEYU( FRICOUNT ) = CSRCALL

                CASE( 'SATURDAY' )
                    IF ( .NOT.DAYFLAG(6) ) CYCLE

                    SATCOUNT = SATCOUNT + 1
                    SATIDU ( SATCOUNT ) = FIELD(9)
                    SATKEYU( SATCOUNT ) = CSRCALL

                CASE( 'SUNDAY' )
                    IF ( .NOT.DAYFLAG(7) ) CYCLE

                    SUNCOUNT = SUNCOUNT + 1
                    SUNIDU ( SUNCOUNT ) = FIELD(9)
                    SUNKEYU( SUNCOUNT ) = CSRCALL

                CASE( 'WEEKEND' )
                    IF ( .NOT.WKEFLAG ) CYCLE

                    IF ( DAYFLAG( 6 ) ) THEN
                        SATCOUNT = SATCOUNT + 1
                        SATIDU ( SATCOUNT ) = FIELD(9)
                        SATKEYU( SATCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 7 ) ) THEN
                        SUNCOUNT = SUNCOUNT + 1
                        SUNIDU ( SUNCOUNT ) = FIELD(9)
                        SUNKEYU( SUNCOUNT ) = CSRCALL
                    END IF

                CASE( 'WEEKDAY' )

                    IF ( .NOT. WKDFLAG ) CYCLE

                    IF ( DAYFLAG( 1 ) ) THEN
                        MONCOUNT = MONCOUNT + 1
                        MONIDU ( MONCOUNT ) = FIELD(9)
                        MONKEYU( MONCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 2 ) ) THEN
                        TUECOUNT = TUECOUNT + 1
                        TUEIDU ( TUECOUNT ) = FIELD(9)
                        TUEKEYU( TUECOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 3 ) ) THEN
                        WEDCOUNT = WEDCOUNT + 1
                        WEDIDU ( WEDCOUNT ) = FIELD(9)
                        WEDKEYU( WEDCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 4 ) ) THEN
                        THUCOUNT = THUCOUNT + 1
                        THUIDU ( THUCOUNT ) = FIELD(9)
                        THUKEYU( THUCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 5 ) ) THEN
                        FRICOUNT = FRICOUNT + 1
                        FRIIDU ( FRICOUNT ) = FIELD(9)
                        FRIKEYU( FRICOUNT ) = CSRCALL
                    END IF

                CASE( 'ALLDAY' )

                    IF ( DAYFLAG( 1 ) ) THEN
                        MONCOUNT = MONCOUNT + 1
                        MONIDU ( MONCOUNT ) = FIELD(9)
                        MONKEYU( MONCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 2 ) ) THEN
                        TUECOUNT = TUECOUNT + 1
                        TUEIDU ( TUECOUNT ) = FIELD(9)
                        TUEKEYU( TUECOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 3 ) ) THEN
                        WEDCOUNT = WEDCOUNT + 1
                        WEDIDU ( WEDCOUNT ) = FIELD(9)
                        WEDKEYU( WEDCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 4 ) ) THEN
                        THUCOUNT = THUCOUNT + 1
                        THUIDU ( THUCOUNT ) = FIELD(9)
                        THUKEYU( THUCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 5 ) ) THEN
                        FRICOUNT = FRICOUNT + 1
                        FRIIDU ( FRICOUNT ) = FIELD(9)
                        FRIKEYU( FRICOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 6 ) ) THEN
                        SATCOUNT = SATCOUNT + 1
                        SATIDU ( SATCOUNT ) = FIELD(9)
                        SATKEYU( SATCOUNT ) = CSRCALL
                    END IF

                    IF ( DAYFLAG( 7 ) ) THEN
                        SUNCOUNT = SUNCOUNT + 1
                        SUNIDU ( SUNCOUNT ) = FIELD(9)
                        SUNKEYU( SUNCOUNT ) = CSRCALL
                    END IF

                CASE( 'HOURLY' )

                    METCOUNT = METCOUNT + 1
                    METIDU ( METCOUNT ) = FIELD(9)
                    METKEYU( METCOUNT ) = CSRCALL
                    METREFFLAG( JSPC )  = .TRUE.

                    IF ( LASTMET .EQ. MISSMET ) THEN
                        METPROTYPE = FIELD( 8 )
                        LASTMET    = METPROTYPE
                    ELSE IF ( LASTMET .NE. METPROTYPE ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG, '( A,I10 )' ) 'ERROR:  inconsistent MET PROF TYPE at line', L
                        CALL M3MESG( MESG )
                    END IF

                CASE DEFAULT

                    WRITE( MESG, '( 3A,I10 )' ) 'ERROR:  unknown XREFTYPE="', TRIM( FIELD(8) ), '" at line', L
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    W     = W + 1
                    IF ( W .GT. MXWARN ) EXIT

            END SELECT

        END DO      !  end loop reading and processing TREF file

        IF ( W .GT. MXWARN ) THEN
            CALL M3MESG( 'Maximum number of errors exceeded' )
        END IF


        DAYCOUNT = MONCOUNT + TUECOUNT + WEDCOUNT + THUCOUNT +
     &             FRICOUNT + SATCOUNT + SUNCOUNT

        JJ = MTHCOUNT + WEKCOUNT + DOMCOUNT + DAYCOUNT + METCOUNT
        IF ( JJ .LE. 0 ) THEN
            CALL M3EXIT( PNAME, 0,0, 'ERROR:  No TXREFs matching inventory found', 2 )
        END IF


        IF ( EFLAG ) THEN
            CALL M3EXIT( PNAME,0,0, 'ERROR:  reading XREF file', 2 )
        END IF


C.........   Now sort all these tables and check for duplicates:

        IF ( .NOT. SORTREF( 'MONTH-of-YEAR ', MTHCOUNT, MTHCOUNT, MTHIDU, MTHKEYU, MTHIDS, MTHKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing MONTHLY XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'DAY-of-WEEK   ', WEKCOUNT, WEKCOUNT, WEKIDU, WEKKEYU, WEKIDS, WEKKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing WEEKLY XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'DAY-of-MONTH  ', DOMCOUNT, DOMCOUNT, DOMIDU, DOMKEYU, DOMIDS, DOMKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing DAILY XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'MONDAY        ', MONCOUNT, DAYCOUNT, MONIDU, MONKEYU, MONIDS, MONKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing MONDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'TUESDAY       ', TUECOUNT, TUECOUNT, TUEIDU, TUEKEYU, TUEIDS, TUEKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing TUESDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'WEDNESDAY     ', WEDCOUNT, WEDCOUNT, WEDIDU, WEDKEYU, WEDIDS, WEDKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing WEDNESDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'THURSDAY      ', THUCOUNT, THUCOUNT, THUIDU, THUKEYU, THUIDS, THUKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing THURSDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'FRIDAY        ', FRICOUNT, FRICOUNT, FRIIDU, FRIKEYU, FRIIDS, FRIKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing FRIDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'SATURDAY      ', SATCOUNT, SATCOUNT, SATIDU, SATKEYU, SATIDS, SATKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing SATURDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'SUNDAY        ', SUNCOUNT, SUNCOUNT, SUNIDU, SUNKEYU, SUNIDS, SUNKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing SUNDAY DIURNAL XREFs' )
        END IF

        IF ( .NOT. SORTREF( 'MET-BASED     ', METCOUNT, METCOUNT, METIDU, METKEYU, METIDS, METKEYS ) ) THEN
            EFLAG = .TRUE.
            CALL M3MESG( 'ERROR:  processing MET-BASED HOURLY XREFs' )
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( PNAME,0,0, 'ERROR:  sorting/processing XREFs', 2 )
        END IF

        DEALLOCATE( MTHIDU,
     &              WEKIDU,
     &              DOMIDU,
     &              MONIDU,
     &              TUEIDU,
     &              WEDIDU,
     &              THUIDU,
     &              FRIIDU,
     &              SATIDU,
     &              SUNIDU,
     &              METIDU,
     &             MTHKEYU,
     &             WEKKEYU,
     &             DOMKEYU,
     &             MONKEYU,
     &             TUEKEYU,
     &             WEDKEYU,
     &             THUKEYU,
     &             FRIKEYU,
     &             SATKEYU,
     &             SUNKEYU,
     &             METKEYU )

C.........   Read in the relevant profile-tables:

        IF ( MTHCOUNT .GT. 0 ) THEN         !  month-of-year

            MTHNAME = CATEGORY(1:1) //  'TPRO_MONTHLY'
            NMON  = CSVPROF( MTHNAME, 12,
     &                       MTHIDP, MONFAC,
     &                       MTHCOUNT, MTHIDS  )

        ELSE

            NMON  = NULLPROF( 12, MTHIDP, MONFAC )

        END IF      !  if mthcount > 0


        IF ( WEKCOUNT .GT. 0 ) THEN         !  day-of-week

            WEKNAME = CATEGORY(1:1) //  'TPRO_WEEKLY'
            NWEK  = CSVPROF( WEKNAME, 7,
     &                       WEKIDP, WEKFAC,
     &                       WEKCOUNT, WEKIDS  )

        ELSE

            NWEK  = NULLPROF( 7, WEKIDP, WEKFAC )

        END IF      !  if wekcount > 0


C.........  hour-of-day:  all these use TPROF_HOURLY:

        IF ( DAYCOUNT .GT. 0 ) THEN

            !....  Accumulate all these IDs into MONIDS(:)

            N = MONCOUNT
            DO I = 1, TUECOUNT
                N = N + 1
                MONIDS( N ) = TUEIDS( I )
            END DO
            DO I = 1, WEDCOUNT
                N = N + 1
                MONIDS( N ) = WEDIDS( I )
            END DO
            DO I = 1, THUCOUNT
                N = N + 1
                MONIDS( N ) = THUIDS( I )
            END DO
            DO I = 1, FRICOUNT
                N = N + 1
                MONIDS( N ) = FRIIDS( I )
            END DO
            DO I = 1, SATCOUNT
                N = N + 1
                MONIDS( N ) = SATIDS( I )
            END DO
            DO I = 1, SUNCOUNT
                N = N + 1
                MONIDS( N ) = SUNIDS( I )
            END DO

            DIUNAME = CATEGORY(1:1) //  'TPRO_HOURLY'
            NHRL  = CSVPROF( DIUNAME, 24,
     &                       HRLIDP, HRLFAC,
     &                       DAYCOUNT, MONIDS  )

        ELSE

            NHRL  = NULLPROF( 24, HRLIDP, HRLFAC )

        END IF      !  if DAYCOUNT > 0


        IF ( DOMCOUNT .GT. 0 ) THEN         !  day-of-month

            DOMNAME = CATEGORY(1:1) //  'TPRO_DAILY'
            NDOM  = CSVDOMP( DOMNAME, SDATE, EDATE,
     &                       DOMIDP, DOMFAC,
     &                       DOMCOUNT, DOMIDS  )

        ELSE

            NDOM = NULLDOMP( DOMIDP, DOMFAC )

        END IF      !  if domcount > 0


        IF ( METCOUNT .GT. 0 ) THEN       !  met based

C.............  Determine which hourly profiles to apply
            MESG = 'Specifies the basis of hourly profiles: [YEAR|MONTH]'
            CALL ENVSTR( 'HOURLY_TPROF_BASE', MESG, ' ',HOUR_TPROF, IOS )
            CALL UPCASE( HOUR_TPROF )

            IF( .NOT. ( HOUR_TPROF=='MONTH' .OR. HOUR_TPROF=='YEAR' .OR.
     &                  HOUR_TPROF=='DAY' ) ) THEN
                MESG = 'ERROR: MUST define the basis of hourly profiles '//
     &                 'for a correct hourly conversion.'
     &                 //CRLF()//BLANK10//':: Define HOURLY_TPROF_BASE to '//
     &                 '[YEAR, MONTH, or DAY]'
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

            METFLAG = .TRUE.
            METNAME = CATEGORY(1:1) //  'TPRO_HOURLY_NCF'

            IF ( .NOT.OPEN3( METNAME, FSREAD3, PNAME ) ) THEN
                MESG = 'Could not open "' // TRIM( METNAME ) // '"'
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            ELSE  IF( .NOT. DESC3( METNAME ) ) THEN
                MESG = 'Could not DESC3(' // TRIM( METNAME ) // ')'
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            ELSE

                NMETPROF = NROWS3D
                ALLOCATE( COUNTIES( NMETPROF ), STAT = IOS )
                CALL CHECKMEM( IOS, 'COUNTIES', 'PROCTPRO' )

                IF ( .NOT.READ3( METNAME, 'COUNTIES', 1, SDATE, STIME, COUNTIES ) ) THEN
                    MESG = 'Could not READ3("TPRO_HOUR","COUNTIES",...)'
                    CALL M3EXIT( PNAME,SDATE,STIME, MESG, 2 )
                END IF

                DO V = 1, NIPPA

                    IF ( .NOT.METREFFLAG( V ) )  CYCLE

                    WRITE( CPOS, '(I5)' ) V
                    CALL PADZERO( CPOS )

                    DO S = 1, NSRC

                        CFIP = CIFIP( S )
                        TSCC = CSCC( S )
                        CALL PADZERO( CFIP )
                        CALL PADZERO( TSCC )
                        CALL BLDCSRC( CFIP, TSCC, BLANK,
     &                                BLANK, BLANK, BLANK,
     &                                BLANK, CPOS, CSRCALL )

C.........................  If this source in XREF and this FIP in METFIPS
C.........................  use subscript into (unsorted) COUNTIES:

                        I = FINDC( CSRCALL, METCOUNT, METKEYS )         !  index in sorted XREF, or 0
                        IF ( I .GT. 0 )  THEN
                            K = FIND1( STR2INT( CFIP ), NMETPROF, COUNTIES )  !  index into sorted list,   or 0
                            METPROF( S,V ) = K                          !  index into unsorted list, or 0
                        END IF

                    END DO      !  end loop on sources S

                END DO          !  end loop on pollutants V

                DEALLOCATE( COUNTIES, METKEYS )

            END IF              !  if not open3(); else if not desc3(); else...

        ELSE

            NMETPROF = 0
            METNAME  = CMISS
            METFLAG  = .FALSE.

        END IF      ! if met xrefs and profiles used


C.........  Map cross references into profile-subscripts:
C.........  Use zero-subscript as a sentinel-value for "none exist"

        ALLOCATE( MTHPDEX( 0:MTHCOUNT ),
     &            WEKPDEX( 0:WEKCOUNT ),
     &            DOMPDEX( 0:DOMCOUNT ),
     &            MONPDEX( 0:MONCOUNT ),
     &            TUEPDEX( 0:TUECOUNT ),
     &            WEDPDEX( 0:WEDCOUNT ),
     &            THUPDEX( 0:THUCOUNT ),
     &            FRIPDEX( 0:FRICOUNT ),
     &            SATPDEX( 0:SATCOUNT ),
     &            SUNPDEX( 0:SUNCOUNT ), STAT = IOS )
        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG, '( A, I10 )' )
     &           'ERROR:  prof-index allocation failure.  STAT=', IOS
            CALL M3EXIT( PNAME, 0,0, MESG, 2 )
        END IF

        MTHPDEX( 0 ) = 0
        WEKPDEX( 0 ) = 0
        DOMPDEX( 0 ) = 0
        MONPDEX( 0 ) = 0
        TUEPDEX( 0 ) = 0
        WEDPDEX( 0 ) = 0
        THUPDEX( 0 ) = 0
        FRIPDEX( 0 ) = 0
        SATPDEX( 0 ) = 0
        SUNPDEX( 0 ) = 0


C.........  NOTE:  must search N+1 elements of zero-based arrays *IDP(0:N)
C.........         and then adjust by (-1) to skip element 0

        DO N = 1, MTHCOUNT
            MTHPDEX( N ) = FINDC( MTHIDS( N ), NMON+1, MTHIDP ) - 1
            IF ( MTHPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for month-of-year XREF profile-ID '// MTHIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, WEKCOUNT
            WEKPDEX( N ) = FINDC( WEKIDS( N ), NWEK+1, WEKIDP ) - 1
            IF ( WEKPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for day-of-week XREF profile-ID '// WEKIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, DOMCOUNT
            DOMPDEX( N ) = FINDC( DOMIDS( N ), NDOM+1, DOMIDP ) - 1
            IF ( DOMPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for day-of-month XREF profile-ID '// DOMIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, MONCOUNT
            MONPDEX( N ) = FINDC( MONIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( MONPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Monday hour-of-day XREF profile-ID '// MONIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, TUECOUNT
            TUEPDEX( N ) = FINDC( TUEIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( TUEPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Tuesday hour-of-day XREF profile-ID '// TUEIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, WEDCOUNT
            WEDPDEX( N ) = FINDC( WEDIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( WEDPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Wednesday hour-of-day XREF profile-ID '// MTHIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, THUCOUNT
            THUPDEX( N ) = FINDC( THUIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( THUPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Thursday hour-of-day XREF profile-ID '// THUIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, FRICOUNT
            FRIPDEX( N ) = FINDC( FRIIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( FRIPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Friday hour-of-day XREF profile-ID '// FRIIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, SATCOUNT
            SATPDEX( N ) = FINDC( SATIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( SATPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Saturday hour-of-day XREF profile-ID '// SATIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        DO N = 1, SUNCOUNT
            SUNPDEX( N ) = FINDC( SUNIDS( N ), NHRL+1, HRLIDP ) - 1
            IF ( SUNPDEX( N ) .LE. 0 ) THEN
                MESG = 'No profile for Sunday hour-of-day XREF profile-ID '// SUNIDS( N )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        IF ( EFLAG ) THEN
            CALL M3EXIT( PNAME, 0,0, 'Unable to resolve profile(s) for XREF(s)', 2 )
        END IF

C.........  Normalize profiles:

        CALL NORMTPRO()


C******************  INTERNAL SUBPROGRAMS  *****************************


      CONTAINS


        SUBROUTINE TIMEFLAGS()

C Set active-month/active day-of-week flags for the time period SDATE:EDATE.
C Compute SDATE, EDATE for this set of episodes.

C           INTEGER, EXTERNAL :: JSTEP3, WKDAY

            INTEGER, PARAMETER :: DAYSTEP = 240000

            INTEGER     I, J, IT, N, D, JDATE, JTIME

            INTEGER     MAXZONE

C.................  begin body  .....................

            SDATE = ITDATE( 1 )
            EDATE = ITDATE( 1 )

C.............  Find time zone min, max

            MAXZONE = MAXVAL( TZONES )

C.............  Initialize flag-arrays:

            MONFLAG = .FALSE.
            DAYFLAG = .TRUE.

            DO IT = 1, SIZE( ITDATE )        !  loop over episodes for this run

C.................  because of local-time corrections may need one day on each side
C                   of run-start, run-end:

                JDATE = ITDATE( IT )
                JTIME = 0
                CALL NEXTIME( JDATE, JTIME, -MAXZONE*10000 )    !  earliest local starting date for this episode
                STIME = JTIME
                N = RUNLEN( IT ) / 10000        ! hours
                N = 1 + (N+23) / 24             ! days [ (N+23)/24 "rounds up" the days for runlen]

                DO D = 0, N

                    CALL DAYMON( JDATE, I, J )
                    MONFLAG( I ) = .TRUE.

                    IF ( SECSDIFF( SDATE,0, JDATE,0 ) .LT. 0 )  SDATE = JDATE
                    CALL NEXTIME( JDATE, JTIME, DAYSTEP )
                    IF ( SECSDIFF( EDATE,0, JDATE,0 ) .GT. 0 )  EDATE = JDATE

                END DO

            END DO                          !  end loop over episodes M

            WKDFLAG =  (    DAYFLAG(1) .OR. DAYFLAG(2) .OR. DAYFLAG(3)
     &                 .OR. DAYFLAG(4) .OR. DAYFLAG(5) )

            WKEFLAG = ( DAYFLAG(6) .OR. DAYFLAG(7) )

            RETURN

        END SUBROUTINE TIMEFLAGS


C-----------------------------------------------------------------------------------
C   Generate "null" profile TFAC( NFIELDS, 0:0 )
C-----------------------------------------------------------------------------------

        INTEGER FUNCTION  NULLPROF( NFIELDS, IDSTR, TFAC  )

C.............  Arguments:

            INTEGER      ,              INTENT(IN   ) :: NFIELDS
            CHARACTER(16), ALLOCATABLE, INTENT(  OUT) :: IDSTR( : )         !  profile IDs
            REAL         , ALLOCATABLE, INTENT(  OUT) :: TFAC ( :,: )       !  factors in profile

            CHARACTER(24), PARAMETER :: PNAME = 'PROCTPRO/NULLPROF'

C.............  Local variables:

            INTEGER     ISTAT

            !.............  body of function NULLPROF()  .....................

            ALLOCATE(  IDSTR( 0:0 ),
     &          TFAC( NFIELDS,0:0 ), STAT = ISTAT )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' )
     &            'ERROR:  Allocation failure:  STAT=', ISTAT
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            TFAC( 1:NFIELDS,0 ) = 1.0
            IDSTR( 0 )        =  '0'
            NULLPROF          =   0

            RETURN

        END  FUNCTION  NULLPROF


C-----------------------------------------------------------------------------------


        INTEGER FUNCTION  CSVPROF( FNAME, NFIELDS, IDSTR, TFAC,
     &                             IDCNT, IDLIST )

C  Open and count the CSV-profile file FNAME.
C  Allocate  both arguments and local arrays.
C  Read FNAME, and sort data onto output arguments,
C  filtering out data only for ID's in IDLIST
C
C  Lines must be ASCII CSV of the form
C      <character-string ID>, TFAC(1), ..., TFAC(NFIELDS) [comment...]
C------------------------------------------------------------------------------

C.............  Arguments:

            CHARACTER(*) ,              INTENT(IN   ) :: FNAME
            INTEGER      ,              INTENT(IN   ) :: NFIELDS
            CHARACTER(16), ALLOCATABLE, INTENT(  OUT) :: IDSTR( : )         !  profile IDs
            REAL         , ALLOCATABLE, INTENT(  OUT) :: TFAC ( :,: )       !  factors in profile
            INTEGER      ,              INTENT(IN   ) :: IDCNT
            CHARACTER(16),              INTENT(IN   ) :: IDLIST( IDCNT )

C.............  Parameters:

            CHARACTER(24), PARAMETER :: PNAME = 'PROCTPRO/CSVPROF'

C.............  Local variables:

            INTEGER     I, J, K, L, M, N, W, ISTAT
            INTEGER     IDINDX( IDCNT )
            INTEGER     FDEV, NSORT, NLINES, NDATA

            LOGICAL     EFLAG

            CHARACTER(  16) :: THISID, LASTID, AKEY
            CHARACTER(  16) :: IDSORT( IDCNT )
            CHARACTER(  32) :: FIELDS( NFIELDS+2 )
            CHARACTER( 256) :: MESG
            CHARACTER(1024) :: LINE

            CHARACTER(20), ALLOCATABLE :: CKEY( : )
            INTEGER      , ALLOCATABLE :: INDX( : )
            REAL         , ALLOCATABLE :: FACS( :,: )

C.............  body of function CSVPROF()  .....................
C.............  Create sorted-unique list of input IDs:

            DO N = 1, IDCNT
                IDINDX( N ) = N
            END DO

            CALL SORTIC( IDCNT, IDINDX, IDLIST )

            NSORT  = 0
            LASTID = CMISS
            DO I = 1, IDCNT
                THISID = IDLIST( IDINDX( I ) )
                IF ( THISID .NE. LASTID ) THEN
                    NSORT           = NSORT + 1
                    IDSORT( NSORT ) = THISID
                    LASTID          = THISID
                END IF
            END DO


C.............  Open and count FNAME

            FDEV = CSVOPEN( FNAME, NLINES, NDATA )

            ALLOCATE(    CKEY( NDATA ),
     &         FACS( NFIELDS,  NDATA ),
     &                IDSTR( 0:NDATA ),
     &         TFAC( NFIELDS,0:NDATA ), STAT = ISTAT )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' )
     &            'ERROR:  Allocation failure:  STAT=', ISTAT
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            IDSTR( 0 )        =  '0'
            TFAC( 1:NFIELDS,0 ) = 1.0

C.............  Read file:  CKEY and FACS

            EFLAG = .FALSE.

            M = 0
            DO L = 1, NLINES

                READ( FDEV, '( A )', IOSTAT=ISTAT ) LINE
                IF ( ISTAT .NE. 0 ) THEN
                    WRITE( MESG, '( 3 A, I10, 1X, A, I10 )' )
     &                'ERROR:: reading "', TRIM( FNAME ),
     &                '" at line', L,  '"--IOSTAT=', ISTAT
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    W     = W + 1
                    IF ( W .GT. MXWARN ) EXIT
                    CYCLE
                ELSE IF ( BLKORCMT( LINE ) ) THEN
                    CYCLE
                END IF

                CALL PARSLINE( LINE, NFIELDS+2, FIELDS )
                AKEY = FIELDS( 1 )

                IF ( FINDC( AKEY, NSORT, IDSORT ) .LE. 0 ) CYCLE        ! ID does not show up in XREF

                M = M + 1
                CKEY( M ) = AKEY

                DO J = 1, NFIELDS
                    READ( FIELDS( J+1 ), *, IOSTAT=ISTAT ) FACS( J,M )
                    IF ( ISTAT .NE. 0 ) THEN
                        WRITE( MESG, '( 3 A, I10, 1X, A, I10 )' )
     &                  'ERROR:: reading FACTORS from"', TRIM( FNAME ),
     &                 '" at line', L,  '"--IOSTAT=', ISTAT
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                        W     = W + 1
                        IF ( W .GT. MXWARN ) EXIT
                        CYCLE
                    END IF
                END DO

            END DO

99          CLOSE( FDEV )       !  completed input of this file

            IF ( EFLAG ) THEN
                MESG = 'ERROR:  Fatal error(s) reading ' // FNAME
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            CALL  SORTPRO( M, NFIELDS, CKEY, FACS, IDSTR, TFAC )

            WRITE( MESG, '( 4 A, I10, 2X, A, I4, 2X, A )' )
     &        PNAME, ':  file "', TRIM( FNAME ),
     &        '" processed:', M, ' active data-rows', NFIELDS, 'fields'
            CALL M3MESG( MESG )
            DEALLOCATE( CKEY, FACS )

            CSVPROF = M
            RETURN

        END  FUNCTION  CSVPROF


C------------------------------------------------------------------------------
C   Generate "null" profile DMFAC( 31,12,0:0 )
C-----------------------------------------------------------------------------------

        INTEGER FUNCTION  NULLDOMP( IDSTR, DMFAC )

C.............  Arguments:

            CHARACTER(16), ALLOCATABLE, INTENT(  OUT) :: IDSTR( : )         !  profile IDs
            REAL,          ALLOCATABLE, INTENT(  OUT) :: DMFAC( :,:,: )     !  (31,12,0:NDOM==0)

C.............  Local variables:

            INTEGER     ISTAT

            !.............  body of function NULLDOMP()  .....................

            ALLOCATE(       IDSTR( 0:0 ),
     &                DMFAC( 31,12,0:0 ), STAT = ISTAT )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' )
     &            'ERROR:  sorted-profile allocation failure:  STAT=', ISTAT
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            IDSTR( 0 )     =  '0'
            DMFAC( :,:,0 ) = 1.0
            NULLDOMP       =   0

            RETURN

        END  FUNCTION  NULLDOMP


C------------------------------------------------------------------------------


        INTEGER FUNCTION  CSVDOMP( FNAME, SDATE, EDATE, IDSTR, DMFAC,
     &                             IDCNT, IDLIST )

C  Open and count the CSV-day-profile file FNAME.
C  Allocate  both arguments and local arrays.
C  Read FNAME, and sort data onto output arguments,
C  filtering out data only for ID's in IDLIST
C
C  Lines must be ASCII CSV of the form
C      <character-string ID>, MONTH, TFAC(1), ..., TFAC(MON_DAYS(MONTH)) [comment...]
C
C  Output ID's in IDSTR are TRIM(ID)//MM where MM is the 2-digit month
C  for the indicated profile-line
C------------------------------------------------------------------------------

C.............  Arguments:

            CHARACTER(*) ,              INTENT(IN   ) :: FNAME
            INTEGER      ,              INTENT(IN   ) :: SDATE, EDATE
            CHARACTER(16), ALLOCATABLE, INTENT(  OUT) :: IDSTR( : )         !  profile IDs
            REAL         , ALLOCATABLE, INTENT(  OUT) :: DMFAC( :,:,: )     !  (31,12,NDOM)
            INTEGER      ,              INTENT(IN   ) :: IDCNT
            CHARACTER(16),              INTENT(IN   ) :: IDLIST( IDCNT )

C.............  Parameters:

            INTEGER      , PARAMETER :: DAYSTEP = 240000
            CHARACTER(24), PARAMETER :: PNAME   = 'PROCTPRO/CSVDOMP'

C           INTEGER, EXTERNAL :: JSTEP3

C.............  Local variables:

            INTEGER     I, J, K, L, M, N, W, ISTAT, IMON
            INTEGER     IDINDX( IDCNT )
            INTEGER     NSORT, NLINES, NDATA, FDEV
            INTEGER     JDATE, JTIME, NRECS, MON, MDAY
            LOGICAL     LEAPYEAR
            LOGICAL  :: EFLAG = .FALSE.

            CHARACTER(  16) :: THISID, LASTID, AKEY, MISS
            CHARACTER(  16) :: IDSORT( IDCNT )
            CHARACTER(  32) :: FIELDS( 33 )
            CHARACTER( 256) :: MESG
            CHARACTER(1136) :: LINE

            CHARACTER(16), ALLOCATABLE :: CIDU( : )
            REAL         , ALLOCATABLE :: FACS( :,:,: )

C..................  body of function CSVDOMP()  .....................
C.............  Create sorted-unique list of input IDs:

            MISS = CMISS        !  truncated to same length as AKEY

            DO N = 1, IDCNT
                IDINDX( N ) = N
            END DO

            CALL SORTIC( IDCNT, IDINDX, IDLIST )

            NSORT  = 0
            LASTID = CMISS
            DO I = 1, IDCNT
                THISID = IDLIST( IDINDX( I ) )
                IF ( THISID .NE. LASTID ) THEN
                    NSORT           = NSORT + 1
                    IDSORT( NSORT ) = THISID
                    LASTID          = THISID
                END IF
            END DO


C.............  Check if run is impacted by leap year
            LEAPYEAR = .FALSE.
            JTIME = 0
            JDATE = SDATE
            NRECS = JSTEP3( EDATE, JTIME, SDATE, JTIME, DAYSTEP )
            DO I = 1, NRECS
                CALL DAYMON( JDATE, MON, MDAY )
                IF ( MON == 2 .AND. ISLEAP( JDATE ) ) THEN
                    LEAPYEAR = .TRUE.
                    EXIT
                END IF
                CALL NEXTIME( JDATE, JTIME, DAYSTEP )
            END DO

C.............  Open and count FNAME

            FDEV = CSVOPEN( FNAME, NLINES, NDATA )

            ALLOCATE(       CIDU( NSORT ),
     &              FACS( 31,12,0:NSORT ), STAT = ISTAT )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' )
     &            'ERROR:  unsorted-profile allocation failure:  STAT=', ISTAT
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            FACS  = -9999.9      !  array assignments
            FACS( :,:,0 ) = 1.0

C.............  Read file:  CKEY = ID//MONTH, and FACS

            W = 0
            DO L = 1, NLINES

                READ( FDEV, '( A )', END=99, IOSTAT=ISTAT ) LINE
                IF ( ISTAT .NE. 0 ) THEN
                    WRITE( MESG, '( 3 A, I10, 1X, A, I10 )' )
     &                'ERROR: reading "', TRIM( FNAME ),
     &                '" at line', L,  '"--IOSTAT=', ISTAT
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    W     = W + 1
                    IF ( W .GT. MXWARN ) EXIT
                    CYCLE
                ELSE IF ( BLKORCMT( LINE ) ) THEN
                    CYCLE
                END IF

                CALL PARSLINE( LINE, 33, FIELDS )
                AKEY = FIELDS( 1 )

                N    = FINDC( AKEY, NSORT, IDSORT )
                IF ( N .LE. 0 ) CYCLE        ! ID does not show up in XREF

                CIDU( N  ) = AKEY

                READ( FIELDS( 2 ), *, IOSTAT=ISTAT ) IMON
                IF ( ISTAT .NE. 0 ) THEN
                    WRITE( MESG, '( 3 A, I10, 1X, A, I10 )' )
     &                'ERROR: reading MONTH from"', TRIM( FNAME ),
     &                 '" at line', L,  '"--IOSTAT=', ISTAT
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    W     = W + 1
                    IF ( W .GT. MXWARN ) EXIT
                    CYCLE
                ELSE IF ( IMON .LT. 1 .OR. IMON .GT. 12 ) THEN
                    WRITE( MESG, '( 3 A, I10, 1X, A, I10 )' )
     &                'ERROR: invalid MONTH from"', TRIM( FNAME ),
     &                 '" at line', L,  '"--IOSTAT=', ISTAT
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    W     = W + 1
                    IF ( W .GT. MXWARN ) EXIT
                    CYCLE
                ELSE IF ( .NOT.MONFLAG( IMON ) ) THEN
                    CYCLE
                END IF

                IF ( IMON .NE. 2 ) THEN
                    NN = MON_DAYS( IMON )
                ELSE IF ( LEAPYEAR ) THEN
                    NN = 29     ! = MON_DAYS( IMON ) + 1
                ELSE
                    NN = 28     ! = MON_DAYS( IMON )         ! = 28 for normal february
                END IF

                DO J = 1, NN
                    READ( FIELDS( J+2 ), *, IOSTAT=ISTAT ) FACS( J,IMON,N )
                    IF ( ISTAT .NE. 0 ) THEN
                        WRITE( MESG, '( 3 A, I10, 1X, A, I4, 1X, A, I10 )' )
     &                  'ERROR:: reading FACTORS from"', TRIM( FNAME ),
     &                  '" at line', L,  'field', J+2,  '"--IOSTAT=', ISTAT
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                        W     = W + 1
                        IF ( W .GT. MXWARN ) EXIT
                        CYCLE
                    END IF
                END DO

            END DO

99          CLOSE( FDEV )       !  completed input of this file

            IF ( EFLAG ) THEN
                MESG = 'ERROR:  Fatal error(s) reading ' // FNAME
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

C.............  Sort/process IDs: Count actually-occurring IDs:

            M = 0
            DO I = 1, NSORT
                IF ( CIDU( I ) .NE. MISS )  M = M + 1
                IDINDX( I ) = I
            END DO

            CALL SORTIC( NSORT, IDINDX, CIDU )

            ALLOCATE( IDSTR( 0:M ),
     &                DMFAC( 31,12,0:M ), STAT = ISTAT )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' )
     &            'ERROR:  sorted-profile allocation failure:  STAT=', ISTAT
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            IDSTR( 0 )     = '0'
            DMFAC( :,:,0 ) = 1.0

            IDSTR = BLANK
            N      = 0
            LASTID = CMISS
            DO I = 1, NSORT
                K = IDINDX( I )
                IF ( CIDU( K ) .NE. LASTID .AND.
     &               CIDU( K ) .NE. BLANK ) THEN
                    N = N + 1
                    LASTID     = CIDU( K )
                    IDSTR( N ) = CIDU( K )
                    DMFAC( :,:,N ) = FACS( :,:,K )
                END IF
            END DO

            DEALLOCATE( FACS, CIDU )
            CSVDOMP = N

            RETURN

        END  FUNCTION  CSVDOMP


C------------------------------------------------------------------------------


        INTEGER FUNCTION  CSVOPEN( FNAME, NLINES, NDATA )

C  Open, and count lines and data-lines in the CSV-profile file FNAME.
C  Rewind after completion.
C------------------------------------------------------------------------------

C.............  Arguments:

            CHARACTER(*), INTENT(IN   ) :: FNAME        !  logical file name
            INTEGER     , INTENT(  OUT) :: NLINES       !  total number of lines
            INTEGER     , INTENT(  OUT) :: NDATA        !  number of non-comment lines

C.............  Externals and Parameters:

C           INTEGER, EXTERNAL  :: GETEFILE

            CHARACTER(24), PARAMETER :: PNAME = 'PROCTPRO/CSVOPEN'

C.............  Local variables:

            INTEGER     F, I, J, K, L, M, N, W
            INTEGER     FDEV

            LOGICAL     EFLAG

            CHARACTER(  1) :: CBUF
            CHARACTER(256) :: LINE, MESG

C..................  body of function CSVOPEN()  .....................

            CALL NAMEVAL(  FNAME, LINE )
            CALL UPCASE( LINE )
            IF ( LINE .EQ. 'NONE' ) THEN
                MESG = 'File turned off: "' // TRIM( FNAME ) //'"=="NONE"'
                CALL M3MESG( MESG )
                NLINES  = 0
                NDATA   = 0
                CSVOPEN = 0
                RETURN
            END IF

            FDEV = GETEFILE( FNAME, .TRUE., .TRUE., PNAME )
            IF ( FDEV .LT. 0 ) THEN
                MESG = 'ERROR:  Could not open "' // FNAME
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

C.................  Count entries in FDEV

            EFLAG = .FALSE.
            L     = 0
            W     = 0
            M     = 0
            DO

                L = L + 1
                READ( FDEV, '( A )', END=11, IOSTAT=ISTAT ) CBUF
                IF ( ISTAT .NE. 0 ) THEN
                    WRITE( MESG, '( 3 A, I10, 1X, A, 1X, I10 )' )
     &                'ERROR:: counting "', TRIM( FNAME ),
     &                'at line', L,  '"--IOSTAT=', ISTAT
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    W     = W + 1
                    IF ( W .GT. MXWARN ) EXIT
                    CYCLE
                ELSE IF ( BLKORCMT( LINE ) ) THEN
                    CYCLE
                ELSE
                    M = M + 1
                END IF

            END DO

11          CONTINUE
            IF ( EFLAG ) THEN
                MESG = 'ERROR:  Fatal error(s) counting ' // FNAME
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            NLINES  = L-1
            NDATA   = M
            CSVOPEN = FDEV

            REWIND( FDEV )

            RETURN

        END  FUNCTION  CSVOPEN


C------------------------------------------------------------------------------
C  Sort profile data into output data-structures
C------------------------------------------------------------------------------


        SUBROUTINE  SORTPRO( NROWS, NFIELDS, CKEY, FACS, CSRT, FSRT )

C.............  Arguments:

            INTEGER     , INTENT(IN   ) :: NROWS, NFIELDS
            CHARACTER(*), INTENT(IN   ) ::  CKEY( NROWS )           !  input IDs
            REAL        , INTENT(IN   ) ::  FACS( NFIELDS,  NROWS ) !  input coeffs
            CHARACTER(*), INTENT(  OUT) ::  CSRT( 0:NROWS )         !  sorted IDs
            REAL        , INTENT(INOUT) ::  FSRT( NFIELDS,0:NROWS ) !  output coeffs

C.............  Local variables:

            INTEGER     I, J, K, L, M, N, W
            INTEGER     INDX( NROWS )

C..................  body of function SORTPRO()  .....................

            DO L = 1, NROWS
                INDX( L ) = L
            END DO

            CALL SORTIC( NROWS, INDX, CKEY )

            DO L = 1, NROWS
                K = INDX( L )
                CSRT (   L ) = CKEY( K )
                FSRT ( :,L ) = FACS( :,K )
            END DO

        END  SUBROUTINE  SORTPRO


C------------------------------------------------------------------------------
C  Sort XREF data into output data-structures
C------------------------------------------------------------------------------

        LOGICAL FUNCTION  SORTREF( TYPE, COUNT, NDIM, IDLIST, KEYLIST, IDSORT, KEYSORT )

            CHARACTER(*),                    INTENT(IN   ) :: TYPE
            INTEGER,                         INTENT(IN   ) :: COUNT
            INTEGER,                         INTENT(IN   ) :: NDIM      ! usually same as COUNT
            CHARACTER(16),                   INTENT(IN   ) ::  IDLIST( COUNT )
            CHARACTER(ALLLEN3),              INTENT(IN   ) :: KEYLIST( COUNT )
            CHARACTER(16),      ALLOCATABLE, INTENT(  OUT) ::  IDSORT( : )
            CHARACTER(ALLLEN3), ALLOCATABLE, INTENT(  OUT) :: KEYSORT( : )

            INTEGER             INDX( COUNT )
            INTEGER             I, K, ISTAT
            LOGICAL             EFLAG
            CHARACTER(ALLLEN3)  LASTKEY
            CHARACTER(256)      MESG

C..............  Always need to allocate tables, even if size is zero:

            WRITE( MESG, '( 3A, I9 )' )
     &            'PROCTPRO:  TREF type "', TYPE,
     &            '" active reference count=', COUNT
            CALL M3MESG( MESG )

            ALLOCATE( IDSORT( NDIM ), KEYSORT( NDIM ),STAT=ISTAT )
            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A,I10 )' )
     &                'ERROR:  allocating XREF tables.  STAT=', ISTAT
                CALL M3MESG( MESG )
                SORTREF = .FALSE.
                RETURN
            END IF

            IF ( COUNT .EQ. 0 ) THEN
                SORTREF = .TRUE.
                RETURN
            END IF

            DO I = 1, COUNT
                INDX( I ) = I
            END DO

            CALL SORTIC( COUNT, INDX, KEYLIST )

            LASTKEY = CMISS
            EFLAG   = .FALSE.
            DO I = 1, COUNT
                K = INDX( I )
                IDSORT ( I ) = IDLIST( K )
                KEYSORT( I ) = KEYLIST( K )
                IF ( LASTKEY .EQ. KEYLIST( K ) ) THEN
                    EFLAG = .TRUE.
                ELSE
                    LASTKEY = KEYLIST( K )
                END IF
            END DO

            IF ( EFLAG ) THEN
                CALL M3MESG( 'ERROR:  duplicate keys in XREF' )
            END IF

            SORTREF = ( .NOT.EFLAG )
            RETURN

        END  FUNCTION  SORTREF


C------------------------------------------------------------------------------


        LOGICAL FUNCTION  ISLEAP( JDATE )

C............  begin body  ........................

            INTEGER, INTENT(IN   ) :: JDATE

            INTEGER     YEAR

            YEAR = JDATE / 1000

            IF ( MOD( YEAR,4 ) .NE. 0 ) THEN            !  2001,2002, etc... normal years
                ISLEAP = .FALSE.
            ELSE IF ( MOD( YEAR,100 ) .NE. 0 ) THEN     !  "normal" leap-years
                ISLEAP = .TRUE.
            ELSE IF ( MOD( YEAR,400 ) .NE. 0 ) THEN     !  1800,1900,2100 ... century nonleap years by Gregory's rule
                ISLEAP = .FALSE.
            ELSE                                        !  1600,2000,2400 ... century leap years by Gregory's rule
                ISLEAP = .TRUE.
            END IF

             RETURN

        END FUNCTION ISLEAP


      END SUBROUTINE PROCTPRO
