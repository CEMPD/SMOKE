
        SUBROUTINE RDMEDSPD( FDEV, TZONE, TYPNAM, LASTFLAG )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads the day- or hour-specific emissions in
C      MEDS format. It appends the records to the global storage from the MODDAYHR
C
C  PRECONDITIONS REQUIRED:
C      Must complete processing annual/avg inventory
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created by B.H. Baek on 12/2013
C
C      01/2026 by tranhuy:  Use M3UTILIO
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
C.........  MODULES for public variables
        USE M3UTILIO

        USE MODSOURC, ONLY: CIFIP, CSOURC, NMEDGAI, COABDST

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, ITFACA, SCASIDX,
     &                      ITKEEPA, SORTCAS, SCASIDX, NUNIQCAS,
     &                      UCASNPOL, UNIQCAS, UCASIDX, UCASNKEP

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NIPPA, NSRC, EANAM, NCHARS

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, LPDSRC, IDXSRC, PDEMOUT

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, USEDAYLT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
C        CHARACTER(2) CRLF
C        INTEGER      ENVINT
C        LOGICAL      ENVYN, CHKINT
C        INTEGER      FIND1
C        INTEGER      FINDC
C        INTEGER      INDEX1
C        INTEGER      JULIAN
C        INTEGER      SECSDIFF
C        INTEGER      STR2INT
C        REAL         STR2REAL
C        REAL         YR2DAY
C        INTEGER      YEAR4
C        INTEGER      GETTZONE
C        LOGICAL      ISDSTIME

C        EXTERNAL     CRLF, ENVINT, ENVYN, FIND1, FINDC, INDEX1, JULIAN, 
C     &               SECSDIFF, STR2INT, STR2REAL, YEAR4, YR2DAY, CHKINT,
C     &               GETTZONE, ISDSTIME
        INTEGER, EXTERNAL :: GETTZONE


C.........  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN)  :: FDEV           ! file unit no.
        INTEGER,      INTENT (IN)  :: TZONE          ! output time zone
        CHARACTER(*), INTENT (IN)  :: TYPNAM         ! 'day' or 'hour'
        LOGICAL,      INTENT (IN)  :: LASTFLAG       ! process last inv input file

C...........   Local list of bad sources to prevent duplicate writing of error
C              messages
        CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )

C...........   Local parameters
        INTEGER, PARAMETER :: MXSEG = 60     ! max no of segments

C...........   Local segment arrays
        INTEGER            :: SPSTAT( MXSPDAT ) = 0    ! true: special data variable used
        INTEGER            :: EASTAT( NIPPA )    ! > 0 poll present in data
        INTEGER            :: EAIDX ( NIPPA )    ! index of poll present in data

C...........   Local list of FIPS start/end positions to facilitate
C              faster lookups
        INTEGER, ALLOCATABLE, SAVE :: STARTSRC( : )
        INTEGER, ALLOCATABLE, SAVE :: ENDSRC( : )

C...........   Local list of arrays for warning handling
        LOGICAL, ALLOCATABLE, SAVE :: WARNKEEP( : ) ! true: write warning for Keep = N
        LOGICAL, ALLOCATABLE, SAVE :: WARNMULT( : ) ! true: write warning for Multiple pollutants from a single pollutant in Inventory Table

C...........   Temporary read arrays
        REAL          :: TDAT = 0.0       ! temporary data values

C...........   Other local variables
        INTEGER          D, H, HS, I, J, L, L1, L2, NP, V, S, T    ! counters and indices
        INTEGER          SC, EC        ! start/end col no for inv poll
        INTEGER          ES, NS, SS    ! end src, tmp no. src, start sourc

        INTEGER          CIDX             ! tmp data index
        INTEGER          COD              ! data index
        INTEGER          DAY              ! tmp day of month
        INTEGER, SAVE :: ICC = 0          ! tmp country code from header
        INTEGER          IOS              ! i/o status
        INTEGER          IREC             ! record counter
        INTEGER       :: INVFMT = 12      ! MEDS fomrat index
        INTEGER          JDATE, TDATE     ! tmp Julian date
        INTEGER          JTIME, TTIME     ! tmp HHMMSS time
        INTEGER       :: NHOUR = 1        ! processing hr (day=24, hr=1)
        INTEGER, SAVE :: SDATE = 0        ! tmp starting Julian date
        INTEGER, SAVE :: STIME = 0        ! tmp starting HHMMSS time
        INTEGER       :: OUTSTEP = 10000  ! output time step
        INTEGER, SAVE :: MXWARN       	  ! max no. warnings
        INTEGER, SAVE :: NBADSRC = 0      ! no. bad sources
        INTEGER, SAVE :: NFIELD = 1       ! number of data fields
        INTEGER, SAVE :: NPOA   = 6       ! no of pol in MEDS
        INTEGER, SAVE :: NSTEPS = 0       ! number of time steps
        INTEGER, SAVE :: NWARN( 5 )       ! warnings counter
        INTEGER          YR4              ! MEDS inv year
        INTEGER          ZONE             ! source time zones
        INTEGER       :: RDEV = 0         !  unit no. for REPINVEN file

        REAL             CONVFAC          ! tmp conversion factor from Inventory Table

        LOGICAL       :: EFLAG  = .FALSE. ! TRUE iff ERROR
        LOGICAL       :: WARNOUT = .FALSE.! true: then output warnings
        LOGICAL, SAVE :: DAYFLAG = .FALSE.! true: processing day-specific inv
        LOGICAL, SAVE :: INITIAL  = .TRUE.! true: initialize sdate/stime for one time
        LOGICAL, SAVE :: ONETIME  = .TRUE.! true: one time routine called
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called
        LOGICAL, SAVE :: TFLAG  = .FALSE. ! true: use SCCs for matching with inv

        CHARACTER( 3 ) :: STA = '006'     ! state code for CA (=006)
        CHARACTER(100) :: BUFFER = ' '    ! src description buffer 
        CHARACTER(1920):: LINE   = ' '    ! line buffer 
        CHARACTER(512) :: MESG   = ' '    ! message buffer
        CHARACTER(CHRLEN3) ::GAI = ' '    ! GAI lookup code
 
        CHARACTER( 3 )    ARBN, CNTY
        CHARACTER(FIPLEN3) CFIP      ! tmp co/st/cy code
        CHARACTER(FIPLEN3) LFIP      ! previous st/co FIPS code
        CHARACTER(CASLEN3) CDAT      ! tmp Inventory data (input) name
        CHARACTER(IOVLEN3) CNAM      ! tmp SMOKE name
        CHARACTER(PLTLEN3) FCID      ! tmp facility ID
        CHARACTER(CHRLEN3) SKID      ! tmp stack ID
        CHARACTER(CHRLEN3) DVID      ! tmp device ID
        CHARACTER(CHRLEN3) PRID      ! tmp process ID
        CHARACTER(SCCLEN3) TSCC      ! tmp source category code
        CHARACTER(ALLLEN3) CSRC      ! tmp source string
        CHARACTER(NAMLEN3) OUTIDX    ! name for integer source index
        CHARACTER(NAMLEN3) ONAME     ! output logical filename
 
        CHARACTER(IOVLEN3) POLNAM( 6 )

        CHARACTER(16) :: PROGNAME = 'RDMEDSPD' !  program name

C***********************************************************************
C   begin body of program RDMEDSPD
C.........  First time routine called
        IF( FIRSTIME ) THEN

C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET , ' ', 100, I )

C.............  Allocate memory for bad source storage
            ALLOCATE( BADSRC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

C.............  list of MEDS pol names
            POLNAM( 1 ) = 'CO'
            POLNAM( 2 ) = 'NOX'
            POLNAM( 3 ) = 'SOX'
            POLNAM( 4 ) = 'TOG'
            POLNAM( 5 ) = 'PM'
            POLNAM( 6 ) = 'NH3'

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

C.............  Initialize local variables
            NWARN = 0  ! array
            INVFMT = 12
            SPSTAT = 0
            FIRSTIME = .FALSE.

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

            NIPPA = NPOA   ! same no of poll (6) in MEDS (CO,NOX,SOX,TOG,PM,NH3)

C.............  Determine if file is day- or hour-specific by the length of the
C               lines. Make sure day- and hour-specific data are not in the
C               same file.
C.............  If the file is hourly but the only the daily is to be read, then
C               behave as if it is a daily file.

C.............  Set Julian day from MMDDYY8 SAS format
            IF ( TYPNAM == 'DAY' .OR. TYPNAM == 'day' ) THEN
                DAYFLAG = .TRUE.
                OUTIDX  = 'INDXD'
                IF( LINE( 64:64 ) /= '-' ) THEN
                    MESG = 'ERROR: CAN NOT process Hourly MEDS inventories '//
     &                     'when DAY_SPECIFIC_YN is set to Y'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            ELSE IF( TYPNAM == 'HOUR' .OR. TYPNAM == 'hour' ) THEN
                DAYFLAG = .FALSE.
                OUTIDX  = 'INDXH'
                IF( STR2INT( LINE( 64:65 ) ) < 0 ) THEN
                    MESG = 'ERROR: CAN NOT process Daily MEDS inventories '//
     &                     'when HOUR_SPECIFIC_YN is set to Y'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            ELSE 
                MESG = 'INTERNAL ERROR: Do not know type ' // TYPNAM
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            YR4 = 2000 + STR2INT( LINE( 59:60 ) )
            JDATE = STR2INT( ADJUSTL( LINE( 61:63 ) ) )
            JDATE = YR4 * 1000 + JDATE

            JTIME = 00000
            IF( .NOT. DAYFLAG ) JTIME = STR2INT( LINE( 64:65 ) ) * 10000

C.............  Search for time zone for current county
            GAI = ADJUSTL( LINE( 71:73 ) )  ! GAI lookup code

            IF( .NOT. ALLOCATED( COABDST ) ) THEN
                 MESG='ERROR: MUST provide GAI_LOOKUP_TABLE file '
     &              //'for pregridded MEDS-formatted inventory'
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = INDEX1( GAI, NMEDGAI, COABDST( :,1 ) )
            IF( I < 1 ) THEN
                ARBN = ADJUSTR( LINE( 68:70 ) )
                CALL PADZERO( ARBN )
                WRITE( CNTY, '(I3.3)' ) STR2INT( LINE( 57:58 ) )
                CFIP = ARBN // STA // CNTY // '000'
            ELSE 
                CFIP = COABDST( I,2 )      ! FIPS code
            END IF

            I = FINDC( CFIP, NCOUNTY, CNTYCOD )

C.............  If time zone name is not found, thenoutput error
            IF( I .LE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find time zone for county: '//
     &                 CFIP // ' from GEOCODE4 file'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Set time zone number
            ZONE = GETTZONE( CFIP )
 
C.............  If daily emissions are not in the output time zone, print 
C               warning
            IF( WARNOUT .AND. .NOT. DAYFLAG .AND. ZONE .NE. TZONE .AND.
     &          NWARN( 1 ) .LE. MXWARN ) THEN
                WRITE( MESG,94010 ) 
     &                'WARNING: Time zone ', ZONE, 'in hourly-specific ' //
     &                'file at line', IREC, CRLF() // BLANK10 //  
     &                'does not match output time zone', TZONE
                CALL M3MESG( MESG )
                NWARN( 1 ) = NWARN( 1 ) + 1
            END IF

C.............  Check if date is in daylight time, if local zone has
C               already been converted, and if this FIPS code is
C               exempt from daylight time or not.
            IF( ISDSTIME( JDATE ) .AND. USEDAYLT( I ) ) THEN
                ZONE = ZONE - 1
            END IF

C.............  Convert date and time to output time zone.
            CALL NEXTIME( JDATE, JTIME, ( ZONE - TZONE ) * 10000 )

C.............  Intialize SDATE and STIME for one time
            IF( INITIAL ) THEN
                SDATE = JDATE
                STIME = JTIME
                INITIAL = .FALSE.
            END IF

C.............  Write emissions for this time step
            IF( SDATE .NE. JDATE .OR. STIME .NE. JTIME ) THEN

                TDATE = SDATE
                TTIME = STIME
                IF( DAYFLAG ) NHOUR = 24

                DO T = 1, NHOUR

C.....................  Open day-specific or hour-specific output file
                    IF( ONETIME ) THEN
                        CALL OPENPDOUT( NSRC, NPOA, TZONE, TDATE, TTIME, OUTSTEP,
     &                              INVFMT, TYPNAM, .FALSE., EAIDX, SPSTAT,
     &                              ONAME, RDEV )
                        ONETIME = .FALSE.
                    END IF

C.....................   Output daily or hour emissions to PDAY output file
                    IF ( .NOT. WRITE3( ONAME, OUTIDX, TDATE, TTIME, IDXSRC ) ) THEN
                         L2   = LEN_TRIM( ONAME )
                         MESG= 'Error writing output file "' // ONAME(1:L2) // '"'
                         CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                    END IF

                    DO V = 1, NPOA 
                        IF ( .NOT. WRITE3( ONAME, POLNAM( V ), TDATE, TTIME, PDEMOUT( :,V ) ) ) THEN
                             L2   = LEN_TRIM( ONAME )
                             MESG= 'Error writing output file "' // ONAME(1:L2) // '"'
                             CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                        END IF
                    END DO

                    IF( DAYFLAG )  CALL NEXTIME( TDATE, TTIME, 10000 )

                END DO

            END IF

C.............  Set key for searching sources
            FCID = ADJUSTl( LINE( 43:51 ) )  ! platn/facility ID (=FCID)
            SKID = ADJUSTL( LINE( 37:39 ) )  ! Column ID (=pointID, SKID)
            DVID = ADJUSTL( LINE( 52:56 ) )  ! stack ID (=DeviceID), DVID)
            PRID = ADJUSTL( LINE( 40:42 ) )  ! Row ID (=ProcessID, PRID)

            TSCC = ADJUSTR( LINE(  9:22 ) )     ! SCC from MEDS format
            CALL PADZERO( TSCC )

C.............  If FIPS code is not the same as last time, then
C               look it up and get indidies
            IF( CFIP .NE. LFIP ) THEN
                J = FINDC( CFIP, NINVIFIP, INVCFIP )
                IF( J .LE. 0 ) THEN
                    MESG = 'INTERNAL ERROR: Could not find FIPS code ' // 
     &                     CFIP // ' in internal list.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                END IF

                SS = STARTSRC( J )
                ES = ENDSRC( J )
                NS = ES - SS + 1
                LFIP = CFIP
            END IF

C.............  Build source characteristics field for searching inventory
            CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,
     &                    TSCC, CHRBLNK3, POLBLNK3, CSRC )
 
C.............  Search for this record in sources
            J = FINDC( CSRC, NS, CSOURC( SS ) )

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
                        WRITE( MESG,94010 ) 'WARNING: Period-specific record at line',
     &                       IREC, ' does not match inventory sources: '//
     &                         CRLF() // BLANK10 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN( 3 ) = NWARN( 3 ) + 1
                    END IF

                END IF

                CYCLE               !  to head of read loop

C.............  Otherwise, update master list of sources in the inventory
            ELSE
                S = SS - 1 + J         ! calculate source number
                LPDSRC( S ) = .TRUE.

            END IF

            SC = 79     ! starting col no for poll inv

C.............  Loop over no of pollutant for the source
            DO NP = 1, NPOA
          
C.............  Check pollutant code and set index I
                CDAT = POLNAM( NP )         ! CO, NOX, SOX, TOG, PM, and NH3 (in order)
                CDAT = ADJUSTL( CDAT ) 
                CALL UPCASE( CDAT ) 

C.................  Look up pollutant name in unique sorted array of
C                   Inventory pollutant names
                CIDX  = FINDC( CDAT, NUNIQCAS, UNIQCAS )

C.................  Check to see if data name is in inventory list
                COD  = INDEX1( CDAT, NIPPA, EANAM )
            
                IF( COD < 1 ) THEN
                    WRITE( MESG,94010 )
     &                   'WARNING: Skipping pollutant "'// TRIM(CDAT)//
     &                   '" at line', IREC, '- not in annual inventory'
                        CALL M3MESG( MESG )
                ELSE IF( CIDX < 1 ) THEN
                    WRITE( MESG,94010 )
     &                   'WARNING: Skipping pollutant "'// TRIM(CDAT)//
     &                   '" at line', IREC, '- not in inventory table'
                        CALL M3MESG( MESG )
                ELSE
                    EASTAT( COD ) = CIDX 
                    EAIDX ( NP  ) = COD

                END IF

C.................  Set conversion factor from Inventory Table. Default is 1.0
                CONVFAC = ITFACA( SCASIDX( UCASIDX( CIDX ) ) )
                
                EC = SC + 9
C.................  Convert kg to tons unit
                TDAT = CONVFAC * STR2REAL( LINE( SC:EC ) ) * 1000.0 * GM2TON

                PDEMOUT( S,NP ) = PDEMOUT( S,NP ) + TDAT

                SC = EC + 1

            END DO

            SDATE = JDATE
            STIME = JTIME

        END DO

299     CONTINUE   ! Exit from read loop

C.........  Write emissions for last date/hour input file
        IF( LASTFLAG ) THEN

            TDATE = SDATE
            TTIME = STIME
            IF( DAYFLAG ) NHOUR = 24

            DO T = 1, NHOUR

C.................  Open day-specific or hour-specific output file
                IF( ONETIME ) THEN
                    CALL OPENPDOUT( NSRC, NPOA, TZONE, TDATE, TTIME, OUTSTEP,
     &                              INVFMT, TYPNAM, .FALSE., EAIDX, SPSTAT,
     &                              ONAME, RDEV )
                    ONETIME = .FALSE.
                END IF

C.................   Output daily or hour emissions to PDAY output file
                IF ( .NOT. WRITE3( ONAME, OUTIDX, TDATE, TTIME, IDXSRC ) ) THEN
                     L2   = LEN_TRIM( ONAME )
                     MESG= 'Error writing output file "' // ONAME(1:L2) // '"'
                     CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                END IF

                DO V = 1, NPOA 
                    IF ( .NOT. WRITE3( ONAME, POLNAM( V ), TDATE, TTIME, PDEMOUT( :,V ) ) ) THEN
                         L2   = LEN_TRIM( ONAME )
                         MESG= 'Error writing output file "' // ONAME(1:L2) // '"'
                         CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                    END IF
                END DO

                IF( DAYFLAG )  CALL NEXTIME( TDATE, TTIME, 10000 )

            END DO

        END IF

C.........  Abort if error found while reading file
        IF( EFLAG ) THEN
            MESG = 'Problem processing day- or hour-specific data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( I6.6 )

        END SUBROUTINE RDMEDSPD
