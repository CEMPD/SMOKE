
        SUBROUTINE WRIDAOUT( DDEV, VDEV, TDEV, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Write IDA output file based on source characteristics in memory and
C      temporary pollutant and activity files.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C...........   This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        INTEGER         FIND1
        INTEGER         GETIFDSC
        INTEGER         INDEX1

        EXTERNAL        FIND1, GETIFDSC, INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER      , INTENT (IN) :: DDEV           ! emissions unit no.
        INTEGER      , INTENT (IN) :: VDEV           ! activity unit no.
        INTEGER      , INTENT (IN) :: TDEV( NIPPA )  ! tmp files unit nos.
        INTEGER      , INTENT(OUT) :: STATUS         ! exit status

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: PLANTLEN = 15 ! IDA fmt plant ID length
        INTEGER, PARAMETER :: POINTLEN = 15 ! IDA fmt point ID length 
        INTEGER, PARAMETER :: STACKLEN = 12 ! IDA fmt stack ID length
        INTEGER, PARAMETER :: BOILRLEN = 5  ! IDA fmt boiler ID length
        INTEGER, PARAMETER :: SEGMTLEN = 2  ! IDA fmt segment ID length
        INTEGER, PARAMETER :: PDESCLEN = 40 ! IDA fmt plant description length
        INTEGER, PARAMETER :: SCCLEN   = 10 ! IDA fmt SCC field length

C.........  Parameters for output IDA formats
        INTEGER    , PARAMETER :: AROLEN ( 7 ) = ( /10,10,11,7,3,6,0/ )
        INTEGER    , PARAMETER :: MBOPLEN( 7 ) = ( /10,10,0,0,0,0,0/ )
        INTEGER    , PARAMETER :: MBOALEN( 7 ) = ( /14,0,0,0,0,0,0/ )
        INTEGER    , PARAMETER :: PTOLEN ( 7 ) = ( /13,13,7,3,10,3,3/ )

        CHARACTER*6, PARAMETER :: AROFMT ( 7 )  = 
     &                        ( / ' F10.4',' F10.4',' F11.4','  F7.2',
     &                            '    I3','  F6.2','      ' / )
        CHARACTER*6, PARAMETER :: MBOPFMT( 7 ) = 
     &                        ( / ' F10.4',' F10.4','      ','      ',
     &                            '      ','      ','      ' / )
        CHARACTER*6, PARAMETER :: MBOAFMT( 7 ) = 
     &                        ( / ' F16.3','      ','      ','      ',
     &                            '      ','      ','      ' / )
        CHARACTER*6, PARAMETER :: PTOFMT ( 7 ) = 
     &                        ( / ' F13.4',' F13.4','  F7.2','    I3',
     &                            ' F10.5','    I3','    I3' / )

C...........   Local allocatable arrays
        REAL        , ALLOCATABLE :: DATARECS( : )  ! pol/act data
        CHARACTER*52, ALLOCATABLE :: COUTRECS( : )    ! formatted pol/act data

C...........   Local fixed arrays
        INTEGER         OUTTYPE( 7 )

C...........   IDA output variables (names same as IDA format description)

        INTEGER         STID, CYID, BEGYR, ENDYR, HOURS, START, NDAY
        INTEGER         WEEKS, SIC

        REAL            STKHGT, STKDIAM, STKTEMP, STKFLOW, STKVEL
        REAL            BOILCAP, WINTHRU, SPRTHRU, SUMTHRU, FALTHRU
        REAL            THRUPUT, MAXRATE, HEATCON, SULFCON, ASHCON
        REAL            NETDC, LATC, LONC

        CHARACTER(LEN=1)        CAPUNITS, OFFSHORE  
        CHARACTER(LEN=LNKLEN3 ) CLNK  
        CHARACTER(LEN=PLANTLEN) PLANTID  
        CHARACTER(LEN=POINTLEN) POINTID
        CHARACTER(LEN=STACKLEN) STACKID
        CHARACTER(LEN=BOILRLEN) BLRID
        CHARACTER(LEN=SEGMTLEN) SEGMENT
        CHARACTER(LEN=ORSLEN3 ) ORISID
        CHARACTER(LEN=PDESCLEN) PLNTDESC
        CHARACTER(LEN=SCCLEN  ) SCC

C...........   Other local variables

        INTEGER         I, J, L, L1, L2, K, N, S  ! counters and indices

        INTEGER         COID     ! tmp country code
        INTEGER         FDEV     ! tmp unit no.
        INTEGER         FIP      ! tmp FIPS state and county code
        INTEGER         IOS      ! i/o status
        INTEGER         LCOID    ! previous country ID in loop
        INTEGER         LYEAR    ! previous year in loop
        INTEGER         NCHAR    ! number of strings returned from PARSCSRC
        INTEGER         YEAR     ! tmp 4-digit year

        REAL            CTOF     ! celius to farenheit
        REAL            M2FT     ! meters to feet

        LOGICAL IDACOLS( 7 )
        DATA    IDACOLS / 6*.TRUE., .FALSE. /

        CHARACTER*4   CYEAR          !  character 4-digit year
        CHARACTER*100 CHARS( 7 )     !  source fields for output
        CHARACTER*200 ACTBUF         !  activity list buffer
        CHARACTER*200 POLBUF         !  pollutant list buffer
        CHARACTER*200 UNTBUF         !  activity units buffer
        CHARACTER*300 ACTVFMT        !  output format buffer for activities
        CHARACTER*300 EMISFMT        !  output format buffer for emissions
        CHARACTER*300 MESG           !  message buffer

        CHARACTER*16 :: PROGNAME = 'WRIDAOUT' ! program name

C***********************************************************************
C   begin body of subroutine WRIDAOUT

C.........  Initializations
        STATUS = 0
        LCOID = -9
        LYEAR = -9

C.........  Rewind temporary files
        DO I = 1, NIPPA
            REWIND( TDEV( I ) )
        END DO

C.........  Write buffers for pollutants and activities to use in headers
        IF( NIPOL .GT. 0 ) THEN
            POLBUF = EANAM( 1 )

            DO I = 2, NIPOL
                L = LEN_TRIM( POLBUF )
                POLBUF = POLBUF( 1:L ) // ' ' // EANAM( I )
            END DO

        END IF

        IF( NIACT .GT. 0 ) THEN
            ACTBUF = ACTVTY( 1 )

            J = INDEX1( ACTVTY( 1 ), NIPPA, EANAM )
            L2 = LEN_TRIM( EAUNIT( J ) ) 
            UNTBUF = '"' // EAUNIT( J )( 1:L2 ) // '"'

            DO I = 2, NIACT
                L = LEN_TRIM( ACTBUF )
                ACTBUF = ACTBUF( 1:L ) // ' ' // ACTVTY( I )

                L = LEN_TRIM( UNTBUF )
                J = INDEX1( ACTVTY( I ), NIPPA, EANAM )
                L2 = LEN_TRIM( EAUNIT( J ) ) 
                UNTBUF = UNTBUF( 1:L )// ' "'// EAUNIT(J)( 1:L2 )// '"'

            END DO

        END IF

C.........  Allocate local memory
        N = MAX( NPPOL, NPACT )
        ALLOCATE( DATARECS( N ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DATARECS', PROGNAME )
        ALLOCATE( COUTRECS( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'COUTRECS', PROGNAME )

C.........  Read data from temporary files and write each IDA record
        SELECT CASE( CATEGORY )

C............................................................................
C.........  For area sources...
        CASE ( 'AREA' )

C.............  Create output format
            WRITE( EMISFMT, 93320 ) NIPOL

C.............  Determine output types (1=real, 0=int)
            OUTTYPE = 1    ! array
            DO I = 1, NPPOL
                K = INDEX( AROFMT( I ), 'I' )
                IF( K .GT. 0 ) THEN
                    OUTTYPE( I ) = 0
                END IF
            END DO

C.............  Write initial header
            L = LEN_TRIM( POLBUF )
            WRITE( DDEV, 93000 ) 
     &             '#TYPE     Area Source Emission Inventory',
     &             '#DESC     Output from SMOKE',
     &             '#POLID    ' // POLBUF( 1:L )

C.............  Write area-source characteristics to output file
            DO S = 1, NSRC

C.................  Store others in temporary variables
                COID = IFIP( S ) / 100000
                FIP  = IFIP( S ) - COID * 100000
                STID = FIP / 1000 
                CYID = FIP - STID * 1000
                SCC  = CSCC ( S )
                YEAR = INVYR( S )

C.................  Write out header
                CALL WRITE_IDA_HEADER( DDEV, IOS )
                IF( IOS .GT. 0 ) CYCLE

C.................  Read data from temporary files
                DO J = 1, NIPOL
                    FDEV = TDEV( J )
                    READ( FDEV, * ) ( DATARECS( I ), I = 1, NPPOL )

                    CALL FORMAT_OUTREC( NPPOL, DATARECS, AROLEN, 
     &                                  AROFMT, COUTRECS( J ) )
                END DO

C.................  Convert data records to ASCII strings to filter out missing
C                   values

C.................  Write out main entries
                WRITE( DDEV, EMISFMT ) STID, CYID, SCC,
     &               ( COUTRECS( J ), J = 1, NIPOL )

                LCOID = COID
                LYEAR = YEAR

            END DO  ! loop through sources

C............................................................................
C.........  For mobile sources...
        CASE ( 'MOBILE' )

C.............  Create output formats
            WRITE( EMISFMT, 93330 ) NIPOL
            WRITE( ACTVFMT, 93335 ) NIACT

C.............  Set output types (1=real, 0=int)
            OUTTYPE = 1    ! array

C.............  Write initial header
            L = LEN_TRIM( POLBUF )
            IF( DDEV .GT. 0 ) WRITE( DDEV, 93000 ) 
     &             '#TYPE     Motor Vehicle Emission Inventory',
     &             '#DESC     Output from SMOKE',
     &             '#POLID    ' // POLBUF( 1:L )

            L1 = LEN_TRIM( ACTBUF )
            L2 = LEN_TRIM( UNTBUF )
            IF( VDEV .GT. 0 ) WRITE( VDEV, 93000 ) 
     &             '#TYPE     Motor Vehicle Activity Inventory',
     &             '#DESC     Output from SMOKE',
     &             '#DATA     ' // ACTBUF( 1:L1 ),
     &             '#UNITS    ' // UNTBUF( 1:L2 )

C.............  Write mobile-source characteristics to output file
            DO S = 1, NSRC

C.................  Store others in temporary variables
                COID  = IFIP( S ) / 100000
                FIP   = IFIP( S ) - COID * 100000
                STID  = FIP / 1000 
                CYID  = FIP - STID * 1000
                CLNK  = CLINK( S )
                SCC   = CSCC ( S )
                YEAR  = INVYR( S )

C.................  Set link to zero if blank
                IF ( CLNK .EQ. ' ' ) CLNK = ADJUSTR( '0' )

C.................  Read emissions from temporary files
                K = 0
                DO J = 1, NIPOL
                    K = K + 1
                    FDEV = TDEV( K )
                    READ( FDEV, * ) ( DATARECS( I ), I = 1, NPPOL )

                    CALL FORMAT_OUTREC( NPPOL, DATARECS, MBOPLEN,
     &                                  MBOPFMT, COUTRECS( J ) )
                END DO

C.................  Write out emissions (fixed formatted)
                IF( DDEV .GT. 0 ) THEN

C.....................  Write out header
                    CALL WRITE_IDA_HEADER( DDEV, IOS )
                    IF( IOS .GT. 0 ) CYCLE            

                    WRITE( DDEV, EMISFMT ) STID, CYID, CLNK, SCC,
     &                   ( COUTRECS( J ), J = 1, NIPOL )

                END IF

C.................  Read activities from temporary files
                DO J = 1, NIACT
                    K = K + 1
                    FDEV = TDEV( K )
                    READ( FDEV, * ) ( DATARECS( I ), I = 1, NPACT )

                    CALL FORMAT_OUTREC( NPACT, DATARECS, MBOALEN, 
     &                                  MBOAFMT, COUTRECS( J ) )
                END DO

C.................  Write out activities (list formatted)
                IF( VDEV .GT. 0 ) THEN

C.....................  Write out header
                    CALL WRITE_IDA_HEADER( VDEV, IOS )
                    IF( IOS .GT. 0 ) CYCLE

                    WRITE( VDEV, ACTVFMT ) STID, CYID, CLNK, SCC,
     &                   ( COUTRECS( J ), J = 1, NIACT )

                END IF                    

                LCOID = COID
                LYEAR = YEAR

            END DO  ! loop through sources

C............................................................................
C.........  For point sources...
        CASE ( 'POINT' )

C.............  Set values of variables that SMOKE does not (yet) support
            BEGYR    = -9
            ENDYR    = -9
            BOILCAP  = -9.
            CAPUNITS = ' '
            WINTHRU  = 0.
            SPRTHRU  = 0.
            SUMTHRU  = 0.
            FALTHRU  = 0.
            HOURS    = 0
            START    = 0
            NDAY     = 0
            WEEKS    = 0
            THRUPUT  = -9.
            MAXRATE  = -9.
            HEATCON  = -9.
            SULFCON  = -9.
            ASHCON   = -9.
            NETDC    = -9.
            OFFSHORE = ' '

C.............  Compute conversion constants
            M2FT  = 1./FT2M
            CTOF  = 1./FTOC

C.............  Create output format
            WRITE( EMISFMT, 93340 ) NIPOL

C.............  Determine output types (1=real, 0=int)
            OUTTYPE = 1    ! array
            DO I = 1, NPPOL
                K = INDEX( PTOFMT( I ), 'I' )
                IF( K .GT. 0 ) THEN
                    OUTTYPE( I ) = 0
                END IF
            END DO

C.............  Write initial header
            L = LEN_TRIM( POLBUF )
            WRITE( DDEV, 93000 ) 
     &             '#TYPE     Point Source Emission Inventory',
     &             '#DESC     Output from SMOKE',
     &             '#POLID    ' // POLBUF( 1:L )

C.............  Write point-source characteristics to output file
            DO S = 1, NSRC

                CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP,
     &                         IDACOLS, NCHAR, CHARS )

C.................  Truncate character string variables
                PLANTID  = CHARS( 2 ) 
                POINTID  = CHARS( 3 )
                STACKID  = CHARS( 4 )
                SEGMENT  = CHARS( 5 )
                SCC      = CHARS( 6 )

                BLRID    = CBLRID( S )
                PLNTDESC = CPDESC( S )

C.................  Convert units of stack parameters
                STKHGT  = STKHT( S ) * M2FT
                STKDIAM = STKDM( S ) * M2FT
                STKTEMP = ( STKTK( S ) - CTOK ) / CTOF + 32.
                STKVEL  = STKVE( S ) * M2FT
                STKFLOW = STKVEL * 0.25 * PI * STKDIAM * STKDIAM

C.................  Store others in temporary variables
                COID = IFIP( S ) / 100000
                FIP  = IFIP( S ) - COID * 100000
                STID = FIP / 1000 
                CYID = FIP - STID * 1000

                SIC    = ISIC ( S )                
                YEAR   = INVYR( S )                
                LATC   = YLOCA( S )                
                LONC   = XLOCA( S )                
                ORISID = CORIS( S )

C.................  Write out header
                CALL WRITE_IDA_HEADER( DDEV, IOS )
                IF( IOS .GT. 0 ) CYCLE

C.................  Read data from temporary files
                DO J = 1, NIPOL
                    FDEV = TDEV( J )
                    READ( FDEV, * ) ( DATARECS( I ), I = 1, NPPOL )

                    CALL FORMAT_OUTREC( NPPOL, DATARECS, PTOLEN, 
     &                                  PTOFMT, COUTRECS( J ) )
                END DO

C.................  Write out main entries
                WRITE( DDEV, EMISFMT ) STID, CYID, PLANTID, POINTID,
     &                 STACKID, ORISID, BLRID, SEGMENT, PLNTDESC, SCC,
     &                 BEGYR, ENDYR, STKHGT, STKDIAM, STKTEMP, STKFLOW, 
     &                 STKVEL, BOILCAP, CAPUNITS, WINTHRU, SPRTHRU, 
     &                 SUMTHRU, FALTHRU, HOURS, START, NDAY, WEEKS, 
     &                 THRUPUT, MAXRATE, HEATCON, SULFCON, ASHCON, 
     &                 NETDC, SIC, LATC, LONC, OFFSHORE, 
     &                 ( COUTRECS( J ), J = 1, NIPOL )

                LCOID = COID
                LYEAR = YEAR

            END DO  ! loop through sources

        END SELECT  ! for source category

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93320   FORMAT( '(I2, I3, A10, ', I3.3, '(A47) )' )

93330   FORMAT( '(I2, I3, A10, A10, ', I3.3, '(A20) )' )           ! mb emis

93335   FORMAT( '(I2, 1X, I3, 1X, A10, A10, 1X,', I3.3, '(A17) )' ) ! mb activity

93340   FORMAT( '(I2, I3, A15, A15, A12, A6, A5, A2, A40, A10, I4, I4,', ! point
     &          'F4.0, F6.2, F4.0, F10.2, F9.2, F8.2, A1, 4F2.0, 2I2, ',
     &          'I1, I2, F11.1, F12.3, F8.2, F5.2, F5.2, F9.3, I4,',
     &          '2F9.4, A1, ', I3.3,'(A52) )' )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for writing the COUNTRY and YEAR
C               parts of the IDA header fields, when the COUNTRY or YEAR are
C               inconsistent with the previous country or year.
            SUBROUTINE WRITE_IDA_HEADER( FDEV, LOCSTAT )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: FDEV    ! file unit no.
            INTEGER, INTENT(OUT) :: LOCSTAT ! exit status

C.............  Subprogram local variables
            INTEGER     K

C-------------------------------------------------------------------------

            LOCSTAT = 0

            IF( COID .NE. LCOID .OR. YEAR .NE. LYEAR ) THEN

                WRITE( CYEAR, '(I4)' ) YEAR

                K = FIND1( COID, NCOUNTRY, CTRYCOD )

                IF( K .GT. 0 ) THEN

                    WRITE( FDEV, 93000 ) 
     &                 '#COUNTRY  ' // CTRYNAM( K ),
     &                 '#YEAR     ' // CYEAR

                ELSE
                    STATUS = 1
                    LOCSTAT = 1
                    WRITE( MESG,94010 ) 'Invalid country code', K,
     &                     'found at source', S
                    CALL M3MESG( MESG )

                END IF

            END IF

C---------------------  FORMAT  STATEMENTS   -----------------------------

C...............   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE WRITE_IDA_HEADER

C-------------------------------------------------------------------------
C-------------------------------------------------------------------------

C.............  This internal subprogram is for formatting the pollutant- and
C               activity-specific data for output to IDA format.  It screens
C               for missing values.
            SUBROUTINE FORMAT_OUTREC( NPPOL, DATAVALS, FMTLEN, 
     &                                FMTBUF, BUFFER )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: NPPOL             ! number of data vals
            REAL        , INTENT (IN) :: DATAVALS( NPPOL ) ! data values
            INTEGER     , INTENT (IN) :: FMTLEN  ( NPPOL ) ! format lengths
            CHARACTER(*), INTENT (IN) :: FMTBUF  ( NPPOL ) ! format info
            CHARACTER(*), INTENT(OUT) :: BUFFER            ! fmtted output

C.............  Subprogram local variables
            INTEGER       I, L

            CHARACTER*10  FMT
            CHARACTER*20  BUFLOC

C-------------------------------------------------------------------------

            BUFFER = ' '

C.............  Build format for first field
            FMT = '(' // FMTBUF( 1 ) // ')'

C.............  Initialize output buffer with first value if it is not missing
            IF( DATAVALS( 1 ) .GT. BADVAL3 ) 
     &          WRITE( BUFFER, FMT ) DATAVALS( 1 )
            L = FMTLEN( 1 )

C.............  For remaining fields, build formats and add value to buffer
C               for non-missing values
            DO I = 2, NPPOL

                FMT = '(' // FMTBUF( I ) // ')'

                IF( OUTTYPE ( I ) .EQ. 0      .AND.
     &              DATAVALS( I ) .GT. IMISS3       ) THEN

                    WRITE( BUFLOC, FMT ) INT( DATAVALS( I ) )
                    BUFFER = BUFFER( 1:L ) // BUFLOC

                ELSE IF ( OUTTYPE ( I ) .NE. 0       .AND.
     &                    DATAVALS( I ) .GT. BADVAL3       ) THEN
                    WRITE( BUFLOC, FMT ) DATAVALS( I )
                    BUFFER = BUFFER( 1:L ) // BUFLOC

                END IF

                L = L + FMTLEN( I )

            END DO

C---------------------  FORMAT  STATEMENTS   -----------------------------

C...............   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE FORMAT_OUTREC

        END SUBROUTINE WRIDAOUT

