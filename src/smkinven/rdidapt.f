
        SUBROUTINE RDIDAPT( FDEV, NRAWIN, NRAWBP, WKSET,
     &                      NRAWOUT, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 186
C
C  DESCRIPTION:
C      This subroutine reads the IDA format point-source inventory
C      files.  It call be called multiple times for multiple files.
C
C  PRECONDITIONS REQUIRED:
C      Files must be opened and their unit numbers stored in FDEV() in the
C      order listed in the description.  OTHER?
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      created by M. Houyoux (04/99) 
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the point source inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
         INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
      	LOGICAL                CHKINT
        LOGICAL                CHKREAL
        CHARACTER*2            CRLF
        INTEGER                ENVINT
        LOGICAL                ENVYN
        INTEGER                GETFLINE
        INTEGER                GETNLIST
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY 

        EXTERNAL    CRLF, ENVYN, GETFLINE, GETNLIST, INDEX1, 
     &              STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file
        INTEGER     , INTENT (IN) :: NRAWIN ! total raw record-count 
        INTEGER     , INTENT (IN) :: NRAWBP ! total raw record times pols
        INTEGER     , INTENT (IN) :: WKSET  ! weekly profile interpretation
        INTEGER     , INTENT(OUT) :: NRAWOUT! outgoing source * pollutants
        LOGICAL     , INTENT(OUT) :: EFLAG  ! outgoing error flag
        INTEGER     ,INTENT(INOUT):: NDROP  !  number of records dropped
        REAL        ,INTENT(INOUT):: EDROP( MXIDAT )  ! emis dropped per pol

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: BLIDLEN  = 6   ! width of boiler field
        INTEGER, PARAMETER :: DESCLEN  = 40  ! width of plant description field
        INTEGER, PARAMETER :: MXPOLFIL = 53  ! maximum pollutants in file
        INTEGER, PARAMETER :: PTOTWIDE = 52  ! total width of all pol fields
        INTEGER, PARAMETER :: PTNONPWD = 249 ! width of non-pol fields

C...........   Local parameters, dependent
        INTEGER, PARAMETER :: LINSIZ  = PTNONPWD + MXPOLFIL * PTOTWIDE

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER, PARAMETER :: ISINIT( NPTPPOL3 ) = 
     &                              ( / 250,263,276,283,286,296,299 / )

        INTEGER, PARAMETER :: IEINIT( NPTPPOL3 ) = 
     &                              ( / 262,275,282,285,295,298,301 / )

C...........   Local arrays
        INTEGER          IS( NPTPPOL3 )  ! start position for each pol char
        INTEGER          IE( NPTPPOL3 )  ! end position for each pol char

C...........   Counters of total number of input records
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count
        INTEGER, SAVE :: NSRCPOL = 0 ! cumulative source x pollutants count

C.........  Temporary variables for storing source characteristics.  These
C           variables must be the width of the fields for global source
C           characteristics definition for use in BLDCSRC.
        CHARACTER(LEN=PLTLEN3) FCID  ! tmp plant ID
        CHARACTER(LEN=CHRLEN3) PTID  ! tmp point ID
        CHARACTER(LEN=CHRLEN3) SKID  ! tmp stack ID
        CHARACTER(LEN=CHRLEN3) SGID  ! tmp segment ID

C.........  Local allocatable arrays

C...........   Other local variables
        INTEGER         I, J, K, L, V  ! counters and indices

        INTEGER         CNY     !  county code
        INTEGER         COD     !  tmp pollutant position in INVDNAM
        INTEGER         CPRI    !  tmp primary control device code
        INTEGER         CSEC    !  tmp secondary control device code
        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP, SCC, SIC  ! tmp fip, scc, sic
        INTEGER         ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER, SAVE:: MXWARN  !  maximum number of warnings
        INTEGER         NPOL    !  number of pollutants in file
        INTEGER, SAVE:: NWARN =0!  number of warnings in this routine
        INTEGER         SS      !  counter for sources
        INTEGER         STA     !  state code
        INTEGER         TPF     !  tmp temporal adjustments setting

        REAL            CEFF    !  tmp control effectiveness
        REAL            DAY2YR  ! factor to convert from daily data to annual
        REAL            EANN    !  tmp annual-ave emission value
        REAL            EMFC    !  tmp emission factor
        REAL            EOZN    !  tmp ozone-season-ave emission value
        REAL            LAT     !  tmp Y-coordinate 
        REAL            LON     !  tmp X-coordinate
        REAL            RBUF    !  tmp real value
        REAL            REFF    !  tmp rule effectiveness
        REAL            DM, HT, FL, TK, VL  ! Temporary stack parms

        LOGICAL, SAVE :: CFLAG              ! true: recalc vel w/ flow & diam
        LOGICAL, SAVE :: FFLAG    = .FALSE. ! true: fill in 0. annual with seasonal
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: 1st time routine called
        LOGICAL, SAVE :: WFLAG    = .FALSE. ! true: all lat-lons to western hemi

        CHARACTER*20    CNTRY   !  country name
        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=BLIDLEN)  BLID  ! tmp boiler ID
        CHARACTER(LEN=IOVLEN3)  CBUF  ! tmp pollutant name
        CHARACTER(LEN=POLLEN3)  CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3)  CFIP  ! character FIP code
        CHARACTER(LEN=CHRLEN3)  CHAR4 ! tmp characteristic 4
        CHARACTER(LEN=ORSLEN3)  CORS  ! tmp DOE plant ID
        CHARACTER(LEN=IOVLEN3)  CPOL  ! tmp pollutant code
        CHARACTER(LEN=DESCLEN)  DESC  ! tmp plant description
        CHARACTER(LEN=PTNONPWD) LINPT1 ! non-emissions part of format
        CHARACTER(LEN=PTOTWIDE) LINEMS ! emissions part of format
        CHARACTER(LEN=SCCLEN3)  TSCC  ! tmp scc

        CHARACTER*16 :: PROGNAME = 'RDIDAPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDIDAPT

C.........  Set up settings the first time the subroutine is called
        IF( FIRSTIME ) THEN

C.............  Get settings from the environment
            CFLAG = ENVYN( 'VELOC_RECALC', 
     &                 'Flag for recalculating velocity', .FALSE., IOS )

            WFLAG = ENVYN( 'WEST_HSPHERE',
     &                 'Western hemisphere flag', .TRUE., IOS )

            MESG = 'Fill in 0. annual data based on seasonal data.'
            FFLAG = ENVYN( 'FILL_ANN_WSEAS', MESG, .FALSE., IOS )

            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

        ENDIF

C.........  Reinitialize for multiple subroutine calls
        EFLAG = .FALSE.
        ICC   = -9
        INY   = 0
        NPOL  = 0

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        SS   = NSRCSAV
        ES   = NSRCPOL
        IREC = 0
        DO

C.............  Read a line of IDA file as a character string
C.............  If line is a header line or blank, it will advance anyway
            READ( FDEV, 93010, END=199, IOSTAT=IOS, ADVANCE="NO", 
     &            EOR=1001 ) 
     &            LINPT1
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file before column ', PTNONPWD,
     &              'at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            L = LEN_TRIM( LINPT1 )  ! store width of line and check

C.............  Skip blank lines
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set 
C               properly
            CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
     &                   LINPT1, ICC, INY, NPOL, IOS )

C.............  Interpret error status
            IF( IOS .EQ. 4 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Maximum allowed data variables ' //
     &                 '(MXPOLFIL=', MXPOLFIL, CRLF() // BLANK10 //
     &                 ') exceeded in input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Read state and county code
            CALL READ_INTEGER( 2, IREC, .FALSE.,  LINPT1( 1:2 ), 
     &                         'state code' , STA, EFLAG )
            STA = MAX( STA, 0 )

            CALL READ_INTEGER( 3, IREC, .FALSE., LINPT1( 3:5 ), 
     &                         'county code', CNY, EFLAG )
            CNY = MAX( CNY, 0 )

C.............  Read stack height and convert units
            CALL READ_REAL( 4, IREC, .TRUE., LINPT1( 120:123 ), 
     &                      'stack height' , HT, EFLAG )
            IF( HT .LT. 0. ) HT = 0.
            HT   = HT * FT2M                 ! ft to m

C.............  Read stack diameter and convert units
            CALL READ_REAL( 6, IREC, .TRUE., LINPT1( 124:129 ), 
     &                      'stack diameter' , DM, EFLAG )
            IF( DM .LT. 0. ) DM = 0.
            DM   = DM * FT2M                 ! ft to m

C.............  Read exit temperature and convert units
            CALL READ_REAL( 4, IREC, .TRUE., LINPT1( 130:133 ), 
     &                      'stack exit temperature' , TK, EFLAG )
            IF( TK .LT. 0. ) THEN
                TK = 0.
            ELSE
                TK   = ( TK - 32 ) * FTOC + CTOK ! F to K
            END IF

C.............  Read exit flow rate and convert units
            CALL READ_REAL( 10, IREC, .TRUE., LINPT1( 134:143 ), 
     &                      'stack exit temperature' , FL, EFLAG )
            IF( FL .LT. 0. ) FL = 0.
            FL   = FL * FT2M3                ! ft^3/s to m^3/s

C.............  Read exit velocity and convert units
            CALL READ_REAL( 7, IREC, .TRUE., LINPT1( 144:152 ), 
     &                      'stack exit velocity' , VL, EFLAG )

C.............  Recompute velocity if that option has been set
            IF ( CFLAG .OR. VL .EQ. 0. ) THEN
                RBUF = ( 0.25 * PI * DM * DM )
                IF ( RBUF .GT. 0 ) VL = FL / RBUF

C.............  Otherwise, convert units
            ELSE
                VL = VL * FT2M   ! ft/s to m/s
            END IF

C.............  Read Standard Industrial Code
            CALL READ_INTEGER( 4, IREC, .TRUE., LINPT1( 227:230 ), 
     &                         'SIC', SIC, EFLAG )

C.............  If SIC is missing, set to zero
            SIC = MAX( SIC, 0 )

C.............  Read stack latitude
            CALL READ_REAL( 9, IREC, .FALSE., LINPT1( 231:239 ), 
     &                      'latitude' , LAT, EFLAG )

C.............  Read stack longitude and correct hemisphere if necessary
            CALL READ_REAL( 9, IREC, .FALSE., LINPT1( 240:248 ), 
     &                      'longtiude' , LON, EFLAG )
            IF( WFLAG .AND. LON .GT. 0 ) LON = -LON
       
            IF( FIRSTIME ) DAY2YR  = 1. / YR2DAY( INY )
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields that remain

            FIP  = ICC * 100000 + 1000 * STA + CNY
            FCID = ADJUSTL( LINPT1(   6:20  ) )  ! plant ID
            PTID = ADJUSTL( LINPT1(  21:35  ) )  ! point ID
            SKID = ADJUSTL( LINPT1(  36:47  ) )  ! stack ID
            CORS = ADJUSTL( LINPT1(  48:53  ) )  ! DOE plant ID
            BLID = ADJUSTL( LINPT1(  54:59  ) )  ! boiler ID
            SGID = ADJUSTL( LINPT1(  60:61  ) )  ! segment ID
            DESC = ADJUSTL( LINPT1(  62:101 ) )  ! plant description
            TSCC = ADJUSTL( LINPT1( 102:111 ) )  ! SCC code

C.............  Set the default temporal resolution of the data
            TPF  = MTPRFAC * WKSET

C.............  Increment source number
            SS = SS + 1

C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Loop through pollutants and store data so that there is one
C               record for each pollutant.  This will be consistent with
C               the other reader routines.
            DO V = 1, NPOL

C.................  Non-advancing read for all but the last pollutant and
C                   advancing read for the last pollutant
                IF ( V .LE. NPOL ) THEN
                    READ( FDEV, 93020, END=199, IOSTAT=IOS, 
     &                    ADVANCE="NO", EOR = 1003 ) LINEMS
                ELSE
                    READ( FDEV, 93020, END=199, IOSTAT=IOS ) 
     &                  LINEMS

                END IF

                CBUF = TMPNAM( V )
                L = LEN_TRIM( CBUF )

C.................  Read annual emissions for pollutant V
                CALL READ_REAL( 13, IREC, .TRUE., LINEMS( 1:13 ), 
     &                          CBUF( 1:L ) // ' annual emissions', 
     &                          EANN, EFLAG )

C.................  Read seasonal emissions for pollutant V
                CALL READ_REAL( 13, IREC, .TRUE., LINEMS( 14:26 ), 
     &                          CBUF( 1:L ) // ' seasonal emissions', 
     &                          EOZN, EFLAG )

C.................  Read control efficiency for pollutant V
                CALL READ_REAL( 7, IREC, .TRUE., LINEMS( 27:33 ), 
     &                          CBUF( 1:L ) // ' control efficiency', 
     &                          CEFF, EFLAG )

C.................  Read rule effectiveness for pollutant V
                CALL READ_REAL( 3, IREC, .TRUE., LINEMS( 34:36 ), 
     &                          CBUF( 1:L ) // ' rule effectiveness', 
     &                          REFF, EFLAG )

C.................  Read emission factor for pollutant V
                CALL READ_REAL( 10, IREC, .TRUE., LINEMS( 37:46 ), 
     &                          CBUF( 1:L ) // ' emission factor', 
     &                          EMFC, EFLAG )

C.................  Read primary control code for pollutant V
                CALL READ_INTEGER( 3, IREC, .TRUE., LINEMS( 47:49 ), 
     &                             CBUF(1:L)// ' primary control code', 
     &                             CPRI, EFLAG )

C.................  Read secondary control code for pollutant V
                CALL READ_INTEGER( 3, IREC, .TRUE., LINEMS( 50:52 ), 
     &                             CBUF(1:L)//' secondary control code', 
     &                             CSEC, EFLAG )

C.................  If there has been an error, do not try to finish storing 
C                   the records.  Instead  go to next line of file.
                IF( EFLAG ) CYCLE

C.................  Warning if missing both annual and ozone-season information
                IF( EANN  .LT. 0.     .AND.
     &              EOZN  .LT. 0.     .AND.
     &              NWARN .LT. MXWARN       ) THEN

                    WRITE( MESG,94010 ) 'WARNING: Missing annual '//
     &                     'AND seasonal emissions at ' //
     &                     'line', IREC, 'for ' // CPOL
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
                END IF

C.................  Replace annual data with ozone-season information if
C                   user option is set
                IF( FFLAG        .AND. 
     &              EANN .LE. 0. .AND.
     &              EOZN .GT. 0.       ) THEN

                    WRITE( MESG,94010 ) 'NOTE: Using seasonal ' //
     &                     'emissions to fill in annual emissions' //
     &                     CRLF() // BLANK10 // 'at line', IREC,
     &                     'for "' // CBUF( 1:L ) // '"'
                    CALL M3MESG( MESG )

                    EANN = EOZN * DAY2YR

C.....................  Remove monthly factors for this source. Note that this
C                       will impact ALL pollutants, even if only one pollutant
C                       gets filled.  This is necessary unless TPF is changed
C                       to be pollutant-dependent.
                    TPF  = WKSET

                END IF

C.................  Store data in final arrays if there is enough memory
                ES = ES + 1

                IF ( ES .LE. NRAWBP ) THEN

                    J = DATPOS( V )
                    INDEXA ( ES     ) = ES
                    INRECA ( ES     ) = SS                    
                    POLVLA ( ES,NEM ) = INVDCNV( J ) * EANN
                    POLVLA ( ES,NOZ ) = EOZN
                    POLVLA ( ES,NCE ) = CEFF
                    POLVLA ( ES,NRE ) = REFF
                    POLVLA ( ES,NEF ) = EMFC
                    POLVLA ( ES,NC1 ) = REAL( CPRI ) ! store as real for now
                    POLVLA ( ES,NC2 ) = REAL( CSEC ) ! store as real for now
                    
                    WRITE( CCOD,94125 ) J

                    CHAR4 = TSCC 
                    CALL BLDCSRC( CFIP, FCID, PTID, SKID, SGID, 
     &                            CHAR4, CHRBLNK3, CCOD, CSOURCA( ES ) )

                END IF  !  if ES in range

            END DO      !  end of loop through pollutants

C.............  Store source characteristics if dimension is okay
C.............  This is done after pollutants stored because TPF might change
C               in pollutants loop
            IF( SS .LE. NRAWIN ) THEN

                IFIPA  ( SS ) = FIP
                ISICA  ( SS ) = SIC
                TPFLGA ( SS ) = TPF
                INVYRA ( SS ) = INY
                STKHTA ( SS ) = HT
                STKDMA ( SS ) = DM
                STKTKA ( SS ) = TK
                STKVEA ( SS ) = VL
                XLOCAA ( SS ) = LON
                YLOCAA ( SS ) = LAT
                CSCCA  ( SS ) = TSCC
                CORISA ( SS ) = CORS
                CBLRIDA( SS ) = BLID
                CPDESCA( SS ) = DESC

            END IF 

        END DO          !  to head of FDEV-read loop

199     CONTINUE        !  exit from the FDEV-read loop

        CLOSE( FDEV )

        WRITE( MESG,94010 ) 
     &         'IDA FILE processed:'  // CRLF() // BLANK10 //
     &              'This-file source-count', SS - NSRCSAV,
     &         CRLF() // BLANK10 //
     &              'Cumulative source-count', SS,
     &         CRLF() // BLANK10 //
     &              'This-file source*pollutant-count', ES - NSRCPOL,
     &         CRLF() // BLANK10 //
     &              'Cumulative source*pollutant-count', ES

        CALL M3MSG2( MESG )

C.........  Update saved cumulative counts
        NSRCSAV = SS        !  source
        NSRCPOL = ES        !  source*pollutant

C.........  Write message if overflow occurred
        IF( NSRCSAV .GT. NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source memory allocation ' //
     &             'insufficient for IDA inventory'
            CALL M3MSG2( MESG )

        END IF

        IF( NSRCPOL .GT. NRAWBP ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source by pollutant memory ' //
     &             'allocation insufficient for IDA inventory'
            CALL M3MSG2( MESG )

        ELSE
            NRAWOUT = NSRCPOL

        END IF

C.........  Deallocate local allocatable arrays 
        DEALLOCATE( TMPNAM, DATPOS )

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

1001    WRITE( MESG,94010 ) 'End of line encountered ' //
     &                      'unexpectedly at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

1003    WRITE( MESG,94010 ) 'End of line encountered ' //
     &                      'unexpectedly for variable', V, 
     &                      'at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A )  

93020   FORMAT( A52 )    ! must match value of PTOTWIDE

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 ) 

        CONTAINS

            SUBROUTINE READ_INTEGER( LENGTH, IREC, OFLAG, STRING, DESC, 
     &                               INTVAL, EFLAG )

            INTEGER      , INTENT ( IN )  :: LENGTH
            INTEGER      , INTENT ( IN )  :: IREC
            LOGICAL      , INTENT ( IN )  :: OFLAG  !  true: field is optional 
            CHARACTER(LEN=LENGTH), INTENT( IN ) :: STRING
            CHARACTER*(*), INTENT ( IN )  :: DESC
            INTEGER      , INTENT ( OUT ) :: INTVAL
            LOGICAL      , INTENT ( OUT ) :: EFLAG 

C.............  Check to see if the field is blank
            IF ( STRING .EQ. ' ' ) THEN

C.................  If field is blank and optional, then set to missing
                IF ( OFLAG ) THEN
                    INTVAL = IMISS3

C.................  If field is blank and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' is blank at line', IREC
                    CALL M3MESG( MESG )
                    INTVAL = IMISS3
                    RETURN
                END IF
            END IF

C.............  Try to read value and see what the error status is
            READ( STRING, *, IOSTAT = IOS ) INTVAL

C.............  If error, then write message and continue
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: ' // DESC //
     &                 ' has non-integer value "' // STRING //
     &                 '"' // CRLF() // BLANK 10 // 'at line', IREC
                CALL M3MESG( MESG )
                INTVAL = IMISS3

C.............  Check if missing value
            ELSE IF ( INTVAL .EQ. -9 ) THEN

C.................  If field is missing and optional, then set to missing
                IF ( OFLAG ) THEN
                    INTVAL = IMISS3

C.................  If field is missing and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' has -9 missing value at line', IREC
                    CALL M3MESG( MESG )
                    INTVAL = IMISS3

                END IF

            END IF

            RETURN

C............................................................................

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_INTEGER

C--------------------------------------------------------------------------
C--------------------------------------------------------------------------

            SUBROUTINE READ_REAL( LENGTH, IREC, OFLAG, STRING, DESC, 
     &                            REALVAL, EFLAG )

            INTEGER      , INTENT ( IN )  :: LENGTH
            INTEGER      , INTENT ( IN )  :: IREC
            LOGICAL      , INTENT ( IN )  :: OFLAG  !  true: field is optional 
            CHARACTER(LEN=LENGTH), INTENT( IN ) :: STRING
            CHARACTER*(*), INTENT ( IN )  :: DESC
            REAL         , INTENT ( OUT ) :: REALVAL
            LOGICAL      , INTENT ( OUT ) :: EFLAG 

C.............  Check to see if the field is blank
            IF ( STRING .EQ. ' ' ) THEN

C.................  If field is blank and optional, then set to missing
                IF ( OFLAG ) THEN
                    REALVAL = BADVAL3

C.................  If field is blank and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' is blank at line', IREC
                    CALL M3MESG( MESG )
                    REALVAL = BADVAL3
                END IF

                RETURN

            END IF

C.............  Try to read value and see what the error status is
            READ( STRING, *, IOSTAT = IOS ) REALVAL

C.............  If error, then write message and continue
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: ' // DESC //
     &                 ' has non-readable value "' // STRING //
     &                 '"' // CRLF() // BLANK 10 // 'at line', IREC
                CALL M3MESG( MESG )
                REALVAL = BADVAL3

C.............  Check if missing value
            ELSE IF ( REALVAL .EQ. -9 ) THEN

C.................  If field is missing and optional, then set to zero
                IF ( OFLAG ) THEN
                    REALVAL = BADVAL3

C.................  If field is missing and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' has -9 missing value at line', IREC
                    CALL M3MESG( MESG )
                    REALVAL = BADVAL3

                END IF

            END IF

            RETURN

C............................................................................

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_REAL


        END SUBROUTINE RDIDAPT
