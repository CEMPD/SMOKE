
        SUBROUTINE RDIDAAR( FDEV, NRAWBP, WKSET,
     &                      CURREC, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine reads the IDA format area-source inventory
C      files.  It can read multiple IDA files, if needed.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      created by M. Houyoux (04/99) 
C
C**************************************************************************
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
         INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

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

        EXTERNAL    CRLF, ENVINT, ENVYN, GETFLINE, GETNLIST, INDEX1, 
     &              STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file
        INTEGER     , INTENT (IN) :: NRAWBP ! total raw record times pols
        INTEGER     , INTENT (IN) :: WKSET  ! weekly profile interpretation
        INTEGER     ,INTENT(INOUT):: CURREC ! current no. source * pollutants
        LOGICAL     , INTENT(OUT) :: EFLAG  ! outgoing error flag
        INTEGER     ,INTENT(INOUT):: NDROP  !  number of records dropped
        REAL        ,INTENT(INOUT):: EDROP( MXIDAT )  ! emis dropped per pol

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 63  ! maximum pollutants in file
        INTEGER, PARAMETER :: AROTWIDE = 47  ! total width of all pol fields
        INTEGER, PARAMETER :: ARNONPWD = 15  ! width of non-pol fields

C...........   Local parameters, dependent
        INTEGER, PARAMETER :: LINSIZ  = ARNONPWD + MXPOLFIL * AROTWIDE

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER, PARAMETER :: ISINIT( NARPPOL3 ) = 
     &                              ( / 16,26,36,47,54,57 / )

        INTEGER, PARAMETER :: IEINIT( NARPPOL3 ) = 
     &                              ( / 25,35,46,53,56,62 / )

C...........   Local arrays
        INTEGER          IS( NARPPOL3 )  ! start position for each pol char
        INTEGER          IE( NARPPOL3 )  ! end position for each pol char

C...........   Other local variables
        INTEGER         I, J, K, L, N, V  ! counters and indices

        INTEGER         COD     !  tmp pollutant position in INVDNAM
        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP     !  tmp FIPS code
        INTEGER         ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER         LDEV    !  unit no. for log file
        INTEGER, SAVE:: MXWARN  !  maximum number of warnings
        INTEGER         NPOL    !  number of pollutants in file
        INTEGER, SAVE:: NWARN =0!  number of warnings in this routine
        INTEGER         NWRLINE !  number of lines of file written to log
        INTEGER         TPF     !  tmp temporal adjustments setting

        REAL            CEFF    !  tmp control effectiveness
        REAL            DAY2YR  !  factor to convert from daily data to annual
        REAL            EANN    !  tmp annual-ave emission value
        REAL            EMFC    !  tmp emission factor
        REAL            EOZN    !  tmp ozone-season-ave emission value
        REAL            REFF    !  tmp rule effectiveness
        REAL            RPEN    !  tmp rule penetration

        LOGICAL, SAVE:: FFLAG    = .FALSE. ! true: fill in 0. annual with seasonal
        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called

        CHARACTER*20    CNTRY   !  country name
        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIP code
        CHARACTER(LEN=IOVLEN3) CBUF  ! tmp pollutant code
        CHARACTER(LEN=LINSIZ)  LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp scc
        
        CHARACTER(LEN=300)     TENLINES( 10 ) ! first ten lines of inventory file

        CHARACTER*16 :: PROGNAME = 'RDIDAAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDIDAAR

        IF ( FIRSTIME ) THEN

            MESG = 'Fill in 0. annual data based on seasonal data.'
            FFLAG = ENVYN( 'FILL_ANN_WSEAS', MESG, .FALSE., IOS )

            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

        END IF

C.........  Reinitialize for multiple subroutine calls
        EFLAG = .FALSE.
        ICC   = -9
        INY   = 0
        NPOL  = 0
        
C.........  Get log file number for reports
        LDEV = INIT3()        
        NWRLINE = 0

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        ES   = CURREC
        IREC = 0
        DO

C.............  Read a line of IDA file as a character string
            READ( FDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            L = LEN_TRIM( LINE )  ! store width of line and check

C.............  Skip blank lines
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set 
C               properly
            CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
     &                   LINE, ICC, INY, NPOL, IOS )

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

C.............  Write first ten lines to log file
            IF( NWRLINE < 10 ) THEN
            	NWRLINE = NWRLINE + 1
            	TENLINES( NWRLINE ) = TRIM( LINE )
            
                IF( NWRLINE == 10 ) THEN
                    MESG = 'First 10 lines of IDA area inventory:'
                    WRITE( LDEV,* ) TRIM( MESG )
             
                    DO I = 1,NWRLINE
                        WRITE( LDEV,* ) TRIM( TENLINES( I ) )
                    END DO
                END IF
            END IF
            
C.............  Make sure that all of the needed integer values are integers...

C.............  Check state/county codes, error for missing
            IF( .NOT. CHKINT( LINE( 1:2 ) ) .OR. 
     &          .NOT. CHKINT( LINE( 3:5 ) ) .OR.
     &          LINE( 1:2 ) .EQ. ' '        .OR.
     &          LINE( 3:5 ) .EQ. ' '             ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: State and/or county ' //
     &                 'code is non-integer or missing at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Initialize start and end positions
            IS = ISINIT - AROTWIDE  ! array
            IE = IEINIT - AROTWIDE  ! array

C.............  Make sure that all of the needed real values are real...

C.............  Emissions and associated data
            DO V = 1, NPOL

C.................  Update start and end positions
                DO K = 1, NARPPOL3
                    IS( K ) = IS( K ) + AROTWIDE
                    IE( K ) = IE( K ) + AROTWIDE
                END DO

                IF( .NOT. CHKREAL( LINE( IS(1):IE(1) ) ) .OR.
     &              .NOT. CHKREAL( LINE( IS(2):IE(2) ) ) .OR.
     &              .NOT. CHKREAL( LINE( IS(3):IE(3) ) ) .OR.
     &              .NOT. CHKREAL( LINE( IS(4):IE(4) ) ) .OR.
     &              .NOT. CHKREAL( LINE( IS(5):IE(5) ) )      ) THEN

                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM( V ) )
                    WRITE( MESG,94010 ) 'ERROR: Emission data, ' //
     &                     'control percentages, and/or emission ' //
     &                     CRLF() // BLANK10 // 'factor for' //
     &                     TMPNAM( V )( 1:L ) // ' are not a number ' //
     &                     'or have bad formatting at line', IREC
                    CALL M3MESG( MESG )
                    NDROP = NDROP + 1

                END IF

                IF( LINE( IS(1):IE(1) ) .EQ. ' ' .AND.
     &              LINE( IS(2):IE(2) ) .EQ. ' '       ) THEN

                    L = LEN_TRIM( TMPNAM( V ) )
                    WRITE( MESG,94010 ) 'WARNING: All emissions ' //
     &                     'data for ' // TMPNAM( V )( 1:L ) //  
     &                     ' are missing at line', IREC
                    CALL M3MESG( MESG )
                    LINE( IS(1):IE(1) ) = '0.'
                    LINE( IS(2):IE(2) ) = '0.'

                END IF

            END DO

C.............  If there has been an error, do not try to store any of the
C               records.  Instead  go to next line of file.
            IF( EFLAG ) CYCLE

            IF( FIRSTIME ) DAY2YR  = 1. / YR2DAY( INY )
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields...

            FIP  = ICC * 100000 + 1000 * STR2INT( LINE( 1:2 ) ) +
     &             STR2INT( LINE( 3:5 ) )

            TSCC = LINE( 6:15 )    ! SCC code

C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Initialize start and end positions
            IS = ISINIT - AROTWIDE  ! array
            IE = IEINIT - AROTWIDE  ! array

C.............  Loop through pollutants and store data so that there is one
C               record for each pollutant.  This will be consistent with
C               the other reader routines.
            DO V = 1, NPOL

C.................  Set the default temporal resolution of the data
                TPF  = MTPRFAC * WKSET
            
                CBUF = TMPNAM( V )
                L = LEN_TRIM( CBUF )

C.................  Update start and end positions
                DO K = 1, NARPPOL3
                    IS( K ) = IS( K ) + AROTWIDE
                    IE( K ) = IE( K ) + AROTWIDE
                END DO

C.................  Read annual emissions for pollutant V
                N = IE(1) - IS(1) + 1
                CALL READ_REAL( N, IREC, .TRUE., LINE( IS(1):IE(1) ), 
     &                          CBUF( 1:L ) // ' annual emissions', 
     &                          EANN, EFLAG )

C.................  Read seasonal emissions for pollutant V
                N = IE(2) - IS(2) + 1
                CALL READ_REAL( N, IREC, .TRUE., LINE( IS(2):IE(2) ), 
     &                          CBUF( 1:L ) // ' seasonal emissions', 
     &                          EOZN, EFLAG )

C.................  Read emission factor for pollutant V
                N = IE(3) - IS(3) + 1
                CALL READ_REAL( 10, IREC, .TRUE., LINE( IS(3):IE(3) ), 
     &                          CBUF( 1:L ) // ' emission factor', 
     &                          EMFC, EFLAG )

C.................  Read control efficiency for pollutant V
                N = IE(4) - IS(4) + 1
                CALL READ_REAL( N, IREC, .TRUE., LINE( IS(4):IE(4) ), 
     &                          CBUF( 1:L ) // ' control efficiency', 
     &                          CEFF, EFLAG )

C.................  Read rule effectiveness for pollutant V
                N = IE(5) - IS(5) + 1
                CALL READ_REAL( N, IREC, .TRUE., LINE( IS(5):IE(5) ), 
     &                          CBUF( 1:L ) // ' rule effectiveness', 
     &                          REFF, EFLAG )

C.................  Read rule penetration for pollutant V
                N = IE(6) - IS(6) + 1
                CALL READ_REAL( N, IREC, .TRUE., LINE( IS(6):IE(6) ), 
     &                          CBUF( 1:L ) // ' rule penetration', 
     &                          RPEN, EFLAG )

                IF( EANN .LT. AMISS3 ) THEN

                    IF( NWARN .LT. MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Missing annual '//
     &                         'emissions at line', IREC, 'for ' // CBUF
                        NWARN = NWARN + 1
                        CALL M3MESG( MESG )
                    END IF

                    IF ( EOZN .LT. AMISS3 ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Missing annual '//
     &                         'AND seasonal emissions at ' //
     &                         'line', IREC, 'for ' // CBUF
                        NWARN = NWARN + 1
                        CALL M3MESG( MESG )
                    END IF

                ELSE IF ( EOZN .LT. AMISS3 ) THEN

                    IF( NWARN .LT. MXWARN ) THEN
                        WRITE( MESG,94010 ) 
     &                         'WARNING: Missing seasonal '//
     &                         'emissions at line', IREC, 'for ' // CBUF
                        NWARN = NWARN + 1
                        CALL M3MESG( MESG )
                    END IF

                END IF

                IF( CEFF .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing control ' //
     &                     'efficiency at line', IREC, 'for ' // CBUF 
                    IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF
                
                IF( REFF .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                     'effectiveness at line', IREC, 'for '// CBUF
                    IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF

                IF( RPEN .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                     'pentration at line', IREC, 'for ' // CBUF 
                    IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF

C.................  Replace annual data with ozone-season information if
C                   user option is set
                IF( FFLAG        .AND. 
     &              EANN .LE. 0. .AND.
     &              EOZN .GT. 0.       ) THEN

                    WRITE( MESG,94010 ) 'NOTE: Using seasonal ' //
     &                     'emissions to fill in annual emissions' //
     &                     CRLF() // BLANK10 // 'at line', IREC,
     &                     'for ' // CBUF
                    CALL M3MESG( MESG )

                    EANN = EOZN * DAY2YR

C.....................  Remove monthly factors for this source.
                    TPF  = WKSET

                END IF

C.................  Store data in final arrays if there is enough memory
                ES = ES + 1

                IF ( ES .LE. NRAWBP ) THEN

                    J = DATPOS( V )
                    IFIPA  ( ES     ) = FIP
                    TPFLGA ( ES     ) = TPF
                    INVYRA ( ES     ) = INY
                    CSCCA  ( ES     ) = TSCC                     
                    POLVLA ( ES,NEM ) = INVDCNV( J ) * EANN
                    POLVLA ( ES,NOZ ) = EOZN
                    POLVLA ( ES,NEF ) = EMFC
                    POLVLA ( ES,NCE ) = CEFF
                    POLVLA ( ES,NRE ) = REFF
                    POLVLA ( ES,NRP ) = RPEN

                    WRITE( CCOD,94125 ) J

                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                            CCOD, CSOURCA( ES ) )

                END IF  !  if ES in range

            END DO      !  end of loop through pollutants

        END DO          !  to head of FDEV-read loop

199     CONTINUE        !  exit from the FDEV-read loop

        CLOSE( FDEV )

        WRITE( MESG,94010 ) 
     &         'IDA FILE processed:'  // CRLF() // BLANK10 //
     &              'This-file source*pollutant-count', ES - CURREC,
     &         CRLF() // BLANK10 //
     &              'Cumulative source*pollutant-count', ES

        CALL M3MSG2( MESG )

C.........  Write message if overflow occurred
        IF( ES .GT. NRAWBP ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source by pollutant memory ' //
     &             'allocation insufficient for IDA inventory'
            CALL M3MSG2( MESG )

        ELSE
            CURREC = ES

        END IF

C.........  Deallocate local allocatable arrays 
        DEALLOCATE( TMPNAM, DATPOS )

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        CONTAINS

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
                    RETURN
                END IF
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


        END SUBROUTINE RDIDAAR
