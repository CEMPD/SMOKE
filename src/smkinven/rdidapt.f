
        SUBROUTINE RDIDAPT( FDEV, NRAWIN, NRAWBP, WKSET,
     &                      NRAWOUT, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 
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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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

        EXTERNAL    CRLF, ENVYN, GETFLINE, GETNLIST, INDEX1, 
     &              STR2INT, STR2REAL

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

C...........   Other local variables
        INTEGER         I, J, K, L, V  ! counters and indices

        INTEGER         COD     !  tmp pollutant position in INVDNAM
        INTEGER         CPRI    !  tmp primary control device code
        INTEGER         CSEC    !  tmp secondary control device code
        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP, SCC, SIC  ! tmp fip, scc, sic
        INTEGER         ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER         MXWARN  !  maximum number of warnings
        INTEGER         NPOL    !  number of pollutants in file
        INTEGER, SAVE:: NWARN =0!  number of warnings in this routine
        INTEGER         SS      !  counter for sources
        INTEGER         TPF     !  tmp temporal adjustments setting

        REAL            CEFF    !  tmp control effectiveness
        REAL            EANN    !  tmp annual-ave emission value
        REAL            EMFC    !  tmp emission factor
        REAL            EOZN    !  tmp ozone-season-ave emission value
        REAL            LAT     !  tmp Y-coordinate 
        REAL            LON     !  tmp X-coordinate
        REAL            REFF    !  tmp rule effectiveness
        REAL            DM, HT, FL, TK, VL  ! Temporary stack parms

        LOGICAL, SAVE :: CFLAG              ! true: recalc vel w/ flow & diam
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: 1st time routine called
        LOGICAL, SAVE :: WFLAG    = .FALSE. ! true: all lat-lons to western hemi

        CHARACTER*20    CNTRY   !  country name
        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=BLIDLEN) BLID  ! tmp boiler ID
        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIP code
        CHARACTER(LEN=ORSLEN3) CORS  ! tmp DOE plant ID
        CHARACTER(LEN=IOVLEN3) CPOL  ! tmp pollutant code
        CHARACTER(LEN=DESCLEN) DESC  ! tmp plant description
        CHARACTER(LEN=LINSIZ)  LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp scc

        CHARACTER*16 :: PROGNAME = 'RDIDAPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDIDAPT

        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Set up settings the first time the subroutine is called
        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

C.............  Get settings from the environment
            CFLAG = ENVYN( 'VELOC_RECALC', 
     &                 'Flag for recalculating velocity', .FALSE., IOS )

            WFLAG = ENVYN( 'WEST_HSPHERE',
     &                 'Western hemisphere flag', .TRUE., IOS )

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
        TPF  = MTPRFAC * WKSET
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
		NDROP = NDROP + 1
            END IF

C.............  Check SIC code, warning for missing
            IF( .NOT. CHKINT( LINE( 227:230 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: SIC code is non-' //
     &                 'integer at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( NWARN .LT. MXWARN .AND. 
     &                LINE( 227:230 ) .EQ. ' ' ) THEN
                WRITE( MESG,94010 ) 'WARNING: Missing SIC code at ' //
     &                 'line', IREC, '. Default 0000 will be used.'
                NWARN = NWARN + 1
                CALL M3MESG( MESG )

            END IF

C.............  Initialize start and end positions
            IS = ISINIT - PTOTWIDE  ! array
            IE = IEINIT - PTOTWIDE  ! array

C.............  Check primary and secondary control equipment codes
            DO V = 1, NPOL

                DO K = 6,7
                    IS( K ) = IS( K ) + PTOTWIDE
                    IE( K ) = IE( K ) + PTOTWIDE
                END DO

C.................  Check for bad values, missing is okay
                IF( .NOT. CHKINT( LINE( IS(6):IE(6) ) ) .OR.
     &              .NOT. CHKINT( LINE( IS(7):IE(7) ) )     ) THEN

                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM( V ) )
                    WRITE( MESG,94010 ) 'ERROR: Primary and/or ' //
     &                     'secondary control codes for ' //
     &                     TMPNAM( V )( 1:L ) // ' are non-integer ' //
     &                     'at line', IREC
                    CALL M3MSG2( MESG )

                END IF

            END DO

C.............  Make sure that all of the needed real values are real...

C.............  Stack height, diam, exit temperature, flow, & velocity
            IF( .NOT. CHKREAL( LINE( 120:123 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 124:129 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 130:133 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 134:143 ) )      ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Stack parameters are ' //
     &                 'not numbers or have bad formatting' // CRLF() //
     &                 BLANK10 // 'at line', IREC
                CALL M3MSG2( MESG )

            END IF

C.............  Check stack coordinates, missing is an error
            IF( .NOT. CHKREAL( LINE( 231:239 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 240:248 ) )      ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: latitude and/or ' //
     &                 'longitude are not numbers or have ' // CRLF() //
     &                 BLANK10 // 'bad formatting at line', IREC
                CALL M3MESG( MESG )

            ELSE IF( LINE( 231:239 ) .EQ. ' ' .OR.
     &               LINE( 240:248 ) .EQ. ' '      ) THEN
      
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: latitude and/or ' //
     &                 'longitude are missing at line', IREC
                CALL M3MESG( MESG )

            END IF

C.............  Emissions and associated data
            DO V = 1, NPOL

C.................  Update start and end positions
                DO K = 1, NPTPPOL3
                    IS( K ) = IS( K ) + PTOTWIDE
                    IE( K ) = IE( K ) + PTOTWIDE
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
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields...

            FIP  = ICC * 100000 + 1000 * STR2INT( LINE( 1:2 ) ) +
     &             STR2INT( LINE( 3:5 ) )

            FCID = ADJUSTL( LINE(   6:20  ) )  ! plant ID
            PTID = ADJUSTL( LINE(  21:35  ) )  ! point ID
            SKID = ADJUSTL( LINE(  36:47  ) )  ! stack ID
            CORS = ADJUSTL( LINE(  48:53  ) )  ! DOE plant ID
            BLID = ADJUSTL( LINE(  54:59  ) )  ! boiler ID
            SGID = ADJUSTL( LINE(  60:61  ) )  ! segment ID
            DESC = ADJUSTL( LINE(  62:101 ) )  ! plant description
            TSCC = ADJUSTL( LINE( 102:111 ) )  ! SCC code

            SIC  = MAX( STR2INT ( LINE( 227:230 ) ), 0 )
            HT   = STR2REAL( LINE( 120:123 ) ) * FT2M  ! ft to m
            DM   = STR2REAL( LINE( 124:129 ) ) * FT2M  ! ft to m
            TK   = ( STR2REAL( LINE( 130:133 ) ) - 32 ) * FTOC + CTOK ! F to K
            FL   = STR2REAL( LINE( 134:143 ) ) * FT2M3 ! ft^3/s to m^3/s
            LAT  = STR2REAL( LINE( 231:239 ) )
            LON  = STR2REAL( LINE( 240:248 ) )

            IF( WFLAG .AND. LON .GT. 0 ) LON = -LON

C.............  Recalculate velocity if it is bad or if flag is set
            IF( CFLAG .OR. .NOT. CHKREAL( LINE( 144:152 ) ) ) THEN
                VL = FL / ( 0.25 * PI * DM * DM )
            ELSE
                VL = STR2REAL( LINE( 144:152 ) ) * FT2M   ! ft/s to m/s
            END IF

C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Store source characteristics if dimension is okay
            SS = SS + 1

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

C.............  Initialize start and end positions
            IS = ISINIT - PTOTWIDE  ! array
            IE = IEINIT - PTOTWIDE  ! array

C.............  Loop through pollutants and store data so that there is one
C               record for each pollutant.  This will be consistent with
C               the other reader routines.
            DO V = 1, NPOL

                CPOL = TMPNAM( V )

C.................  Update start and end positions
                DO K = 1, NPTPPOL3
                    IS( K ) = IS( K ) + PTOTWIDE
                    IE( K ) = IE( K ) + PTOTWIDE
                END DO

                EANN = STR2REAL( LINE( IS(1):IE(1) ) )
                EOZN = STR2REAL( LINE( IS(2):IE(2) ) )
                CEFF = STR2REAL( LINE( IS(3):IE(3) ) )
                REFF = STR2REAL( LINE( IS(4):IE(4) ) )
                EMFC = STR2REAL( LINE( IS(5):IE(5) ) )
                CPRI = MAX( IMISS3, STR2INT ( LINE( IS(6):IE(6) ) ))
                CSEC = MAX( IMISS3, STR2INT ( LINE( IS(7):IE(7) ) ))

                IF( NWARN .LT. MXWARN .AND. 
     &              EANN  .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing annual ' //
     &                 'emissions for at line', IREC, 'for ' // CPOL
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
                END IF

                IF( NWARN .LT. MXWARN .AND. 
     &              EOZN  .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing ozone ' //
     &                 'season emissions at line', IREC, 'for ' // CPOL
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
                END IF

                IF( NWARN .LT. MXWARN .AND. 
     &              CEFF  .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing control ' //
     &                     'efficiency at line', IREC, 'for ' // CPOL
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
                END IF
                
                IF( NWARN .LT. MXWARN .AND. 
     &              REFF .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                     'effectiveness at line', IREC, 'for '//CPOL
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
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
 
                    CALL BLDCSRC( CFIP, FCID, PTID, SKID, SGID, 
     &                            TSCC, CHRBLNK3, CCOD, CSOURCA( ES ) )

                END IF  !  if ES in range

            END DO      !  end of loop through pollutants

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

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDIDAPT
