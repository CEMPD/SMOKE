        SUBROUTINE RDIDAAR( FDEV, NRAWIN, NRAWBP, MXIPOL, 
     &                      INVPNAM, NRAWOUT, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line XXX
C
C  DESCRIPTION:
C      This subroutine reads the IDA format area-source inventory
C      files. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      created by M. Houyoux (04/99) 
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CONST3.EXT'    !  physical constants
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

        EXTERNAL    CRLF, ENVINT, ENVYN, GETFLINE, GETNLIST, INDEX1, 
     &              STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file
        INTEGER     , INTENT (IN) :: NRAWIN ! total raw record-count 
        INTEGER     , INTENT (IN) :: NRAWBP ! total raw record times pols
        INTEGER     , INTENT (IN) :: MXIPOL ! max no of inventory pols
        CHARACTER(*), INTENT (IN) :: INVPNAM( MXIPOL ) ! inv pol names
        INTEGER     , INTENT(OUT) :: NRAWOUT! outgoing source * pollutants
        LOGICAL     , INTENT(OUT) :: EFLAG  ! outgoing error flag
        INTEGER     , INTENT(OUT) :: NDROP  !  number of records dropped
        REAL        , INTENT(OUT) :: EDROP( MXIPOL )  ! emis dropped per pol

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

C...........   Local allocatable arrays
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: TMPNAM( : )! pol names for file
        INTEGER               , ALLOCATABLE :: POLPOS( : )! pol pos in INVPNAM

C...........   Local arrays
        INTEGER          IS( NARPPOL3 )  ! start position for each pol char
        INTEGER          IE( NARPPOL3 )  ! end position for each pol char

C...........   Counters of total number of input records
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count
        INTEGER, SAVE :: NSRCPOL = 0 ! cumulative source x pollutants count

C...........   Other local variables
        INTEGER         I, J, K, L, V  ! counters and indices

        INTEGER         COD     !  tmp pollutant position in INVPNAM
        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP     !  tmp FIPS code
        INTEGER         ICC     !  position of CNTRY in CNTRYNM3
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
        REAL            REFF    !  tmp rule effectiveness
        REAL            RPEN    !  tmp rule penetration

        CHARACTER*20    CNTRY   !  country name
        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVPNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIP code
        CHARACTER(LEN=IOVLEN3) CPOL  ! tmp pollutant code
        CHARACTER(LEN=LINSIZ)  LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp scc

        CHARACTER*16 :: PROGNAME = 'RDIDAAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDIDAAR

        MXWARN = ENVINT( WARNSET, ' ', 99999999, IOS )

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
        TPF  = MTPRFAC * WTPRFAC
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

C.............  Scan for header lines
            IF( LINE( 1:1 ) .EQ. '#' ) THEN

                IF ( LINE(2:8) .EQ. 'COUNTRY' ) THEN  ! read in country-name
                    CNTRY = ADJUSTL( LINE( 9:L ) )
                    ICC   = INDEX1( CNTRY, MXCNTRY3, CNTRYNM3 )

                    IF ( ICC .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( CNTRY )
                        WRITE( MESG, 94010 )
     &                         'Unknown country name "' // CNTRY( 1:L )
     &                         // '" encountered at line', IREC,
     &                         '. Must edit EMCNST3.EXT and recompile.'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF

                    ICC   = CNTRYCD3( ICC )
          
                ELSEIF ( LINE(2:5) .EQ. 'YEAR' ) THEN ! read in inventory year
                    INY = STR2INT( LINE( 6:L ) )
                    IF ( INY .LT. 1971 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 ) 'Invalid year ', INY, 
     &                         'encountered at line ', IREC
                        CALL M3MESG( MESG )
                    END IF

                ELSEIF ( LINE(2:6) .EQ. 'POLID' ) THEN ! read in pollutants

C.....................  Allocate memory for current file for reading pol names
C                       and storing positions in master list
                    IF( .NOT. ALLOCATED( TMPNAM ) ) THEN
                        LINE = LINE( 7:L )
                        L = LEN_TRIM( LINE )
                        NPOL = GETNLIST( L, LINE )
                        ALLOCATE( TMPNAM( NPOL ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'TMPNAM', PROGNAME )
                        ALLOCATE( POLPOS( NPOL ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'POLPOS', PROGNAME )
                    END IF

                    IF( NPOL .GT. MXPOLFIL ) THEN
                        WRITE( MESG,94010 ) 'INTERNAL ERROR: Maximum '//
     &                         'pollutants allowed (MXPOLFIL) is', 
     &                         MXPOLFIL, 'but file has', NPOL
                        CALL M3MSG2( MESG )
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

C.....................  Parse the header line into the pollutant names
                    CALL PARSLINE( LINE, NPOL, TMPNAM )

C.....................  Store the position in master list of each pollutant
C.....................  Write error if pollutant is not found.
                    DO V = 1, NPOL

                        CPOL = TMPNAM( V )
                        L = LEN_TRIM( CPOL )
                        COD = INDEX1( CPOL, MXIPOL, INVPNAM )
                        IF( COD .LE. 0 ) THEN
                            EFLAG = .TRUE.
                            MESG = 'ERROR: Pollutant "' // CPOL( 1:L )//
     &                             '" not in master pollutant list!'
                            CALL M3MSG2( MESG )
                        ELSE
                            POLPOS( V ) = COD
                        END IF

                    END DO

                END IF

                CYCLE   ! to next iteration

            END IF

C.............  If the line is not a header line, make sure that all of the 
C               important header lines have been read in...

C.............  Check for country header
            IF( ICC .LT. 0 ) THEN
                EFLAG = .TRUE.
                ICC = 9         ! to turn off error message
                MESG = 'ERROR: Country name was not set with ' //
     &                 '#COUNTRY header before first data line.'
                CALL M3MSG2( MESG )
            END IF

C.............  Check for inventory year header
            IF( INY .EQ. 0 ) THEN
                EFLAG = .TRUE.
                INY = 1       ! to turn off error message
                MESG = 'ERROR: Inventory year was not set with ' //
     &                 '#YEAR header before first data line.'
                CALL M3MSG2( MESG )
            END IF

C.............  Check for pollutant names header
            IF( NPOL .EQ. 0 ) THEN
                EFLAG = .TRUE. 
                NPOL = -1       ! to turn off error message
                MESG = 'ERROR: Pollutants were not set with ' //
     &                 '#POLID header before first data line.'
                CALL M3MESG( MESG )
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

                END IF

                IF( LINE( IS(1):IE(1) ) .EQ. ' ' .AND.
     &              LINE( IS(2):IE(2) ) .EQ. ' '       ) THEN

                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM( V ) )
                    WRITE( MESG,94010 ) 'ERROR: All emissions ' //
     &                     'data for ' // TMPNAM( V )( 1:L ) //  
     &                     ' are missing at line', IREC
                    CALL M3MESG( MESG )

                END IF

            END DO

C.............  If there has been an error, do not try to store any of the
C               records.  Instead  go to next line of file.
            IF( EFLAG ) CYCLE
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields...

            FIP  = ICC * 100000 + 1000 * STR2INT( LINE( 1:2 ) ) +
     &             STR2INT( LINE( 3:5 ) )

            TSCC = LINE( 6:15 )    ! SCC code

C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Store source characteristics if dimension is okay
            SS = SS + 1

            IF( SS .LE. NRAWIN ) THEN

                IFIPA  ( SS ) = FIP
                TPFLGA ( SS ) = TPF
                INVYRA ( SS ) = INY
                CSCCA  ( SS ) = TSCC
 
            END IF 

C.............  Initialize start and end positions
            IS = ISINIT - AROTWIDE  ! array
            IE = IEINIT - AROTWIDE  ! array

C.............  Loop through pollutants and store data so that there is one
C               record for each pollutant.  This will be consistent with
C               the other reader routines.
            DO V = 1, NPOL

C.................  Update start and end positions
                DO K = 1, NARPPOL3
                    IS( K ) = IS( K ) + AROTWIDE
                    IE( K ) = IE( K ) + AROTWIDE
                END DO

                EANN = STR2REAL( LINE( IS(1):IE(1) ) )
                EOZN = STR2REAL( LINE( IS(2):IE(2) ) )
                EMFC = STR2REAL( LINE( IS(3):IE(3) ) )
                CEFF = STR2REAL( LINE( IS(4):IE(4) ) )
                REFF = STR2REAL( LINE( IS(5):IE(5) ) )
                RPEN = STR2REAL( LINE( IS(6):IE(6) ) )

                IF( NWARN .LT. MXWARN .AND. EANN .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing annual ' //
     &                 'emissions at line', IREC
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
                END IF

                IF( NWARN .LT. MXWARN  .AND. EOZN .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing ozone ' //
     &                 'season emissions at line', IREC
                    NWARN = NWARN + 1
                    CALL M3MESG( MESG )
                END IF

                IF( CEFF .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing control ' //
     &                 'efficiency at line', IREC, '.' // CRLF() // 
     &                 BLANK10 // 'Default of 100. will be used.'
                    IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                    CEFF = 100.
                END IF
                
                IF( REFF .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                 'effectiveness at line', IREC, '.' // CRLF() // 
     &                 BLANK10 // 'Default of 100. will be used.'
                    IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                    REFF = 100.
                END IF

                IF( RPEN .LT. AMISS3 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                 'pentration at line', IREC, '.' // CRLF() // 
     &                 BLANK10 // 'Default of 100. will be used.'
                    IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                    RPEN = 100.
                END IF

C.................  Store data in final arrays if there is enough memory
                ES = ES + 1

                IF ( ES .LE. NRAWBP ) THEN

                    INDEXA ( ES     ) = ES                    
                    POLVLA ( ES,NEM ) = EANN
                    POLVLA ( ES,NOZ ) = EOZN
                    POLVLA ( ES,NEF ) = EMFC
                    POLVLA ( ES,NCE ) = CEFF
                    POLVLA ( ES,NRE ) = REFF
                    POLVLA ( ES,NRP ) = RPEN

                    WRITE( CCOD,94125 ) POLPOS( V )

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

        DEALLOCATE( TMPNAM, POLPOS )

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDIDAAR
