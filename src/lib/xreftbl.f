
        SUBROUTINE PXREFTBL( OPTYPE, NXREF, MXCHRS, NCHARS, JSCC, 
     &                       NIPOL, EINAM )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine will do the logic for grouping the point sources cross 
C      references.  In the logic phase, it will merely store the group number 
C      and the count in the group for each x-ref entry. Then, depending on
C      the type of operation (temporal or speciation), it will allocate memory 
C      for an populate the appropriate grouped tables.  The unsorted x-ref
C      data and grouped tables will be passed by a module for cross-reference
C      data.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   This module is for cross reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         STR2INT

        EXTERNAL   CRLF, STR2INT 

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: OPTYPE ! operation type (tmprl or spec)
        INTEGER     , INTENT (IN) :: NXREF  ! no. unsorted x-ref entries
        INTEGER     , INTENT (IN) :: MXCHRS ! maximum number of source chars
        INTEGER     , INTENT (IN) :: NCHARS ! actual number of source chars
        INTEGER     , INTENT (IN) :: JSCC   ! pos of SCC in src def'n if ne 0
        INTEGER     , INTENT (IN) :: NIPOL  ! number of inventory pollutants
        CHARACTER(*), INTENT (IN) :: EINAM( NIPOL ) ! pollutant name

C...........   Local parameters
        INTEGER, PARAMETER :: SNFLEN3 = SRCLEN3 - FIPLEN3

C...........   Arrays for intial pass through x-ref to determine degree of 
C              each record and store the type and count of that type
        INTEGER         XTYPE( NXREF )  ! 
        INTEGER         XTCNT( NXREF )  !  

C...........  Arrays for counting number of x-ref records in each matching
C             degree and comparing current record with previous record as part
C             of the counting. Note that although all arrays below are 
C             available for each degree of matching, only the array elements 
C             that are appropriate for a given degree are actually populated. 
        INTEGER                N   ( 0:NXTYPES )  ! cnt for degree of matching
        INTEGER                PIFIP ( NXTYPES )  ! previous co/st/cy code

        CHARACTER(LEN=SCCLEN3) PTSCC ( NXTYPES )  ! previous SCC
        CHARACTER(LEN=SRCLEN3) PCSRC ( NXTYPES )  ! previous CSRC
        CHARACTER(LEN=SS5LEN3) PCSSC ( NXTYPES )  ! previous CSRC(part) // SCC

C...........   Array of source characeristics
        CHARACTER*300           CHARS( MXCHRS )

C...........   Other local variables
        INTEGER       I, J, J1, J2, K, L, T     ! counter and indices

        INTEGER       ICYID            ! temporary county code
        INTEGER       IFIP             ! temporary FIPS code
        INTEGER       IMON, IWEK, IDIU ! temporary temporal profile codes
        INTEGER       IDIUPS           ! tmp w/ pollutant-specific indicator
        INTEGER       IOS              ! i/o status
        INTEGER       ISP              ! temporary pollutant position in EINAM
        INTEGER       LOPT             ! length of OPTYPE
        INTEGER       LSA              ! in SCC, position at end of 1st part
        INTEGER       LSB              ! in SCC, position at start of 2nd part
        INTEGER       NT               ! code for specificity of x-ref entry
        INTEGER    :: PISP = IMISS3    ! previous iteration ISP

        LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred
        LOGICAL    :: DEFAULT( NIPOL ) ! true: if default entry in x-ref
        LOGICAL    :: POLDFLT          ! true: okay to have pol-spec defaults
        LOGICAL    :: SFLAG = .FALSE.  ! true: operation type is speciation
        LOGICAL    :: TFLAG = .FALSE.  ! true: operation type is temporal

        CHARACTER*300          BUFFER        ! source definition buffer
        CHARACTER*300          MESG          ! message buffer

        CHARACTER(LEN=STALEN3) CSTA          ! temporary (character) state code
        CHARACTER(LEN=SCLLEN3) SCCL5         ! left digits of TSCC
        CHARACTER(LEN=SCRLEN3) SCCR5         ! right 5 digits of TSCC
        CHARACTER(LEN=SCRLEN3) SCRZERO       ! buf for 0 right 5 digits of TSCC
        CHARACTER(LEN=SNFLEN3) CNFIP         ! characterstics without FIPS code
        CHARACTER(LEN=SRCLEN3) CSRC          ! temporary source char string
        CHARACTER(LEN=FIPLEN3) CFIP          ! temporary (character) FIPS code
        CHARACTER(LEN=FIPLEN3) FIPZERO       ! buffer for zero FIPS code
        CHARACTER(LEN=SCCLEN3) PSCC          ! previous SCC
        CHARACTER(LEN=SCCLEN3) TSCC          ! temporary SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO       ! buffer for zero SCC
        CHARACTER(LEN=SS5LEN3) CSRCSCC       ! buffer for source // SCC
        CHARACTER(LEN=SPNLEN3) SPCODE        ! tmp for speciation profile code

        CHARACTER*16 :: PROGNAME = 'PXREFTBL' ! program name

C***********************************************************************
C   begin body of subroutine PXREFTBL

        LOPT = LEN_TRIM( OPTYPE )

C.........  Check for valid operation type
        SELECT CASE( OPTYPE )
        CASE( 'SPECIATION' )
            POLDFLT = .TRUE.
            SFLAG   = .TRUE.
        CASE( 'TEMPORAL' ) 
            POLDFLT = .FALSE.
            TFLAG   = .TRUE.
        CASE DEFAULT

            MESG = 'INTERNAL ERROR: Operation type "' // 
     &             OPTYPE( 1:LOPT ) //
     &             '" not known in subroutine ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT
        
C.........  Set up zero strings for FIPS code of zero and SCC code of zero
        FIPZERO = REPEAT( '0', FIPLEN3 )
        SCCZERO = REPEAT( '0', SCCLEN3 )
        SCRZERO = REPEAT( '0', SCRLEN3 )

C.........  Initialize default array
        DEFAULT = .FALSE.   ! array

C.........  Initialize arrays for counting number of x-ref records in each
C           degree of matching
        N      = 0   ! arrays
        PIFIP  = 0
        PSCC   = ' '
        PCSRC  = ' '
        PCSSC  = ' '

C.........  Set up for parsing SCCs
        LSB = SCCLEN3 - SCRLEN3 + 1
        LSA = LSB - 1

C.........  Initialize source characteristics
        CHARS = ' '  ! array

C.........  Loop through and count entries of each type. Store type.
C.........  For CSRC, don't include pollutant for grouping.
        DO I = 1, NXREF

            J = INDXTA( I )

            CSRC    = CSRCTA( J )( 1:SRCLEN3 )
            ISP     = ISPTA ( J )
            TSCC    = CSCCTA( J )

            DO J = 1, NCHARS
                CHARS( J ) = CSRC( PTBEGL3( J ):PTENDL3( J ) )
            ENDDO

C.............  Rearrange CHARS if SCC is a part of the source definition
C               because we have now stored SCC separately.  It will be
C               much easier to group the source characteristics this way.
C               Still go back to using original definition when storing
C               in tables (i.e., use CSRC)
            IF( JSCC .GT. 0 ) THEN
                DO J = JSCC, NCHARS - 1
                    CHARS( J ) = CHARS( J + 1 )
                ENDDO
                CHARS( NCHARS ) = ' '
            ENDIF

C.............  Set up partial strings for checking
            CFIP    = CHARS( 1 )
            IFIP    = STR2INT( CFIP )         ! For checking previous
            CSTA    = CFIP( 1:STALEN3 )
            ICYID   = IFIP - STR2INT( CSTA ) * 1000
            SCCL5   = TSCC(   1:LSA     )
            SCCR5   = TSCC( LSB:SCCLEN3 )
            CSRCSCC = CSRC // TSCC
            CNFIP   = CSRC( PLTPOS3:SRCLEN3 )
  
C.............  Select cases
C.............  Note that since these are sorted in order of increasing FIPS
C               code, SCC, pollutant index, etc., that the entries with zero for
C               these characteristics will appear earlier in the sorted list
            IF( CHARS( 1 ) .EQ. FIPZERO ) THEN           ! FIPS code is default

                IF( TSCC .EQ. SCCZERO ) THEN                   ! SCC is default

                    IF( POLDFLT ) THEN           ! Pollutant-specific permitted

                        IF( ISP .EQ. 0 ) THEN
                            NT = 1

                        ELSEIF( ISP .NE. 0 .AND. .NOT.    ! Pollutant specified
     &                          DEFAULT( ISP )         ) THEN

                            DEFAULT( ISP ) = .TRUE.
                            NT = 1

                        ELSE                                  ! Report and skip
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE                     ! Pollutant-specific not permitted

                        IF( ISP .EQ. 0 .AND. .NOT.    ! Pollutant not specified
     &                      DEFAULT( 1 )           ) THEN 

                            DEFAULT( 1 ) = .TRUE.
                            NT = 1

                        ELSEIF( ISP .NE. 0 ) THEN             ! Report and skip
                            MESG = 'Cannot use pollutant-specific ' //
     &                             'ultimate-default'
                            CALL REPORT_INVALID_XREF( MESG )
                            NT = 0

                        ELSEIF( DEFAULT( 1 ) ) THEN           ! Report and skip
                            CALL REPORT_DUP_XREF
                            NT = 0
                        ENDIF

                    ENDIF

                ELSEIF( SCCR5 .EQ. SCRZERO ) THEN        ! 5 or 3-digit SCC

                    NT = 2
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    ENDIF

                ELSE                                         ! Complete SCC

                    NT = 3
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    ENDIF

                ENDIF

            ELSEIF( ICYID .EQ. 0 ) THEN            ! County code is default

                IF( TSCC .EQ. SCCZERO ) THEN         ! SCC code is default

                    IF( ISP .EQ. 0 ) THEN         ! Pollutant not specified
                        NT = 4
                        IF( IFIP .NE. PIFIP( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP( NT ) = IFIP

                        ELSE
                            CALL REPORT_DUP_XREF
                            NT = 0
                        ENDIF

                    ELSE                                  ! Report and skip
                        MESG = 'Cannot use pollutant-specific ' //
     &                         'Country/State-default'
                        CALL REPORT_INVALID_XREF( MESG )
                        NT = 0

                    ENDIF

                ELSEIF( SCCR5 .EQ. SCRZERO ) THEN        ! 5 or 3-digit SCC

                    NT = 5
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC
                       
                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    ENDIF

                ELSE                                         ! Complete SCC

                    NT = 6
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    ENDIF

                ENDIF

            ELSEIF( CNFIP .EQ. ' ' ) THEN    ! Country/St/Co code is complete

                IF( TSCC .EQ. SCCZERO ) THEN            ! SCC code is default

                    IF( ISP .EQ. 0 ) THEN !   & non-FIP field blank
                        NT = 7
                        IF( IFIP .NE. PIFIP( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP( NT ) = IFIP

                        ELSE
                            CALL REPORT_DUP_XREF
                            NT = 0
                        ENDIF

                    ELSE                                  ! Report and skip
                        MESG = 'Cannot use pollutant-specific ' //
     &                         'Country/State/County-default'
                        CALL REPORT_INVALID_XREF( MESG )
                        NT = 0

                   ENDIF

                ELSEIF( SCCR5 .EQ. SCRZERO ) THEN        ! 5 or 3-digit SCC

                    NT = 8
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    ENDIF

                ELSE                                         ! Complete SCC

                    NT = 9
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    ENDIF

                ENDIF                                              ! End SCC 

            ELSE                                        ! Plant is specified

                IF( TSCC .EQ. SCCZERO ) THEN           ! SCC code is default

C.....................  Loop through plant-specific characteristics. Only the
C                       plant is permitted to not have an SCC not specified
                    NT = 15
                    DO J = MXCHRS, 2, -1

                        IF( NT .EQ. 10 .AND. CHARS( J ) .NE. ' ' ) THEN

                            IF( CSRC .NE. PCSRC( NT ) ) THEN
                                N( NT ) = N( NT ) + 1
                                PCSRC( NT ) = CSRC
                                EXIT                      ! End loop with NT

                            ELSEIF( ISP .NE. PISP ) THEN
                                EXIT                      ! End loop with NT

                            ELSE
                                CALL REPORT_DUP_XREF
                                NT = 0
                                EXIT                      ! End loop with NT
                           ENDIF

                        ELSEIF( NT         .GT. 10  .AND. 
     &                          CHARS( J ) .NE. ' '       ) THEN

                            MESG = 'SCC is not specified'
                            CALL REPORT_INVALID_XREF( MESG )
                            NT = 0
                            EXIT                      ! End loop with NT

                        ELSEIF( NT .EQ. 10 ) THEN
                            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )

                            MESG = 'INTERNAL ERROR: Check PXREFTBL ' //
     &                             'for processing record: ' //
     &                             CRLF() // BLANK10 // BUFFER( 1:L ) //
     &                             ' POL:' // EINAM( ISP )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                        ENDIF

                        NT = NT - 1

                    ENDDO           ! End loop on plant characteristics

                ELSEIF( SCCR5 .EQ. SCRZERO ) THEN         ! 5 or 3-digit SCC

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )

                    MESG = 'Partial SCC "' // TSCC // '" is given ' //
     &                     'instead of full SCC'
                    CALL REPORT_INVALID_XREF( MESG )
                    NT = 0

                ELSE                                          ! Complete SCC
C.....................  Loop through plant-specific characteristics,
C                       and store the most specific entries first.
C.....................  Process NT 16 through 11
                    NT = 16
                    DO J = MXCHRS, 2, -1

                        IF( CHARS( J ) .NE. ' ' ) THEN

                            IF( CSRCSCC .NE. PCSSC( NT ) ) THEN
                                N( NT ) = N( NT ) + 1
                                PCSSC( NT ) = CSRCSCC
                                EXIT                      ! End loop with NT

                            ELSEIF( ISP .NE. PISP ) THEN
                                EXIT                      ! End loop with NT

                            ELSE
                                CALL REPORT_DUP_XREF
                                NT = 0
                                EXIT                      ! End loop with NT
                            ENDIF

                        ENDIF

                        NT = NT - 1

                    ENDDO           ! End loop on plant characteristics

                ENDIF ! End SCC

            ENDIF ! End degree of Country/State/County code, or plant specified

            PISP = ISP

            XTYPE( I ) = NT
            XTCNT( I ) = N( NT )         ! Dimensioned from 0 so NT can = 0

        ENDDO                            ! End Loop on sorted x-ref entries

C.........  Allocate memory for the tables, depending on the operation type
C.........  Note that types 4 and 7 are not pollutant-specific, so there is
C           no NIPOL dimension for those tables.
C.........  Also, initialize profile codes
        IF( TFLAG ) THEN
            ALLOCATE( MPRT01( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MPRT01', PROGNAME )
            ALLOCATE( WPRT01( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WPRT01', PROGNAME )
            ALLOCATE( DPRT01( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DPRT01', PROGNAME )
            MPRT01 = IMISS3 ! arrays
            WPRT01 = IMISS3
            DPRT01 = IMISS3
        ELSEIF( SFLAG ) THEN
            ALLOCATE( CSPT01( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSPT01', PROGNAME )
            CSPT01 = EMCMISS3
        ENDIF

        J = N( 2 )                                       ! SCC=left, FIP=0
        IF( J .GT. 0 ) THEN
            ALLOCATE( CHRT02( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT02', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT02( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT02', PROGNAME )
                ALLOCATE( WPRT02( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT02', PROGNAME )
                ALLOCATE( DPRT02( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT02', PROGNAME )
                MPRT02 = IMISS3 ! arrays
                WPRT02 = IMISS3
                DPRT02 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT02( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT02', PROGNAME )
                CSPT02 = EMCMISS3

            ENDIF
        ENDIF

        J = N( 3 )                                   ! SCC=all, FIP=0
        IF( J .GT. 0 ) THEN
            ALLOCATE( CHRT03( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT03', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT03( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT03', PROGNAME )
                ALLOCATE( WPRT03( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT03', PROGNAME )
                ALLOCATE( DPRT03( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT03', PROGNAME )
                MPRT03 = IMISS3 ! arrays
                WPRT03 = IMISS3
                DPRT03 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT03( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT03', PROGNAME )
                CSPT03 = EMCMISS3

            ENDIF
        ENDIF
                
        J = N( 4 )                                 ! SCC=0, FIP=state
        IF( J .GT. 0 ) THEN
            ALLOCATE( CHRT04( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT04', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT04( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT04', PROGNAME )
                ALLOCATE( WPRT04( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT04', PROGNAME )
                ALLOCATE( DPRT04( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT04', PROGNAME )
                MPRT04 = IMISS3 ! arrays
                WPRT04 = IMISS3
                DPRT04 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT04( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT04', PROGNAME )
                CSPT04 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 5 )                                 ! SCC=left, FIP=state
        IF( J .GT. 0 ) THEN
            ALLOCATE( CHRT05( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT05', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT05( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT05', PROGNAME )
                ALLOCATE( WPRT05( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT05', PROGNAME )
                ALLOCATE( DPRT05( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT05', PROGNAME )
                MPRT05 = IMISS3 ! arrays
                WPRT05 = IMISS3
                DPRT05 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT05( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT05', PROGNAME )
                CSPT05 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 6 )  
        IF( J .GT. 0 ) THEN                        ! SCC=all, FIP=state
            ALLOCATE( CHRT06( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT06', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT06( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT06', PROGNAME )
                ALLOCATE( WPRT06( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT06', PROGNAME )
                ALLOCATE( DPRT06( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT06', PROGNAME )
                MPRT06 = IMISS3 ! arrays
                WPRT06 = IMISS3
                DPRT06 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT06( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT06', PROGNAME )
                CSPT06 = EMCMISS3

            ENDIF
        ENDIF
                        
        J = N( 7 )   
        IF( J .GT. 0 ) THEN                          ! SCC=0, FIP=all
            ALLOCATE( CHRT07( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT07', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT07( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT07', PROGNAME )
                ALLOCATE( WPRT07( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT07', PROGNAME )
                ALLOCATE( DPRT07( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT07', PROGNAME )
                MPRT07 = IMISS3 ! arrays
                WPRT07 = IMISS3
                DPRT07 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT07( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT07', PROGNAME )
                CSPT07 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 8 )
        IF( J .GT. 0 ) THEN                         ! SCC=left, FIP=all
            ALLOCATE( CHRT08( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT08', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT08( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT08', PROGNAME )
                ALLOCATE( WPRT08( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT08', PROGNAME )
                ALLOCATE( DPRT08( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT08', PROGNAME )
                MPRT08 = IMISS3 ! arrays
                WPRT08 = IMISS3
                DPRT08 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT08( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT08', PROGNAME )
                CSPT08 = EMCMISS3

            ENDIF
        ENDIF
                        
        J = N( 9 )
        IF( J .GT. 0 ) THEN                          ! SCC=all, FIP=all
            ALLOCATE( CHRT09( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT09', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT09( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT09', PROGNAME )
                ALLOCATE( WPRT09( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT09', PROGNAME )
                ALLOCATE( DPRT09( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT09', PROGNAME )
                MPRT09 = IMISS3 ! arrays
                WPRT09 = IMISS3
                DPRT09 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT09( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT09', PROGNAME )
                CSPT09 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 10 )
        IF( J .GT. 0 ) THEN                       ! PLANT=non-blank, SCC=0
            ALLOCATE( CHRT10( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT10', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT10( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT10', PROGNAME )
                ALLOCATE( WPRT10( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT10', PROGNAME )
                ALLOCATE( DPRT10( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT10', PROGNAME )

                MPRT10 = IMISS3 ! arrays
                WPRT10 = IMISS3
                DPRT10 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT10( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT10', PROGNAME )
                CSPT10 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 11 )         
        IF( J .GT. 0 ) THEN                      ! PLANT=non-blank, SCC=all
            ALLOCATE( CHRT11( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT11', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT11( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT11', PROGNAME )
                ALLOCATE( WPRT11( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT11', PROGNAME )
                ALLOCATE( DPRT11( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT11', PROGNAME )
                MPRT11 = IMISS3 ! arrays
                WPRT11 = IMISS3
                DPRT11 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT11( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT11', PROGNAME )
                CSPT11 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 12 )        
        IF( J .GT. 0 ) THEN                      ! CHAR1=non-blank, SCC=all
            ALLOCATE( CHRT12( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT12', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT12( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT12', PROGNAME )
                ALLOCATE( WPRT12( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT12', PROGNAME )
                ALLOCATE( DPRT12( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT12', PROGNAME )
                MPRT12 = IMISS3 ! arrays
                WPRT12 = IMISS3
                DPRT12 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT12( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT12', PROGNAME )
                CSPT12 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 13 )  
        IF( J .GT. 0 ) THEN                      ! CHAR2=non-blank, SCC=all
            ALLOCATE( CHRT13( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT13', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT13( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT13', PROGNAME )
                ALLOCATE( WPRT13( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT13', PROGNAME )
                ALLOCATE( DPRT13( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT13', PROGNAME )
                MPRT13 = IMISS3 ! arrays
                WPRT13 = IMISS3
                DPRT13 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT13( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT13', PROGNAME )
                CSPT13 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 14 )
        IF( J .GT. 0 ) THEN                      ! CHAR3=non-blank, SCC=all
            ALLOCATE( CHRT14( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT14', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT14( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT14', PROGNAME )
                ALLOCATE( WPRT14( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT14', PROGNAME )
                ALLOCATE( DPRT14( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT14', PROGNAME )
                MPRT14 = IMISS3 ! arrays
                WPRT14 = IMISS3
                DPRT14 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT14( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT14', PROGNAME )
                CSPT14 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 15 )
        IF( J .GT. 0 ) THEN                      ! CHAR4=non-blank, SCC=all
            ALLOCATE( CHRT15( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT15', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT15( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT15', PROGNAME )
                ALLOCATE( WPRT15( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT15', PROGNAME )
                ALLOCATE( DPRT15( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT15', PROGNAME )
                MPRT15 = IMISS3 ! arrays
                WPRT15 = IMISS3
                DPRT15 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT15( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT15', PROGNAME )
                CSPT15 = EMCMISS3

            ENDIF
        ENDIF
            
        J = N( 16 )
        IF( J .GT. 0 ) THEN                      ! CHAR5=non-blank, SCC=all
            ALLOCATE( CHRT16( J ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT16', PROGNAME )

            IF( TFLAG ) THEN
                ALLOCATE( MPRT16( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRT16', PROGNAME )
                ALLOCATE( WPRT16( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'WPRT16', PROGNAME )
                ALLOCATE( DPRT16( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DPRT16', PROGNAME )
                MPRT16 = IMISS3 ! arrays
                WPRT16 = IMISS3
                DPRT16 = IMISS3

            ELSEIF( SFLAG ) THEN
                ALLOCATE( CSPT16( J,NIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSPT16', PROGNAME )
                CSPT16 = EMCMISS3

            ENDIF
        ENDIF
            
C.........  Now store the tables, depending on the operation type
C.........  First store the characteristics for each profile, then store the
C           operation-specific profile numbers.
        DO I = 1, NXREF
            J      = INDXTA( I )
            CSRC   = CSRCTA( J )
            ISP    = ISPTA ( J )
            TSCC   = CSCCTA( J )

C.............  Set up partial strings for saving
            SCCL5  = TSCC( 1:LSA )
            CFIP   = CSRC( 1:FIPLEN3 ) 
            CSTA   = CSRC( 1:STALEN3 )

            T      = XTYPE ( I )
            K      = XTCNT ( I )

            SELECT CASE ( T )

            CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate
            CASE( 1 )  ! Skip this because no source-characteristics
            CASE( 2 )
                CHRT02( K ) = SCCL5
            CASE( 3 )
                CHRT03( K ) = TSCC
            CASE( 4 )                   ! NOTE:  pol-specific was excluded
                CHRT04( K ) = CSTA
            CASE( 5 )
                CHRT05( K ) = CSTA // SCCL5
            CASE( 6 )
                CHRT06( K ) = CSTA // TSCC
            CASE( 7 )                   ! NOTE:  pol-specific was excluded
                CHRT07( K ) = CFIP
            CASE( 8 )
                CHRT08( K ) = CFIP // SCCL5
            CASE( 9 )
                CHRT09( K ) = CFIP // TSCC
            CASE( 10 )
                CHRT10( K ) = CSRC( 1:PTENDL3( 2 ) )
            CASE( 11 )
                CHRT11( K ) = CSRC( 1:PTENDL3( 2 ) ) // TSCC
            CASE( 12 )
                CHRT12( K ) = CSRC( 1:PTENDL3( 3 ) ) // TSCC
            CASE( 13 )
                CHRT13( K ) = CSRC( 1:PTENDL3( 4 ) ) // TSCC
            CASE( 14 )
                CHRT14( K ) = CSRC( 1:PTENDL3( 5 ) ) // TSCC
            CASE( 15 )
                CHRT15( K ) = CSRC( 1:PTENDL3( 6 ) ) // TSCC
            CASE( 16 )
                CHRT16( K ) = CSRC( 1:PTENDL3( 7 ) ) // TSCC
            CASE DEFAULT

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'INTERNAL ERROR: Point source cross-' // 
     &                 'reference category', T, 
     &                 'not known in subroutine ' // PROGNAME
                CALL M3MESG( MESG )

            END SELECT

        ENDDO                            ! End Loop on sorted x-ref entries

        IF( TFLAG ) THEN

            DO I = 1, NXREF

                J      = INDXTA( I )
                ISP    = ISPTA ( J )
                IMON   = MPRNA ( J )
                IWEK   = WPRNA ( J )
                IDIU   = DPRNA ( J )
                IDIUPS = IDIU + ADDPS  ! for pollutant specific

                T      = XTYPE ( I )
                K      = XTCNT ( I )
C
C.................  Populate tables depending on type. Note that the pollutant-
C                   specific entries are assumed to always come after the
C                   non-specific ones (based on the previous sorting).
C.................  The temporal pol-specific entries are stored by adding 90000
C                   to the monthly profile number (which has a maximum of 3 
C                   digits) so that the pol-specific can be identified later
                SELECT CASE ( T )

                CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate

                CASE( 1 )
                    MPRT01 = IMON  ! Arrays
                    WPRT01 = IWEK
                    DPRT01 = IDIU

                CASE( 2 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT02( K,: ) = IMON
                        WPRT02( K,: ) = IWEK
                        DPRT02( K,: ) = IDIU

                    ELSE
                        MPRT02( K,ISP ) = IMON
                        WPRT02( K,ISP ) = IWEK
                        DPRT02( K,ISP ) = IDIUPS

                    ENDIF

                CASE( 3 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT03( K,: ) = IMON
                        WPRT03( K,: ) = IWEK
                        DPRT03( K,: ) = IDIU

                    ELSE
                        MPRT03( K,ISP ) = IMON
                        WPRT03( K,ISP ) = IWEK
                        DPRT03( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 4 )                   ! NOTE:  pol-specific was excluded
                    MPRT04( K ) = IMON
                    WPRT04( K ) = IWEK
                    DPRT04( K ) = IDIU

                CASE( 5 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT05( K,: ) = IMON
                        WPRT05( K,: ) = IWEK
                        DPRT05( K,: ) = IDIU

                    ELSE
                        MPRT05( K,ISP ) = IMON
                        WPRT05( K,ISP ) = IWEK
                        DPRT05( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 6 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT06( K,: ) = IMON
                        WPRT06( K,: ) = IWEK
                        DPRT06( K,: ) = IDIU

                    ELSE
                        MPRT06( K,ISP ) = IMON
                        WPRT06( K,ISP ) = IWEK
                        DPRT06( K,ISP ) = IDIUPS

                    ENDIF

                CASE( 7 )                   ! NOTE:  pol-specific was excluded
                    MPRT07( K ) = IMON
                    WPRT07( K ) = IWEK
                    DPRT07( K ) = IDIU

                CASE( 8 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT08( K,: ) = IMON
                        WPRT08( K,: ) = IWEK
                        DPRT08( K,: ) = IDIU

                    ELSE
                        MPRT08( K,ISP ) = IMON
                        WPRT08( K,ISP ) = IWEK
                        DPRT08( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 9 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT09( K,: ) = IMON
                        WPRT09( K,: ) = IWEK
                        DPRT09( K,: ) = IDIU

                    ELSE
                        MPRT09( K,ISP ) = IMON
                        WPRT09( K,ISP ) = IWEK
                        DPRT09( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 10 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT10( K,: ) = IMON
                        WPRT10( K,: ) = IWEK
                        DPRT10( K,: ) = IDIU

                    ELSE
                        MPRT10( K,ISP ) = IMON
                        WPRT10( K,ISP ) = IWEK
                        DPRT10( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 11 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT11( K,: ) = IMON
                        WPRT11( K,: ) = IWEK
                        DPRT11( K,: ) = IDIU

                    ELSE
                        MPRT11( K,ISP ) = IMON
                        WPRT11( K,ISP ) = IWEK
                        DPRT11( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 12 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT12( K,: ) = IMON
                        WPRT12( K,: ) = IWEK
                        DPRT12( K,: ) = IDIU

                    ELSE
                        MPRT12( K,ISP ) = IMON
                        WPRT12( K,ISP ) = IWEK
                        DPRT12( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 13 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT13( K,: ) = IMON
                        WPRT13( K,: ) = IWEK
                        DPRT13( K,: ) = IDIU

                    ELSE
                        MPRT13( K,ISP ) = IMON
                        WPRT13( K,ISP ) = IWEK
                        DPRT13( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 14 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT14( K,: ) = IMON
                        WPRT14( K,: ) = IWEK
                        DPRT14( K,: ) = IDIU

                    ELSE
                        MPRT14( K,ISP ) = IMON
                        WPRT14( K,ISP ) = IWEK
                        DPRT14( K,ISP ) = IDIUPS

                    ENDIF
                    
                CASE( 15 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT15( K,: ) = IMON
                        WPRT15( K,: ) = IWEK
                        DPRT15( K,: ) = IDIU

                    ELSE
                        MPRT15( K,ISP ) = IMON
                        WPRT15( K,ISP ) = IWEK
                        DPRT15( K,ISP ) = IDIUPS	

                    ENDIF
                                        
                CASE( 16 )
                    IF( ISP .EQ. 0 ) THEN
                        MPRT16( K,: ) = IMON
                        WPRT16( K,: ) = IWEK
                        DPRT16( K,: ) = IDIU

                    ELSE
                        MPRT16( K,ISP ) = IMON
                        WPRT16( K,ISP ) = IWEK
                        DPRT16( K,ISP ) = IDIUPS	

                    ENDIF
                                        
                CASE DEFAULT

                END SELECT

            ENDDO                            ! End Loop on sorted x-ref entries

        ELSEIF( SFLAG ) THEN

            DO I = 1, NXREF

                J      = INDXTA( I )
                ISP    = ISPTA ( J )
                SPCODE = CSPRNA( J )

                T      = XTYPE ( I )
                K      = XTCNT ( I )

                IF( ISP .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'INTERNAL ERROR: Pollutant code was 0 ' //
     &                     'for speciation cross-reference entry ' //
     &                     CRLF()// BLANK10// 'in program ' // PROGNAME
                    CALL M3MESG( MESG )
                ENDIF

C.................  Populate tables depending on type. Note that all profiles 
C                   are, by definition, pollutant-specific.
                SELECT CASE ( T )

                CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate

                CASE( 1 )
                    CSPT01( ISP ) = SPCODE  ! Arrays

                CASE( 2 )
                    CSPT02( K,ISP ) = SPCODE

                CASE( 3 )
                    CSPT03( K,ISP ) = SPCODE

                CASE( 4 )                   ! NOTE:  pol-specific was excluded
                    CSPT04( K,ISP ) = SPCODE

                CASE( 5 )
                    CSPT05( K,ISP ) = SPCODE

                CASE( 6 )
                    CSPT06( K,ISP ) = SPCODE

                CASE( 7 )                   ! NOTE:  pol-specific was excluded
                    CSPT07( K,ISP ) = SPCODE

                CASE( 8 )
                    CSPT08( K,ISP ) = SPCODE

                CASE( 9 )
                    CSPT09( K,ISP ) = SPCODE
                    
                CASE( 10 )
                    CSPT10( K,ISP ) = SPCODE

                CASE( 11 )
                    CSPT11( K,ISP ) = SPCODE

                CASE( 12 )
                    CSPT12( K,ISP ) = SPCODE

                CASE( 13 )
                    CSPT13( K,ISP ) = SPCODE

                CASE( 14 )
                    CSPT14( K,ISP ) = SPCODE

                CASE( 15 )
                    CSPT15( K,ISP ) = SPCODE

                CASE( 16 )
                    CSPT16( K,ISP ) = SPCODE

                CASE DEFAULT

                END SELECT

            ENDDO                            ! End Loop on sorted x-ref entries

        END IF   ! End of different storage for different operation types

        IF( EFLAG ) THEN
            MESG = 'Problem processing cross-reference records.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store count of records in each group in final variable
        DO I = 1, NXTYPES
            TXCNT( I ) = N( I )
        ENDDO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram writes a warning message for 
C               duplicate entries in the cross-reference file.
            SUBROUTINE REPORT_DUP_XREF

C.............  Local variables
            INTEGER       :: L1
            INTEGER       :: L2
            CHARACTER*300 :: BUFFER
            CHARACTER*300 :: MESG

C......................................................................

            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

            MESG = 'WARNING: Duplicate entry in ' // OPTYPE( 1:LOPT ) //
     &             ' x-ref file:' // CRLF() // BLANK10 //
     &             BUFFER( 1:L2 )

            IF( ISP .GT. 0 ) THEN
                L1 = LEN_TRIM( MESG )
                MESG = MESG( 1:L1 ) // ' POL:' // EINAM( ISP )
            END IF

            CALL M3MSG2( MESG )

            END SUBROUTINE REPORT_DUP_XREF

C----------------------------------------------------------------------

C.............  This internal subprogram writes a warning message for 
C               invalid cross-reference entries
            SUBROUTINE REPORT_INVALID_XREF( INMESG )

C.............  Subprogram arguments
            CHARACTER(*) INMESG       ! Input message to include on line 2

C.............  Local parameters
            CHARACTER*48, PARAMETER :: PART1 = 
     &         'WARNING: Skipping invalid cross-reference entry.'

C.............  Local variables
            INTEGER       :: L1
            INTEGER       :: L2
            CHARACTER*300 :: BUFFER
            CHARACTER*300 :: MESG

C......................................................................

            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L1 )

            L1 = LEN_TRIM( INMESG )
            L2 = LEN_TRIM( BUFFER )

            MESG = PART1  // 
     &             CRLF() // BLANK10 // INMESG( 1:L1 ) // ':' //
     &             CRLF() // BLANK10 // BUFFER( 1:L2 )
     &                
            IF( ISP .GT. 0 ) THEN
                L1 = LEN_TRIM( MESG )
                MESG = MESG( 1:L1 ) // ' POL:' // EINAM( ISP )
            END IF

            CALL M3MESG( MESG )

            END SUBROUTINE REPORT_INVALID_XREF

C----------------------------------------------------------------------

        END SUBROUTINE PXREFTBL
