
        SUBROUTINE XREFTBL( OPTYPE, NXREF )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine will do the logic for grouping the cross references for
C      any source category.  In the logic phase, it will merely store the group 
C      number and the count in the group for each x-ref entry. Then, depending on
C      the type of operation (e.g., temporal, speciation, reactivity packet), 
C      it will allocate memory for an populate the appropriate grouped 
C      tables.  The unsorted x-ref data and grouped tables will be passed by a
C      module for cross-reference data.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
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

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         STR2INT

        EXTERNAL   CRLF, STR2INT 

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: OPTYPE ! operation type (tmprl,spec,ctg...)
        INTEGER     , INTENT (IN) :: NXREF  ! no. ungrouped x-ref entries

C...........   Local parameters
        INTEGER, PARAMETER :: SNFLEN3 = SRCLEN3 - FIPLEN3

C...........   Arrays for intial pass through x-ref to determine degree of 
C              each record and store the type and count of that type
        INTEGER         XTYPE( NXREF )  ! group number of x-ref entry
        INTEGER         XTCNT( NXREF )  ! position in group of x-ref entry

C...........  Arrays for counting number of x-ref records in each matching
C             degree and comparing current record with previous record as part
C             of the counting. Note that although all arrays below are 
C             available for each degree of matching, only the array elements 
C             that are appropriate for a given degree are actually populated. 
        INTEGER                N     ( 0:NXTYPES ) ! cnt for degree of matching
        INTEGER                PIFIP (   NXTYPES ) ! previous co/st/cy code

        CHARACTER(LEN=SCCLEN3) PTSCC ( NXTYPES )  ! previous SCC
        CHARACTER(LEN=SRCLEN3) PCSRC ( NXTYPES )  ! previous CSRC
        CHARACTER(LEN=SS5LEN3) PCSSC ( NXTYPES )  ! previous CSRC(part) // SCC

C...........   Array of source characeristics
        CHARACTER*300           CHARS( MXCHRS )

C...........   Other local variables
        INTEGER       I, J, J1, J2, L  ! counter and indices

        INTEGER       ICYID            ! temporary county code
        INTEGER       IDUM             ! dummy integer
        INTEGER       IFIP             ! temporary FIPS code
        INTEGER       IOS              ! i/o status
        INTEGER       ISP              ! temporary pollutant position in EANAM
        INTEGER       LOPT             ! length of OPTYPE
        INTEGER       L1R, L2R, L3R    ! SCC level 1, 2, and 3 right positions
        INTEGER       NT               ! code for specificity of x-ref entry
        INTEGER    :: PISP = IMISS3    ! previous iteration ISP

        LOGICAL    :: CFLAG = .FALSE.  ! true: operation type is control cntls
        LOGICAL    :: DEFAULT( NIPPA ) ! true: if default entry in x-ref
        LOGICAL    :: DFLAG = .FALSE.  ! true: operation type is additive cntls
        LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred
        LOGICAL    :: FFLAG = .FALSE.  ! true: operation type is emis. factors
        LOGICAL    :: GFLAG = .FALSE.  ! true: operation type is ctg cntls
        LOGICAL    :: IFLAG = .FALSE.  ! true: operation type is gridding
        LOGICAL    :: JFLAG = .FALSE.  ! true: operation type is projection
        LOGICAL    :: LFLAG = .FALSE.  ! true: operation type is allowable cntls
        LOGICAL    :: MFLAG = .FALSE.  ! true: operation type is VMT mix
        LOGICAL    :: NFLAG = .TRUE.   ! true: has pol or activity-specific
        LOGICAL    :: OFLAG = .FALSE.  ! true: any control packet
        LOGICAL    :: PFLAG = .FALSE.  ! true: operation type is speeds
        LOGICAL    :: POADFLT          ! true: okay to have pol/act-spec dfaults
        LOGICAL    :: RFLAG = .FALSE.  ! true: operation type is reactivty cntls
        LOGICAL    :: SFLAG = .FALSE.  ! true: operation type is speciation
        LOGICAL    :: TFLAG = .FALSE.  ! true: operation type is temporal

        CHARACTER*1            CDUM          ! dummy character string
        CHARACTER*300          BUFFER        ! source definition buffer
        CHARACTER*300          MESG          ! message buffer

        CHARACTER(LEN=STALEN3) CSTA          ! temporary (character) state code
        CHARACTER(LEN=SCCLEN3) SCCL          ! left digits of TSCC
        CHARACTER(LEN=SCCLEN3) SCCR          ! RHS digits of TSCC (old method)
        CHARACTER(LEN=SCCLEN3) SCCR_A        ! RHS for level 1 of SCC
        CHARACTER(LEN=SCCLEN3) SCCR_B        ! RHS for level 2 of SCC
        CHARACTER(LEN=SCCLEN3) SCCR_C        ! RHS for level 3 of SCC
        CHARACTER(LEN=SCCLEN3) SCRZERO       ! buf for 0 right 5 digits of TSCC
        CHARACTER(LEN=SNFLEN3) CNFIP         ! characterstics without FIPS code
        CHARACTER(LEN=SRCLEN3) CSRC          ! temporary source char string
        CHARACTER(LEN=FIPLEN3) CFIP          ! temporary (character) FIPS code
        CHARACTER(LEN=FIPLEN3) FIPZERO       ! buffer for zero FIPS code
        CHARACTER(LEN=SCCLEN3) PSCC          ! previous SCC
        CHARACTER(LEN=SCCLEN3) TSCC          ! temporary SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO       ! buffer for zero SCC
        CHARACTER(LEN=SCCLEN3) SCCZ_A        ! buffer for level 1 zero string
        CHARACTER(LEN=SCCLEN3) SCCZ_B        ! buffer for level 2 zero string
        CHARACTER(LEN=SCCLEN3) SCCZ_C        ! buffer for level 3 zero string
        CHARACTER(LEN=SS5LEN3) CSRCSCC       ! buffer for source // SCC

        CHARACTER*16 :: PROGNAME = 'XREFTBL' ! program name

C***********************************************************************
C   begin body of subroutine XREFTBL

        LOPT = LEN_TRIM( OPTYPE )

C.........  Check for valid operation type
        SELECT CASE( OPTYPE )

        CASE( 'ADD' )
            POADFLT = .TRUE.
            DFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'ALLOWABLE' )
            POADFLT = .TRUE.
            LFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'CONTROL' )
            POADFLT = .TRUE.
            CFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'EMS_CONTROL' )
            POADFLT = .TRUE.
            CFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'CTG' )
            POADFLT = .TRUE.
            GFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'EMISFACS' )
            POADFLT = .TRUE.
            FFLAG   = .TRUE.
        CASE( 'GRIDDING' )
            POADFLT = .FALSE.
            IFLAG   = .TRUE.
            NFLAG   = .FALSE.
        CASE( 'PROJECTION' )
            POADFLT = .FALSE.
            JFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'REACTIVITY' )
            POADFLT = .FALSE.
            RFLAG   = .TRUE.
            OFLAG   = .TRUE.
        CASE( 'SPEED' ) 
            POADFLT = .FALSE.
            PFLAG   = .TRUE.
            NFLAG   = .FALSE.
        CASE( 'SPECIATION' )
            POADFLT = .TRUE.
            SFLAG   = .TRUE.
        CASE( 'TEMPORAL' ) 
            POADFLT = .TRUE.
            TFLAG   = .TRUE.
        CASE( 'VMTMIX' ) 
            POADFLT = .FALSE.
            MFLAG   = .TRUE.
            NFLAG   = .FALSE.

        CASE DEFAULT

            MESG = 'INTERNAL ERROR: Operation type "' // 
     &             OPTYPE( 1:LOPT ) //
     &             '" not known in subroutine ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT
        
C.........  Set up zero strings for FIPS code of zero and SCC code of zero
        L1R = SCCLEV1 + 1
        L2R = SCCLEV2 + 1
        L3R = SCCLEV3 + 1
        FIPZERO = REPEAT( '0', FIPLEN3 )
        SCCZERO = REPEAT( '0', SCCLEN3 )
        SCRZERO = REPEAT( '0', SCCLEN3 - LSCCEND )
        SCCZ_A = REPEAT( '0', SCCLEN3 - SCCLEV1 )
        SCCZ_B = REPEAT( '0', SCCLEN3 - SCCLEV2 )
        SCCZ_C = REPEAT( '0', SCCLEN3 - SCCLEV3 )

C.........  Initialize default array
        DEFAULT = .FALSE.   ! array

C.........  Initialize arrays for counting number of x-ref records in each
C           degree of matching
        N      = 0    ! arrays
        PIFIP  = 0
        PCSRC  = ' '
        PCSSC  = ' '
        PTSCC  = ' '

C.........  Initialize source characteristics
        CHARS = ' '  ! array

C.........  Loop through and count entries of each type. Store type.
C.........  For CSRC, don't include pollutant for grouping.
        ISP    = 0
        PSCC   = ' '
        DO I = 1, NXREF

            J = INDXTA( I )

            CSRC    = CSRCTA( J )( 1:SC_ENDP( NCHARS ) )
            TSCC    = CSCCTA( J )
            IF( NFLAG ) ISP = ISPTA ( J )  ! no pollutants for gridding

            DO J = 1, NCHARS
                CHARS( J ) = CSRC( SC_BEGP( J ):SC_ENDP( J ) )
            END DO

C.............  Rearrange CHARS if SCC is a part of the source definition
C               because we have now stored SCC separately.  It will be
C               much easier to group the source characteristics this way.
C               Still go back to using original definition when storing
C               in tables (i.e., use CSRC)
            IF( JSCC .GT. 0 ) THEN
                DO J = JSCC, NCHARS - 1
                    CHARS( J ) = CHARS( J + 1 )
                END DO
                CHARS( NCHARS ) = ' '
            END IF

C.............  Set up partial strings for checking
            CFIP    = CHARS( 1 )
            IFIP    = STR2INT( CFIP )         ! For checking previous
            CSTA    = CFIP( 1:STALEN3 )
            ICYID   = IFIP - STR2INT( CSTA ) * 1000
            SCCL    = TSCC(       1:LSCCEND )
            SCCR    = TSCC( RSCCBEG:SCCLEN3 )
            CSRCSCC = CSRC // TSCC
            CNFIP   = CSRC( SC_BEGP( PLTIDX ):SC_ENDP( NCHARS ) )

C.............  More partial strings for special SCC levels matching
            SCCR_A = TSCC( L1R:SCCLEN3 )
            SCCR_B = TSCC( L2R:SCCLEN3 )
            SCCR_C = TSCC( L3R:SCCLEN3 )

C.............  Select cases
C.............  Note that since these are sorted in order of increasing FIPS
C               code, SCC, pollutant index, etc., that the entries with zero for
C               these characteristics will appear earlier in the sorted list
            IF( IFIP .EQ. 0 ) THEN                       ! FIPS code is default

                IF( TSCC .EQ. SCCZERO ) THEN                   ! SCC is default

                    IF( POADFLT ) THEN           ! Pollutant-specific permitted

                        IF( ISP .EQ. 0 ) THEN
                            NT = 1
                            N( NT ) = N( NT ) + 1

                        ELSEIF( ISP .NE. 0 .AND. .NOT.    ! Pollutant specified
     &                          DEFAULT( ISP )         ) THEN

                            DEFAULT( ISP ) = .TRUE.
                            NT = 1
                            N( NT ) = N( NT ) + 1

                        ELSE                                  ! Report and skip
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE                     ! Pollutant-specific not permitted

                        IF( ISP .EQ. 0 .AND. .NOT.    ! Pollutant not specified
     &                      DEFAULT( 1 )           ) THEN 

                            DEFAULT( 1 ) = .TRUE.
                            NT = 1
                            N( NT ) = N( NT ) + 1

                        ELSE IF( ISP .NE. 0 ) THEN            ! Report and skip
                            MESG = 'Cannot use pollutant-specific ' //
     &                             'ultimate-default'
                            CALL REPORT_INVALID_XREF( MESG )
                            NT = 0

                        ELSE IF( DEFAULT( 1 ) ) THEN          ! Report and skip
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    END IF

                ELSE IF( SCCR .EQ. SCRZERO ) THEN        ! Left SCC

                    NT = 2
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE                                         ! Complete SCC

                    NT = 3
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

C.................  Section for special SCC levels for controls. This is not an
C                   efficient way to implement this, but it's needed so long
C                   as the old Right-left method is still needed.
                IF ( OFLAG .AND. NT .NE. 1 ) THEN

                    IF( SCCR_A .EQ. SCCZ_A ) THEN         !  Level-1 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 17
                        IF( TSCC .NE. PTSCC( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PTSCC( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE IF( SCCR_B .EQ. SCCZ_B ) THEN    !  Level-2 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 18
                        IF( TSCC .NE. PTSCC( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PTSCC( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE IF( SCCR_C .EQ. SCCZ_C ) THEN    !  Level-3 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 19
                        IF( TSCC .NE. PTSCC( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PTSCC( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    END IF 

                END IF 

            ELSEIF( ICYID .EQ. 0 ) THEN            ! County code is default

                IF( TSCC .EQ. SCCZERO ) THEN         ! SCC code is default

                    NT = 4
                    IF( IFIP .NE. PIFIP( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP( NT ) = IFIP

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

c                    ELSE                                  ! Report and skip
c                        MESG = 'Cannot use pollutant-specific ' //
c     &                         'Country/State-default'
c                        CALL REPORT_INVALID_XREF( MESG )
c                        NT = 0
c                    END IF

                ELSEIF( SCCR .EQ. SCRZERO ) THEN         ! left SCC

                    NT = 5
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC
                       
                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE                                     ! Complete SCC

                    NT = 6
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

C.................  Section for special SCC levels for controls. This is not an
C                   efficient way to implement this, but it's needed so long
C                   as the old Right-left method is still needed.
                IF ( OFLAG .AND. NT .NE. 4 ) THEN

                    IF( SCCR_A .EQ. SCCZ_A ) THEN         !  State/Level-1 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 20
                        IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                      TSCC .NE. PTSCC( NT )      ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP ( NT ) = IFIP
                            PTSCC ( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE IF( SCCR_B .EQ. SCCZ_B ) THEN    !  State/Level-2 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 21
                        IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                      TSCC .NE. PTSCC( NT )      ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP ( NT ) = IFIP
                            PTSCC ( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE IF( SCCR_C .EQ. SCCZ_C ) THEN    !  State/Level-3 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 22
                        IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                      TSCC .NE. PTSCC( NT )      ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP ( NT ) = IFIP
                            PTSCC ( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    END IF 

                END IF 

            ELSEIF( CNFIP .EQ. ' ' .OR.
     &              CNFIP .EQ. TSCC     ) THEN  ! Country/St/Co code is complete

                IF( TSCC .EQ. SCCZERO ) THEN            ! SCC code is default

                    NT = 7
                    IF( IFIP .NE. PIFIP( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP( NT ) = IFIP

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

c                    ELSE                                  ! Report and skip
c                        MESG = 'Cannot use pollutant-specific ' //
c     &                         'Country/State/County-default'
c                        CALL REPORT_INVALID_XREF( MESG )
c                        NT = 0
c                    END IF

                ELSEIF( SCCR .EQ. SCRZERO ) THEN        ! Left SCC

                    NT = 8
                    IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                  TSCC .NE. PTSCC( NT )      ) THEN
                        N( NT ) = N( NT ) + 1
                        PIFIP ( NT ) = IFIP
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

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
                    END IF

                END IF                                              ! End SCC 

C.................  Section for special SCC levels for controls. This is not an
C                   efficient way to implement this, but it's needed so long
C                   as the old Right-left method is still needed.
                IF ( OFLAG .AND. NT .NE. 7 ) THEN

                    IF( SCCR_A .EQ. SCCZ_A ) THEN         !  FIPS/Level-1 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 23
                        IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                      TSCC .NE. PTSCC( NT )      ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP ( NT ) = IFIP
                            PTSCC ( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE IF( SCCR_B .EQ. SCCZ_B ) THEN    !  FIPS/Level-2 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 24
                        IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                      TSCC .NE. PTSCC( NT )      ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP ( NT ) = IFIP
                            PTSCC ( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    ELSE IF( SCCR_C .EQ. SCCZ_C ) THEN    !  FIPS/Level-3 SCC

                        N( NT ) = N( NT ) - 1
                        NT = 25
                        IF( IFIP .NE. PIFIP( NT ) .OR. 
     &                      TSCC .NE. PTSCC( NT )      ) THEN
                            N( NT ) = N( NT ) + 1
                            PIFIP ( NT ) = IFIP
                            PTSCC ( NT ) = TSCC

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF

                    END IF 

                END IF 

            ELSE                                        ! Plant is specified

                IF( TSCC .EQ. SCCZERO ) THEN           ! SCC code is default

C.....................  Loop through plant-specific characteristics. Only the
C                       plant is permitted to not have an SCC not specified
                    NT = 9 + MXCHRS - 1
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
                           END IF

                        ELSEIF( NT         .GT. 10  .AND. 
     &                          CHARS( J ) .NE. ' '       ) THEN

                            MESG = 'SCC is not specified'
                            CALL REPORT_INVALID_XREF( MESG )
                            NT = 0
                            EXIT                      ! End loop with NT

                        ELSEIF( NT .EQ. 10 ) THEN
                            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )

                            MESG = 'INTERNAL ERROR: Check XREFTBL ' //
     &                             'for processing record: ' //
     &                             CRLF() // BLANK10 // BUFFER( 1:L ) //
     &                             ' POA:' // EANAM( ISP )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                        END IF

                        NT = NT - 1

                    ENDDO           ! End loop on plant characteristics

                ELSEIF( SCCR .EQ. SCRZERO ) THEN         ! Left SCC

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )

                    MESG = 'Partial SCC "' // TSCC // '" is given ' //
     &                     'instead of full SCC'
                    CALL REPORT_INVALID_XREF( MESG )
                    NT = 0

                ELSE                                          ! Complete SCC
C.....................  Loop through plant-specific characteristics,
C                       and store the most specific entries first.
C.....................  Process NT 16 through 11
                    NT = 9 + MXCHRS
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
                            END IF

                        END IF

                        NT = NT - 1

                    ENDDO           ! End loop on plant characteristics

                END IF ! End SCC

            END IF ! End degree of Country/State/County code, or plant specified

            PISP = ISP

            XTYPE( I ) = NT
            XTCNT( I ) = N( NT )         ! Dimensioned from 0 so NT can = 0

        ENDDO                            ! End Loop on sorted x-ref entries

C.........  Allocate the memory for the source-characteristics portion of the
C           grouped cross-reference tables
        CALL ALOCCHRT( N( 1 ) )

C.........  Populate the grouped tables of cross-reference source-characteristics 
        CALL FILLCHRT( NXREF, XTYPE, XTCNT( 1 ) )

C.........  Depending on the operation type, first allocate memory for the tables
C           and initialize profile codes where needed.
C.........  Then, populate the tables from the sorted and post-processed 
C           cross-reference tables
 
C.........  Temporal x-ref tables
        IF( TFLAG ) THEN

            CALL ALOCTTBL( NIPPA, N( 1 ) )
            CALL FILLTTBL( NIPPA, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) ) 

C.........  Speciation x-ref tables
        ELSE IF( SFLAG ) THEN 

            CALL ALOCSTBL( NIPPA, N( 1 ) )
            CALL FILLSTBL( NIPPA, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) ) 

C.........  Gridding x-ref tables
        ELSE IF( IFLAG ) THEN
            CALL ALOCGTBL( N( 1 ) )
            CALL FILLGTBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) ) 

C.........  Emission factor x-ref tables
        ELSE IF( FFLAG ) THEN
            CALL ALOCETBL( NIACT, N( 1 ) )
            CALL FILLETBL( NIACT, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) ) 

C.........  Vehicle mix
        ELSE IF( MFLAG ) THEN
             CALL ALOCMTBL( N( 1 ) )
             CALL FILLMTBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

C.........  Speeds
        ELSE IF( PFLAG ) THEN
             CALL ALOCPTBL( N( 1 ) )
             CALL FILLPTBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

C.........  All control x-ref tables
        ELSE

            CALL ALOCCTBL( NIPPA, N( 1 ) )
            CALL FILLCTBL( NIPPA, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

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

            L1 = LEN_TRIM( MESG )
            IF( TSCC .NE. SCCZERO ) MESG = MESG( 1:L1 ) // 
     &                                     ' TSCC: ' // TSCC
            IF( ISP .GT. 0 ) THEN
                L1 = LEN_TRIM( MESG )
                MESG = MESG( 1:L1 ) // ' POL:' // EANAM( ISP )
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
                MESG = MESG( 1:L1 ) // ' POA:' // EANAM( ISP )
            END IF

            CALL M3MESG( MESG )

            END SUBROUTINE REPORT_INVALID_XREF

C----------------------------------------------------------------------

        END SUBROUTINE XREFTBL
