
        SUBROUTINE FILLCHRT( NXREF, XTYPE, XTCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine populates the source-characteristics part of the grouped
C      cross-reference tables. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
C
C****************************************************************************/
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
        USE MODXREF, ONLY: CHRT02, CHRT03, CHRT04, CHRT05,
     &          CHRT06, CHRT07, CHRT08, CHRT09, CHRT10,
     &          CHRT11, CHRT12, CHRT13, CHRT14, CHRT15, CHRT16,
     &          CHRT02A, CHRT02B, CHRT02C,
     &          CHRT05A, CHRT05B, CHRT05C,
     &          CHRT08A, CHRT08B, CHRT08C,
     &          CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31,
     &          CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37,
     &          INDXTA, CSRCTA, CSCCTA, CMACTA 

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, LSCCEND, SCCLEV1, SCCLEV2, SCCLEV3

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         SETSCCTYPE
        
        EXTERNAL        SETSCCTYPE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
        INTEGER     , INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
        INTEGER     , INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

C...........   Local field position array
        INTEGER, ALLOCATABLE :: ENDLEN( : ) 

C...........   Other local variables
        INTEGER       I, J, K, L, T    ! counter and indices

        INTEGER       IOS              ! i/o status

        LOGICAL    :: EFLAG = .FALSE.  ! true: error was detected
        LOGICAL       SCCFLAG          ! true: SCC type is different from previous

        CHARACTER*300          MESG    ! message buffer

        CHARACTER(LEN=STALEN3) CSTA    ! temporary (character) state code
        CHARACTER(LEN=SCCLEN3) SCCL    ! left digits of TSCC
        CHARACTER(LEN=FIPLEN3) CFIP    ! temporary (character) FIPS code
        CHARACTER(LEN=SRCLEN3) CSRC    ! temporary source char string
        CHARACTER(LEN=SCCLEN3) TSCC    ! temporary SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO ! buffer for zero SCC
        CHARACTER(LEN=SICLEN3) CSIC    ! buffer for SIC
        CHARACTER(LEN=SICLEN3) CSICL   ! buffer for left 2-digit SIC
        CHARACTER(LEN=MACLEN3) CMCT    ! buffer for MACT code

        CHARACTER*16 :: PROGNAME = 'FILLCHRT' ! program name

C***********************************************************************
C   begin body of subroutine FILLCHRT

C.........  Set up zero string for SCC code of zero
        SCCZERO = REPEAT( '0', SCCLEN3 )

C.........  Set the local field position array based on the source category
        SELECT CASE ( CATEGORY )
        CASE( 'AREA' ) 
            ALLOCATE( ENDLEN( MXARCHR3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ENDLEN', PROGNAME )
            ENDLEN = 1  ! array
            ENDLEN( 1:MXARCHR3 ) = ARENDL3( 1:MXARCHR3 )

        CASE( 'MOBILE' )
            ALLOCATE( ENDLEN( MXMBCHR3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ENDLEN', PROGNAME )
            ENDLEN = 1  ! array
            ENDLEN( 1:MXMBCHR3 ) = MBENDL3( 1:MXMBCHR3 )

        CASE( 'POINT' )
            ALLOCATE( ENDLEN( MXPTCHR3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ENDLEN', PROGNAME )
            ENDLEN = 1  ! array
            ENDLEN( 1:MXPTCHR3 ) = PTENDL3( 1:MXPTCHR3 )

        END SELECT

C.........  Store the source characteristics for each x-ref entry, depending
C           on the group (XTYPE) and the position in that group (XTCNT)

        DO I = 1, NXREF
            J      = INDXTA( I )
            CSRC   = CSRCTA( J )
            TSCC   = CSCCTA( J )
            
            IF( ALLOCATED( CMACTA ) ) THEN
                CMCT = CMACTA( J )
            ELSE
                CMCT = ' '
            END IF

C.............  Set up partial strings for country/state/county
            CFIP   = CSRC( 1:FIPLEN3 ) 
            CSTA   = CSRC( 1:STALEN3 )

C.............  Determine whether SIC is imbedded in SCC field
            K = INDEX( TSCC, SICNOTE )

C.............  If SIC imbedded, setup SIC fields
            CSIC = ' '
            IF( K .GT. 0 ) THEN
                L = K + LEN( SICNOTE )
                CSIC  = TSCC( L : L + SICLEN3 - 1 )
                TSCC  = SCCZERO
                CSICL = CSIC( 1:2 )
            END IF

C.............  Set up partial SCC strings for saving

C.............  Set type of SCC                
            SCCFLAG = SETSCCTYPE( TSCC )
            SCCL   = TSCC( 1:LSCCEND )

            T      = XTYPE ( I )  ! extract what group this entry is in
            K      = XTCNT ( I )  ! extract position in that group

            SELECT CASE ( T )

            CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate
            CASE( 1 )  ! Skip this because no source-characteristics
            CASE( 2 )
                CHRT02( K ) = SCCL
            CASE( 3 )
                CHRT03( K ) = TSCC
            CASE( 4 )
                CHRT04( K ) = CSTA
            CASE( 5 )
                CHRT05( K ) = CSTA // SCCL
            CASE( 6 )
                CHRT06( K ) = CSTA // TSCC
            CASE( 7 )
                CHRT07( K ) = CFIP
            CASE( 8 )
                CHRT08( K ) = CFIP // SCCL
            CASE( 9 )
                CHRT09( K ) = CFIP // TSCC
            CASE( 10 )
                CHRT10( K ) = CSRC( 1:ENDLEN( 2 ) )
            CASE( 11 )
                CHRT11( K ) = CSRC( 1:ENDLEN( 2 ) ) // TSCC
            CASE( 12 )
                IF( TSCC .EQ. SCCZERO ) TSCC = ' '
                CHRT12( K ) = CSRC( 1:ENDLEN( 3 ) ) // TSCC
            CASE( 13 )
                IF( TSCC .EQ. SCCZERO ) TSCC = ' '
                CHRT13( K ) = CSRC( 1:ENDLEN( 4 ) ) // TSCC
            CASE( 14 )
                IF( TSCC .EQ. SCCZERO ) TSCC = ' '
                CHRT14( K ) = CSRC( 1:ENDLEN( 5 ) ) // TSCC
            CASE( 15 )
                IF( TSCC .EQ. SCCZERO ) TSCC = ' '
                CHRT15( K ) = CSRC( 1:ENDLEN( 6 ) ) // TSCC
            CASE( 16 )
                IF( TSCC .EQ. SCCZERO ) TSCC = ' '
                CHRT16( K ) = CSRC( 1:ENDLEN( 7 ) ) // TSCC

C.............  NOTE- Cases added in version 1.4 (initially) for cntl/proj only
            CASE( 17 )
                CHRT02A( K ) = TSCC( 1:SCCLEV1 )
            CASE( 18 )
                CHRT02B( K ) = TSCC( 1:SCCLEV2 )
            CASE( 19 )
                CHRT02C( K ) = TSCC( 1:SCCLEV3 )
            CASE( 20 )
                CHRT05A( K ) = CSTA // TSCC( 1:SCCLEV1 )
            CASE( 21 )
                CHRT05B( K ) = CSTA // TSCC( 1:SCCLEV2 )
            CASE( 22 )
                CHRT05C( K ) = CSTA // TSCC( 1:SCCLEV3 )
            CASE( 23 )
                CHRT08A( K ) = CFIP // TSCC( 1:SCCLEV1 )
            CASE( 24 )
                CHRT08B( K ) = CFIP // TSCC( 1:SCCLEV2 )
            CASE( 25 )
                CHRT08C( K ) = CFIP // TSCC( 1:SCCLEV3 )
            CASE( 26 )
                CHRT26( K ) = CSIC( 1:2 ) 
            CASE( 27 )
                CHRT27( K ) = CSIC 
            CASE( 28 )
                CHRT28( K ) = CSTA // CSIC( 1:2 ) 
            CASE( 29 )
                CHRT29( K ) = CSTA // CSIC 
            CASE( 30 )
                CHRT30( K ) = CFIP // CSIC( 1:2 ) 
            CASE( 31 )
                CHRT31( K ) = CFIP // CSIC 
                
C.............  MACT based cases
            CASE( 32 )
                CHRT32( K ) = CMCT
            CASE( 33 )
                CHRT33( K ) = TSCC // CMCT
            CASE( 34 )
                CHRT34( K ) = CSTA // CMCT
            CASE( 35 )
                CHRT35( K ) = CSTA // TSCC // CMCT
            CASE( 36 )
                CHRT36( K ) = CFIP // CMCT
            CASE( 37 )
                CHRT37( K ) = CFIP // TSCC // CMCT
                
            CASE DEFAULT

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'INTERNAL ERROR: Cross-' // 
     &                 'reference category', T, 
     &                 'not known in subroutine ' // PROGNAME
                CALL M3MESG( MESG )

            END SELECT

        ENDDO                            ! End Loop on sorted x-ref entries

C.........  Deallocate local memory
        DEALLOCATE( ENDLEN )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE FILLCHRT
