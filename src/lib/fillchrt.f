
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
        USE MODXREF

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
        INTEGER     , INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
        INTEGER     , INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

C...........   Local field position array
        INTEGER, ALLOCATABLE :: ENDLEN( : ) 

C...........   Other local variables
        INTEGER       I, J, K, T       ! counter and indices

        INTEGER       IOS              ! i/o status

        LOGICAL    :: EFLAG = .FALSE.  ! true: error was detected

        CHARACTER*300          MESG    ! message buffer

        CHARACTER(LEN=STALEN3) CSTA    ! temporary (character) state code
        CHARACTER(LEN=SCCLEN3) SCCL    ! left digits of TSCC
        CHARACTER(LEN=FIPLEN3) CFIP    ! temporary (character) FIPS code
        CHARACTER(LEN=SRCLEN3) CSRC    ! temporary source char string
        CHARACTER(LEN=SCCLEN3) TSCC    ! temporary SCC

        CHARACTER*16 :: PROGNAME = 'FILLCHRT' ! program name

C***********************************************************************
C   begin body of subroutine FILLCHRT

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

C.............  Set up partial strings for saving
            SCCL   = TSCC( 1:LSCCEND )
            CFIP   = CSRC( 1:FIPLEN3 ) 
            CSTA   = CSRC( 1:STALEN3 )

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
                CHRT12( K ) = CSRC( 1:ENDLEN( 3 ) ) // TSCC
            CASE( 13 )
                CHRT13( K ) = CSRC( 1:ENDLEN( 4 ) ) // TSCC
            CASE( 14 )
                CHRT14( K ) = CSRC( 1:ENDLEN( 5 ) ) // TSCC
            CASE( 15 )
                CHRT15( K ) = CSRC( 1:ENDLEN( 6 ) ) // TSCC
            CASE( 16 )
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
