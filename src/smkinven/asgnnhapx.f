
        SUBROUTINE ASGNNHAPX( NRAWBP )

C***********************************************************************
C  subroutine body starts at line 107
C
C  DESCRIPTION:
C      For each source, find the most specific match to determine
C      if the source should be excluded from the calculation of NONHAPVOC
C      pollutant in the inventory.  The algorithm allows for all four
C      levels of SCC matching.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/02 by M. Houyoux
C
C************************************************************************
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
C...........   This module contains the source ararys
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDC

        EXTERNAL        CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NRAWBP  ! no. raw records by pollutant

C.........  Other local variables
        INTEGER          I, J, LS, S    !  counters and indices

        INTEGER          IOS         ! i/o status
        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          F7, F8, F9, F10, F11        ! tmp find indices

        LOGICAL       :: EFLAG    = .FALSE.

        CHARACTER*256             MESG     ! message buffer
        CHARACTER(LEN=STALEN3) CSTA     ! tmp Country/state code
        CHARACTER(LEN=STSLEN3) CSTAS_A  ! tmp Country/state code // level 1 SCC
        CHARACTER(LEN=STSLEN3) CSTAS_B  ! tmp Country/state code // level 2 SCC
        CHARACTER(LEN=STSLEN3) CSTAS_C  ! tmp Country/state code // level 3 SCC
        CHARACTER(LEN=STSLEN3) CSTAS_D  ! tmp Country/state code // all SCC
        CHARACTER(LEN=SCCLEN3) TSCC_A   ! tmp level 1 of SCC
        CHARACTER(LEN=SCCLEN3) TSCC_B   ! tmp level 2 of SCC
        CHARACTER(LEN=SCCLEN3) TSCC_C   ! tmp level 3 of SCC
        CHARACTER(LEN=SCCLEN3) TSCC_D   ! tmp level 4 (all) of SCC
        CHARACTER(LEN=SRCLEN3) CSRC     ! tmp source chars string
        CHARACTER(LEN=SS0LEN3) CSSC0    ! tmp FIPS // Plant // SCC
        CHARACTER(LEN=SS1LEN3) CSSC1    ! tmp source chars through char1 // SCC
        CHARACTER(LEN=SS2LEN3) CSSC2    ! tmp source chars through char2 // SCC
        CHARACTER(LEN=SS3LEN3) CSSC3    ! tmp source chars through char3 // SCC
        CHARACTER(LEN=SS4LEN3) CSSC4    ! tmp source chars through char4 // SCC
        CHARACTER(LEN=SS5LEN3) CSSC5    ! tmp source chars through char5 // SCC
        CHARACTER(LEN=FIPLEN3) CFIP     ! tmp (character) FIPS code
        CHARACTER(LEN=FPLLEN3) CFIPPLT  ! tmp FIPS code // plant id
        CHARACTER(LEN=FPSLEN3) CFIPS_A  ! tmp FIPS code // level 1 of SCC
        CHARACTER(LEN=FPSLEN3) CFIPS_B  ! tmp FIPS code // level 2 of SCC
        CHARACTER(LEN=FPSLEN3) CFIPS_C  ! tmp FIPS code // level 3 of SCC
        CHARACTER(LEN=FPSLEN3) CFIPS_D  ! tmp FIPS code // level 4 (all) of SCC

        CHARACTER*16 :: PROGNAME = 'ASGNNHAPX' ! program name

C***********************************************************************
C   begin body of subroutine ASGNNHAPX

C.........  Allocate memory for source-based array from MODLISTS
        ALLOCATE( LNONHAP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LNONHAP', PROGNAME )
        LNONHAP = .TRUE.   ! array

C.........  Loop through the unsorted sources
        S  = 0
        DO I = 1, NRAWBP

C............  Set source number from previous iteration
            LS = S

C............  Get sorted position and source number for this source
            J = INDEXA( I )
            S = SRCIDA( I )

C............  If record has the same source as the previous iteration, 
C              skip it
            IF ( S .EQ. LS ) CYCLE

C.............  Create selection 
            SELECT CASE ( CATEGORY )

            CASE ( 'AREA', 'MOBILE' )
                CSRC    = CSOURCA( J )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                TSCC_D  = CSCCA( J )            ! Level 4 (all)
                TSCC_A  = TSCC_D( 1:SCCLEV1 )   ! Level 1
                TSCC_B  = TSCC_D( 1:SCCLEV2 )   ! Level 2
                TSCC_C  = TSCC_D( 1:SCCLEV3 )   ! Level 3

                CFIPS_A = CFIP // TSCC_A
                CFIPS_B = CFIP // TSCC_B
                CFIPS_C = CFIP // TSCC_C
                CFIPS_D = CFIP // TSCC_D

                CSTAS_A = CSTA // TSCC_A
                CSTAS_B = CSTA // TSCC_B
                CSTAS_C = CSTA // TSCC_C
                CSTAS_D = CSTA // TSCC_D

            END SELECT

C.............  Try for any FIPS code & SCC matches; then any Cy/st code & 
C               SCC matches
            F11= FINDC( CFIPS_D, TXCNT( 9 ), CHRT09 ) 
            F10= FINDC( CFIPS_C, TXCNT( 25), CHRT08C ) 
            F9 = FINDC( CFIPS_B, TXCNT( 24), CHRT08B ) 
            F8 = FINDC( CFIPS_A, TXCNT( 23), CHRT08A ) 
            F7 = FINDC( CSTAS_D, TXCNT( 6 ), CHRT06 ) 
            F6 = FINDC( CSTAS_C, TXCNT( 22), CHRT05C ) 
            F5 = FINDC( CSTAS_B, TXCNT( 21), CHRT05B ) 
            F4 = FINDC( CSTAS_A, TXCNT( 20), CHRT05A ) 

            IF( F11 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F10 .GT. 0 )THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F9 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F8 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F7 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F6 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F5 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F4 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Now try for county or state match prior to SCC-only match,
C               which is different from other SMOKE assignment priorities.
C               This is different because it is likely that an entire 
C               state or county must be excluded from the NONHAP calculation.

C.............  Try for any FIPS code match
            F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

            IF( F0 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Try for any country/state code match 
            F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

            IF( F0 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Try for SCC matches
            F3 = FINDC( TSCC_D , TXCNT( 3 ), CHRT03 ) 
            F2 = FINDC( TSCC_C , TXCNT( 19), CHRT02C ) 
            F1 = FINDC( TSCC_B , TXCNT( 18), CHRT02B ) 
            F0 = FINDC( TSCC_A , TXCNT( 17), CHRT02A ) 

            IF( F3 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F2 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F1 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            ELSE IF( F0 .GT. 0 ) THEN
                LNONHAP( S ) = .FALSE.
                CYCLE                       !  to end of sources-loop
            END IF

        END DO        !  end loop on source x pollutants

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ASGNNHAPX
