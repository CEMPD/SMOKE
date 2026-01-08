
        SUBROUTINE ASGNNHAPX

C***********************************************************************
C  subroutine body starts at line 112
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
C     09/2025 by HT UNC-IE:  Use M3UTILIO
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************
        USE M3UTILIO

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE MODSOURC, ONLY: CSOURC, CSCC

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: NHAP_EXCL, LNONHAP, TXCNT,
     &                     CHRT02A, CHRT02B, CHRT02C,
     &                     CHRT03, CHRT04,
     &                     CHRT05A, CHRT05B, CHRT05C,
     &                     CHRT06, CHRT07,
     &                     CHRT08A, CHRT08B, CHRT08C,
     &                     CHRT09, CHRT10, CHRT11, CHRT12,
     &                     CHRT13, CHRT14, CHRT15, CHRT16

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NSRC, SCCLEV1, SCCLEV2, SCCLEV3

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
c       CHARACTER(2)    CRLF
c       INTEGER         FINDC
c       LOGICAL         SETSCCTYPE

c       EXTERNAL        CRLF, FINDC, SETSCCTYPE
        LOGICAL, EXTERNAL :: SETSCCTYPE

C...........   SUBROUTINE ARGUMENTS

C.........  Other local variables
        INTEGER          I, J, LS, S    !  counters and indices

        INTEGER          IOS         ! i/o status
        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          F7, F8, F9, F10, F11        ! tmp find indices

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL          SCCFLAG           ! true: SCC type is different from previous

        CHARACTER(256)     MESG     ! message buffer
        CHARACTER(STALEN3) CSTA     ! tmp Country/state code
        CHARACTER(STSLEN3) CSTAS_A  ! tmp Country/state code // level 1 SCC
        CHARACTER(STSLEN3) CSTAS_B  ! tmp Country/state code // level 2 SCC
        CHARACTER(STSLEN3) CSTAS_C  ! tmp Country/state code // level 3 SCC
        CHARACTER(STSLEN3) CSTAS_D  ! tmp Country/state code // all SCC
        CHARACTER(SCCLEN3) TSCC_A   ! tmp level 1 of SCC
        CHARACTER(SCCLEN3) TSCC_B   ! tmp level 2 of SCC
        CHARACTER(SCCLEN3) TSCC_C   ! tmp level 3 of SCC
        CHARACTER(SCCLEN3) TSCC_D   ! tmp level 4 (all) of SCC
        CHARACTER(SRCLEN3) CSRC     ! tmp source chars string
        CHARACTER(SS0LEN3) CSSC0    ! tmp FIPS // Plant // SCC
        CHARACTER(SS1LEN3) CSSC1    ! tmp source chars through char1 // SCC
        CHARACTER(SS2LEN3) CSSC2    ! tmp source chars through char2 // SCC
        CHARACTER(SS3LEN3) CSSC3    ! tmp source chars through char3 // SCC
        CHARACTER(SS4LEN3) CSSC4    ! tmp source chars through char4 // SCC
        CHARACTER(SS5LEN3) CSSC5    ! tmp source chars through char5 // SCC
        CHARACTER(FIPLEN3) CFIP     ! tmp (character) FIPS code
        CHARACTER(FPLLEN3) CFIPPLT  ! tmp FIPS code // plant id
        CHARACTER(FPSLEN3) CFIPS_A  ! tmp FIPS code // level 1 of SCC
        CHARACTER(FPSLEN3) CFIPS_B  ! tmp FIPS code // level 2 of SCC
        CHARACTER(FPSLEN3) CFIPS_C  ! tmp FIPS code // level 3 of SCC
        CHARACTER(FPSLEN3) CFIPS_D  ! tmp FIPS code // level 4 (all) of SCC
        CHARACTER(SS5LEN3):: CSRC5=' '! tmp source chars through char5
        CHARACTER(SS4LEN3):: CSRC4=' '! tmp source chars through char4
        CHARACTER(SS3LEN3):: CSRC3=' '! tmp source chars through char3
        CHARACTER(SS2LEN3):: CSRC2=' '! tmp source chars through char2
        CHARACTER(SS1LEN3):: CSRC1=' '! tmp source chars through char1
        CHARACTER(SS5LEN3):: CHK16=' '! tmp source chars through char5// SCC
        CHARACTER(SS4LEN3):: CHK15=' '! tmp source chars through char4// SCC
        CHARACTER(SS3LEN3):: CHK14=' '! tmp source chars through char3// SCC
        CHARACTER(SS2LEN3):: CHK13=' '! tmp source chars through char2// SCC
        CHARACTER(SS1LEN3):: CHK12=' '! tmp source chars through char1// SCC
        CHARACTER(SS0LEN3):: CHK11=' '! tmp FIPS // Plant // SCC
        CHARACTER(FPLLEN3):: CHK10=' '! tmp FIPS code // plant id

        CHARACTER(16) :: PROGNAME = 'ASGNNHAPX' ! program name

C***********************************************************************
C   begin body of subroutine ASGNNHAPX

C.........  Allocate memory for source-based array from MODLISTS
        ALLOCATE( LNONHAP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LNONHAP', PROGNAME )
        LNONHAP = .TRUE.   ! array

C.........  Loop through the sorted sources
        DO S = 1, NSRC

C.............  Create selection 
            CSRC    = CSOURC( S )
            CFIP    = CSRC( 1:FIPLEN3 )
            CSTA    = CFIP( 1:STALEN3 )
                
C.............  Set type of SCC                
            SCCFLAG = SETSCCTYPE ( CSCC( S ) )
            TSCC_D  = CSCC( S )             ! Level 4 (all)
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

            IF ( CATEGORY == 'POINT' ) THEN

                CHK16   = CSRC( 1:PTENDL3( 7 ) ) // TSCC_D
                CHK15   = CSRC( 1:PTENDL3( 6 ) ) // TSCC_D
                CHK14   = CSRC( 1:PTENDL3( 5 ) ) // TSCC_D
                CHK13   = CSRC( 1:PTENDL3( 4 ) ) // TSCC_D
                CHK12   = CSRC( 1:PTENDL3( 3 ) ) // TSCC_D
                CHK11   = CSRC( 1:PTENDL3( 2 ) ) // TSCC_D
                CHK10   = CSRC( 1:PTENDL3( 2 ) )

                CSRC5   = CSRC( 1:PTENDL3( 7 ) ) 
                CSRC4   = CSRC( 1:PTENDL3( 6 ) ) 
                CSRC3   = CSRC( 1:PTENDL3( 5 ) ) 
                CSRC2   = CSRC( 1:PTENDL3( 4 ) ) 
                CSRC1   = CSRC( 1:PTENDL3( 3 ) ) 

C.................  Look for plant/char/scc assignments at various levels
                F6 = FINDC( CHK16, TXCNT( 16 ), CHRT16 )
                F5 = FINDC( CHK15, TXCNT( 15 ), CHRT15 )
                F4 = FINDC( CHK14, TXCNT( 14 ), CHRT14 )
                F3 = FINDC( CHK13, TXCNT( 13 ), CHRT13 )
                F2 = FINDC( CHK12, TXCNT( 12 ), CHRT12 )

C.................  Look for plant/char assignments at various levels
                IF( F6 .LE. 0 ) F6 = FINDC( CSRC5, TXCNT( 16 ), CHRT16 )
                IF( F5 .LE. 0 ) F5 = FINDC( CSRC4, TXCNT( 15 ), CHRT15 ) 
                IF( F4 .LE. 0 ) F4 = FINDC( CSRC3, TXCNT( 14 ), CHRT14 ) 
                IF( F3 .LE. 0 ) F3 = FINDC( CSRC2, TXCNT( 13 ), CHRT13 ) 
                IF( F2 .LE. 0 ) F2 = FINDC( CSRC1, TXCNT( 12 ), CHRT12 ) 

C.................  Look for plant/SCC and plant assignments
                F1 = FINDC( CHK11, TXCNT( 11 ), CHRT11 ) 
                F0 = FINDC( CHK10, TXCNT( 10 ), CHRT10 )

C.................  Choose the most specific assignment first
                IF( F6 .GT. 0 ) THEN
                    LNONHAP( S ) = .FALSE.
                    CYCLE                       !  to end of sources-loop
                ELSE IF( F5 .GT. 0 ) THEN
                    LNONHAP( S ) = .FALSE.
                    CYCLE                       !  to end of sources-loop
                ELSE IF( F4 .GT. 0 ) THEN
                    LNONHAP( S ) = .FALSE.
                    CYCLE                       !  to end of sources-loop
                ELSE IF( F3 .GT. 0 ) THEN
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
                  
            END IF

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

C.........  Reverse non-HAP exclusions to inclusions if the header is
C           defined as '/INCLUDE/' in NHAPEXCLUDE file

        IF( .NOT. NHAP_EXCL ) THEN
            DO S = 1, NSRC
                IF( LNONHAP( S ) ) THEN
                    LNONHAP( S ) = .FALSE.   ! non-HAP exlusions
                ELSE
                    LNONHAP( S ) = .TRUE.    ! non-HAP inclusions
                ENDIF
            END DO
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ASGNNHAPX
