
        SUBROUTINE PREPLM( EMLAYS, HMIX, HTS, PSFC, TS, DDZH, QV, TA, 
     &                     UW, VW, ZH, ZF, ZSTK, PRES, LSTK, LPBL, TSTK, 
     &                     WSTK, DTHDZ, WSPD, ZZF )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C    Computes the values needed for the PLMRIS subroutine from the 
C    meteorology data.
C
C  PRECONDITIONS REQUIRED:
C    Interpolated (to the location of a source) meteorology data as input,
C    vertical grid structure.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       I/O API 
C
C  REVISION  HISTORY:
C       Copied from preplm.f v 1.2 in DAQM-V2 Emissions Preprocessor by
C           M. Houyoux 3/99
C
C***********************************************************************
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
C***********************************************************************
 
        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    ! physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        REAL          POLY

        EXTERNAL      POLY

C...........   SUBROUTINE ARGUMENTS (NOTE: All met parms are per-source)
        INTEGER, INTENT (IN) :: EMLAYS          ! no. emissions layers
        REAL   , INTENT (IN) :: HMIX            ! mixing height
        REAL   , INTENT (IN) :: HTS             ! stack height
        REAL   , INTENT (IN) :: PSFC            ! surface pressure
        REAL   , INTENT (IN) :: TS              ! surface temperature
        REAL   , INTENT (IN) :: DDZH ( EMLAYS ) ! 1/( zh(l) - zh(l-1) )
        REAL   , INTENT (IN) :: QV   ( EMLAYS ) ! mixing ratio
        REAL   , INTENT (IN) :: TA   ( EMLAYS ) ! absolute temperature
        REAL   , INTENT (IN) :: UW   ( EMLAYS ) ! x-direction winds
        REAL   , INTENT (IN) :: VW   ( EMLAYS ) ! y-direction winds
        REAL   , INTENT (IN) :: ZF   ( EMLAYS ) ! layer surface height (m)
        REAL   , INTENT (IN) :: ZH   ( EMLAYS ) ! layer center  height (m)
        REAL   , INTENT (IN) :: ZSTK ( EMLAYS ) ! zf( l,s ) - stkht(s) (m)
        REAL, INTENT(IN OUT) :: PRES ( EMLAYS ) ! pessure
        INTEGER, INTENT(OUT) :: LSTK            ! first L: ZF(L) > STKHT
        INTEGER, INTENT(OUT) :: LPBL            ! first L: ZF(L) > mixing layer
        REAL   , INTENT(OUT) :: TSTK            ! tmptr at top of stack (K)
        REAL   , INTENT(OUT) :: WSTK            ! wind speed @ top of stack(m/s)
        REAL   , INTENT(OUT) :: DTHDZ( EMLAYS ) ! gradient of THETV
        REAL   , INTENT(OUT) :: WSPD ( EMLAYS ) ! wind speed (m/s)
        REAL   , INTENT(OUT) :: ZZF  ( 0:EMLAYS )! elevation at full-levels

C...........   Local variables

        INTEGER       L, M

        REAL          ES
        REAL          QSFC
        REAL          TVSFC
        REAL          THETG
        REAL          THV1
        REAL          THVK
        REAL          TV( EMLAYS )     ! virtual temperature
        REAL          P, Q

C***********************************************************************
C   begin body of subroutine PREPLM

C...........   Convert pressure to millibars from pascals, compute wind speed,
C              and compute virtual temperature
 
        DO L = 1, EMLAYS
 
            PRES( L ) = 1.0E-2 * PRES( L )
 
            P = UW( L )
            Q = VW( L )
            WSPD( L ) = SQRT( P * P  +  Q * Q )
            TV  ( L ) = TA( L ) * 
     &                  ( 1. + 0.622 * ( QV( L ) / ( 1. + QV( L ) ) ) )
 
        END DO

        ES    = 6.1078 * EXP( 5384.21 / CTOK - 5384.21 / TS )
        QSFC  = 0.622  * ES / ( PSFC - ES )
        TVSFC = TS   * ( 1.0 + 0.6077 * QSFC )
        THETG = TVSFC  * ( 1000.0 / PSFC )**0.286
        THV1  = TV( 1 )*( 1000.0 / PRES( 1 ) )**0.286
        DTHDZ( 1 ) = ( THV1 - THETG ) / ZH( 1 )

        IF ( HMIX .LE. ZF( 1 ) ) THEN
            LPBL = 1
        END IF
        IF ( HTS .LE. ZF( 1 ) ) THEN
            LSTK = 1
        END IF

        ZZF( 0 ) = 0.0
        ZZF( 1 ) = ZF( 1 )
        DO L = 2, EMLAYS
 
            IF ( HMIX  .GT. ZF( L-1 ) )  LPBL = L
            IF ( HTS .GT. ZF( L-1 ) )  LSTK = L
 
            THVK = TV( L ) * ( 1000.0 / PRES( L ) )**0.286
            DTHDZ( L ) = DDZH( L-1 ) * ( THVK - THV1 )
            THV1 = THVK
 
            ZZF( L ) = ZF( L )
 
        END DO
 
        M    = MAX( 1, LSTK - 2 )
        TSTK =      POLY( HTS, ZH( M ), TA( M ), 3 )
        WSTK = MAX( POLY( HTS, ZH( M ), WSPD( M ), 3 ), 0.1 )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

	END SUBROUTINE PREPLM
