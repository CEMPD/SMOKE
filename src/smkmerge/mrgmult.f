
        SUBROUTINE MRGMULT( NSRC, NG, NL, NMAT1, NMAT2, 
     &                      KEY1, KEY2, KEY3, KEY4, EMSRC, 
     &                      CUMATX, CAMATX, SMATX,
     &                      NX, IX, GMATX, GOUT1, GOUT2 )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine multiplies a source-emissions vector with a gridding
C      matrix and optionally a speciation array and multiplicative control
C      array. An additive control array canbe added to the emissions.  Which
C      matrices are applied depend on the setting of the keys in the 
C      subroutine call.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC        ! number of source
        INTEGER     , INTENT (IN) :: NG          ! (local) number of grid cells
        INTEGER     , INTENT (IN) :: NL          ! (local) number of layers
        INTEGER     , INTENT (IN) :: NMAT1       ! dim 1 for gridding matrix
        INTEGER     , INTENT (IN) :: NMAT2       ! dim 2 for gridding matrix
        INTEGER     , INTENT (IN) :: KEY1        ! inven emissions index
        INTEGER     , INTENT (IN) :: KEY2        ! mult controls index
        INTEGER     , INTENT (IN) :: KEY3        ! additive controls index
        INTEGER     , INTENT (IN) :: KEY4        ! speciation index
        REAL        , INTENT (IN) :: EMSRC ( NSRC,* ) ! source-based emissions
        REAL        , INTENT (IN) :: CUMATX( NSRC,* ) ! mult control factors
        REAL        , INTENT (IN) :: CAMATX( NSRC,* ) ! additive control factors
        REAL        , INTENT (IN) :: SMATX ( NSRC,* ) ! speciation factors
        INTEGER     , INTENT (IN) :: NX    ( NG )     ! no. of sources per cell
        INTEGER     , INTENT (IN) :: IX    ( NMAT1 )  ! list of sources per cell
        REAL        , INTENT (IN) :: GMATX ( NMAT2 )  ! gridding coefficients
        REAL        , INTENT(OUT) :: GOUT1 ( NG, NL ) ! one-time gridded emis
        REAL     , INTENT(IN OUT) :: GOUT2 ( NG, NL ) ! cumulative gridded emis

C.........  Other local variables
        INTEGER         C, J, K, L, S   ! counters and indicies   
        REAL*8          SUM1            ! sum for GOUT1   
        REAL*8          SUM2            ! sum for GOUT2  
        REAL*8          VAL             ! tmp value  

        CHARACTER*16 :: PROGNAME = 'MRGMULT' ! program name

C***********************************************************************
C   begin body of subroutine MRGMULT

C.........  Check if this is a valid inventory pollutant for this call, and
C           if the number of layers is one.
        IF( NL .EQ. 1 .AND. KEY1 .GT. 0 ) THEN

C............. If multiplicative controls, additive controls, and speciation
            IF( KEY2 .GT. 0 .AND. KEY3 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = CAMATX( S,KEY3 ) + 
     &                        EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                        CUMATX( S,KEY2 ) * SMATX( S,KEY4 )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C............. If multiplicative controls & additive controls
            ELSE IF( KEY2 .GT. 0 .AND. KEY3 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = CAMATX( S,KEY3 ) + 
     &                        EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                        CUMATX( S,KEY2 )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C............. If multiplicative controls & speciation
            ELSE IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                        CUMATX( S,KEY2 ) * SMATX( S,KEY4 )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C............. If additive controls & speciation
            ELSE IF( KEY3 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = CAMATX( S,KEY3 ) + 
     &                        EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                        SMATX ( S,KEY4 )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C............. If multiplicative controls only
            ELSE IF( KEY2 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                        CUMATX( S,KEY2 )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C............. If additive controls only
            ELSE IF( KEY3 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = CAMATX( S,KEY3 ) + 
     &                        EMSRC ( S,KEY1 ) * GMATX( K )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C.............  If speciation only
            ELSE IF( KEY4 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = EMSRC( S,KEY1 ) * GMATX( K ) * 
     &                        SMATX( S,KEY4 )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO


C.............  If inventory pollutant only
            ELSE 
                K = 0
                DO C = 1, NG

                    SUM1 = 0.
                    SUM2 = GOUT2( C,1 )
    
                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        VAL = EMSRC( S,KEY1 ) * GMATX( K )
                        SUM1 = SUM1 + VAL
                        SUM2 = SUM2 + VAL
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

            END IF  ! End which of controls and speciation

C.........  If we need to use layer fractions...

        ELSE IF( NL .GT. 1 .AND. KEY1 .GT. 0 ) THEN

C............. If multiplicative controls, additive controls, and speciation
            IF( KEY2 .GT. 0 .AND. KEY3 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = CAMATX( S,KEY3 ) + LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                            CUMATX( S,KEY2 ) * SMATX( S,KEY4 )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If multiplicative controls & additive controls
            ELSE IF( KEY2 .GT. 0 .AND. KEY3 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = CAMATX( S,KEY3 ) + LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                            CUMATX( S,KEY2 )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If multiplicative controls & speciation
            ELSE IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                            CUMATX( S,KEY2 ) * SMATX( S,KEY4 )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If additive controls & speciation
            ELSE IF( KEY3 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = CAMATX( S,KEY3 ) + LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                            SMATX( S,KEY4 )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If multiplicative controls only
            ELSE IF( KEY2 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                            CUMATX( S,KEY2 )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If additive controls only
            ELSE IF( KEY3 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = CAMATX( S,KEY3 ) + LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If speciation only
            ELSE IF( KEY4 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * 
     &                            SMATX( S,KEY4 )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C.............  If inventory pollutant and layer fractions
            ELSE 

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = 0.
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            VAL = LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K )
                            SUM1 = SUM1 + VAL
                            SUM2 = SUM2 + VAL
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

            END IF  ! End which of controls and speciation

        END IF      ! End if no inventory emissions or L > 1

        RETURN

        END SUBROUTINE MRGMULT
