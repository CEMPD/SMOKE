
        SUBROUTINE MRGMULT( NSRC, NG, NL, NMAT1, NMAT2, 
     &                      KEY1, KEY2, KEY4, ISPC, FG, FT, 
     &                      EMSRC, RINFO, CUMATX, SMATX,
     &                      NX, IX, GMATX, ICNY, GOUT1, GOUT2,
     &                      COUT1, COUT2, COUT4, COUT5 )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine multiplies a source-emissions vector with a gridding
C      matrix and optionally a speciation array and multiplicative control
C      array.  Which matrices are applied depend on the setting of the 
C      keys in the subroutine call.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: NPSRC, ELEVFLAG, ELEVADJ, LFRAC, PINGFLAG,
     &       INLINEFLAG, SRCGRPFLAG, ISRCGRP, EMGGRD, SUBSECFLAG, IGRPNUM

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: ELEVFLTR

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY

        IMPLICIT NONE

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC        ! number of source
        INTEGER     , INTENT (IN) :: NG          ! (local) number of grid cells
        INTEGER     , INTENT (IN) :: NL          ! (local) number of layers
        INTEGER     , INTENT (IN) :: NMAT1       ! dim 1 for gridding matrix
        INTEGER     , INTENT (IN) :: NMAT2       ! dim 2 for gridding matrix
        INTEGER     , INTENT (IN) :: KEY1        ! inven emissions index
        INTEGER     , INTENT (IN) :: KEY2        ! mult controls index
        INTEGER     , INTENT (IN) :: KEY4        ! speciation matrix index
        INTEGER     , INTENT (IN) :: ISPC        ! output species index
        REAL        , INTENT (IN) :: FG          ! gridded output units conv
        REAL        , INTENT (IN) :: FT          ! st/co tot output units conv
        REAL        , INTENT (IN) :: EMSRC ( NSRC,* ) ! source-based emissions
        REAL        , INTENT (IN) :: RINFO ( NSRC,2 ) ! reactivity information
        REAL        , INTENT (IN) :: CUMATX( NSRC,* ) ! mult control factors
        REAL        , INTENT (IN) :: SMATX ( NSRC,* ) ! speciation factors
        INTEGER     , INTENT (IN) :: NX    ( NG )     ! no. of sources per cell
        INTEGER     , INTENT (IN) :: IX    ( NMAT1 )  ! list of sources per cell
        REAL        , INTENT (IN) :: GMATX ( NMAT2 )  ! gridding coefficients
        INTEGER     , INTENT (IN) :: ICNY  ( NSRC )   ! county index by source
        REAL        , INTENT(OUT) :: GOUT1 ( NG, NL ) ! one-time gridded emis
        REAL     , INTENT(IN OUT) :: GOUT2 ( NG, NL ) ! cumulative gridded emis
        REAL     , INTENT(IN OUT) :: COUT1 ( NCOUNTY, * )! no control county
        REAL     , INTENT(IN OUT) :: COUT2 ( NCOUNTY, * )! multiplv cntl county
        REAL     , INTENT(IN OUT) :: COUT4 ( NCOUNTY, * )! reactivity cntl cnty
        REAL     , INTENT(IN OUT) :: COUT5 ( NCOUNTY, * )! all control cntl cnty

C.........  Other local variables
        INTEGER         C, J, K, L, S   ! counters and indicies
        INTEGER         IDX             ! index to list of counties in grid 
        INTEGER         GIDX            ! index to source group
        REAL            GFAC            ! tmp gridding factor
        REAL            FG0             ! gridding conv fac div. totals conv fac
        REAL(8)         SUM1            ! sum for GOUT1   
        REAL(8)         SUM2            ! sum for GOUT2
        REAL(8)         MULT            ! tmp value with multiplictv controls
        REAL(8)         REAC            ! tmp value with reactivity controls
        REAL(8)         VAL             ! tmp value  
        REAL(8)         VMP             ! tmp market penetration value  

        CHARACTER(16) :: PROGNAME = 'MRGMULT' ! program name

C***********************************************************************
C   begin body of subroutine MRGMULT

        FG0 = FG / FT

C.........  If the sources are point sources and there are elevated sources,
C           transfer ELEVFLTR to ELEVADJ
        IF( NSRC .EQ. NPSRC .AND. (ELEVFLAG .OR. PINGFLAG .OR. INLINEFLAG) ) THEN

            ELEVADJ( 1:NPSRC ) = ELEVFLTR( 1:NPSRC )

        ELSE

            ELEVADJ = 0.

        END IF

C.........  Check if this is a valid inventory pollutant for this call, and
C           if the number of layers is one.
        IF( NL .LE. 1 .AND. KEY1 .GT. 0 ) THEN

C............. If multiplicative controls & speciation
            IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,1 )
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX  = ICNY( S )
                        GFAC = GMATX( K ) * FT

                        VAL = EMSRC( S,KEY1 )* SMATX( S,KEY4 )* GFAC
                        COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL

                        MULT = VAL * CUMATX( S,KEY2 )
                        COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + MULT

                        VMP  = RINFO( S,2 )
                        REAC = ( VAL * (1.-VMP) + RINFO( S,1 ) * VMP ) *
     &                         GFAC
                        COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + REAC
              
                        VAL = MULT * (1.-VMP) + RINFO( S,1 )* VMP* GFAC

                        COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0
                        
                        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                            GIDX = ISRCGRP( S )
                            IF( SUBSECFLAG ) GIDX = IGRPNUM( ISRCGRP( S ) )
                            EMGGRD( C,GIDX ) = 
     &                          EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C............. If multiplicative controls only
            ELSE IF( KEY2 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,1 )
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX  = ICNY( S )
                        GFAC = GMATX( K ) * FT

                        VAL = EMSRC( S,KEY1 )*GFAC
                        COUT1( IDX,KEY1 ) = COUT1( IDX,KEY1 ) + VAL

                        MULT = VAL * CUMATX( S,KEY2 )
                        COUT2( IDX,KEY1 ) = COUT2( IDX,KEY1 ) + MULT
                        COUT4( IDX,KEY1 ) = COUT4( IDX,KEY1 ) + VAL
                        COUT5( IDX,KEY1 ) = COUT5( IDX,KEY1 ) + MULT

                        VAL  = MULT
                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0
                        
                        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                            GIDX = ISRCGRP( S )
                            IF( SUBSECFLAG ) GIDX = IGRPNUM( ISRCGRP( S ) )
                            EMGGRD( C,GIDX ) = 
     &                          EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C.............  If speciation only
            ELSE IF( KEY4 .GT. 0 ) THEN

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,1 )
                    SUM2 = GOUT2( C,1 )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX  = ICNY( S )
                        GFAC = GMATX( K ) * FT

                        VAL = EMSRC( S,KEY1 ) * SMATX( S,KEY4 ) * GFAC
                        COUT1( IDX,ISPC )= COUT1( IDX,ISPC ) + VAL
                        COUT2( IDX,ISPC )= COUT2( IDX,ISPC ) + VAL

                        VMP = RINFO( S,2 )
                        VAL = VAL * (1.-VMP) + RINFO( S,1 ) * VMP * GFAC

                        COUT4( IDX,ISPC )= COUT4( IDX,ISPC ) + VAL
                        COUT5( IDX,ISPC )= COUT5( IDX,ISPC ) + VAL

                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0
                        
                        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                            GIDX = ISRCGRP( S )
                            IF( SUBSECFLAG ) GIDX = IGRPNUM( ISRCGRP( S ) )
                            EMGGRD( C,GIDX ) = 
     &                          EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

C.............  If inventory pollutant only
            ELSE 
                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,1 )
                    SUM2 = GOUT2( C,1 )
    
                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX = ICNY( S )

                        VAL = EMSRC( S,KEY1 ) * GMATX( K ) * FT
                        COUT1( IDX,KEY1 )= COUT1( IDX,KEY1 ) + VAL
                        COUT2( IDX,KEY1 )= COUT2( IDX,KEY1 ) + VAL
                        COUT4( IDX,KEY1 )= COUT4( IDX,KEY1 ) + VAL
                        COUT5( IDX,KEY1 )= COUT5( IDX,KEY1 ) + VAL

                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0
                        
                        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                            GIDX = ISRCGRP( S )
                            IF( SUBSECFLAG ) GIDX = IGRPNUM( ISRCGRP( S ) )
                            EMGGRD( C,GIDX ) = 
     &                          EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END DO

                    GOUT1( C,1 ) = SUM1
                    GOUT2( C,1 ) = SUM2

                END DO

            END IF  ! End which of controls and speciation

C.........  If we need to use layer fractions...

        ELSE IF( KEY1 .GT. 0 ) THEN

C............. If multiplicative controls & speciation & layer fractions
            IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = GOUT1( C,L )
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            IDX  = ICNY( S )
                            GFAC = GMATX( K ) * LFRAC( S,L ) * FT

                            VAL = EMSRC( S,KEY1 )* SMATX( S,KEY4 )* GFAC
                            COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL

                            MULT = VAL * CUMATX( S,KEY2 )
                            COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + MULT

                            VMP  = RINFO( S,2 )
                            REAC = ( VAL * (1.-VMP) + RINFO(S,1)* VMP )*
     &                             GFAC
                            COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + REAC

                            VAL = MULT* (1.-VMP) + RINFO(S,1)* VMP* GFAC

                            COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

                            VAL  = VAL * ( 1.-ELEVADJ( S ) )
                            SUM1 = SUM1 + VAL * FG0
                            SUM2 = SUM2 + VAL * FG0

                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If multiplicative controls and layer fractions
            ELSE IF( KEY2 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = GOUT1( C,L )
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            IDX  = ICNY( S )
                            GFAC = GMATX( K ) * LFRAC( S,L ) * FT

                            VAL = EMSRC( S,KEY1 )*GFAC
                            COUT1( IDX,KEY1 ) = COUT1( IDX,KEY1 ) + VAL

                            MULT = VAL * CUMATX( S,KEY2 )
                            COUT2( IDX,KEY1 ) = COUT2( IDX,KEY1 ) + MULT
                            COUT4( IDX,KEY1 ) = COUT4( IDX,KEY1 ) + VAL
                            COUT5( IDX,KEY1 ) = COUT5( IDX,KEY1 ) + MULT

                            VAL  = MULT
                            VAL  = VAL * ( 1.-ELEVADJ( S ) )
                            SUM1 = SUM1 + VAL * FG0
                            SUM2 = SUM2 + VAL * FG0

                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

C............. If speciation and layer fraction
            ELSE IF( KEY4 .GT. 0 ) THEN

                DO L = 1, NL

                    K = 0
                    DO C = 1, NG

                        SUM1 = GOUT1( C,L )
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            IDX  = ICNY( S )
                            GFAC = GMATX( K ) * LFRAC( S,L ) * FT

                            VAL = EMSRC( S,KEY1 )* SMATX( S,KEY4 )* GFAC
                            COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL
                            COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + VAL

                            VMP  = RINFO( S,2 )
                            VAL  = ( VAL*(1.-VMP) + 
     &                               RINFO( S,1 ) * VMP * GFAC )

                            COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + VAL
                            COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

                            VAL  = VAL * ( 1.-ELEVADJ( S ) )
                            SUM1 = SUM1 + VAL * FG0
                            SUM2 = SUM2 + VAL * FG0
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

                        SUM1 = GOUT1( C,L )
                        SUM2 = GOUT2( C,L )

                        DO J = 1, NX( C )
                            K = K + 1
                            S = IX( K )
                            IDX = ICNY( S )

                            VAL = LFRAC( S,L ) *
     &                            EMSRC ( S,KEY1 ) * GMATX( K ) * FT
                            COUT1( IDX,KEY1 ) = COUT1( IDX,KEY1 ) + VAL
                            COUT2( IDX,KEY1 ) = COUT2( IDX,KEY1 ) + VAL
                            COUT4( IDX,KEY1 ) = COUT4( IDX,KEY1 ) + VAL
                            COUT5( IDX,KEY1 ) = COUT5( IDX,KEY1 ) + VAL

                            VAL  = VAL * ( 1.-ELEVADJ( S ) )
                            SUM1 = SUM1 + VAL * FG0
                            SUM2 = SUM2 + VAL * FG0
                        END DO

                        GOUT1( C,L ) = SUM1
                        GOUT2( C,L ) = SUM2

                    END DO
                END DO

            END IF  ! End which of controls and speciation

        END IF      ! End if no inventory emissions or L > 1

        RETURN

        END SUBROUTINE MRGMULT
