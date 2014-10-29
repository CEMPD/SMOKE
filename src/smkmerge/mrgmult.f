
        SUBROUTINE MRGMULTP( NSRC, NG, NL, NMAT1, NMAT2,
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
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       New optimized, parallel version Carlie J. Coats, Jr., 2014
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
C.........  MODMERGE: major data structure and control flags
C.........  MODELEV: arrays for plume-in-grid and major sources
C.........  MODSTCY: arrays for state and county summaries

        USE MODMERGE, ONLY: NPSRC, ELEVFLAG, ELEVADJ, LFRAC, PINGFLAG,
     &                      INLINEFLAG, SRCGRPFLAG, ISRCGRP, EMGGRD

        USE MODELEV, ONLY: ELEVFLTR

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
        INTEGER     , INTENT (IN) :: NX    ( 0:NG )   ! cumulative no. of sources per cell
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
        REAL            MULT            ! tmp value with multiplictv controls
        REAL            REAC            ! tmp value with reactivity controls
        REAL            REML            ! tmp value with reac&mult  controls
        REAL            VAL             ! tmp value
        REAL            VTMP            ! tmp market penetration value
        REAL            VBAR

        REAL            AVAL( NSRC )            ! tmp value
    
        REAL(8)         SUM             ! sum for GOUT*

        CHARACTER(16), PARAMETER :: PROGNAME = 'MRGMULT' ! program name

C***********************************************************************
C   begin body of subroutine MRGMULT

        IF( KEY1 .LE. 0 .OR. NL .LE. 0 ) RETURN

C.........  If the sources are point sources and there are elevated sources,
C           transfer ELEVFLTR to ELEVADJ

        IF ( NSRC .EQ. NPSRC .AND. 
     &       ( ELEVFLAG .OR. PINGFLAG .OR. INLINEFLAG ) ) THEN

            ELEVADJ( 1:NPSRC ) = ELEVFLTR( 1:NPSRC )

        ELSE

            ELEVADJ( 1:NPSRC ) = 0.0

        END IF

C.........  Check if this is a valid inventory pollutant for this call
C............. Per-source computations:
C............. If multiplicative controls & speciation

        IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

!$OMP       PARALLEL DO
!$OMP&         DEFAULT( NONE ),
!$OMP&          SHARED( NSRC, ICNY, EMSRC, SMATX, CUMATX, RINFO,
!$OMP&                  ELEVADJ, FT, ISPC, KEY1, KEY2, KEY4, 
!$OMP&                  AVAL, COUT1, COUT2, COUT4, COUT5 ),
!$OMP&         PRIVATE( S, IDX, VAL, MULT, VTMP, VBAR, REAC, REML )

            DO S = 1, NSRC

                IDX  = ICNY( S )
                VAL  = EMSRC( S,KEY1 ) * SMATX( S,KEY4 ) * FT
                MULT = VAL * CUMATX( S,KEY2 )
                VBAR = 1.0 - RINFO( S,2 )
                VTMP = RINFO( S,1 ) * RINFO( S,2 ) * FT
                REAC = VAL  * VBAR + VTMP
                REML = MULT * VBAR + VTMP

                AVAL( S ) = REML * ( 1.0 - ELEVADJ( S ) )

                COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL
                COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + MULT
                COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + REAC
                COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + REML

            END DO

C............. If multiplicative controls only

        ELSE IF( KEY2 .GT. 0 ) THEN

!$OMP       PARALLEL DO
!$OMP&         DEFAULT( NONE ),
!$OMP&          SHARED( NSRC, ICNY, EMSRC, CUMATX, FT,
!$OMP&                  ELEVADJ, ISPC, KEY1, KEY2, KEY4, 
!$OMP&                  AVAL, COUT1, COUT2, COUT4, COUT5 ),
!$OMP&         PRIVATE( S, IDX, VAL, MULT )

            DO S = 1, NSRC

                IDX  = ICNY( S )
                VAL  = EMSRC( S,KEY1 ) * FT
                MULT = VAL * CUMATX( S,KEY2 )

                AVAL( S ) = MULT * ( 1.0 - ELEVADJ( S ) )

                COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL
                COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + MULT
                COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + VAL
                COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + MULT

            END DO

C.............  If speciation only

        ELSE IF( KEY4 .GT. 0 ) THEN

!$OMP       PARALLEL DO
!$OMP&         DEFAULT( NONE ),
!$OMP&          SHARED( NSRC, ICNY, EMSRC, SMATX, FT,
!$OMP&                  ELEVADJ, ISPC, KEY1, KEY4, 
!$OMP&                  AVAL, COUT1, COUT2, COUT4, COUT5 ),
!$OMP&         PRIVATE( S, IDX, VAL )

            DO S = 1, NSRC

                IDX  = ICNY( S )
                VAL  = EMSRC( S,KEY1 ) * SMATX( S,KEY4 ) * FT

                AVAL( S ) = VAL * ( 1.0 - ELEVADJ( S ) )

                COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL
                COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + VAL
                COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + VAL
                COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

            END DO

C.............  If inventory pollutant only

        ELSE

!$OMP       PARALLEL DO
!$OMP&         DEFAULT( NONE ),
!$OMP&          SHARED( NSRC, ICNY, EMSRC, FT,
!$OMP&                  ELEVADJ, ISPC, KEY1,
!$OMP&                  AVAL, COUT1, COUT2, COUT4, COUT5 ),
!$OMP&         PRIVATE( S, IDX, VAL )

            DO S = 1, NSRC

                IDX  = ICNY( S )
                VAL  = EMSRC( S,KEY1 ) * FT

                AVAL( S ) = VAL * ( 1.0 - ELEVADJ( S ) )

                COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL
                COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + VAL
                COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + VAL
                COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

            END DO

        END IF  ! End which of controls and speciation

C............. Grid computations:  apply gridding matrix, etc.:

        FG0 = FG / FT

        IF ( NL .GT. 1  ) THEN

!$OMP       PARALLEL DO
!$OMP&        DEFAULT( NONE ),
!$OMP&         SHARED( NL, NG, NX, IX, GMATX, AVAL, LFRAC, FG0, 
!$OMP&                 GOUT1, GOUT2 ),
!$OMP&        PRIVATE( L, C, SUM, K, S, VAL )

            DO L = 1, NL
            DO C = 1, NG

                SUM = 0.0D0

                DO  K = NX( C-1 ) + 1, NX( C )
                    S   = IX( K )
                    VAL = GMATX( K )*AVAL( S )*LFRAC( S,L )
                    SUM = SUM + VAL

                END DO

                GOUT1( C,L ) = GOUT1( C,L ) + SUM * FG0
                GOUT2( C,L ) = GOUT1( C,L ) + SUM * FG0

            END DO
            END DO
            
        ELSE        !!  else nl=1

!$OMP       PARALLEL DO
!$OMP&        DEFAULT( NONE ),
!$OMP&         SHARED( NG, NX, IX, GMATX, AVAL, FG0, GOUT1, GOUT2 ),
!$OMP&        PRIVATE( C, SUM, K, S, VAL, GIDX )

            DO C = 1, NG

                SUM = 0.0D0

                DO K = NX( C-1 ) + 1, NX( C )
                    S   = IX( K )
                    VAL = GMATX( K ) * AVAL( S )
                    SUM = SUM + VAL
                END DO

                GOUT1( C,1 ) = GOUT1( C,1 ) + SUM * FG0
                GOUT2( C,1 ) = GOUT1( C,1 ) + SUM * FG0

            END DO

        END IF      ! end if  nl > 1 or not
        
        IF ( SRCGRPFLAG ) THEN

!$OMP       PARALLEL DO
!$OMP&        DEFAULT( NONE ),
!$OMP&         SHARED( NG, NX, IX, GMATX, AVAL, ISRCGRP, FG0, EMGGRD ),
!$OMP&        PRIVATE( C, K, S, VAL, GIDX )

            DO C = 1, NG

                DO K = NX( C-1 ) + 1, NX( C )
                    S    = IX( K )
                    VAL  = GMATX( K ) * AVAL( S ) * FG0
                    GIDX = ISRCGRP( S )
                    EMGGRD( C,GIDX ) = EMGGRD( C,GIDX ) + VAL
                END DO

            END DO

        END IF      ! end if srcgrpflag, or not

        RETURN

        END SUBROUTINE MRGMULTP
