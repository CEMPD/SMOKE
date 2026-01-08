
        SUBROUTINE APPLUMAT( NSRC, NMATX, VAL, NU, IU, CU, VALBYSRC )

C***********************************************************************
C  subroutine APPLUMAT body starts at line 72
C
C  DESCRIPTION:
C      Applies the "ungridding" matrix to gridded data to compute a per-source
C      value of the data.  If the ungridding matrix has no factors for a source,
C      a missing value is returned for that source (i.e., BADVAL3)
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************

C.........  This module contains the global variables for the 3-d grid
        USE M3UTILIO

        USE MODGRID, ONLY: NCOLS, YOFF, XDIFF, XOFF

        IMPLICIT NONE

C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC          ! no. sources
        INTEGER     , INTENT (IN) :: NMATX         ! ugridding factors dimension
        REAL        , INTENT (IN) :: VAL( * )      ! gridded data
        INTEGER     , INTENT (IN) :: NU( NSRC  )   ! no. cells per source
        INTEGER     , INTENT (IN) :: IU( NMATX )   ! grid cell IDs 
        REAL        , INTENT (IN) :: CU( NMATX )   ! ungridding coefficients
        REAL        , INTENT(OUT) :: VALBYSRC( NSRC ) ! ungridded data

C...........   Other local variables
        INTEGER     C, J, K, S     ! counters and indices
        INTEGER     COL            ! subgrid column number
        INTEGER     ROW            ! subgrid row number

        REAL        RDUM        ! tmp value for summing over cells

        CHARACTER(16) :: PROGNAME = 'APPLUMAT' ! program name

C***********************************************************************
C   begin body of subroutine APPLUMAT

C.........  Apply ungridding matrix from a (possible) subgrid to data on base 
C           grid.  If no subgrid, then XOFF and YOFF will be 1 and no problem.
        K = 0
        DO S = 1, NSRC

            IF( NU( S ) .GT. 0 ) THEN
                RDUM = 0.0
            ELSE
                RDUM = BADVAL3
            END IF

            DO J = 1, NU( S )
                K = K + 1

C.................  Get column and row from subgrid
                C = IU( K )
                ROW = C / NCOLS          ! note: integer math
                IF( MOD( C, NCOLS ) .GT. 0. ) ROW = ROW + 1
                COL = C - ( ROW-1 ) * NCOLS

C.................  Compute cell number of base grid
                C = ( ROW + YOFF - 1 ) * ( NCOLS + XDIFF ) + COL + XOFF

                RDUM = RDUM + VAL( C ) * CU( K )
            END DO

            VALBYSRC( S ) = RDUM

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE APPLUMAT
