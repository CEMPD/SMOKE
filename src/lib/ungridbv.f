
        SUBROUTINE UNGRIDBV( NC, NR, XREFS, YREFS, 
     &                       NPTS, XLOC, YLOC, NU, CU )

C***************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine computes "ungridding" matrices to be used by the
C      I/O API routine BMATVEC() based on the grid cell coordinates in
C      a variable grid GRIDCRO2D file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen 7/04 based on ungridb.f
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
C***************************************************************************

        IMPLICIT NONE
        
C.........  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: NC           ! number of columns in grid
        INTEGER,      INTENT (IN) :: NR           ! number of rows in grid
        REAL,         INTENT (IN) :: XREFS( NC,NR ) ! grid cell center x coordinates
        REAL,         INTENT (IN) :: YREFS( NC,NR ) ! grid cell center y coordinates
        INTEGER,      INTENT (IN) :: NPTS         ! number of point-source locations
        REAL*8,       INTENT (IN) :: XLOC( NPTS ) ! X point coordinates
        REAL*8,       INTENT (IN) :: YLOC( NPTS ) ! Y point coordinates
        INTEGER,      INTENT(OUT) :: NU( 4,NPTS ) ! surrounding grid cells for each source
        REAL,         INTENT(OUT) :: CU( 4,NPTS ) ! fraction of each grid cell
        
C.........  Local variables
        INTEGER          I, K, S     ! counters and indices
        INTEGER          COL         ! column for current point
        INTEGER          ROW         ! row for current point
        
        REAL             LEFT        ! left edge of matching grid cell
        REAL             RIGHT       ! right edge of matching grid cell
        REAL             BOTTOM      ! bottom edge of matching grid cell
        REAL             TOP         ! top edge of matching grid cell
        REAL             X, Y        ! normalized difference between point location
                                     !   and grid cell center
        REAL             P, Q        ! fractions used for calculating coefficients

        CHARACTER(16) :: PROGNAME = 'UNGRIDBV' ! program name

C***********************************************************************
C   begin body of subroutine UNGRIDBV

C.........  Loop through point-source locations
        DO S = 1, NPTS

C.............  Find column

C.............  Check if point is to the left of the grid
            IF( XLOC( S ) < XREFS( 1,1 ) ) THEN
                LEFT = 1
                RIGHT = 1
                
                X = 0.
                
C.............  Check if point is to the right of the grid
            ELSE IF( XLOC( S ) > XREFS( NC, 1 ) ) THEN
                LEFT = NC
                RIGHT = NC
                
                X = 0.

C.............  Loop through columns to find the correct one                
            ELSE
                DO I = 2, NC
                    IF( XLOC( S ) < XREFS( I,1 ) ) THEN
                        LEFT = I - 1
                        RIGHT = I
                        
                        X = ( XLOC( S ) - XREFS( I-1,1 ) ) /
     &                      ( XREFS( I,1 ) - XREFS( I-1,1 ) )
                        EXIT
                    END IF
                END DO
            END IF
            
C.............  Find row

C.............  Check if point is below the grid
            IF( YLOC( S ) < YREFS( 1,1 ) ) THEN
                BOTTOM = 1
                TOP = 1
                
                Y = 0.
            
C.............  Check if point is above the grid
            ELSE IF( YLOC( S ) > YREFS( 1,NR ) ) THEN
                BOTTOM = NR
                TOP = NR
                
                Y = 0.

C.............  Loop through rows to find the correct one
            ELSE
                DO I = 2, NR
                    IF( YLOC( S ) < YREFS( 1,I ) ) THEN
                        BOTTOM = I - 1
                        TOP = I
                        
                        Y = ( YLOC( S ) - YREFS( 1,I-1 ) ) /
     &                      ( YREFS( 1,I ) - YREFS( 1,I-1 ) )
                        EXIT
                    END IF
                END DO
            END IF

C.............  Set grid cells surrounding point location
            NU( 1,S ) = ( BOTTOM - 1 ) * NC + LEFT
            NU( 2,S ) = ( BOTTOM - 1 ) * NC + RIGHT
            NU( 3,S ) = ( TOP - 1 ) * NC + LEFT
            NU( 4,S ) = ( TOP - 1 ) * NC + RIGHT

C.............  Calculate fractions for each surrounding grid cell
            P = 1. - X
            Q = 1. - Y
            CU( 1,S ) = P * Q
            CU( 2,S ) = X * Q
            CU( 3,S ) = P * Y
            CU( 4,S ) = X * Y
        
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE UNGRIDBV
