
        SUBROUTINE UNGRIDBV( NC, NR, LAT, LON, 
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
        REAL,         INTENT (IN) :: LAT( NC,NR ) ! grid cell center latitudes
        REAL,         INTENT (IN) :: LON( NC,NR ) ! grid cell center longitudes
        INTEGER,      INTENT (IN) :: NPTS         ! number of point-source locations
        REAL,         INTENT (IN) :: XLOC( NPTS ) ! X point coordinates
        REAL,         INTENT (IN) :: YLOC( NPTS ) ! Y point coordinates
        INTEGER,      INTENT(OUT) :: NU( 4,NPTS ) ! surrounding grid cells for each source
        REAL,         INTENT(OUT) :: CU( 4,NPTS ) ! fraction of each grid cell
        
C.........  Local variables
        INTEGER          I, K, S     ! counters and indices
        INTEGER          COL         ! column for current point
        INTEGER          ROW         ! row for current point
        
        REAL             X, Y        ! normalized difference between point location
                                     !   and grid cell center
        REAL             P, Q        ! fractions used for calculating coefficients

        CHARACTER(16) :: PROGNAME = 'UNGRIDBV' ! program name

C***********************************************************************
C   begin body of subroutine UNGRIDBV

C.........  Loop through point-source locations
        DO S = 1, NPTS

C.............  Find column
            DO I = 1, NC
            
C.................  Check if point is to the left of the grid
                IF( I == 1 .AND. XLOC( S ) < LON( I,1 ) ) THEN
                    COL = 0
                    X = 0.
                    EXIT
                END IF
               
C.................  Check if point is to the right of the grid
                IF( I == NC .AND. XLOC( S ) > LON( I,1 ) ) THEN
                    COL = NC
                    X = 1.
                    EXIT
                END IF

C.................  Check if point is between current and next location               
                IF( XLOC( S ) >= LON( I,1 ) .AND. 
     &              XLOC( S ) <= LON( I+1,1 ) ) THEN
                    COL = I
                    X = ( XLOC( S ) - LON( I,1 ) ) / 
     &                  ( LON( I+1,1 ) - LON( I,1 ) )
                    EXIT
                END IF            
            END DO
            
C.............  Find row
            DO I = 1, NR
                
C.................  Check if point is below the grid
                IF( I == 1 .AND. YLOC( S ) < LAT( 1,I ) ) THEN
                    ROW = 0
                    Y = 0.
                    EXIT
                END IF
                
C.................  Check if point is above the grid
                IF( I == NR .AND. YLOC( S ) > LAT( 1,I ) ) THEN
                    ROW = NR
                    Y = 1.
                    EXIT
                END IF

C.................  Check if point is between current and next locations                
                IF( YLOC( S ) >= LAT( 1,I ) .AND. 
     &              YLOC( S ) <= LAT( 1,I+1 ) ) THEN
                    ROW = I
                    Y = ( YLOC( S ) - LAT( 1,I ) ) /
     &                  ( LAT( 1,I+1 ) - LAT( 1,I ) )
                    EXIT
                END IF
            END DO

C.............  Set grid cells surrounding point location
            IF( ROW == 0 ) THEN               ! row below grid
                IF( COL == 0 ) THEN           ! column to the left of grid
                    K = 1
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                ELSE IF( COL == NC ) THEN     ! column to the right of grid
                    K = NC
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                ELSE                          ! column inside grid
                    K = COL
                    NU( 1,S ) = K
                    NU( 2,S ) = K + 1
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                END IF
            ELSE IF( ROW == NR ) THEN         ! row above grid
                IF( COL == 0 ) THEN           ! column to the left of grid
                    K = ( NR - 1 ) * NC + 1
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                ELSE IF( COL == NC ) THEN     ! column to the right of grid
                    K = NR * NC
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                ELSE                          ! column inside grid
                    K = ( NR - 1 ) * NC + COL
                    NU( 1,S ) = K
                    NU( 2,S ) = K + 1
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                END IF
            ELSE                              ! row inside grid
                IF( COL == 0 ) THEN           ! column to the left of grid
                    K = ( ROW - 1 ) * NC + 1
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K + NC
                    NU( 4,S ) = K + NC
                ELSE IF( COL == NC ) THEN     ! column to the right of grid
                    K = ROW * NC
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K + NC
                    NU( 4,S ) = K + NC
                ELSE                          ! column inside grid
                    K = ( ROW - 1 ) * NC + COL
                    NU( 1,S ) = K
                    NU( 2,S ) = K + 1
                    NU( 3,S ) = K + NC
                    NU( 4,S ) = K + NC + 1
                END IF
            END IF

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
