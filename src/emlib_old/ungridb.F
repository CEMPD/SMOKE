
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  UNGRIDB( NCOLS, NROWS, XORIG, YORIG, XCELL, YCELL,
     &                       NPTS, XLOC, YLOC, NU, CU )

C***********************************************************************
C  subroutine body starts at line  51
C
C  FUNCTION:
C 	computes "ungridding" matrices to be used by BMATVEC()
C	for program LAYPOINT, etc., to perform bilinear interpolation
C	from a grid to a set of locations { <XLOC(S),YLOC(S)>, S=1:NPTS }
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C	prototype 12/95 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        INTEGER		NCOLS, NROWS	!  number of grid columns, rows
        REAL*8		XORIG, YORIG	!  X,Y coords of LL grid corner
        REAL*8		XCELL, YCELL	!  X,Y direction cell size
        INTEGER		NPTS	        !  number of (point-source) locations
        REAL		XLOC( NPTS ) 	!  X point coordinates
        REAL		YLOC( NPTS ) 	!  Y point coordinates
        INTEGER		NU( 4,NPTS )    !  single-indexed subscripts into grid
        REAL            CU( 4,NPTS )    !  coefficients

C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER		S		!  source counter
        INTEGER		C, R		!  subscripts into doubly-indexed grid
        INTEGER		K		!  subscript  into singly-indexed grid
        REAL		DDX, DDY	!  inverse cell size
        REAL		XD0, YD0	!  center of LL cell
        REAL		X, Y		!  grid-normal coords of point
        REAL		P, Q 		!  linear-interpolation coeffs


C***********************************************************************
C   begin body of subroutine  UNGRIDB

        DDX = 1.0 / XCELL
        DDY = 1.0 / YCELL
        
        XD0 = XORIG + 0.5D0 * XCELL
        YD0 = YORIG + 0.5D0 * YCELL

        DO  11  S = 1, NPTS
            
            X = DDX * ( XLOC( S ) - XD0 )	!  normalized grid coords
            Y = DDY * ( YLOC( S ) - YD0 )	!  normalized grid coords
            C = 1 + INT( X )        ! truncated to integer
            R = 1 + INT( Y )        ! truncated to integer
            X = MOD( X, 1.0 )   ! trapped between 0 and 1
            Y = MOD( Y, 1.0 )   ! trapped between 0 and 1

            IF ( R .LT. 1 ) THEN                  !  r below    grid

                IF ( C .LT. 1 ) THEN              !  c left of  grid

                    K = 1
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                    CU( 1,S ) = 1.0
                    CU( 2,S ) = 0.0
                    CU( 3,S ) = 0.0
                    CU( 4,S ) = 0.0

                ELSE IF ( C .GE. NCOLS - 1 ) THEN !  c right of grid

                    K = NCOLS
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                    CU( 1,S ) =  1.0
                    CU( 2,S ) =  0.0
                    CU( 3,S ) =  0.0
                    CU( 4,S ) =  0.0

                ELSE                              !  c "in" the grid

                    K = C
                    NU( 1,S ) = K
                    NU( 2,S ) = K + 1
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                    CU( 1,S ) = 1.0 - X
                    CU( 2,S ) = X
                    CU( 3,S ) = 0.0 
                    CU( 4,S ) = 0.0

                END IF

            ELSE IF ( R .GT. NROWS - 1 ) THEN     !  r above    grid

                IF ( C .LT. 1 ) THEN                !  c left of  grid

                    K = ( NROWS - 1 ) * NCOLS
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                    CU( 1,S ) = 1.0
                    CU( 2,S ) = 0.0
                    CU( 3,S ) = 0.0
                    CU( 4,S ) = 0.0

                ELSE IF ( C .GE. NCOLS - 1 ) THEN !  c right of grid

                    K = NROWS * NCOLS
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                    CU( 1,S ) = 1.0
                    CU( 2,S ) = 0.0
                    CU( 3,S ) = 0.0
                    CU( 4,S ) = 0.0

                ELSE                                !  c "in" the grid

                    K = ( NROWS - 1 ) * NCOLS  +  C
                    NU( 1,S ) = K
                    NU( 2,S ) = K + 1
                    NU( 3,S ) = K
                    NU( 4,S ) = K
                    CU( 1,S ) = 1.0 - X
                    CU( 2,S ) = X
                    CU( 3,S ) = 0.0
                    CU( 4,S ) = 0.0

                END IF

            ELSE                                  !  r "in" the grid

                IF ( C .LT. 1 ) THEN              !  c left of  grid

                    K = ( R - 1 ) * NCOLS
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K + NCOLS
                    NU( 4,S ) = K + NCOLS
                    CU( 1,S ) =  1.0 - Y
                    CU( 2,S ) =  0.0
                    CU( 3,S ) =  Y
                    CU( 4,S ) =  0.0

                ELSE IF ( C .GE. NCOLS - 1 ) THEN !  c right of grid

                    K = R * NCOLS
                    NU( 1,S ) = K
                    NU( 2,S ) = K
                    NU( 3,S ) = K + NCOLS
                    NU( 4,S ) = K + NCOLS
                    CU( 1,S ) =  1.0 - Y
                    CU( 2,S ) =  0.0
                    CU( 3,S ) =  Y
                    CU( 4,S ) =  0.0

                ELSE                              !  c "in" the grid

                    K = ( R - 1 ) * NCOLS  +  C
                    NU( 1,S ) = K
                    NU( 2,S ) = K + 1
                    NU( 3,S ) = K + NCOLS
                    NU( 4,S ) = K + NCOLS + 1
                    P = 1.0 - X
                    Q = 1.0 - Y
                    CU( 1,S ) =  P * Q
                    CU( 2,S ) =  X * Q
                    CU( 3,S ) =  P * Y
                    CU( 4,S ) =  X * Y

                END IF

            END IF      !  end computing bilinear interpolation matrix

11      CONTINUE        !  end matrix computation loop on point sources

        RETURN
        END

