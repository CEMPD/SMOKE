
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  LNK2GRD( XBEG, YBEG, XEND, YEND, 
     &                       NCEL, ACEL, AFAC, ALEN )

C***********************************************************************
C  subroutine body starts at line
C
C  FUNCTION:  Given a link with end points XBEG, YBEG, XEND, YEND:
C	Compute the number NCEL of cells intersected by the link,
C	the cell-numbers (in terms of storage order col + (row-1)*ncols),
C	the total length ALEN of the link, and the lengths
C	AFAC( I ) of cell:link intersections.
C
C  PRECONDITIONS REQUIRED:
C	Grid description set into FDESC3.EXT prior to call
C	right handed coord system (XCELL3D and YCELL3D positive)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C	none
C
C  REVISION  HISTORY:
C	prototype 1/96 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'CHDIMS3.EXT'   !  chemistry dimensioning parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters


C...........   ARGUMENTS and their descriptions:

        REAL        XBEG, YBEG, XEND, YEND	!  end points of link
        INTEGER     NCEL			!  number of intersections
        INTEGER     ACEL( MXCFIP )		!  cell #:  col + (row-1)*ncols
        REAL        AFAC( MXCFIP )		!  frac of link in cell
        REAL        ALEN			!  link length


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        REAL 	XA, YA, XB, YB
        REAL 	XTOP, YTOP
        REAL 	DXDY, DYDX
        REAL 	DDX, DDY, DD
        REAL 	DXLNK, DYLNK
        REAL    XDUM( MXCFIP + 1 )
        REAL    YDUM( MXCFIP + 1 )
        INTEGER XCOL( MXCFIP + 1 )
        INTEGER YROW( MXCFIP + 1 )
        REAL    A, Z, AA, ZZ
        INTEGER JA, JZ, I, J, K, L, M
        INTEGER NX, NY
        INTEGER ROW, COL


C***********************************************************************
C   begin body of subroutine  LNK2GRD

        DDX   = 1.0 / XCELL3D
        DDY   = 1.0 / YCELL3D
        DXLNK = ( XEND - XBEG )
        DYLNK = ( YEND - YBEG )
        ALEN  = SQRT( DXLNK * DXLNK  +  DYLNK * DYLNK )

        IF ( ALEN .LT. 1.0E-4 ) THEN	!  10 cm threshold for length :-)
            NCEL = 0
            RETURN
        END IF
        
C.......   Compute intersections of link with {X = XORIG3D + k * XCELL3D}

        IF ( XBEG .LT. XEND ) THEN
            A = XBEG
            Z = XEND
        ELSE 
            A = XEND
            Z = XBEG
        END IF

        XTOP  = XORIG3D  + XCELL3D * FLOAT( NCOLS3D )

        IF ( A .GT. XTOP  .OR.  Z .LT. XORIG3D ) THEN
            NCEL = 0
            RETURN
        ELSE
            IF ( A .GT. XORIG3D ) THEN
                JA = -INT( DDX * ( XORIG3D - A ) )
                J  = 1
                XDUM( 1 ) = A
                XCOL( 1 ) = JA
            ELSE
                JA  = 0
                J   = 0
            END IF
            IF ( Z .LT. XTOP ) THEN
                JZ = -INT( DDX * ( XORIG3D - Z ) )
                DO  11   I = JA, JZ
                    J = J + 1
                    XDUM( J ) = XORIG3D + XCELL3D * FLOAT( I )
                    XCOL( J ) = I + 1
11              CONTINUE
                J = J + 1
                XDUM( J ) = Z
            ELSE
                JZ = NCOLS3D
                DO  22  I = JA, JZ
                    J = J + 1
                    XDUM( J ) = XORIG3D + XCELL3D * FLOAT( I )
                    XCOL( J ) = I + 1
22              CONTINUE
            END IF
            NX = J

        END IF

C.......   Compute intersections of link with {Y = YORIG3D + k * YCELL3D}

        IF ( YBEG .LT. YEND ) THEN
            AA = YBEG
            ZZ = YEND
        ELSE 
            AA = YEND
            ZZ = YBEG
        END IF

        YTOP  = YORIG3D  + YCELL3D * FLOAT( NCOLS3D )

        IF ( AA .GT. YTOP  .OR.  ZZ .LT. YORIG3D ) THEN
            NCEL = 0
            RETURN
        ELSE
            IF ( AA .GT. YORIG3D ) THEN
                JA = -INT( DDY * ( YORIG3D - AA ) )
                J  = 1
                YDUM( 1 ) = AA
                YROW( 1 ) = JA
            ELSE
                JA  = 0
                J   = 0
            END IF
            IF ( Z .LT. YTOP ) THEN
                JZ = -INT( DDY * ( YORIG3D - ZZ ) )
                DO  33   I = JA, JZ
                    J = J + 1
                    YDUM( J ) = YORIG3D + YCELL3D * FLOAT( I )
                    YROW( J ) = I + 1
33              CONTINUE
                J = J + 1
                YDUM( J ) = ZZ
            ELSE
                JZ = NCOLS3D
                DO  44  I = JA, JZ
                    J = J + 1
                    YDUM( J ) = YORIG3D + YCELL3D * FLOAT( I )
                    YROW( J ) = I + 1
44              CONTINUE
            END IF
            NY = J

        END IF


C...........   Put sentinels on the ends of both lists:

        XDUM( NX + 1 ) = 9.9E36
        YDUM( NY + 1 ) = 9.9E36

C...........   Merge the two intersection-lists to compute fractions:

        XA  = XDUM( 1 )
        YA  = YDUM( 1 )
        ROW = YROW( 1 )
        COL = XCOL( 1 )

        IF ( DXLNK .GT. 0.0 ) THEN

            IF ( DYLNK .GT. 0.0 ) THEN

                DD   = 1.0 / DYLNK
                DXDY = DXLNK * DD
                DYDX = DYLNK / DXLNK
                DD   = ALEN * DD
                K    = 0
                L    = 2
                M    = 2

111             CONTINUE	!  merge-loop

                    K = K + 1
                    ACEL( K ) = COL + NCOLS3D * ( ROW - 1 )

                    XB = A   +  DXDY * ( YDUM( L ) - AA )
                    YB = AA  +  DYDX * ( XDUM( M ) -  A )

                    IF ( XB .LT. XDUM( M ) ) THEN
                        YB = YDUM( L )
                        AFAC( K ) = DD * ( YB - YA )
                        ROW = YROW( L )
                        L   = L + 1
                        GO TO 111
                    ELSE IF ( YB .LT. YDUM( L ) ) THEN
                        XB = XDUM( M )
                        AFAC( K ) = DD * ( YB - YA )
                        COL = XCOL( M )
                        M   = M + 1
                        GO TO 111
                    ELSE
                        XB = XDUM( M )
                        YB = YDUM( L )
                        AFAC( K ) = DD * ( YB - YA )
                        ROW = YROW( L )
                        COL = XCOL( M )
                        L   = L + 1
                        M   = M + 1
                        GO TO 111
                    END IF

            ELSE 			!  dylnk zero:  problem in X only

                ROW = NCOLS3D * ( YROW( 1 ) - 1 )
                DO  122  K = 2, NX
                    ACEL( K ) = XCOL( I-1 ) + ROW
                    AFAC( K ) = XDUM( I ) - A
                    A         = XDUM( I )
122             CONTINUE

            END IF

        ELSE IF ( DYLNK .GT. 0.0 ) THEN		!  problem in Y only

                COL = XCOL( 1 )
                DO  133  K = 2, NX
                    ACEL( K ) = NCOLS3D * ( XCOL( I-1 ) - 1 ) + COL
                    AFAC( K ) = YDUM( I ) - AA
                    AA        = YDUM( I )
133             CONTINUE
        END IF


        RETURN
        END

