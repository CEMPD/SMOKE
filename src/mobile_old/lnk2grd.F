
        SUBROUTINE  LNK2GRD( XBEG, YBEG, XEND, YEND, 
     &                       NCEL, ACEL, AFAC, ALEN )

C***********************************************************************
C  subroutine body starts at line  94
C
C  FUNCTION:  Given a link with end points XBEG, YBEG, XEND, YEND:
C       Compute the number NCEL of cells intersected by the link,
C       the cell-numbers (in terms of storage order col + (row-1)*ncols),
C       the total length ALEN of the link, and the lengths
C       AFAC( I ) of cell:link intersections.
C
C  PRECONDITIONS REQUIRED:
C       Grid description set into FDESC3.EXT prior to call
C       right handed coord system (XCELL3D and YCELL3D positive)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       none
C
C  REVISION  HISTORY:
C       prototype 1/96 by CJC
C       New algorithm 5/96 by CJC  (removes termination check bug)
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'MBDIMS3.EXT'   ! mobile-source dimensioning parameters
        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations


C...........   ARGUMENTS and their descriptions:

        REAL        XBEG, YBEG, XEND, YEND      !  end points of link
        INTEGER     NCEL                        !  number of intersections
        INTEGER     ACEL( MXCFIP )              !  cell #:  col + (row-1)*ncols
        REAL        AFAC( MXCFIP )              !  frac of link in cell
        REAL        ALEN                        !  link length


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        REAL     XA, YA, XX, YY
        REAL     DDX, DDY                !  1/cellsize
        REAL     FF, GG                  !  scratch factor vbles
        REAL     FA, FX, FY, GX, GY      !  scratch factor vbles
        REAL     XLNK, YLNK              !  link projection in x(y) coord.
        REAL     DXLNK, DYLNK            !  1/link projection in x(y) coord.
        REAL     XFAC( 0:MXCFIP + 1 )    !  frac of link traversed
        REAL     YFAC( 0:MXCFIP + 1 )    !  frac of link traversed
        INTEGER  XROW( MXCFIP + 1 )      !  subscript for this grid intersection
        INTEGER  XCOL( MXCFIP + 1 )      !  subscript for this grid intersection
        INTEGER  YROW( MXCFIP + 1 )      !  subscript for this grid intersection
        INTEGER  YCOL( MXCFIP + 1 )      !  subscript for this grid intersection
        REAL     FAC                     !  fraction of link traversed so far
        INTEGER  J, IR, IC, JZ, CCC, RRR
        INTEGER  NX, NY
        INTEGER  IX, IY
        INTEGER  XINC, YINC
        INTEGER  ROW, COL


C***********************************************************************
C   begin body of subroutine  LNK2GRD

C...........   Initializations

        XLNK = XEND - XBEG
        YLNK = YEND - YBEG
                
        DDX   = 1.0 / XCELL3D
        DDY   = 1.0 / YCELL3D

        COL  = 1 + INT( DDX * ( XBEG - XORIG3D ) )	!  starting cell
        ROW  = 1 + INT( DDY * ( YBEG - YORIG3D ) )

        CCC  = 1 + INT( DDX * ( XEND - XORIG3D ) )	!  ending cell
        RRR  = 1 + INT( DDY * ( YEND - YORIG3D ) )

        ALEN = SQRT( XLNK**2 + YLNK**2 )		!  link length


C...............   Check for 0-D, 1-D cases, using tolerance of
C...............   1.0E-5 * cellsize (on the order of 0.5-10.0 meters).
C...........  0-D case, within-a-cell case:

        IF ( ( CCC .EQ. COL  .AND. RRR .EQ. ROW ) .OR.
     &       ( ALEN .LT. 1.0E-5 * XCELL3D )         ) THEN

            IF ( COL .LT. 1  .OR.  COL .GT. NCOLS3D .OR.
     &           ROW .LT. 1  .OR.  ROW .GT. NCOLS3D ) THEN

                NCEL = 0

            ELSE

                NCEL = 1
                ACEL( 1 ) = COL + NCOLS3D * ( ROW - 1 )
                AFAC( 1 ) = ALEN

            END IF

            RETURN

        END IF             !  O-D or within-a-cell case


C...........   1-D problems:
C...............   (Near)constant X -- 1-D problem in Y only

        IF ( ( COL .EQ. CCC )  .OR.
     &       ( ABS( XLNK ) .LT. 1.0E-5 * XCELL3D ) ) THEN
 
            IF ( COL .LT. 1 .OR. COL .GT. NCOLS3D ) THEN
                NCEL = 0
                RETURN
            END IF
 
            DYLNK = ALEN / ABS( YLNK )  !  safe "/" -- 0-D case already done.
 
            IF ( ROW .LT. RRR ) THEN             ! initialization, termination

                IF ( RRR .LT. 1  .OR.  ROW .GT. NROWS3D) THEN
                    NCEL = 0
                    RETURN
                ELSE IF ( ROW .LT. 1 ) THEN
                    IY = 1
                    IF ( RRR .LE. NROWS3D ) THEN
                        J  = 1
                        NY = RRR - 1
                        YY = YORIG3D + YCELL3D * FLOAT( RRR - 1 )
                        ACEL( 1 ) = COL + NCOLS3D * ( RRR - 1 )
                        AFAC( 1 ) = DYLNK * ( YEND - YY )
                    ELSE
                        J  = 0
                        NY = NCOLS3D
                    END IF
                ELSE            !  row >= 1 
                    IY = ROW + 1
                    YA = YORIG3D + YCELL3D * FLOAT( ROW )
                    ACEL( 1 ) = COL + NCOLS3D * ( ROW - 1 )
                    AFAC( 1 ) = DYLNK * ( YA - YBEG )
                    IF ( RRR .LE. NROWS3D ) THEN
                        J  = 2
                        NY = RRR  - 1
                        YY = YORIG3D + YCELL3D * FLOAT( RRR - 1 )
                        ACEL( 2 ) = COL + NCOLS3D * ( RRR - 1 )
                        AFAC( 2 ) = DYLNK * ( YEND - YY )
                    ELSE
                        J  = 1
                        NY = NROWS3D
                    END IF
                END IF

            ELSE                        !  row > rr, by construction

                IF ( ROW .LT. 1  .OR.  RRR .GT. NROWS3D ) THEN
                    NCEL = 0
                    RETURN
                ELSE IF ( RRR .LT. 1 ) THEN
                    IY = 1
                    IF ( ROW .LE. NROWS3D ) THEN
                        J  = 1
                        NY = ROW - 1
                        YY = YORIG3D + YCELL3D * FLOAT( ROW - 1 )
                        ACEL( 1 ) = COL + NCOLS3D * ( ROW - 1 )
                        AFAC( 1 ) = DYLNK * ( YEND - YY )
                    ELSE
                        J  = 0
                        NY = NCOLS3D
                    END IF
                ELSE                    !  row >= 1 
                    IY = RRR + 1
                    YA = YORIG3D + YCELL3D * FLOAT( RRR )
                    ACEL( 1 ) = COL + NCOLS3D * ( RRR - 1 )
                    AFAC( 1 ) = DYLNK * ( YA - YBEG )
                    IF ( ROW .LE. NROWS3D ) THEN
                        J  = 2
                        NY = ROW  - 1
                        YY = YORIG3D + YCELL3D * FLOAT( ROW - 1 )
                        ACEL( 2 ) = COL + NCOLS3D * ( ROW - 1 )
                        AFAC( 2 ) = DYLNK * ( YEND - YY )
                    ELSE
                        J  = 1
                        NY = NROWS3D
                    END IF
                END IF

            END IF                      ! row <,> rrr (unequal, by construction)
 
            FF = DYLNK * YCELL3D
            DO  11  IR = IY, NY
                J  = J + 1
                ACEL( J ) = COL + NCOLS3D * ( IR - 1 )
                AFAC( J ) = FF
11          CONTINUE

            NCEL = J
            RETURN

        END IF          !  if problem in Y only


C...............   (Near)-constant Y -- 1-D problem in X only:

        IF ( ( ROW .EQ. RRR ) .OR. 
     &       ( ABS( YLNK ) .LT. 1.0E-5 * YCELL3D ) ) THEN

            IF ( ROW .LT. 1  .OR. ROW .GT. NROWS3D ) THEN
                NCEL = 0
                RETURN
            END IF

            DXLNK = ALEN / ABS( XLNK )  !  safe "/" -- 0-D case already done.
 
            IF ( COL .LT. CCC ) THEN             !  initialization, termination

                IF ( CCC .LT. 1  .OR.  COL .GT. NCOLS3D ) THEN
                    NCEL = 0
                    RETURN
                ELSE IF ( COL .LT. 1 ) THEN
                    IX = 1
                    IF ( CCC .LE. NCOLS3D ) THEN
                        J  = 1
                        NX = CCC - 1
                        XX = XORIG3D + XCELL3D * FLOAT( CCC - 1 )
                        ACEL( 1 ) = CCC + NCOLS3D * ( ROW - 1 )
                        AFAC( 1 ) = DXLNK * ( XEND - XX )
                    ELSE
                        J  = 0
                        NX = NCOLS3D
                    END IF
                ELSE            !  col >= 1 
                    IX = COL + 1
                    XA = XORIG3D + XCELL3D * FLOAT( COL )
                    ACEL( 1 ) = COL + NCOLS3D * ( ROW - 1 )
                    AFAC( 1 ) = DXLNK * ( XA - XBEG )
                    IF ( CCC .LE. NCOLS3D ) THEN
                        J  = 2
                        NX = CCC  - 1
                        XX = XORIG3D + XCELL3D * FLOAT( CCC - 1 )
                        ACEL( 2 ) = CCC + NCOLS3D * ( ROW - 1 )
                        AFAC( 2 ) = DXLNK * ( XEND - XX )
                    ELSE
                        J  = 1
                        NX = NCOLS3D
                    END IF
                END IF

            ELSE                        !  col > ccc, by construction

                IF ( COL .LT. 1  .OR.  CCC .GT. NCOLS3D ) THEN
                    NCEL = 0
                    RETURN
                ELSE IF ( CCC .LT. 1 ) THEN
                    IX = 1
                    IF ( COL .LE. NCOLS3D ) THEN
                        J  = 1
                        NX = COL - 1
                        XX = XORIG3D + XCELL3D * FLOAT( COL - 1 )
                        ACEL( 1 ) = COL + NCOLS3D * ( ROW - 1 )
                        AFAC( 1 ) = DXLNK * ( XEND - XX )
                    ELSE
                        J  = 0
                        NX = NCOLS3D
                    END IF
                ELSE                    !  ccc >= 1 
                    IX = CCC + 1
                    XA = XORIG3D + XCELL3D * FLOAT( CCC )
                    ACEL( 1 ) = CCC + NCOLS3D * ( ROW - 1 )
                    AFAC( 1 ) = DXLNK * ( XA - XBEG )
                    IF ( COL .LE. NCOLS3D ) THEN
                        J  = 2
                        NX = COL  - 1
                        XX = XORIG3D + XCELL3D * FLOAT( COL - 1 )
                        ACEL( 2 ) = COL + NCOLS3D * ( ROW - 1 )
                        AFAC( 2 ) = DXLNK * ( XEND - XX )
                    ELSE
                        J  = 1
                        NX = NCOLS3D
                    END IF
                END IF

            END IF                      ! col <,> ccc (unequal, by construction)
 
            FF = DXLNK * XCELL3D
            JZ = NCOLS3D * ( ROW - 1 )
            DO  22  IC = IX, NX
                J  = J + 1
                ACEL( J ) = IC + JZ
                AFAC( J ) = FF
22          CONTINUE

            NCEL = J
            RETURN

        END IF          !  if problem in X only


C...............   Else 2-D case: 
C.......   (by construction, we know that | XLNK |, | YLNK | > 1.0E5*CELL3D)
C.......   Find the intersections of link with the grid, starting
C.......   from XBEG, YBEG:

        IF ( COL .LT. CCC ) THEN         !  establish loop bounds
            IX = MAX( COL, 1 )
            NX = MIN( CCC, NCOLS3D )
        ELSE
            IX = MAX( CCC, 1 )
            NX = MIN( COL, NCOLS3D )
        END IF

        YFAC( 0 ) = 9.0E36      !  sentinel on end-of-list
        XFAC( 0 ) = 9.0E36      !  sentinel on end-of-list
        DXLNK = 1.0 / XLNK
        FA    = DXLNK * ( XORIG3D - XBEG )      !  starting fraction
        FX    = XCELL3D * DXLNK
        FY    = DDY * ( YBEG - YORIG3D )
        GY    = DDY * YLNK
        J     = 0
        DO  111  IC = IX, NX
            GG = FA + FX * FLOAT( IC )  !  link fraction from XBEG thru col IC
            IF ( GG .GE. 0.0  .AND. GG .LE. 1.0 ) THEN
                RRR = 1 + INT( FY + GY * GG )
                IF ( RRR .GE. 1  .AND. RRR .LE. NROWS3D ) THEN
                    J = J + 1
                    XFAC( J ) = GG
                    XCOL( J ) = IC
                    XROW( J ) = RRR
                END IF
            END IF
111     CONTINUE

        IF ( J .EQ. 0 ) THEN
            NCEL = 0
            RETURN
        END IF

        NX = J
        XFAC( J+1 ) = 9.0E36    !  sentinel on end-of-list


C.......   Find the intersections of link with the grid, starting
C.......   from XBEG, YBEG:

        IF ( ROW .LT. RRR ) THEN         !  establish loop bounds
            IY = MAX( ROW, 1 )
            NY = MIN( RRR, NROWS3D )
        ELSE
            IY = MAX( RRR, 1 )
            NY = MIN( ROW, NROWS3D )
        END IF

        DYLNK = 1.0 / YLNK
        FA    = DYLNK * ( YORIG3D - YBEG )      !  starting fraction
        FY    = YCELL3D * DYLNK
        FX    = DDX * ( XBEG - XORIG3D )
        GX    = DDX * XLNK
        J     = 0
        DO  122  IR = IY, NY
            GG = FA + FY * FLOAT( IR )  !  link fraction from XBEG thru row IR
            IF ( GG .GE. 0.0  .AND. GG .LE. 1.0 ) THEN
                CCC = 1 + INT( FX + GX * GG )
                IF ( CCC .GE. 1  .AND. CCC .LE. NCOLS3D ) THEN
                    J = J + 1
                    YFAC( J ) = GG
                    YROW( J ) = IR
                    YCOL( J ) = CCC
                END IF
            END IF
122     CONTINUE

        IF ( J .EQ. 0 ) THEN
            NCEL = 0
            RETURN
        END IF

        NY = J
        YFAC( J+1 ) = 9.0E36    !  sentinel on end-of-list


C...........   Now merge the two intersection lists:
C.......   Get the starting point:

        IF ( XBEG .LT. XEND ) THEN
            XINC = 1
            IX   = 1
        ELSE
            XINC = -1
            IX   = NX
        END IF

        IF ( YBEG .LT. YEND ) THEN
            YINC = 1
            IY   = 1
        ELSE
            YINC = -1
            IY   = NY
        END IF

        IF ( COL .LT. 1  .OR.  COL .GT. NCOLS3D .OR.
     &       ROW .LT. 1  .OR.  ROW .GT. NCOLS3D ) THEN

            FF = 0.0    !  start inside the grid

        ELSE IF ( XFAC( IX ) .GT. YFAC( IY ) ) THEN

            FF  = MIN( 1.0, MAX( 0.0, XFAC( IX ) - FX ) )
            COL = XCOL( IX )
            ROW = XROW( IX )

        ELSE    !  else xfac(1) <= yfac(1 )

            FF  = MIN( 1.0, MAX( 0.0, YFAC( IY ) - GY ) )
            COL = YCOL( IY )
            ROW = YROW( IY )

        END IF

        J = 0
133     CONTINUE        !  merge-loop (terminates when both lists hit sentinels)
            
            IF ( XFAC( IX ) .LT. YFAC( IY ) ) THEN
                FAC = XFAC( IX ) - FF
                FF  = XFAC( IX )
                COL = IX
                IX  = IX + XINC
            ELSE
                FAC = YFAC( IY ) - FF
                FF  = YFAC( IY )
                ROW = IY
                IY  = IY + YINC
            END IF
            
            IF ( FF .LE. 2.0 ) THEN
                J = J + 1
                AFAC( J ) = ALEN * MIN( 1.0, MAX( 0.0, FAC ) )
                ACEL( J ) = COL  +  NCOLS3D * ( ROW - 1 )
                GO TO  133              !  to head of loop
            END IF


C...........   Merge complete.  Return the number of cells found:

        NCEL = J
        RETURN
        END

