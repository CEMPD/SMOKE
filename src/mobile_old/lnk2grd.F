
        SUBROUTINE  LNK2GRD( XBEG, YBEG, XEND, YEND, 
     &                       NCEL, ACEL, AFRAC, ALEN )

C***********************************************************************
C  subroutine body starts at line  94
C
C  FUNCTION:  Given a link with end points XBEG, YBEG, XEND, YEND:
C       Compute the number NCEL of cells intersected by the link,
C       the cell-numbers (in terms of storage order col + (row-1)*ncols),
C       the total length ALEN of the link, and the fraction
C       AFRAC( I ) of each cell:link intersections.
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

C...........   COMMON in FDESC3.EXT provides XCELL3D, YCELL3D, etc.

C...........   ARGUMENTS and their descriptions:

        REAL        XBEG, YBEG, XEND, YEND      !  end points of link
        INTEGER     NCEL                        !  number of intersections
        INTEGER     ACEL ( MXCFIP )             !  cell #:  col + (row-1)*ncols
        REAL        AFRAC( MXCFIP )             !  length of link in cell
        REAL        ALEN                        !  link length


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        REAL     DDX, DDY                !  1/cellsize
        REAL     DXLNK                   !  One over x-dir link length
        REAL     DYLNK                   !  One over y-dir link length
        REAL     END                     !  scratch ending coordinate
        REAL     FF                      !  scratch factor vbles
        REAL     FAC                     !  fraction of link traversed so far
        REAL     RSAV
        REAL     START                   !  start starting coordinate
        REAL     XX, YY
        REAL     XLNK                    !  link length in x direction
        REAL     YLNK                    !  link length in y direction
        REAL     XFAC( 0:MXCFIP + 1 )    !  frac of link traversed
        REAL     YFAC( 0:MXCFIP + 1 )    !  frac of link traversed

        INTEGER  CCC                     !  ending   x-cell of link
        INTEGER  COL                     !  starting x-cell of link
        INTEGER  IX, IY
        INTEGER  ISAV
        INTEGER  J, IR, IC, JZ  
        INTEGER  NX, NY
        INTEGER  ROW                     !  starting y-cell of link
        INTEGER  RRR                     !  ending   y-cell of link
        INTEGER  XINC, YINC              !  merge counters increments
        INTEGER  XCOL( MXCFIP + 1 )      !  subscript for this grid intersection
        INTEGER  YROW( MXCFIP + 1 )      !  subscript for this grid intersection
        INTEGER  XB, XE                  !  pntr to calc fracs @ beg and end of link
        INTEGER  YB, YE                  !  pntr to calc fracs @ beg and end of link

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

        IF( ALEN .EQ. 0.0 ) THEN 

            NCEL = -1
            RETURN

        ELSEIF ( ( CCC .EQ. COL  .AND. RRR .EQ. ROW ) .OR.
     &           ( ALEN .LT. 1.0E-5 * XCELL3D )            ) THEN

C.............  Case: link is outside domain
            IF ( COL .LT. 1  .OR.  COL .GT. NCOLS3D .OR.
     &           ROW .LT. 1  .OR.  ROW .GT. NCOLS3D      ) THEN

                NCEL = 0

C.............  Case: link is completely in a single cell
            ELSE

                NCEL = 1
                ACEL ( 1 ) = COL + NCOLS3D * ( ROW - 1 )
                AFRAC( 1 ) = 1.0

            END IF

            RETURN

        END IF             !  O-D or within-a-cell case

C...........   Calculate once and only once for remaining cases

        DXLNK = 1.0 / ABS( XLNK )  !  safe "/" -- 0-D case already done.
        DYLNK = 1.0 / ABS( YLNK )  !  safe "/" -- 0-D case already done.

C...........   1-D problems:

C.........  Case:  (Near)constant X -- 1-D problem in Y only

        IF ( ( COL .EQ. CCC )  .OR.
     &       ( ABS( XLNK ) .LT. 1.0E-5 * XCELL3D ) ) THEN
 
C.............  Case:  Link is outside domain
            IF ( COL .LT. 1 .OR. COL .GT. NCOLS3D ) THEN
                NCEL = 0
                RETURN
            END IF
 
C.............  For cases in which YBEG is not located in a grid cell lower
C.............  than YEND, then flip the coordinates and use the same algorithm.

            IF ( ROW .GT. RRR ) THEN

                ISAV = RRR
                RRR  = ROW
                ROW  = ISAV
 
                RSAV = YEND
                YEND = YBEG
                YBEG = RSAV

            ENDIF 

C.............  Process intersections moving from bottom to top of domain

C.............  Case: Link outside of domain
            IF ( RRR .LT. 1  .OR.  ROW .GT. NROWS3D ) THEN
                NCEL = 0
                RETURN

C.............  Case: Link crosses domain at some point
            ELSE

C.................  Initialize for all cells (inside domain) intersecting link
                J  = 0
                FF = DYLNK * YCELL3D
                DO  11  IR = ROW, RRR

                    IF( IR .GE. 1 .AND. IR. LE. NROWS3D ) THEN
                        J  = J + 1
                        ACEL ( J ) = COL + NCOLS3D * ( IR - 1 )
                        AFRAC( J ) = FF
                    ENDIF

11              CONTINUE

                NCEL = J

C.................  Reset first fraction if it starts on interior of domain
                IF( YBEG .GT. YORIG3D ) THEN
                    YY = YORIG3D + YCELL3D * FLOAT( ROW )
                    AFRAC( 1 ) = DYLNK * ( YY - YBEG )
                ENDIF

C.................  Reset last  fraction if it ends   on interior of domain
                IF( YEND .LT. YORIG3D + YCELL3D * NROWS3D ) THEN
                    YY = YORIG3D + YCELL3D * FLOAT( RRR - 1 )
                    AFRAC( NCEL ) = DYLNK * ( YEND - YY )
                ENDIF

                RETURN

            END IF      !  if link crosses over domain at all

        END IF          !  if link crosses cells in Y direction only


C.........  Case: (Near)-constant Y -- 1-D problem in X only:

        IF ( ( ROW .EQ. RRR ) .OR. 
     &       ( ABS( YLNK ) .LT. 1.0E-5 * YCELL3D ) ) THEN

C.............  Case:  Link is outside domain
            IF ( ROW .LT. 1  .OR. ROW .GT. NROWS3D ) THEN
                NCEL = 0
                RETURN
            END IF

C.............  For cases in which XBEG is not located in a grid cell to the left
C.............  of XEND, then flip the coordinates and use the same algorithm.

            IF ( COL .GT. CCC ) THEN

                ISAV = CCC
                CCC  = COL
                COL  = ISAV

                RSAV = XEND
                XEND = XBEG
                XBEG = RSAV

            ENDIF

C.................  Case:  Link is outside domain
            IF ( CCC .LT. 1  .OR.  COL .GT. NCOLS3D ) THEN
                NCEL = 0
                RETURN

C.............  Case: Link crosses domain at some point
            ELSE

C.................  Initialize for all cells (inside domain) intersecting link
                J  = 0
                FF = DXLNK * XCELL3D
                JZ = NCOLS3D * ( ROW - 1 )
                DO  22  IC = COL, CCC
                    J  = J + 1
                    ACEL ( J ) = IC + JZ
                    AFRAC( J ) = FF
22              CONTINUE

                NCEL = J

C.................  Reset first fraction if it starts on interior of domain
                IF( XBEG .GT. XORIG3D ) THEN
                    XX = XORIG3D + XCELL3D * FLOAT( COL )
                    AFRAC( 1 ) = DXLNK * ( XX - XBEG )
                ENDIF

C.................  Reset last  fraction if it ends   on interior of domain
                IF( XEND .LT. XORIG3D + XCELL3D * NCOLS3D ) THEN
                    XX = XORIG3D + XCELL3D * FLOAT( CCC - 1 )
                    AFRAC( NCEL ) = DXLNK * ( XEND - XX )
                ENDIF

                RETURN

            END IF      !  if link crosses over domain at all
 
        END IF          !  if link crosses cells in X direction only


C...................................................................
C.........  2-D CASE: 
C.........  (by construction, we know that | XLNK |, | YLNK | > 1.0E5*CELL3D)
C.........  Find the intersections of link with the grid, starting
C.........  from XBEG, YBEG:
C...................................................................

C...................................................................
C.........  Set fractions of link between vertical grid lines by
C.........  using the X fractions 
C...................................................................

C.........  Establish loop bounds and starting and ending coordinates
C.........  Also, precalculate pointers (XB & XE) for adjusting fractions 
C.........  based on orientation of link.  NOTE: XINC is used for 
C.........  stepping through loop _AND_ for changing sign in the XFAC calc
        IX    = COL
        NX    = CCC
        START = XBEG
        END   = XEND

        IF ( COL .LT. CCC ) THEN
            XINC = 1
            XB   = IX
            XE   = NX - 1
        ELSE
            XINC = -1
            XB   = IX - 1
            XE   = NX
        END IF

        J         = 0         !  cell counter
        XFAC( J ) = 9.0E36    !  sentinel on front end-of-list

C.........  Initialize for all cells (inside domain) intersecting link
        FF = DXLNK * XCELL3D
        DO  111  IC = IX, NX, XINC

C.............  First X-cell of link is inside domain
            IF( IC .EQ. IX .AND. START .GE. XORIG3D .AND.
     &          IC .GT. 1  .AND. IC    .LE. NCOLS3D       ) THEN
                J         = J + 1
                XX        = XORIG3D + XCELL3D * FLOAT( XB )
                XCOL( J ) = IC
                XFAC( J ) = DXLNK * FLOAT( XINC ) * ( XX - START )

C.............  Initialize first cell when link starts outside domain
            ELSEIF( IC .EQ. 1 ) THEN
                J         = J + 1
                XCOL( J ) = IC
                XFAC( J ) = FF

C.............  Last X-cell of link is inside domain
            ELSEIF( IC  .EQ. NX .AND. 
     &              END .LE. XORIG3D + XCELL3D * NCOLS3D .AND.
     &              IC  .GT. 1  .AND. IC .LE. NCOLS3D          ) THEN
                J         = J + 1
                XX        = XORIG3D + XCELL3D * FLOAT( XE )
                XCOL( J ) = IC
                XFAC( J ) = XFAC( J-1 ) + 
     &                      DXLNK * FLOAT( XINC ) * ( END - XX )

C.............  Set fractions for interior of domain
            ELSEIF( IC .GT. 1 .AND. IC .LE. NCOLS3D ) THEN
                J         = J + 1
                XCOL( J ) = IC
                XFAC( J ) = XFAC( J-1 ) + FF

            ENDIF

111     CONTINUE

        NX           = J         !  total number of columns intersected
        XFAC( NX+1 ) = 9.0E36    !  sentinel on tail end-of-list

C.........  Case: Link is outside domain
        IF ( NX .EQ. 0 ) THEN
            NCEL = 0
            RETURN
        END IF


C...................................................................
C.........  Set fractions of link between horizontal grid lines by
C.........  using the Y fractions 
C...................................................................

C.........  Establish loop bounds and starting and ending coordinates
C.........  Also, precalculate pointers (YB & YE) for adjusting fractions 
C.........  based on orientation of link.  NOTE: YINC is used for 
C.........  stepping through loop _AND_ for changing sign in the YFAC calc
        IY    = ROW
        NY    = RRR
        START = YBEG
        END   = YEND

        IF ( ROW .LT. RRR ) THEN
            YINC = 1
            YB   = IY
            YE   = NY - 1
        ELSE
            YINC = -1
            YB   = IY - 1
            YE   = NY
        END IF

        J         = 0         !  cell counter
        YFAC( J ) = 9.0E36    !  sentinel on front end-of-list

C.........  Calculate fractions for cells (inside domain) intersecting link
        FF = DYLNK * YCELL3D
        DO  222  IR = IY, NY, YINC

C.............  First Y-cell of link is inside domain
            IF( IR .EQ. IY .AND. START .GE. YORIG3D .AND.
     &          IR .GT. 1  .AND. IR    .LE. NROWS3D      ) THEN
                J         = J + 1
                YY        = YORIG3D + YCELL3D * FLOAT( YB )
                YROW( J ) = IR
                YFAC( J ) = DYLNK * FLOAT( YINC ) * ( YY - START )

C.............  Initialize first cell when link starts outside domain
            ELSEIF( IR .EQ. 1 ) THEN
                J         = J + 1
                YROW( J ) = IR
                YFAC( J ) = FF

C.............  Last Y-cell of link is inside domain
            ELSEIF( IR  .EQ. NY .AND. 
     &              END .LE. YORIG3D + YCELL3D * NROWS3D .AND.
     &              IR  .GT. 1 .AND. IR .LE. NROWS3D           ) THEN
                J         = J + 1
                YY        = YORIG3D + YCELL3D * FLOAT( YE )
                YROW( J ) = IR
                YFAC( J ) = YFAC( J-1 ) + 
     &                      DYLNK * FLOAT( YINC ) * ( END - YY )

C.............  Set fractions for interior of domain
            ELSEIF( IR .GT. 1 .AND. IR .LE. NROWS3D ) THEN
                J         = J + 1
                YROW( J ) = IR
                YFAC( J ) = YFAC( J-1 ) + FF

            ENDIF

222     CONTINUE

        NY           = J         !  total number of columns intersected
        YFAC( NY+1 ) = 9.0E36    !  sentinel on tail end-of-list

C.........  Case: Link is outside domain
        IF ( NY .EQ. 0 ) THEN
            NCEL = 0
            RETURN
        END IF


C...................................................................
C........   Now merge the two intersection lists:
C........   This algorithm starts at one end of the link,
C........   and continues to the other.  The XFAC and YFAC variables
C........   contain the total fraction of the link up to the point 
C........   where the link intersects column IX (for XFAC) and/or 
C........   row IY (for YFAC).   The next intersection produces the
C........   next AFRAC for cell (IC,IR).  If the link passes through
C........   a cell corner, then both the IX counter and IY counter
C........   are incremented by 1.
C...................................................................

C.........  Perform the merge.  Loop terminates when both lists hit
C........  sentinels on either end of fractions arrays

        J  = 0
        FF = 0
        IX = 1
        IY = 1
        XINC = 1
        YINC = 1
133     CONTINUE

C.............  Intersect new column next           
            IF ( XFAC( IX ) .LT. YFAC( IY ) ) THEN
                FAC = XFAC( IX ) - FF
                FF  = XFAC( IX )
                IC  = XCOL( IX )
                IR  = YROW( IY )
                IX  = IX + XINC

C.............  Intersect new column and new row at same time
C.............  or both arrays have reached end of lists
            ELSEIF( XFAC( IX ) .EQ. YFAC( IY ) ) THEN
                FAC = XFAC( IX ) - FF
                FF  = XFAC( IX )
                IC  = XCOL( IX )
                IR  = YROW( IY )
                IX  = IX + XINC
                IY  = IY + YINC

C.............  Intersect new row next
            ELSE
                FAC = YFAC( IY ) - FF
                FF  = YFAC( IY )
                IC  = XCOL( IX )
                IR  = YROW( IY )
                IY  = IY + YINC
            END IF
            
            IF ( FF .LE. 2.0 ) THEN   ! Check for sentinel
                J = J + 1
                ACEL ( J ) = IC +  NCOLS3D * ( IR - 1 )
                AFRAC( J ) = MIN( 1.0, MAX( 0.0, FAC ) )
                GO TO  133              !  to head of loop
            END IF


C...........   Merge complete.  Return the number of cells found:

        NCEL = J
        RETURN
        END

