
        SUBROUTINE LNK2GRD( NDIM, XBEGIN, YBEGIN, XENDIN, YENDIN, 
     &                      NCEL, ACEL, AFRAC, ALEN, EFLAG )

C***********************************************************************
C  subroutine body starts at line 111
C
C  FUNCTION:  Given a link with end points XBEG, YBEG, XEND, YEND:
C       Compute the number NCEL of cells intersected by the link,
C       the cell-numbers (in terms of storage order col + (row-1)*ncols),
C       the total length ALEN of the link, and the fraction
C       AFRAC( I ) of each cell:link intersections.
C
C  PRECONDITIONS REQUIRED:
C       Grid description stored in MODGRID prior to call
C       right handed coord system (XCELL and YCELL positive)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Copied from lnk2grd.F 3.2 by mhouyoux
C       Updated for new SMOKE in 5/99
C
C***********************************************************************
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

C...........   MODULES for public variables
C.........  This module contains the global variables for the 3-d grid
        USE M3UTILIO

        USE MODGRID, ONLY: XCELL, YCELL, XORIG, YORIG, NCOLS, NROWS

        IMPLICIT NONE

C...........   INCLUDES:

C        INCLUDE 'PARMS3.EXT'    ! I/O API constants
C        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
C        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations

C...........   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN) :: NDIM          !  max number of intersections
        REAL   , INTENT (IN) :: XBEGIN        !  end points of link
        REAL   , INTENT (IN) :: YBEGIN        !  end points of link
        REAL   , INTENT (IN) :: XENDIN        !  end points of link
        REAL   , INTENT (IN) :: YENDIN        !  end points of link
        INTEGER, INTENT(OUT) :: NCEL          !  number of intersections
        INTEGER, INTENT(OUT) :: ACEL ( NDIM ) !  cell #:  col + (row-1)*ncols
        REAL   , INTENT(OUT) :: AFRAC( NDIM ) !  length of link in cell
        REAL   , INTENT(OUT) :: ALEN          !  link length
        LOGICAL, INTENT(OUT) :: EFLAG         !  error flag

C.........  Local arrays dimensioned with subroutine arguments
        INTEGER, ALLOCATABLE, SAVE :: XCOL( : )    !  subscript for this grid intersection
        INTEGER, ALLOCATABLE, SAVE :: YROW( : )    !  subscript for this grid intersection
        REAL, ALLOCATABLE, SAVE    :: XFAC( : )    !  frac of link traversed
        REAL, ALLOCATABLE, SAVE    :: YFAC( : )    !  frac of link traversed

C...........   LOCAL VARIABLES and their descriptions:

        REAL     DDX, DDY                !  1/cellsize
        REAL     DXLNK                   !  One over x-dir link length
        REAL     DYLNK                   !  One over y-dir link length
        REAL     ENDS                    !  scratch ending coordinate
        REAL     FF                      !  scratch factor vbles
        REAL     FAC                     !  fraction of link traversed so far
        REAL     RSAV
        REAL     START                   !  start starting coordinate
        REAL     XX, YY
        REAL     XBEG, YBEG              !  local beginning coords
        REAL     XEND, YEND              !  local ending coords
        REAL     XLNK                    !  link length in x direction
        REAL     YLNK                    !  link length in y direction

        INTEGER  CCC                     !  ending   x-cell of link
        INTEGER  COL                     !  starting x-cell of link
        INTEGER  IX, IY
        INTEGER  IOS                     !  i/o status
        INTEGER  ISAV
        INTEGER  J, IR, IC, JZ  
        INTEGER  NX, NY
        INTEGER  ROW                     !  starting y-cell of link
        INTEGER  RRR                     !  ending   y-cell of link
        INTEGER  XINC, YINC              !  merge counters increments
        INTEGER  XB, XE                  !  pntr to calc fracs @ link beg & end
        INTEGER  YB, YE                  !  pntr to calc fracs @ link beg & end

        LOGICAL :: FIRSTIME = .TRUE.     !  true: first time routine called

        CHARACTER(16) :: PROGNAME = 'LNK2GRD'  ! progname name

C***********************************************************************
C   begin body of subroutine  LNK2GRD

C...........   Allocate memory first time routine is called
        IF( FIRSTIME ) THEN
            ALLOCATE( XCOL( NDIM+1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XCOL', PROGNAME )
            ALLOCATE( YROW( NDIM+1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YROW', PROGNAME )
            ALLOCATE( XFAC( 0:NDIM+1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XFAC', PROGNAME )
            ALLOCATE( YFAC( 0:NDIM+1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YFAC', PROGNAME )
            FIRSTIME = .FALSE.
        END IF

C...........   Initializations
        XBEG = XBEGIN
        XEND = XENDIN
        YBEG = YBEGIN
        YEND = YENDIN

        XCOL = 0   ! array
        YROW = 0   ! array
        XFAC = 0.  ! array
        YFAC = 0.  ! array

        XLNK = XEND - XBEG
        YLNK = YEND - YBEG
                
        DDX   = 1.0 / XCELL
        DDY   = 1.0 / YCELL

        COL  = 1 + INT( DDX * ( XBEG - XORIG ) )      !  starting cell
        ROW  = 1 + INT( DDY * ( YBEG - YORIG ) )

        CCC  = 1 + INT( DDX * ( XEND - XORIG ) )      !  ending cell
        RRR  = 1 + INT( DDY * ( YEND - YORIG ) )

        ALEN = SQRT( XLNK**2 + YLNK**2 )                !  link length

C...........  Initialize the number of cells to zero
        NCEL = 0
        
C...............   Check for 0-D, 1-D cases, using tolerance of
C...............   1.0E-5 * cellsize (on the order of 0.5-10.0 meters).
C...............   Within a cell case:
        IF ( ( CCC .EQ. COL  .AND. RRR .EQ. ROW ) .OR.
     &           ( ALEN .LT. 1.0E-5 * XCELL )            ) THEN

C.............  Make sure 1-cell link is inside domain
            IF ( COL .GE. 1  .AND.  COL .LE. NCOLS .AND.
     &           ROW .GE. 1  .AND.  ROW .LE. NROWS      ) THEN

                NCEL = 1
                ACEL ( 1 ) = COL + NCOLS * ( ROW - 1 )
                AFRAC( 1 ) = 1.0

            END IF

            RETURN

        END IF             !  O-D or within-a-cell case

C...........   Calculate once and only once for remaining cases
        IF( XLNK == 0.0 ) XLNK = 1.0
        IF( YLNK == 0.0 ) YLNK = 1.0
        DXLNK = 1.0 / ABS( XLNK )  !  safe "/" -- 0-D case already done.
        DYLNK = 1.0 / ABS( YLNK )  !  safe "/" -- 0-D case already done.

C...........   1-D problems:

C.........  Case:  (Near)constant X -- 1-D problem in Y only

        IF ( ( COL .EQ. CCC )  .OR.
     &       ( ABS( XLNK ) .LT. 1.0E-5 * XCELL ) ) THEN
 
C.............  Case:  Link is outside domain
            IF ( COL .LT. 1 .OR. COL .GT. NCOLS ) THEN
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
            IF ( RRR .LT. 1  .OR.  ROW .GT. NROWS ) THEN
                NCEL = 0
                RETURN

C.............  Case: Link crosses domain at some point
            ELSE

C.................  Initialize for all cells (inside domain) intersecting link
                J  = 0
                FF = DYLNK * YCELL
                DO  11  IR = MAX( 1,ROW ) , MIN( RRR,NROWS )

                    IF( IR .GE. 1 .AND. IR. LE. NROWS ) THEN
                        J  = J + 1
                        ACEL ( J ) = COL + NCOLS * ( IR - 1 )
                        AFRAC( J ) = FF
                    ENDIF

11              CONTINUE

                NCEL = J

C.................  Reset first fraction if it starts on interior of domain
                IF( YBEG .GT. YORIG ) THEN
                    YY = YORIG + YCELL * FLOAT( ROW )
                    AFRAC( 1 ) = DYLNK * ( YY - YBEG )
                    IF( AFRAC( 1 ) < 0.0 ) AFRAC( 1 ) = 0.0
                ENDIF

                IF( NCEL .GT. NDIM ) CALL REPORT_BAD_CELL

C.................  Reset last  fraction if it ends   on interior of domain
                IF( YEND .LT. YORIG + YCELL * NROWS ) THEN
                    YY = YORIG + YCELL * FLOAT( RRR - 1 )
                    AFRAC( NCEL ) = DYLNK * ( YEND - YY )
                    IF( AFRAC( NCEL ) < 0.0 ) AFRAC( NCEL ) = 0.0
                ENDIF

                RETURN

            END IF      !  if link crosses over domain at all

        END IF          !  if link crosses cells in Y direction only


C.........  Case: (Near)-constant Y -- 1-D problem in X only:

        IF ( ( ROW .EQ. RRR ) .OR. 
     &       ( ABS( YLNK ) .LT. 1.0E-5 * YCELL ) ) THEN

C.............  Case:  Link is outside domain
            IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN
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
            IF ( CCC .LT. 1  .OR.  COL .GT. NCOLS ) THEN
                NCEL = 0
                RETURN

C.............  Case: Link crosses domain at some point
            ELSE

C.................  Initialize for all cells (inside domain) intersecting link
                J  = 0
                FF = DXLNK * XCELL
                JZ = NCOLS * ( ROW - 1 )
                DO  22  IC = MAX( 1,COL ) , MIN( CCC,NCOLS )
                    J  = J + 1
                    ACEL ( J ) = IC + JZ
                    AFRAC( J ) = FF
22              CONTINUE
                NCEL = J

                IF( NCEL .GT. NDIM ) CALL REPORT_BAD_CELL

C.................  Reset first fraction if it starts on interior of domain
                IF( XBEG .GT. XORIG ) THEN
                    XX = XORIG + XCELL * FLOAT( COL )
                    AFRAC( 1 ) = DXLNK * ( XX - XBEG )
                    IF( AFRAC( 1 ) < 0.0 ) AFRAC( 1 ) = 0.0
                ENDIF

C.................  Reset last  fraction if it ends   on interior of domain
                IF( XEND .LT. XORIG + XCELL * NCOLS ) THEN
                    XX = XORIG + XCELL * FLOAT( CCC - 1 )
                    AFRAC( NCEL ) = DXLNK * ( XEND - XX )
                    IF( AFRAC( NCEL ) < 0.0 ) AFRAC( NCEL ) = 0.0
                ENDIF

                RETURN

            END IF      !  if link crosses over domain at all
 
        END IF          !  if link crosses cells in X direction only

C...................................................................
C.........  2-D CASE: 
C.........  (by construction, we know that | XLNK |, | YLNK | > 1.0E5*CELL)
C.........  Find the intersections of link with the grid, starting
C.........  from XBEG, YBEG:
C...................................................................

C...................................................................
C.........  Set fractions of link between vertical grid lines by
C.........  using the X fractions 
C...................................................................

C.........  Establish loop bounds and starting and ending coordinates
C.........  Also, precalculate pointers (XB & XE) for adjusting fractions 
C.........  based on orientation of link. 
        IX    = COL
        NX    = CCC
        START = XBEG
        ENDS  = XEND

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
        FF = DXLNK * XCELL

        DO  IC = IX, NX, XINC

C.............  Skip cell ID is outside domain
            IF( IC < 1 ) CYCLE

C.............  First X-cell of link is inside domain
            IF( IC .EQ. IX .AND. START .GE. XORIG .AND.
     &          IC .GE. 1  .AND. IC    .LE. NCOLS       ) THEN
                J         = J + 1
                XX        = XORIG + XCELL * FLOAT( XB )
                XCOL( J ) = IC
                XFAC( J ) = DXLNK * FLOAT( XINC ) * ( XX - START )

C.............  Initialize first cell when link starts outside domain
            ELSEIF( IC .EQ. 1 ) THEN
                J         = J + 1
                XCOL( J ) = IC
                XFAC( J ) = FF

C.............  Last X-cell of link is inside domain
            ELSEIF( IC   .EQ. NX .AND. 
     &              ENDS .LE. XORIG + XCELL * NCOLS .AND.
     &              IC   .GE. 1  .AND. IC .LE. NCOLS          ) THEN
                J         = J + 1
                XX        = XORIG + XCELL * FLOAT( XE )
                XCOL( J ) = IC
                XFAC( J ) = XFAC( J-1 ) + 
     &                      DXLNK * FLOAT( XINC ) * ( ENDS - XX )

C.............  Set fractions for interior of domain and non-end of link
            ELSEIF( IC .GT. 1 .AND. IC .LE. NCOLS ) THEN
                J         = J + 1
                XCOL( J ) = IC
                XFAC( J ) = XFAC( J-1 ) + FF

            ENDIF

            IF( XFAC( J ) < 0.0 ) XFAC( J ) = 0.0

        END DO

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
C.........  based on orientation of link.  Note that YINC is used for 
C.........  stepping through loop _AND_ for changing sign in the YFAC calc
        IY    = ROW
        NY    = RRR
        START = YBEG
        ENDS  = YEND

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
        FF = DYLNK * YCELL

        DO  IR = IY, NY, YINC

C.............  Skip cell ID is outside domain
            IF( IR < 1 ) CYCLE

C.............  First Y-cell of link is inside domain
            IF( IR .EQ. IY .AND. START .GE. YORIG .AND.
     &          IR .GE. 1  .AND. IR    .LE. NROWS      ) THEN
                J         = J + 1
                YY        = YORIG + YCELL * FLOAT( YB )
                YROW( J ) = IR
                YFAC( J ) = DYLNK * FLOAT( YINC ) * ( YY - START )

C.............  Initialize first cell when link starts outside domain
            ELSEIF( IR .EQ. 1 ) THEN
                J         = J + 1
                YROW( J ) = IR
                YFAC( J ) = FF

C.............  Last Y-cell of link is inside domain
            ELSEIF( IR   .EQ. NY .AND. 
     &              ENDS .LE. YORIG + YCELL * NROWS .AND.
     &              IR   .GE. 1 .AND. IR .LE. NROWS           ) THEN
                J         = J + 1
                YY        = YORIG + YCELL * FLOAT( YE )
                YROW( J ) = IR
                YFAC( J ) = YFAC( J-1 ) + 
     &                      DYLNK * FLOAT( YINC ) * ( ENDS - YY )

C.............  Set fractions for interior of domain and non-end of link
            ELSEIF( IR .GT. 1 .AND. IR .LE. NROWS ) THEN
                J         = J + 1
                YROW( J ) = IR
                YFAC( J ) = YFAC( J-1 ) + FF

            ENDIF

            IF( YFAC( J ) < 0.0 ) YFAC( J ) = 0.0

        END DO

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

C.............  Intersect new column and new row at same time
C.............  or both arrays have reached end of lists
            IF( ABS( XFAC( IX ) - YFAC( IY ) ) .LT. 1.0E-05 ) THEN
                FAC = XFAC( IX ) - FF
                FF  = XFAC( IX )
                IC  = XCOL( IX )
                IR  = YROW( IY )
                IX  = IX + XINC
                IY  = IY + YINC

C.............  Intersect new column next           
            ELSE IF ( XFAC( IX ) .LT. YFAC( IY ) ) THEN
                FAC = XFAC( IX ) - FF
                FF  = XFAC( IX )
                IC  = XCOL( IX )
                IR  = YROW( IY )
                IX  = IX + XINC

C.............  Intersect new row next
            ELSE
                FAC = YFAC( IY ) - FF
                FF  = YFAC( IY )
                IC  = XCOL( IX )
                IR  = YROW( IY )
                IY  = IY + YINC
            END IF
            IF( IR > 0 .AND. IC > 0 ) THEN
                IF ( FF .LE. 2.0 ) THEN   ! Check for sentinel
                    J = J + 1
                    ACEL ( J ) = IC +  NCOLS * ( IR - 1 )
                    AFRAC( J ) = MIN( 1.0, MAX( 0.0, FAC ) )
                    GO TO  133              !  to head of loop
                END IF
            END IF

C...........   Merge complete.  Return the number of cells found:

        NCEL = J
        RETURN

C ********************** INTERNAL SUBPROGRAMS ****************************

        CONTAINS

            SUBROUTINE REPORT_BAD_CELL

C.............  Local variables
            CHARACTER(300) MESG 

C..........................................................................

            WRITE( MESG,94010 ) 'INTERNAL ERROR: Bad number of cells ', 
     &             NCEL, 'computed in ' // PROGNAME
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

            END SUBROUTINE REPORT_BAD_CELL

        END SUBROUTINE LNK2GRD


