
        LOGICAL FUNCTION INGRID( XX, YY, NCOLS, NROWS, COL, ROW )

C***************************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function determines if coordinates are inside the grid or not.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 9/2000
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
C***************************************************************************
        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL    DSCM3GRD
        EXTERNAL   DSCM3GRD

C...........   SUBROUTINE ARGUMENTS
        REAL   , INTENT (IN) :: XX    ! X-coordinate in projection
        REAL   , INTENT (IN) :: YY    ! Y-ccordinate in projection
        INTEGER, INTENT(OUT) :: NCOLS ! no. of columns
        INTEGER, INTENT(OUT) :: NROWS ! no. of rows
        INTEGER, INTENT(OUT) :: COL   ! column number
        INTEGER, INTENT(OUT) :: ROW   ! row number

C...........   Local variables

        INTEGER          C, I  !  indices and counters.
        INTEGER, SAVE :: NC    !  saved no. columns
        INTEGER, SAVE :: NR    !  saved no. rows

        REAL             XDUM, YDUM ! tmp X and Y coordinates
        REAL, SAVE    :: XX0, YY0   ! X and Y origin in coordinates of grid
        REAL, SAVE    :: XX1, YY1   ! X and Y upper right position of grid
        REAL, SAVE    :: DDX, DDY   ! Inverse cell length in X and Y directions

        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
                                 
        CHARACTER(16)   COORD     !  coordinate system name
        CHARACTER(16)   COORUNIT  !  coordinate system projection units
        CHARACTER(16)   GRDNM     !  grid name
        CHARACTER(80)   GDESC     !  grid description
        CHARACTER(300)  MESG      !  message buffer 

        CHARACTER(16) :: PROGNAME = 'INGRID' ! program name

C***********************************************************************
C   begin body of function INGRID

C.........  Get grid name from the environment and read grid parameters
        IF( FIRSTIME ) THEN

            IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT, 
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

                MESG = 'Could not get Models-3 grid description'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.............  Store grid parameters for later processing
            ELSE

                XX0 = SNGL( XORIG3D )
                YY0 = SNGL( YORIG3D )
                XX1 = XX0 + FLOAT( NCOLS3D ) * SNGL( XCELL3D )
                YY1 = YY0 + FLOAT( NROWS3D ) * SNGL( YCELL3D )
                DDX = 1.0 / SNGL( XCELL3D )
                DDY = 1.0 / SNGL( YCELL3D )
                NR  = NROWS3D
                NC  = NCOLS3D

            END IF

            FIRSTIME = .FALSE.

        END IF

        NROWS = NR
        NCOLS = NC

C.........  Initialize function
        INGRID = .TRUE.

C.........  Check X starting coordinate            
        XDUM = DDX * ( XX - XX0 )
        IF ( XDUM .LT. 0.0  )  THEN
            INGRID = .FALSE.
            RETURN
        END IF

C.........  Check X cell            
        COL = 1 + INT( XDUM )
        IF ( COL .GT. NCOLS .OR. COL .LE. 0 )  THEN
            INGRID = .FALSE.
            RETURN
        END IF

C.........  Check Y starting coordinate            
        YDUM = DDY * ( YY - YY0 )
        IF ( YDUM .LT. 0.0  )  THEN
            INGRID = .FALSE.
            RETURN
        END IF

C.........  Check Y cell            
        ROW = 1 + INT( YDUM )
        IF ( ROW .GT. NROWS .OR. ROW .LE. 0 )  THEN
            INGRID = .FALSE.
            RETURN
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END FUNCTION INGRID

