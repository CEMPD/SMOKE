
        SUBROUTINE GENPTCEL( NRECS, NGRID, XLOCA, YLOCA, NEXCLD,
     &                       NX, INDX, GN, SN )

C***************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine creates the scratch gridding matrix, which contains
C      the cell IDs for each source.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 8/99
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL    DSCM3GRD
        EXTERNAL   DSCM3GRD

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NRECS          ! no. records w/ coordinates
        INTEGER, INTENT (IN) :: NGRID          ! no. grid cells
        REAL   , INTENT (IN) :: XLOCA( NRECS ) ! X-coordinate in projection
        REAL   , INTENT (IN) :: YLOCA( NRECS ) ! Y-ccordinate in projection
        INTEGER, INTENT(OUT) :: NEXCLD         ! no. sources excluded
        INTEGER, INTENT(OUT) :: NX   ( NGRID ) ! no. of sources per cell
        INTEGER, INTENT(OUT) :: INDX ( NRECS ) ! sorting index
        INTEGER, INTENT(OUT) :: GN   ( NRECS ) ! cell number
        INTEGER, INTENT(OUT) :: SN   ( NRECS ) ! record number

C...........   Local variables

        INTEGER         C, I  !  indices and counters.

        INTEGER         COL   ! tmp column number
        INTEGER         ROW   ! tmp row number
        INTEGER         CSAV  ! saved value of C
        INTEGER         NROWS, NCOLS  ! No. of rows, columns, and cells

        REAL            XX, YY, XDUM, YDUM ! tmp X and Y coordinates
        REAL            XX0, YY0  ! X and Y origin in coordinates of grid
        REAL            XX1, YY1  ! X and Y upper right position of grid
        REAL            DDX, DDY  ! Inverse cell length in X and Y directions
                                 
        CHARACTER*16    COORD     !  coordinate system name
        CHARACTER*16    COORUNIT  !  coordinate system projection units
        CHARACTER*16    GRDNM     !  grid name
        CHARACTER*80    GDESC     !  grid description
        CHARACTER*300   MESG      !  message buffer 

        CHARACTER*16 :: PROGNAME = 'GENPTCEL' ! program name

C***********************************************************************
C   begin body of subroutine GENPTCEL

C.........  Get grid name from the environment and read grid parameters
        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT, 
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Store grid parameters for later processing
        ELSE

            XX0   = SNGL( XORIG3D )
            YY0   = SNGL( YORIG3D )
            XX1   = XX0 + FLOAT( NCOLS3D ) * SNGL( XCELL3D )
            YY1   = YY0 + FLOAT( NROWS3D ) * SNGL( YCELL3D )
            DDX   = 1.0 / SNGL( XCELL3D )
            DDY   = 1.0 / SNGL( YCELL3D )
            NROWS = NROWS3D
            NCOLS = NCOLS3D

            IF( NCOLS * NROWS .NE. NGRID ) THEN
 
                MESG = 'INTERNAL ERROR: Number of cells in "' //
     &                 PROGNAME( 1:16 ) // '" are inconsistent '//
     &                 'with calling program'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF

        END IF

C.........  Initialize number of sources per cell
        NX = 0   ! array

C.........  Initialize scratch gridding matrix - before sparse storage
        NEXCLD = 0
        CSAV = 0
        DO I = 1, NRECS
            GN  ( I ) = 0
            SN  ( I ) = 0
            INDX( I ) = I

            XX = XLOCA( I )
            YY = YLOCA( I ) 
            
            XDUM = DDX * ( XX - XX0 )
            IF ( XDUM .LT. 0.0  )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            END IF

            COL = 1 + INT( XDUM )
            IF ( COL .GT. NCOLS .OR. COL .LE. 0 )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            END IF

            YDUM = DDY * ( YY - YY0 )
            IF ( YDUM .LT. 0.0  )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            END IF

            ROW = 1 + INT( YDUM )
            IF ( ROW .GT. NROWS .OR. ROW .LE. 0 )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            END IF
            
            C = COL  +  NCOLS * ( ROW - 1 )

            IF( C .LE. NGRID ) THEN
                NX( C ) = NX( C ) + 1

            ELSEIF( C .GT. CSAV ) THEN
                CSAV = C

            END IF

            GN( I ) = C
            SN( I ) = I

        END DO        !  end loop on sources I, computing gridding matrix.
        
        IF ( CSAV .GT. NGRID ) THEN
            WRITE( MESG,94010 ) 
     &             'INTERNAL ERROR: Number of grid cells C=', CSAV,
     &             'exceeds NGRID=', NGRID
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 3 ) 
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GENPTCEL

