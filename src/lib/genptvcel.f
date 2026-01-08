
        SUBROUTINE GENPTVCEL( NRECS, NGRID, XLOCA, YLOCA, NEXCLD,
     &                        NX, INDX, GN, SN )

C***************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine creates the scratch gridding matrix, which contains
C      the cell IDs for each source using a variable grid.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen 7/04 based on genptcel.f
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

C.........  This module contains the global variables for the 3-d grid
        USE M3UTILIO

        USE MODGRID, ONLY: GRDNM, GDTYP, NCOLS, NROWS, 
     &                     P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                     XORIG, YORIG, XCELL, YCELL

        USE MODGRDLIB
 
        IMPLICIT NONE

C...........   INCLUDES
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(16), EXTERNAL :: PROMPTMFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NRECS          ! no. records w/ coordinates
        INTEGER, INTENT (IN) :: NGRID          ! no. grid cells
        REAL(8), INTENT (IN) :: XLOCA( NRECS ) ! X-coordinate in projection
        REAL(8), INTENT (IN) :: YLOCA( NRECS ) ! Y-ccordinate in projection
        INTEGER, INTENT(OUT) :: NEXCLD         ! no. sources excluded
        INTEGER, INTENT(OUT) :: NX   ( NGRID ) ! no. of sources per cell
        INTEGER, INTENT(OUT) :: INDX ( NRECS ) ! sorting index
        INTEGER, INTENT(OUT) :: GN   ( NRECS ) ! cell number
        INTEGER, INTENT(OUT) :: SN   ( NRECS ) ! record number

C...........   Allocatable arrays
        REAL*8, ALLOCATABLE :: XVALS( :,: ) ! x values for grid cell boundaries
        REAL*8, ALLOCATABLE :: YVALS( :,: ) ! y values for grid cell boundaries

C...........   Local variables

        INTEGER         C, I, J  !  indices and counters

        INTEGER         COL   ! tmp column number
        INTEGER         ROW   ! tmp row number
        INTEGER         CSAV  ! saved value of C
        INTEGER         IOS   ! i/o status
        INTEGER         NCDOT ! number of columns in dot file
        INTEGER         NRDOT ! number of rows in dot file

        REAL            STEP      ! distance to edge of column or row
        REAL*8          XX, YY    ! tmp X and Y coordinates
        REAL            XMIN      ! minimum location for column
        REAL            XMAX      ! maximum location for column
        REAL            YMIN      ! minimum location for row
        REAL            YMAX      ! maximum location for row

        LOGICAL ::      EFLAG = .FALSE. ! true: an error has occured

        CHARACTER(16)   GNAME     !  logical file name for GRID_DOT_2D file
        CHARACTER(300)  MESG      !  message buffer 

        CHARACTER(16) :: PROGNAME = 'GENPTVCEL' ! program name

C***********************************************************************
C   begin body of subroutine GENPTVCEL

C.........  Open 2-D grid parameters file to get cell coordinates
C           Use the GRID_DOT_2D file to get cell boundaries
        GNAME = PROMPTMFILE(
     &          'Enter name for DOT-POINT SURFACE GRID file',
     &          FSREAD3, 'GRID_DOT_2D', PROGNAME )
    
        IF( .NOT. DESC3( GNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( GNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check grid against previously set grid description
        CALL CHKGRID( GNAME, 'DOT', 0, EFLAG )
        
        IF( EFLAG ) THEN
            MESG = TRIM( GNAME ) // ' does not match previously set ' //
     &             'grid description'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        NCDOT = NCOLS + 1
        NRDOT = NROWS + 1
        
C.........  Allocate memory for grid cell coordinates
        ALLOCATE( XVALS( NCDOT, NRDOT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XVALS', PROGNAME )
        ALLOCATE( YVALS( NCDOT, NRDOT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YVALS', PROGNAME )
        
C.........  Read grid cell coordinates
        IF( .NOT. READ3( GNAME, 'LON', 1, 0, 0, XVALS ) ) THEN
            MESG = 'Could not read LON from file "' //
     &             TRIM( GNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        IF( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, YVALS ) ) THEN
            MESG = 'Could not read LAT from file "' //
     &             TRIM( GNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Convert coordinates to map projection units
        CALL CONVRTXY( NCDOT, NRDOT, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
     &                 XCENT, YCENT, XVALS, YVALS )

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

            IF( .NOT. INGRID( XLOCA(I), YLOCA(I), NCOLS, NROWS, COL, ROW ) ) THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            END IF

C.............  Find correct column for point
            DO J = 1, NCOLS
            
                IF( XX >= XVALS( J,1 ) .AND. XX <= XVALS( J+1,1 ) ) EXIT

            END DO

            IF( J <= NCOLS ) THEN 
                COL = J
            END IF
            
C.............  Find correct row for point
            DO J = 1, NROWS
            
                IF( YY >= YVALS( 1,J ) .AND. YY <= YVALS( 1,J+1 ) ) EXIT
                
            END DO

            IF( J <= NROWS ) THEN 
                ROW = J
            END IF

C.............  Compute grid cell number based on column and row            
            C = COL + NCOLS * ( ROW - 1 )

            IF( C .LE. NGRID ) THEN
                NX( C ) = NX( C ) + 1

            ELSE IF( C .GT. CSAV ) THEN
                CSAV = C

            END IF

            GN( I ) = C
            SN( I ) = I

        END DO        !  end loop on sources I, computing gridding matrix.
        
        IF( NCOLS * NROWS .NE. NGRID ) THEN

            MESG = 'INTERNAL ERROR: Number of cells in "' //
     &                 PROGNAME( 1:16 ) // '" are inconsistent '//
     &                 'with calling program'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

        IF ( CSAV .GT. NGRID ) THEN
            WRITE( MESG,94010 ) 
     &             'INTERNAL ERROR: Number of grid cells C=', CSAV,
     &             'exceeds NGRID=', NGRID
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 ) 
        END IF

C......... Deallocate local memory
    	DEALLOCATE( XVALS, YVALS )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GENPTVCEL

