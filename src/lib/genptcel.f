
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

        USE MODGRDLIB
        USE M3UTILIO

        IMPLICIT NONE

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

C...........   Local variables

        INTEGER         C, I  !  indices and counters.

        INTEGER         COL   ! tmp column number
        INTEGER         ROW   ! tmp row number
        INTEGER         CSAV  ! saved value of C
        INTEGER      :: NROWS = 0  ! No. of rows, columns, and cells
        INTEGER      :: NCOLS = 0  ! No. of rows, columns, and cells

        CHARACTER(16)   COORD     !  coordinate system name
        CHARACTER(16)   COORUNIT  !  coordinate system projection units
        CHARACTER(16)   GRDNM     !  grid name
        CHARACTER(80)   GDESC     !  grid description
        CHARACTER(300)  MESG      !  message buffer 

        CHARACTER(16) :: PROGNAME = 'GENPTCEL' ! program name

C***********************************************************************
C   begin body of subroutine GENPTCEL

C.........  Initialize number of sources per cell
        NX = 0   ! array

C.........  Initialize scratch gridding matrix - before sparse storage
        NEXCLD = 0
        CSAV = 0
        DO I = 1, NRECS
            GN  ( I ) = 0
            SN  ( I ) = 0
            INDX( I ) = I

            IF( .NOT. INGRID( XLOCA( I ), YLOCA( I ), NCOLS, NROWS, COL, ROW ) ) THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            END IF
                         
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
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 3 ) 
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GENPTCEL

