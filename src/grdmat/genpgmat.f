
        SUBROUTINE GENPGMAT( FNAME, NPSRC, NGRID, XLOCA, YLOCA,
     &                       NX, IX, NCOEF, CMAX, CMIN )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine creates the point source gridding matrix and
C      writes it to a NetCDF file.  This subroutine is created to be
C      consisent with those for area and mobile sources, which are needed
C      to ensure that the sparse matrix is stored contiguously in memory. 
C
C  PRECONDITIONS REQUIRED:
C      File must be opened and its logical name input through the FNAME
C      argument.  Memory for NX and IX must be allocated prior to the 
C      subroutine call to ensure that it is contiguous. The x and y coordinates
C      must already be converted to the output grid.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
        LOGICAL         DSCM3GRD
        INTEGER         TRIMLEN

        EXTERNAL DSCM3GRD, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=NAMLEN3) FNAME ! matrix output inventory logical name
        INTEGER       NPSRC          ! Actual source count
        INTEGER       NGRID          ! Actual grid cell count
        REAL          XLOCA( NPSRC ) ! X-coordinate in proper coordinate system
        REAL          YLOCA( NPSRC ) ! Y-coordinate in proper coordinate system
        INTEGER       NX( NGRID )    ! Number of sources per cell
        INTEGER       IX( NPSRC )    ! Source number, w/ order based on NX
        INTEGER       NCOEF          ! Number of coefficients
        INTEGER       CMAX           ! Maximum number of sources per cell
        INTEGER       CMIN           ! Minimum number of sources per cell

C...........   Scratch Gridding Matrix

        INTEGER      INDX( NPSRC )
        INTEGER      GN  ( NPSRC )
        INTEGER      SN  ( NPSRC )

C...........   Other local variables

        INTEGER         C, I, J, K, K2, R !  indices and counters.

        INTEGER         COL   ! tmp column number
        INTEGER         ROW   ! tmp row number
        INTEGER         CSAV  ! saved value of C
        INTEGER         IOS   ! i/o status
        INTEGER         KSAV  ! saved value of K
        INTEGER         K2SAV ! saved value of K2
        INTEGER         NEXCLD! number of sources not in grid.  Will
                              !    appear as blanks in SN( ) and GN( )
        INTEGER         NROWS, NCOLS  ! No. of rows, columns, and cells

        REAL            XX, YY, XDUM, YDUM ! tmp X and Y coordinates
        REAL            XX0, YY0  ! X and Y origin in coordinates of grid
        REAL            XX1, YY1  ! X and Y upper right position of grid
        REAL            DDX, DDY  ! Inverse cell length in X and Y directions
                                 
        LOGICAL         GFLAG     !  generate output gridding matrix if true

        CHARACTER*16    COORD     !  coordinate system name
        CHARACTER*16    COORUNIT  !  coordinate system projection units
        CHARACTER*16    GRDNM     !  grid name
        CHARACTER*80    GDESC     !  grid description
        CHARACTER*300   MESG      !  message buffer 

        CHARACTER*16 :: PROGNAME = 'GENPGMAT' ! program name

C***********************************************************************
C   begin body of subroutine GENPGMAT

C.........  Get grid name from the environment and read grid parameters
        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT, 
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'ERROR: Could not get Models-3 grid description'
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

            ENDIF

        ENDIF

C.........  Initialize number of sources per cell
        DO C = 1, NGRID
            NX( C ) = 0
        ENDDO

C.........  Initialize temporary storage table
        NEXCLD = 0
        CSAV = 0
        DO I = 1, NPSRC
            GN  ( I ) = 0
            SN  ( I ) = 0
            INDX( I ) = I

            XX = XLOCA( I )
            YY = YLOCA( I ) 
            
            XDUM = DDX * ( XX - XX0 )
            IF ( XDUM .LT. 0.0  )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            ENDIF

            COL = 1 + INT( XDUM )
            IF ( COL .GT. NCOLS )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            ENDIF

            YDUM = DDY * ( YY - YY0 )
            IF ( YDUM .LT. 0.0  )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            ENDIF

            ROW = 1 + INT( YDUM )
            IF ( ROW .GT. NROWS )  THEN
                NEXCLD = NEXCLD + 1
                CYCLE  ! To end of loop
            ENDIF
            
            C = COL  +  NCOLS * ( ROW - 1 )

            IF( C .LE. NGRID ) THEN
                NX( C ) = NX( C ) + 1

            ELSEIF( C .GT. CSAV ) THEN
                CSAV = C

            ENDIF

            GN( I ) = C
            SN( I ) = I

        ENDDO        !  end loop on sources I, computing gridding matrix.
        
        IF ( CSAV .GT. NGRID ) THEN
            WRITE( MESG,94010 ) 
     &             'INTERNAL ERROR: Number of grid cells C=', CSAV,
     &             'exceeds NGRID=', NGRID
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 3 ) 
        END IF

        CMAX = NX( 1 )
        CMIN = CMAX

C.........  Sort temporary gridding matrix
        CALL SORTI2( NPSRC, INDX, GN, SN )
            
        GFLAG = ( FNAME( 1:5 ) .NE. 'NONE ' )        

C...........   Compute gridding matrix and/or compute statistics

        K     = 0
        KSAV  = 0
        K2SAV = 0
        IF ( GFLAG ) THEN          ! write gridding matrix to file

            CALL M3MSG2( 'Computing gridding matrix and statistics...' )

C.............   Compress matrix into I/O representation from scratch
C                representation.
C.............   Compute statistics

            DO R = 1, NGRID
            
                J = NX( R )
                   
                IF( J .GT. CMAX ) THEN
                    CMAX = J
                ELSEIF( J .LT. CMIN ) THEN
                    CMIN = J
                ENDIF

                IF( J .NE. 0 ) THEN
                   
                    DO C = 1, J
                        K  = K + 1
                        K2 = K + NEXCLD

                        IF( K .LE. NPSRC .AND. K2 .LE. NPSRC )
     &                      IX( K ) = SN( INDX( K2 ) )

                        IF( K  .GT. KSAV )  KSAV  = K
                        IF( K2 .GT. K2SAV ) K2SAV = K2

                    ENDDO

                ENDIF
                   
            ENDDO    !  end of loop on cells K for this FIP
            
C.............  Give error(s) if memory allocation exceeded 
            IF( K .GT. NPSRC ) THEN
                WRITE( MESG,94010 ) 
     &                  'INTERNAL ERROR: Number of gridding ' //
     &                  'coefficients K=', K, 
     &                  'exceeds NPSRC=', NPSRC
                CALL M3MSG2( MESG ) 
            ENDIF
                   
            IF( K2 .GT. NPSRC ) THEN
                WRITE( MESG,94010 ) 
     &                  'INTERNAL ERROR: Number of gridding ' //
     &                  'coefficients K2=', K2, 
     &                  'exceeds NPSRC=', NPSRC
                CALL M3MSG2( MESG ) 
            ENDIF

C.............  Stop program if memory exceedance errors were written
            IF( K .GT. NPSRC .OR. K2 .GT. NPSRC ) 
     &          CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 ) 

            CALL M3MSG2( 'Writing out GRIDDING MATRIX file...' )

            IF( .NOT. WRITE3( FNAME, 'ALL', 0, 0, NX ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 
     &              'Error writing GRIDDING MATRIX file.', 2 )
            ENDIF

        ELSE            !....  not GFLAG:  report matrix stats only
                            
            CALL M3MSG2( 'Computing gridding statistics...' )

            DO R = 1, NGRID  ! Must start loop at 1 to get correct K value
            
                J = NX( R )
                K = K + J
                   
                IF( J .GT. CMAX ) THEN
                    CMAX = J
                ELSEIF ( J .LT. CMIN ) THEN
                    CMIN = J
                ENDIF

            ENDDO    !  end of loop on cells K for this FIP

        ENDIF  !  if gflag, or not

        NCOEF = K

        WRITE( MESG,94010 ) 
     &      'NOTE: Number of sources excluded from grid was', NEXCLD

        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END

