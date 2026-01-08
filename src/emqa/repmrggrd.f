
        SUBROUTINE REPMRGGRD( RCNT, NX, IX, CX, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       The REPMRGGRD routine assigns grid cells and gridding factors to each 
C       of the sources still selected by the program. If there is a subgrid,
C       it will be used to subselect the grid cells.  The final OUTREC array
C       is created by this routine before grouping, summing, and writing
C       emissions data.  If the report contains normalization by grid cell
C       these factors will be included in the gridding factors.
C 
C  PRECONDITIONS REQUIRED:
C       Routine is only called if gridding is being used
C       Gridding matrix is allocated and populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: INDEXA

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: NMATX, LSUBGRID, ALLRPT, NSUBGRID, SUBGNAM,
     &                      VALIDCEL

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: OUTSRC, OUTBIN, OUTCELL, OUTGFAC,
     &                      NOUTREC, NSRCDROP

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, GDTYP, XCELL, YCELL, NCOLS, YORIG

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

C...........  EXTERNAL FUNCTIONS and their descriptions:
C       INTEGER     INDEX1

C        EXTERNAL    INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: RCNT        ! current report number
        INTEGER, INTENT (IN) :: NX( NGRID ) ! no. srcs per cell
        INTEGER, INTENT (IN) :: IX( NMATX ) ! src IDs
        REAL   , INTENT (IN) :: CX( NMATX ) ! gridding coefficients
        LOGICAL, INTENT(OUT) :: EFLAG       ! true: error found

C...........  Local allocatable arrays
        REAL, ALLOCATABLE, SAVE :: NORMFAC( : )

C...........  Local variables
        INTEGER         C          ! indices and counters
        INTEGER         IOS        ! i/o status
        INTEGER         ROW        ! tmp row number

        REAL            FAC        ! temporary factor
        REAL            Y          ! tmp Y cell center
        REAL            Y0         ! y origin less 0.5*dy

        LOGICAL      :: FIRSTFLAG = .TRUE.  ! True: no normalize by cell area yet found

        CHARACTER(300)  MESG       ! message buffer

        CHARACTER(16) :: PROGNAME = 'REPMRGGRD' ! program name

C***********************************************************************
C   begin body of subroutine REPMRGGRD

C.........  For first time a normalized grid is encountered by this routine,
C           compute the cell-area normalization factors for full grid.
        IF ( ALLRPT( RCNT )%NORMCELL .AND. FIRSTFLAG ) THEN

            ALLOCATE( NORMFAC( NGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NORMFAC', PROGNAME )

            SELECT CASE( GDTYP )
            CASE( LAMGRD3, UTMGRD3 )
                FAC   = 1. / ( XCELL * YCELL )
                NORMFAC = FAC   ! array

            CASE( LATGRD3 )

                Y0 = YORIG - 0.5 * YCELL
                FAC = ( PI * REARTH / 180. )**2 * XCELL * YCELL
                DO C = 1, NGRID
                    ROW = 1 + INT( (C-1) / NCOLS )
                    Y = Y0 + REAL( ROW ) * YCELL
                    NORMFAC( C ) = 1. / ( COS( Y ) * FAC )
                END DO

            CASE DEFAULT

                WRITE( MESG,94010 ) 'INTERNAL ERROR: Grid type', GDTYP,
     &                 'not recognized in ' // PROGNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END SELECT

            FIRSTFLAG = .FALSE.

        END IF

C.........  Set report-specific local settings
        LSUBGRID = ( ALLRPT( RCNT )%SUBGNAM .NE. ' ' )

C.........  Deallocate source list and bins before reallocating for source-cell 
C           intersections.  Valid sources are identified with INDEXA != 0.
        IF( ALLOCATED( OUTSRC ) ) DEALLOCATE( OUTSRC, OUTBIN )

C.........  Deallocate cell and gridding information if they have been
C           allocated for a previous report
        IF( ALLOCATED( OUTCELL ) ) DEALLOCATE( OUTCELL, OUTGFAC )

C.........  Count cell-source intersections...
        CALL SOURCE_CELL_X( 'COUNT', NX, IX, CX )

C.........  Warning if no source-cell intersections
        IF( NOUTREC .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = BLANK5 // 'ERROR: No source-cell intersections '//
     &             'for report.  Output will be empty.'
            CALL M3MSG2( MESG )

        END IF

C.........  Allocate memory for output record arrays
        ALLOCATE( OUTSRC( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTSRC', PROGNAME )
        ALLOCATE( OUTBIN( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTBIN', PROGNAME )
        ALLOCATE( OUTCELL( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTCELL', PROGNAME )
        ALLOCATE( OUTGFAC( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTGFAC', PROGNAME )

C.........  Store sources, cells, and gridding factors
        CALL SOURCE_CELL_X( 'STORE', NX, IX, CX )        

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal function counts or computes the source-cell
C               intersections for a whole grid or subgrid, and for specifically-
C               selected sources from the inventory.
            SUBROUTINE SOURCE_CELL_X( COMMAND, NX, IX, CX )

C.............  Subprogram arguments (note - array dimensions from MODREPRT)
            CHARACTER(*), INTENT (IN) :: COMMAND     ! "COUNT" or "STORE"
            INTEGER     , INTENT (IN) :: NX( NGRID ) ! no. srcs per cell
            INTEGER     , INTENT (IN) :: IX( NMATX ) ! src IDs
            REAL        , INTENT (IN) :: CX( NMATX ) ! gridding coefficients

C.............  Local variables
            INTEGER   C, IC, IG, J, K, N, S

C----------------------------------------------------------------------

C.............  Special case of no subgrid and all sources selected
C.............  Have this case explicitly because it will not be slowed down
C               by uncessary counting and conditionals.
            IF( NSRCDROP .EQ. 0 .AND. .NOT. LSUBGRID ) THEN

                NOUTREC = NMATX

                IF( COMMAND .EQ. 'COUNT' ) RETURN

                K = 0 
                DO C = 1, NGRID

                    DO J = 1, NX( C )

                        K = K + 1
                        OUTCELL( K ) = C
                        OUTSRC ( K ) = IX( K )
                        OUTGFAC( K ) = CX( K )

                    END DO

                END DO

C.............  Special case of no subgrid and partial sources
            ELSE IF ( .NOT. LSUBGRID ) THEN

                K = 0
                N = 0 
                DO C = 1, NGRID

                    DO J = 1, NX( C )

                        K = K + 1
                        S = IX( K )

C.........................  Skip source if it has not been selected
                        IF( INDEXA( S ) .LE. 0 ) CYCLE

C.........................  Increment count of output records
                        N = N + 1

C.........................  If only counting, then skip the rest of the loop
                        IF( COMMAND .EQ. 'COUNT' ) CYCLE

C.........................  Store the factors for the current cell and source
                        OUTCELL( N ) = C
                        OUTSRC ( N ) = S
                        OUTGFAC( N ) = CX( K )

                    END DO

                END DO

                NOUTREC = N     ! Store output count of source-cell interstns

C.............  Process for a subgrid and/or partial source list
            ELSE

C.................  Determine subgrid index based on name (the validity has 
C                   already been confirmed elsewhere).
                IG = INDEX1( ALLRPT( RCNT )%SUBGNAM, NSUBGRID, SUBGNAM )

                K  = 0 
                N  = 0
                IC = 1
                DO C = 1, NGRID

C.....................  If current cell is valid for the subgrid...
                    IF( C .EQ. VALIDCEL( IC,IG ) ) THEN

C........................  Loop through sources for current cell
                       DO J = 1, NX( C )

                            K = K + 1
                            S = IX( K )

C.............................  Skip source if it has not been selected
                            IF( INDEXA( S ) .LE. 0 ) CYCLE

C.............................  Increment count of output records
                            N = N + 1

C.............................  If only counting, then skip the rest of the loop
                            IF( COMMAND .EQ. 'COUNT' ) CYCLE

C.............................  Store the factors for the current cell and src
                            OUTCELL( N ) = C
                            OUTSRC ( N ) = S
                            OUTGFAC( N ) = CX( K )

                        END DO  ! End loop over sources for current cell

                        IC = IC + 1

C.....................  Otherwise, step through gridding matrix
                    ELSE

                        K = K + NX( C )

                    END IF      ! Check if valid cell

                END DO          ! End loop over cells

                NOUTREC = N     ! Store output count of source-cell interstns

            END IF              ! If ( all sources and full grid ) or not

C.............  If needed, update all gridding factors with division by cell
C               area.
            IF ( ALLRPT( RCNT )%NORMCELL ) THEN

                DO N = 1, NOUTREC
                    C = OUTCELL( N )
                    OUTGFAC( N ) = OUTGFAC( N ) * NORMFAC( C )
                ENDDO

            END IF

            RETURN
 
            END SUBROUTINE SOURCE_CELL_X

        END SUBROUTINE REPMRGGRD


