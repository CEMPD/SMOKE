
        SUBROUTINE REPMRGGRD( RCNT, NX, IX, CX, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       The REPMRGGRD routine assigns grid cells and gridding factors to each 
C       of the sources still selected by the program. If there is a subgrid,
C       it will be used to subselect the grid cells.  The final OUTREC array
C       is created by this routine before grouping, summing, and writing
C       emissions data.
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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     INDEX1

        EXTERNAL    INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: RCNT        ! current report number
        INTEGER, INTENT (IN) :: NX( NGRID ) ! no. srcs per cell
        INTEGER, INTENT (IN) :: IX( NMATX ) ! src IDs
        REAL   , INTENT (IN) :: CX( NMATX ) ! gridding coefficients
        LOGICAL, INTENT(OUT) :: EFLAG       ! true: error found

C...........  Local variables

        INTEGER         IOS        ! i/o status

        CHARACTER*300   MESG       ! message buffer

        CHARACTER*16 :: PROGNAME = 'REPMRGGRD' ! program name

C***********************************************************************
C   begin body of subroutine REPMRGGRD

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

                    IF( C .EQ. VALIDCEL( IC,IG ) ) THEN

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

                    END IF      ! Check if valid cell

                END DO          ! End loop over cells

                NOUTREC = N     ! Store output count of source-cell interstns

            END IF              ! If ( all sources and full grid ) or not

            RETURN
 
            END SUBROUTINE SOURCE_CELL_X

        END SUBROUTINE REPMRGGRD


