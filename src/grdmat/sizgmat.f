
        SUBROUTINE SIZGMAT( CATEGORY, NSRC, MXSCEL, MXCSRC, NMATX )

C***********************************************************************
C  subroutine body starts at line 102
C
C  DESCRIPTION:
C      This subroutine determines the sizes needed for creating the for the 
C      area and mobile gridding matrices
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 5/99
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C...........   This module contains the gridding surrogates tables
        USE MODSURG

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL   INGRID
        EXTERNAL  INGRID

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  !  source category
        INTEGER     , INTENT (IN) :: NSRC      !  local number of sources
        INTEGER     , INTENT(OUT) :: MXSCEL    !  max sources per cell   
        INTEGER     , INTENT(OUT) :: MXCSRC    !  max cells per source   
        INTEGER     , INTENT(OUT) :: NMATX     !  no. src-cell intersections   

C...........   Local arrays dimensioned by subroutine arguments
C...........   Note that the NGRID dimension could conceivably be too small if 
C              a link winds through the whole domain, but this is a case that
C              is not worth going to extra trouble for since it is not realistic
        INTEGER         NX  ( NGRID )    ! number of srcs per cell  
        INTEGER         ACEL( NGRID )    ! number of cell intersections per src
        REAL            AFAC( NGRID )    ! fraction of link in cell

C...........   Other local variables
        INTEGER         C, F, J, K, I, N, S          ! counters and indices

        INTEGER         CCNT             ! counters for no. non-zero-surg cells
        INTEGER      :: CELLSRC = 0      ! cell number as source char
        INTEGER         COL              ! tmp column
        INTEGER         ID1, ID2         ! primary and 2ndary surg codes
        INTEGER         ISIDX            ! tmp surrogate ID code index
        INTEGER         NCEL             ! tmp number of cells
        INTEGER         ROW              ! tmp row

        REAL            ALEN        ! link length

        LOGICAL      :: EFLAG = .FALSE. ! true: error flag
        LOGICAL      :: LFLAG = .FALSE. ! true: location data available
        LOGICAL      :: XYSET = .FALSE. ! true: X/Y available for src

        CHARACTER*256   MESG        ! message buffer

        CHARACTER(LEN=LNKLEN3) :: CLNK = ' '   ! tmp link ID

        CHARACTER*16 :: PROGNAME = 'SIZGMAT' ! program name

C***********************************************************************
C   begin body of subroutine SIZGMAT

C.........  Print status message
        MESG = 'Computing gridding matrix size...'
        CALL M3MSG2( MESG )

C.........  Set flag to indicate that XLOCA/YLOCA are available
        LFLAG = ALLOCATED( XLOCA )

C.........  Initialize the count of sources per cell
        NX = 0   ! array

C.........  Loop through sources
        MXCSRC  = 0
        CELLSRC = 0
        DO S = 1, NSRC

            IF( CATEGORY .EQ. 'AREA' ) CELLSRC = CELLID( S )
            IF( CATEGORY .EQ. 'MOBILE' ) CLNK = CLINK( S )

C.............  Determine if x/y location is available
            XYSET = .FALSE.
            IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

C.............  If cell-specific source...
            IF ( CELLSRC .GT. 0 ) THEN
                NCEL = 1
                ACEL( 1 ) = CELLID( S )
                AFAC( 1 ) = 1.

C.............  Check if source has been converted to point src
            ELSE IF( XYSET ) THEN

C................  If source is in the domain....
                IF( INGRID( XLOCA( S ), YLOCA( S ), 
     &                      NCOLS, NROWS, COL, ROW  ) ) THEN

C....................  Set as 1 cell and get the cell number
                    NCEL = 1
                    ACEL( 1 ) = ( ROW-1 ) * NCOLS + COL
                    AFAC( 1 ) = 1.

                END IF

C............  If area/non-link source...
            ELSE IF( CLNK .EQ. ' ' ) THEN

C.................  Retrieve the index to the surrogates cy/st/co list
                ISIDX = SRGIDPOS( S )
                F     = SGFIPPOS( S )

C.................  Retrieve the cell intersection info from the
C                   surrogates tables from MODSURG
                IF ( F .GT. 0 ) THEN

                    NCEL = NCELLS( F )
                    ACEL( 1:NCEL ) = FIPCELL( 1:NCEL, F )      ! arrays

                    DO K = 1, NCEL
                        CALL SETFRAC( S, ISIDX, K, F, 1, .FALSE., 
     &                                ' ', ID1, ID2, AFAC( K ) )
                    END DO

C.................  Otherwise, skip this source because it's outside the grid
                ELSE  
                    NCEL = 0

                END IF

C.............  If link source, determine the number of cells for this source
            ELSE

                CALL LNK2GRD( NGRID, XLOC1( S ), YLOC1( S ), XLOC2( S ), 
     &                        YLOC2( S ), NCEL, ACEL, AFAC, ALEN, EFLAG)

C.................  Make sure that there was enough storage 
                IF ( EFLAG ) THEN
                    WRITE( MESG,94010 )
     &                  'INTERNAL ERROR: Overflow for source', S
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

            END IF
C
C.............  Loop through the cells for this source and increment the number
C               of sources per cell. 
            CCNT = 0
            DO N = 1, NCEL

                IF( AFAC( N ) .GT. 0. ) THEN
                    C = ACEL( N )
                    NX( C ) = NX( C ) + 1
                    CCNT = CCNT + 1
                END IF

            END DO    ! End loop on cells for this source

C.............  Update the maximum number of cells per source
            IF( CCNT .GT. MXCSRC ) MXCSRC = CCNT
 
        END DO        ! End loop on sources

C.........  Abort if error
        IF( EFLAG ) THEN
            MESG = 'Problem determining memory for gridding matrix.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Determine maximum number of sources per cell
C.........  And determine the total number of source-cell intersections
        MXSCEL = NX( 1 )
        NMATX  = NX( 1 )          
        DO C = 2, NGRID
            
            J = NX( C )
            IF( J .GT. MXSCEL ) THEN
                MXSCEL = J
            END IF

            NMATX = NMATX + J

        END DO        ! End loop on cells

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SIZGMAT
