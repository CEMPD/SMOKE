
        SUBROUTINE BLDSRCCELL( NSRC, NGRID, NMAT, NX, IX, CX )

C***********************************************************************
C  subroutine BLDSRCCELL body starts at line
C
C  DESCRIPTION:
C      This subroutine uses the gridding matrix to build a list of grid
C      cells and associated fractions for each source.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 4/10 by C. Seppanen
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: NSRCCELLS, SRCCELLS, SRCCELLFRACS

        IMPLICIT NONE

C.........  INCLUDES:

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NSRC        ! number of sources
        INTEGER, INTENT(IN) :: NGRID       ! number of grid cells
        INTEGER, INTENT(IN) :: NMAT        ! dimension for matrixes
        INTEGER, INTENT(IN) :: NX( NGRID ) ! number of sources per cell
        INTEGER, INTENT(IN) :: IX( NMAT )  ! list of sources per cell
        REAL,    INTENT(IN) :: CX( NMAT )  ! list of source fractions per cell

C.........  LOCAL VARIABLES and their descriptions:

C.........  Other local variables
        INTEGER   INDX, NG, NS       ! indexes and counters
        INTEGER   MXNSRCCELLS        ! max. no. cells per source
        INTEGER   CELL               ! cell number
        INTEGER   SRC                ! source number
        INTEGER   IOS                ! error status

        CHARACTER(16) :: PROGNAME = 'BLDSRCCELL' ! program name

C***********************************************************************
C   begin body of subroutine BLDSRCCELL

C.........  Determine maximum number of cells per source
        ALLOCATE( NSRCCELLS( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NSRCCELLS', PROGNAME )
        NSRCCELLS = 0   ! array

        MXNSRCCELLS = 0
        INDX = 0
        DO NG = 1, NGRID

C.............  Loop over sources for current grid cell
            DO NS = 1, NX( NG )
            
                INDX = INDX + 1
                SRC = IX( INDX )
                NSRCCELLS( SRC ) = NSRCCELLS( SRC ) + 1
                IF( NSRCCELLS( SRC ) > MXNSRCCELLS ) THEN
                    MXNSRCCELLS = NSRCCELLS( SRC )
                END IF
                
            END DO
        
        END DO

C.........  Store list of cells and fractions for each source        
        ALLOCATE( SRCCELLS( MXNSRCCELLS, NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCCELLS', PROGNAME )
        ALLOCATE( SRCCELLFRACS( MXNSRCCELLS, NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCCELLFRACS', PROGNAME )
        NSRCCELLS = 0       ! array
        SRCCELLS = 0        ! array
        SRCCELLFRACS = 0.0  ! array

        INDX = 0
        DO NG = 1, NGRID
        
            DO NS = 1, NX( NG )

                INDX = INDX + 1
                SRC = IX( INDX )
                NSRCCELLS( SRC ) = NSRCCELLS( SRC ) + 1
                
                CELL = NSRCCELLS( SRC )
                SRCCELLS    ( CELL, SRC ) = NG
                SRCCELLFRACS( CELL, SRC ) = CX( INDX )
            
            END DO
        
        END DO

        END SUBROUTINE BLDSRCCELL
