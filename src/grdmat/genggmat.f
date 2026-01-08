
        SUBROUTINE GENGGMAT( GNAME, FDEV, MXSCEL, NASRC, NMATX, 
     &                       NX, IX, CX, NCOEF, CMAX, CMIN )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by B. Baek 2/09
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C************************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: XLOCA, YLOCA, CIFIP, CELLID, CSOURC

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NCHARS, NSRC

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS 
C       CHARACTER(2)    CRLF
C       INTEGER         FIND1
C       INTEGER         FINDC
        LOGICAL         DSCM3GRD

C        EXTERNAL        CRLF, FIND1, FINDC, DSCM3GRD
        EXTERNAL     DSCM3GRD

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: GNAME         ! gridding mtx logical name
        INTEGER     , INTENT (IN) :: FDEV          ! surg codes report file
        INTEGER     , INTENT (IN) :: MXSCEL        ! max sources per cell
        INTEGER     , INTENT (IN) :: NASRC         ! no. mobile sources
        INTEGER     , INTENT (IN) :: NMATX         ! no. source-cell intersects
        INTEGER     , INTENT(OUT) :: NX  ( NGRID ) ! no. srcs per cell
        INTEGER     , INTENT(OUT) :: IX  ( NMATX ) ! src IDs 
        REAL        , INTENT(OUT) :: CX  ( NMATX ) ! gridding coefficients
        INTEGER     , INTENT(OUT) :: NCOEF         ! no. of gridding coeffs
        INTEGER     , INTENT(OUT) :: CMAX          ! max no. of sources per cell
        INTEGER     , INTENT(OUT) :: CMIN          ! min no. of sources per cell

C...........   Local allocatable arrays...

C...........   Scratch Gridding Matrix (subscripted by source-within-cell, cell)
        INTEGER, ALLOCATABLE :: IS ( :,: ) ! source IDs for each cell
        REAL   , ALLOCATABLE :: CS ( :,: ) ! factors 

C...........   Temporary array for flagging sources that are outside the
C              domain and for flagging sources with all zero surrogate fractions
        LOGICAL, ALLOCATABLE :: INDOMAIN( : ) ! true: src is in the domain

C...........   Temporary arrays for storing surrogate codes to use
        INTEGER, ALLOCATABLE :: SURGID1( : ) ! primary surrogate code
        INTEGER, ALLOCATABLE :: SURGID2( : ) ! secondary surrogate code

C...........   Other local variables
        INTEGER         C, F, I, J, K, N, S !  indices and counters.

        INTEGER         COL     ! tmp column
        INTEGER         ID1,ID2 ! tmp primary and secondary surg codes
        INTEGER         IOS     ! i/o status
        INTEGER         ISIDX   ! tmp surrogate ID code index
        INTEGER         JMAX    ! counter for storing correct max dimensions
        INTEGER         L2      ! string length
        INTEGER         NCEL    ! tmp number of cells 
        INTEGER         NNOSRG  ! no. of cy/st/co codes with no surrogates
        INTEGER         ROW     ! tmp row

        REAL            FRAC    ! tmp surrogate fraction

        LOGICAL      :: EFLAG = .FALSE.  !  true: error detected
        LOGICAL      :: LFLAG = .FALSE.  !  true: location data available
        LOGICAL      :: XYSET = .FALSE. ! true: X/Y available for src

        CHARACTER(FIPLEN3) CFIP   !  tmp country/state/county code
        CHARACTER(FIPLEN3) LFIP   !  cy/st/co code from previous iteration
        CHARACTER(16)   COORUNIT  !  coordinate system projection units
        CHARACTER(80)   GDESC     !  grid description
        CHARACTER(256)  MESG      !  message buffer 

        CHARACTER(SRCLEN3)    CSRC  ! tmp source chars string

        CHARACTER(16) :: PROGNAME = 'GENGGMAT' ! program name

C***********************************************************************
C   begin body of subroutine GENGGMAT

C.........  Initialize number of sources per cell counter
        NX = 0   ! Array

C.........  Allocate memory for temporary gridding matrix and other
        ALLOCATE( IS( MXSCEL, NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IS', PROGNAME )
        ALLOCATE( CS( MXSCEL, NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CS', PROGNAME )
        ALLOCATE( INDOMAIN( NASRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDOMAIN', PROGNAME )
        ALLOCATE( SURGID1( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SURGID1', PROGNAME )
        ALLOCATE( SURGID2( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SURGID2', PROGNAME )
        SURGID1 = 0   ! array
        SURGID2 = 0   ! array

C.........  Set flag to indicate that XLOCA/YLOCA are available
        LFLAG = ALLOCATED( XLOCA )

C.......   Compute gridding matrix:
C.......       First case:   explicit link (ILINK > 0
C.......       Second case:  some LNKDEF entry applies
C.......       Third case:   FIP+roadtype cross-reference match
C.......       Fourth case:  state+roadtype cross-reference match
C.......       fifth case:   roadtype cross-reference match
C.......       sixth case:   fallback default

        MESG = 'Computing gridding matrix and statistics...'
        CALL M3MSG2( MESG )

        LFIP  = ' '
        NNOSRG   = 0
        JMAX  = -1

        DO S = 1, NSRC
            
            CFIP = CIFIP ( S )
            C    = CELLID( S )
            CSRC = CSOURC( S )

C.............  Initialize sources as being in the domain
            INDOMAIN( S ) = .TRUE.

C.............  Special case for source already assigned a grid cell
            IF( C .GT. 0 ) THEN
            
                J = NX( C )

C.................  Check that the maximum number of sources per cell is ok                
                IF( J .LT. MXSCEL ) THEN
                    J = J + 1
                    IS( 1,C ) = S
                    CS( 1,C ) = 1.0
                    
                ELSE
                    IF( J+1 .GT. JMAX ) JMAX = J+1
                END IF
                
                NX( C ) = J
                
                CYCLE           ! To head of loop over sources
            END IF

        END DO        !  end loop on sources S, computing gridding matrix.

C.........  Abort if overflow occurred
        IF ( JMAX .GT. MXSCEL ) THEN   
            
            WRITE( MESG,94010 )
     &       'INTERNAL ERROR: Gridding matrix not ' //
     &       'written.' // CRLF() // BLANK10 //
     &       'Arrays would have overflowed.' 
     &       // CRLF() // BLANK10 // 
     &       'Current maximum sources per cell (MXSCEL) =', MXSCEL, '.'
     &       // CRLF() // BLANK10 // 
     &       'Actual  maximum sources per cell          =', JMAX  , '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

C.........  Compress matrix into I/O representation from scratch 
C.........  representation and compute statistics.
                   
        K    = 0
        CMAX = 0
        CMIN = 99999999
        DO C = 1, NGRID  ! Loop through cells
            
            J = NX( C )
                   
            IF (      J .GT. CMAX ) THEN
                CMAX = J
            ELSE IF ( J .GT. 0 .AND. J .LT. CMIN ) THEN
                CMIN = J
            END IF
                   
            DO N = 1, J  ! Loop through sources in this cell
                K = K + 1
                IF ( K .LE. NMATX ) THEN
                    S       = IS( N,C )
                    IX( K ) = S
                    CX( K ) = CS( N,C )
                END IF
            END DO
                   
        END DO    !  end of loop on cells C for this FIP

        NCOEF = K

C.........  Write gridding matrix
        MESG = 'Writing out GRIDDING MATRIX file...'
        CALL M3MSG2( MESG )

        IF( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
            MESG = 'Error writing GRIDDING MATRIX file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write output surrogates codes
        DO S = 1, NSRC
            WRITE( FDEV,93360 ) S, SURGID1( S ), SURGID2( S )
        END DO

C.........  Dellallocate locally allocated memory
        DEALLOCATE( IS, CS, INDOMAIN, SURGID1, SURGID2 )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93360   FORMAT( I8, 1X, I8, 1X, I8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GENGGMAT

