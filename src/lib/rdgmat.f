
        SUBROUTINE RDGMAT( FNAME, NGRID, NMAT1, NMAT2, NX, IX, CX )

C***************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads a gridding matrix for any source category
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  External functions
C       CHARACTER(2) CRLF
C        EXTERNAL     CRLF

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME       ! gridding matrix name
        INTEGER     , INTENT (IN) :: NGRID       ! number of grid cells
        INTEGER     , INTENT (IN) :: NMAT1       ! dim 1 for matrix
        INTEGER     , INTENT (IN) :: NMAT2       ! dim 2 for matrix
        INTEGER     , INTENT(OUT) :: NX( NGRID ) ! number of sources per cell
        INTEGER                   :: IX( NMAT1 ) ! list of sources per cell
        REAL                      :: CX( NMAT2 ) ! coefficients for sources

C.........  Other local variables
        INTEGER         C       !  tmp cell number
        INTEGER         NSUM    !  count of gridding matrix size

        CHARACTER(300)  MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDGMAT' ! program name

C***********************************************************************
C   begin body of subroutine RDGMAT

C.........  Read matrix
        IF ( .NOT. READ3( FNAME, 'ALL', 1, 0, 0, NX ) ) THEN

            MESG = 'Could not read gridding matrix from file "' //
     &             FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF      !  if read3() failed for gridding matrix

C.........  Check to make sure that data are consistent with header
        NSUM = 0
        DO C = 1, NGRID
            NSUM = NSUM + NX( C )
        ENDDO

        IF( NSUM .GT. NMAT1 ) THEN

            MESG = 'ERROR: Gridding matrix dimension is inconsistent '//
     &             'with records count!' // CRLF() // '          ' //
     &             'Delete gridding matrix and recreate it.'

            CALL M3MSG2( MESG )

            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

        RETURN

        END SUBROUTINE RDGMAT
