
        SUBROUTINE RDGMAT( FNAME, NGRID, NMAT1, NMAT2, NX, IX, CX )

C***********************************************************************
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
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME       ! gridding matrix name
        INTEGER     , INTENT (IN) :: NGRID       ! number of grid cells
        INTEGER     , INTENT (IN) :: NMAT1       ! dim 1 for matrix
        INTEGER     , INTENT (IN) :: NMAT2       ! dim 2 for matrix
        INTEGER     , INTENT(OUT) :: NX( NGRID ) ! number of sources per cell
        INTEGER     , INTENT(OUT) :: IX( NMAT1 ) ! list of sources per cell
        REAL        , INTENT(OUT) :: CX( NMAT2 ) ! coefficients for sources

C.........  Other local variables
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDGMAT' ! program name

C***********************************************************************
C   begin body of subroutine RDGMAT

        IF ( .NOT. READ3( FNAME, 'ALL', 1, 0, 0, NX ) ) THEN

            MESG = 'Could not read gridding matrix from file "' //
     &             FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF      !  if read3() failed for gridding matrix

        RETURN

        END SUBROUTINE RDGMAT
