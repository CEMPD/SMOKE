
        SUBROUTINE RDUMAT( FNAME, NSRC, NMAT1, NMAT2, NU, IU, CU )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads an ungridding matrix for any source category
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

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME       ! ungridding matrix name
        INTEGER     , INTENT (IN) :: NSRC        ! number of sources
        INTEGER     , INTENT (IN) :: NMAT1       ! dim 1 for matrix
        INTEGER     , INTENT (IN) :: NMAT2       ! dim 2 for matrix
        INTEGER     , INTENT(OUT) :: NU( NSRC  ) ! number of cells per source
        INTEGER                   :: IU( NMAT1 ) ! list of cells per source
        REAL                      :: CU( NMAT2 ) ! coefficients for cells

C.........  Other local variables
        CHARACTER(300)  MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDUMAT' ! program name

C***********************************************************************
C   begin body of subroutine RDUMAT

        MESG = 'Reading ungridding matrix...'
        CALL M3MSG2( MESG )

        IF ( .NOT. READ3( FNAME, 'ALL', 1, 0, 0, NU ) ) THEN

            MESG = 'Could not read ungridding matrix from file "' //
     &             FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF      !  if read3() failed for ungridding matrix

        RETURN

        END SUBROUTINE RDUMAT
