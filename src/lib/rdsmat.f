
        SUBROUTINE RDSMAT( FNAME, VDESC, SMAT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads a speciation matrix for any source category.
C      The variable description is provided as the name for reading, and
C      this routine finds which variable name to read and reads in the data.
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
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         INDEX1

        EXTERNAL        CRLF, INDEX1

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! speciation matrix file name
        CHARACTER(*), INTENT (IN) :: VDESC     ! variable description to read
        REAL        , INTENT(OUT) :: SMAT( * ) ! coeffs for sources

C.........  Other local variables
        INTEGER         J, L, L1, L2, V       !  counters and indices

        CHARACTER(LEN=IODLEN3 ) DBUF    !  variable description buffer
        CHARACTER(LEN=IOVLEN3 ) VBUF    !  variable name buffer
        CHARACTER*300           MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDSMAT' ! program name

C***********************************************************************
C   begin body of subroutine RDSMAT

C.........  Retrieve file header
        IF ( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description of file ' // FNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Make sure variable name has proper spacing that is the same as
C           file header
        L1 = INDEX( VDESC, '_' )
        L2 = LEN_TRIM( VDESC )
        DBUF = VDESC( 1:L1-1 )
        DBUF = DBUF( 1:IOVLEN3 ) // '_' // VDESC( L1+1:L2 )

C.........  Find variable description in list of descriptions
        J = INDEX1( DBUF, NVARS3D, VDESC3D )

        IF( J .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Speciation variable description "'//
     &             DBUF // '" not found in speciation matrix.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Read variable and print nice error message if cannot
        VBUF = VNAME3D( J )
        IF ( .NOT. READ3( FNAME, VBUF, 1, 0, 0, SMAT ) ) THEN

            L  = LEN_TRIM( FNAME )
            L1 = LEN_TRIM( VBUF )
            MESG = 'Could not read speciation matrix from file "' //
     &             FNAME( 1:L ) // '"' // CRLF() // BLANK10 //
     &             'for variable "' // VBUF( 1:L1 ) //
     &             '" with description "' // VDESC( 1:L2 ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF    !  if read3() failed for speciation matrix

        RETURN

        END SUBROUTINE RDSMAT
