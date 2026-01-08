
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
C       Created ??/???? by ???
C       09/2025 by HT UNC-IE:  Use M3UTILIO
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
C***************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
c       CHARACTER(2)    CRLF
c       INTEGER         INDEX1

c       EXTERNAL        CRLF, INDEX1

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! speciation matrix file name
        CHARACTER(*), INTENT (IN) :: VDESC     ! variable description to read
        REAL        , INTENT(OUT) :: SMAT( * ) ! coeffs for sources

C.........  Other local variables
        INTEGER         J, L, L1, L2, V       !  counters and indices

        CHARACTER(IODLEN3) DBUF    !  variable description buffer
        CHARACTER(IOVLEN3) VBUF    !  variable name buffer
        CHARACTER(300)     MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDSMAT' ! program name

C***********************************************************************
C   begin body of subroutine RDSMAT

C.........  Retrieve file header
        IF ( .NOT. DESCSET( FNAME, ALLFILES ) ) THEN
            MESG = 'Could not get description of file ' // FNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Make sure variable name has proper spacing that is the same as
C           file header
        L1 = INDEX( VDESC, SPJOIN )
        L2 = LEN_TRIM( VDESC )
        DBUF = VDESC( 1:L1-1 )
        DBUF = DBUF( 1:IOVLEN3 ) // SPJOIN // VDESC( L1+1:L2 )

C.........  Find variable description in list of descriptions
        J = INDEX1( DBUF, NVARSET, VDESCSET )

        IF( J .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Speciation variable description "'//
     &             DBUF // '" not found in speciation matrix.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Read variable and print nice error message if cannot
        VBUF = VNAMESET( J )
        IF ( .NOT. READSET( FNAME, VBUF, 1, ALLFILES, 
     &                      0, 0, SMAT ) ) THEN

            L  = LEN_TRIM( FNAME )
            L1 = LEN_TRIM( VBUF )
            MESG = 'Could not read speciation matrix from file "' //
     &             FNAME( 1:L ) // '"' // CRLF() // BLANK10 //
     &             'for variable "' // VBUF( 1:L1 ) //
     &             '" with description "' // VDESC( 1:L2 ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF    !  if readset() failed for speciation matrix

        RETURN

        END SUBROUTINE RDSMAT
