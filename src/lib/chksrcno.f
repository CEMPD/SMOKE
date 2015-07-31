
        SUBROUTINE CHKSRCNO( CATDESC, FILNAM, NTEST, NSRC, EFLAG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine compares the number of rows in an I/O API header 
C      (NROWS3D) with the input number of sources for the purpose of checking
C      the number of sources between two files
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
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

        IMPLICIT NONE

C.........  INCLUDE FILES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  EXTERNAL FUNCTIONS
        CHARACTER(2)  CRLF
        EXTERNAL      CRLF

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*)    CATDESC    ! source category description
        CHARACTER(*)    FILNAM     ! name of file being checked
        INTEGER         NTEST      ! number of sources to check
        INTEGER         NSRC       ! number of sources to compare against
        LOGICAL         EFLAG      ! error flag

C.........  Other local variables

        CHARACTER(300)  MESG 

        CHARACTER(16) :: PROGNAME = 'CHKSRCNO' ! program name

C***********************************************************************
C   begin body of subroutine CHKSRCNO

C.........  Check the number of sources
        IF( NTEST .NE. NSRC ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'ERROR: Dimension mismatch. Source number in ' //
     &             FILNAM // ' file: ', NTEST, ', ' // CRLF() //
     &             BLANK10 // 'but in ' // CATDESC // 
     &             ' source inventory file: ', NSRC
            CALL M3MSG2( MESG )

        END IF 

        RETURN
C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE CHKSRCNO

