
        SUBROUTINE WRCHRSCC( CSCC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes character string SCC codes
C
C  PRECONDITIONS REQUIRED:
C      Index sorted in order of increasing ISCC value and ISCC populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
C
C************************************************************************
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

C.........  MODULES for public variables
C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVSCC, INVSCC 

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CRL, NSRC

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     PROMPTFFILE
        EXTERNAL    PROMPTFFILE

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CSCC( NSRC )  !  unsorted SCCs

C...........   Other local variables
        INTEGER       J, S
        INTEGER       FDEV                !  output file unit number

        CHARACTER(300) MESG                !  message buffer
        CHARACTER(SCCLEN3) TSCC, LSCC !  current and previous 10-digit SCC

        CHARACTER(16) :: PROGNAME = 'WRCHRSCC' !  program name

C***********************************************************************
C   begin body of subroutine WRCHRSCC

C.............  Get unique-SCC output file
        FDEV = PROMPTFFILE( 
     &          'Enter the name of the ACTUAL SCC output file',
     &          .FALSE., .TRUE., CRL // 'SCC', PROGNAME )

C.........  Write out SCCs list (this subroutine expect the data structure
C           that is being provided.

        CALL M3MSG2( 'Writing out ACTUAL SCC file...' )

C.........  Call the routine that generates the unique SCC list 
        CALL GENUSLST

C.........  Write unique SCCs list  
        DO J = 1, NINVSCC

            WRITE( FDEV, '(A)', ERR=6001 ) INVSCC( J )

        END DO 

        RETURN

6001    MESG = 'Problem writing SCC file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE WRCHRSCC
