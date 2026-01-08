
        SUBROUTINE WRPTREF( NPSRC, IDIU, IWEK, IMON )

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C      This subroutine writes a temporal x-ref file in SMOKE format. It is
C      only used for the EMS-95 input files.
C
C  PRECONDITIONS REQUIRED:
C      Input arrays populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
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

C.........  MODULES for public variables
        USE M3UTILIO

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CRL

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C       INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS and their descriptions:
C        INTEGER     PROMPTFFILE
C        EXTERNAL    PROMPTFFILE

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NPSRC          !  actual source count
        INTEGER, INTENT (IN) :: IDIU( NPSRC )  !  source FIPS (county) ID
        INTEGER, INTENT (IN) :: IWEK( NPSRC )  !  source SCC
        INTEGER, INTENT (IN) :: IMON( NPSRC )  !  source SIC

C...........   Other local variables
        INTEGER         S

        INTEGER         FDEV             !  output file unit number

        CHARACTER(300)  MESG             !  message buffer

        CHARACTER(16) :: PROGNAME = 'WRPTREF' !  program name

C***********************************************************************
C   begin body of program WRPTREF

        MESG = 'Enter the name of the TEMPORAL X-REF output file'
        FDEV = PROMPTFFILE( MESG, .FALSE., .TRUE.,
     &                      CRL // 'TREF_ALT', PROGNAME )

        MESG = 'Writing out TEMPORAL CROSS-REFERENCE file...'
        CALL M3MSG2( MESG )

        DO S = 1, NPSRC

            WRITE( FDEV, 93040, ERR=6001 ) 1, IWEK( S ), IDIU( S )

        ENDDO 

        RETURN

6001    MESG = 'ERROR writing temporal x-ref file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93040   FORMAT( I5, ',' ,I5, ',' ,I5 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE WRPTREF
