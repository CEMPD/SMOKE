
        SUBROUTINE WRPTREF( FDEV, NPSRC, IDIU, IWEK, IMON )

C***********************************************************************
C  subroutine body starts at line 70
C
C  DESCRIPTION:
C      This subroutine writes a temporal x-ref file in SMOKE format. It is
C      only used for the EMS-95 input files.
C
C  PRECONDITIONS REQUIRED:
C      File opened with unit number FDEV
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
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  SUBROUTINE ARGUMENTS
        INTEGER      FDEV           !  ASCII file unit
        INTEGER      NPSRC          !  actual source count
        INTEGER      IDIU( NPSRC )  !  source FIPS (county) ID
        INTEGER      IWEK( NPSRC )  !  source SCC
        INTEGER      IMON( NPSRC )  !  source SIC

C...........   Other local variables
        INTEGER         S

        CHARACTER*300   MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'WRPTREF' !  program name

C***********************************************************************
C   begin body of program WRPTREF

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

        END
