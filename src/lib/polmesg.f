
        SUBROUTINE POLMESG( NPOL, NAMES )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes out a message stating that the pollutants in the
C      argument list are being processed.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C
C  REVISION  HISTORY:
C      Created 3/99 by M. Houyoux
C
C************************************************************************
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

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NPOL           !  number of pollutants
        CHARACTER(*), INTENT (IN) :: NAMES( NPOL )  !  pollutant names

C...........   Other local variables
        INTEGER       I, L1, L2

        CHARACTER*300 MESG           !  message buffer

        CHARACTER*16 :: PROGNAME = 'POLMESG' !  program name

C***********************************************************************
C   begin body of subroutine POLMESG

        MESG = 'Processing pollutants: '

        DO I = 1, NPOL

            L1 = LEN_TRIM( MESG )
            IF( NAMES( I ) .EQ. ' ' ) EXIT

            L2 = LEN_TRIM( NAMES( I ) )
            MESG = MESG( 1:L1 ) // ', "' // NAMES( I )( 1:L2 ) // '"'

        ENDDO

        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE POLMESG
