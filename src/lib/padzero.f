
        SUBROUTINE PADZERO( STRING )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine right-justifies a string and replaces the leading spaces
C      with zeros.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
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
C Pathname: @(#)$Id$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*), INTENT(IN OUT) :: STRING ! character string to adjust

        INTEGER         I, L

        CHARACTER*16 :: PROGNAME = 'PADZERO' ! program name

C***********************************************************************
C   begin body of subroutine PADZERO

        STRING = ADJUSTR( STRING )

        L = LEN( STRING )

        DO I = 1, L

            IF( STRING( I:I ) .EQ. ' ' ) THEN

                STRING( I:I ) = '0'

            ELSE

                EXIT

            ENDIF

        ENDDO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )   

        END
