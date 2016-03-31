
        SUBROUTINE PADNZERO( N, STRING )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine right-justifies a string and replaces the leading spaces
C      with zeros for a string of size N.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 3/02 by G. Cano
C
C****************************************************************************/
C
C Project Title: EDSS Tools Library
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
C Pathname: @(#)$Id$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT(IN)     :: N      ! length of the string to pad
        CHARACTER(N), INTENT(IN OUT) :: STRING ! character string to adjust

        INTEGER         I, L

        CHARACTER(16) :: PROGNAME = 'PADNZERO' ! program name

C***********************************************************************
C   begin body of subroutine PADZERO

        STRING = ADJUSTR( STRING )

        DO I = 1, N

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
