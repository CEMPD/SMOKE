
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  TOTV( VNAME, N, VEC )

C***********************************************************************
C  subroutine body starts at line  41
C
C  FUNCTION:
C
C
C  PRECONDITIONS REQUIRED:
C
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        CHARACTER*16	VNAME
        INTEGER		N
        REAL		VEC( N )


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER		S
        REAL		TOT


C***********************************************************************
C   begin body of subroutine  TOTV

        TOT = 0.0
        DO 111 S = 1, N
            TOT = TOT + VEC( S )
111     CONTINUE

        WRITE( *,92010 ) VNAME, TOT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92010	FORMAT( 5X, 'Species "', A, '" total = ', 1PG12.6 )

        END

