
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  MSMATVEC( N, M, A, V, C )

C***********************************************************************
C  subroutine body starts at line  43
C
C  FUNCTION:  apply a layered diagonal matrix to a layered vector
C             and sum to existing output vector
C
C  PRECONDITIONS REQUIRED: Initialization of C input vector
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 2/95 by CJC
C       modified  6/95 by M Houyoux to add to exsiting C
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         N               ! length of input vector
        INTEGER         M               ! number of "layers" to aggregate
        REAL            A( N, M )       ! diagonal coeff matrix
        REAL            V( N, M )       ! input  vector
        REAL            C( N    )       ! output vector


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         L, R
        REAL            S


C***********************************************************************
C   begin body of subroutine  MSMATVEC

        DO  22  R = 1, N
            
            S = 0.0
            DO  11  L = 1, M
                S = S  +  A( R,L ) * V( R,L )
11          CONTINUE
            C( R ) = C( R ) + S
        
22      CONTINUE

        RETURN
        END

