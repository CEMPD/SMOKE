
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  PCOEF( N, X, Y, C )

C***********************************************************************
C  subroutine body starts at line
C
C  FUNCTION:  returns array C of polynomial coefficients
C	for the polynomial fitting the data 
C       Y( K ) = P( X( K ) ) = sum_I
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C	Prototype 11/95 by CJC adapted from _Numerical_Recipes_,
C	section 3.5 (pp. 92-94), routine POLCOE
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        INTEGER		N
        REAL		X( N )
        REAL		Y( N )
        REAL		C( N )


C...........   PARAMETERS and their descriptions:

        INTEGER		NMAX
        PARAMETER     ( NMAX = 15 )


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER		I, J
        REAL		P, F, B
        REAL		S( NMAX )
        REAL		T( NMAX )
        CHARACTER*80    MESG

C***********************************************************************
C   begin body of subroutine  PCOEF

        IF ( N .GT. NMAX ) THEN
            WRITE( MESG, '(A, I5, 2X, A )' )
     &          'Requested degree', N, 'exceeds max=15'
            CALL M3EXIT( 'PCOEF', 0, 0, MESG, 2 )
        END IF

        DO  11  I = 1, N

            S( I ) = 0.0
            C( I ) = 0.0
            T( I ) = FLOAT( I )

11      CONTINUE
        
        S( N ) = -X( 1 )
        
        DO  22  I = 2, N
        DO  21  J = N - I + 1, N - 1

            S( J ) = S( J ) - X( I ) * S( J+1 )

21      CONTINUE
22      CONTINUE

        DO  33  I = 1, N

            P = T( N )

            DO  31  J = N - 1, 1, -1

                P = T( J ) * S( J+1 ) + P * X( I )

31          CONTINUE		!  end loop on J

            F = Y( I ) / P
            B = 1.0

            DO  32  J = N, 1, -1

                C( J ) = C( J ) + B * F
                B      = S( J ) + B * X( I )

32          CONTINUE	!  end loop on J

33      CONTINUE	!  end loop on I

        RETURN

        END

