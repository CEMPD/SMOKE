
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  SMATVEC( NCOLS, NROWS, NCOFF, N, I, M, U, V )

C***********************************************************************
C  subroutine body starts at line  43
C
C  FUNCTION:  multiply a sparse matrix <N,I,C> by a vector U and 
C             return the result V
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 2/95 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         NCOLS           ! length of input vector
        INTEGER         NROWS           ! length of output vector
        INTEGER         NCOFF           ! max number of coefficients
        
        INTEGER         N( NROWS )      ! # of entries per row
        INTEGER         I( NCOFF )      ! columns list
        REAL            M( NCOFF )      ! coefficient array
        
        REAL            U( NCOLS )      !  input vector
        REAL            V( NROWS )      ! output vector


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         R, C, K
        REAL            SUM


C***********************************************************************
C   begin body of subroutine  SMATVEC

        K = 0
        DO  22  R = 1, NROWS
            
            SUM = 0.0
            
            DO  11  C = 1, N( R )
                K = K + 1
                SUM = SUM  +  M( K ) * U( I( K ) )
11          CONTINUE
            
            V( R ) = SUM
        
22      CONTINUE

      RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

91000   FORMAT ( //5X , '*** ERROR ABORT in subroutine SMATVEC ***',
     &            /5X , A ,
     &           // )        !  generic error message format


C...........   Informational (LOG) message formats... 92xxx


C...........   Formatted file I/O formats............ 93xxx


C...........   Internal buffering formats............ 94xxx


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

