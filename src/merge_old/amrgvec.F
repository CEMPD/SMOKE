
        SUBROUTINE  AMRGVEC( NCOLA, A, NCOFF, NA, IA, MA, NROWS, V )

C***********************************************************************
C  subroutine body starts at line  63
C
C  FUNCTION:  multiply vectors A and P by sparse matrices <NA,IA,CA> 
C	and add the result to V, for CSGMERGE program.
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 2/96 by CJC
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         NCOLA           ! length of input vector A
        INTEGER         NCOFF           ! max number of A-coefficients
        INTEGER         NROWS           ! length of output vector
        REAL            A ( NCOLA )     !  input vector A
        INTEGER         NA( NROWS )     ! # of entries per row
        INTEGER         IA( NCOFF )     ! columns list
        REAL            MA( NCOFF )     ! coefficient array
        REAL            V ( NROWS )     ! output vector


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         R, C, KA
        REAL            SUM


C***********************************************************************
C   begin body of subroutine  AMRGVEC

        KA = 0

        DO  33  R = 1, NROWS
            
            SUM = V( R )
            
            DO  11  C = 1, NA( R )
                KA = KA + 1
                SUM = SUM  +  MA( KA ) * A( IA( KA ) )
11          CONTINUE

            V( R ) = SUM
        
33      CONTINUE

        RETURN
        END

