
        SUBROUTINE  AGRDSUM( NCOLS, NROWS, NCOFF, N, I, M, 
     &                       PTR, U, V )

C***********************************************************************
C  subroutine body starts at line  67
C
C  DESCRIPTION:  
C       multiply a sparse matrix <N,I,C> by a vector U and then process the 
C       result based on a sorted key to return V with a value for each 
C       different entry in the key
C
C  PRECONDITIONS REQUIRED:  
C       argument 'V' must be initialized 
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 4/97 by M Houyoux
C
C***********************************************************************
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
C****************************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         NCOLS           ! length of input vector
        INTEGER         NROWS           ! length of output vector
        INTEGER         NCOFF           ! max number of coefficients
        
        INTEGER         N( NROWS )      ! # of entries per row of sparse matrix
        INTEGER         I( NCOFF )      ! columns list of sparse matrix
        REAL            M( NCOFF )      ! coefficient array of sparse matrix
     
        INTEGER         PTR ( NCOLS )   !  summation pointer
        REAL            U   ( NCOLS )   !  input vector
        REAL            V   ( *     )   !  output vector

C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         R, C, K, L, P

C***********************************************************************
C   begin body of subroutine  AGRDSUM

        K = 0
        DO  22  R = 1, NROWS

            DO  11  C = 1, N( R )
                K = K + 1
                L =   I( K )
                P = PTR( L )
                V( P ) = V( P ) +  M( K ) * U( L )
11          CONTINUE
            
22      CONTINUE

      RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END

