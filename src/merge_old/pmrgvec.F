
        SUBROUTINE  PMRGVEC( NCOLP, P, NP, IP, NROWS, V )

C***********************************************************************
C  subroutine body starts at line  62
C
C  DESCRIPTION:  
C       multiply vectors A and P by sparse incidence matrix
C	<NP,IP> and add the result to V for CSGMERGE program.
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C	Prototype 2/96 by CJC
c
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
        
        INTEGER         NCOLP           ! length of input vector P
        INTEGER         NROWS           ! length of output vector
        REAL            P ( NCOLP )     !  input vector P
        INTEGER         NP( NROWS )     ! # of entries per row
        INTEGER         IP( NCOLP )     ! columns list
        REAL            V( NROWS )      ! output vector


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         R, C, KP
        REAL            SUM


C***********************************************************************
C   begin body of subroutine  PMRGVEC

        KP = 0

        DO  33  R = 1, NROWS
            
            SUM = V( R )
            
            DO  22  C = 1, NP( R )
                KP = KP + 1
                SUM = SUM  +  P( IP( KP ) )
22          CONTINUE

            V( R ) = SUM
        
33      CONTINUE

        RETURN
        END

