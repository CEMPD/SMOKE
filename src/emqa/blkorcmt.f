
        LOGICAL FUNCTION BLKORCMT( LINE )

C***********************************************************************
C  function body starts at line 49
C
C  DESCRIPTION:
C    The BLKORCMT routine will return TRUE if line is blank or a comment
C    line.
C
C  PRECONDITIONS REQUIRED:
C    LINE is defined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***********************************************************************

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: LINE   ! line of data

C***********************************************************************
C   begin body of function BLKORCMT

        BLKORCMT = .FALSE.

        IF( LINE        .EQ. ' ' .OR.
     &      LINE( 1:1 ) .EQ. '#'      ) BLKORCMT = .TRUE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END FUNCTION BLKORCMT
        
 
