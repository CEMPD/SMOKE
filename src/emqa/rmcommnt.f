
        SUBROUTINE RMCOMMNT( CMT_DELIM, LINE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The RMCOMMNT routine will delete the comment part of a line, if it is
C    there.
C
C  PRECONDITIONS REQUIRED:
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
        CHARACTER(*), INTENT (IN)     :: CMT_DELIM   ! comment delimeter
        CHARACTER(*), INTENT (IN OUT) :: LINE        ! line of data

C...........   Local variables
        INTEGER      J

C***********************************************************************
C   begin body of subroutine RMCOMMNT

        J = INDEX( LINE, CMT_DELIM )

C.........  If comment delimeter is found, then remove remainder of line
        IF( J .GT. 0 ) THEN

            J = J + LEN_TRIM( CMT_DELIM ) - 1

            LINE = LINE( 1:J )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RMCOMMNT
        
 
