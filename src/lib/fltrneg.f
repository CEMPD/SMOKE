
        SUBROUTINE FLTRNEG( STRING )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine filters the valid "missing" terms that can be in a 
C      cross-reference file changes the argument value to blank when it is
C      one of the missing values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     
C
C  REVISION  HISTORY:
C      Started 6/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN OUT) :: STRING   ! string for filter '-9'

C***********************************************************************
C   Begin body of subroutine FLTRNEG

        IF( INDEX( STRING, '-9' ) .GT.  0  .OR.
     &      STRING                .EQ. '0'      ) STRING = ' '

        RETURN

        END SUBROUTINE FLTRNEG
