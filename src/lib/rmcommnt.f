
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
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
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
C           MRH update 1/2/09: remove the comment field too and allow for
C           length of comment delimiter greater than 1 character (e.g.,
C           Smkreport's ## in-line comment delimiter.
        IF( J .GT. 0 ) THEN

            J = J - 1
            LINE = LINE( 1:J )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RMCOMMNT
        
 
