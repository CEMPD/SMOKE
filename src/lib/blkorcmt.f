
        LOGICAL FUNCTION BLKORCMT( LINE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The BLKORCMT routine will return true of the buffer given it is a blank
C    or a comment.
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

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: LINE        ! line of data

C...........   Local variables
        INTEGER      J

C***********************************************************************
C   begin body of subroutine BLKORCMT


        IF( LINE .EQ. ' ' .OR. LINE( 1:1 ) .EQ. CINVHDR ) THEN
            BLKORCMT = .TRUE.

        ELSE
            BLKORCMT = .FALSE.

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END FUNCTION BLKORCMT
        
 
