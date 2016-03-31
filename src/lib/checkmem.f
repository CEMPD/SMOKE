
        SUBROUTINE CHECKMEM( MSTATUS, ONVAR, CALLER )
 
C***********************************************************************
C  subroutine body starts at line  105
C
C  DESCRIPTION:
C       Reports an error and exits if memory status flag is non-zero.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C***********************************************************************
C
C Project Title: EDSS Tools Library
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
C***************************************************************************
 
      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

       INTEGER     , INTENT( IN ) ::  MSTATUS !  ALLOCATE function exit status
       CHARACTER(*), INTENT( IN ) ::  ONVAR   !  Variable name(s) of previous ALLOCATE statement
       CHARACTER(*), INTENT( IN ) ::  CALLER  !  Name of calling program

C...........   Local variables

       CHARACTER(256)  MESG
       CHARACTER(32)   PNAME

C***********************************************************************
C   begin body of function CHECKMEM

C.........  Abort if memory status is non-zero

        IF( MSTATUS .GT. 0 ) THEN            
            PNAME = TRIM( CALLER ) // ':' // 'CHECKMEM'
            WRITE( MESG, '( 3 A, I10 )' )
     &             'Failure allocating memory for "', TRIM( ONVAR ),
     &             '":  STATUS=', MSTATUS
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

        END SUBROUTINE CHECKMEM

