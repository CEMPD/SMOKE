
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
C       Adapted 10/98 by M Houyoux
C
C***********************************************************************
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
 
      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

       INTEGER       MSTATUS !  ALLOCATE function exit status
       CHARACTER*(*) ONVAR   !  Variable name of previous ALLOCATE statement
       CHARACTER*(*) CALLER  !  Name of calling program

C...........   ARGUMENTS and their descriptions:
       INTEGER      TRIMLEN
       EXTERNAL     TRIMLEN

C...........   Local variables

       INTEGER         L1
       INTEGER         L2
       CHARACTER*256   MESG

       CHARACTER*16 :: PROGNAME = 'CHECKMEM' ! program name

C***********************************************************************
C   begin body of function CHECKMEM

C.........  Get lengths of input character strings
        L1 = TRIMLEN( ONVAR )
        L2 = TRIMLEN( CALLER )

C.........  Abort if memory status is non-zero

        IF( MSTATUS .GT. 0 ) THEN            
            MESG = 'Failure allocating memory for "' // ONVAR( 1:L1 ) //
     &             '" variable'
            CALL M3EXIT( CALLER( 1:L2 ), 0, 0, MESG, 2 )
        ENDIF

        RETURN

        END

