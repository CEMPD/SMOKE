C Version @(#)$Id$ $Source$ $Date$ 

C..  Created by M. Houyoux on 12/16.  Although LNBLNK is referenced many
C..  times throughout MOBILE5, I cannot find any such function.  I believe the 
C..  purpose of the function is to give the last non-blank field number in
C..  a left-justified character variable.  As such, I hae copied TRIMLEN from
C..  the Model3 stuff and modified it to serve at LNBLNK.

        INTEGER  FUNCTION  LNBLNK ( STRING )

C***********************************************************************
C  function body starts at line 36
C
C  FUNCTION:  return the effective length of argument CHARACTER*(*) STRING,
C             after trailing blanks have been trimmed.
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:  
C             Prototype 8/91 by CJC
C             Version 2/93 for CRAY by CJC
C
C***********************************************************************

      IMPLICIT NONE


C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*)	STRING


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER  	L, K


C***********************************************************************
C   begin body of function  LNBLNK

        L = LEN( STRING )
        DO  11  K = L, 1, -1
            IF ( STRING( K:K ) .NE. ' ' ) THEN
                GO TO  12
            END IF
11      CONTINUE

        K = 1

12      CONTINUE

        LNBLNK = K

        RETURN

        END

