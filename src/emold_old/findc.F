
C Version "@(#)$Id$ $Source$ $Date$ 

        INTEGER FUNCTION FINDC( KEY, N, LIST )

C***********************************************************************
C  function body starts at line 44
C
C  RETURNS:
C       subscript at which the targeted key appears, or 
C       -1 in case of failure.
C       EXAMPLE:  search for FIP in table of FIP values
C
C  PRECONDITIONS REQUIRED:
C       Sorted table <N,LIST> for searching
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 8/1/96 by CJC
C
C***********************************************************************

      IMPLICIT NONE


C...........   ARGUMENTS and their descriptions:
        
        CHARACTER*(*)  KEY           !  key
        INTEGER        N             !  table size
        CHARACTER*(*)  LIST( N )     !  table to search for KEY


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER  LO
        INTEGER  HI
        INTEGER  M


C***********************************************************************
C   begin body of function  FINDC

        LO = 1
        HI = N
        
11      CONTINUE
            
            IF ( LO .GT. HI ) THEN
            
                FINDC = -1
                RETURN
                
            END IF
           
            M = ( LO + HI ) / 2
            IF ( KEY .GT. LIST( M ) ) THEN
                LO = M + 1
                GO TO  11
            ELSE IF ( KEY .LT. LIST( M ) ) THEN
                HI = M - 1
                GO TO  11
            END IF          !  end of bin search loop for this K
        
        
        FINDC = M
        RETURN
        END

