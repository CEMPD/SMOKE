
C Version "@(#)$Id$ $Source$ $Date$ 

        INTEGER FUNCTION LOCAT1( K1, N, LIST1 )

C***********************************************************************
C  function body starts at line 48
C
C  RETURNS:
C       subscript at which the targeted key should be inserted
C       or -1 if the key is found.
C       EXAMPLE:  search for <FIP> in table of FIP values.
C
C  PRECONDITIONS REQUIRED:
C       Sorted table <N,LIST1> for searching
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 12/95 by MRH copied from FIND1 of CJC
C
C***********************************************************************

      IMPLICIT NONE


C...........   ARGUMENTS and their descriptions:
        
        INTEGER  K1             !  first  key
        INTEGER  N              !  table size
        INTEGER  LIST1( N )     !  table to search for K1


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER  LO
        INTEGER  HI
        INTEGER  M


C***********************************************************************
C   begin body of function  LOCAT1

        LO = 1
        HI = N

        IF( N .EQ. 0 ) THEN

            LOCAT1 = -1
            RETURN

        ENDIF 

11      CONTINUE
           
            IF ( LO .GT. HI ) THEN
            
                LOCAT1 = LO
                RETURN
                
            END IF
           
            M = ( LO + HI ) / 2
            IF ( K1 .GT. LIST1( M ) ) THEN
                LO = M + 1
                GO TO  11
            ELSE IF ( K1 .LT. LIST1( M ) ) THEN
                HI = M - 1
                GO TO  11
            END IF          !  end of bin search loop for this <K1>
        
        
        LOCAT1 = -1         ! key found
        RETURN
        END

