
C Version "@(#)$Id$ $Source$ $Date$ 

        INTEGER FUNCTION LOCAT3( K1, K2, K3, N, LIST1, LIST2, LIST3 )

C***********************************************************************
C  function body starts at line 48
C
C  RETURNS:
C       subscript at which the targeted key-triple should be inserted
C       or -1 if the key-triple is found.
C       EXAMPLE:  search for <FIP,ID7,ID3> in table of FIP and ASC values
C       (where ID7 is leading-7 digits and ID3 is trailing 3 digits
C       in a 10-digit Area Source Code).
C
C  PRECONDITIONS REQUIRED:
C       Sorted table <N,LIST1,LIST3,LIST3> for searching
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 12/95 by MRH copied from FIND3 of CJC
C
C***********************************************************************

      IMPLICIT NONE


C...........   ARGUMENTS and their descriptions:
        
        INTEGER  K1             !  first  key
        INTEGER  K2             !  second key
        INTEGER  K3             !  third  key
        INTEGER  N              !  table size
        INTEGER  LIST1( N )     !  table to search for K1
        INTEGER  LIST2( N )     !  table to search for K2
        INTEGER  LIST3( N )     !  table to search for K3


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER  LO
        INTEGER  HI
        INTEGER  M


C***********************************************************************
C   begin body of function  LOCAT3

        LO = 1
        HI = N

        IF( N .EQ. 0 ) THEN

            LOCAT3 = -1
            RETURN

        ENDIF 

11      CONTINUE
           
            IF ( LO .GT. HI ) THEN
            
                LOCAT3 = LO
                RETURN
                
            END IF
           
            M = ( LO + HI ) / 2
            IF ( K1 .GT. LIST1( M ) ) THEN
                LO = M + 1
                GO TO  11
            ELSE IF ( K1 .LT. LIST1( M ) ) THEN
                HI = M - 1
                GO TO  11
            ELSE IF ( K2 .GT. LIST2( M ) ) THEN
                LO = M + 1
                GO TO  11
            ELSE IF ( K2 .LT. LIST2( M ) ) THEN
                HI = M - 1
                GO TO  11
            ELSE IF ( K3 .GT. LIST3( M ) ) THEN
                LO = M + 1
                GO TO  11
            ELSE IF ( K3 .LT. LIST3( M ) ) THEN
                HI = M - 1
                GO TO  11
            END IF          !  end of bin search loop for this <K1,K2,K3>
        
        
        LOCAT3 = -1         ! key-tuple found
        RETURN
        END

