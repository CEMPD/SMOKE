      LOGICAL FUNCTION PREFIX(STR1, STR2) 
C 
C  Decide whether STR1 is a prefix of STR2 -- that is, whether STR2 
C  begins with a group of characters that is identical to the group of 
C  non-blank characters contained in STR1. Each string is assumed to 
C  to contain only one left-justified group of non-blank characters. 
C 
C  Called by CHKLAB. 
C 
C  Changes: (Last change first) 
C 
C  18 Dec 98 @DynTel-ZK 2-000 Zk Explicit Typing 
C 
C  Input on call: 
C 
C    STR1 - The input string which might be a prefix of STR2. 
C    STR2 - The input string which might contain STR1 as a prefix. 
C 
C  Output on return: 
C 
C    FUNCTION PREFIX is TRUE if STR1 is a prefix of STR2. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  PTR      I    Pointer to current position ins STR1, STR2. 
C 
C  Notes: 
C 
C  PREFIX was added for MOBILE6. 
C 
C 
      IMPLICIT NONE 
      CHARACTER*(*)  STR1 
      CHARACTER*(*)  STR2 
C 
      INTEGER        PTR 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      PREFIX = .FALSE. 
      PTR = 1 
C 
  110 IF (STR1(PTR:PTR).EQ.' '.AND.PTR.GT.1) THEN 
        PREFIX = .TRUE. 
      ELSEIF (STR1(PTR:PTR).NE.STR2(PTR:PTR)) THEN 
        PREFIX = .FALSE. 
      ELSE IF (PTR.EQ.LEN(STR1)) THEN 
        PREFIX = .TRUE. 
      ELSE IF (PTR.EQ.LEN(STR2)) THEN 
        IF (STR1(PTR+1:PTR+1).EQ.' ') THEN   
          PREFIX = .TRUE. 
        ELSE 
          PREFIX = .FALSE. 
        ENDIF 
      ELSE 
        PTR = PTR + 1 
        GOTO 110  
      ENDIF 
C 
      RETURN 
      END 
