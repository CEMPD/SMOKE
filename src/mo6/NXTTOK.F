      SUBROUTINE NXTTOK(STR1, STR2, PTR1, PTR2) 
C 
C  This subroutine returns the next blank-delimited token from a string. 
C 
C  Called by CHKLAB. 
C 
C  Changes: (Last change first) 
C 
C 30 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 27 Oct 98 @DynTel-MLA 2-000  Blanks are returned (instead of a special 
C           sentinel character) when the end of STR1 is reached. The 
C           subroutine can be called repeatedly after the end of STR1 is 
C           reached and it will return blanks without causing an error. 
C 
C  Input on call: 
C 
C    STR1 - The input string to be parsed. 
C    PTR2 - The last character of the previous token returned from STR1. 
C 
C  Output on return: 
C 
C    STR2 - The next blank-delimited string of characters from STR1. 
C    PTR1 - Position in STR1 of the first character of the token returned  
C           in STR2. 
C    PTR2 - Position in STR1 of the last character of the token returned  
C           in STR2 (or the following blank character). 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  BUFFER   C    Local working storage. 
C  ERROR    C    Error string returned in STR2 for all error conditions. 
C  I        I    Loop control variable. 
C 
C  Notes: 
C 
C  NXTTOK was added for MOBILE6. 
C 
      IMPLICIT NONE 
      INTEGER        PTR1 
      INTEGER        PTR2 
      CHARACTER*(*)  STR1 
      CHARACTER*(*)  STR2 
C 
      CHARACTER      BUFFER*200 
      CHARACTER*13   ERROR 
      INTEGER        I 
C 
      PARAMETER      (ERROR = '*** ERROR ***') 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     PTR2 should be pointing to the last character in the last token 
C     found in STR1. (It should be set to zero when this function is 
C     called to return the first token in a string.) The initial value 
C     of PTR1 is irrelevant. It is set by this function. 
C 
  100 IF ((LEN(STR1)+2).GT.LEN(BUFFER)) THEN 
        GOTO 180  
      ENDIF 
C 
      BUFFER = STR1 
      BUFFER(LEN(BUFFER):LEN(BUFFER)) = ' ' 
C 
      IF (PTR2.GE.LEN(BUFFER)) THEN 
        PTR1 = LEN(BUFFER) 
      ELSE 
        PTR1 = PTR2 + 1 
      ENDIF 
C 
      DO WHILE (PTR1.LT.LEN(BUFFER) .AND. BUFFER(PTR1:PTR1).EQ.' ' ) 
          PTR1 = PTR1 + 1 
      END DO 
C 
      PTR2 = PTR1 
  120 DO WHILE (PTR2.LT.LEN(BUFFER) .AND. BUFFER(PTR2:PTR2).NE.' ' ) 
          PTR2 = PTR2 + 1 
      END DO 
C 
      IF ((PTR2-PTR1).GE.LEN(STR2)) THEN 
        GOTO 180  
      ENDIF 
C 
      STR2 = BUFFER(PTR1:PTR2) 
C 
      GOTO 190  
C 
C     For all error conditions, return the error string in STR2. 
C 
  180 IF (LEN(STR2).GE.LEN(ERROR)) THEN 
        STR2 = ERROR 
      ELSE 
        DO 185 I = 1, LEN(STR2) 
          STR2(I:I) = ERROR(I:I) 
  185   CONTINUE 
      ENDIF 
C 
  190 RETURN 
      END 
