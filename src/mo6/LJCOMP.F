       SUBROUTINE LJCOMP(STR1, STR2) 
C 
C  This subroutine left justifies a string, compressing multiple embedded  
C  blanks to single blanks, but does not convert the string to upper case. 
C 
C  Called by INISCN. 
C 
C  Changes: (Last change first) 
C 
C  19 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    STR1 - The input string to be converted. 
C 
C  Output on return: 
C 
C    STR2 - A left-justified copy of STR1, blanks compressed. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CH       C    The next character from STR1. 
C  PREV     C    The last character added to STR2. 
C  PTR1     I    Current character position in STR1. 
C  PTR2     I    Current character position in STR2. 
C 
C  Notes: 
C 
C    LJCOMP was added for MOBILE6.                      
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      CHARACTER*(*)  STR1 
      CHARACTER*(*)  STR2 
C 
      CHARACTER      CH 
      CHARACTER      PREV 
      INTEGER        PTR1 
      INTEGER        PTR2 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     Clear string 2. Initialize so leading blanks are not copied. 
C     The first nonblank character goes into position 1 in string 2. 
C 
  001 STR2 = ' ' 
      PREV = ' ' 
      PTR2 = 1 
C 
C     Copy characters from string 1 to string 2. Skip leading blanks and  
C     compress multiple embedded blanks to single blanks. 
C 
      DO 20 PTR1 = 1,LEN(STR1) 
C 
        CH = STR1(PTR1:PTR1) 
        IF (CH.EQ.' '.AND.PREV.EQ.' ') THEN 
          CYCLE 
        ENDIF 
C 
        STR2(PTR2:PTR2) = CH 
        PREV = CH 
        PTR2 = PTR2 + 1 
        IF (PTR2.GT.LEN(STR2)) GOTO 990 
C 
   20 CONTINUE 
 
  990 RETURN 
      END 
 
 
 
 
 
 
 
 
 
 
