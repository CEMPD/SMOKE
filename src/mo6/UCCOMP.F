      SUBROUTINE UCCOMP(STR1, STR2) 
C 
C  This subroutine converts a character string to upper case and left 
C  justifies it, compressing multiple embedded blanks to single blanks. 
C 
C  Called by DSPTCH,GETALH,GETERQ,GETRFG,GETRJN,GETSPR,GETTAM, 
C            GETWIN,FNDLAB, and OUTHDR. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  13 January 1999 @DynTel-zk 2-000 Explicit Typing 
C 
c  Input on call: 
C 
C    STR1 - The input string to be converted. 
C 
C  Output on return: 
C 
C    STR2 - An uppercase, left-justified copy of STR1, blanks compressed. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CH       C    The next character from STR1. 
C  LOWCASE  C    String containing the lower case letters. 
C  PREV     C    The last character added to STR2. 
C  PTR1     I    Current character position in STR1. 
C  PTR2     I    Current character position in STR2. 
C  UCPTR    I    Pointer to parallel arrays of upper/lower case letters. 
C  UPCASE   C    String containing the upper case letters. 
C 
C  Notes: 
C 
C    UCCOMP was added for MOBILE6.                      
C 
C  Notes: 
      IMPLICIT NONE 
C 
      CHARACTER*(*)  STR1 
      CHARACTER*(*)  STR2 
C 
      CHARACTER      CH 
      CHARACTER      PREV 
      INTEGER        PTR1 
      INTEGER        PTR2 
      INTEGER        UCPTR 
      CHARACTER(LEN=26) :: UPCASE='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
      CHARACTER(LEN=26) :: LOWCASE='abcdefghijklmnopqrstuvwxyz' 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     Clear string 2. Initialize so leading blanks are not copied. 
C     The first nonblank character goes into position 1 in string 2. 
C 
  001 STR2 = ' ' 
      PREV = ' ' 
      PTR2 = 1 
C 
C     Copy characters from string 1 to string 2. Convert alphabetic 
C     characters to upper case. Skip leading and multiple blanks. 
C 
      DO 20 PTR1 = 1,LEN(STR1) 
C 
        CH = STR1(PTR1:PTR1) 
        IF (CH.EQ.' '.AND.PREV.EQ.' ') THEN 
          CYCLE 
        ENDIF 
C 
        DO 10 UCPTR = 1,26 
          IF (CH.EQ.LOWCASE(UCPTR:UCPTR)) THEN
             CH = UPCASE(UCPTR:UCPTR) 
            EXIT 
          ENDIF 
   10   CONTINUE 
C 
        STR2(PTR2:PTR2) = CH 
        PREV = CH 
        PTR2 = PTR2 + 1 
        IF (PTR2.GT.LEN(STR2)) EXIT 
C 
   20 CONTINUE 
 
  990 RETURN 
      END 
