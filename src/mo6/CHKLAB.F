      SUBROUTINE CHKLAB(INLAB,STDLAB,LABTYPE,SECTION) 
C 
C  Interpret the label text and return the full standard text for the 
C  label, as well as a token number corresponding to the label type. 
C 
C  Called by MOBILE. 
C 
C  Calls NXTTOK, PREFIX 
C 
C  Input on call: 
C 
C    INLAB   - Label entered by the user in the Mobile6 input file. 
C 
C  Output on return: 
C 
C    LABTYPE - An integer that identifies the label type. 
C    SECTION - Position for this label in the input (global data input, 
C              scenario input, either global or scenario, or header). 
C    STDLAB  - Full (unabbreviated, uppercase, left justified) label text. 
C 
C  Changes : (Last change first) 
C 
C  27 Jul 99 @DynTel-MLA 2-000 Removed CHKLAB from the module LABELS and 
C            made it an ordinary external procedure (not a module proc.). 
C 
C  Module variables: 
C 
C    LABLIST - Array of derived type LABELDEF (see module LABELS). 
C    NLABS   - Dimension of the LABELDEF array (see module LABELS). 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  TOK1     C    Next 'word' from the input record label. 
C  TOK2     C    Next 'word' from the MINLAB string. 
C  TOK3     C    Next 'word' from the FULLAB string. 
C  T1END    I    Ending position of TOK1 in the input record string. 
C  T1ST     I    Starting position of TOK1 in the input record string. 
C  T2END    I    Ending position of TOK2 in the MINLAB string. 
C  T2ST     I    Starting position of TOK2 in the MINLAB string. 
C  T3END    I    Ending position of TOK3 in the FULLAB string. 
C  T3ST     I    Starting position of TOK3 in the FULLAB string. 
C 
C  Notes: 
C 
C  CHKLAB was added for MOBILE6. 
C 
      USE LABELS 
      IMPLICIT NONE 
 
      CHARACTER(LEN = *), INTENT(IN)   :: INLAB 
      INTEGER, INTENT(OUT)             :: LABTYPE 
      INTEGER, INTENT(OUT)             :: SECTION 
      CHARACTER(LEN = *), INTENT(OUT)  :: STDLAB 
C 
      INTEGER             ::  I 
      CHARACTER (LEN=40)  ::  TOK1 
      CHARACTER (LEN=40)  ::  TOK2 
      CHARACTER (LEN=40)  ::  TOK3 
      INTEGER             ::  T1END 
      INTEGER             ::  T1ST 
      INTEGER             ::  T2END 
      INTEGER             ::  T2ST 
      INTEGER             ::  T3END 
      INTEGER             ::  T3ST 
C 
C 
      LOGICAL  ::  PREFIX 
C 
C. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     Loop through the labels. Look for one that matches the label on 
C     the record read from the input file. 
C 
  100 DO 120 I = 1,NLABS 
C 
        T1ST = 0 
        T1END = 0 
        T2ST = 0 
        T2END = 0 
        T3ST = 0 
        T3END = 0 
C 
C       Compare the user input wuth the MIN and FULL labels. 
C       Return the token number if there is a match... (or zero). 
C 
  110   CALL NXTTOK(INLAB,TOK1,T1ST,T1END) 
        CALL NXTTOK(LABLIST(I)%MINLAB,TOK2,T2ST,T2END) 
        CALL NXTTOK(LABLIST(I)%FULLAB,TOK3,T3ST,T3END) 
C 
        IF (TOK1.EQ.' '.AND.TOK2.EQ.' ') THEN 
          LABTYPE = LABLIST(I)%LABID 
          STDLAB = LABLIST(I)%FULLAB 
          SECTION = LABLIST(I)%SECTION 
          GOTO 190 
        ELSEIF (TOK1.EQ.' '.OR.TOK2.EQ.' ') THEN 
          CYCLE 
        ELSEIF (PREFIX(TOK2,TOK1).AND.PREFIX(TOK1,TOK3)) THEN 
          GOTO 110 
        ELSE 
          CYCLE 
        ENDIF 
C 
  120 CONTINUE 
C 
C     If we fall out of the loop without finding any match, set LABTYPE 
C     to zero and return. 
C 
      LABTYPE = 0 
      SECTION = 0 
      STDLAB = ' *** ERROR ***' 
C 
  190 RETURN 
      END SUBROUTINE CHKLAB 
