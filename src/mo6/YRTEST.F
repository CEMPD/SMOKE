      SUBROUTINE YRTEST(IY,NERR,IY1,IY2,INERR) 
C 
C  YRTEST converts 2-digit input year to 4-digit year. 
C 
C  Called by GETATP, GETERQ, GETIMP, GETNAT, PARSEC and RDCITY. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  12 Oct 00 AIR Task 02: Bug Fix 279: Allow 2 and 4 digit input years 
C  15 January 1999 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IY,NERR,IY1,IY2,INERR 
C    common blocks: 
C    /YEARS4/ IYSTART 
C 
C  Output on return: 
C 
C    argument list: IY,INERR 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  II       I    IY converted to 4 digit year 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'YEARS4.I' 
C 
C     declare parameter list 
C 
      INTEGER INERR 
      INTEGER IY 
      INTEGER IY1 
      INTEGER IY2 
      INTEGER NERR 
C 
C     declare local variables 
C 
      INTEGER II 
C 
      II=IY 
      IF(IY.LT.100) II=IY+1900 
      IF(II.LT.IYSTART) II=II+100 
      IF(IY.LT.0.OR.II.LT.IY1.OR.II.GT.IY2) 
     *   CALL QUITER(0.,IY,NERR,INERR) 
      IY=II 
      RETURN 
      END 
