      FUNCTION IRNLPT(MY,IVTL) 
C 
C  IRNLPT selects the pointer to the running loss base emission rate for 
C  the given model year and vehicle class. 
C 
C  Called by RNLOSS. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  23 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IVTL 
C    common blocks: 
C    /RUNLS3/ MYGRNL,MAXRNL 
C 
C  Output on return: 
C 
C  function: IRNLPT 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE  
C 
      INCLUDE 'RUNLS3.I' 
C 
      INTEGER IGTL, MY, IVTL, IRNLPT 
C 
      DO 10 IGTL=1,MAXRNL 
      IF(MY.LE.MYGRNL(IGTL,IVTL)) GOTO 20 
   10 CONTINUE 
C 
C  Set default. 
C 
      IGTL=1 
      CALL QUITER(0.,3,61,MY) 
   20 IRNLPT=IGTL 
C 
      RETURN 
      END 
