      REAL FUNCTION CRANKC(MY,IV) 
C 
C  CRANKC selects the yearly value for crankcase HC. 
C 
C  Called by EVP_HRLY. 
C 
C  Calls QUITER. 
C 
C  Changes : 
C 
C 11 Oct 2000 @DynTel-BG separate model year variable, QTR_MY, for 
C        fatal error call to QUITER 
C 12 Dec 1999 @DynTel-BG 1-006 defined function arguments with INTENT 
C             statements 
C  3 Sep 1998 @DynTel-ZK 2-000 Explicit Typing  
C 
C  Input on call: 
C 
C    argument list: MY,IV 
C    common blocks: 
C    /CRACOM/ MAXCRA,CRAVAL,MYGCRA 
C    /REGION/ IREJN 
C 
C  Output on return: 
C 
C    function: CRANKC 
C 
C  Notes. 
C 
      IMPLICIT NONE  
      INCLUDE 'CRACOM.I' 
      INCLUDE 'REGION.I' 
C 
C 
C     Declare parameters. 
C 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
      INTEGER :: QTR_MY 
C 
      INTEGER IG 
C 
      QTR_MY=MY 
      DO 10 IG=1,MAXCRA 
      IF(MY.LE.MYGCRA(IG,IV,IREJN)) GOTO 20 
   10 CONTINUE 
C 
C  Set default. 
C 
      IG=1 
      CALL QUITER(CRAVAL(IG,IV,IREJN),5,60,QTR_MY) 
   20 CRANKC=CRAVAL(IG,IV,IREJN) 
C 
      RETURN 
      END 
