C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION IGLPTR(MY,IVA)
C
C  IGLPTR selects the pointer to the fraction of vehicles with extra load.
C
C  Called by BIGALH.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IVA
C    common blocks:
C    /ALUPLT/ MYGXLC
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: IGLPTR
C
C  Notes:
C
C  None.
C
C
      COMMON /ALUPLT/ MYGXLC(6,3,2),MYGTWC(6,3,2)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      DO 10 IG=1,6
      IF(MY.LE.MYGXLC(IG,IVA,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,8,61,MY)
   20 IGLPTR=IG
C
      RETURN
      END
