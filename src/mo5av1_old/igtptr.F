C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION IGTPTR(MY,IVA)
C
C  IGTPTR selects the pointer to the fraction of vehicles with trailer in tow.
C
C  Called by BIGALH.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IVA
C    common blocks:
C    /ALUPLT/ MYGTWC
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: IGTPTR
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
      IF(MY.LE.MYGTWC(IG,IVA,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,9,61,MY)
   20 IGTPTR=IG
C
      RETURN
      END
