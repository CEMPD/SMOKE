C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION IGEPTR(MY,IVA)
C
C  IGEPTR selects the pointer to the fraction of vehicles equipped with a/c.
C
C  Called by BIGALH.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IVA
C    common blocks:
C    /ALUPAC/ MYGACE
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: IGEPTR
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'ALUPAC.I'
      INCLUDE 'REGION.I'
C
      DO 10 IG=1,8
      IF(MY.LE.MYGACE(IG,IVA,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,6,61,MY)
   20 IGEPTR=IG
C
      RETURN
      END
