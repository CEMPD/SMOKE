C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION IGSFPT(MY,IV)
C
C  IGSFPT selects the pointer to the gas/diesel sales fractions.
C
C  Called by GSFIDX.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IV
C    common blocks:
C    /GSFCOM/ MAXGSF,MYGSF
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: IGSFPT
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'GSFCOM.I'
      INCLUDE 'REGION.I'
C
      DO 10 IG=1,MAXGSF
      IF(MY.LE.MYGSF(IG,IV,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default for case where MY exceeds highest yr in MYGSF, given IV & IREJN.
C
      IG=1
      CALL QUITER(0.,16,61,MY)
   20 IGSFPT=IG
C
      RETURN
      END
