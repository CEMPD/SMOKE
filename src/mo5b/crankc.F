C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION CRANKC(MY,IV)
C
C  CRANKC selects the yearly value for crankcase HC.
C
C  Called by CCEVRT.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IV
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
C
      INCLUDE 'CRACOM.I'
      INCLUDE 'REGION.I'
C
      DO 10 IG=1,MAXCRA
      IF(MY.LE.MYGCRA(IG,IV,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(CRAVAL(IG,IV,IREJN),5,60,MY)
   20 CRANKC=CRAVAL(IG,IV,IREJN)
C
      RETURN
      END
