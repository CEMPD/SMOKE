C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION IDLPTR(MY,IP,IV)
C
C  IDLPTR selects the model year group pointer into the slope & intercept
C  arrays for the idle emission rate equation.
C
C  Called by BEFIDL.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IP,IV
C    common blocks:
C    /IDLEQ2/ MYGIDL,MAXIDL
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: IDLPTR
C
C  Notes:
C
C  None.
C
C
      COMMON /IDLEQ2/ MYGIDL(15,3,8,2),MAXIDL
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      DO 10 IG=1,MAXIDL
      IF(MY.LE.MYGIDL(IG,IP,IV,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,17,61,MY)
C
   20 IDLPTR=IG
      RETURN
      END
