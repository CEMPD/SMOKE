C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE EXMPTR(MY,IV)
C
C  EXMPTR selects the pointer to the non-methane fractions.
C
C  Called by BIGCFX.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IV
C    common blocks:
C    /NMETH2/ MYGMTH,MTHMYG
C    /REGION/ IREJN
C
C  Output on return:
C
C
C  Notes:
C
C  EXMPTR was modified for MOBILE4.1.
C
C
      COMMON /NMETH2/ MYGMTH(7,8,2),MTHMYG
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      DO 10 IG=1,7
      IF(MY.LE.MYGMTH(IG,IV,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,14,61,MY)
C
   20 MTHMYG=IG
C
      RETURN
      END
