C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION ITCPTR(MY,IP,IVTEMP)
C
C  ITCPTR selects the pointer to the temperature correction coefficients.
C
C  Called by BIGCFX.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IP,IVTEMP
C    common block: /TEMPC2/ MAXTCF,MYGTCF
C
C  Output on return:
C
C    function: ITCPTR
C
C  Notes:
C
C  None.
C
C
      COMMON /TEMPC2/ MAXTCF,MYGTCF(9,3,5),MYGTFD(3),ISHIFT
C
      DO 10 IG=1,MAXTCF
      IF(MY.LE.MYGTCF(IG,IP,IVTEMP)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,15,61,MY)
   20 ITCPTR=IG
C
      RETURN
      END
