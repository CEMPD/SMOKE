C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION IRULPT(MY,IVTL)
C
C  IRULPT selects the pointer to the running loss base emission rate for
C  the given model year and vehicle class.
C
C  Called by RULOSS.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IVTL
C    common blocks:
C    /RUNLS3/ MYGRUL,MAXRUL
C
C  Output on return:
C
C  function: IRULPT
C
C  Notes:
C
C  None.
C
C
      COMMON /RUNLS3/ MYGRUL(4,4),RULRVP(4),RULTEM(4),RULSPD(3),MAXRUL
C
      DO 10 IGTL=1,MAXRUL
      IF(MY.LE.MYGRUL(IGTL,IVTL)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IGTL=1
      CALL QUITER(0.,3,61,MY)
   20 IRULPT=IGTL
C
      RETURN
      END
