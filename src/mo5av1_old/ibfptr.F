C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION IBFPTR(MY,IP,IVLM)
C
C  IBFPTR selects the pointer to the normalized bag fractions in BFSTEP and
C  BFR50K (if 1981+ LDGV) or in BFRCOF (all other MY x IVLM cases).
C
C  Called by BIGIN2 and BEFIDL.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IP,IVLM
C    common blocks:
C    /OPMOD2/ MYSTEP,MAXBFR,MYGBFR
C
C  Output on return:
C
C    function: IBFPTR
C
C  Notes:
C
C  None.
C
C
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,6)
C
      IF(MYSTEP.EQ.2) GOTO 20
C
      DO 10 IG=1,MAXBFR
      IF(MY.LE.MYGBFR(IG,IP,IVLM)) GOTO 30
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,10,61,MY)
      GOTO 30
C
C  1981+ LDGV case: coefficients change every year 1981 through 1992.
C  1984+ LDGT1/2 case: coefficients change every year 1984 through 1992.
C  1993 - 2020 are same as 1992's.  Figure index accordingly.
C
   20 IF(MY.LE.1992) IG=MY-1981+1
      IF(MY.GE.1993) IG=12
C
   30 IBFPTR=IG
C
      RETURN
      END
