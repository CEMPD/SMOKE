C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION DISCAL(LIM,LKINK,ID,IY,ICY)
C
C  DISCAL computes the disablement rate for the given vehicle type, model
C  year interval, I/M program (or lack thereof), mileage distribution, and
C  disablement class.
C
C  Called by DISIMP.
C
C  Calls ITAMPT.
C
C  Input on call:
C
C    parameter list: LIM,LKINK,ID
C    common blocks:
C    /LOOKUP/ IVTAM
C    /MYCODE/ MY
C    /TAMEQ1/ TAMZML,TAMDR,TAMA50
C    /TAMEQ2/ F50K
C    /TAMEQ3/ CUMALL,CUMNIM,CUMIM,GT50KA,GT50KN,GT50KI
C
C  Output on return:
C
C    function: DISCAL
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C                         tampering class subject to 50,000 mile DR kink?"
C  M5KINK    I    switch: 1 = no & 2 = yes to query "50,000 kink applies?"
C                         (myg 2 + not misfueling + 50K+ + 50K+ slope steeper)
C
C  Notes:
C
C  M5KINK logic assumes TAMDR slopes are all positive.
C  DISCAL was modified for MOVILE4.1 v4 &' v5 to change the over 50K cases
C  from a set of rates to an adjustment (TAMA50) to the under 50K rates.
C
C
C
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMEQ1/ TAMZML(9,4,3,2),TAMDR(9,4,3,2),TAMA50(9,4,3,2)
      COMMON /TAMEQ2/ MYGTAM(3,4),F50K
      COMMON /TAMEQ3/ CUMCAP,CUMALL,CUMNIM,CUMIM,GT50KA,GT50KN,GT50KI
C
C  The slope for MYG 2 changes at 50,000 miles, except for the 3 misfueling
C  categories.  The over 50K slopes are stored as in TAMA50.
C
      IGD=ITAMPT(6)
C
      M5KINK=1
      IF(LKINK.EQ.2) GOTO 10
C
C  Case 1: a single rate line is used for the entire model year interval:
C          (1) all mileage is on the non-I/M rate line, or
C          (2) all mileage is on the     I/M rate line.
C
      IF(GT50KA.GT.0.0) M5KINK=2
      DISCAL=TAMZML(ID,IVTAM,IGD,LIM)
C
      IF(M5KINK.EQ.1) DISCAL=DISCAL
     *    +TAMDR(ID,IVTAM,IGD,LIM)*CUMALL
C
      IF(M5KINK.EQ.2) DISCAL=DISCAL
     *    +TAMDR(ID,IVTAM,IGD,LIM)*F50K
     *    +TAMDR(ID,IVTAM,IGD,LIM)*TAMA50(ID,IVTAM,IGD,LIM)*GT50KA
C
      RETURN
C
C  Case 2: An I/M program starts within the interval, covering the model year
C          and vehicle class being evaluated.  Both the non-I/M and I/M rate
C          lines are used - the I/M "kink": non-I/M line intercept, first
C          mileage segment times the non-I/M slope & second times the I/M slope.
C          The resulting rate may be replaced with other estimates, depending
C          on the 2 lines and at what point they are evaluated.
C
   10 IF(ICY.LE.1990.AND.GT50KI.GT.0.0) THEN
        M5KINK=2
      ELSEIF(ICY.EQ.1991.AND.IY.EQ.2.AND.GT50KI.GT.0.0) THEN
        M5KINK=2
      ELSEIF((GT50KI.GE.0.0)
     *      .AND.(ICY.GE.1992.OR.(IY.EQ.1.AND.1991.EQ.ICY))) THEN
        M5KINK=2
      ENDIF
      DISCAL=TAMZML(ID,IVTAM,IGD,1)
C
      IF(M5KINK.EQ.1) DISCAL=DISCAL
     *    +TAMDR(ID,IVTAM,IGD,1)*CUMNIM
     *    +TAMDR(ID,IVTAM,IGD,2)*CUMIM
C
      IF(M5KINK.EQ.2.AND.GT50KN.GT.0.0) DISCAL=DISCAL
     *    +TAMDR(ID,IVTAM,IGD,1)*F50K
     *    +TAMDR(ID,IVTAM,IGD,1)*TAMA50(ID,IVTAM,IGD,1)*GT50KN
     *    +TAMDR(ID,IVTAM,IGD,2)*TAMA50(ID,IVTAM,IGD,2)*CUMIM
C
      IF(M5KINK.EQ.2.AND.GT50KN.EQ.0.0) DISCAL=DISCAL
     *    +TAMDR(ID,IVTAM,IGD,1)*CUMNIM
     *    +TAMDR(ID,IVTAM,IGD,2)*(F50K-CUMNIM)
     *    +TAMDR(ID,IVTAM,IGD,2)*TAMA50(ID,IVTAM,IGD,2)*GT50KI
C
      RETURN
      END
