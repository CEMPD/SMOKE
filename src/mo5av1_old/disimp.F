C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE DISIMP(LJDX,LAY,LIM,LKINK,IY,ICY)
C
C  DISIMP takes the passed in ATP case, figures out what I/M case applies,
C  calls DISCAL to calculate it and if necessary substitutes Non-I/M / I/M
C  line limit values.
C
C  Called by DISATP.
C
C  Calls AMAX1 and AMIN1 (FORTRAN library function) and DISCAL.
C
C  Input on call:
C
C    parameter list: LJDX,LAY,LIM,LKINK
C    common blocks:
C    /CUMCOM/ CUMMIL
C    /LOOKUP/ IVTAM
C    /MYCODE/ JDX,IMDXSY
C    /TAMEQ2/ F50K
C    /TAMEQ3/ CUMCAP
C
C  Output on return:
C
C    common blocks:
C    /TAMEQ1/ F50K
C    /TAMEQ2/ BTR
C    /TAMEQ3/ CUMALL,CUMNIM,CUMIM,GT50KA,GT50KN,GT50KI
C    /TAMEQ4/ BTR
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  DALLIM    R    disablement rate for the all I/M mileage segment case
C  DKINKI    R    disablement rate for the kinked I/M mileage segment case
C  DNONIM    R    disablement rate for the all non-I/M mileage segment case
C
C  Notes:
C
C  CUMCAP was added for MOBILE4.1 Version 02.
C
C
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMEQ1/ TAMZML(9,4,3,2),TAMDR(9,4,3,2),TAMA50(9,4,3,2)
      COMMON /TAMEQ2/ MYGTAM(3,4),F50K
      COMMON /TAMEQ3/ CUMCAP,CUMALL,CUMNIM,CUMIM,GT50KA,GT50KN,GT50KI
      COMMON /TAMEQ4/ BTR(9,2)
C
C  Lookup and scale the needed cumulative mileage segments.  Based on these
C  segments, compute the mileage on each segment greater than 50,000 miles.
C
      CUMALL=AMIN1(CUMCAP,CUMMIL(LJDX,IVTAM)*.0001)
      GT50KA=AMAX1(0.0,CUMALL-F50K)
      IF(ICY.EQ.1991.AND.IY.EQ.1) THEN
        CUMNIM=CUMALL
        CUMIM=0.0
        GT50KI=0.0
        GT50KN=0.0
      ELSE
        IF(ICY.GE.1992) THEN
          CUMNIM=CUMALL
          CUMIM=0.0
          GT50KI=0.0
          GT50KN=0.0
        ENDIF
      ENDIF
      IF(LKINK.EQ.1) GOTO 10
C
      CUMNIM=AMIN1(CUMCAP,CUMMIL(JDX-IMDXSY+1,IVTAM)*.0001)
      IF(ICY.EQ.1991.AND.IY.EQ.1.AND.CUMNIM.GT.CUMALL) THEN
        CUMNIM=CUMALL
      ELSE
        IF(ICY.GE.1992.AND.CUMNIM.GT.CUMALL) CUMNIM=CUMALL
      ENDIF
      CUMIM=CUMALL-CUMNIM
      GT50KN=AMAX1(0.0,CUMNIM-F50K)
      GT50KI=AMAX1(0.0,GT50KA-GT50KN)
C
   10 DO 40 ID=1,9
C
C  Case 1: No I/M program applies during the model year interval.
C
      DNONIM=DISCAL(1,1,ID,IY,ICY)
      IF(LIM.EQ.2) GOTO 20
      BTR(ID,LAY)=DNONIM
      GOTO 40
C
C  Case 2: I/M program applies to all years of the model year interval.
C
   20 DALLIM=DISCAL(2,1,ID,IY,ICY)
      IF(LKINK.EQ.2) GOTO 30
      BTR(ID,LAY)=DALLIM
      IF(DALLIM.GT.DNONIM) BTR(ID,LAY)=DNONIM
      GOTO 40
C
C  Case 3: I/M program starts to apply during interval -> kink in line's slope.
C
   30 DKINKI=DISCAL(LIM,LKINK,ID,IY,ICY)
      BTR(ID,LAY)=DKINKI
      IF(DKINKI.LT.DALLIM) BTR(ID,LAY)=DALLIM
      IF(DKINKI.GT.DNONIM) BTR(ID,LAY)=DNONIM
C
   40 CONTINUE
C
      RETURN
      END
