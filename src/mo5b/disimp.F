C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE DISIMP(LJDX,LAY,LIM,LKINK)
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
C  Sep-20-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar
C  year removed IY and ICY from the parameter list and thus from calls
C  to DISCAL
C
C
      INCLUDE 'CUMCOM.I'
      INCLUDE 'LOOKUP.I'
      INCLUDE 'MYCODE.I'
      INCLUDE 'TAMEQ1.I'
      INCLUDE 'TAMEQ2.I'
      INCLUDE 'TAMEQ3.I'
      INCLUDE 'TAMEQ4.I'
C
C  Lookup and scale the needed cumulative mileage segments.  Based on these
C  segments, compute the mileage on each segment greater than 50,000 miles.
C
      CUMALL=AMIN1(CUMCAP,CUMMIL(LJDX,IVTAM)*.0001)
      GT50KA=AMAX1(0.0,CUMALL-F50K)
      CUMNIM=CUMALL
      CUMIM=0.0
      GT50KI=0.0
      GT50KN=0.0
      IF(LKINK.EQ.1) GOTO 10
C
      CUMNIM=AMIN1(CUMCAP,CUMMIL(JDX-IMDXSY+1,IVTAM)*.0001)
      IF(CUMNIM.GT.CUMALL) CUMNIM=CUMALL
      CUMIM=CUMALL-CUMNIM
      GT50KN=AMAX1(0.0,CUMNIM-F50K)
      GT50KI=AMAX1(0.0,GT50KA-GT50KN)
C
   10 DO 40 ID=1,9
C
C  Case 1: No I/M program applies during the model year interval.
C
      DNONIM=DISCAL(1,1,ID)
      IF(LIM.EQ.2) GOTO 20
      BTR(ID,LAY)=DNONIM
      GOTO 40
C
C  Case 2: I/M program applies to all years of the model year interval.
C
   20 DALLIM=DISCAL(2,1,ID)
      IF(LKINK.EQ.2) GOTO 30
      BTR(ID,LAY)=DALLIM
      IF(DALLIM.GT.DNONIM) BTR(ID,LAY)=DNONIM
      GOTO 40
C
C  Case 3: I/M program starts to apply during interval -> kink in line's slope.
C
   30 DKINKI=DISCAL(LIM,LKINK,ID)
      BTR(ID,LAY)=DKINKI
      IF(DKINKI.LT.DALLIM) BTR(ID,LAY)=DALLIM
      IF(DKINKI.GT.DNONIM) BTR(ID,LAY)=DNONIM
C
   40 CONTINUE
C
      RETURN
      END
