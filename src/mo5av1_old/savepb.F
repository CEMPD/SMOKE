C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE SAVEPB
C
C  SAVEPB computes and saves the misfueling category BTRs for all
C  IVTAM x IDX cases, using LPOD = the last year of leaded gasoline, not ICY.
C  FINDPB will then use the resulting table to get the LPOD misfueling rates
C  with which to replace the corresponding ICY rates.  The latter are in
C  error because DISATP in using ICY assumes misfueling continues after LPOD.
C
C  Called by ONESEC.
C
C  Calls DISATP, and MYRIMT.
C
C  Input on call:
C
C    common blocks:
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG
C    /IMPAR1/ ICYIM
C    /IMPAR2/ ITEST
C    /MAXIMA/ MAXYRS
C    /TAMPB1/ LPOD
C    /YEARS4/ IY2020
C
C  Output on return:
C
C    common blocks:
C    /LOOKUP/ IVTAM
C    /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IMDXSY,IMKINK
C    /TAMPB1/ PBBTR
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IM       I    flag relating run's I/M program to relative model year
C                IDX:  1 = no I/M   2 = I/M for this IDX, given ICY,
C                IVTAM, and I/M parameter.
C  IMXSY    I    Same definition as in TAMPER.
C  IY       I    Same definition as in TAMPER.
C  LDX1ST   I    JDX order pointer to first model year covered by ATP
C  LDXIV    I    Same definition as in TAMPER.
C  LDXLST   I    JDX order pointer to last model year covered by ATP
C  LMYR     I    flag relating model year to ATP parameters
C                  1 = model year not covered by an ATP
C                  2 = JDX is covered by an ATP (but IVTAM may not be,
C                      hence, check for inclusion to set LMYRVT)
C  LMYRIM   I    flag relating model year, vehicle and calendar year to
C                (LPOD) I/M programs(s). 1=I/M program 1 and 2=I/M
C                program 2 while a 0=no I/M program for MYR,IPOD,IVTAM.
C
C
C  Notes:
C
C  For MOBILE4.1 SAVEPB was created.  SAVEPB is identical to
C  TAMPER, except that LPOD is used instead of ICY and SAVEPB is
C  identical to TAMPER, except that LPOD is used instead of ICY
C  and the save to PBBTR replaces the calls to EEFGRP, EMIRAT,
C  BAGEME, and EVPEME.
C  For MOBILE5 SAVEPB was modified to include the IM240 effect on
C  tampering.  Needed to add a call to use the MYRIMT function.
C  March 10, 1993 IY=1 was added to the parameter list to the
C  call to DISATP to guarantee the 1990 inventory emission
C  remain the same in newer versions of MOBILE as in MOBILE5.
C
C
      COMMON /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(4)
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMEQ4/ BTR(9,2)
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,4)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C  Relate start year of anti-tampering program (ATP) and 1st and
C  last model years covered to calendar year.  Put them in JDX
C  order.  For example:
C
C      LPOD = LAPSY => LDXSY = 1
C      LPOD = LAPSY + 1 => LDXSY = 2
C          .
C          .
C      LPOD = LAPSY + 24 => LDXSY = 25
C
C  JDX order is used because these values are used to access the JDX ordered
C  CUMMIL.  Note that LDXSY <= 0 => LPOD < LAPSY => no ATP was in effect during
C  any part of the calendar year's 25 model year window.
C
C  The ATP has to be in effect at least 1 year prior to LPOD for it to affect
C  LPOD's emissions.  The mileage accumulated through LAPSY all goes to the no
C  program case.
C
C  Only calculate these values once per scenario => compute them here before
C  the tamper loops begin.
C
      LDXSY=LPOD-LAPSY+1
      LDX1ST=LPOD-LAP1ST+1
      LDXLST=LPOD-LAPLST+1
C
C  For each tampering vehicle type and each model year in the 24 year window
C  from LPOD-23 thru LPOD (25+ year case is always zeroed), compute the bag and
C  evaporative tampering offsets.
C
      DO 90 IVTAM=1,4
      DO 80 IDX=2,MAXYRS
C
C  MY pointer info using IDX & LPOD may be done here, before sub calls.
C  Identify this loop's model year for use in ITAMPT lookups.
C
      MY=LPOD-(MAXYRS-IDX)
      LMYRIM=MYRIMT(MY,IVTAM,LPOD)
      IMDXSY=LPOD-IY2020+1
      IF(LMYRIM.NE.0) IMDXSY=LPOD-ICYIM(LMYRIM)+1
C
C  Assign the reverse order (LPOD -> 1, increases with age) my index for use in
C  apportioning mileage weights.
C
      JDX=MAXYRS-IDX+1
C
C  Reset (turn off) the ATP my range switch.  Turn on if there is an ATP
C  program, it starts before the calendar year, the model year JDX is in the
C  affected range and the vehicle type has been selected by the user for
C  coverage.
C
C      LMYRVT = 1 => switch is off.  ATP not a factor in my JDX's tampering
C      LMYRVT = 2 => switch is ON.  ATP is a factor in my JDX's tampering
C
C  Note: ATPFLG = 1 => LAPSY = 2020 => LDXSY < 1 => 1st 2 criteria can be
C        checked using LDXSY.
C
      LDXSY=LPOD-LAPSY+1
      LDX1ST=LPOD-LAP1ST+1
      LDXLST=LPOD-LAPLST+1
      LDXIV=LVTFLG(IVTAM)
C
C  Expand ATP parameters to include IM240 coverage.
C
      IF(LMYRIM.EQ.0) GOTO 60
      IF(ITEST(LMYRIM).NE.4) GOTO 60
      IMXSY=LPOD-ICYIM(LMYRIM)+1
      IF(LDXSY.LT.IMXSY) LDXSY=IMXSY
      IF(JDX.LT.LDX1ST.OR.LDX1ST.LT.1)  LDX1ST=JDX
      IF(JDX.GT.LDXLST.OR.LDXLST.GT.25) LDXLST=JDX
      IF(LVTFLG(IVTAM).LT.ILDT(IVTAM,LMYRIM)) LDXIV=2
C
C  Determine if this IVTAM/JDX is ATP covered.
C
   60 LMYR=1
      LMYRVT=1
      IF(LDXSY.GT.1.AND.LDX1ST.GE.JDX.AND.JDX.GE.LDXLST)  LMYR=2
      IF(LMYR.EQ.2.AND.LDXIV.EQ.2) LMYRVT=2
C
C
C  Set kinked disablement rate curves flag.
C
      IM=LMYRIM+1
      IF(IM.NE.1) IM=2
      IMKINK=1
      IF(IM.EQ.2.AND.JDX.GT.IMDXSY) IMKINK=2
C
C  Figure the base tampering rates for both exhaust and non-exhaust emissions.
C
      IY=1
      CALL DISATP(IM,IY)
C
C  Now save (only) the misfueling categories BTRs for later use by FINDPB.
C
      DO 20 IAY=1,2
      PBBTR(1,IAY,IDX,IVTAM)=BTR(3,IAY)
      PBBTR(2,IAY,IDX,IVTAM)=BTR(4,IAY)
      PBBTR(3,IAY,IDX,IVTAM)=BTR(9,IAY)
   20 CONTINUE
C
   80 CONTINUE
   90 CONTINUE
C
      RETURN
      END
