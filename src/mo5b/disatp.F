C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE DISATP(IM)
C
C  DISATP computes the disablement rates for the given vehicle
C  type, model year and mileage distribution.  If the start year
C  of an anti-tampering program precedes the calendar year, then
C  2 sets of rates are produced.  The first set covers through
C  the ATP start year (the "previous" interval) and the second
C  set both "previous" and "subsequent" (to start year)
C  intervals.
C
C  Called by TAMPER, and SAVEPB.
C
C  Calls DISIMP.
C
C  Input on call:
C
C    parameters list: IM
C    common blocks:
C    /MYCODE/ JDX,LDXSY,LMYRVT,IMDXSY,IMKINK
C
C  Output on return:
C
C    common blocks:
C    /MYCODE/ IAY
C    /TAMEQ4/ BTR
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C
C
C  Notes:
C
C  For MOBILE5 DISATP was modified to pass IM and eliminate the
C  common IMPAR3.
C  March 9, 1993 IY was added to the parameter list to guarantee
C  the 1990 inventory emission factors remain the same in newer
C  versions of MOBILE as in MOBILE5.  And ICY is needed.
C  Sep-19-1994 @ CSC - ked (for bsg) Request 438 Fix 1990 and earlier calendar
C  year removed IY from the parameter list, and the calculation of ICY is no 
C  longer needed and thus removed IY and ICY from the calls to DISIMP
C
C
      INCLUDE 'MYCODE.I'
      INCLUDE 'TAMEQ4.I'
C
C  Setup section.
C
      DO 10 IAY=1,2
      DO 10 ID=1,9
      BTR(ID,IAY)=0.0
   10 CONTINUE
C
C  Calculate mileage accumulated up through start year of anti-tampering
C  program (ATP) and mileage since.  Then compute the base tampering rates.
C  There are 3 cases:
C
C  (1) There is no ATP or the ATP is not in effect until the calendar year
C      or later or the model year (JDX) being checked on this call is not in
C      the range included in the ATP or the VT is not covered by the ATP:
C      Check if ATP-applies switch is on.  Off => all mileage is "previous"
C      (to any ATP) => IAY = 1 for all tampering.
C
      IF(LMYRVT.EQ.1) CALL DISIMP(JDX,1,IM,IMKINK)
C
C  (2) There is an ATP and it starts before the model year and JDX is in the
C      ATP affected MY range.  All mileage is "subsequent" (to LAPSY) =>
C      IAY = 2 for all tampering.
C
      IF(LMYRVT.EQ.2.AND.LDXSY.GT.JDX)
     * CALL DISIMP(JDX,2,IM,IMKINK)
C
C  (3) There is an ATP and it starts before ICY but at or after my JDX and
C      JDX is in the ATP affected MY range.  This case => split mileage into
C      "previous" and overall (= total = "previous" + "subsequent") parts.
C      "Subsequent" category sizes will be computed in subroutine EFFGRP as
C      the difference between the corrected overall and "previous" sizes.
C
      IF(LMYRVT.EQ.1.OR.(LMYRVT.EQ.2.AND.LDXSY.GT.JDX)) GOTO 20
C
C  Calculate the "previous" segment.
C
      LJDX=JDX-LDXSY+1
      LKINK=1
      IF(IMKINK.EQ.1) LIM=IM
      IF(IMKINK.EQ.2.AND.IMDXSY.GT.LDXSY) LKINK=2
      IF(IMKINK.EQ.2) LIM=LKINK
      CALL DISIMP(LJDX,1,LIM,LKINK)
C
C  Calculate the overall segment.
C
      CALL DISIMP(JDX,2,IM,IMKINK)
C
C  Range checks constraining all cases:
C
C  (1) Cannot tamper with < 0% or > 100% of the fleet.
C
   20 DO 40 ID=1,9
      DO 30 IAY=1,2
      IF(BTR(ID,IAY).LT.0.0) BTR(ID,IAY)=0.0
      IF(BTR(ID,IAY).GT.1.0) BTR(ID,IAY)=1.0
   30 CONTINUE
   40 CONTINUE
C
C  (2) Note that (1) insures overall misfueling rate (ID=9) is <= 100%.
C      Now make the sum of the 2 components of misfueling be = the overall
C      misfueling, both insuring the sum is <= 100% and preventing a kink
C      up in the slope of the sum's curve, if either the NCKS or TNKS misfueling
C      curve has been forced up to zero by (1).  The tradeoff in the latter
C      case is a kink down in the slope of the fuel inlet (NCKS) or other (TNKS)
C      that was not forced up to zero.  Finally, the 3 constraint checks'
C      conditions are mutually exclusive, so the check order does not matter.
C      The 3 checks are the only use of the overall misfueling rates (ID=9).
C      The overall misfueling rates are not used hereafter.
C
      DO 50 IAY=1,2
      IF((BTR(3,IAY)+BTR(4,IAY)).GT.1.0) BTR(4,IAY)=1.0-BTR(3,IAY)
      IF(BTR(3,IAY).GT.BTR(9,IAY)) BTR(3,IAY)=BTR(9,IAY)
      IF(BTR(4,IAY).GT.BTR(9,IAY)) BTR(4,IAY)=BTR(9,IAY)
   50 CONTINUE
C
C  Check for negative growth
C
      IF (LMYRVT.EQ.1.OR.(LMYRVT.EQ.2.AND.LDXSY.GT.JDX)) GOTO 56
      DO 55 ID=1,9
        IF(BTR(ID,2)-BTR(ID,1) .LT. 0.0) BTR(ID,1)=BTR(ID,2)
   55 CONTINUE
C
C  (3) Make sure that Canister Only (ID=6) is <= Canister+Gas Cap (ID=8)
C
   56 IF(BTR(6,1).GT.BTR(8,1)) BTR(6,1)=BTR(8,1)
      IF(BTR(6,2).GT.BTR(8,2)) BTR(6,2)=BTR(8,2)
C
C  Calculate subsequent rates for exhaust tampering.
C
      IF(LMYRVT.NE.1) THEN
          DO 60 ID=1,5
            BTR(ID,2)=BTR(ID,2)-BTR(ID,1)
   60     CONTINUE
      ENDIF
      RETURN
      END
