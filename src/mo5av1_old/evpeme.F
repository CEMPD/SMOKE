C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE EVPEME(IY)
C
C  EVPEME serves the same functions for evaporative (evap = hot soak + diurnal)
C  and crankcase (cc) HC as EFFGRP, EMIRAT and BAGEME do for exhaust IP,
C  producing hs & du HC tampering rates & cc HC offset.  The task is simplified
C  into one subroutine because there is no overlap categorization or
C  weighting.  The appropriate numbers are looked up and multiplied
C  together with only anti-tampering effectiveness reductions (if any) applied.
C
C  Called by TAMPER.
C
C  Calls ITAMPT,OTCALC.
C
C  Input on call:
C
C    parameters list: IY
C    common blocks:
C    /EVPDAT/ EVPAER,EVPTGS,CCEMI
C    /LOOKUP/ IVTAM
C    /MYCODE/ IDX,JDX,LDXSY,LMYRVT
C    /TAMEQ4/ BTR
C
C  Output on return:
C
C    common blocks:
C    /IM240P/ DSIZE
C    /LOOKUP/ IHG
C    /TAMOUT/ THS,TDU,TCC
C    /TAMOU2/ TOB,GCONLY
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ICY      I    See parameter dictionary
C  OTR      R    overall (previous + subsequent my ranges) tampering rate
C  OTRCAR   R    carbureted hot soak affecting overall tampering rate
C  OTRINJ   R    fuel injected hot soak affecting overall tampering rate
C  PCTEQP   R    percentage of IVTAM fleet equipped with IHG technology,
C                given IG3
C
C  Notes:
C
C  For MOBILE4.1 EVPEME was modified to store diurnal and
C  crankcase tampering rates in DSIZE.
C  March 9, 1993 IY was added to the parameter list to guarantee
C  the 1990 inventory emission factors remain the same in newer
C  versions of MOBILE as in MOBILE5.  And ICY is needed.
C
C
      COMMON /EVPDAT/ EVPAER(3,2),EVPTGS(3,2,4),CCEMI(7,4)
      COMMON /IM240P/ DSIZE(25,4)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMEQ4/ BTR(9,2)
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TAMOU2/ TOB(25,4),GCONLY(25,4)
C
C  For MOBILE4, crankcase (IH = 3) processing does not change.
C  Diurnal processing only changes in its selection of the new evaporative
C  cannister + gas cap tampering (ID = 8) instead of the previously used
C  evaporative cannister only (ID = 6).  Since aside from data point selection,
C  the procedures for calculating the diurnal and crankcase offsets are the
C  same, a loop structure remains the convenient way to program the procedure
C  Hot soak, on the other hand, has become more complicated and has therefore
C  been pulled from the old IH loop and treated separately, after the diurnal
C  and crankcase offsets have been computed.
C
C  MOBILE4 defers HS and DU tampering offset calculation to CCEVRT,
C  returning only the % equipped * overall tampering rates instead.
C
C  Calculate ICY
C
      ICY=JDX-1+MY
C
C  Store the non-ATP rates, but first check ICY and IY to
C  guarantee the 1990 inventory emission factors remain the same
C  in newer versions of MOBILE as in MOBILE5.
C
      IF(ICY.GE.1992) THEN
        DSIZE(IDX,IVTAM)=BTR(8,2)
        IF(LMYRVT.EQ.1) DSIZE(IDX,IVTAM)=BTR(8,1)
      ENDIF
      IF(ICY.EQ.1991.AND.IY.EQ.1) THEN
        DSIZE(IDX,IVTAM)=BTR(8,2)
        IF(LMYRVT.EQ.1) DSIZE(IDX,IVTAM)=BTR(8,1)
      ENDIF
C
      IF(ICY.EQ.1991.AND.IY.EQ.2)
     *DSIZE(IDX,IVTAM)=BTR(8,1)+BTR(8,2)
      IF(ICY.LE.1990) DSIZE(IDX,IVTAM)=BTR(8,1)+BTR(8,2)
C
      DO 10 IH=2,3
C
C  There are only two evaporative cases, not three, when indexing the
C  anti-tampering effectiveness rates, evaporative technology group sizes,
C  and the model year groups for the group sizes and the evaporative emission
C  impact rates:  for these arrays, the hot soak and diurnal numbers are
C  identical.  For convenience we use IHG instead of IH in these arrays.
C
      IHG=IH-1
C
C  Assign the tampering index (see paragraph preceding DO statement).
C
      ID=9-IHG
C
C  Get the overall tampering rate.  As in DISATP, there are 3 cases:
C
C    (1) No ATP or its start year is after ICY or IDX is not in affected my
C        range: use the "previous" base rate without applying the
C        effectiveness rates.
C
      IF(LMYRVT.EQ.1) OTR=BTR(ID,1)
C
C    (2) There is an ATP and its start year is before IDX and IDX is in the
C        affected my range: use the "subsequent" rate * the effectiveness
C        rate.
C
      IF(LMYRVT.EQ.2.AND.LDXSY.GT.JDX.AND.ID.EQ.7)
     *  OTR=OTCALC(ID,2,IHG,2)
      IF(LMYRVT.EQ.2.AND.LDXSY.GT.JDX.AND.ID.EQ.8)
     *  OTR=OTCALC(ID-2,2,IHG,2)+OTCALC(ID,2,3,2)-OTCALC(ID-2,2,3,2)
C
C    (3) There is an ATP and its start year is at or after IDX but before ICY
C        and IDX is in the affected my range: sum the rates * their
C        corresponding effectiveness rates.  Note that in this case, BTR holds
C        the "previous" and overall rates => subtract 1st from 2nd to get the
C        "subsequent" rates.
C
      IF(LMYRVT.EQ.2.AND.LDXSY.LE.JDX.AND.ID.EQ.7)
     *  OTR=OTCALC(ID,1,IHG,1)+OTCALC(ID,2,IHG,2)-OTCALC(ID,1,IHG,2)
      IF(LMYRVT.EQ.2.AND.LDXSY.LE.JDX.AND.ID.EQ.8)
     *  OTR=OTCALC(ID-2,1,IHG,1)+OTCALC(ID,1,3,1)-OTCALC(ID-2,1,3,1)
     *     +OTCALC(ID-2,2,IHG,2)+OTCALC(ID,2,3,2)-OTCALC(ID-2,2,3,2)
     *    -(OTCALC(ID-2,1,IHG,2)+OTCALC(ID,1,3,2)-OTCALC(ID-2,1,3,2))
C
C  Lookup model year group pointers.
C
      IG3=ITAMPT(3)
C
      IF(IH.EQ.3) IG4=ITAMPT(4)
C
C  Lookup and apply the percentage equipped and emission impact factors.
C
      PCTEQP=EVPTGS(IG3,IHG,IVTAM)
      IF(IH.EQ.2) TDU(IDX,IVTAM)=OTR*PCTEQP
      IF(IH.EQ.3) TCC(IDX,IVTAM)=OTR*PCTEQP*CCEMI(IG4,IVTAM)
   10 CONTINUE
C
C  Now calculate the hot soak offset.  This case is complicated by both the
C  usage of 2 tampering types (evaporative cannister only and evaporative
C  cannister plus gas caps) and the weighting together of fuel injected and
C  carbureted emission impact rates.
C
      IH=1
      IHG=IH
C
C  Assign the carbureted tampering index = evaporative cannister only, since
C  carbureted hot soak HC is not affected by gas cap tampering.  The fuel
C  injected index is ID+2 = 8 = cannister + gas caps.
C
      ID=6
C
C  Get the overall tampering rates.  As with diurnal and crankcase, there are
C  3 cases, but now 2 rates are computed.
C
C    (1) No ATP or its start year is after ICY or IDX is not in affected my
C        range: use the "previous" base rate without applying the
C        effectiveness rates.
C
      IF(LMYRVT.NE.1) GOTO 20
      OTRCAR=BTR(ID,1)
      OTRINJ=BTR(ID+2,1)
      GOTO 40
C
C    (2) There is an ATP and its start year is before IDX and IDX is in the
C        affected my range: use the "subsequent" rate * the effectiveness
C        rate.
C
   20 IF(LDXSY.LE.JDX) GOTO 30
      OTRCAR=OTCALC(ID,2,IHG,2)
      OTRINJ=OTCALC(ID,2,IHG,2)+OTCALC(ID+2,2,3,2)-OTCALC(ID,2,3,2)
      GOTO 40
C
C    (3) There is an ATP and its start year is at or after IDX but before ICY
C        and IDX is in the affected my range: sum the rates * their
C        corresponding effectiveness rates.  Note that in this case, BTR holds
C        the "previous" and overall rates => subtract 1st from 2nd to get the
C        "subsequent" rates.
C
   30 OTRCAR=OTCALC(ID,1,IHG,1)+OTCALC(ID,2,IHG,2)-OTCALC(ID,1,IHG,2)
      OTRINJ=OTCALC(ID,1,IHG,1)+OTCALC(ID+2,1,3,1)-OTCALC(ID,1,3,1)
     *      +OTCALC(ID,2,IHG,2)+OTCALC(ID+2,2,3,2)-OTCALC(ID,2,3,2)
     *     -(OTCALC(ID,1,IHG,2)+OTCALC(ID+2,1,3,2)-OTCALC(ID,1,3,2))
C
C  Save the unweighted (by PCTEQP) Evaporative Cannister (EC) Only tampering
C  rate for figuring the On Board (OB) EC tampering impact on refueling losses.
C  Lookup model year group pointer for the percentage equipped table.
C
   40 TOB(IDX,IVTAM)=OTRCAR
      IG3=ITAMPT(3)
C
C  Lookup the percentage equipped with canisters factor and apply it to
C  the otr's to get the tampered percentage of carbureted and fuel injected
C  cases for hot soak.
C
      PCTEQP=EVPTGS(IG3,IHG,IVTAM)
      THS(1,IDX,IVTAM)=PCTEQP*OTRCAR
      THS(2,IDX,IVTAM)=PCTEQP*OTRINJ
C
C  Store the gas cap only tampering values for running loss weighting.
C
      GCONLY(IDX,IVTAM)=OTRINJ-OTRCAR
      IF(GCONLY(IDX,IVTAM).LT.0.0) GCONLY(IDX,IVTAM)=0.0
C
      RETURN
      END
