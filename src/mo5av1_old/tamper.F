C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE TAMPER(ICY,IY)
C
C  TAMPER computes tampering effects in grams/mile on bag 1, 2 & 3
C  HC, CO, & NOX and on hotsoak, diurnal & crankcase HC.
C
C  Called by EFCALX.
C
C  Calls ATPAER, BAGEME, CALAER, DISATP, EFFGRP, EMIRAT, EVPEME,
C  FINDPB, INITAP, MYRIMT
C
C  Input on call:
C
C    parameter list: ICY,IY
C    common blocks:
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP
C    /FLAGS3/ ATPFLG
C    /FLAGS4/ PRTFLG
C    /IMPAR1/ ICYIM
C    /IMPAR2/ ITEST
C    /IMPAR6/ IFREQ,INTYP
C    /IMPAR8/ ICALC
C    /MAXIMA/ MAXYRS
C    /TAMPB1/ LPOD
C    /YEARS4/ IY2020
C
C  Output on return:
C
C    common blocks:
C    /IMPAR8/ ICALC
C    /LOOKUP/ IVTAM
C    /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IMDXSY,IMKINK
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CRAIMP   R    IM240 I/M Compliance rate
C  D240DI   I    Array of disablements to use if IM240 applies
C  DNODIS   I    Array of 8*1 for initializing disablements to none
C  FQTIMP   I    IM240 I/M frequency 1=annual 2=biennial
C  ID       I    loop counting variable for disablement index
C  IM       I    flag relating run's I/M program to relative model year
C                IDX:  1 = no I/M   2 = I/M for this IDX, given ICY,
C                IVTAM, and I/M parameter
C  IMXSY    I    Local version of IMDXSY.
C  LDX1ST   I    JDX order pointer to first model year covered by ATP
C  LDXIV    I    Local ATP vehicle coverage.
C  LDXLST   I    JDX order pointer to last model year covered by ATP
C  LMYR     I    flag relating model year to ATP parameters
C                  1 = model year not covered by an ATP
C                  2 = JDX is covered by an ATP (but IVTAM may not be,
C                      hence, check for inclusion to set LMYRVT)
C  LMYRIM   I    flag relating model year, vehicle and calendar year to
C                I/M program(s). 1=I/M program 1 and 2=I/M program 2
C                0=no I/M program for MYR,CY,IVTAM
C  PGMIMP   I    IM240 I/M program type 1=test only 2=test and repair
C
C
C  Notes:
C
C  For MOBILE4.1 TAMPER was changed to call FINDPB
C  For MOBILE5 TAMPER was modified to allow checking for IM240
C  catalyst removal.
C  March 9, 1993 IY was added to the parameter list to guarantee
C  the 1990 inventory emission factors remain the same in newer
C  versions of MOBILE as in MOBILE5.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG,D240DI(8),DNODIS(8)
      INTEGER ATPPGM,ATPFQT,PGMIMP,FQTIMP,DISTYP
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
C
      COMMON /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(4)
      COMMON /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP(8)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /IMPAR8/ ICALC
      COMMON /IM240P/ DSIZE(25,4)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMEQ4/ BTR(9,2)
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,4)
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C     Putting a 3 in CAN, PCV, or CAP has no effect unless ATPFLG=2
C     Also IM240 effects CAT, ECK, and TNK.
C
      DATA D240DI/1,3,3,3,1,1,1,1/
      DATA DNODIS/8*1/
C
C  Relate start year of anti-tampering program (ATP) and 1st and last model
C  years covered to calendar year.  Put them in JDX order.  For example:
C
C      ICY = LAPSY => LDXSY = 1
C      ICY = LAPSY + 1 => LDXSY = 2
C          .
C          .
C      ICY = LAPSY + 24 => LDXSY = 25
C
C  JDX order is used because these values are used to access the
C  JDX ordered CUMMIL.  Note that LDXSY <= 0 => ICY < LAPSY => no
C  ATP was in effect during any part of the calendar year's 25
C  model year window.
C
C  The ATP has to be in effect at least 1 year prior to ICY for it
C  to affect ICY's emissions.  The mileage accumulated through
C  LAPSY all goes to the no ICY's emissions.  The mileage
C  accumulated through LAPSY all goes to the no program case.
C
C  Only calculate these values once per scenario => compute them
C  here before the tamper loops begin.
C
      LDXSY=ICY-LAPSY+1
      LDX1ST=ICY-LAP1ST+1
      LDXLST=ICY-LAPLST+1
C
C  For each tampering vehicle type and each model year in the 24 year window
C  from ICY-23 thru ICY (25+ year case is always zeroed), compute the bag and
C  evaporative tampering offsets.
C
      ICALC = 1
      DO 90 IVTAM=1,4
      DO 80 IDX=2,MAXYRS
C
C  MY pointer info using IDX & ICY may be done here, before sub calls.
C  Identify this loop's model year for use in ITAMPT lookups.
C
      MY=ICY-(MAXYRS-IDX)
      LMYRIM=MYRIMT(MY,IVTAM,ICY)
C
C  Relate I/M start year to calendar year, using JDX ordering.  This pointer
C  will be used within IDX loop to set kinked tampering rate curves flag.
C  Is this convoluteted you bet!
C
      IMDXSY=ICY-IY2020+1
      IF(LMYRIM.NE.0) IMDXSY=ICY-ICYIM(LMYRIM)+1
C
      IF(LMYRIM.EQ.0) THEN
          IF(ICALC.EQ.1) GOTO 50
          ICALC=1
          IF(ATPFLG.EQ.2) CALL ATPAER(1,ATPPGM,ATPFQT,CRATP,DISTYP)
          IF(ATPFLG.EQ.1) CALL ATPAER(1,1,1,0.0,DNODIS)
          GOTO 40
      ENDIF
      IF(ITEST(LMYRIM).EQ.4) THEN
          IF(ICALC.EQ.2) GOTO 50
          ICALC=2
          PGMIMP=INTYP(LMYRIM)
          FQTIMP=IFREQ(LMYRIM)
          CRAIMP=CRIM(LMYRIM)
          IF(ATPFLG.EQ.1) GOTO 30
          DO 20 ID=1,8
            IF(D240DI(ID).EQ.1.AND.DISTYP(ID).EQ.2) D240DI(ID)=2
   20     CONTINUE
   30     CALL ATPAER(2,PGMIMP,FQTIMP,CRAIMP,D240DI)
          GOTO 40
      ENDIF
          IF(ICALC.EQ.1) GOTO 50
          ICALC=1
          IF(ATPFLG.EQ.2) CALL ATPAER(1,ATPPGM,ATPFQT,CRATP,DISTYP)
          IF(ATPFLG.EQ.1) CALL ATPAER(1,1,1,0.0,DNODIS)
   40 CALL INITAP
      CALL CALAER
C
C  Assign the reverse order (ICY -> 1, increases with age) my index for use in
C  apportioning mileage weights.
C
   50 JDX=MAXYRS-IDX+1
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
      LDXSY=ICY-LAPSY+1
      LDX1ST=ICY-LAP1ST+1
      LDXLST=ICY-LAPLST+1
      LDXIV=LVTFLG(IVTAM)
C
C  Expand ATP parameters to include IM240 coverage
C
      IF(LMYRIM.EQ.0) GOTO 60
      IF(ITEST(LMYRIM).NE.4) GOTO 60
      IMXSY=ICY-ICYIM(LMYRIM)+1
      IF(LDXSY.LT.IMXSY) LDXSY=IMXSY
      IF(JDX.LT.LDX1ST.OR.LDX1ST.LT.1) LDX1ST=JDX
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
C  Set kinked disablement rate curves flag.
C
   70 IM=LMYRIM+1
      IF(IM.NE.1) IM=2
      IMKINK=1
      IF(IM.EQ.2.AND.JDX.GT.IMDXSY) IMKINK=2
C
C  Figure the base tampering rates for both exhaust and
C  non-exhaust emissions.
C
      CALL DISATP(IM,IY)
C
C  If ICY is beyond LPOD, use the tampering rates as of LPOD for
C  the misfueling categories (ID=3,4,9) for vehicles in ICY's
C  fleet produced thru LPOD.
C
      IF(ICY.GT.LPOD) CALL FINDPB(ICY)
C
C  The next 3 subroutines compute the tampering bag emission
C  additives.
C
      CALL EFFGRP(IY)
      CALL EMIRAT
      CALL BAGEME
C
C  The non-exhaust HC emission additives calculation procedure is
C  much simpler because there are no overlaps.  EVPEME does all
C  the work.
C
      IF(PRTFLG.EQ.1.OR.PRTFLG.EQ.4) CALL EVPEME(IY)
C
   80 CONTINUE
   90 CONTINUE
C
      RETURN
      END
