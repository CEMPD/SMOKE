      SUBROUTINE ATTAMPER(ICY,ITTT) 
C 
C  ATTAMPER computes Air Toxics tampering effects in grams/mile on 
C  bag 1, 2 & 3 HC and on hotsoak, diurnal & crankcase HC. 
C 
C  Called by EFCALX. 
C 
C  Calls ATPAER, BAGEME, CALAER, DISATP, EFFGRP, ATEMIRAT, EVPEME, 
C  FINDIM, FINDPB, INITAP 
C 
C  Changes: (Last change first) 
C 
C  26 Feb 02 New MOBILE6.3 routine, based on 27 Nov 00 TAMPER.FOR 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: ICY 
C 
C    common blocks: 
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG 
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP 
C    /IMPAR1/ ICYIM,MODYR 
C    /IMPAR2/ ITEST 
C    /IMPAR6/ IFREQ,INTYP 
C    /IMPAR8/ ICALC 
C    /MAXIMA/ MAXYRS 
C    /TAMPB1/ LPOD 
C    /VVNAMS/ VVGAS  (IVTYPE.I) 
C    /YEARS4/ IYEND 
C 
C  Output on return: 
C 
C    common blocks: 
C    /IMPAR8/ ICALC 
C    /LOOKUP/ IVTAM 
C    /MYCODE/ LDXSY,LMYRVT,IMDXSY,IMKINK 
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
C  ICALC    I    indicates the ATPAER functions status 
C                0 - has not been called or has been called to reset 
C                1 - ATP only using DISTYP 
C                2 - IM240 only using D240DI 
C                3 - both IM240 and ATP cover use both DISTYP and D240DI 
C  IIM240   I    flag intificating IM240 is present 1=no 2=yes 
C  IM       I    flag relating run's I/M program to relative model year 
C                IDX:  1 = no I/M   2 = I/M for this IDX, given ICY, 
C                IVTAM, and I/M parameter 
C  ITTT     I    Tambag and Tambas Toggle 
C  LATPTP   I    ATP type 1=no ATP 2=ATP without evap 3=ATP with evap 
C  LDX1ST   I    JDX order pointer to first model year covered by ATP 
C  LDXIV    I    Local ATP vehicle coverage. 
C  LDXLST   I    JDX order pointer to last model year covered by ATP 
C  LIM240   I    flag relating model year, vehicle and calendar year to 
C                I/M program(s). 1=I/M program 1 and 2=I/M program 2 
C                0=no I/M program for MYR,CY,IVTAM 
C  LMYR     I    flag relating model year to ATP parameters 
C                  1 = model year not covered by an ATP 
C                  2 = JDX is covered by an ATP (but IVTAM may not be, 
C                      hence, check for inclusion to set LMYRVT) 
C  LMYRIM   I    flag relating model year, vehicle and calendar year to 
C                I/M program(s). 1=I/M program 1 and 2=I/M program 2 
C                0=no I/M program for MYR,CY,IVTAM 
C  PGMIMP   I    IM240 I/M program type 1=test only 2=test and repair 
C  S240DI   I    Saves array of disablements to use if IM240 applies 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATPAR1.I' 
      INCLUDE 'ATPAR2.I' 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IM240P.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'TAMEQ4.I' 
      INCLUDE 'TAMPB1.I' 
      INCLUDE 'TAMOUT.I' 
      INCLUDE 'YEARS4.I' 
C 
      INTEGER FINDIM 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER FQTIMP,IV,PGMIMP 
      INTEGER CHKATP 
      INTEGER ITTT 
      INTEGER IDX 
      INTEGER JDX 
      INTEGER ACTUAL_MY 
      INTEGER MY 
      INTEGER  LDXLST, ICALC, ID, LDXIV, LMYR, IM, LATPTP, LDX1ST, 
     *         LIM240, LMYRIM, IIM240 
      INTEGER, DIMENSION(8) :: D240DI 
      INTEGER, DIMENSION(8) :: S240DI=(/1,3,3,3,1,1,1,1/) 
      INTEGER, DIMENSION(8) :: DNODIS=1 
      INTEGER, SAVE :: IEXHAUST = 1 
      REAL CRAIMP 
C 
C     Putting a 3 in CAN, PCV, or CAP has no effect unless ATPFLG=2 
C     Also IM240 effects CAT, ECK, and TNK. 
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
      IF(LDXSY.GT.25) LDXSY=25 
      IF(LDXSY.LT.1) LDXSY=1 
      IF(LDX1ST.GT.25) LDX1ST=25 
      IF(LDX1ST.LT.1) LDX1ST=1 
      IF(LDXLST.GT.25) LDX1ST=25 
      IF(LDXLST.LT.1) LDXLST=1 
C 
      IF(ITTT.EQ.2) THEN 
        LDXSY  = 1 
        LDX1ST = 1 
        LDXLST = 1 
      ENDIF 
C 
C  For each tampering vehicle type and each model year in the 24 year window 
C  from ICY-23 thru ICY (25+ year case is always zeroed), compute the bag and 
C  evaporative tampering offsets. 
C 
C  Initialize 
C 
      ICALC = 0 
      IVTAM = 0 
      DO 90 IV = 1,MAXVEH 
C 
        IF (VVGAS(IV).EQ.0)  CYCLE 
        IVTAM = IVTAM + 1 
C 
        DO 80 IDX=2,MAXYRS 
C 
C  MY pointer info using IDX & ICY may be done here, before sub calls. 
C  Identify this loop's model year for use in ITAMPT lookups. 
C 
          ACTUAL_MY=ICY+IDX-MAXYRS 
          MY=ACTUAL_MY 
          IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
            MY = MIN(CAAMY,MY) 
          ENDIF 
C 
          LIM240=FINDIM(ICY,MY,IV,IEXHAUST) 
          LATPTP=CHKATP(MY,IV,ICY) 
C 
          IF(ITTT.EQ.2) THEN 
            LIM240 = 0 
            LATPTP = 1 
          ENDIF 
C 
          LMYRIM=LIM240 
C 
C  Relate I/M start year to calendar year, using JDX ordering.  This pointer 
C  will be used within IDX loop to set kinked tampering rate curves flag. 
C 
          IMDXSY=0 
          IIM240=1 
          IF(LMYRIM.NE.0) THEN 
C 
            IMDXSY=ICY-ICYIM(1,LMYRIM)+1 
C 
C  When IM240 is the applicable test type set tampering program parameters 
C 
            IF(ITEST(LIM240).EQ.4) IIM240=2 
C 
            IF(IIM240.EQ.2) THEN 
              PGMIMP=INTYP(LIM240) 
              IF (PGMIMP.GT.2) PGMIMP = 2 
              FQTIMP=IFREQ(LIM240) 
              CRAIMP=CRIM(LIM240) 
            ENDIF 
C 
          ENDIF 
C 
C  ATP is off, IM240 is off, check if this calculation has happened already 
C 
          IF(LATPTP.EQ.1) THEN 
            IF(IIM240.EQ.1) THEN 
              IF(ICALC.EQ.0) GOTO 70 
C 
              CALL INITAP 
              CALL ATPAER(MY,IDX,1,1,1,0.0,DNODIS) 
              CALL CALAER 
C 
C  Set indicator that the calculation has happened already 
C 
              ICALC=0 
            ENDIF 
C 
C  ATP is off, IM240 is on, check if this calculation has happened already 
C 
            IF(IIM240.NE.1) THEN 
              IF(ICALC.EQ.2) GOTO 70 
              CALL INITAP 
              CALL ATPAER(MY,IDX,2,PGMIMP,FQTIMP,CRAIMP,S240DI) 
              CALL CALAER 
C 
C  Set indicator that the calculation has happened already 
C 
              ICALC=2 
            ENDIF 
          ENDIF 
C 
C  Check presence of ATP 
C 
          IF(LATPTP.NE.1) THEN 
C 
C  Check presence of IM240 
C 
            IF(IIM240.EQ.1) THEN 
C 
C  ATP is on, IM240 is off, check if this calculation has happened already 
C 
              IF(ICALC.EQ.1) GOTO 70 
              CALL INITAP 
              CALL ATPAER(MY,IDX,1,ATPPGM,ATPFQT,CRATP,DISTYP) 
              CALL CALAER 
C 
C  Set indicator that the calculation has happened already 
C 
              ICALC=1 
            ENDIF 
C 
C  ATP is on, IM240 is on, check if this calculation has happened already 
C 
            IF(IIM240.NE.1) THEN 
              IF(ICALC.EQ.3) GOTO 70 
              CALL INITAP 
              DO 10 ID=1,8 
                D240DI(ID)=S240DI(ID) 
   10         CONTINUE 
              DO 20 ID=1,8 
                IF(D240DI(ID).EQ.1.AND.DISTYP(ID).EQ.2) D240DI(ID)=2 
   20         CONTINUE 
              CALL ATPAER(MY,IDX,2,PGMIMP,FQTIMP,CRAIMP,D240DI) 
              CALL CALAER 
C 
C  Set indicator that the calculation has happened already 
C 
              ICALC=3 
            ENDIF 
C 
          ENDIF 
C 
C  Assign the reverse order (ICY -> 1, increases with age) my index for use in 
C  apportioning mileage weights. 
C 
   70 JDX=MAXYRS-IDX+1 
C 
C  Reset (turn off) the ATP my range switch.  Turn on if there is an ATP 
C  program, it starts before the calendar year, the model year JDX is in the 
C  affected range and the vehicle type has been selected by the user for 
C  coverage. 
C 
C  LMYRVT = 1 => switch is off.  ATP not a factor in my JDX's tampering 
C  LMYRVT = 2 => switch is ON.  ATP is a factor in my JDX's tampering 
C 
C  Note: ATPFLG = 1 => LAPSY = 2050 => LDXSY < 1 => 1st 2 criteria can be 
C  checked using LDXSY. 
C 
          LDXSY=ICY-LAPSY+1 
          LDX1ST=ICY-LAP1ST+1 
          LDXLST=ICY-LAPLST+1 
          IF(LDXSY.GT.25) LDXSY=25 
          IF(LDXSY.LT.1) LDXSY=1 
          IF(LDX1ST.GT.25) LDX1ST=25 
          IF(LDX1ST.LT.1) LDX1ST=1 
          IF(LDXLST.GT.25) LDX1ST=25 
          IF(LDXLST.LT.1) LDXLST=1 
          LDXIV=LVTFLG(IVTAM) 
C 
          IF(ITTT.EQ.2) THEN 
            LDXSY  = 1 
            LDX1ST = 1 
            LDXLST = 1 
          ENDIF 
C 
C  Expand ATP parameters to include IM240 coverage 
C 
          IF(IIM240.EQ.2) THEN 
C 
C  If IM240 I/M starts before the ATP, use the IM240 start date to 
C  determine ATP effects. 
C 
          IF(LDXSY.LT.ICY-ICYIM(1,LIM240)+1) 
     *      LDXSY=ICY-ICYIM(1,LIM240)+1 
C 
C  If IM240 I/M model year coverage is greater than ATP MYR coverage, 
C  expand ATP coverage to match I/M program 
C 
            IF(ICY-MODYR(1,LIM240)+1.GT.LDX1ST) 
     *        LDX1ST=ICY-MODYR(1,LIM240)+1 
            IF(ICY-MODYR(2,LIM240)+1.LT.LDXLST) 
     *        LDXLST=ICY-MODYR(2,LIM240)+1 
            IF(LVTFLG(IVTAM).LT.ILDT(IVTAM,LIM240)) LDXIV=2 
          ENDIF 
C 
C  Determine if this IVTAM/JDX is ATP covered. 
C 
          LMYR=1 
          LMYRVT=1 
          IF(LDXSY.GT.1.AND.LDX1ST.GE.JDX.AND.JDX.GE.LDXLST)  LMYR=2 
          IF(LMYR.EQ.2.AND.LDXIV.EQ.2) LMYRVT=2 
C 
C  Set kinked disablement rate curves flag. 
C 
          IM=LMYRIM+1 
          IF(IM.NE.1) IM=2 
          IMKINK=1 
          IF(IM.EQ.2.AND.JDX.GT.IMDXSY) IMKINK=2 
C 
C  Figure the base tampering rates for both exhaust and 
C  non-exhaust emissions. 
C 
          CALL DISATP(MY,JDX,IM) 
C 
C  If ICY is beyond LPOD, use the tampering rates as of LPOD for 
C  the misfueling categories (ID=3,4,9) for vehicles in ICY's 
C  fleet produced thru LPOD. 
C 
          IF(ICY.GT.LPOD) CALL FINDPB(ICY,ACTUAL_MY,IDX) 
C 
C  The next 3 subroutines compute the tampering bag emission 
C  additives. 
C 
          CALL EFFGRP(MY) 
          CALL ATEMIRAT(MY) ! AT version 
          CALL BAGEME(IDX,ITTT) 
C 
C  The non-exhaust HC emission additives calculation procedure is 
C  much simpler because there are no overlaps.  EVPEME does all 
C  the work. 
C 
          CALL EVPEME(MY,IDX,JDX) 
C 
   80   CONTINUE 
   90 CONTINUE 
C 
      RETURN 
      END 
