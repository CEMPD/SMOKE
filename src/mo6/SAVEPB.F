      SUBROUTINE SAVEPB 
C 
C  SAVEPB computes and saves the misfueling category BTRs for all 
C  IVTAM x IDX cases, using LPOD = the last year of leaded gasoline, not ICY. 
C  FINDPB will then use the resulting table to get the LPOD misfueling rates 
C  with which to replace the corresponding ICY rates.  The latter are in 
C  error because DISATP in using ICY assumes misfueling continues after LPOD. 
C 
C  Called by CHKINP. 
C 
C  Calls DISATP, EVPEME, FINDIM. 
C 
C  Changes: (Last change first) 
C 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I; pass MY, JDX, and IDX through parameter lists 
C  20 Sep 00 @EPA- BG Added constant IEXHAUST=1 for call to FINDIM 
C 15 Sep 00 AIR Task 08: Updated FINDIM calls to select exhaust I/M. 
C 28 Jan 99 @DynTel-MLA 2-663  Replaced calls to MYRATP and MYRIMT with 
C           calls to FINDIM. MYRATP and MYRIMT resolved overlapping I/M 
C           programs, a situation that is not allowed in Mobile6. FINDIM 
C           simply returns the number of the effective I/M program, if any. 
C 06 Jan 99 @DynTel-ZK 2-000 Explicit Typing 
C 
C  6 Nov 98 @DynTel-MLA 2-663  A new dimension was added to ICYIM to allow 
C           the array to hold program end dates as well as start dates. All 
C           references to ICYIM in existing code were updated accordingly. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    common blocks: 
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG 
C    /IMPAR1/ ICYIM,MODYR 
C    /IMPAR2/ ITEST 
C    /MAXIMA/ MAXYRS 
C    /TAMPB1/ LPOD 
C    /VVNAMS/ VVGAS  (IVTYPE.I) 
C    /YEARS4/ IYEND 
C 
C  Output on return: 
C 
C    common blocks: 
C    /LOOKUP/ IVTAM 
C    /MYCODE/ LDXSY,LMYRVT,IMDXSY,IMKINK 
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
C  LIM240   I    flag relating model year, vehicle and calendar year to 
C                I/M program(s). 1=I/M program 1 and 2=I/M program 2 
C                0=no I/M program for MYR,CY,IVTAM 
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
C  For MOBILE4.1 SAVEPB was created.  SAVEPB is 
C  identical to TAMPER, except that LPOD is used instead of ICY 
C  and the save to PBBTR replaces the calls to EEFGRP, EMIRAT, 
C  BAGEME, and EVPEME. 
C  For MOBILE5 SAVEPB was modified to include the IM240 effect on 
C  tampering.  Needed to add a call to use the MYRIMT function. 
C  March 10, 1993 IY=1 was added to the argument list to the 
C  call to DISATP to guarantee the 1990 inventory emission 
C  remain the same in newer versions of MOBILE as in MOBILE5. 
C  15-April-1993 @ ARC-bsg for Subtask 2-213 
C  SAVEPB was changed to modify the I/M deterrence effect 
C  This is accomplished by using the MYRATP function 
C  instead of the MYRIMT function. 
C  June 1, 1993 @ ARC-bsg Subtask 238 
C  Sep-22-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar year 
C  removed IY from the argument list 
C  16 Jan 96 @ DynTel-MLA 2-610 'SAVEPB Inconsistency'  Added new variable  
C    LIM240 and call to function MYRIMT. LIM240 is used in place of LMYRIM  
C    in some cases. This aligns SAVEPV with TAMPER - the characteristics of 
C    the second (governing) I/M program are used to determine whether ATP 
C    benefits should be extended back to start year of the first I/M prog. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. Added variable IV to 
C            map certain values of IVTAM. 
C 
      IMPLICIT NONE 
      INCLUDE 'ATPAR1.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'TAMEQ4.I' 
      INCLUDE 'TAMPB1.I' 
      INCLUDE 'YEARS4.I' 
C 
      INTEGER FINDIM 
      INTEGER LMYRIM 
      INTEGER LIM240 
      INTEGER IIM240 
      INTEGER LDXLST 
      INTEGER LDXIV 
      INTEGER LMYR 
      INTEGER IM 
      INTEGER LDX1ST 
      INTEGER IV 
      INTEGER IDX 
      INTEGER JDX 
      INTEGER MY 
      INTEGER, SAVE :: IEXHAUST=1 
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
      IF(LDXSY.GT.25) LDXSY=25 
      IF(LDXSY.LT.1) LDXSY=1 
      IF(LDX1ST.GT.25) LDX1ST=25 
      IF(LDX1ST.LT.1) LDX1ST=1 
      IF(LDXLST.GT.25) LDX1ST=25 
      IF(LDXLST.LT.1) LDXLST=1 
C 
C  For each tampering vehicle type and each model year in the 24 year window 
C  from LPOD-23 thru LPOD (25+ year case is always zeroed), compute the bag and 
C  evaporative tampering offsets. 
C 
      IVTAM = 0 
      DO 90 IV = 1,MAXVEH 
C 
        IF (VVGAS(IV).EQ.0)  CYCLE 
        IVTAM = IVTAM + 1 
C 
        DO 80 IDX=2,MAXYRS 
C 
C  MY pointer info using IDX & LPOD may be done here, before sub calls. 
C  Identify this loop's model year for use in ITAMPT lookups. 
C 
          MY=LPOD-(MAXYRS-IDX) 
          LMYRIM=FINDIM(LPOD,MY,IV,IEXHAUST) 
          LIM240=LMYRIM 
C 
C  Relate I/M start year to calendar year, using JDX ordering.  This pointer 
C  will be used within IDX loop to set kinked tampering rate curves flag. 
C 
          IMDXSY=0 
          IIM240=1 
          IF(LMYRIM.NE.0) THEN 
            IMDXSY=LPOD-ICYIM(1,LMYRIM)+1 
            IF(ITEST(LIM240).EQ.4) IIM240=2 
          ENDIF 
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
C  LMYRVT = 1 => switch is off.  ATP not a factor in my JDX's tampering 
C  LMYRVT = 2 => switch is ON.  ATP is a factor in my JDX's tampering 
C 
C  Note: ATPFLG = 1 => LAPSY = 2050 => LDXSY < 1 => 1st 2 criteria can be 
C  checked using LDXSY. 
C 
          LDXSY=LPOD-LAPSY+1 
          LDX1ST=LPOD-LAP1ST+1 
          LDXLST=LPOD-LAPLST+1 
          IF(LDXSY.GT.25) LDXSY=25 
          IF(LDXSY.LT.1) LDXSY=1 
          IF(LDX1ST.GT.25) LDX1ST=25 
          IF(LDX1ST.LT.1) LDX1ST=1 
          IF(LDXLST.GT.25) LDX1ST=25 
          IF(LDXLST.LT.1) LDXLST=1 
          LDXIV=LVTFLG(IVTAM) 
C 
C  Expand ATP parameters to include IM240 coverage 
C 
          IF(IIM240.EQ.2) THEN 
C 
C  If IM240 I/M starts before the ATP, use the IM240 start date to 
C  determine ATP effects. 
C 
            IF(LDXSY.LT.LPOD-ICYIM(1,LIM240)+1) 
     *        LDXSY=LPOD-ICYIM(1,LIM240)+1 
C 
C  If IM240 I/M model year coverage is greater than ATP MYR coverage, 
C  expand ATP coverage to match I/M program 
C 
            IF(LPOD-MODYR(1,LIM240)+1.GT.LDX1ST) 
     *        LDX1ST=LPOD-MODYR(1,LIM240)+1 
            IF(LPOD-MODYR(2,LIM240)+1.LT.LDXLST) 
     *        LDXLST=LPOD-MODYR(2,LIM240)+1 
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
C  Now save (only) the misfueling categories BTRs for later use by FINDPB. 
C 
          DO 20 IAY=1,2 
            PBBTR(1,IAY,IDX,IVTAM)=BTR(3,IAY) 
            PBBTR(2,IAY,IDX,IVTAM)=BTR(4,IAY) 
            PBBTR(3,IAY,IDX,IVTAM)=BTR(9,IAY) 
   20     CONTINUE 
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
