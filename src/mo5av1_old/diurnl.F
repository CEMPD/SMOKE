C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION DIURNL(MY,IDX,IVGAS,IY)
C
C  DIURNL selects the yearly value for diurnal HC.
C  DIURNL selects the finj dependent failed offset.
C
C  Called by CCEVRT.
C   4087
C  Calls EVMAIN, ETPCAL.
C
C  Input on call:
C
C    parameter list: MY,IDX,IVGAS
C    common blocks:
C    /CITCIN/ UDI,IUDI,PFRATE,IPFCU
C    /CITRV1/ RVPICY
C    /EVAPAR/ ISTDC,ISTDT,ISTDP,FINJ
C    /FLAGS3/ TEMFLG
C    /FLAGS4/ HCFLAG
C    /MAXIMA/ MAXYRS
C    /TAMOUT/ TDU
C    /TEMPS/  AMBT,TEMMIN,TEMMAX,TEMEVP
C
C  Output on return:
C
C    function: DIURNL
C    /EVAPGR/ EVP,PARD
C    /TECHFL/ CONSTF
C
C  Local array subscripts
C
C  AGWT (5)   - AGWT (IDU)
C  CARRFL (3) - CARRFL (IFDS)
C  DU (5)     - DU (IDU)
C
C  Local variable / array dictionary
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CARRFL   R    Failed offset to be indexed on finj type and used
C                inside EVUNC to calculate the failed (diurnal)
C  CFMDU    R    multiple diurnal correction factor (dep. on RVPICY)
C  AGEWT    R    age dependent diurnal component weighting factors
C  DU       R    diurnal component (1=full, 2=partial: 8AM-11AM,
C                3=partial: 10AM-3PM, 4=partial: 8AM-2PM, 5=multiple)
C  DUPF     R    the EVMAIN diurnal loss for Pass, Failed Purge, Failed
C                Pressure and controlled and uncontrolled
C  DUPASS   R    Pass emission rate to insure FPressure and FPurge
C                rates never fall below the pass rate
C  DUPURG   R    Saves purge failure emission rate for MDU check
C  DUUNC    R    diurnal loss for uncontrolled (failed purge) vehicles
C  FDU      R    weighted full diurnal loss
C  IDU      I    index for DU indicating full, 8AM-11AM, 10AM-3PM,
C                8AM-2PM, and multiple
C  IVTEMP   I    IVGAS
C  PDU      R    weighted partial diurnal loss
C  MDU      R    weighted multiple diurnal loss
C  RSTLS    R    incremental resting loss HC rate
C
C
C  Notes:
C
C  Certification fuels support from MOBILE4 was removed in MOBILE4.1.
C  Modified for MOBILE4.1 to add support for pressure/purge fail.
C  DIURNL was modified for MOBILE5v13 to include the changes to
C  reflect the differences between a 1 hour and an 8 hour measurement
C  for diurnal emissions.  Added the local array for referencing
C  the constant for diurnal fail for passing through common.
C
C
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /CITRV1/ RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX
      COMMON /EVAPGR/ EVP(4),GREVP(4,9),VGREVP(4),PARD(3),GRPD(3,9)
      COMMON /EVAPAR/ IFDS,ISTDC,ISTDT,ISTDP,FINJ(2)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TECHFL/ CONSTF
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION DU(5),AGEWT(5),CARRFL(3)
      REAL MDU
      DATA CARRFL /8.193,5.6626,3.6121/
      ICY=MY-IDX+MAXYRS
C
C  Diurnal is zero if any evaporative temperature is 40F or less,
C  or if minimum temperature is 25F or less.
C
      DIURNL=0.0
      EVP(2)=0.0
      EVP(3)=0.0
      IF(TEMMAX   .LE.40.0 .OR.
     *   AMBT     .LE.40.0.AND.TEMFLG.EQ.2 .OR.
     *   TEMEVP(1).LE.40.0 .OR.
     *   TEMEVP(2).LE.40.0 .OR.
     *   TEMMIN   .LE.25.0) RETURN
C
C  If the max and min temperatures are the same, warn the user that
C  the diurnals (partial and full) are zero.
C
      IF(TEMMAX-TEMMIN.LT.1.) THEN
        IF(IDX.EQ.1.AND.IVGAS.EQ.1) CALL QUITER(0.,0,146,INERR)
        RETURN
      ENDIF
C
      KEYEVP=2
C
C  Grab resting loss rate first
C
      IVTEMP=IVGAS
      IF(IVGAS.EQ.5) IVTEMP=8
      RSTLS=RSTLOS(MY,IVTEMP,TEMEVP(1))
C
C  Compute diurnals for full, partial 8AM-11AM, partial 10AM-3PM, and
C  partial 8AM-2PM conditions.
C  Compute the multiple diurnals for CARB, TBI, and PFI and weight by
C  sales fractions. Multiple diurnals use the in-use adjustment factor.
C  HDGV and pre-1981 LDGV/T do not have fuel injection. MC do not
C  have multiple diurnals.
C
C  Determine Multiple Diurnal Correction Factor for pass vehicles
C
      CFMDU=EXP(4.8491-0.33424*RVPICY)
C
      DU(5)=0.0
C
C  Initialize Full Diurnal for Evap. Test Procedure adjustment.
C
      DUFULL=0.0
C
      DO 40 IDU=1,4
C
      DU(IDU)=0.0
C
C  HDGV and MC only have full diurnals, so do not calculate partials.
C
      IF(IVGAS.GE.4.AND.IDU.GE.2) GOTO 40
C
C  Set UDI index and compute RIVAL.
C
      IUDI=IDU+1
      RIVAL=UDI(IUDI)/UDI(1)
C
C  All HDGV and pre-1981 LDGV/T are 100% carbureted, but are affected by
C  pass/fail.  For these vehicle classes, figure a carbureted number
C  accordingly and then skip the fuel injected (TBI and PBI) rate cases.
C  Motorcycles are 100% carbureted and not affected by pass/fail.
C
C  Loop and combine for Passed, Failed Purge and Failed Pressure
C  and for carbureted and fuel injected vehicles.
C  Limit RIVAL for full diurnals, fuel-injected vehicles only.
C
      DUPASS=0.0
C
      MAXKEQ=3
      IF(IVGAS.EQ.5) MAXKEQ=1
C
      DO 30 KEYEQ=1,MAXKEQ
C
C  Loop and combine for Passed, Failed Purge and Failed Pressure
C
C     IFDS = 1 : Carbureted
C            2 : Throttle-body fuel injection
C            3 : Ported or Multi-point fuel injection
C
      MAXFDS=3
      IF(IVGAS.GE.4.OR.MY.LT.1981) MAXFDS=1
C
      DO 20 IFDS=1,MAXFDS
C
      IF(IFDS.EQ.1) ISTUSE=ISTDC
      IF(IFDS.EQ.2) ISTUSE=ISTDT
      IF(IFDS.EQ.3) ISTUSE=ISTDP
C
      CONSTF=CARRFL(IFDS)
      IF(IDU.EQ.1.AND.IFDS.GT.1.AND.RIVAL.LT.0.5) RIVAL=0.5
C
      CALL EVMAIN(MY,KEYEVP,KEYEQ,IVGAS,ISTUSE,RIVAL,DUPF,IY,ICY)
      IF(KEYEQ.EQ.1.AND.MY.GE.1981.AND.IVGAS.EQ.1) THEN
         IF(IFDS.EQ.1) DUPF=DUPF*2.3297
         IF(IFDS.EQ.2) DUPF=DUPF*2.0301
         IF(IFDS.EQ.3) DUPF=DUPF*1.8502
      ENDIF
      IF(KEYEQ.EQ.1.AND.MY.GE.1981.AND.
     *  (IVGAS.EQ.2.OR.IVGAS.EQ.3)) THEN
         IF(IFDS.EQ.1) DUPF=DUPF*2.3297
         IF(IFDS.EQ.2) DUPF=DUPF*1.9246
         IF(IFDS.EQ.3) DUPF=DUPF*1.9246
      ENDIF
C
C  Subtract off the resting loss from each diurnal.  Limit 0.<DU.
C
      DUPF=AMAX1(DUPF-RSTLS,0.0)
C
C  Don't allow purge or pressure rates to fall below pass rates.
C
      IF(KEYEQ.GT.1.AND.DUPF.LT.DUPASS) DUPF=DUPASS
C
C  Save the purge rate for future multiple check.
C
      IF(KEYEQ.EQ.2) DUPURG=DUPF
C
C  All HDGV and pre-1981 LDGV/T and MC are 100% carbureted
C  so skip the fuel injection weighting.
C
      IF(MAXFDS.EQ.1) GO TO 10
C
C  Weight carbureted and fuel injected rates by corresponding sales fractions
C  and add results together to obtain the fleet diurnal rate.
C
      IF(IFDS.EQ.1) DUPF=(1.0-FINJ(1)-FINJ(2))*DUPF
      IF(IFDS.GT.1) DUPF=DUPF*FINJ(IFDS-1)
C
C  Save the composite pass rate for future purge/pressure checks.
C
   10 IF(KEYEQ.EQ.1) DUPASS=DUPASS+DUPF
C
C  Apply Pass, Failed Purge, and Failed Pressure weighting.
C
      DUPF=DUPF*PFRATE(KEYEQ,IPFCU)
C
C  Evap. Test Procedure adjustment.
C
      IF(IDU.EQ.1.AND.IVGAS.LE.4)
     *     DUFULL=DUFULL+DUPF*(1.0-ETPCAL(2,MY,KEYEQ,0.,0.))
C
      IF(IDU.GT.1.AND.IVGAS.LE.4)
     *     DUPF=DUPF*(1.0-ETPCAL(5,MY,KEYEQ,0.,0.))
C
C  Sum across technologies to form a fleet value.
C
      DU(IDU)=DU(IDU)+DUPF
C
   20 CONTINUE
C
C  On full diurnal iteration, handle multiple diurnal.
C  Pass portion gets the correction factor, on fail purge iteration,
C  check that pass portion is no greater than fail purge portion.
C
      IF(IDU.EQ.1) THEN
        IF(KEYEQ.EQ.1) DU(5)=DU(1)*CFMDU
        IF(KEYEQ.EQ.2.AND.DU(5).GT.DUPURG)
     *    DU(5)=DUPURG
C
C  Evap. Test Procedure adjustment to weighted multiples.
C  (Since value of pass rate is not known until purge is calculated,
C  must wait until KEYEQ=2 before applying pass ETP.)
C
        IF(KEYEQ.EQ.2.AND.IVGAS.LE.4)
     *     DU(5)=DU(5)*(1.0-((0.078+0.037)*ETPCAL(3,MY,1,0.,0.)
     *              +(0.021+0.014+0.008+0.003)*ETPCAL(4,MY,1,0.,0.))
     *                    /(0.078+0.037+0.021+0.014+0.008+0.003))
C
        IF(KEYEQ.GT.1) THEN
          IF(IVGAS.LE.4)
     *      DUPF=DUPF*(1.0-((0.078+0.037)*ETPCAL(3,MY,KEYEQ,0.,0.)
     *              +(0.021+0.014+0.008+0.003)*ETPCAL(4,MY,KEYEQ,0.,0.))
     *                   /(0.078+0.037+0.021+0.014+0.008+0.003))
          DU(5)=DU(5)+DUPF
        ENDIF
      ENDIF
C
   30 CONTINUE
C
C  Evap. Test Procedure adjustment.
C
      IF(IDU.EQ.1.AND.IVGAS.LE.4) DU(IDU)=DUFULL
C
   40 CONTINUE
C
C  Store the In-Use diurnal rate
C
      EVP(2)=DU(1)
C
C  Store the unmodified multiple diurnal rate
C
      IF(DU(5).LT.DU(1)) DU(5)=DU(1)
      EVP(3)=DU(5)
C
C  Compute the age weighted full, partial, and multiple diurnal emission
C  components. Weighting factors are vehicle dependent.
C  If HCFLAG=6, use average weighting factors instead of age
C  dependent ones for LDGV/T.
C
      GOTO (50,50,50,60,70),IVGAS
C
C  LDGV/T, with age
C
   50 IF(HCFLAG.NE.6) THEN
C
      AGE=FLOAT(MAXYRS+1-IDX)
      AGEWT(1)=(25.58-0.5*AGE+6.31+0.84*AGE)/100.
      AGEWT(2)=(37.0-0.73*AGE)/100.
      AGEWT(3)=(8.11-0.16*AGE)/100.
      AGEWT(4)=(4.22-0.08*AGE)/100.
      AGEWT(5)=(4.45-0.09*AGE+6.62+0.88*AGE)/100.
C
      FDU=DU(1)*AGEWT(1)
      PDU=DU(2)*AGEWT(2)
     *   +DU(3)*AGEWT(3)
     *   +DU(4)*AGEWT(4)
      MDU=DU(5)*AGEWT(5)
      GOTO 80
      ENDIF
C LDGV/T, average
      IF(HCFLAG.EQ.6) THEN
      FDU=DU(1)*0.3396
      PDU=DU(2)*0.3236+DU(3)*0.0709+DU(4)*0.0369
      MDU=DU(5)*0.161
      GOTO 80
      ENDIF
C HDGV
   60 FDU=DU(1)*0.498
      PDU=0.0
      MDU=(DU(5)*0.5*0.137)+(DU(5)*0.094)
      GOTO 80
C MC
   70 FDU=DU(1)
      PDU=0.0
      MDU=0.0
C
C  Combine the emission components to obtain the total diurnal value.
C
   80 DIURNL=FDU+PDU+MDU
C
C  Store the partial diurnals
C
      DO 90 IDU=2,4
      IPARTL=IDU-1
      PARD(IPARTL)=DU(IDU)
   90 CONTINUE
C
C  For HCFLAG=3 or 5, replace the full diurnal with weighted full
C  plus partials. Since HDGV and MC don't have partials, skip.
C
      IF((HCFLAG.EQ.3.OR.HCFLAG.EQ.5).AND.IVGAS.LE.3) THEN
C
      AGESUM=0.0
      DO 100 IDU=1,4
      AGESUM=AGESUM+AGEWT(IDU)
  100 CONTINUE
C
      EVP(2)=0.0
      DO 110 IDU=1,4
      EVP(2)=EVP(2)+DU(IDU)*AGEWT(IDU)/AGESUM
  110 CONTINUE
      ENDIF
C
      RETURN
      END
