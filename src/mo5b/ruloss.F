C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION RULOSS(MY,IVTL,RVPHS)
C
C  RULOSS computes the Running Loss value for each model year and gas
C  vehicle class at the given RVP and running Loss temperature using
C  4-point linear interpolation.
C
C  Called by HCCALX.
C
C  Calls ETPCAL.
C
C  Input on call:
C
C    parameter list: MY,IVTL,RVPHS
C    common blocks:
C    /CITCIN/ PFRATE,IPFCU
C    /FLAGS3/ TEMFLG
C    /RUNLS1/ TLEMI,TLVMT,TLVMTU,MAXTL
C    /RUNLS2/ TLFAIL
C    /RUNLS3/ RULRVP,RULTEM,RULSPD,MAXRUL
C    /TAMOUT/ THS
C    /TEMPS/  TEMEVP
C
C  Output on return:
C
C    function: RULOSS
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  EM11     R    running loss emission (or offset) at lower boundary
C                temperature and lower boundary RVP
C  EM12     R    running loss emission (or offset) at lower boundary
C                temperature and upper boundary RVP
C  EM21     R    running loss emission (or offset) at upper boundary
C                temperature and lower boundary RVP
C  EM22     R    running loss emission (or offset) at upper boundary
C                temperature and upper boundary RVP
C  EMT1     R    interpolated running loss emission (or offset)
C                at the input RVP and the lower boundary temperature
C  EMT2     R    interpolated running loss emission (or offset)
C                at the input RVP and the upper boundary temperature
C  FACTOR   R    misc. intermediate factor for interpolations and such
C  RSTL     R    resting loss rate to be subtracted off
C  RUFAIL   R    running loss emission offset for canister pressure/
C                purge fail vehicles
C  RUPASS   R    running loss emission rate for pass vehicles
C  RVP1     R    lower boundary RVP
C  RVP2     R    upper boundary RVP
C  RVPINT   R    running loss RVP used for interpolation/extrapolation
C  RVPRL    R    running loss RVP
C  TEMRL    R    running loss temperature
C  TMP1     R    lower boundary temperature
C  TMP2     R    upper boundary temperature
C
C  Notes:
C
C  Oxygenated and cert. fuels support from MOBILE4 was removed in MOBILE4.1.
C  RULOSS was modified in Version 4.1 to activate the use of running loss
C  values by trip length rather than using trip-length weighted values
C  and to interpolate running loss emission by vehicle speed.
C  The tampering effects calculated in CAPOFF and CANOFF were replaced
C  with pressure/purge failed vehicle information in MOBILE4.1
C
C
      INCLUDE 'CITCIN.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'RUNLS1.I'
      INCLUDE 'RUNLS2.I'
      INCLUDE 'RUNLS3.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'TAMOUT.I'
      INCLUDE 'TEMPS.I'
C
C  Initialize
C
      RULOSS=0.0
      IF(TLVMTU(1).GE.0.0) GO TO 20
C
C  First initialize TLVMTU array to built-in values from TLVMT
C  unless some value was already read in (or RULOSS already called)
C
      DO 10 ITL=1,MAXTL
      TLVMTU(ITL)=TLVMT(ITL)
   10 CONTINUE
C
C  Running loss is zero if any evaporative temperature is 40F or less,
C  or if minimum temperature is 25F or less.
C
   20 IF(TEMMAX   .LE.40.0 .OR.
     *   AMBT     .LE.40.0.AND.TEMFLG.EQ.2 .OR.
     *   TEMEVP(1).LE.40.0 .OR.
     *   TEMEVP(2).LE.40.0 .OR.
     *   TEMMIN   .LE.25.0) RETURN
C
C  Choice of RVP depends on whether or not cert fuel change applies to the
C  current MY.
C
      RVPRL=RVPHS
C
C  Reset RVP and Running Loss Temperature if above running loss limits
C
      IF(RVPRL.LT.RULRVP(1)) RVPRL=RULRVP(1)
      IF(RVPRL.GT.RULRVP(4)) RVPRL=RULRVP(4)
C
      TEMRL=TEMEVP(2)
      IF(TEMRL.GT.RULTEM(4)) TEMRL=RULTEM(4)
C
C  Set the interpolation RVP.
C  (Allow for extrapolation 0.5 RVP below RULRVP lower limit.)
C
      RVPINT=RVPHS
      IF(RVPINT.LT.RULRVP(1)-0.5) RVPINT=RULRVP(1)-0.5
      IF(RVPINT.GT.RULRVP(4)) RVPINT=RULRVP(4)
C
C  Find and assign the bracketing temperature values.
C
      DO 30 ITLT=1,3
      IF(TEMRL.LE.RULTEM(ITLT)) GOTO 40
   30 CONTINUE
C
      ITLT=4
C
   40 TMP1=40.0
      IF(ITLT.GT.1) TMP1=RULTEM(ITLT-1)
      TMP2=RULTEM(ITLT)
C
C  Find and assign the bracketing RVP values.
C
      DO 50 ITLR=2,3
      IF(RVPRL.LE.RULRVP(ITLR)) GOTO 60
   50 CONTINUE
C
      ITLR=4
C
   60 RVP1=RULRVP(ITLR-1)
      RVP2=RULRVP(ITLR)
C
C  Find and assign the bracketing SPD values (TLVMT groups)
C
      SPDRL=SPD(IVTL)
      ISPD1=1
      ISPD2=1
      IF(SPDRL.GT.RULSPD(1)) ISPD2=2
      IF(SPDRL.GE.RULSPD(2)) ISPD1=2
      IF(SPDRL.GT.RULSPD(2)) ISPD2=3
      IF(SPDRL.GE.RULSPD(3)) ISPD1=3
      FACTOR=0.
      IF(ISPD1.NE.ISPD2)
     * FACTOR=(SPDRL-RULSPD(ISPD1))/(RULSPD(ISPD2)-RULSPD(ISPD1))
C
C  Grab the resting loss rate for this temperature.
C
      RSTL=RSTLOS(MY,IVTL,TEMRL)
C
C  Determine model year pointer (Pre-Controlled groups use TLFAIL)
C
      IG=IRULPT(MY,IVTL)
      IF(IG.EQ.1) GO TO 80
C
C  Reduce IG by one to use correct TLEMI components for pass IV
C
      IGTL=IG-1
C
C  Assign bracketing running loss emission values
C  Interpolate on vehicle speed as the values are added up
C  Subtract off resting losses for each iteration.
C  Convert to g/mi and take the bag fraction to subtract off.
C
      EM11=0.0
      EM12=0.0
      EM21=0.0
      EM22=0.0
      RSTL=RSTL/(MAXTL*SPDRL)
C
      DO 70 ITL=1,MAXTL
      IF(ITLT.GT.1) EM11=EM11+(TLEMI(ITL,ISPD1,ITLR-1,ITLT-1,IGTL,IVTL)
     *   + FACTOR * (TLEMI(ITL,ISPD2,ITLR-1,ITLT-1,IGTL,IVTL) -
     *          TLEMI(ITL,ISPD1,ITLR-1,ITLT-1,IGTL,IVTL)) - RSTL) *
     *               TLVMTU(ITL)
      IF(ITLT.GT.1) EM12=EM12+(TLEMI(ITL,ISPD1,ITLR,ITLT-1,IGTL,IVTL)
     *   + FACTOR * (TLEMI(ITL,ISPD2,ITLR,ITLT-1,IGTL,IVTL) -
     *          TLEMI(ITL,ISPD1,ITLR,ITLT-1,IGTL,IVTL)) - RSTL) *
     *               TLVMTU(ITL)
      EM21=EM21+(TLEMI(ITL,ISPD1,ITLR-1,ITLT,IGTL,IVTL)
     *   + FACTOR * (TLEMI(ITL,ISPD2,ITLR-1,ITLT,IGTL,IVTL) -
     *          TLEMI(ITL,ISPD1,ITLR-1,ITLT,IGTL,IVTL)) - RSTL) *
     *               TLVMTU(ITL)
      EM22=EM22+(TLEMI(ITL,ISPD1,ITLR,ITLT,IGTL,IVTL)
     *   + FACTOR * (TLEMI(ITL,ISPD2,ITLR,ITLT,IGTL,IVTL) -
     *          TLEMI(ITL,ISPD1,ITLR,ITLT,IGTL,IVTL)) - RSTL) *
     *               TLVMTU(ITL)
   70 CONTINUE
C
C  Interpolate first on RVP, then temperature
C
      EMT1=RVPINT*(EM12-EM11)/(RVP2-RVP1)
     *    +(RVP2*EM11-RVP1*EM12)/(RVP2-RVP1)
      EMT2=RVPINT*(EM22-EM21)/(RVP2-RVP1)
     *    +(RVP2*EM21-RVP1*EM22)/(RVP2-RVP1)
C
      RUPASS=TEMRL*(EMT2-EMT1)/(TMP2-TMP1)
     *      +(TMP2*EMT1-TMP1*EMT2)/(TMP2-TMP1)
C
C  Compute the running loss for canister pressure/purge fail vehicles.
C
C  Assign bracketing running loss pressure/purge fail offset values
C  Don't forget to subtract off fraction of resting loss rate (g/hr).
C
   80 EM11=0.0
      EM12=0.0
      EM21=0.0
      EM22=0.0
      RSTL=RSTL/FLOAT(MAXTL)
C
      DO 90 ITL=1,MAXTL
C
C  Figure conversion factor for g-->g/mi (plus include the %TLVMT)
C
      FACTOR=TLVMTU(ITL)*FLOAT(MAXTL)/(SPD(IVTL)*ITL)
C
      IF(ITLT.GT.1) EM11=EM11+(TLFAIL(ITL,ITLR-1,ITLT-1)-RSTL)*FACTOR
      IF(ITLT.GT.1) EM12=EM12+(TLFAIL(ITL,ITLR,ITLT-1)-RSTL)*FACTOR
      EM21=EM21+(TLFAIL(ITL,ITLR-1,ITLT)-RSTL)*FACTOR
      EM22=EM22+(TLFAIL(ITL,ITLR,ITLT)-RSTL)*FACTOR
   90 CONTINUE
C
C  Interpolate first on RVP, then temperature
C
      EMT1=RVPINT*(EM12-EM11)/(RVP2-RVP1)
     *    +(RVP2*EM11-RVP1*EM12)/(RVP2-RVP1)
      EMT2=RVPINT*(EM22-EM21)/(RVP2-RVP1)
     *    +(RVP2*EM21-RVP1*EM22)/(RVP2-RVP1)
C
      RUFAIL=TEMRL*(EMT2-EMT1)/(TMP2-TMP1)
     *      +(TMP2*EMT1-TMP1*EMT2)/(TMP2-TMP1)
      IF(IG.EQ.1) RUPASS=RUFAIL
C
C  Reset fail if less than pass.
C
      IF(RUFAIL.LT.RUPASS) RUFAIL=RUPASS
C
C  Evap. Test Procedure adjustment.
C
      IF(IVTL.LE.4) THEN
      RUPASS=RUPASS*(1.0-ETPCAL(6,MY,1,RVPRL,TEMRL))
      RUFAIL=RUFAIL*(1.0-(ETPCAL(6,MY,2,RVPRL,TEMRL)
     *                   +ETPCAL(6,MY,3,RVPRL,TEMRL))/2.)
      ENDIF
C
C  Multiply pass rate by percentage of passed IV
C  Multiply fail rate by percentage of pressure and purge failed IV
C
      RULOSS=RUPASS*PFRATE(1,IPFCU)
     *      +RUFAIL*(PFRATE(2,IPFCU)+PFRATE(3,IPFCU))
C
C  If speed is below 7 or above 48, extrapolate the running loss rate
C
      IF(SPDRL.LT.RULSPD(1)) RULOSS=RULOSS*RULSPD(1)/SPDRL
      IF(SPDRL.GT.RULSPD(3)) RULOSS=RULOSS*RULSPD(3)/SPDRL
C
      RETURN
      END
