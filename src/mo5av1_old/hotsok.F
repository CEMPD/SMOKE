C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION HOTSOK(MY,IDX,IVGAS,IY)
C
C  HOTSOK selects the yearly value for hot soak HC.
C
C  Called by CCEVRT.
C
C  Calls EVMAIN, ETPCAL.
C
C  Input on call:
C
C    parameter list: MY,IDX,IVGAS
C    common blocks:
C    /CITCIN/ PFRATE,IPFCU
C    /CITRV2/ RVPHS
C    /EVAHS1/ HSEQ
C    /EVAPAR/ ISTDC,ISTDT,ISTDP,FINJ
C    /EVAPHD/ IVVMAP
C    /FLAGS3/ TEMFLG
C    /FLAGS4/ HCFLAG
C    /MAXIMA/ MAXYRS
C    /TAMOUT/ THS
C    /TEMPS/  AMBT,TEMMIN,TEMMAX,TEMEVP
C
C  Output on return:
C
C    function: HOTSOK
C    /EVAPGR/ EVP
C
C  Local array subscripts
C
C  HSTCF  (3) - HSTCF  (IFDS)
C
C  Local variable / array dictionary
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  HSPASS   R    Pass emission rate to insure FPressure and FPurge
C                never fall below the pass rate
C  HSPF     R    the EVMAIN hot soak loss for Pass, Failed Purge, Failed
C                Pressure and controlled and uncontrolled
C  HSTCF    R    hot soak temperature correction factor for CARB, TBI,
C                and PBI
C  HSTEMP   R    temperature used in the hot soak temperature correction
C                factor equation
C  IVV      I    vehicle (IVV=1) or truck (IVV=2) classification switch
C  RIVAL    R    hot soak RVP adjusted for weathering
C  FULADJ   R    FTP vs In-Use fuel level adjustment factor for fuel
C                injected vehicles
C  TRIPFR   R    fraction of vehicle fleet with daily trips
C
C  Notes:
C
C  Certification fuels support (KEYEQ=3) from MOBILE4 was removed in MOBILE4.1.
C  Modified for MOBILE4.1 to add support for pressure/purge fail.
C
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /EVAHS1/ HSEQ(3,2,9),EFHS(3,2,5,5)
      COMMON /EVAPGR/ EVP(4),GREVP(4,9),VGREVP(4),PARD(3),GRPD(3,9)
      COMMON /EVAPAR/ IFDS,ISTDC,ISTDT,ISTDP,FINJ(2)
      COMMON /EVAPHD/ HDSAL,HDWGT(2),IVVMAP(5)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION HSTCF(3)
      DATA FULADJ/0.88/
      ICY=MY-IDX+MAXYRS
C
C  Hot soak is zero if any evaporative temperature is 40F or less,
C  or if minimum temperature is 25F or less.
C
      HOTSOK=0.0
      EVP(1)=0.0
      IF(TEMMAX   .LE.40.0 .OR.
     *   AMBT     .LE.40.0.AND.TEMFLG.EQ.2 .OR.
     *   TEMEVP(1).LE.40.0 .OR.
     *   TEMEVP(2).LE.40.0 .OR.
     *   TEMMIN   .LE.25.0) RETURN
C
      KEYEVP=1
      RIVAL=RVPHS
C
C  Lower bound on Hot Soak RVP is 5.0, reset if necessary.
C
      IF(RIVAL.LT.5.0) RIVAL=5.0
C
C  Select the base temperature to be used for the hot soak correction
C  factor. Use the calculated hot soak temperature as default. However,
C  if that temperature is below 70F, reset to 70F, then interpolate later.
C
      HSTEMP=TEMEVP(1)
      IF(HSTEMP.LT.70.0) HSTEMP=70.0
C
C  Compute the Hot Soak Temperature Correction Factor for each fuel
C  delivery system.
C
      IVV=IVVMAP(IVGAS)
C
      DO 20 IFDS=1,3
C
      HSTCF(IFDS)=(HSEQ(IFDS,IVV,6)
     *            +HSEQ(IFDS,IVV,7)*HSTEMP
     *            +HSEQ(IFDS,IVV,8)*HSTEMP**2)
     *            /HSEQ(IFDS,IVV,9)
C
   20 CONTINUE
C
C  All HDGV AND pre-1981 LDGV/T are 100% carbureted, but are affected by
C  pass/fail.  For these vehicle classes, figure a carbureted number
C  accordingly and then skip the fuel injected rate case.
C  Motorcycles are 100% carbureted and not affected by pass/fail.
C
C  Loop and combine for Passed, Failed Purge and Failed Pressure
C  and for carbureted and fuel injected vehicles.
C  (MC do not have purge or pressure checks.)
C
      HSPASS=0.0
C
      MAXKEQ=3
      IF(IVGAS.EQ.5) MAXKEQ=1
C
      DO 50 KEYEQ=1,MAXKEQ
C
      MAXFDS=3
      IF(IVGAS.GE.4.OR.MY.LT.1981.OR.KEYEQ.GT.1) MAXFDS=1
C
C  Loop through fuel delivery systems
C
C     IFDS = 1 : Carbureted
C            2 : Throttle-body fuel injection
C            3 : Ported or Multi-point fuel injection
C
      DO 40 IFDS=1,MAXFDS
C
      IF(IFDS.EQ.1) ISTUSE=ISTDC
      IF(IFDS.EQ.2) ISTUSE=ISTDT
      IF(IFDS.EQ.3) ISTUSE=ISTDP
C
      CALL EVMAIN(MY,KEYEVP,KEYEQ,IVGAS,ISTUSE,RIVAL,HSPF,IY,ICY)
      JFDS=IFDS
      IF(IVGAS.GE.4.OR.MY.LT.1981) JFDS=1
C
C  Skip the TCF for Failed Purge and Failed Pressure and only if
C  the last standard group
C
      IF(KEYEQ.GT.1) THEN
        IF(IVGAS.EQ.1.AND.ISTUSE.GE.5) GOTO 60
        IF(IVGAS.EQ.2.AND.ISTUSE.GE.5) GOTO 60
        IF(IVGAS.EQ.3.AND.ISTUSE.GE.3) GOTO 60
        IF(IVGAS.EQ.4.AND.ISTUSE.GE.2) GOTO 60
      ENDIF
      HSPF=HSPF*HSTCF(JFDS)
C
C  Subtract off resting losses
C
   60 IVTEMP=IVGAS
      IF(IVGAS.EQ.5) IVTEMP=8
      HSPF=AMAX1(HSPF-RSTLOS(MY,IVTEMP,HSTEMP),0.0)
C
C  Don't allow purge or pressure rates to fall below pass rates.
C
      IF(KEYEQ.GT.1.AND.HSPF.LT.HSPASS) HSPF=HSPASS
C
C  All HDGV and pre-1981 LDGV/T and MC are 100% carbureted
C  so skip the fuel injection weighting.
C
      IF(MAXFDS.EQ.1) GO TO 30
C
C  Weight carbureted and fuel injected rates by corresponding sales fractions
C  For fuel injected vehices, adjust for the difference between in-use
C  and FTP fuel tank levels
C
      IF(IFDS.EQ.1) HSPF=(1.0-FINJ(1)-FINJ(2))*HSPF
      IF(IFDS.GT.1) HSPF=HSPF*FINJ(IFDS-1)*FULADJ
C
C  Save the composite pass rate for future purge/pressure checks.
C
   30 IF(KEYEQ.EQ.1) HSPASS=HSPASS+HSPF
C
C  Apply Pass, Failed Purge, and Failed Pressure weighting.
C
      HSPF=HSPF*PFRATE(KEYEQ,IPFCU)
C
C  If temperature is below 70F, interpolate with assumption
C  that HS=0 AT 40F.
C
      IF(TEMEVP(1).LT.70.0) HSPF=HSPF*(TEMEVP(1)-40.0)/(70.0-40.0)
C
C  Evap. Test Procedure adjustment calculation.
C
      IF(IVGAS.LE.4) HSPF=HSPF*(1.0-ETPCAL(1,MY,KEYEQ,RVPHS,HSTEMP))
C
C  Combine the hot soak components together to form a fleet value.
C
      EVP(1)=EVP(1)+HSPF
C
   40 CONTINUE
C
   50 CONTINUE
C
C  Multiply fleet value by the fraction of vehicles with daily trips.
C  HDGV and MC daily trip fraction is 1.0, all other vehicles use age
C  (age=MAXYRS+1-IDX) equation (or average value of 0.7615, if HCFLAG=6).
C
      IF(HCFLAG.NE.6) TRIPFR=(87.07-1.72*FLOAT(MAXYRS+1-IDX))/100.
      IF(HCFLAG.EQ.6) TRIPFR=0.7615
      IF(IVGAS.GE.4) TRIPFR=1.0
C
      HOTSOK=EVP(1)*TRIPFR
C
      RETURN
      END
