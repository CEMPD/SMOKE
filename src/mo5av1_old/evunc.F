C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION EVUNC(KEYEVP,KEYEQ,IVGAS,RIVAL)
C
C  EVUNC calculates a hot soak or diurnal for an 1981+ failed
C  (KEYEQ = 2 or 3) emission equation, given the emission standard
C  context (veh type, MY, etc) and the fuel RVP.
C
C  Called by EVMAIN, EVPRE, EV85, EVSTD.
C
C  Input on call:
C
C    parameter list: KEYEVP,KEYEQ,IVGAS,RIVAL
C    common blocks:
C    /CITRV2/ RVP090
C    /CITCIN/ UDI,IUDI
C    /REGION/ IREJN
C    /TECHFL/ CONSTF
C    /TEMPS/ TEMEVP
C
C  Output on return:
C
C    function: EVUNC
C
C
C  Notes.
C
C  Modified for MOBILE4.1 to add support for pressure/purge fail.
C  Modified for MOBILE5v13 to reflect the differences between
C  measurement for 1 hour of diurnal as opposed to 8 hours of
C  diurnal emission.  Added common CONSTF for failed constant.
C
C
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /EVADU2/ HIDU(3,2),R2DU(2,2)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /TECHFL/ CONSTF
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      IF(KEYEVP.EQ.2) GOTO 30
C
      GOTO(10,20),KEYEQ-1
C
C  Hot Soak, Failed Purge
C
   10 EVUNC=EXP(1.76223+0.4443*(RIVAL-RVP090)
     1      +0.05114*(TEMEVP(1)-82.0))
      RETURN
C
C  Hot Soak, Failed Pressure
C
   20 EVUNC=EXP(1.774+0.4443*(RIVAL-RVP090)
     1      +0.05114*(TEMEVP(1)-82.0))
      RETURN
C
C  Diurnal
C
   30 GOTO(40,50),KEYEQ-1
C
C  Diurnal, Failed Purge
C
   40 EVUNC=-6.1152+11.6952*RIVAL+CONSTF
      IF(IVGAS.GE.2) EVUNC=EVUNC*1.17
      EVUNC=AMAX1(EVUNC,0.0)
      RETURN
C
C  Diurnal, Failed Pressure
C
   50 EVUNC=2.8997+9.0903*RIVAL+CONSTF
      IF(IVGAS.GE.2) EVUNC=EVUNC*1.17
      EVUNC=AMAX1(EVUNC,0.0)
      RETURN
C
      END
