C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION EV85(KEYEVP,KEYEQ,IVGAS,RIVAL,IY,ICY)
C
C  EV85 calculates a hot soak or diurnal rate for the 1985+ HDGV emission
C  standard, given the equation group and the fuel RVP or Wade Index,
C  respectively.
C
C  Called by EVMAIN.
C
C  Calls EV81 and EVUNC.
C
C  Input on call:
C
C    parameter list: KEYEVP,KEYEQ,KEYEQ,IVGAS,RIVAL
C
C    common blocks:
C    /EVAPHD/ HDSAL,HDWGT
C
C  Output on return:
C
C    function: EV85
C
C
C  Notes:
C
C  Certification fuels support (KEYEQ=3) from MOBILE4 was removed
C  in MOBILE4.1.  EV85 was modified for MOBILE4.1 to include
C  pass/fail.  EV85 was modified for MOBILE5v13 to reflect the
C  changes in diurnal, the changes account for the difference
C  between a 1 hour measurement for diurnal and an 8 hour
C  measurement of diurnal emissions.
C
C
      COMMON /EVAPHD/ HDSAL,HDWGT(2),IVVMAP(5)
C
      IF(KEYEQ.GE.2)GOTO 10
C
C  Passed.
C
      EV85=EV81(KEYEVP,1,IVGAS,RIVAL,IY,ICY)*HDWGT(1)*HDSAL
     *    +EV81(KEYEVP,1,IVGAS,RIVAL,IY,ICY)*HDWGT(2)*(1.0-HDSAL)
      IF(KEYEVP.EQ.2) EV85=EV85*2.3297
      RETURN
C
C  Failed Purge or Pressure.
C
   10 EV85=EVUNC(KEYEVP,KEYEQ,IVGAS,RIVAL)*HDWGT(1)*HDSAL
     *    +EVUNC(KEYEVP,KEYEQ,IVGAS,RIVAL)*HDWGT(2)*(1.0-HDSAL)
      RETURN
      END
