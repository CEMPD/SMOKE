C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION EVPRE(KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL)
C
C  EVPRE calculates a hot soak or diurnal rate for a pre-control emission
C  standard, given the vehicle type and fuel RVP or Wade Index, respectively.
C
C  Called by EVMAIN.
C
C  Input on call:
C
C    parameter list: KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL
C    common blocks:
C    /CITRV2/ RVP090,RVP115
C    /EVADU1/ EFDU
C    /EVADU2/ R2DU
C    /EVAHS1/ EFHS
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: EVPRE
C
C  Notes:
C
C  JSTD should be 1 for all vehicle types except MC.  MC can be 1 or 2 or 3.
C  EVPRE was modified for MOBILE4.1 to include pass/fail.
C
C
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /EVADU1/ DUEQ(3,2,3),EFDU(3,2,5,5)
      COMMON /EVADU2/ HIDU(3,2),R2DU(2,2)
      COMMON /EVAHS1/ HSEQ(3,2,9),EFHS(3,2,5,5)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      XRIVAL=RIVAL
C
      IF(KEYEVP.EQ.2) GOTO 10
C
C  Set up for Hot Soak call for Pass, Failed Purge, Failed Pressure
C  Reset RIVAL to lower bound of 7.2 RVP)
C
      R1=RVP090
      R2=RVP115
C
      EF090=EFHS(KEYEQ,1,JSTD,IVGAS)
      EF115=EFHS(KEYEQ,2,JSTD,IVGAS)
C
      IF(XRIVAL.LT.7.2) XRIVAL=7.2
C
      GOTO 20
C
C  Set up for Diurnal rate call for Pass, Failed Purge, Failed Pressure
C  (Reset RIVAL to lower bound of 0.75)
C
   10 KSTD=1
C
      R1=1.0
      R2=R2DU(KSTD,IREJN)
C
      EF090=EFDU(KEYEQ,1,JSTD,IVGAS)
      EF115=EFDU(KEYEQ,2,JSTD,IVGAS)
C
      IF(XRIVAL.LT.0.75) XRIVAL=0.75
C
C  Set up linear interpolation or ratio extrapolation.
C  Since pre-control, all calls to EVUNC use KEYEQ=3.
C
   20 Y2=EVUNC(KEYEVP,3,IVGAS,XRIVAL)
C
      Y1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Y1=EVUNC(KEYEVP,3,IVGAS,R1)
C
      IF(XRIVAL.LT.R1) X2=EVUNC(KEYEVP,3,IVGAS,R1)
      IF(XRIVAL.GE.R1) X2=EVUNC(KEYEVP,3,IVGAS,R2)
C
      X1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) X1=Y1
C
      IF(XRIVAL.LT.R1) Z2=EF090
      IF(XRIVAL.GE.R1) Z2=EF115
C
      Z1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Z1=EF090
C
C  Determine EVPRE
C
      EVPRE=Z1+(Z2-Z1)*(Y2-Y1)/(X2-X1)
C
      RETURN
C
      END
