C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION EVSTD(KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL,IY,ICY)
C
C  EVSTD calculates the hot soak or diurnal rate for a standard in the
C  controlled to pre-1981 LDGV and LDGT1/2 classes, given the fuel RVP
C  or Wade Index, respectively.
C
C  Called by EVMAIN.
C
C  Calls EV81.
C
C  Input on call:
C
C    parameter list: KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL
C    common blocks:
C    /CITRV2/ RVP090
C    /EVADU1/ EFDU
C    /EVADU2/ R2DU
C    /EVAHS1/ EFHS
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: EVSTD
C
C  Local array subscripts:
C
C  R2HS(2)    -  R2HS ( KSTD )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  R1       R    low RVP or diurnal Index based on low RVP for EV81 call
C  R2       R    high RVP or diurnal Index based on high RVP for EV81 call
C  R2HS     R    hot soak values for R2 based on standard
C
C  Notes:
C
C  The emission rates in this case are to be fitted through the 1981+ LDGV
C  or LDGT1/2 curves.
C  EVSTD was modified for MOBILE4.1 to include pass/fail.
C
C
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /EVADU1/ DUEQ(3,2,3),EFDU(3,2,5,5)
      COMMON /EVADU2/ HIDU(3,2),R2DU(2,2)
      COMMON /EVAHS1/ HSEQ(3,2,9),EFHS(3,2,5,5)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      DIMENSION R2HS(2)
      DATA R2HS/11.5,11.7/
C
      IFDS=1
C
      KSTD=1
      IF(JSTD.EQ.4) KSTD=2
C
      XRIVAL=RIVAL
C
      IF(KEYEVP.EQ.2) GOTO 10
C
C  Set up for Hot Soak calls.
C
      R1=RVP090
      R2=R2HS(KSTD)
C
      EF090=EFHS(KEYEQ,1,JSTD,IVGAS)
      EF115=EFHS(KEYEQ,2,JSTD,IVGAS)
C
      GOTO 20
C
C  Set up for Diurnal calls.
C  (Reset XRIVAL to lower bound of 0.75 for fails)
C
   10 R1=1.0
      R2=R2DU(KSTD,IREJN)
C
      EF090=EFDU(KEYEQ,1,JSTD,IVGAS)
      EF115=EFDU(KEYEQ,2,JSTD,IVGAS)
C
      IF(KEYEQ.GT.1.AND.XRIVAL.LT.0.75) XRIVAL=0.75
C
C  Set up linear interpolation or ratio extrapolation.
C  Handle pass and failed purge/pressure accordingly.
C
   20 IF(KEYEQ.GT.1) GOTO 30
C
C  Pass
C
      Y2=EV81(KEYEVP,IFDS,IVGAS,XRIVAL,IY,ICY)
C
      Y1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2)
     *Y1=EV81(KEYEVP,IFDS,IVGAS,R1,IY,ICY)
C
      IF(XRIVAL.LT.R1) X2=EV81(KEYEVP,IFDS,IVGAS,R1,IY,ICY)
      IF(XRIVAL.GE.R1) X2=EV81(KEYEVP,IFDS,IVGAS,R2,IY,ICY)
C
      GOTO 40
C
C  Failed Purge and Pressure
C
   30 Y2=EVUNC(KEYEVP,KEYEQ,IVGAS,XRIVAL)
C
      Y1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Y1=EVUNC(KEYEVP,KEYEQ,IVGAS,R1)
C
      IF(XRIVAL.LT.R1) X2=EVUNC(KEYEVP,KEYEQ,IVGAS,R1)
      IF(XRIVAL.GE.R1) X2=EVUNC(KEYEVP,KEYEQ,IVGAS,R2)
C
C  Complete remainder of set up.
C
   40 X1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) X1=Y1
C
      IF(XRIVAL.LT.R1) Z2=EF090
      IF(XRIVAL.GE.R1) Z2=EF115
C
      Z1=0.0
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Z1=EF090
C
C  Determine EVSTD
C
      EVSTD=Z1+(Z2-Z1)*(Y2-Y1)/(X2-X1)
C
      RETURN
C
      END
