C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION EVHI(KEYEVP,KEYEQ,RIVAL)
C
C  EVHI calculates the high altitude, 1972-76 model year, hot soak or
C  diurnal for LDGVs and LDGT1s given the fuel RVP or Wade Index,
C  respectvely.
C
C  Called by EVMAIN.
C
C  Input on call:
C
C    parameter list: KEYEVP,RIVAL
C    common blocks:
C    /CITRV2/ RVP090,RVP115
C    /EVADU2/ HIDU,R2DU
C    /EVAHIA/ HIHS
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: EVHI
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  R1       R    low RVP or diurnal Index
C  R2       R    high RVP or diurnal Index
C
C  Notes:
C
C  Modified for MOBILE4.1 to add support for pressure/purge fail.
C
C
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /EVADU2/ HIDU(3,2),R2DU(2,2)
      COMMON /EVAHIA/ HIHS(3,2),HIADJ
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      IF(KEYEVP.EQ.2) GOTO 10
C
C  Hot soak rate call
C
      R1=RVP090
      R2=RVP115
      EVHI=(HIHS(KEYEQ,2)*(RIVAL-R1)-HIHS(KEYEQ,1)*(RIVAL-R2))/(R2-R1)
      RETURN
C
C  Diurnal rate call (Note: EVHI currently only called when IREJN = 2.)
C
   10 R1=1.0
      KSTD=1
      R2=R2DU(KSTD,IREJN)
      EVHI=(HIDU(KEYEQ,2)*(RIVAL-R1)-HIDU(KEYEQ,1)*(RIVAL-R2))/(R2-R1)
      RETURN
C
      END
