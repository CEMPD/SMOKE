      REAL FUNCTION EVPRE(IH,KEYEQ,JSTD,IVGAS,RIVAL) 
C 
C  EVPRE calculates a hot soak rate for a pre-control emission standard, 
C  given the vehicle type and fuel RVP or Wade Index, respectively. 
C 
C  Called by EVMAIN. 
C 
C  Changes: (Last change first) 
C 
C  20 Jun 00 @EPA-ddj Removed IVGAS from parameter list of EVUNC. 
C   8 Mar  1999 @Dyntel-BG 2-664 added hourly index to parameter list 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  22 January 99 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IH,KEYEQ,JSTD,IVGAS,RIVAL 
C    common blocks: 
C    /CITRV2/ RVP090,RVP115 
C    /EVAHS1/ EFHS 
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
      IMPLICIT NONE 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'EVAHS1.I' 
C 
C     Declare parameter list 
C 
      INTEGER IH 
      INTEGER KEYEQ 
      INTEGER JSTD 
      INTEGER IVGAS 
      REAL    RIVAL 
C 
C     Declare called functions 
C 
      REAL EVUNC 
C 
C     Declare local variables 
C 
      REAL EF090 
      REAL EF115 
      REAL R1 
      REAL R2 
      REAL X1 
      REAL X2 
      REAL XRIVAL 
      REAL Y1 
      REAL Y2 
      REAL Z1 
      REAL Z2 
C 
      XRIVAL=RIVAL 
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
C  Set up linear interpolation or ratio extrapolation. 
C  Since pre-control, all calls to EVUNC use KEYEQ=3. 
C 
      Y2=EVUNC(IH,3,XRIVAL) 
C 
      Y1=0.0 
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Y1=EVUNC(IH,3,R1) 
C 
      IF(XRIVAL.LT.R1) X2=EVUNC(IH,3,R1) 
      IF(XRIVAL.GE.R1) X2=EVUNC(IH,3,R2) 
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
