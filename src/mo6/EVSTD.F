      REAL FUNCTION EVSTD(IH,KEYEQ,JSTD,IVGAS,RIVAL,MY) 
C 
C  EVSTD calculates the hot soak for a standard in the controlled 
C  to pre-1981 LDGV and LDGT classes, given the fuel RVP or Wade 
C  Index, respectively. 
C 
C  Called by EVMAIN. 
C 
C  Calls EV81. 
C 
C  Changes : (Last change first) 
C 
C  20 Jun 00 @EPA-ddj Removed IVGAS from parameter list of EVUNC. 
C   2 Jun 00 @EPA-ddj Bug 86, Added SAVE to R2HS type statement. 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  28 Sep  1999 @EPA-EG          Updated EVUNC function to reflect Arcadis Hot Soak Analysis. 
C   8 Mar  1999 @Dyntel-BG 2-664 added hourly index to parameter list 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  25 Sept 1998 @DynTel-Zk 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IH,KEYEQ,JSTD,IVGAS,RIVAL 
C    common blocks: 
C    /CITRV2/ RVP090 
C    /EVAHS1/ EFHS 
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
C  or LDGT curves. 
C  EVSTD was modified for MOBILE4.1 to include pass/fail. 
C  Sep-20-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar year 
C  removed IY and ICY from the argument list and thus from the calls 
C  to EV81 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE  
      INCLUDE 'CITRV2.I' 
      INCLUDE 'EVAHS1.I' 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: KEYEQ 
      INTEGER, INTENT(IN) :: JSTD 
      INTEGER, INTENT(IN) :: IVGAS 
      INTEGER, INTENT(IN) :: MY 
      REAL,    INTENT(IN) :: RIVAL 
      INTEGER IFDS, KSTD 
      REAL XRIVAL,R1,R2,EF090,EF115,EV81,EVUNC,X1,Y1,Y2,Z1,X2,Z2 
      REAL, DIMENSION(2), SAVE :: R2HS=(/11.5,11.7/) 
C 
      IFDS=1 
C 
      KSTD=1 
      IF(JSTD.EQ.4) KSTD=2 
C 
      XRIVAL=RIVAL 
C 
C  Set up for Hot Soak calls. 
C 
      R1=RVP090 
      R2=R2HS(KSTD) 
C 
      EF090=EFHS(KEYEQ,1,JSTD,IVGAS) 
      EF115=EFHS(KEYEQ,2,JSTD,IVGAS) 
C 
C  Set up linear interpolation or ratio extrapolation. 
C  Handle pass and failed purge/pressure accordingly. 
C 
      IF(KEYEQ.GT.1) GOTO 30 
C 
C  Pass 
C 
      Y2=EV81(IFDS,IVGAS,XRIVAL,MY) 
C 
      Y1=0.0 
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Y1=EV81(IFDS,IVGAS,R1,MY) 
C 
      IF(XRIVAL.LT.R1) X2=EV81(IFDS,IVGAS,R1,MY) 
      IF(XRIVAL.GE.R1) X2=EV81(IFDS,IVGAS,R2,MY) 
C 
      GOTO 40 
C 
C  Failed Purge and Pressure 
C 
   30 Y2=EVUNC(IH,KEYEQ,XRIVAL) 
C 
      Y1=0.0 
      IF(XRIVAL.GE.R1.AND.XRIVAL.LE.R2) Y1=EVUNC(IH,KEYEQ,R1) 
C 
      IF(XRIVAL.LT.R1) X2=EVUNC(IH,KEYEQ,R1) 
      IF(XRIVAL.GE.R1) X2=EVUNC(IH,KEYEQ,R2) 
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
