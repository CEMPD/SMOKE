      FUNCTION EV85(IH,KEYEQ,IVGAS,RIVAL,MY,IFDS) 
C 
C  EV85 calculates a hot soak or diurnal rate for the 1985+ HDGV emission 
C  standard, given the equation group and the fuel RVP or Wade Index, 
C  respectively. 
C 
C  Called by EVMAIN. 
C 
C  Calls EV81 and EVUNC. 
C 
C  Changes: (Last change first) 
C 
C  24 Aug 00 @EPA-djb Bug 233, Passed IFDS from EVMAIN. Use IFDS for 
C            EV81 call, instead of 1 (carbureted only).  Eliminate HDSAL. 
C  20 Jun 00 @EPA-ddj Removed IVGAS from parameter list of EVUNC. 
C  28 Sep  1999 @EPA-EG          Updated EVUNC function to reflect Arcadis Hot Soak Analysis. 
C   8 Mar  1999 @Dyntel-BG 2-664 added hourly index to parameter list 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  23 Sept 1998 @DynTel-ZK Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IH,KEYEQ,KEYEQ,IVGAS,RIVAL 
C 
C    common blocks: 
C    /EVAPHD/ HDWGT 
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
C  Sep-20-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar 
C  year removed IY and ICY from the argument list and thus from the calls to 
C  EV81 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAPHD.I' 
      REAL EV85, EV81, EVUNC 
      REAL, INTENT(IN)    :: RIVAL 
      INTEGER, INTENT(IN) :: KEYEQ 
      INTEGER, INTENT(IN) :: IVGAS 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: MY 
C 
C  Heavy duty gasoline vehicles with weight less than 14,000 lbs GVWR (2b & 3) 
C  use the hot soak of LDT times a 1.5 factor.  Heavier HDGV use a 2.0 factor. 
C 
C  Passing vehicles. (KEYEQ = 1) 
C 
      IF(KEYEQ.EQ.1) THEN 
        IF(IVGAS.EQ.6 .OR. IVGAS.EQ.7) THEN 
          EV85=EV81(IFDS,IVGAS,RIVAL,MY)*HDWGT(1) !1.5 factor 
        ELSE  !IV = 8 to 13 and 25 
          EV85=EV81(IFDS,IVGAS,RIVAL,MY)*HDWGT(2) !2.0 factor 
        ENDIF 
C 
C  Failed Purge or Pressure. (KEYEQ = 2 or 3) 
C 
      ELSE 
        IF(IVGAS.EQ.6 .OR. IVGAS.EQ.7) THEN 
          EV85=EVUNC(IH,KEYEQ,RIVAL)*HDWGT(1) 
        ELSE 
          EV85=EVUNC(IH,KEYEQ,RIVAL)*HDWGT(2) 
        ENDIF 
      ENDIF 
      RETURN 
      END 
