      REAL FUNCTION EVUNC(IH,KEYEQ,RIVAL) 
C 
C  EVUNC calculates a hot soak or diurnal for an 1981+ failed 
C  (KEYEQ = 2 or 3) emission equation, given the emission standard 
C  context (veh type, MY, etc) and the fuel RVP. 
C 
C  Called by EVMAIN, EVPRE, EV85, EVSTD. 
C 
C  Changes : (Last change first) 
C 
C  20 Jun 00  @EPA-bag  Added INTENT qualifier to parameter declaration; removed 
C                       IVGAS from parameter list 
C  08 Jun 00: Removed type declarations for intrinsic functions. 
C  28 Sep  1999 @EPA-EG          Updated EVUNC function to reflect Arcadis Hot Soak Analysis. 
C   8 Mar  1999 @Dyntel-BG 2-664 Added hourly index to parameter list 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  25 Sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IH,KEYEQ,RIVAL 
C    common blocks: 
C    /CITRV2/ RVP090 
C    /TEMPS/ TEMHRLY 
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
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'TEMPS.I' 
      INTEGER, INTENT(IN) :: KEYEQ 
      INTEGER, INTENT(IN) :: IH 
      REAL,    INTENT(IN) :: RIVAL 
C 
C 
      IF(KEYEQ.EQ.2) THEN 
C 
C  Hot Soak, Failed Purge 
C 
        IF(RIVAL.GT.RVP090) THEN 
          EVUNC=EXP(1.76223+0.444300*(RIVAL-RVP090) 
     *          +0.05114*(TEMHRLY(IH)-82.0)) 
        ELSE 
          EVUNC=EXP(1.76223+0.552175*(RIVAL-RVP090) 
     *        +0.05114*(TEMHRLY(IH)-82.0)) 
        ENDIF 
C 
      ELSEIF(KEYEQ.EQ.3) THEN 
C 
C  Hot Soak, Failed Pressure 
C 
        IF(RIVAL.GT.RVP090) THEN 
          EVUNC=EXP(1.774+0.444300*(RIVAL-RVP090) 
     *          +0.05114*(TEMHRLY(IH)-82.0)) 
C 
        ELSE 
          EVUNC=EXP(1.774+0.413356*(RIVAL-RVP090) 
     *        +0.05114*(TEMHRLY(IH)-82.0)) 
C 
        ENDIF 
C 
      ENDIF 
C 
C 
      RETURN 
      END 
