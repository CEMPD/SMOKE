      REAL FUNCTION 
     *     IMTEST(MY,ACTUAL_MY,JDX,VAGE_LSTIM,IP,IV,ISR,IFDS,IMPGM,ISTD, 
     *     ICY,IDX,IH,IROAD,CAT,INERR) 
C 
C  IMTEST returns the correction factor for an inspection and 
C  maintenance program test group. This include OBD I/M. 
C  The correction is applicable to the 
C  1981+ model years for light duty gasoline vehicles and trucks. In 
C  periodic IM programs there may be two integer values of VAGE_LSTIM for each 
C  model year. The fraction of vehicles corresponding to each value of 
C  VAGE_LSTIM is calculated in IMPROG where the total additive IM correction 
C  factor is also computed. 
C 
C  Called by IMPROG 
C 
C  Calls ARBEF, EF_HI_FRCN, EF81_93, EF_POST93, IDRATE, METHANE 
C 
C  Changes: 
C 
C  09 May 01 @EPA-elg Added Exhaust I/M Factor for 1996+ Model Years. 
C  04 Dec 00 @EPA-bag Task X6, added ACTUAL_MY to parameter list 
C  06 Sep 00 @EPA-elg BUG247 Added code to properly calculate after I/M 
C                            repair levels for 1994-95 MYs. 
C  29 Aug 00 @EPA-elg BUG241 Replaced the calls to METHANE with calls to IMCF. 
C  25 Jul 00 @EPA-djb Added methane speciation and added a special 
C            1994-1995 model year case to the IMTEST calculation. 
C            Changed IMTT from REAL to INTEGER value. 
C  26 Jun 00 @EPA-elg   Added code to apply I/M Discount. 
C  14 Jun 00 @EPA-elg   Replaced EF_POST with EF_CASE. 
C  14 Jun 00 @EPA-elg   Added T2CRT50K to repair effects calculation. 
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list 
C  08 Jun 00 @EPA-djb Bug113 Added ISTD to IDRATE call for 
C            the EF_HI_FRCN call. 
C  15 May 00 @EPA-elg   Added condition to prevent an OBD I/M program on pre-1996 MYs 
C                       or an Exhaust I/M test on post-1996 MYs. 
C  24 Apr 00 AIR Task 02: Included INERR in IMTEST, EF_HI_FRCN and IDRATE 
C                         parameter lists. 
C  09 SEP 99 @DynTel-BG Renamed functions HIEMIT.FOR, HNBEF.FOR, and HN94B.FOR 
C            to EF_HI_FRCN, EF81_93, and EF_POST93, respectively 
C  26 AUG 99 @DynTel-HXQ Deleted the non-used flag IGM 
C  16 Aug 99 @DynTel-bg  2-663  All logical and arithmetic expressions 
C            evaluated with the same variable types, simplified logic 
C            for determining whether a correction need to be calculated, 
C            CASE statement for model year ranges. 
C  29 July 99 @Dyntel-JWR 2-663 Redefined the variable IMPFIX to  
C            reflect a correction in the waived vehicle term in the 
C            I/M credit model (3rd term in eqn 12a of M6IM001). As 
C            a result IMPFIX no longer represents the number of  
C            repaired vehicles  
C  28 May 99 @DynTel-HXQ 2-650 Removed units conversions for start 
C            emissions 
C  18 May 99 @DynTel-HXQ 2-650 Changing the code of converting 
C                              start emissions from g/start to g/mile. 
C  28 Apr 99 @Dyntel-ddj 2-694, Corrected conditional on IF statement 
C            from IMTT.NE.12 to IMTT.LT.12.  Previous version prevented 
C            execution of IMTT.GE.12 ELSE statement.  Also added 
C            statement function AR94B. 
C  17 Mar 99 @Dyntel-RJD 2-694 Extended the scope of IMTEST so that it 
C            also returns the I/M benefits for 1994+ model year LDGV/Ts. 
C  13 Oct 98 @Dyntel-RJD 2-663 Created IMTEST to facilitate I/M benefits 
C            calculations in IMPROG. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,JDX,VAGE_LSTIM,IP,IV,ISR,IFDS,IMPGM,ISTD,INERR 
C    common blocks: 
C    /IMPAR1/ CRIM,WAV81PLUS 
C    /IMPAR2/ ITEST 
C    /IMPAR6/ IFREQ 
C    /MAXIMA/ MAXYRS 
C    /MYRSAV/ AMAR 
C    /LDGBER/ CRT050K 
C    /LDGIMC/ WAVRDC,NONCOM,ASEBZML,IDR,VMTIM 
C    /LDGOBD/ HIFB3, HIHB, OBDREP, R_S_CONV,HIGH3,HIGH4 
C    /TPDCOM/ TPDCO 
C 
C  Output on return: 
C 
C    function: IMTEST 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  ----------------------------------------------------- 
C   HIGH_FRC    R    Fraction of all vehicles that are high emitters 
C   H3          R    Fraction of high emitters after correction for OBD 
C   H4          R    Fraction of high emitters after correction for OBD and 
C                    OBD I/M. 
C   HIGH        I    Integer whose value is 2, indicating high emitter class. 
C   HN94_TMP    R    Local variable for 1994+ after repair emissions. 
C   IAIM        I    The integer part of the age VAGE_LSTIM. 
C   IMPFIX      R    Fraction of vehicles that are essentially repaired. 
C   IMTT        I    I/M test type. 
C   MEF         R    Local value for basic emission rate. 
C   MEFARBEF    R    Local value for after repair emissions. 
C   MEFHIGH     R    Local value for high emitter emission rate. 
C   NORMAL      I    Integer whose value is 1, indicating normal emitter class. 
C   NONCOM      R    None compliance rate. 
C   RAIM        R    For modeling purposes, in selective IM programs such as 
C                    RSD and COIM vehicles are assumed to be eligible for an 
C                    IM test every six months following the date of sale of 
C                    the vehicle into the fleet. RAIM is the difference 
C                    between the age of the vehicle AIM and its integer part 
C                    IAIM. RAIM is always 0 or 0.5. 
C   VMTAIM      R    Vehicle mileage in 10000 mile units at ge VAGE_LSTIM. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'LDGBER.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'STDIST.I' 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'TIER2.I' 
      INCLUDE 'WEEKLY.I' 
C 
C  Declare external functions 
C 
      REAL    EF81_93 
      REAL    EF_CASE 
      REAL    EF_HI_FRCN 
      REAL    IDRATE 
      REAL    ARBEF 
C 
C  Declare parameter list 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: IMPGM 
      INTEGER, INTENT(IN) :: ISTD 
      INTEGER, INTENT(INOUT) :: INERR 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IROAD 
C 
      REAL,    INTENT(IN) :: CAT 
      REAL,    INTENT(IN) :: VAGE_LSTIM 
C 
C  Declare local variables 
C 
      REAL    HIGH_FRC 
      REAL    HI3 
      REAL    HI4 
      REAL    HN94_TMP 
      REAL    IMPFIX 
      REAL    MEF 
      REAL    MEFARBEF 
      REAL    MEFHIGH 
      REAL    NONCOM 
      REAL    RAIM 
      REAL    VMTAIM 
      INTEGER IAIM 
      INTEGER IMTT 
      INTEGER, SAVE :: NORMAL=1,HIGH=2 
C 
C 
      IMTEST = 0.0 
C 
C     Determine whether or not to compute a correction factor 
      IF (VAGE_LSTIM.GT.0.49.AND.IMPGM.NE.0) THEN 
C 
        MEF  = 0.0 
        MEFHIGH = 0.0 
        MEFARBEF = 0.0 
C 
        IAIM = INT(VAGE_LSTIM) 
        RAIM = VAGE_LSTIM-REAL(IAIM) 
        IMTT = ITEST(IMPGM) 
C 
        IF(ITEST(IMPGM).LT.12 .AND.MY.GE.1996) IMTT = 12 
C 
        IF (IAIM.EQ.0) THEN 
          VMTAIM = 10.0*RAIM*AMAR(1,IV) 
        ELSE 
          VMTAIM = VMTIM(IAIM,IV) + 10.0*RAIM*AMAR(IAIM+1,IV) 
        END IF 
C 
C       Determine the post '94 after repair rates 
C 
        IF(ISTD.NE.4) THEN 
          AR94B = OBDREP*CRT050K(ISTD,IP,IV)* 
     &                        R_S_CONV(ISTD,IP,ISR) 
        ELSEIF(ISTD.EQ.4) THEN                    !If Tier2 use T2CRT50K 
          AR94B = OBDREP*T2CRT50K(IV,IP,JDX)* 
     &                        R_S_CONV(ISTD,IP,ISR) 
        ENDIF 
C 
C       After repair emissions are either Normals or capped at 1.5 times 
C       Cert Standards. 
C 
        HN94_TMP = EF_CASE(MY,JDX,IP,IV,ISR,VMTAIM,ISTD,NORMAL,IFDS) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *            IH,IROAD,CAT,HN94_TMP) 
        MEF = AR94B 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *            IH,IROAD,CAT,MEF) 
        AR94B = MEF 
        AR94B=AMIN1(AR94B,HN94_TMP) 
C 
C Calculate the correction factor corresponding to the program 
C number, IMPGM. 
C 
C 
C       determine the fraction of vehicles which are I&M non-compliant 
        NONCOM =  1.0 - (CRIM(IMPGM) / 100.0) 
C 
C       determine the fraction of vehicles which have been repaired 
        IMPFIX =  1.0 - WAV81PLUS(IMPGM)*WAVRDC(IMPGM) - NONCOM 
C 
C       determine the fraction of vehicles which are high emitters 
        HIGH_FRC = EF_HI_FRCN(MY,IP,IV,VMTAIM,IFDS,ISTD,INERR) 
C 
C         Test type indices greater than or equal to 12 correspond to 
C         OBD I&M programs (see also GETIPR.FOR) : 
C 
C          1 = idle test 
C          2 = 2500 rpm idle test 
C          3 = loaded idle test 
C          4 = IM240 test 
C          5 = ASM test procedure 
C          6 = ASM 2525 PHASE-IN 
C          7 = ASM 2525 FINAL 
C          8 = ASM 5015 PHASE-IN 
C          9 = ASM 5015 FINAL 
C         10 = ASM 2525/5015 PHASE-IN 
C         11 = ASM 2525/5015 FINAL 
C         12 = OBD I/M 
C 
        IF(MY.GE.1981 .AND. MY.LE.1993) THEN     !Tier 0 Case 
          MEFHIGH  = EF81_93(MY,IP,IV,VMTAIM,IFDS,ISR,HIGH) 
          CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,2, 
     *            IH,IROAD,CAT,MEFHIGH) 
          MEFARBEF = 
     *        ARBEF(MY,JDX,VAGE_LSTIM,IP,IV,VMTAIM,IFDS,ISR,IMPGM,1) 
          CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *            IH,IROAD,CAT,MEFARBEF) 
          IMTEST = -HIGH_FRC* 
     *       IDRATE(MY,JDX,VAGE_LSTIM,IP,IV,ISR,IFDS,ISTD,IMPGM,INERR) 
     *         * ( (1.0 - WAV81PLUS(IMPGM)*WAVRDC(IMPGM) - NONCOM) 
     *         *  MEFHIGH - IMPFIX *  MEFARBEF ) 
C 
C 
        ELSEIF(MY.EQ.1994 .OR. MY.EQ.1995) THEN  !Non-OBD Case 
            IF(ISTD.EQ.1) THEN 
              MEFARBEF = 
     *        ARBEF(MY,JDX,VAGE_LSTIM,IP,IV,VMTAIM,IFDS,ISR,IMPGM,1) 
              CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *                  IH,IROAD,CAT,MEFARBEF) 
            ELSEIF(ISTD.GT.1) THEN 
              MEFARBEF = AR94B 
            ENDIF 
C 
            MEFHIGH = EF_CASE(MY,JDX,IP,IV,ISR,VMTAIM,ISTD,HIGH,IFDS) 
            CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,2, 
     *            IH,IROAD,CAT,MEFHIGH) 
          IMTEST = -HIGH_FRC* 
     *       IDRATE(MY,JDX,VAGE_LSTIM,IP,IV,ISR,IFDS,ISTD,IMPGM,INERR) 
     *         * ( (1.0 - WAV81PLUS(IMPGM)*WAVRDC(IMPGM) - NONCOM) 
     *         *  MEFHIGH - IMPFIX *  MEFARBEF ) 
C 
C 
        ELSEIF(MY.GE.1996) THEN                  !OBD Case 
          IF(IMTT.GE.12) THEN                    !OBD I/M Cases 
            MEFHIGH=EF_CASE(MY,JDX,IP,IV,ISR,VMTAIM,ISTD,HIGH,IFDS) 
            CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,2, 
     *            IH,IROAD,CAT,MEFHIGH) 
            HI3 = HIGH3(JDX,IP,IV) 
            HI4 = HIGH4(JDX,IP,IV) 
            IF (VMTAIM.GT.3.6) THEN 
              IMTEST = (HI4 - HI3) 
     *         * ((1.0 - WAV81PLUS(IMPGM) * WAVRDC(IMPGM) - NONCOM) 
     *         * MEFHIGH - IMPFIX * AR94B) 
            ELSE    !Under warranty: 100% compliance & no waivers. 
              IMTEST = (HI4 - HI3) * (MEFHIGH-AR94B) 
            END IF  !VMTAIM Condition 
          ELSE                                   !Non-OBD I/M should not occur. 
            IMTEST = 0.0 
          ENDIF     !IMTT Condition 
C 
        ELSE                                     !Pre-1981 MY should not occur. 
            IMTEST = 0.0 
        ENDIF       !Model Year Condition 
C 
C  Apply I/M program discount. 
C 
        IMTEST = IMTEST * DISCNT(IP,INTYP(IMPGM)) 
C 
      ENDIF 
C 
      RETURN 
      END 
