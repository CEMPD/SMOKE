      REAL FUNCTION IMPROG(ICY,MY,ACTUAL_MY,JDX,IP,IV,VMTAGE,ISR,IFDS, 
     *                     ISTD,INERR,IDX,IH,IROAD,CAT) 
C 
C  IMPROG returns the additive correction factor for an inspection and 
C  maintenance program. The correction is applicable to 1981+ light 
C  gasoline vehicles and trucks. In general, the vehicle population 
C  corresponding to any given vehicle type and model year can be 
C  partitioned into groups depending on the age VAGE_LSTIM of the vehicle on 
C  the day it is tested. The fraction of vehicles corresponding to each 
C  value of VAGE_LSTIM is calculated. The total additive IM correction factor 
C  is then computed as a weighted average over all the IM test groups. 
C 
C  Called by HRLOOP. 
C 
C  Calls AMIN1, AMAX1, FINDIM, EF_HI_FRCN, EF_POST93 and IMTEST 
C 
C  Changes: 
C 
C  09 May 01 @EPA-elg Added Exhaust I/M Factor for 1996+ Model Years. 
C  03 Jan 01 @EPA-djb Bug307, Removed IMPGM loop from REDPROG calculation. 
C  04 Dec 00 @EPA-bag Task X6, added ACTUAL_MY to parameter list 
C  02 Nov 00 @EPA-djb Bug299, skip calculations for evap I/M cases. 
C  20 Sep 00 @EPA- BG Added constant IEXHAUST=1 for call to FINDIM 
C  15 Sep 00 AIR Task 08: Updated FINDIM calls to select exhaust I/M. 
C  29 Aug 00 @EPA-djb bug241 Removed call to IMPROG_NP 
C  25 Jul 00 @EPA-djb Added methane speciation, removed non-I/M IMREMAIN 
C            calculations and altered/added some comments. 
C  14 Jun 00 @EPA-elg Replaced EF_POST with EF_CASE. 
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list 
C  08 Jun 00 @EPA-djb Bug113 Added ISTD to EF_HI_FRCN call. 
C  17 May 00 @EPA-bag bug#88 replaced HI3 with its intented value, i.e., 
C                         HIGH3(JDX,IP,IV,) from obd_hi_frcn.for 
C  24 Apr 00 AIR Task 02: Included INERR into IMPROG, IMTEST, IMPROG_NP, 
C                         and EF_HI_FRCN parameter lists. 
C  11 Apr 00 @EPA-eg            hrloop revisions 
C  21 Mar 00 @EPA-eg            hrloop revisions 
C  30 Dec 99 @DynTel-ddj 1-021 Added 1.5 x Tier2 Standard cap on 
C            after repair emissions. 
C  21 Oct 99 @Dyntel - JWR- Introduced scaling factor REDPROG to  
C            correct the January and July credits for the case when  
C            the program start year equals the calendar year, by 
C            removing the spurious Oct -Jan contribution. Fixed OBD credits 
C            for 1994 and 1995 model years overwriting all other I&M 
C            credits. 
C  13 SEP 99 @DynTel-BG Added OBD_IM constant parameter; non-I&M OBD credit for 
C            1994 and 1995 needed a sign reversal (OBD for these years increased 
C            emissions rather than decreasing them). 
C  09 SEP 99 @DynTel-BG Renamed functions HIEMIT.FOR and HN94B.FOR 
C            to EF_HI_FRCN and EF_POST93, respectively. 
C  11 Aug 99 @DynTel-bag 2-663 Variable names changed, logical and arithmetic 
C            expressions now evaluated with same variable types, replaced constants 
C            with descriptive parameters, CASE statement for separation of COIM and 
C            COIM and RSD non-periodic I&M programs from the other periodic programs, 
C            re-introduced max. and min. values for post 1993 after repair rates, and 
C            AR94B, determined in IMPROG 
C            from IMPROG_NP 
C   8 Jul 99 @DynTel-ddj 2-705 Removed bypass around I/M calculation 
C            when FTP output option selected. 
C   1 Jun 99 @DynTel-ddj 2-700 Changed implementation logic of FTP_FLAG 
C            to be consistent with other flags. 
C  21 May 99 @DynTel-ddj 2-700, Added code for FTP output option 
C            (FTP_FLAG). 
C  29 Apr 99 @DynTel-ddj 2-694, Added local variable AR94B. 
C  16 Mar 99 @Dyntel-RJD 2-694 Modified IMPROG so that it also returns 
C            the I/M additive correction factor for 1994+ model year 
C            vehicles equipped with OBD. 
C  13 Oct 98 @Dyntel-RJD 2-663 Created IMTEST to correct the 1981-1993 
C            light duty gasoline vehicles and trucks basic emission 
C            factor for the benefits of IM programs. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,JDX,IP,IV,VMTAGE,ISR,IFDS,ISTD,INERR 
C    common blocks: 
C    /EVAL/   FRC, FRN 
C    /IMPAR1/ ICYIM,MODYR,EXEMPTAGE 
C    /IMPAR2/ ILDT,ITEST 
C    /LDGBER/ CRT050K 
C    /LDGIMC/ AGEIM,WAIVERF, WAIVERP, NONCOM 
C    /LDGOBD/ HIGH3,HIGH4 
C    /LDGSIM/ IMFLG_NP,IMPN_NP 
C    /MAXIMA/ MAXYRS 
C 
C  Output on return: 
C 
C    function: IMPROG 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C  AIM1          I    Age index: IVAGEIM(JAN_SEP). 
C  AIM2          I    Age index: IVAGEIM(OCT_DEC). 
C  AVEREF        R    Combined Normal/High emission rate. 
C  BEGIN_DATE    I    Integer whose value is 1, indicating the I/M program start. 
C  CIMPN         I    Local I/M program number. 
C  DFACTn        R    Deterioration factor that reduces I/M benefits as a 
C                     function of time since testing stopped for model year 
C                     segment n vehicles. 
C  END_DATE      I    Integer whose value is 2, indicating the I/M program end. 
C  HI3           R    Fraction of high emitters in the fleet with OBD. 
C  HIFRCN        R    Fraction of high emitters in the fleet. 
C  HIGH          I    Integer whose value is 2, indicating high emitter class. 
C  HIGHAVE       R    Local high emission rate. 
C  IMCRD         R    Fractional reduction in basic emissions from I/M (always negative). 
C  IMPGM         I    I/M program number. The user may enter more than one 
C                     program for each run. Each program is specified for 
C                     a range of vehicle types, calender and model years. 
C  IMPN          I    Array with the I/M program number for each segment (jan/sep or oct/dec). 
C  IMREMAIN      R    The fraction of the fleet average emissions that remain after I/M. 
C  ISEG          I    Integer index of the model year segment. 
C  IMLOG_NP      I    RSD and COIM corrections do not apply unless IMLOG_NP>0. 
C  IVAGEIM       I    integer age of a vehicle on the date of its previous IM test. 
C  IY            I    Age looping variable. 
C  JAN_SEP       I    Integer whose value is 1, indicating the January to September interval. 
C  JTEST         i    Local value of ITEST. 
C  NO _IM_PRGRM  I    Integer whose value is 0, indicating that there is no I/M program. 
C  NORMAL        I    Integer whose value is 1, indicating normal emitter class. 
C  NORMAVE       R    Local normal emission rate. 
C  NOT_REQUIRED  I    Integer whose value is 0, indicating that there is no I/M program. 
C  OCT_DEC       I    Integer whose value is 2, indicating the October to December interval. 
C  REDPROG       R    Scaling factor applied to IMPROG at the end of the I/M 
C                     calculations that appropriately reduces January and July 
C                     results for the case when the I/M program start year equals  
C                     calender year of evaluation. 
C  REPAIR_X      R    Holding variable for AR94B. 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LDGBER.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'LDGSIM.I' 
      INCLUDE 'MAXIMA.I' 
C 
C     Declare external functions. 
C 
      INTEGER FINDIM 
      REAL    EF_CASE 
      REAL    EF81_93 
      REAL    IDRATE 
      REAL    IMTEST 
      REAL    EF_HI_FRCN 
C 
C     Declare parameters. 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IROAD 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: ISTD 
      INTEGER, INTENT(INOUT) :: INERR 
      REAL,    INTENT(IN) :: CAT 
      REAL,    INTENT(IN) :: VMTAGE 
C 
C     Declare local variables. 
C 
      REAL    AVEREF(2) 
      REAL    DFACT1 
      REAL    DFACT2 
      REAL    HI3 
      REAL    HIFRCN 
      REAL    HIGHAVE 
      REAL    IDSUM 
      REAL    IMCRD 
      REAL    IMREMAIN 
      REAL    NORMAVE 
      REAL    OBSUM 
      REAL    REDPROG 
      REAL    REPAIR_X 
      INTEGER AIM1 
      INTEGER AIM2 
      INTEGER IDX 
      INTEGER IMPGM 
      INTEGER ISEG,IY,IMLOG_NP 
      INTEGER IVAGEIM(2),IMPN(2),CIMPN 
      INTEGER JTEST 
      INTEGER, SAVE :: NORMAL=1,HIGH=2 
      INTEGER, SAVE :: JAN_SEP=1,OCT_DEC=2 
      INTEGER, SAVE :: BEGIN_DATE=1,END_DATE=2 
      INTEGER, SAVE :: NOT_REQUIRED=0 
      INTEGER, SAVE :: NO_IM_PRGRM=0 
      INTEGER, SAVE :: IEXHAUST = 1 
C 
C 
      IMPROG = 0.0 
      IDX  = (MAXYRS + 1) - JDX 
C 
C Search algorithm to find I/M program numbers for the previous tests 
C taken by both segment 1 and segment 2 vehicles 
C 
      IMPN(JAN_SEP)  = NO_IM_PRGRM 
      IMPN(OCT_DEC)  = NO_IM_PRGRM 
      IMLOG_NP       = NO_IM_PRGRM 
      IVAGEIM(JAN_SEP) = NOT_REQUIRED 
      IVAGEIM(OCT_DEC) = NOT_REQUIRED 
C 
C Search algorithm to find I/M program numbers for the emissions 
C evaluation date by model year and vehicle type 
C 
      CIMPN = FINDIM(ICY,MY,IV,IEXHAUST) 
C 
C     The models years are partitioned because sales of new model year, 
C     MY, vehicles start in Oct. of the previous calendar year, MY-1 
C     and vehicle mileages start in Jan. of a given model year, MY. So, 
C     find appropriate I&M programs for these two different parts of 
C     the calendar year. (See pages 78-80 of M6.IM.001) 
      DO ISEG  = JAN_SEP,OCT_DEC 
C 
C       loop over the I&M programs 
        DO IMPGM = MAXIMPGM, 1, -1 
C 
C         Do not set these variables for evap I/M cases. 
          IF(ITEST(IMPGM).GE.20) CYCLE 
C 
          IF( 
C 
C            check that the current model year vehicle types have the last 
C            I&M date within the user designated I&M program years 
     *       ICYIM(BEGIN_DATE,IMPGM).LE.MY+AGEIM(IMPGM,JDX,ISEG).AND. 
     *       ICYIM(END_DATE,IMPGM).GE.MY+AGEIM(IMPGM,JDX,ISEG)  .AND. 
C 
C            check that the vehicle model year is in the I&M model year range 
     *       MODYR(BEGIN_DATE,IMPGM).LE.MY.AND. 
     *       MODYR(END_DATE,IMPGM).GE.MY  .AND. 
C 
     *       ILDT(IV,IMPGM).EQ.2) THEN 
C 
C            the two model year segments, JAN_SEP & OCT_DEC, can have 
C            different I&M programs depending on the user designated 
C            program years 
             IMPN(ISEG) = IMPGM 
             IVAGEIM(ISEG) = AGEIM(IMPGM,JDX,ISEG) 
             EXIT 
C 
          END IF 
        END DO 
C 
C       vehicle types with model years and/or I&M inspection dates which 
C       are out of range for user designated program (IMPGM) and 
C       the years the program is in effect 
        IF( IVAGEIM(ISEG) .EQ. NOT_REQUIRED ) CYCLE 
C 
C Search algorithm to find I&M program numbers for the selective or non- 
C periodic I/M tests, i.e., COIM or OBD 
C 
        DO IY = IVAGEIM(ISEG),JDX-1 
           DO IMPGM = MAXIMPGM, 1, -1 
C 
C            Do not set these variables for evap I/M cases. 
             IF(ITEST(IMPGM).GE.20) CYCLE 
C 
             IF (ICYIM(BEGIN_DATE,IMPGM).LE.MY+IY .AND. 
     *           ICYIM(END_DATE,IMPGM)  .GE.MY+IY .AND. 
C 
     *           MODYR(BEGIN_DATE,IMPGM).LE.MY   .AND. 
     *           MODYR(END_DATE,IMPGM)  .GE.MY   .AND. 
C 
     *           ILDT(IV,IMPGM).EQ.2) THEN 
C 
C                non-periodic/COIM or OBD program number holder by 
C                by year and month intervals, ISEG. 
                 IMPN_NP(ISEG,IY+1) = IMPGM 
                 IMPN_NP(ISEG,1) = IMPGM 
                 EXIT 
              END IF 
           END DO 
           IF(IMPN_NP(ISEG,IY+1).GT.NO_IM_PRGRM) 
C 
C  from TESTIM.FOR : 
C             IMFLG_NP = 0 no COIM or OBD 
C                        1 RSD only 
C                        2 COIM only 
C                        3 both RSD and OBD 
C 
     *       IMLOG_NP = IMLOG_NP + IMFLG_NP(IMPN_NP(ISEG,IY+1)) 
        END DO 
      END DO 
C 
      AIM1 = IVAGEIM(JAN_SEP) 
      AIM2 = IVAGEIM(OCT_DEC) 
C 
C Calculate the benefits deterioration factor in the event I&M testing 
C ceases for a given model year. I&M testing stops if the end year for 
C the program is reached before the start date of any new program or if 
C the model year moves out of an I&M moving window. 
C 
      DFACT1 = 1.0 
      DFACT2 = 1.0 
C 
      IF (IMPN(JAN_SEP).NE.NO_IM_PRGRM) THEN 
        IF (CIMPN.EQ.NO_IM_PRGRM .AND. AIM1.GE.1 ) THEN 
           DFACT1 = 1.0 - 0.33*REAL(JDX-AIM1) 
        END IF 
        IF (JDX.GT.EXEMPTAGE(IMPN(JAN_SEP))) 
     *     DFACT1 = 1.0 - 0.33*REAL(JDX-EXEMPTAGE(IMPN(JAN_SEP))) 
      END IF 
C 
      IF (IMPN(OCT_DEC).NE.NO_IM_PRGRM) THEN 
        IF (CIMPN.EQ.NO_IM_PRGRM .AND. AIM2.GE.1 ) THEN 
           DFACT2 = 1.0 - 0.33*REAL(JDX-AIM2) 
        END IF 
        IF (JDX.GT.EXEMPTAGE(IMPN(OCT_DEC))) 
     *     DFACT2 = 1.0 - 0.33*REAL(JDX-EXEMPTAGE(IMPN(OCT_DEC))) 
      END IF 
C 
      IF (DFACT1.LE.0.0) DFACT1 = 0.0 
      IF (DFACT2.LE.0.0) DFACT2 = 0.0 
C 
C In general there are two different values of AGEIM correspondinding to 
C each model year. The I/M benefits must be calculated separately for each 
C of the two ages and then weighted together.(See pages 78-84 of M6.IM.001) 
C 
        IMPROG = DFACT1 * FRN(IV) 
     *     * IMTEST(MY,ACTUAL_MY,JDX,REAL(AIM1),IP,IV,ISR,IFDS, 
     *              IMPN(JAN_SEP),ISTD,ICY,IDX,IH,IROAD,CAT,INERR) 
     *     + DFACT2 * FRC(IV) 
     *     * IMTEST(MY,ACTUAL_MY,JDX,REAL(AIM2),IP,IV,ISR,IFDS, 
     *              IMPN(OCT_DEC),ISTD,ICY,IDX,IH,IROAD,CAT,INERR) 
C 
C Apply the reduction factor REDPROG to fix the cases of 
C program start year = calender  year (annual testing) 
C as well as program start year = cal yr +1 (for biennial 
C testing) by scaling down IMPROG. The reduction factor 
C corrects for the unwanted Oct-Jan term in these cases 
C 
C     Initialize REDPROG 
      REDPROG=1. 
C 
C     Apply REDPROG only to exhaust I/M cases (ITEST<20). 
      JTEST=20 
      DO IMPGM = MAXIMPGM, 1, -1 
        IF(IMPGM.EQ.CIMPN) JTEST=ITEST(IMPGM) 
      END DO 
      IF(JTEST.LT.20) THEN 
C 
C Consider the case where program start year = calendar year 
C for both annual and biennial testing frequencies 
C For evaluation month January we expect no I/M credit 
C while for eval month = July we take 2/3 of the I/M credit 
C because we only want six months worth of I/M instead 
C of nine months 
C 
        IF (ICYIM(BEGIN_DATE,CIMPN).EQ.ICY) THEN 
          IF (MEVAL.EQ.1) THEN 
            REDPROG=0.0             !January 
          ELSE IF (MEVAL.EQ.7) THEN 
            REDPROG=0.67            !July 
          ENDIF 
        ENDIF 
C 
C For CY=Prog Start year + 1 we scale only the biennial credits 
C by 4/5 for january (because we want 12 months not 15) 
C and 6/7 for july (because we want 18 months worth  
C instead of 21)  
C 
        IF (ICYIM(BEGIN_DATE,CIMPN)+1.EQ.ICY) THEN 
          IF (IFREQ(CIMPN).EQ.2) THEN 
            IF (MEVAL.EQ.1) THEN 
              REDPROG=0.8             !January 
            ELSE IF (MEVAL.EQ.7) THEN 
              REDPROG=0.8571          !July 
            ENDIF 
          ENDIF 
        ENDIF 
C 
C Apply the reduction factor (REDPROG) to IMPROG for this I/M program 
        IMPROG=IMPROG*REDPROG 
C 
      ENDIF  !ITEST conditional 
C 
C Calculate high emitter fraction and correct 1994 and later model years 
C emissions for OBD if applicable. 
C 
      HIFRCN = EF_HI_FRCN(MY,IP,IV,VMTAGE,IFDS,ISTD,INERR) 
      HI3    = HIGH3(JDX,IP,IV) 
C 
C 
C Calculate the combined Normal/High base emission rate. 
C 
      IF(MY.GE.1981 .AND. MY.LE.1993) THEN      !Tier 0 Case 
        NORMAVE   = EF81_93(MY,IP,IV,VMTAGE,IFDS,ISR,1) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *                IH,IROAD,CAT,NORMAVE) 
        HIGHAVE   = EF81_93(MY,IP,IV,VMTAGE,IFDS,ISR,2) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,2, 
     *                IH,IROAD,CAT,HIGHAVE) 
        AVEREF(ISR) = NORMAVE*(1.0-HIFRCN) + HIGHAVE*HIFRCN 
C 
      ELSEIF(MY.EQ.1994 .OR. MY.EQ.1995) THEN   !Non-OBD Case 
        NORMAVE = EF_CASE(MY,JDX,IP,IV,ISR,VMTAGE,ISTD,NORMAL,IFDS) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *                IH,IROAD,CAT,NORMAVE) 
        HIGHAVE = EF_CASE(MY,JDX,IP,IV,ISR,VMTAGE,ISTD,HIGH,IFDS) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,2, 
     *                IH,IROAD,CAT,HIGHAVE) 
        AVEREF(ISR)  = NORMAVE*(1.0-HIFRCN) + HIGHAVE*HIFRCN 
C 
      ELSEIF(MY.GE.1996) THEN                   !OBD Case 
        NORMAVE = EF_CASE(MY,JDX,IP,IV,ISR,VMTAGE,ISTD,NORMAL,IFDS) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *                IH,IROAD,CAT,NORMAVE) 
        HIGHAVE = EF_CASE(MY,JDX,IP,IV,ISR,VMTAGE,ISTD,HIGH,IFDS) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,2, 
     *                IH,IROAD,CAT,HIGHAVE) 
        REPAIR_X = AR94B 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,1, 
     *                IH,IROAD,CAT,REPAIR_X) 
        AVEREF(ISR)  = NORMAVE*(1.0-HIFRCN) + HIGHAVE*HI3 + 
     *                 REPAIR_X*(HIFRCN-HI3) 
C 
      ELSE                                      !Pre-1981 case should not occur. 
        NORMAVE = 0.0 
        HIGHAVE = 0.0 
        AVEREF(ISR) = 1.0 
      END IF 
C 
C     Allows for 1996+ Exhaust I/M tests.  It ratios the Exhaust test 
C     Identification rate to the OBD I/M Identification rate. 
C 
      IF(IMPN(JAN_SEP).GT.0) THEN 
       IF(MY.GE.1996 .AND. ITEST(IMPN(JAN_SEP)).LT.12) THEN    !Exhaust I/M After 1996 MY. 
C 
        IF(WMIL.EQ.0.0 .OR. OBDIM.EQ.0.0) THEN 
          IDSUM = 1.0 
        ELSE 
          IDSUM=IDRATE(MY,JDX,REAL(AIM1),IP,IV,ISR,IFDS,ISTD,   !Exhaust I/M ID Rate 
     *                 IMPN(JAN_SEP),INERR) 
          OBSUM= WMIL * OBDIM                                   !OBD I/M ID Rate 
          IDSUM = IDSUM / OBSUM 
 
          IF(IDSUM.GT.1.0) IDSUM = 1.0 
          IF(IDSUM.LT.0.0) IDSUM = 0.0 
        ENDIF 
C 
        IMPROG = IMPROG * IDSUM 
       ENDIF 
      ENDIF 
C 
C     IMREMAIN is the fraction of the fleet average emissions that remain 
C     after the I/M credit (IMCRD) has been subtracted.  It is multiplied 
C     times the basic exhaust emission rate in HRLOOP. 
C     IMCRD (the emission benefit) is a negative number. 
C 
      IF(AVEREF(ISR).LE.0.0) THEN 
        IMCRD = 0.0 
      ELSE 
        IMCRD = IMPROG / AVEREF(ISR) 
      ENDIF 
      IF(MY.LE.1980) IMCRD = 0.0  !Pre-1981 case should not occur. 
C 
      IMREMAIN = 1.0 + IMCRD 
      IMPROG = IMREMAIN 
C 
      RETURN 
      END 
