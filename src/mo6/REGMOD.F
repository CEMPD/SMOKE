      SUBROUTINE REGMOD(ICY,INERR) 
C 
C  REGMOD adjusts model year mileage/registration arrays. 
C 
C  Called by MOBILE. 
C 
C  Calls CALCUM, CEVMPD, DSFIDX, GETCUM, TFCALX and VEHCOUNT. 
C 
C  Changes: (Last change first) 
C 
C  14 May 01 @EPA-djb Bug 323, Added JDX=1 to the loop which calculates 
C            AGED for motorcycles and diesel transit busses (IVC=15 and 16). 
C  26 Feb 01 @EPA-bag bug315 added checks for division by JANSUM(=0.) 
C  26 Sep 00 AIR Task 07: Bug Fix 188. Added call to GETCUM to initialize 
C            CUMMIL arrays. 
C  07 Jul 00 AIR Task 05: Removed unused ICY in CALCUM call. 
C            Removed unused ICY in TFCALX call. 
C  07 Apr 00 @EPA-djb Changes in use of DSFIDX 
C  25 Feb 00 @DynTel-ddj 1-032, Added GSFVCT(28). 
C  16 Feb 00 @DynTel-MLA 1-011  Removed adjustment to year 12 of MEVMYR 
C            for motorcycles (after consultation with EPA). 
C  22 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C  02 Sep 98 DynTel@RJD 0-632  Set the maximum value of the model 
C            year index to 12 for motorcycles. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C  02 Mar 98 DynTel@RJD 2-204 Introduced call to the CALCUM subroutine to  
C            avoid having the same source code in two places. 
C  02 Mar 98 DynTel@RJD 2-204 (after gjr) Changed calculations for annual  
C            mileage and average accumulation rate when calculating July  
C            emissions. July emissions will be calculated by aging the VMT 
C            mix by 6 months rather than interpolating July emissions from  
C            two consecutive Januaries. 
C  02 Mar 98 DynTel@RJD 2-204 Changed calculations to allow for scrapage in  
C            estimating January 1st model year registrations from July 1st  
C            data. Mobile5b only backdated July 1st data for model year  
C            index 1 to take account of new car sales. Consequently, model  
C            year registration mixes changed abruptly on July 1. The revised 
C            code assumes instead that vehicles are scrapped at a uniform  
C            rate throughout the year. 
C  03 Mar 98 DynTel@RJD 2-204  Changed the name of the array storing the Model 
C            Year Registration fractions from JANMYR to MEVMYR to reflect the 
C            fact that January is no longer the only month of evaluation that 
C            the Mobile model can handle. 
C 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: ICY,INERR 
C 
C    common blocks: 
C    /EVAL/   MEVAL, FRC, FRN 
C    /FLAGS2/ MYMRFG 
C    /FLAGS3/ 
C    /DSFCOM/ 
C    /MAXIMA/ MAXYRS,MCYRS 
C    /MYRSAV/ AMAR,JULMYR 
C    /VVNAMS/ VVHEAVY  (IVTYPE.I) 
C    /VVTYPS/ VTMC     (IVTYPE.I) 
C    /VMXCOM/ VALL 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /MYRCAL/ XMYM,MEVMYR,TFMYM, AGED 
C    /REGISF/ DSFVCT 
C 
C  Local array subscripts: 
C 
C   AGED(25,30) - AGED( JDX , IVD ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  FRYA     R    fraction of a year that is remaining until the next July 1 
C  FRYB     R    fraction of a year that has elapsed since the last July 1 
C  ICY      I    calendar year 
C  INERR    I    error count 
C  IV       I    standard veh class groupings (28) 
C  IVC      I    combined gas/diesel veh class groupings (16) 
C  IVD      I    separate gas/diesel veh class groupings (30) 
C  JANMYR   R    age distribution based on january (rather than july) 
C  JANSUM   R    sum of JANMYR for normalization 
C  JDX      I    vehicle age index; JDX=25-IDX+1 
C  MY       I    see parameter dictionary 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'DSFCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'REGISF.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER   ICY 
      INTEGER   INERR 
      INTEGER   IV 
      INTEGER   IVC 
      INTEGER   IVD 
      INTEGER   JDX 
      INTEGER   MY 
C 
      REAL      DSFIDX   !function 
      REAL      FRYA 
      REAL      FRYB 
      REAL      JANMYR 
      REAL      JANSUM 
C 
      DIMENSION JANMYR(25) 
C 
C  Compute the cumulative mileage arrays (XMYM and CUMMIL) which depend 
C  on evaluation month. 
C 
      CALL CALCUM() 
      CALL GETCUM 
C 
C  Figure the model year mileage accumulation rates, TFMYM 
C  using the fraction of January 1 or July 1 registrations in Ith 
C  year of operation with model year index I (FRN and FRC) calculated in 
C  Subroutine CALCUM. 
C 
      DO IV=1,MAXVEH 
        TFMYM(1,IV)=AMAR(1,IV) 
        IF(VVHEAVY(IV).EQ.1.OR.IV.EQ.VTMC) THEN 
          IF (MEVAL.EQ.1) TFMYM(1,IV)=0.0 
        ENDIF 
        DO JDX=2,MAXYRS 
          TFMYM(JDX,IV)=AMAR(JDX-1,IV)*FRN(IV)+AMAR(JDX,IV)*FRC(IV) 
        END DO 
      END DO !IV 
C 
C  Calculate fraction of a year that has elapsed since July 1st when 
C  vehicles were last registered. Note: this module has only been tested 
C  for MEVAL = 1 and 7. 
C 
      IF(MEVAL.GE.7) THEN  
          FRYA=(19-MEVAL)/12.0 
      ELSE IF (MEVAL.LT.7) THEN 
          FRYA=(7-MEVAL)/12.0 
      END IF  
      FRYB=1.0-FRYA 
C 
C  Compute January 1 or July 1 registration distributions, MEVMYR. 
C 
C  Loop by diesel sales groupings 
C 
      DO IVC=1,16 
C 
C  Combined diesel/gas sales vehicle class groupings are: 
C 
C         IVC                     IVD 
C        -----                   ----- 
C     1 : LDV                 1 : LDGV       17 : LDDV 
C     2 : LDT1                2 : LDGT1      18 : LDDT1 
C     3 : LDT2                3 : LDGT2      19 : LDDT2 
C     4 : LDT3                4 : LDGT3      20 : LDDT3 
C     5 : LDT4                5 : LDGT4      21 : LDDT4 
C     6 : HDV2B               6 : HDGV2B     22 : HDDV2B 
C     7 : HDV3                7 : HDGV3      23 : HDDV3 
C     8 : HDV4                8 : HDGV4      24 : HDDV4 
C     9 : HDV5                9 : HDGV5      25 : HDDV5 
C    10 : HDV6               10 : HDGV6      26 : HDDV6 
C    11 : HDV7               11 : HDGV7      27 : HDDV7 
C    12 : HDV8A              12 : HDGV8A     28 : HDDV8A 
C    13 : HDV8B              13 : HDGV8B     29 : HDDV8B 
C    14 : HDBS (school bus)  14 : HDGBS      30 : HDDBS 
C    15 : HDBT (transit bus) 15 : HDDBT 
C    16 : MC                 16 : MC 
C 
C  Step 1: backup ICY's registrations from 7/1 (July) to 1/1 (Jan) - 
C          LDV/T => my starts 10/1 prior yr => 1/1 sales are 1/3 of 7/1 sales. 
C          HDV & MC => my starts 1/1 => 1/1 sales are 0. 
C 
      IF(MEVAL.EQ.1) THEN 
        JANMYR(1)=JULMYR(1,IVC)/3. 
        IF(IVC.GE.6) JANMYR(1)=0.0 
        JANSUM=JANMYR(1) 
        DO JDX=2,MAXYRS 
          JANMYR(JDX)=JULMYR(JDX,IVC) 
          JANSUM=JANSUM+JANMYR(JDX) 
        END DO  !JDX 
        IF(JANSUM>0.0) THEN 
          DO JDX=1,MAXYRS 
            JANMYR(JDX)=JANMYR(JDX)/JANSUM 
          END DO  !JDX 
        ENDIF 
      ELSE      !Month is July 
        JANSUM=0.0 
        DO JDX=1,MAXYRS 
          JANSUM=JANSUM+JULMYR(JDX,IVC) 
        END DO  !JDX 
        IF(JANSUM>0.0) THEN 
          DO JDX=1,MAXYRS 
            JANMYR(JDX)=JULMYR(JDX,IVC)/JANSUM 
          END DO  !JDX 
        ENDIF 
      ENDIF 
C 
C  Step 2: for each model year of each vehicle type, select and apply the 
C          corresponding GSF to the (gas & diesel combined) July 1 
C          registrations. 
C 
      IF(IVC.LE.14) THEN 
        DO JDX=1,MAXYRS 
          MY=ICY-JDX+1 
          AGED(JDX,IVC)=(1.0-DSFIDX(ICY,MY,IVC)) 
     *                 *JANMYR(JDX) 
          AGED(JDX,IVC+16)=DSFIDX(ICY,MY,IVC) 
     *                 *JANMYR(JDX) 
        END DO !JDX 
C 
C  Motorcycles and diesel transit busses don't have diesel sales fractions 
C 
      ELSE 
        DO JDX=1,MAXYRS 
          AGED(JDX,IVC)=JANMYR(JDX) 
        END DO !JDX 
      ENDIF 
C 
      END DO ! IVC Loop 
C 
C  Step 3: normalize AGED. 
C          the sum DAF will be the fraction of vehicle class 
C          which is that IVD (gas or diesel). 
C 
      DO IVD=1,30 
        DAF(IVD)=0.0 
        DO JDX=1,MAXYRS 
          DAF(IVD)=DAF(IVD)+AGED(JDX,IVD) 
        END DO !JDX 
        IF(DAF(IVD).GT.0.0) THEN 
          DO JDX=1,MAXYRS 
            AGED(JDX,IVD)=AGED(JDX,IVD)/DAF(IVD) 
          END DO !JDX 
        ENDIF 
      END DO !IVD 
C 
C  Determine the vehicle counts by vehicle class (VALL) 
C 
      CALL VEHCOUNT(ICY,FRYA,FRYB) 
C 
C Map the 30 vehicle type registration distributions into the 28 vehicle 
C type registration distributions used by MOBILE6. This will mean combining 
C the LDDT1 & 2 and LDDT3 & 4 by the vehicle count fractions. 
C 
      DO JDX=1,MAXYRS 
        DO IV=1,MAXVEH 
          IF(IV.LE.13) THEN 
           MEVMYR(JDX,IV)=AGED(JDX,IV) 
          ELSEIF(IV.EQ.14) THEN 
            MEVMYR(JDX,IV)=AGED(JDX,17) 
          ELSEIF(IV.EQ.15) THEN 
            IF((VALL(18)+VALL(19)).LE.0.0) THEN 
              MEVMYR(JDX,IV)=0.0 
            ELSE 
              MEVMYR(JDX,IV)=(AGED(JDX,18)*VALL(18) 
     *                   + AGED(JDX,19)*VALL(19)) 
     *                    / (VALL(18)+VALL(19)) 
            ENDIF 
          ELSEIF(IV.GE.16.AND.IV.LE.23) THEN 
            MEVMYR(JDX,IV)=AGED(JDX,IV+6) 
          ELSEIF(IV.EQ.24) THEN 
            MEVMYR(JDX,IV)=AGED(JDX,16) 
          ELSEIF(IV.EQ.25) THEN 
            MEVMYR(JDX,IV)=AGED(JDX,14) 
          ELSEIF(IV.EQ.26) THEN 
            MEVMYR(JDX,IV)=AGED(JDX,15) 
          ELSEIF(IV.EQ.27) THEN 
            MEVMYR(JDX,IV)=AGED(JDX,30) 
          ELSEIF(IV.EQ.28) THEN 
            IF((VALL(20)+VALL(21)).LE.0.0) THEN 
              MEVMYR(JDX,IV)=0.0 
            ELSE 
              MEVMYR(JDX,IV)=(AGED(JDX,20)*VALL(20) 
     *                    + AGED(JDX,21)*VALL(21)) 
     *                    / (VALL(20)+VALL(21)) 
            ENDIF 
          ELSE 
            MEVMYR(JDX,IV)=0.0 
          ENDIF 
        END DO !IV 
      END DO !JDX 
C 
C  Use TFMYM to generate miles per day by model year group from TFMYM. 
C 
      CALL CEVMPD 
C 
C  Use MEVMYR, TFMYM and VALL to generate the vmt mix and travel fractions. 
C 
      CALL TFCALX(INERR) 
C 
      RETURN 
      END 
