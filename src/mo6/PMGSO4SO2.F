      SUBROUTINE PMGSO4SO2(ICY,INERR) 
C 
C  This routine computes the gasoline SO4 and SO2 based on the PART5 model. 
C  The output from this routine assumes 150 ppm sulfur for non-cat and  
C  1 ppm for Cat so that the SO4 can be subtracted from the EMFAC-based EFs. 
C 
C  Called by: PMCALX. 
C 
C Changes: 
C 
C 17 Apr 01 AIR Task 11: Changed code to handle by model year 
C           sulfur correction for SO4. 
C 13 Feb 01 AIR Task 11: New MOBILE6 routine. 
C 13 Feb 01 AIR Task 11: Updated to use MOBILE6 vehicle classes, 
C           fuel economies, and tech fractiosn. Also added fuel 
C           sulfur correction to calculations. 
C 17 Nov 00 AIR Task 06: Removed Dead Code. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by PMCALX 
C 
C  Calls PMTAMPER, PMGASSL, PMFRACA 
C 
C  Input on call: 
C    parameter list: ICY 
C 
C  Common blocks: see Includes below 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  B_PPM     R   Base gasoline fuel sulfur, in ppm 
C  I         I   Index into road type specific average speeds 
C  IDX       I   IDX=(MAXYRS+1)-JDX 
C  IH        I   Hour 
C  IROAD     I   Road Type 
C  IV        I   Vehicle Class 
C  JDX       I   Age 
C  MY        I   Model Year 
C  SO2       R   Gaseous sulfur emission factor 
C  SO4       R   Direct sulfate emission factor for fleet average gasoline vehicle 
C  SUMSO2IH  R   Temporary summation of SO2 
C  SUMSO2JDX R   Temporary summation of SO2 
C  SUMSO4IH  R   Temporary summation of SO4 
C  SUMSO4JDX R   Temporary summation of SO4 
C  XSO2      R   SO2, temporary holding variable 
C  XSO4      R   Direct sulfate (S04), temporary holding variable 
C  SULF      R   Non-cat (1) and Cat (2) EMFAC BER sulfur levels. 
C  SLFNGV    R   NGV sulfur correction factor 
C  ZEVFRAC   R   Fleet ZEV fraction 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'NGV.I'     ! NGVFLG1, NGV_START, NGVPHIN, NGVFRAC 
      INCLUDE 'PART1.I'   ! BERSULF 
      INCLUDE 'PART2.I'   ! BSWGTR,CTLFRC,FECONM,FRACA,GSWGTR,RMIS,SULFCF,TAMFRC 
      INCLUDE 'PART3.I'   ! EFPM, DBPM, VPM 
      INCLUDE 'SPEED9.I'  ! FVMT, HVMT, SVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PCTLEV 
      REAL, EXTERNAL :: PMFLECON 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare local variables/arrays. 
C 
      INTEGER I 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD      
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      INTEGER NGVPT 
      REAL    B_PPM 
      REAL    GSULFCF 
      REAL    SLFNGV 
      REAL    SO2 
      REAL    SO4 
      REAL    SUMSO2IH 
      REAL    SUMSO2JDX 
      REAL    SUMSO4IH 
      REAL    SUMSO4JDX 
      REAL    XSO2 
      REAL    XSO4 
      REAL    ZEVFRAC 
C 
C 
C  The PART5 gasoline sulfate emission factors are based on 340 ppm sulfur (0.034 wt%). 
C 
      B_PPM = 340.  
      BSWGTR = B_PPM / 10000. 
C 
C  Compute the sulfur correction for NGV. Assume the NG sulfur 
C  level to be 0.1 ppm. 
C 
      SLFNGV = 0.1 / B_PPM  !NGV 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  If not gasoline fueled, skip 
C 
      IF(VVGASMC(IV).EQ.0) CYCLE 
C 
C  Loop for all ages 
C 
        DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
          IDX=(MAXYRS+1)-JDX 
C 
C  If travel fraction for this age/vehicle is 0, skip. 
C 
          IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
          MY=ICY-JDX+1 
C 
C  Convert the BER associated fuel sulfur to wt%. (Base is 340 ppm, per PART5, 
C  and is assigned to BERSULF in the BERPAR subroutine.) 
C 
      GSWGTR = BERSULF(JDX,IV) / 10000. 
C 
C  If we assume that sulfate emission are zero when fuel sulfur is zero, then 
C  account for different levels by applying the following sulfur correction factor: 
C 
      GSULFCF = GSWGTR/BSWGTR 
C 
C  Jan-31-1994 @ CSC-bsg Request 390 New Sulfate equations Get the 
C  misfueling rates: 
C 
      CALL PMTAMPER(ICY, MY, IV) 
C 
C  Feb-10-1994 @ CSC-bsg Request 390 Modified calculation of CATFRC to  
C  include RMIS and remove RMIS from sulfate calculation. 
C  Get the catalyst fraction of the gasoline vehicle population. 
C  (Use MOBILE6 computed values instead of PART5's FRACA.) 
C 
      CALL PMFRACA(MY,IV,INERR) 
       
      CTLFRC = FRACA(1) + FRACA(2) + FRACA(3) + FRACA(4) * (1.-RMIS) 
C 
C  For each Vehicle Type get the corresponding fuel economy 
C  (PMFLECON uses MOBILE6 values) 
C 
      FECONM = PMFLECON(MY,IV,INERR) 
C 
C  Compute the NGV fraction 
C 
        IF(NGVFLG1.EQ.1) THEN 
           IF(MY.LT.NGV_START) THEN 
              NGVFRAC = 0.0 
           ELSE 
              NGVPT = MY - NGV_START + 1 
              NGVFRAC = NGVPHIN(NGVPT,IV) 
           END IF 
        ELSE 
           NGVFRAC = 0.0 
        END IF 
C 
C  Determine the ZEV fraction for LDGV, LDGT1-LDGT4. 
C 
        ZEVFRAC = PCTLEV(MY,IV,11) 
C 
C  Combine the NGV and ZEV adjustments with the gasoline sulfur correction 
C  to compute the final sulfur correction factor. (Note: This equation form 
C  assumes that the ZEVs are removed first, and the NGVs are from the remainder.) 
C 
      SULFCF = GSULFCF * (1.0 - ZEVFRAC) * (1.0 - NGVFRAC)  !Gas portion 
     *       + SLFNGV * NGVFRAC                             !NG  portion 
C 
C  Initialize 
C 
      SUMSO4JDX=0. 
      SUMSO2JDX=0. 
C 
C  Loop for each hour 
C 
      DO IH=1,MAXIH 
C 
      SUMSO4IH=0. 
      SUMSO2IH=0. 
C 
C  Loop for each road type 
C 
      DO IROAD = 1,4 
C 
C  Calculate the direct sulfate (SO4) and SO2. Freeways, Arterials & Collectors 
C  loop over speeds from 2.5 to 65 mph. Locals and Ramps have only one speed. 
C  Programming Note: many of the variables used by PMGASSL are passed via PART2.I. 
C 
           SELECT CASE (IROAD) 
C 
           CASE(1,2)       ! Freeways, Arterials and Collectors 
C 
            SO4=0. 
            SO2=0. 
C 
            DO I=1,14 
              CALL PMGASSL(AVESPD(I,IROAD),XSO4,XSO2) 
              SO4 = SO4 + XSO4 * SVMT(IROAD,I,IH) 
              SO2 = SO2 + XSO2 * SVMT(IROAD,I,IH) 
            END DO 
C 
           CASE(3)       ! Local Roadways 
C 
              CALL PMGASSL(AVSPL,SO4,SO2) 
C 
           CASE(4)        !Ramp Roadways 
C 
              CALL PMGASSL(AVSPR,SO4,SO2) 
C 
           END SELECT 
C 
C  Build the database arrays. 
C 
            DBPM(1,IROAD,IH,JDX,IV) = SO4 
            DBPM(6,IROAD,IH,JDX,IV) = SO2 
C 
            SUMSO4IH = SUMSO4IH + SO4 * FVMT(IROAD,IH,IV)   !SO4 
            SUMSO2IH = SUMSO2IH + SO2 * FVMT(IROAD,IH,IV)   !SO2 
C   
            END DO  !IROAD 
C 
            SUMSO4JDX = SUMSO4JDX + SUMSO4IH * HVMT(IH)   !SO4 
            SUMSO2JDX = SUMSO2JDX + SUMSO2IH * HVMT(IH)   !SO2 
C 
          END DO  !IH 
C 
          EFPM(1,IV)=EFPM(1,IV)+SUMSO4JDX*TF(IDX,IV)   !SO4 
          EFPM(6,IV)=EFPM(6,IV)+SUMSO2JDX*TF(IDX,IV)   !SO2 
C 
        END DO  !JDX 
C 
        VPM(1)=VPM(1)+EFPM(1,IV)*VMTMIX(IV)   !SO4 
        VPM(6)=VPM(6)+EFPM(6,IV)*VMTMIX(IV)   !SO2 
C 
      END DO  !IV 
C 
      RETURN 
      END 
