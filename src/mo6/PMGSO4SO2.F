      SUBROUTINE PMGSO4SO2(ICY,INERR) 
C 
C  This routine computes the gasoline SO4 and SO2 based on the PART5 model. 
C  It uses the BERs associated sulfur level from PMEMFAC so that the  
C  associated SO4 can be subtracted from the EMFAC-based EFs. 
C 
C  Called by: PMCALX. 
C 
C Changes: 
C 
C 29 Aug 01 AIR Task 23c: Updated the PART5 sulfur values, thereby changing 
C           the base sulfur level from 340 ppm to 148 ppm. 
C 03 May 01 AIR Task 12: Apply associated BER fuel sulfur. 
C 30 Apr 01 AIR Task 12: Code for Tier2/2007HD Rule. 
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
C  B_PPM     R   Base gasoline fuel sulfur content, in ppm 
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
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS5.I'  ! TIER2_FLAG,HDR2007 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART1.I'   ! BERSULF 
      INCLUDE 'PART2.I'   ! BSWGTR,CTLFRC,FECONM,FRACA,GSWGTR,RMIS,SULFCF,TAMFRC 
      INCLUDE 'PART3.I'   ! EFPM, DBPM, VPM 
      INCLUDE 'SPEED9.I'  ! FVMT, HVMT, SVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD      
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      REAL    B_PPM 
      REAL    SO2 
      REAL    SO4 
      REAL    SUMSO2IH 
      REAL    SUMSO2JDX 
      REAL    SUMSO4IH 
      REAL    SUMSO4JDX 
C 
C 
C  The basic PART5 SO4 values were measured on 148 ppm fuel 
C  (Updated 29Aug01 from 340 ppm). 
C  So use that level as the base, in wt%. 
C 
      B_PPM  = 148. 
      BSWGTR = B_PPM / 10000. 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  If not gasoline fueled, skip 
C 
      IF(VVGASMC(IV).EQ.0) CYCLE 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
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
C  In order for the SO4 to be subtracted from the EMFAC-based or Tier2/ 
C  2007HD Rule EFs, the sulfur levels must agree. Therefore, use the 
C  BER associated sulfur level from the PMEMFAC routine and convert to wt% 
C 
      GSWGTR = BERSULF(JDX,IV) / 10000. 
C 
C  The gasoline sulfate emission factors are based on 148 ppm sulfur (0.0148 wt%). 
C  If we assume that sulfate emission are zero when fuel sulfur is zero, then 
C  account for different levels by applying the following sulfur correction factor: 
C 
      SULFCF = 1.0 + (GSWGTR - BSWGTR)/BSWGTR 
C 
C  Get the misfueling rates: 
C 
      CALL PMTAMPER(ICY, MY, IV) 
C 
C  Get the catalyst fraction of the gasoline vehicle population. 
C  (Use MOBILE6 computed values instead of PART5's FRACA.) 
C 
      CALL PMFRACA(MY,IV,INERR) 
       
      CTLFRC = FRACA(1) + FRACA(2) + FRACA(3) + FRACA(4) * (1.-RMIS) 
C 
C  For each Vehicle Type get the corresponding fuel economy 
C  (PMFLECON uses MOBILE6 values) 
C 
      FECONM = PMFLECON(MY,IV) 
C 
C  Compute the SO4 and SO2. 
C  Programming Note: Many of the values used in PMGASSL are passed 
C  in via PART2.I. 
C 
      CALL PMGASSL(SO4,SO2) 
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
