      SUBROUTINE PMGASPM(ICY) 
C 
C  Computes GASPM, which is EMFAC_PM - SO4 or RCP+SOF. 
C  RCP and SOF are dependent on particle size but 
C  are assumed to be independent of fuel sulfur content. 
C 
C 
C  Changes : 
C 
C  24 Aug 01 AIR Task 23b: Updated to handle EPA-based ZML/DR 
C            results. 
C  17 Apr 01 AIR Task 11: Changed code to handle by model year 
C            sulfur correction for SO4 and prevent negative GASPM. 
C  23 Feb 01 AIR Task 11: New Module 
C 
C  Called by PMCALX 
C 
C  Calls P3WPOX 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  ----------------------------------------------------- 
C  CAT         R   Catalyst fraction 
C  IDX         I   IDX=(MAXYRS+1)-JDX 
C  IH          I   Hour 
C  IROAD       I   Road Type 
C  IV          I    Vehicle Loop 
C  IV          I   Vehicle Class 
C  JDX         I   Age 
C  MY          I   Model year  
C  P3W         R   Fraction of 3-way catalysts 
C  POX         R   Fraction of oxidation catalysts 
C  PSIZECF     R   Temporary particulate size correction factor 
C  RCPSOF      R    RCP+SOF = EMFAC_PM - SO4 
C  RCPSOF      R   RCP+SOF (EMFAC_PM - SO4) 
C  SO2         R   Temporary SO2 
C  SO4         R   Temporary SO4 
C  SUMIH       R   Temporary summation over IH 
C  SUMJDX      R   Temporary summation over JDX 
C  ZEVADJ      R   Adjustment for ZEVs 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I'  ! IOUREP,IOUERR 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART3.I'   ! DBPM, EFPM, VPM 
      INCLUDE 'PSIZE.I'   ! PSL, PSNLCT 
      INCLUDE 'SPEED9.I'  ! HVMT, FVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PCTLEV 
C 
C  Declare local variables/arrays 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      REAL    CAT 
      REAL    P3W 
      REAL    POX 
      REAL    PSIZECF 
      REAL    RCPSOF 
      REAL    SUMIH 
      REAL    SUMJDX 
      REAL    ZEVADJ 
C 
C 
C  Re-iInitialize to hold final value. 
C 
      VPM(4)=0. 
C 
C  GASPM is calculated for gas vehicles only. 
C 
      DO IV=1,MAXVEH 
C 
        IF(VVGASMC(IV).EQ.0) CYCLE 
C 
C  Re-initialize to hold final values. 
C 
      EFPM(4,IV)=0. 
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
C  If the travel fraction for this age/vehicle combination is zero, skip. 
C 
      IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
      MY=ICY-JDX+1 
C 
C  Compute the PART6 particle size adjustment factor, as follows: 
C    If MY1974- and CY1990- use single leaded value from PART6, 
C    otherwise use cat/nocat weighted unleaded value from PART6 
C 
      IF(MY.LE.1974.AND.ICY.LE.1990) THEN 
        PSIZECF=PSL 
      ELSE 
        CALL P3WPOX(MY,IV,P3W,POX) 
        CAT=P3W+POX 
        PSIZECF=(1.0-CAT)*PSNL+CAT*PSNLCT 
      ENDIF 
C 
C  Determine the ZEV adjustment for LDGV, LDGT1-LDGT4. 
C  (NGVs are assummed to emit GASPM like their gasoline counterparts, 
C   so no other adjustments are needed.) 
C 
      ZEVADJ = 1.0 - PCTLEV(MY,IV,11) 
C 
C  Initialize 
C 
      SUMJDX=0. 
C 
C  Loop for hour 
C 
      DO IH=1,MAXIH 
C 
      SUMIH=0. 
C 
C  Loop for road type 
C 
      DO IROAD=1,4 
C 
C  Compute GASPM, which is actually RCP+SOF. The SO4 being subtracted 
C  here is from PMGSO4SO2 and is based on the BER associated fuel 
C  sulfur level, not the in-use sulfur level. 
C 
        RCPSOF = DBPM(10,IROAD,IH,JDX,IV) - DBPM(1,IROAD,IH,JDX,IV) 
C 
C  Produce an error message when negative RCPSOF values occur.  
C  Rest such values to zero. 
C 
        IF(RCPSOF.LT.0.) THEN 
CDFK        WRITE(IOUOUT,10) RCPSOF,ICY,MY,IV 
CDFK        WRITE(IOUERR,10) RCPSOF,ICY,MY,IV 
CDFK  10    FORMAT('*** Error:  a negative Gas RCPSOF value of ',F8.3,/ 
CDFK     &         '    occurred for CY=',I4,' MY=,',I4,' and IV=',I2) 
        RCPSOF=0.0 
        END IF 
C 
C  Adjust for particle size 
C 
        RCPSOF = RCPSOF * PSIZECF 
C 
C  Adjust for ZEVs 
C 
        RCPSOF = RCPSOF * ZEVADJ 
C 
C  Store in arrays and compute composites.      
C 
        DBPM(4,IROAD,IH,JDX,IV) = RCPSOF 
C 
        SUMIH = SUMIH + RCPSOF * FVMT(IROAD,IH,IV) 
C 
        END DO  !IROAD 
C 
        SUMJDX = SUMJDX + SUMIH * HVMT(IH) 
C 
        END DO  !IH 
C 
        EFPM(4,IV)= EFPM(4,IV) + SUMJDX * TF(IDX,IV) 
C 
       END DO ! JDX 
C 
       VPM(4)= VPM(4) + EFPM(4,IV) * VMTMIX(IV) 
C 
      END DO ! IV 
C 
      RETURN 
      END 
