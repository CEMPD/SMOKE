      SUBROUTINE PMGASPM(ICY,INERR) 
C 
C  Computes GASPM, which is RCP+SOF using the PART5 carbon method. 
C  RCP and SOF are dependent on particle size but 
C  are assumed to be independent of fuel sulfur content. 
C  Routine also adjusts for BER limits. 
C 
C 
C  Changes : 
C 
C  18 Apr 01 AIR Task 11: Reverted routine back to PART5 method. 
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
C  IDX         I   IDX=(MAXYRS+1)-JDX 
C  IH          I   Hour 
C  IROAD       I   Road Type 
C  IV          I    Vehicle Loop 
C  IV          I   Vehicle Class 
C  JDX         I   Age 
C  MY          I   Model year  
C  RCPSOF      R   RCP+SOF 
C  SO2         R   Temporary SO2 
C  SO4         R   Temporary SO4 
C  SUMIH       R   Temporary summation over IH 
C  SUMJDX      R   Temporary summation over JDX 
C  ZEVADJ      R   Adjustment for ZEVs 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART3.I'   ! DBPM, EFPM, VPM, BERADJ 
      INCLUDE 'PSIZE.I'   ! PSL, PSNLCT 
      INCLUDE 'SPEED9.I'  ! HVMT, FVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: CEFPM 
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
      REAL    RCPSOF 
      REAL    SUMIH 
      REAL    SUMJDX 
      REAL    ZEVADJ 
C 
C 
C  GASPM is calculated for gas vehicles only. 
C 
      DO IV=1,MAXVEH 
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
C  If the travel fraction for this age/vehicle combination is zero, skip. 
C 
      IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
      MY=ICY-JDX+1 
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
C  Compute GASPM, which is actually RCP+SOF, using the PART5 carbon method. 
C  CEFPM returns with values already adjusted for particle size. 
C 
        RCPSOF = CEFPM(IV,MY,ICY,INERR) 
C 
C  Adjust for ZEVs and BER Limits. 
C 
        RCPSOF = RCPSOF * ZEVADJ * BERADJ(IROAD,IH,JDX,IV) 
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
