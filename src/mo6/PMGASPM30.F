      SUBROUTINE PMGASPM30(ICY,INERR) 
C 
C  Computes GASPM, which is RCP+SOF using the PART5 carbon method. 
C  The values are computed aT PM30 and are independent of fuel 
C  sulfur content. 
C 
C 
C  Changes : 
C 
C  23 Oct 01 AIR Task 25: New Routine. 
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
C  IV          I   Vehicle Class 
C  JDX         I   Age 
C  MY          I   Model year  
C  RCPSOF      R   RCP+SOF 
C  ZEVADJ      R   Adjustment for ZEVs 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART3.I'   ! DBPM, EFPM, VPM 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: CEFPM30 
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
C  Loop for hour 
C 
      DO IH=1,MAXIH 
C 
C  Loop for road type 
C 
      DO IROAD=1,4 
C 
C  Compute GASPM, which is actually RCP+SOF, using the PART5 carbon method. 
C  CEFPM30 returns with values at PM30. 
C 
        RCPSOF = CEFPM30(IV,MY,ICY,INERR)  ! PM30 
C 
C  Adjust for ZEVs. 
C 
        RCPSOF = RCPSOF * ZEVADJ 
C 
C  Store in arrays.      
C 
        DBPM(4,IROAD,IH,JDX,IV) = RCPSOF ! PM30 
C 
        END DO  !IROAD 
C 
        END DO  !IH 
C 
       END DO ! JDX 
C 
      END DO ! IV 
C 
      RETURN 
      END 
