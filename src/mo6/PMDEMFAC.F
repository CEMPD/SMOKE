      SUBROUTINE PMDEMFAC(ICY) 
C 
C  Computes the Total PM (PM = RCP + SOF + SO4) based on the EPA methodology. 
C  The output values are at particle size PM30. 
C 
C  Changes : 
C 
C  16 Oct 01 AIR Task 25: Changed code so only diesels are processed. 
C            Removed VMT limits on diesel vehicles. 
C  11 Sep 01 AIR Task 18: Added HAPS calculations 
C  07 Feb 01 AIR Task 15-23c: Fixed bug where BERSULF wasn't being initialized 
C                          to the HDR 2007 rule properly. 
C  28 Aug 01 AIR Task 23c: Implement DR1, DR2 and PMDRAGE 
C  24 Aug 01 AIR Task 23b: Implement EPA-based ZMLs/DRs and methodology. 
C  06 Jun 01 AIR Task 17: Added CA diesel aromatic control program code. 
C  23 May 01 AIR Task 17: Made No-Speed-Corrections as default. Removed 4-Mode 
C            ZMLs/DRs and made MY start in 1950- in stead of 1975-. 
C  04 May 01 AIR Task 12: Tier2/2007 Rule 
C  17 Apr 01 AIR Task 11: Changed code to handle by model year 
C            sulfur correction for gasoline fueled vehicles and 
C            expanded FVMT dimensions. 
C  01 Feb 01 AIR Task 11: New Module 
C 
C     Called by:  PMCALX 
C 
C  Input on call: 
C 
C    Parameter list: ICY 
C 
C    Common blocks:  See Includes below 
C 
C  Output on return:  
C    /PART3.I/ EFPM,VPM,DBPM 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX4.I'     
      INCLUDE 'ATOX5.I'     
      INCLUDE 'CUMCOM.I'   ! CUMMIL 
      INCLUDE 'IVPCOM.I'   ! IVPTRT 
      INCLUDE 'IVTYPE.I'   ! MAXVEH,VTHDG2B,VTHDG8B,VTHDD2B,VTHDD8B,VTGASBUS,VTSCHOOL,VVLDDSL 
      INCLUDE 'MAXIMA.I'   ! MAXYRS, MAXIH 
      INCLUDE 'MYRCAL.I'   ! TF 
      INCLUDE 'PART1.I'    ! PRTCHK, D_PPM, BERSULF, BERMAX 
      INCLUDE 'PART3.I'    ! EFPM, DBPM, VPM 
      INCLUDE 'PART4.I'    ! PMDZML, PMDDR1, PMDDR2, BSULFLVL 
      INCLUDE 'PART5.I'    ! PMSCFC,PMSPDFLG 
      INCLUDE 'SPEED9.I'   ! FVMT, HVMT 
      INCLUDE 'SULFUR.I'   ! S_PPM 
      INCLUDE 'VMXCOM.I'   ! VMTMIX 
C 
C  Declare external functions. 
C 
      INTEGER, EXTERNAL :: IEVPTR 
      REAL   , EXTERNAL :: HDCPIC 
C 
C  Declare parameter list 
C 
      INTEGER ICY     ! Calendar Year 
C 
C  Declare local variables/arrays 
C 
      INTEGER IDX      ! IDX=(MAXYRS+1)-JDX 
      INTEGER IH       ! Hour index 
      INTEGER IROAD    ! Road index 
      INTEGER IV       ! Vehicle Class 
      INTEGER JDX      ! Age 
      INTEGER MY       ! Model year 
      INTEGER MYPM     ! model year pointer into PM arrays, 1950-=1, 2020+=71 
      REAL    BODOM    ! Base odometer mileage 
      REAL    KODOM    ! Mileage after kink 
      REAL    PMBER    ! Basic emmission rates, g/mi 
      REAL    SUM      ! Temporary summation 
      REAL    SUMIH    ! Temporary summation over IH 
C 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip if not a diesel vehicle 
C 
      IF(VVDSL(IV).EQ.0) CYCLE 
C 
C  Skip if VMTMIX is 0. 
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
C  Determine the base and kink mileage, in 10K. 
C 
      IF(JDX.LE.PMDRAGE(IV)) THEN 
        BODOM=CUMMIL(JDX,IV)/10000. 
        KODOM=0. 
      ELSE 
        IF(PMDRAGE(IV).EQ.0) THEN 
          BODOM=0. 
        ELSE 
          BODOM=CUMMIL(PMDRAGE(IV),IV)/10000. 
        END IF 
        KODOM=CUMMIL(JDX,IV)/10000. - BODOM 
      END IF       
C 
C  Determine the model year pointer into the PM_BER arrays. 
C 
      MYPM=MY-1950+1 
      IF(MYPM.LT. 1) MYPM= 1 
      IF(MYPM.GT.71) MYPM=71 
C 
C  COMPUTE THE BERS 
C 
C  Compute the basic emission rate (BER) in g/mi. 
C 
        PMBER = PMDZML(MYPM,IV) 
     *        + PMDDR1(MYPM,IV) * BODOM 
     *        + PMDDR2(MYPM,IV) * KODOM 
C 
C  Apply diesel BER limits, if applicable. Limits are computed by the 
C  PMBERPAR subroutine. If BERMAX is negative, no limit is applied.  
C 
      IF(BERMAX(JDX,IV).GE.0. .AND. 
     *   PMBER.GT.BERMAX(JDX,IV)) PMBER=BERMAX(JDX,IV)      
C 
C  BUILD BY-MODEL-YEAR AND COMPOSITE EMISSION FACTORS 
C 
      SUMIH = 0. 
C 
C  Loop for each hour 
C 
      DO IH = 1,MAXIH 
C 
       SUM=0. 
C 
C  Loop for each road type 
C 
       DO IROAD = 1,4 
C 
C  Build the Total PM here. 
C 
         DBPM(10,IROAD,IH,JDX,IV) = PMBER  ! PM30 
C 
         SUM = SUM + DBPM(10,IROAD,IH,JDX,IV) * FVMT(IROAD,IH,IV) 
C 
         END DO    !IROAD 
C 
         SUMIH = SUMIH + SUM * HVMT(IH) 
C 
        END DO     !IH 
C 
        EFPM(10,IV)= EFPM(10,IV) + SUMIH * TF(IDX,IV) ! PM30 
C 
       END DO ! JDX 
C 
       VPM(10)= VPM(10) + EFPM(10,IV) * VMTMIX(IV) ! PM30 
C 
      END DO ! IV 
C 
      RETURN 
      END 
