      SUBROUTINE PMTEMFAC 
C 
C  Computes the total PM values.  Also computes the PM-based Air Toxics. 
C  The Total PM and AT values are at the particle size selected by the user. 
C 
C  Changes : 
C 
C  04 Jan 02 AIR Task 25: Removed unused variables 
C  28 Nov 01 AIR Task 28: Added ATHEAP for AT dynamic array support. 
C  20 Oct 01 AIR Task 25: New Module 
C 
C     Called by:  PMCALX 
C 
C  Input on call: 
C 
C    Parameter list: None 
C 
C    Common blocks:  See Includes below 
C 
C  Output on return:  
C    /PART3.I/ EFPM,VPM,DBPM 
C 
      USE ATHEAP 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX4.I'     
      INCLUDE 'ATOX5.I'     
      INCLUDE 'IVTYPE.I'   ! MAXVEH 
      INCLUDE 'MAXIMA.I'   ! MAXYRS, MAXIH 
      INCLUDE 'MYRCAL.I'   ! TF 
      INCLUDE 'PART3.I'    ! EFPM, DBPM, VPM 
      INCLUDE 'SPEED9.I'   ! FVMT, HVMT 
      INCLUDE 'VMXCOM.I'   ! VMTMIX 
C 
C  Declare local variables/arrays 
C 
      INTEGER IAT      ! Air toxic pollutant index 
      INTEGER IDX      ! IDX=(MAXYRS+1)-JDX 
      INTEGER IH       ! Hour index 
      INTEGER IMAP     ! Air toxic mapping index 
      INTEGER IROAD    ! Road index 
      INTEGER IV       ! Vehicle Class 
      INTEGER JDX      ! Age 
      REAL    SUM      ! Temporary summation 
      REAL    SUMIH    ! Temporary summation over IH 
C 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
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
C  Build the Total Gas PM here. 
C 
         IF(VVGASMC(IV).EQ.1) 
     *   DBPM(10,IROAD,IH,JDX,IV) = DBPM(4,IROAD,IH,JDX,IV) ! Carbon, at selected particle size 
     *                            + DBPM(1,IROAD,IH,JDX,IV) ! SO4, at in-use sulfur level 
     *                            + DBPM(5,IROAD,IH,JDX,IV) ! Lead, at selected particle size 
C 
C  Build the Total Diesel PM here. 
C 
         IF(VVDSL(IV).EQ.1)  
     *   DBPM(10,IROAD,IH,JDX,IV) = DBPM(2,IROAD,IH,JDX,IV) ! SOF, at selected particle size 
     *                            + DBPM(1,IROAD,IH,JDX,IV) ! SO4, at in-use sulfur level 
     *                            + DBPM(3,IROAD,IH,JDX,IV) ! RCP, at selected particle size 
C 
         SUM = SUM + DBPM(10,IROAD,IH,JDX,IV) * FVMT(IROAD,IH,IV) 
C 
         END DO    !IROAD 
C 
         SUMIH = SUMIH + SUM * HVMT(IH) 
C 
        END DO     !IH 
C 
        EFPM(10,IV)= EFPM(10,IV) + SUMIH * TF(IDX,IV) 
C 
C  Compute Air Toxics based on Total PM. 
C 
        IF(ATHAPS) THEN  ! Build the PM-based exhaust HAPS here 
          DO IMAP=1,NMAP 
            IF(ATHIV(1,4,IMAP,IV)) THEN 
              IAT = MAXIAT + IMAP 
              EFAT(1,IAT,IV) = EFAT(1,IAT,IV) 
     *                       + SUMIH  
     *                       * ATHBER(1,JDX,1,4,IV,IMAP)  
     *                       * TF(IDX,IV) ! Ratios are already in mg/mi 
            END IF ! ATHIV 
          END DO ! IMAP 
        END IF ! ATHAPS 
C 
       END DO ! JDX 
C 
       VPM(10)= VPM(10) + EFPM(10,IV) * VMTMIX(IV) 
C 
       IF(ATHAPS) THEN ! Build the "All Vehicle" exhaust HAPS here 
         DO IMAP=1,NMAP 
           IF(ATHIV(1,4,IMAP,IV)) THEN 
             IAT = MAXIAT + IMAP 
             VAT(1,IAT) = VAT(1,IAT) + EFAT(1,IAT,IV) * VMTMIX(IV) 
           END IF ! ATHIV 
          END DO ! IMAP 
        END IF ! ATHAPS 
C 
      END DO ! IV 
C 
      RETURN 
      END 
