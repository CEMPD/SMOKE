      SUBROUTINE NGV_EF(IP,IV,MY,IDX,VMTAGE) 
C 
C    Computes the base emission factors for NGV fueled vehicles 
C    and Corrects for HC speciation. NGV includes CNG and LNG vehicles. 
C 
C   Called by EFCALX and HCCALX 
C 
C   CALLs:  HDCPIC 
C 
C   Changes: 
C 
C   07 Jul 00 AIR Task 05: Removed unused ICY from parameter list. 
C   24 Feb 00 @EPA-elg    New Subroutine Created 
C 
C  Input on call: 
C    argument list: IP,IV,MY,IDX,VMTAGE 
C    common blocks: 
C 
C    /FLAGS4/  NMHFLG 
C    /NGVDAT1/ NGVZML, NGVDET, NGVHIGH, 
C              NGVFLG1, NGVFLG2, NGV_START, 
C              NATGAS 
C    /NGVDAT2/ NGV_NZML, NGV_NDET, NGV_HIGH, 
C              NGV_HDZML, NGV_HDDET, NGV_P 
C 
C  Output on return: 
C    common blocks: 
C 
C  Local variable dictionary: 
C 
C   Name    Type              Description 
C  ------   ----  --------------------------------------------------- 
C  CONV      R    GBHP conversion to g/mi units for HDVs. 
C  NATGAS    R    Final NGV Emission Factors passed to HRLOOP. 
C  NGVDET    R    Deterioration rate (10,000K miles) of NGV vehicles. 
C  NGVFLG1   I    Toggle for Alternate NGV Phase-in Fractions 0-DEFAULT 1-Alternate. 
C  NGVFLG2   I    Toggle for Alternate NGV Emission Factors 0-DEFAULT 1-Alternate. 
C  NGVHIGH   R    High Emitter Emission Level for NGVs. 
C  NGV_START I    Earliest NGV Vehicle that can be modeled.  Default = 1994 MY 
C  NGVZML    R    Zero mile Emission Level of NGV vehicles. 
C  NMHC      R    NMHC Conversion factor -  NMHC = 0.93 * NMOG 
C  VOC       R    VOC Conversion factor  -  VOC  = 0.45 * NMOG 
C  THC       R    THC Conversion factor  -  THC  = 11.5 * NMOG 
C 
C   Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'NGV.I' 
C 
C 
C     Declare called functions 
C 
      REAL HDCPIC 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: MY 
      REAL,    INTENT(IN) :: VMTAGE 
C 
C     Declare local variables 
C 
      REAL CONV 
      REAL NMHC 
      REAL FTP 
      REAL THC 
      REAL VOC 
C 
      INTEGER INH,ISR 
C 
      PARAMETER (NMHC = 0.93) 
      PARAMETER (VOC  = 0.45) 
      PARAMETER (THC  = 11.5) 
  
C 
C 
      IF(NGVFLG1.EQ.1 .AND. MY.GE.NGV_START) THEN 
C 
      CONV = HDCPIC(MY,IV) 
      IF(IV.EQ.13 .AND. MY.GT. 1976) CONV = 2.296 
C 
C     Conversion factors for IV=13 and MY>1977 are defined as zero 
C     in HDCFAC array.  For NGV the value of 2.296 will be used. 
C     This is the 1976 value. 
C 
C 
      IF(IV.LE.5) THEN 
       DO ISR = 1,2 
        DO INH = 1,2 
         IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
          IF(INH.EQ.1) THEN 
            FTP = NGV_NZML(IP,IV) + NGV_NDET(IP,IV) * VMTAGE 
          ELSEIF(INH.EQ.2) THEN 
            FTP = NGV_HIGH(IP,IV) 
          ELSE 
            RETURN 
          ENDIF 
          IF(ISR.EQ.1) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = NGV_P(IP,ISR) * FTP 
          ELSEIF(ISR.EQ.2) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = NGV_P(IP,ISR) * FTP 
          ENDIF 
         ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
          IF(INH.EQ.1) THEN 
            NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                               NGVDET(IV,IP,ISR) * VMTAGE 
          ELSEIF(INH.EQ.2) THEN 
            NATGAS(IV,IP,IDX,ISR,INH) = NGVHIGH(IV,IP,ISR) 
          ENDIF 
         ENDIF 
C 
        END DO         !INH LOOP 
       END DO          !ISR LOOP 
C 
C      !IV <= 5 CONDITION BLOCK END 
C 
C 
C 
      ELSEIF (IV.GE.6 .AND. IV.LE.13) THEN 
       INH = 1         !NORMAL ONLY - NO HIGHS 
       ISR = 1         !RUNNING ONLY 
       IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
         IF(IV.EQ.6 .OR. IV.EQ.7) THEN 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,1) + 
     *                              NGV_HDDET(IP,1)*VMTAGE   !LIGHT HEAVY-DUTY 
         ELSEIF(IV.EQ.8  .OR. IV.EQ.9  .OR. IV.EQ.10 .OR. 
     *          IV.EQ.11) THEN 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,2) + 
     *                              NGV_HDDET(IP,2)*VMTAGE   !MEDIUM HEAVY-DUTY 
         ELSEIF(IV.EQ.12  .OR. IV.EQ.13) THEN 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,3) + 
     *                              NGV_HDDET(IP,3)*VMTAGE   !HEAVY HEAVY-DUTY 
         ENDIF 
       ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                            NGVDET(IV,IP,ISR) * VMTAGE 
       ENDIF           !USER DEFAULT EF'S OR USER SPECIFIED 
       NATGAS(IV,IP,IDX,ISR,INH) = NATGAS(IV,IP,IDX,ISR,INH) * CONV 
C 
C      IV >= 6 AND IV <= 13 CONDITION  BLOCK END 
C 
C 
C 
      ELSEIF (IV.EQ.14 .OR. IV.EQ.15) THEN 
       INH = 1         !NORMAL ONLY - NO HIGHS 
       DO ISR = 1,2 
       IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
         IF(IV.EQ.14) THEN            !IV=14 USES IV=1 NGV EF'S 
           FTP = NGV_NZML(IP,1) + NGV_NDET(IP,1) * VMTAGE 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_P(IP,ISR) * FTP 
         ELSEIF(IV.EQ.15) THEN        !IV=15 USES IV=3 NGV EF'S 
           FTP = NGV_NZML(IP,3) + NGV_NDET(IP,3) * VMTAGE 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_P(IP,ISR) * FTP 
         ENDIF 
       ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                            NGVDET(IV,IP,ISR) * VMTAGE 
       ENDIF           !USER DEFAULT EF'S OR USER SPECIFIED 
       END DO 
C 
C      IV = 14 OR IV = 15 CONDITION BLOCK END 
C 
C 
C 
      ELSEIF (IV.GE.16 .AND. IV.LE.23) THEN 
       INH = 1         !NORMAL ONLY - NO HIGHS 
       ISR = 1         !RUNNING ONLY 
       IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
         IF(IV.EQ.16 .OR. IV.EQ.17) THEN 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,1) + 
     *                              NGV_HDDET(IP,1)*VMTAGE   !LIGHT HEAVY-DUTY 
         ELSEIF(IV.EQ.18  .OR. IV.EQ.19  .OR. IV.EQ.20 .OR. 
     *          IV.EQ.21) THEN 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,2) + 
     *                              NGV_HDDET(IP,2)*VMTAGE   !MEDIUM HEAVY-DUTY 
         ELSEIF(IV.EQ.22  .OR. IV.EQ.23) THEN 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,3) + 
     *                              NGV_HDDET(IP,3)*VMTAGE   !HEAVY HEAVY-DUTY 
         ENDIF 
       ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                            NGVDET(IV,IP,ISR) * VMTAGE 
       ENDIF           !USER DEFAULT EF'S OR USER SPECIFIED 
       NATGAS(IV,IP,IDX,ISR,INH) = NATGAS(IV,IP,IDX,ISR,INH) * CONV 
C 
C      IV >= 16 AND IV <= 23 CONDITION BLOCK END 
C 
C 
C 
      ELSEIF (IV.EQ.24) THEN      !MOTORCYCLES USE CAR NGV EF'S 
       INH = 1                    !NORMAL ONLY - NO HIGHS 
       DO ISR = 1,2 
       IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
           FTP = NGV_NZML(IP,1) + NGV_NDET(IP,1) * VMTAGE 
           NATGAS(IV,IP,IDX,ISR,INH) = NGV_P(IP,ISR) * FTP 
       ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                            NGVDET(IV,IP,ISR) * VMTAGE 
       ENDIF           !USER DEFAULT EF'S OR USER SPECIFIED 
       END DO 
C 
C      IV = 24 CONDITION  - MOTORCYCLES  BLOCK END 
C 
C 
      ELSEIF (IV.EQ.25) THEN                 !GASOLINE BUS 
       INH = 1         !NORMAL ONLY - NO HIGHS 
       ISR = 1         !RUNNING ONLY 
       IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,3) + 
     *                              NGV_HDDET(IP,3)*VMTAGE   !HEAVY HEAVY-DUTY 
       ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                            NGVDET(IV,IP,ISR) * VMTAGE 
       ENDIF           !USER DEFAULT EF'S OR USER SPECIFIED 
       NATGAS(IV,IP,IDX,ISR,INH) = NATGAS(IV,IP,IDX,ISR,INH) * CONV 
C 
C 
C      IV = 25 CONDITION  - GAS BUS  BLOCK END 
C 
C 
      ELSEIF (IV.GE.26 .AND. IV.LE.27) THEN  !HEAVY HEAVY-DUTY FOR BUSES 
       INH = 1                               !NORMAL ONLY - NO HIGHS 
       ISR = 1                               !RUNNING ONLY 
       IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGV_HDZML(IP,3) + 
     *                            NGV_HDDET(IP,3)*VMTAGE 
       ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
         NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                            NGVDET(IV,IP,ISR) * VMTAGE 
       ENDIF              !USER DEFAULT EF'S OR USER SPECIFIED 
       NATGAS(IV,IP,IDX,ISR,INH) = NATGAS(IV,IP,IDX,ISR,INH) * CONV 
C 
C 
C 
      ELSEIF (IV.EQ.28) THEN                 !LDDT 3 AND 4 
        INH = 1                              !NORMAL ONLY - NO HIGHS 
        DO ISR = 1,2 
          IF(NGVFLG2.EQ.0 .OR. NGVZML(IV,IP,ISR).LT.0.0) THEN 
            FTP = NGV_NZML(IP,4) + NGV_NDET(IP,4) * VMTAGE 
            NATGAS(IV,IP,IDX,ISR,INH) = NGV_P(IP,ISR) * FTP 
          ELSEIF(NGVZML(IV,IP,ISR).GE.0.0) THEN 
            NATGAS(IV,IP,IDX,ISR,INH) = NGVZML(IV,IP,ISR) + 
     *                             NGVDET(IV,IP,ISR) * VMTAGE 
          ENDIF           !USER DEFAULT EF'S OR USER SPECIFIED 
        END DO 
C 
C      IV = 28  CONDITION BLOCK END 
C 
      ENDIF 
C 
      ENDIF     !NGVFLG1 IF 
C 
C 
C     Correct NGV HC emission factors for the proper speciation as 
C     determined by user input flag NMHFLG. 
C 
      IF(IP.EQ.1) THEN 
        DO ISR = 1,2 
          DO INH = 1,2 
            IF(NMHFLG.EQ.1) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = 
     *                (NATGAS(IV,IP,IDX,ISR,INH) / NMHC) * THC 
            ELSEIF(NMHFLG.EQ.2) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = NATGAS(IV,IP,IDX,ISR,INH) 
            ELSEIF(NMHFLG.EQ.3) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = 
     *                (NATGAS(IV,IP,IDX,ISR,INH) / NMHC) * VOC 
            ELSEIF(NMHFLG.EQ.4) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = 
     *                (NATGAS(IV,IP,IDX,ISR,INH) / NMHC) + 
     *                ((NATGAS(IV,IP,IDX,ISR,INH)/NMHC)*THC) - 
     *                  NATGAS(IV,IP,IDX,ISR,INH) 
            ELSEIF(NMHFLG.EQ.5) THEN 
              NATGAS(IV,IP,IDX,ISR,INH) = 
     *                (NATGAS(IV,IP,IDX,ISR,INH) / NMHC) 
            ENDIF 
          END DO 
        END DO 
      ENDIF 
C 
      RETURN 
      END 
