      SUBROUTINE DIRNLFRC(DIFF_TEMP) 
C 
C This function will calculate the diurnal fraction which 
C will be used in DIRNL_HOUR.FOR. These fractions are the percent 
C of the twenty-four hour diurnal emission (determined in DIRNL24.FOR) 
C produced for each hour, soak time, and vehicle fuel system type. 
C    
C  Called by DIRNLINI 
C 
C  Changes:  (Last change first) 
C 
C  11 Sep 00 @EPA-BAG BUG252&253 Removed partial day diurnal fraction 
C            normalization variable NRMLIZ_ID(IFUEL_SYS,ISOAK) 
C  07 Sep 00 @EPA bug#d4 soak time intervals for full day diurnals 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: DIFF_TEMP 
C    function : 
C    modules  : 
C      DITEMDAT - PTEMPCHANGE,CTEMPCHANGE,TTEMPCHANGE,VPMI,MIN_TEMP_TIME 
C      DIRNLDAT - DIURFRA,DIURFRB,DIURFRC,DIURFRD,DIURFRE,DIURFRF,DIURFRG, 
C                 DIURFRH,DIURFRI,DIURFRJ,DIURFRK,DIURFRL 
C    common blocks: 
C 
C 
C  Output on return: 
C 
C    argument list: 
C    modules  : 
C      DIRNLDAT - DIURNL_FRCTN,NRMLIZ 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name           Type              Description 
C  ------          ----  ------------------------------------------------------- 
C   ISOAK           I     soak time loop index 
C   IHOUR           I     time of day/ hourly loop index 
C   IFUEL_SYS       I     vehicle fuel system type loop index 
C                         1    gross liquid leakers 
C                         2    carbureted passing both purge and pressure 
C                         3    carbureted failing pressure 
C                         4    carbureted failing only purge 
C                         5    FI passing both purge and pressure 
C                         6    FI failing pressure 
C                         7    FI failing only purge 
C  IFULL_DAY       I      through checks with hour loop index and soak 
C                         index determination of whether or not a particular 
C                         IHOUR and ISOAK are a full day diurnal 
C  DIRNL_BGN_TIME  I      time at which the diurnal emissions begin 
C  IYES            I      used to compare flag values; it has the value of 1 
C  INO             I      used to compare flag values; it has the value of 0 
C  ZERO_LIMIT      R      parameter to limit differences which are close to zero 
C 
C 
      USE DITEMDAT, ONLY : PTEMPCHANGE,CTEMPCHANGE,TTEMPCHANGE, 
     *                     VPMI,MIN_TEMP_TIME 
      USE DIRNLDAT, ONLY : DIURNL_FRCTN,DIURFRA,DIURFRB,DIURFRC,DIURFRD, 
     *                     DIURFRE,DIURFRF,DIURFRG,DIURFRH,DIURFRI, 
     *                     DIURFRJ,DIURFRK,DIURFRL,NRMLIZ_FD 
C 
      IMPLICIT NONE 
C 
C     declare parameter list 
C 
      REAL, INTENT(IN) :: DIFF_TEMP 
C 
C     Declare local variables. 
C 
      INTEGER :: ISOAK 
      INTEGER :: IHOUR 
      INTEGER :: IFUEL_SYS 
      INTEGER :: IFULL_DAY 
      INTEGER, PARAMETER :: IYES=1,INO=0 
      REAL, PARAMETER :: ZERO_LIMIT=0.001 
C 
C     INITIALIZATIONS 
C 
      IFULL_DAY = INO 
      DIURNL_FRCTN = 0.0 
      NRMLIZ_FD = 0.0 
C 
      IF(DIFF_TEMP < ZERO_LIMIT ) GO TO 99 ! no diurnals for zero 
C 
C     do loop for IHOUR = 6(5a.m.-6a.m.) to 24(11p.m.-12a.m.) hour 
C     FOR IHOUR 1 to 5 there is no diurnal emission 
      DO IHOUR = 6, 24 
C 
C       Third do loop ISOAK = 1 to 20 hour 
C       The hourly diurnal fraction are all the same 
C       for soak time from 20 hour to 72 hour, so that 
C       we only need to calculate the hourly 
C       diurnal fraction for the first 20 hours. 
        DO ISOAK = 1, 20 
C   
C        *The range for the first full day diurnal 
C            IHOUR-MIN_TEMP_TIME+1 <= ISOAK <= IHOUR-MIN_TEMP_TIME+21 
C 
C        *The range for second full day diurnal 
C            IHOUR-MIN_TEMP_TIME+21 < ISOAK <= IHOUR-MIN_TEMP_TIME+45 
C 
C        *The range for third full day diurnal 
C            IHOUR-MIN_TEMP_TIME+45 < ISOAK 
C 
          IF (ISOAK < IHOUR-MIN_TEMP_TIME + 1 ) THEN 
            IFULL_DAY = INO 
          ELSE 
            IFULL_DAY = IYES 
          ENDIF 
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C 
          IF( TTEMPCHANGE(IHOUR,ISOAK) == -1.0) THEN 
C 
            DIURNL_FRCTN(IHOUR,ISOAK,1) = 0.0 
            DIURNL_FRCTN(IHOUR,ISOAK,2) = 0.0 
            DIURNL_FRCTN(IHOUR,ISOAK,3) = 0.0 
            DIURNL_FRCTN(IHOUR,ISOAK,4) = 0.0 
            DIURNL_FRCTN(IHOUR,ISOAK,5) = 0.0 
            DIURNL_FRCTN(IHOUR,ISOAK,6) = 0.0 
            DIURNL_FRCTN(IHOUR,ISOAK,7) = 0.0 
C 
          ELSE 
C 
C           Calculate the Hourly diurnal fraction 
C     
C         do loop for IFUEL_SYS = 1 to 7 for 7 vehicle fuel system categories 
C         IFUEL_SYS  :  1    carbureted failing pressure 
C                       2    carbureted failing only purge 
C                       3    carbureted passing both purge and pressure 
C                       4    FI failing pressure 
C                       5    FI failing only purge 
C                       6    FI passing both purge and pressure 
C                       7    gross liquid leakers 
C 
            DO IFUEL_SYS = 1, 7 
C 
              DIURNL_FRCTN(IHOUR,ISOAK,IFUEL_SYS) = 
     *           DIURFRA(IFUEL_SYS) 
C 
     *         +(DIURFRB(IFUEL_SYS)+DIURFRC(IFUEL_SYS)*VPMI) 
     *           *TTEMPCHANGE(IHOUR,ISOAK) 
C 
     *         +DIURFRD(IFUEL_SYS)*PTEMPCHANGE(IHOUR) 
C 
     *         +(DIURFRE(IFUEL_SYS)+DIURFRF(IFUEL_SYS)*VPMI) 
     *           *CTEMPCHANGE(IHOUR) 
C 
     *         +(DIURFRG(IFUEL_SYS)+DIURFRH(IFUEL_SYS)*VPMI) 
     *           *PTEMPCHANGE(IHOUR) 
     *           *TTEMPCHANGE(IHOUR,ISOAK) 
C 
     *         +(DIURFRI(IFUEL_SYS)+DIURFRJ(IFUEL_SYS)*VPMI) 
     *           *PTEMPCHANGE(IHOUR) 
     *           *PTEMPCHANGE(IHOUR) 
C 
     *         +(DIURFRK(IFUEL_SYS)+DIURFRL(IFUEL_SYS)*VPMI) 
     *           *CTEMPCHANGE(IHOUR) 
     *           *CTEMPCHANGE(IHOUR) 
C 
              IF(DIURNL_FRCTN(IHOUR,ISOAK,IFUEL_SYS)<0.00001) 
     *           DIURNL_FRCTN(IHOUR,ISOAK,IFUEL_SYS) = 0.0 
C 
              IF(IFULL_DAY == IYES) THEN 
                NRMLIZ_FD(IFUEL_SYS,ISOAK) = NRMLIZ_FD(IFUEL_SYS,ISOAK)+ 
     *                              DIURNL_FRCTN(IHOUR,ISOAK,IFUEL_SYS) 
              ENDIF 
            END DO ! FUEL SYSTEM loop 
C 
          END IF 
C 
          END DO ! SOAK loop 
C 
        END DO ! HOUR loop 
C 
 99   RETURN 
      END SUBROUTINE DIRNLFRC 
 
