      REAL FUNCTION DIRNL_HR(IHOUR,MY,IVGAS) 
C 
C Program to calculate diurnal emission, i.e., evaporative emissions mainly  
C due to fuel/fuel system pressure increases caused by increases in air  
C temperature ( T > 40 F and delta_T > 0 F ) (increases in fuel/fuel system 
C temperature).  
C 
C Called by EVP_HRLY 
C 
C Calls 
C 
C Changes : (Last change first) 
C 
C  31 Jan 01 @EPA-bag moved high altitude adjustment to DIRNLDAY.FOR; no longer 
C            applying the to liquid leaks 
C  01 Nov 00 @EPA-bag BUG296 add IVGAS to DIRNL_HR parameter list; precontrol vehicles 
C            have carbureted pressure fail and liquid leak diurnal fractions only 
C  07 Aug 00 @EPA bug#252&253 soak time intervals determining one day, two day, and 
C            three day diurnals; full day and partial day weighted sums for 1996-98 
C            model years 
C  17 Jul 00 AIR Task 06: Added code to allow for user-supplied Diurnal Soak 
C            Activity Fractions. 
C  07 Aug 00 @EPA-bag bug215 Added DIURNL24B array; they are weighted 1999 portion 
C            of the full day diurnals for the model years 1996-1998; DIURNL24A are 
C            are used for all model years and are also the 1995 portion of the full 
C            day diurnals for the model years 1996-1998 
C  24 Jul 00 @EPA-bag D15 Moved call to multiday diurnals to EVP_HRLY.FOR; multi-day 
C            ratios are calculated once for each MY and IV and then passed through 
C            DIRNLDAT.FOR 
C  28 Jun 00 @EPA-ddj D15 Added code for multi-day diurnals; added IV and MY 
C            to parameter list 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code 
C  03 Mar 99 @DynTel-bag 2-664  Replaced IV with IVGAS in parameter list. 
C            Renamed hourly diurnal variables. 
C  27 Nov 98 @DynTel-bag 2-677  Fixed negative diurnal emissions 
C            calculation for vehicles without emissions controls 
C            (MY < 1972). Calculations are checked with vehicles 
C            with emissions controls which fail the pressure test. 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: IHOUR, MY, IVGAS 
C    modules: 
C    DIRNLDAT - NVHCLES_SOAKING,DIURNL_FRCTN,DIURNL24A, 
C               DIURNL24B,DIU1,DIU2,NRMLIZ_FR,NRMLIZ_ID, 
C               R12, R23, R12Y96A,R12Y96B, R23Y96A, R23Y96B 
C    DITEMDAT - MIN_TEMP_TIME 
C    common blocks: 
C    /TEMPS/  TEMHRLY,TEMMIN,TEMMAX 
C    /EVAPACT/ DSACTUSER,DSACTFLAG 
C 
C 
C  Output on return: 
C 
C    argument list  : 
C    common blocks  : 
C    function value : DIRNL_HR 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name            Type              Description 
C  ------           ----  ------------------------------------------------------- 
C   ALPHA            R     Temporary storage variable for diurnals. It used as 
C                          a check for negative hourly diurnal values. 
C   IHOUR            I     MOBIEL6 standard hourly index (IHOUR=1 mapped to the hour 6a.m.-7a.m.) 
C   IHR              I     Diurnal hourly index (IHR=1 mapped to the hour 12a.m.-1a.m.) 
C   IHR_SOAK         I     Soak time hourly do loop index 
C   ZERO_LMT         R     Reasonable estimate of zero 
C   IFUEL_SYS        I      1    carbureted failing pressure 
C                           2    carbureted failing only purge 
C                           3    carbureted passing both purge and pressure 
C                           4    FI failing pressure 
C                           5    FI failing only purge 
C                           6    FI passing both purge and pressure 
C                           7    gross liquid leakers 
C   C                I     Soak time correction factor equals 48 or 72 
C   PINMY96N98       L     Phase-in model years 1996-98 
C   R12              R     Ratio of the 2nd day's diurnal 
C                          emissions to the 1st day's 
C   R23              R     Ratio of the 3rd day's diurnal 
C                          emissions to the 2nd day's 
C   R12Y96A          R     First ratio of the 2nd day's diurnal 
C                          emissions to the 1st day's for phase-in 
C                          model years 1996-98 
C   R12Y96B          R     Second ratio of the 2nd day's diurnal 
C                          emissions to the 1st day's for phase-in 
C                          model years 1996-98 
C   R23Y96A          R     First ratio of the 3rd day's diurnal 
C                          emissions to the 2nd day's for phase-in 
C                          model years 1996-98 
C   R23Y96B          R     Second ratio of the 3rd day's diurnal 
C                          emissions to the 2nd day's for phase-in 
C                          model years 1996-98 
C   E2DAYD           R     Second day diurnal emission 
C   E3DAYD           R     Third day diurnal emission 
C   PRECNTRL         L     logical variable which checks for pre-control vehicles 
C   FAIL_PRSRE       I     constant(=1) used in indexing diurnal fractions 
C                          for fail pressure 
C   DFRCN            R     local variable for holding and applying diurnal fractions; needed for 
C                          pre-control vehicles 
C   DFRCN_NRML       R     local variable for holding and applying diurnal fraction normalization; 
C                          needed for pre-control vehicles 
C 
      USE DIRNLDAT, ONLY : NVHCLES_SOAKING,DIURNL_FRCTN,DIURNL24A, 
     *                     DIURNL24B,DIU1,DIU2,NRMLIZ_FD, 
     *                     R12, R23, R12Y96A,R12Y96B, R23Y96A, R23Y96B 
      USE DITEMDAT, ONLY : MIN_TEMP_TIME 
C 
      IMPLICIT NONE 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'TEMPS.I' 
      INCLUDE 'EVAPACT.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IHOUR 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IVGAS 
C 
C     Declare local variables 
C 
      INTEGER :: IHR 
      INTEGER :: IHR_SOAK 
      INTEGER :: IHR_SOAK_FRC 
      INTEGER :: IFUEL_SYS 
      INTEGER, SAVE :: FAIL_PRSRE=1 
C 
      REAL :: ALPHA 
      REAL :: DFRCN 
      REAL :: DFRCN_NRML 
      REAL :: E2DAYD, E3DAYD 
      REAL, SAVE :: ZERO_LMT = 0.0001 
C 
      LOGICAL :: PINMY96N98 
      LOGICAL :: PRECNTRL=.FALSE. 
C 
C     Initializations: 
C 
      ALPHA  = 0.0 
      DIRNL_HR = 0.0 
      IHR = IHOUR + 6  !diurnal fractions are determined and indexed according to 
C                       a 24 hour clock which has its first hour, IHR=1, during 
C                       the interval 12a.m. to 1a.m.; the hourly index for the rest 
C                       of the code begins, IHOUR=1, with the hour 6a.m. to 7a.m. 
C 
C     check for pre-control vehicles: 
C     IVGAS = 1,2,3,4,5 : light duty gas vehicles 
C     IVGAS = 6,7,8,9,10,11,12,13,15 : heavy duty gas vehicles and gas buses 
C     IVGAS = 14 : motorcycles 
      PRECNTRL= (IVGAS<=3 .AND. MY<=1971) .OR. 
     &          ((IVGAS==4 .OR. IVGAS==5) .AND. MY<=1978) .OR. 
     &          ((IVGAS>=6 .AND. IVGAS/=14) .AND. MY<=1984) 
C 
C         IFUEL_SYS  :  1    carbureted failing pressure 
C                       2    carbureted failing only purge 
C                       3    carbureted passing both purge and pressure 
C                       4    FI failing pressure 
C                       5    FI failing only purge 
C                       6    FI passing both purge and pressure 
C                       7    gross liquid leakers 
C 
      PINMY96N98 = MY>=1996 .AND. MY<=1998 
C 
C     No Diurnals For These Temperature Limits and Time Intervals 
C 
      IF( (  TEMMAX-TEMMIN <= ZERO_LMT) .OR. 
     &    (         TEMMAX <= 40.0    ) .OR. 
     &    ( TEMHRLY(IHOUR) <= 40.0    ) .OR. 
     &    (            IHR <= 6       )       ) GOTO 99 
C 
C 
C       ACTIVITY DISTRIBUTIONS : 
C       The hourly diurnal emission for a vapor/gaseous leaker is  
C       determined from a sum of the activity weighted and diurnal  
C       temporal fraction weighted 24 hour diurnal emission over a 
C       72 hour soak period : 
C 
C (i)   Diurnal Activity Distribution (hourly diurnal emissions as  
C       a function of daily temperature variations and normalized  
C       to 1.00 over a twenty four hour period) 
C (ii)  Diurnal temporal fractions or interrupted diurnal fractions are  
C       a soak time (up to 72 hours) dependent distribution of percentages 
C       of possible soak lengths and 24 hour diurnal emission throughout 
C       the MOBILE6 18 daily periods (24 hour day) 
C       (The activity weight/fraction, NVHCLES_SOAKING,has been coded 
C        using appendix E of report M6.EVP.002.) 
C 
C       Soak time loop IHR_SOAK = 1 to 72 hours 
        DO IHR_SOAK = 2,72 
C         The hourly diurnal fractions use a different clock (see above) 
          IHR_SOAK_FRC = IHR_SOAK 
C         The hourly diurnal fractions are all the same 
C         for soak hours 20 through 72. 
          IF( IHR_SOAK_FRC > 20 ) IHR_SOAK_FRC = 20 
C 
          DO IFUEL_SYS=1,7 
C 
C           except for liquid leaks pre-control vehicles use 
C           fail pressure diurnal fractions 
            DFRCN=DIURNL_FRCTN(IHR,IHR_SOAK_FRC,IFUEL_SYS) 
            DFRCN_NRML=NRMLIZ_FD(IFUEL_SYS,IHR_SOAK_FRC) 
            IF(PRECNTRL .AND. IFUEL_SYS/=7) THEN 
              DFRCN=DIURNL_FRCTN(IHR,IHR_SOAK_FRC,FAIL_PRSRE) 
              DFRCN_NRML=NRMLIZ_FD(IFUEL_SYS,FAIL_PRSRE) 
            ENDIF 
C 
            IF (IHR_SOAK >= IHR-MIN_TEMP_TIME+1) THEN   ! Multi-day diurnals 
C 
C             3rd day diurnal 
              IF (IHR_SOAK >= IHR-MIN_TEMP_TIME+46) THEN 
                IF (PINMY96N98) THEN   ! Phase-in MY 1996-98 
                  E3DAYD = ( R12Y96A(IFUEL_SYS)*R23Y96A(IFUEL_SYS)* 
     &                       DIURNL24A(IFUEL_SYS) ) + 
     &                     ( R12Y96B(IFUEL_SYS)*R23Y96B(IFUEL_SYS)* 
     &                       DIURNL24B(IFUEL_SYS) ) 
                ELSE 
                  E3DAYD = R12(IFUEL_SYS)*R23(IFUEL_SYS)* 
     &                     DIURNL24A(IFUEL_SYS) 
                END IF 
                ALPHA = DFRCN*E3DAYD 
C 
C             2nd day diurnal 
              ELSEIF(IHR_SOAK >= IHR-MIN_TEMP_TIME+22) THEN 
                IF (PINMY96N98) THEN   ! Phase-in MY 1996-98 
                  E2DAYD = R12Y96A(IFUEL_SYS)*DIURNL24A(IFUEL_SYS) 
     &                   + R12Y96B(IFUEL_SYS)*DIURNL24B(IFUEL_SYS) 
                ELSE 
                  E2DAYD = R12(IFUEL_SYS)*DIURNL24A(IFUEL_SYS) 
                END IF 
                ALPHA = DFRCN*E2DAYD 
C 
C             1rst day diurnal 
              ELSE 
                IF (PINMY96N98) THEN   ! Phase-in MY 1996-98 
                  ALPHA = DFRCN* 
     &                    (DIURNL24A(IFUEL_SYS)+DIURNL24B(IFUEL_SYS)) 
                ELSE 
                  ALPHA = DFRCN*DIURNL24A(IFUEL_SYS) 
                ENDIF 
C 
              ENDIF !end of multi-day conditional 
C 
C             normalization of multi-day diurnal fractions 
              IF (DFRCN_NRML > 0.) ALPHA=ALPHA/DFRCN_NRML 
C 
C           partial day diurnals 
            ELSE 
              IF (PINMY96N98) THEN   ! Phase-in MY 1996-98 
                ALPHA = DFRCN* 
     &                  (DIURNL24A(IFUEL_SYS)+DIURNL24B(IFUEL_SYS)) 
              ELSE 
                ALPHA = DFRCN*DIURNL24A(IFUEL_SYS) 
            ENDIF 
C 
            ENDIF ! finish discriminating between diurnal emission times 
            ALPHA = AMAX1(ALPHA, 0.) 
C 
C       Apply default or user-supplied diurnal soak fractions 
            IF(DSACTFLAG.EQ.0) THEN 
              DIRNL_HR=DIRNL_HR + ALPHA*NVHCLES_SOAKING(IHOUR,IHR_SOAK) ! soak hour index is 1 at the interval 6a.m.-7a.m. 
            ELSE 
              DIRNL_HR=DIRNL_HR + ALPHA*DSACTUSER(IHOUR,IHR_SOAK) ! soak hour index is 1 at the interval 6a.m.-7a.m. 
            ENDIF 
C 
          END DO ! end fuel system loop 
C 
        END DO ! end soak time loop 
C 
 99   RETURN 
      END 
