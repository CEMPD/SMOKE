      SUBROUTINE DIRNLTEM(DIFF_TEMP) 
C 
C This subroutine will calculate the temperature profile and temperature 
C differences, and vapor pressure terms which will be used to determine 
C the diurnal fractions and the daily diurnals 
C 
C  Called by DIRNLINI 
C 
C  Calls VPR_PRS 
C 
C  Changes:  (Last change first) 
C 
C 23 Jun 00 @EPA-bag bug67  Renamed NEW_TEMP to DIRNL_TEMP 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: 
C    function     : 
C    modules      : 
C       DITEMDAT.FOR - LOWTEMP,MIN_TEMP_TIME 
C    common blocks: 
C    /TEMPS/ TEMMIN,TEMMAX 
C 
C  Output on return: 
C 
C    argument list: DIFF_TEMP 
C    common blocks: 
C    modules      : 
C       DITEMDAT.FOR - CTEMPCHANGE,TTEMPCHANGE,PTEMPCHANGE,VPLO,VPHI, 
C                      VPMI,VPPR,LESS_THN10_FCTR 
C 
C 
C  local variable dictionary: 
C 
C     Name           Type              Description 
C    ------          ----  --------------------------------------------- 
C  DIFF_TEMP          R      difference between high and low temperatures 
C  ISOAK              I      soak time index 
C  IHOUR              I      hourly index 
C  SOAK               I      soak time used to determine the diurnal begin time 
C  DIRNL_BGN_TIME     I      time at which the diurnal emissions begin 
C  ZERO_LIMIT         R      parameter to limit differences which are close to zero 
C  MIN_TEMP           R      daily minimum temperature 
C  MAX_TEMP           R      daily maximum temperature 
C  RATIOTEMP          R      ratio between the user input temperature range to 
C                            standard temperature range 
C  LFDAY              R      The range for the first full day diurnal 
C  DIRNL_BGN_TEMP     R      temperature at which the diurnal emissions begin 
C  DIURTEMP           R      maximum of daily low temperature and temperature when 
C                            the diurnal emissions begin, i.e.,MAX( LOWTEMP, DIRNL_BGN_TEMP) 
C  THI                R      daily maximum temperature passed to vapor pressure routine 
C  TLO                R      daily minimum temperature passed to vapor pressure routine 
C 
C 
      USE DITEMDAT, ONLY : LOWTEMP,MIN_TEMP_TIME,CTEMPCHANGE, 
     *                     TTEMPCHANGE,PTEMPCHANGE,VPLO,VPHI, 
     *                     VPMI,VPPR,LESS_THN10_FCTR,DIRNL_TEMP, 
     *                     MAX_TEMPST,MIN_TEMPST,TEMP 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'TEMPS.I' 
C 
C     declare parameter list 
C 
      REAL, INTENT(OUT) :: DIFF_TEMP 
C 
C     declare local variables 
C 
      INTEGER ISOAK 
      INTEGER IHOUR 
      INTEGER SOAK 
      INTEGER DIRNL_BGN_TIME 
      REAL, PARAMETER :: ZERO_LIMIT=0.001 
      REAL :: MIN_TEMP 
      REAL :: MAX_TEMP 
      REAL :: RATIOTEMP 
      REAL :: LFDAY 
      REAL :: DIRNL_BGN_TEMP 
      REAL :: DIURTEMP 
      REAL :: THI 
      REAL :: TLO 
C 
C     initializations 
      THI = TEMMAX 
      TLO = TEMMIN 
      CTEMPCHANGE = -1 
      TTEMPCHANGE = -1 
      PTEMPCHANGE = -1 
C 
C     Determine the daily max. - min. temperature difference 
C 
      LESS_THN10_FCTR = 0.0 
      DIFF_TEMP = TEMMAX-TEMMIN 
      IF(DIFF_TEMP < ZERO_LIMIT ) GO TO 99 ! no diurnals for zero 
                                           ! max. - min temperature 
                                           ! differences 
      LESS_THN10_FCTR = 1.0 
      IF( (DIFF_TEMP - 10) < ZERO_LIMIT ) THEN 
        DIFF_TEMP=10.0 
        LESS_THN10_FCTR = (TEMMAX-TEMMIN)/10.0 
        THI = (TEMMIN + TEMMAX)/2.0 + 5.0 
        TLO = (TEMMIN + TEMMAX)/2.0 - 5.0 
      ENDIF 
C 
C     determine the estimated vapor pressure terms 
      CALL VPR_PRS(THI,TLO) 
C 
      MIN_TEMP=TEMMIN 
      MAX_TEMP=TEMMAX 
C 
C     Calculate the ratio between the user input temperature range to 
C     standard temperature range 
      RATIOTEMP=(MAX_TEMP-MIN_TEMP)/(MAX_TEMPST-MIN_TEMPST) 
C 
C     Calculate the new temperature profile 
C 
      DO IHOUR=1, 24 
        DIRNL_TEMP(IHOUR)=MIN_TEMP + RATIOTEMP*(TEMP(IHOUR)-MIN_TEMPST) 
      END DO 
C 
C     *Temperature change during the previous hour for ihour = 1 
      PTEMPCHANGE(1)=DIRNL_TEMP(24)-DIRNL_TEMP(23) 
C     *Temperature change during the current hour for ihour = 1 
      CTEMPCHANGE(1)=DIRNL_TEMP(1)-DIRNL_TEMP(24) 
C 
C     finish filling these temperature difference arrays 
      DO IHOUR=2, 24 
C 
C       *Temperature change during the current hour 
        CTEMPCHANGE(IHOUR)=DIRNL_TEMP(IHOUR)-DIRNL_TEMP(IHOUR-1) 
C 
C       *Temperature change during the previous hour 
        PTEMPCHANGE(IHOUR)=CTEMPCHANGE(IHOUR-1) 
      END DO 
C 
C     FIRST do loop ISOAK = 1 to 20 hour 
C     The hourly diurnal fraction are all the same 
C     for soak time from 20 hour to 72 hour, so that 
C     we only need to calculate the hourly 
C     diurnal fraction for the first 20 hours. 
C 
      DO ISOAK = 1, 20 
C 
C     Second do loop IHOUR = 6(5a.m.-6a.m.) to 24(11p.m.-12a.m.) hour 
C       FOR IHOUR 1 to 5 there is no diurnal emission 
        DO IHOUR = 6, 24 
C 
C        *The range for the first full day diurnal 
C 
          LFDAY=IHOUR-MIN_TEMP_TIME +2 
C 
C         IHOUR-MIN_TEMP_TIME +22 
C 
C        *The range for second full day diurnal 
C 
C         IHOUR-MIN_TEMP_TIME +23 
C         IHOUR-MIN_TEMP_TIME +46 
C 
C         The range for third full day diurnal 
C 
C         IHOUR-MIN_TEMP_TIME +47 
C 
C        *Total temperature change prior to the previous hour 
C         For the interrupted diurnal 
C 
C         full day diurnal 
          IF (ISOAK < LFDAY ) THEN 
             SOAK = IHOUR - ISOAK-1 
C 
C            DIRNL_BGN_TIME = MAX( MIN_TEMP_TIME, SOAK+2) 
C 
             IF(MIN_TEMP_TIME >= (SOAK+2)) THEN 
                 DIRNL_BGN_TIME = MIN_TEMP_TIME ! standard minimum temperature occurs at MIN_TEMP_TIME 
                                                ! DIRNL_BGN_TIME is the diurnal emission start time 
             ELSE 
                 DIRNL_BGN_TIME = SOAK+2 
             END IF 
C 
             DIRNL_BGN_TEMP = DIRNL_TEMP(DIRNL_BGN_TIME) ! diurnal start temperature/ temperature at DIRNL_BGN_TIME 
C 
C            DIURTEMP = MAX( LOWTEMP, DIRNL_BGN_TEMP) 
C 
             IF( LOWTEMP >= DIRNL_BGN_TEMP) THEN  !LOWTEMP = 40, i.e., diurnal low temperature limit 
                DIURTEMP = LOWTEMP 
             ELSE 
                DIURTEMP = DIRNL_BGN_TEMP 
             END IF 
C 
             IF (IHOUR <= DIRNL_BGN_TIME .OR. 
     *          (DIRNL_TEMP(IHOUR-1) <= DIURTEMP.AND. 
     *           DIRNL_TEMP(IHOUR) <= DIURTEMP))THEN 
               TTEMPCHANGE(IHOUR,ISOAK)=-1.0 
             ELSE 
               TTEMPCHANGE(IHOUR,ISOAK)=DIRNL_TEMP(IHOUR-2)- 
     *                                  DIRNL_BGN_TEMP 
C 
               IF(TTEMPCHANGE(IHOUR,ISOAK) < 0) 
     *                     TTEMPCHANGE(IHOUR,ISOAK)=0.0 
             END IF 
C 
          END IF 
C  
          IF (ISOAK >= LFDAY) THEN 
            DIRNL_BGN_TIME=MIN_TEMP_TIME 
            IF(IHOUR <= DIRNL_BGN_TIME .OR. 
     *          (DIRNL_TEMP(IHOUR-1) <= DIRNL_BGN_TEMP.AND. 
     *           DIRNL_TEMP(IHOUR) <= DIRNL_BGN_TEMP)) THEN 
              TTEMPCHANGE(IHOUR,ISOAK)=-1.0 
            ELSE 
              TTEMPCHANGE(IHOUR,ISOAK)=DIRNL_TEMP(IHOUR-2)-MIN_TEMP !determine the total temperature change 
            ENDIF 
          END IF 
C 
C      WRITE(30, 50) IV,IHOUR,ISOAK,DIURNL_FRCTN(IV,IHOUR,ISOAK) 
C 50   FORMAT(3I6, F10.6) 
C 
        END DO !hourly loop 
C 
C 
      END DO !soak time loop 
C 
C     CLOSE(30) 
C 
 99   RETURN 
      END SUBROUTINE DIRNLTEM 
C 
C 
C 
C 
      SUBROUTINE VPR_PRS(THI,TLO) 
C 
C  Calculates the vapor pressure terms for the diurnal calculations 
C 
C  Called by DIRNLTEM 
C 
C  Calls 
C 
C  Changes:  (Last change first) 
C 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: THI,TLO 
C    function     : 
C    modules      : 
C       DITEMDAT.FOR - VP 
C    common blocks: 
C    /CITRV1/ RVPICY 
C 
C  Output on return: 
C 
C    argument list: 
C    modules      : 
C       DITEMDAT.FOR - VPLO,VPHI,VPMI,VPPR 
C    common blocks: 
C 
C 
C  local variable dictionary: 
C 
C     Name           Type              Description 
C    ------          ----  --------------------------------------------- 
C    THI              R      daily maximum temperature 
C    TLO              R      daily minimum temperature 
C 
 
      USE DITEMDAT, ONLY : VP,VPLO,VPHI,VPMI,VPPR 
      IMPLICIT NONE 
C 
      INCLUDE 'CITRV1.I' 
      REAL, INTENT(IN) :: THI 
      REAL, INTENT(IN) :: TLO 
C 
      VPLO = VP(1)*RVPICY*exp(VP(2)*RVPICY+(VP(3)*RVPICY-VP(4)) 
     *       /(VP(5)+TLO)) 
      VPHI = VP(1)*RVPICY*exp(VP(2)*RVPICY+(VP(3)*RVPICY-VP(4)) 
     *       /(VP(5)+THI)) 
      VPMI = (VPLO+VPHI)/2 
      VPPR = (VPHI-VPLO)*VPMI 
C 
      RETURN 
      END SUBROUTINE VPR_PRS 
