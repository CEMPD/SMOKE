      SUBROUTINE DIRNLINI() 
C 
C  DIRNLINI calls the three subroutines which determine the hourly and 
C  soak time dependent diurnal fractions. These subroutines are independent 
C  of vehicle type and vehicle age, so the call is made in HCCALX before 
C  the vehicle and age loops. 
C 
C  Called by HCCALX 
C 
C  Calls DIRNLTEM,DIRNLFRC,ININVEH 
C 
C  Changes:  (Last change first) 
C 
C  27 Mar 01 @epa-bg    BUG321 put a check on diurnal liquid leaks for 
C            temperature ranges equal to zero 
C  24 Jul 00 @epa-bg    1-005 replaced DIRNL_DAY_LQD with DIU1(7) and DIU2(7) 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: 
C    function : 
C    common blocks: 
C    modules      : 
C       DIRNLDAT.FOR - COGLL 
c 
C  Output on return: 
C 
C    argument list: 
C    common blocks: 
C    modules      : 
C       DIRNLDAT.FOR - DIU1(7),DIU2(7) 
C 
C  local variable dictionary: 
C 
C     Name           Type              Description 
C    ------          ----  --------------------------------------------- 
C  DIFF_TEMP          R      difference between high and low temperatures 
C 
C 
      USE DIRNLDAT, ONLY : DIU1,DIU2,COGLL 
      IMPLICIT NONE 
C 
      REAL :: DIFF_TEMP 
C 
      CALL ININVEH() 
      CALL DIRNLTEM(DIFF_TEMP) 
      CALL DIRNLFRC(DIFF_TEMP) 
C 
      DIU1=0.0 
      DIU2=0.0 
      DIU1(7) = COGLL(1)+COGLL(2)*DIFF_TEMP 
      IF( DIU1(7) < 0.0 .OR. DIFF_TEMP<=0.0 ) DIU1(7) = 0.0 
      DIU2(7) = DIU1(7) 
C 
      RETURN 
      END SUBROUTINE DIRNLINI 
