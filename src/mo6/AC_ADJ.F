      REAL FUNCTION AC_ADJ(IV,IH,ACTUAL_MY,ICY,DUM_ACCF) 
C 
C Calculate the fraction of vehicles with functioning  
C air conditioners. (Assume the repair rates for cars and trucks  
C are the same.) Estimate the fraction using freon as the coolant 
C based on model year. 
C 
C Changes : 
C 
C  04 Dec 00 @EPA-bag Task X6, added ACTUAL_MY to parameter list and substituted 
C            for MY 
C  04 Aug 00 @EPA-bg  bug203 removed absolute humidity from 
C            parameter list of HEAT; used common, ALUHIN.I, 
C            in RELHUM.FOR to pass ABSHUM 
C 10 Feb 00 @DynTel-DDJ, 1-013 New Function. 
C 
C 
C Called by:  HRLOOP 
C 
C Calls:  None. 
C 
C Input on call: 
C 
C   argument list : IV, MY, ICY, DUM_ACCF 
C 
C   common blocks : 
C 
C Output on return: 
C 
C   argument list : AC_ADJ 
C 
C   common blocks : 
C 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C IH               I        Hour of the day based on the 1 = 6:00 AM as 
C                           defined by the EPA for MOBILE6 
C IP               I        Pollutant ( IP = 1(HC), 2(CO), 3(NOx) ) 
C VHCLS_FCNNG_AC   R        Fraction of vehicles with functioning AC 
C CAR_OR_TRUCK     I        Car/truck indicator: 1 for car 
C                                                2 for truck 
C IAC              I        Coolant type: 1 freon 
C                                         2 R-124a 
C IMY              I        Model year index 
C IV_AGE           I        Vehicle age 
C 
C Notes : 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'AC.I' 
      INCLUDE 'TEMPS.I' 
C 
      INTEGER, INTENT(IN)     :: IV 
      INTEGER, INTENT(IN)     :: IH 
      INTEGER, INTENT(IN)     :: ACTUAL_MY 
      INTEGER, INTENT(IN)     :: ICY 
C 
      REAL, INTENT(IN)        :: DUM_ACCF 
C 
      REAL T_FAHR 
      REAL DMND_FCTR_0 
      REAL DMND_FCTR_1 
      REAL DMND_FCTR_MRNG 
      REAL DMND_FCTR_DA 
      REAL DMND_FCTR_PK_SN 
      REAL DMND_FCTR_NT 
      REAL DMND_FCTR_ALL 
      REAL DMND_FCTR 
      REAL VHCLS_FCNNG_AC 
C 
      INTEGER HEAT 
      INTEGER HEAT_INDX 
      INTEGER INERR 
      INTEGER CAR_OR_TRUCK 
      INTEGER IAC 
      INTEGER IMY 
      INTEGER IV_AGE 
C 
      AC_ADJ = 0.0 
C 
      IF(DUM_ACCF .GT. 0.0) THEN 
        T_FAHR = TEMHRLY(IH) 
        HEAT_INDX =  HEAT(T_FAHR) 
C 
        DMND_FCTR_0 = 0. 
        DMND_FCTR_1 = 1. 
        DMND_FCTR_MRNG =  DMND(1,1) + HEAT_INDX*(DMND(2,1) + 
     &                    HEAT_INDX*DMND(3,1) ) 
        DMND_FCTR_PK_SN = DMND(1,2) + HEAT_INDX*(DMND(2,2) + 
     &                    HEAT_INDX*DMND(3,2) ) 
        DMND_FCTR_DA =    DMND(1,3) + HEAT_INDX*(DMND(2,3) + 
     &                    HEAT_INDX*DMND(3,3) ) 
        DMND_FCTR_NT =    DMND(1,4) + HEAT_INDX*(DMND(2,4) + 
     &                    HEAT_INDX*DMND(3,4) ) 
        DMND_FCTR_ALL =   DMND(1,5) + HEAT_INDX*(DMND(2,5) + 
     &                    HEAT_INDX*DMND(3,5) ) 
C 
C MORNING/AFTERNOON 
C Time of day is between SUNRISE and begining of PEAKSUN 
C or between end of PEAKSUN and SUNSET 
C 
        IF( (SUNSET .GT. SUNRISE .AND. 
     &      (IH .GE. SUNRISE .AND. IH .LT. PEAKSUN(1) ) ) .OR. 
     &      (SUNSET. LT. SUNRISE .AND. 
     &      (IH .GE. SUNRISE .OR. IH .LT. PEAKSUN(1) )  ) .OR. 
     &      (IH .GE. PEAKSUN(2) .AND. IH .LT. SUNSET ) )THEN 
C 
          SELECT CASE(HEAT_INDX) 
            CASE(:65) 
              DMND_FCTR = DMND_FCTR_0 
            CASE(66:74) 
              DMND_FCTR = DMND_FCTR_DA 
            CASE(75:96) 
              DMND_FCTR = DMND_FCTR_MRNG 
            CASE(97:110) 
              DMND_FCTR = DMND_FCTR_ALL 
            CASE(111:) 
              DMND_FCTR = DMND_FCTR_1 
          END SELECT 
C 
C PEAK SUN 
C Time of day is during PEAKSUN 
C 
        ELSE IF((IH .GE. PEAKSUN(1)) .AND. (IH .LT. PEAKSUN(2))) THEN 
C 
          SELECT CASE(HEAT_INDX) 
            CASE(:65) 
              DMND_FCTR = DMND_FCTR_0 
            CASE(66:76) 
              DMND_FCTR = DMND_FCTR_DA 
            CASE(77:104) 
              DMND_FCTR = DMND_FCTR_PK_SN 
            CASE(105:110) 
              DMND_FCTR = DMND_FCTR_ALL 
            CASE(111:) 
              DMND_FCTR = DMND_FCTR_1 
          END SELECT 
C 
C NIGHT 
C Time of day is between SUNSET and SUNRISE 
C 
        ELSE IF ( (SUNSET .GT. SUNRISE .AND. 
     &          (IH .GE. SUNSET .OR. IH .LT. SUNRISE ) ) .OR. 
     &          (SUNSET .LT. SUNRISE .AND. 
     &          (IH .GE. SUNSET .AND. IH .LT. SUNRISE) )) THEN 
C 
          SELECT CASE(HEAT_INDX) 
            CASE(:74) 
              DMND_FCTR = DMND_FCTR_0 
            CASE(75:101) 
              DMND_FCTR = DMND_FCTR_NT 
            CASE(102:110) 
              DMND_FCTR = DMND_FCTR_ALL 
            CASE(111:) 
              DMND_FCTR = DMND_FCTR_1 
          END SELECT 
C 
        ENDIF 
C 
C      Interpolate between Nighttime (NT) and the appropriate Daytime 
C      Demand Factor based on Cloud Cover. 
C 
        DMND_FCTR = DMND_FCTR*(1 - CLD_CVR) + DMND_FCTR_NT*CLD_CVR 
C 
C If the demand factor is zero no A/C is being used so exit. 
C 
        IF (DMND_FCTR .GT. 0.00001) THEN   ! ELSE ACCF = 0.0 
C 
          IV_AGE = ICY - ACTUAL_MY + 1 
          IV_AGE = MIN(IV_AGE,20) 
C 
          IF ( IV .EQ. 1 ) THEN 
            CAR_OR_TRUCK = 1 
          ELSE IF ( IV .GE. 2 .AND. IV .LE. 5 ) THEN 
              CAR_OR_TRUCK = 2 
            ELSE IF ( IV .GT. 5 ) THEN 
              CALL QUITER(0.0,IV,172,INERR) 
          END IF 
C 
          IF(ACTUAL_MY .LT. 1972) THEN 
            IMY = 1 
          ELSE IF(ACTUAL_MY .GT. 1998) THEN 
            IMY = 27 
          ELSE 
            IMY = ACTUAL_MY - 1971 
          END IF 
C 
C Estimate the fraction using R-134a as the coolant: 
C 
          IF(ACTUAL_MY .LE. 1995) THEN 
            IAC = 1      ! FREON 
          ELSE 
            IAC = 2      ! R-134a 
          END IF 
C 
          VHCLS_FCNNG_AC = VHCLES_AC_EQPD(CAR_OR_TRUCK,IMY) * 
     &            ( 1 - NTFCNNG_AC_UNITS(IAC,IV_AGE) ) 
C 
C Finally, calculate the net air conditioning pollution emission factor 
C 
          AC_ADJ = DUM_ACCF * DMND_FCTR * VHCLS_FCNNG_AC 
          AC_ADJ = MAX(AC_ADJ, 0.0) 
C 
        END IF     ! Demand Factor Greater than 0.00001 LOOP. 
C 
      END IF       ! DUM_ACCF GT 0.0 
C 
      RETURN 
      END 
