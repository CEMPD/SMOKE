      FUNCTION METH81(IV,MY,ISR,VMT,IFDS) 
C 
C  METH81 computes the START and RUNNING methane emissions for 1981 
C  and newer light duty cars and trucks. 
C 
C  Called by BIGMET 
C 
C  Calls   : NONE 
C 
C  Changes : (Last change first) 
C 
C  14 Aug 00 @EPA-djb Bug 214, Start methane coeffs are in grams already. 
C  26 Jul 00 @EPA-djb Corrected fuel delivery mapping, converted engine 
C            start methane to grams and added comments. 
C  07 Jul 00 AIR Task 05: Removed unused MYR from parameter list. 
C   5 Apr 00 @EPA-eg      methane coefficients updated 
C  15 Mar 00 @EPA-eg            hrloop revisions 
C  28 Sep 99 @DynTel-ddj for EPA - New FUNCTION 
C      
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: IV, MY, ISR, VMT, IFDS 
C 
C  Output on return: 
C 
C      METH81 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ------------------------------------------------ 
C  CH4ZML_CAR    R    Methane zero mile coefficient for passenger cars (running). 
C  CH4S1_CAR     R    Methane initial slope coefficient for passenger cars (running). 
C  CH4S2_CAR     R    Methane final slope coefficient for passenger cars (running). 
C  CH4C1_CAR     R    Methane initial flex point mileage for passenger cars (running). 
C  CH4C2_CAR     R    Methane final flex point mileage for passenger cars (running). 
C  CH4ZMLS_CAR   R    Methane zero mile coefficient for passenger cars (start). 
C  CH4S1S_CAR    R    Methane slope coefficient for passenger cars (start). 
C  CH4ZML_TR     R    Methane zero mile coefficient for light trucks (running). 
C  CH4S1_TR      R    Methane initial slope coefficient for light trucks (running). 
C  CH4S2_TR      R    Methane final slope coefficient for light trucks (running). 
C  CH4C1_TR      R    Methane initial flex point mileage for light trucks (running). 
C  CH4C2_TR      R    Methane final flex point mileage for light trucks (running). 
C  CH4ZMLS_TR    R    Methane zero mile coefficient for light trucks (start). 
C  CH4S1S_TR     R    Methane slope coefficient for light trucks (start). 
C  J             I    Appropriate index to the coefficient array. 
C  IPC           I    Index to coefficient array for carbureted vehicles. 
C  IPFI          I    Index to coefficient array for ported fuel injection vehicles. 
C  ITBI          I    Index to coefficient array for throttle body fuel injection vehicles. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'INJECT.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
      REAL,    INTENT(IN) :: VMT 
C 
C     Declare local variables 
C 
      INTEGER IPC 
      INTEGER IPFI 
      INTEGER ITBI 
      INTEGER J 
C 
      REAL, DIMENSION(7) :: CH4ZML_CAR 
      REAL, DIMENSION(7) :: CH4S1_CAR 
      REAL, DIMENSION(7) :: CH4S2_CAR 
      REAL, DIMENSION(7) :: CH4C1_CAR 
      REAL, DIMENSION(7) :: CH4C2_CAR 
      REAL, DIMENSION(7) :: CH4ZMLS_CAR 
      REAL, DIMENSION(7) :: CH4S1S_CAR 
C 
      REAL, DIMENSION(5) :: CH4ZML_TR 
      REAL, DIMENSION(5) :: CH4S1_TR 
      REAL, DIMENSION(5) :: CH4S2_TR 
      REAL, DIMENSION(5) :: CH4C1_TR 
      REAL, DIMENSION(5) :: CH4C2_TR 
      REAL, DIMENSION(5) :: CH4ZMLS_TR 
      REAL, DIMENSION(5) :: CH4S1S_TR 
C 
C     Declare external functions 
C 
      REAL METH81 
C 
C 
C     Methane Coefficients for Cars (light duty gasoline vehicles) 
C 
C                   81-82  81-82  83-85  86-93  83-87  88-93  88-93 
C                   CARB    FI    CARB   CARB    FI     TBI    PFI 
C 
      CH4ZML_CAR=(/0.0845,0.0271,0.0721,0.0405,0.0365,0.0240,0.0167/) 
      CH4S1_CAR=(/ 0.019, 0.020, 0.002, 0.019, 0.000, 0.004, 0.004/) 
      CH4S2_CAR=(/ 0.009, 0.020, 0.000, 0.009, 0.006, 0.004, 0.005/) 
      CH4C1_CAR=(/ 2.211, 1.392, 2.000, 1.519, 2.121, 3.218, 1.547/) 
      CH4C2_CAR=(/10.000,26.540,10.000, 7.191, 8.129,10.000, 6.789/) 
C 
      CH4ZMLS_CAR=(/0.2003,0.1005,0.1373,0.1023,0.0963,0.0711,0.0936/) 
      CH4S1S_CAR =(/0.021,0.014, 0.018, 0.004, 0.001, 0.000, 0.002/) 
C 
C 
C     Methane Coefficients for Light Duty Gas Trucks 
C 
C                   81-83   81-87   84-93   88-93   88-93 
C                   CARB     FI     CARB     TBI     PFI 
C 
      CH4ZML_TR=(/ 0.1033, 0.0594, 0.1118, 0.0253, 0.0291/) 
      CH4S1_TR=(/  0.003,  0.000,  0.000,  0.004,  0.004/) 
      CH4S2_TR=(/  0.003,  0.004,  0.003,  0.005,  0.009/) 
      CH4C1_TR=(/  1.235,  4.867,  9.072,  1.625,  1.918/) 
      CH4C2_TR=(/  8.035, 10.000, 10.000,  5.446, 10.000/) 
C 
      CH4ZMLS_TR=(/0.4062, 0.0927, 0.2707,  0.1402, 0.0505/) 
      CH4S1S_TR =(/0.010,  0.003,  0.009,   0.000,  0.000/) 
C 
C  Initialize 
C 
      METH81 = 0.0 
C 
C  Only 1981-1993 model years 
C 
      IF (MY.GE.1981 .AND. MY.LE.1993) THEN 
C 
C  Light Duty Gasoline Vehicles (Cars) 
C 
        IF(IV.EQ.1) THEN 
          IF(MY.LE.1982) THEN 
            IPC   = 1 
            IPFI  = 2 
            ITBI  = 2 
          ELSE IF(MY.LE.1985) THEN 
            IPC   = 3 
            IPFI  = 5 
            ITBI  = 5 
          ELSE IF(MY.LE.1987) THEN 
            IPC   = 4 
            IPFI  = 5 
            ITBI  = 5 
          ELSE 
            IPC   = 4 
            ITBI  = 6 
            IPFI  = 7 
          ENDIF 
C 
          IF(IFDS.EQ.1) THEN 
             J = IPC 
          ELSE IF(IFDS.EQ.2) THEN 
             J = ITBI 
          ELSE 
             J = IPFI 
          END IF 
C 
          IF(ISR.EQ.1) THEN                    !Running emissions 
              IF(VMT.LE.CH4C1_CAR(J)) THEN         !less than CH4C1_CAR 
                METH81 = CH4ZML_CAR(J) 
              ELSE IF(VMT.LE.CH4C2_CAR(J)) THEN    !between CH4C1_CAR and CH4C2_CAR 
                METH81 = CH4ZML_CAR(J) + 
     &                    CH4S1_CAR(J)*(VMT-CH4C1_CAR(J)) 
              ELSE                                 !greater than CH4C2_CAR 
                METH81 = CH4ZML_CAR(J) + 
     &                    CH4S1_CAR(J)*(CH4C2_CAR(J) - CH4C1_CAR(J)) + 
     &                    CH4S2_CAR(J)*(VMT-CH4C2_CAR(J)) 
              ENDIF 
          ELSE                                  !Start emissions 
              METH81 = CH4ZMLS_CAR(J) + CH4S1S_CAR(J)*VMT !grams per engine start 
          ENDIF 
C 
C   Light Duty Gasoline Trucks 
C 
        ELSEIF(IV.GE.2 .AND. IV.LE.5) THEN 
          IF(MY.LE.1983) THEN 
            IPC   = 1 
            IPFI  = 2 
            ITBI  = 2 
          ELSE IF(MY.LE.1987) THEN 
            IPC   = 3 
            IPFI  = 2 
            ITBI  = 2 
          ELSE 
            IPC   = 3 
            ITBI  = 4 
            IPFI  = 5 
          ENDIF 
C 
          IF(IFDS.EQ.1) THEN 
             J = IPC 
          ELSE IF(IFDS.EQ.2) THEN 
             J = ITBI 
          ELSE 
             J = IPFI 
          END IF 
C 
          IF(ISR.EQ.1) THEN                  !Running emissions 
              IF(VMT.LE.CH4C1_TR(J)) THEN         !less than CH4C1_TR 
                METH81 = CH4ZML_TR(J) 
              ELSE IF(VMT.LE.CH4C2_TR(J)) THEN    !between CH4C1_TR and CH4C2_TR 
                METH81 = CH4ZML_TR(J) + CH4S1_TR(J)*(VMT-CH4C1_TR(J)) 
              ELSE                                !greater than CH4C2_TR 
                METH81 = CH4ZML_TR(J) + 
     *                    CH4S1_TR(J)*(CH4C2_TR(J) - CH4C1_TR(J)) + 
     *                    CH4S2_TR(J)*(VMT-CH4C2_TR(J)) 
              END IF 
          ELSE                               !Start emissions 
              METH81 = CH4ZMLS_TR(J) + CH4S1S_TR(J)*VMT   !grams per engine start 
          ENDIF 
C 
        END IF !IV Condition 
      END IF   !MY Condition 
C 
      RETURN 
      END 
