      SUBROUTINE DIURATIO(IV, MY) 
C 
C DIURATIO calculates the ratio of the 2nd day's diurnal 
C emissions to the 1st day's (R12) and the ratio of the 3rd 
C day's diurnal emissions to the 2nd day's (R23). 
C For the following four strata, both R12 and R23 will be 
C set equal to 1: 
C   - the stratum of "gross liquid leakers" GLL 
C   - the stratum of "pre-1972" LDVs and LDTs 
C   - the stratum of "pre-1985" HDTs 
C   - the stratum of motorcycles 
C 
C Changes : 
C 
C  21 Sep 00 @EPA bag bug261 motorcycle ratios set to 1 
C  07 Aug 00 @EPA bug252&253 1999 pass both ratios set to 1 rather than 
C            calculated from equations 
C  07 Aug 00 @EPA-bag bug215 removed weighted 24 hour diurnal array DIURNL24 
C  24 Jul 00 @EPA-bag D15    Moved call to multiday diurnals to EVP_HRLY.FOR; 
C            multi-day ratios are calculated once for each MY and IV and then 
C            passed through DIRNLDAT.FOR 
C  21 Jan 00 @DynTel-kml #1-015, New subroutine. 
C 
C Called by:  EVP_HRLY 
C 
C Calls: 
C 
C Input on call: 
C 
C   argument list : IV, MY 
C 
C   common blocks : 
C    DIRNLDAT - DIU1, DIU2, R12COEFF, R23COEFF 
C    DITEMDAT - LESS_THN10_FCTR 
C    /CITRV1/ RVPICY 
C    /TEMPS/ TEMMIN, TEMMAX 
C 
C Output on return: 
C 
C   argument list : 
C 
C   common blocks : 
C     DIRNLDAT - R12, R23, R12Y96A,R12Y96B, R23Y96A, R23Y96B 
C 
C   function      : 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C MIDTEM           R        Midpoint temperature 
C TR               R        Temporary R12/R23 ratio 
C A, B, C, D, E    R        Coefficients of equations to calculate 
C                           ratios R12 and R23 
C IR               I        Column index for arrays R12COEFF 
C                           and R23COEFF: IR=1,4 
C IV               I        Gasoline fueled vehicle type index 
C MY               I        Vehicle model year 
C IFDS             R        Fuel delivery system IFDS= 
C                              1: carb fail pressure 
C                              2: carb fail purge 
C                              3: carb pass both 
C                              4: fi fail pressure 
C                              5: fi fail purge 
C                              6: fi pass both 
C                              7: liquid leakers 
C DIUXMIN          R        Minimum value of diurnal 
c MAXRATIO         R        Maximum value of R12/R23 ratios 
C 
C Notes : 
C 
      USE DITEMDAT, ONLY : LESS_THN10_FCTR 
      USE DIRNLDAT, ONLY : DIU1, DIU2, R12COEFF, R23COEFF, 
     *                     R12, R23, R12Y96A, R12Y96B, R23Y96A, R23Y96B 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'TEMPS.I' 
C 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
C 
      REAL, SAVE :: DIUXMIN = 0.00001 
      REAL, SAVE :: ZERO_LIMIT=0.001 
C 
      REAL MIDTEM 
      REAL TR 
      REAL A 
      REAL B 
      REAL C 
      REAL D 
      REAL E 
C 
      INTEGER IR 
      INTEGER IFDS 
      INTEGER, SAVE :: IMOTRCYCL=24 
C 
      R12=0.0 
      R23=0.0 
      R12Y96A=0.0 
      R12Y96B=0.0 
      R23Y96A=0.0 
      R23Y96B=0.0 
C 
      IF( LESS_THN10_FCTR < ZERO_LIMIT ) GO TO 99 ! no diurnals for zero 
                                                  ! max. - min temperature 
                                                  ! differences 
      DO IFDS=1,7 
C initialize TR 
        TR=0.0 
C 
C Calculate mapping of IFDS into columns of arrays 
C R12COEFF and R23COEFF: 
C 
C        IFDS= 1,2  Carb fail either purge/pressure 
C              3    Carb pass both purge/pressure 
C              4,5  FI fail either purge/pressure 
C              6    FI pass both purge/pressure 
C              7    Gross leaker 
C 
C   IFDS:                            1 2 3 4 5 6 7 
C   IR=column of R12COEFF/R23COEFF:  1 1 2 3 3 4 none 
C 
 
        SELECT CASE (IFDS) 
           CASE (1)      ! Carb fail either purge/pressure 
              IR = IFDS 
           CASE (2:4) 
              IR = IFDS - 1 
           CASE (5:6) 
              IR = IFDS - 2 
           CASE DEFAULT  ! Gross Liquid Leakers (GLL) 
              IR = 0 
        END SELECT 
C 
C Calculate ratios R12/R23 
C 
        IF( (MY.GE.1972 .OR. VVLDGAS(IV) .NE.1) .AND. !not pre-72 LDGV/T 
     &      (MY.GE.1979 .OR. VVLDGT34(IV).NE.1) .AND. !not pre-79 LDGT3/4 
     &      (MY.GE.1985 .OR. VVHDGV(IV)  .NE.1) .AND. !not pre-85 HDGT 
     &      (MY.LT.1999 .OR. VVLDGAS(IV) .NE.1  .OR. 
     &      IR.GE.3) .AND. !not post-98 LDGV/T passing purge and press. 
     &      (VVDSL(IV) .NE. 1)  .AND.     !Not diesel vehicle/truck 
     &      (IFDS .NE. 7)       .AND.     !Not Gross Liquid Leaker (GLL) 
     &      (IV .NE. IMOTRCYCL) )  THEN   !No MOTORCYCLES 
C 
           A = R12COEFF(IR,1)    ! Calculate R12 
           B = R12COEFF(IR,2) 
           C = R12COEFF(IR,3) 
           D = R12COEFF(IR,4) 
           E = R12COEFF(IR,5) 
C 
C Calculate midpoint temperature 
C 
           MIDTEM = (TEMMIN + TEMMAX) / 2. 
C 
           IF (MY.GE.1996 .AND. MY.LE.1998) THEN 
             IF(IFDS<=3) THEN !NO CARBURETED VEHICLES 
                R12Y96B(IFDS) = 0.0 
                R12Y96A(IFDS) = 0.0 
             ELSE 
               IF (DIU1(IFDS) .GE. DIUXMIN) THEN 
                 TR = E + (A + B*MIDTEM + C*RVPICY + D*MIDTEM*RVPICY) 
     &                 / DIU1(IFDS) 
               ELSE 
                 TR = 0 
               END IF 
               R12Y96A(IFDS) = AMAX1(E,TR) 
               IF (DIU2(IFDS) .GE. DIUXMIN .AND. IFDS/=6) THEN 
                 TR =  E + (A + B*MIDTEM + C*RVPICY + D*MIDTEM*RVPICY) 
     &                  / DIU2(IFDS) 
               ELSE 
                 TR = 0 
               END IF 
               R12Y96B(IFDS) = AMAX1(E,TR) 
               IF(IFDS==6) R12Y96B(IFDS) = 1.0 
             ENDIF 
           ELSE 
             IF (DIU1(IFDS) .GE. DIUXMIN) THEN 
               TR = E + (A + B*MIDTEM + C*RVPICY + D*MIDTEM*RVPICY) 
     &               / DIU1(IFDS) 
             ELSE 
               TR = 0 
             END IF 
             R12(IFDS) = AMAX1(E,TR) 
               IF(IFDS==6.AND.MY>=1996) R12(IFDS) = 1.0 
           END IF 
C 
           A = R23COEFF(IR,1)    ! Calculate R23 
           B = R23COEFF(IR,2) 
           C = R23COEFF(IR,3) 
           D = R23COEFF(IR,4) 
           E = R23COEFF(IR,5) 
           IF (MY.GE.1996 .AND. MY.LE.1998) THEN 
             IF(IFDS<=3) THEN !NO CARBURETED VEHICLES 
               R23Y96A(IFDS) = 0.0 
               R23Y96B(IFDS) = 0.0 
             ELSE 
               IF (DIU1(IFDS) .GE. DIUXMIN) THEN 
                 TR = E + (A + B*MIDTEM + C*RVPICY + D*MIDTEM*RVPICY) 
     &                   / DIU1(IFDS) 
               ELSE 
                TR = 0 
               END IF 
               R23Y96A(IFDS) = AMAX1(E,TR) 
               IF (DIU2(IFDS) .GE. DIUXMIN) THEN 
                 TR = E + (A + B*MIDTEM + C*RVPICY + D*MIDTEM*RVPICY) 
     &               / DIU2(IFDS) 
               ELSE 
                 TR = 0 
               END IF 
               R23Y96B(IFDS) = AMAX1(E,TR) 
               IF(IFDS==6) R23Y96B(IFDS) = 1.0 
             ENDIF 
           ELSE 
             IF (DIU1(IFDS) .GE. DIUXMIN) THEN 
               TR = E + (A + B*MIDTEM + C*RVPICY + D*MIDTEM*RVPICY) 
     &               / DIU1(IFDS) 
             ELSE 
               TR = 0 
             END IF 
               R23(IFDS) = AMAX1(E,TR) 
             IF(IFDS==6.AND.MY>=1996) R23(IFDS) = 1.0 
           END IF 
        ELSE 
           R12(IFDS) = 1. 
           R23(IFDS) = 1. 
           R12Y96A(IFDS) = 1. 
           R12Y96B(IFDS) = 1. 
           R23Y96A(IFDS) = 1. 
           R23Y96B(IFDS) = 1. 
        END IF 
C 
      END DO 
C 
 99   RETURN 
      END 
