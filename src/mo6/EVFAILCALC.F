      SUBROUTINE EVFAILCALC(IV) 
C 
C  This subroutine computes the fail rates IL, AZFP, and AZGC ONLY. The 
C  routine duplicates the method found in the speadsheet M6IMFIXC.XLS. 
C 
C  Called by EVAPCALC. 
C 
C  Calls None. 
C 
C  Changes: (Last change first) 
C 
C  15 Sep 00 AIR Task 08: New routine. 
C 
C  
C  Input on call: 
C 
C    argument list : IV 
C    common blocks: 
C     /LDGOBD/  WMIL, OBDF 
C 
C  Output on return: 
C 
C    argument list : None 
C    common blocks: 
C     /EVIM/ FAIL_IL, FAIL_AZFP, FAIL_AZGC 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  ---------------------------------------------------------- 
C  A        R    IL Rate Coefficients 
C  B        R    AZ Pressure Rate Coefficients 
C  C        R    AZ Gas Cap Rate Coefficients 
C  I        I    Do Loop Index, 1=MY95-, 2=MY99+, 3=ETP 
C  AGE      R    Vehicle Age with fractional part 
C  IV       I    Vehicle class 
C  IAGE     I    Vehicle Age, integer portion only 
C  JDX      I    MOBILE6 JDX 
C  X        R    Temporary variable to hold failure rate comparison component 
C  Y        R    Temporary variable to hold no I/M rate comparinson component 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVIM.I' 
      INCLUDE 'LDGOBD.I' 
C 
      INTEGER I 
      INTEGER IAGE 
      INTEGER IV 
      INTEGER JDX 
      REAL    AGE 
      REAL    X 
      REAL    Y 
C 
C  Equation coefficients for IL, AZ pressure, and AZ gas cap. 
C 
      REAL, DIMENSION(3), SAVE :: A=(/30.0     , 1.4462  , 0.8051  /) 
      REAL, DIMENSION(3), SAVE :: B=(/ 0.0115  ,-0.0012  , 0.0005  /) 
      REAL, DIMENSION(3), SAVE :: C=(/ 0.984878,-0.001948,-0.000499/) 
C 
C  Loop for all jdx 
C 
      DO JDX=1,26 
C 
C  Assign age index and age 
C 
      IAGE=JDX-1 
      AGE=IAGE 
C 
C  Compute the IL rates: MY95-, MY99+ and ETP rates 
C  Only LDGV-LDGT2 receive ETP/OBD, HDGV and GasBus do not. 
C  (Note: coefficients result in %, so convert to fractions) 
C 
      FAIL_IL(1,JDX)=(1.0/(1.0/A(1)+A(2)*A(3)**AGE))/100. 
C 
      FAIL_IL(3,JDX)=(1.0/(1.0/A(1)+A(2)*A(3)**(AGE/2.0)))/100. 
C 
      FAIL_IL(2,JDX)=FAIL_IL(3,JDX) 
C 
      IF(IV.LE.5) THEN 
C 
        IF(IAGE.LE.3) 
     *   FAIL_IL(2,JDX)=FAIL_IL(3,JDX)*(1.0-OBDF(1)*WMIL) 
C        
        IF(IAGE.GE.4.AND.IAGE.LE.6) 
     *   FAIL_IL(2,JDX)=FAIL_IL(3,JDX)*(1.0-OBDF(2)*WMIL) 
     *                 -FAIL_IL(3,4)*WMIL*(OBDF(1)-OBDF(2)) 
C 
        IF(IAGE.GE.7) 
     *   FAIL_IL(2,JDX)=FAIL_IL(3,JDX) 
     *                 -WMIL*(FAIL_IL(3,7)*OBDF(2) 
     *                 +FAIL_IL(3,4)*(OBDF(1)-OBDF(2))) 
C 
       ENDIF 
C 
C  Compute the AZ pressure rates, MY95-, MY99+ and ETP rates 
C  Only LDGV-LDGT2 receive ETP and OBD, HDGV and GasBus do not. 
C 
      IF(IAGE.LE.1) THEN 
        FAIL_AZFP(1,JDX)=B(1)+B(2)+B(3) 
        FAIL_AZFP(3,JDX)=B(1)+B(2)*(1./2.)+B(3)*(1./2.)**2 
      ELSE 
        FAIL_AZFP(1,JDX)=B(1)+B(2)*AGE+B(3)*AGE**2 
        FAIL_AZFP(3,JDX)=B(1)+B(2)*(AGE/2.)+B(3)*(AGE/2.)**2 
      ENDIF 
C 
      FAIL_AZFP(2,JDX)=FAIL_AZFP(3,JDX) 
C 
      IF(IV.LE.5) THEN 
C 
        IF(IAGE.LE.3) 
     *   FAIL_AZFP(2,JDX)=FAIL_AZFP(3,JDX)*(1.0-OBDF(1)*WMIL) 
C 
        IF(IAGE.GE.4.AND.IAGE.LE.6) 
     *   FAIL_AZFP(2,JDX)=FAIL_AZFP(3,JDX)*(1.0-OBDF(2)*WMIL) 
     *                   -FAIL_AZFP(3,4)*WMIL*(OBDF(1)-OBDF(2)) 
C 
        IF(IAGE.GE.7) 
     *   FAIL_AZFP(2,JDX)=FAIL_AZFP(3,JDX) 
     *                   -WMIL*(FAIL_AZFP(3,7)*OBDF(2) 
     *                   +FAIL_AZFP(3,4)*(OBDF(1)-OBDF(2))) 
C 
      ENDIF 
C 
C  Compute the AZ gas cap rates, MY95- and MY99+ 
C 
      IF(IAGE.LE.1) THEN 
       FAIL_AZGC(1,JDX)=C(1)+C(2)+C(3) 
      ELSE 
       FAIL_AZGC(1,JDX)=C(1)+C(2)*AGE+C(3)*AGE**2 
      ENDIF 
C 
       FAIL_AZGC(2,JDX)=FAIL_AZGC(1,JDX) 
C 
C  Compare the failure rates to the no I/M values and adjust if necessary 
C   
      DO I=1,2 
        X = FAIL_IL(I,JDX) * FAIL_AZGC(I,JDX) + FAIL_AZFP(I,JDX) 
        Y = EVNOIM(3,JDX,I) + EVNOIM(4,JDX,I) 
C 
        IF(X.GT.Y) THEN 
          FAIL_IL  (I,JDX)=FAIL_IL  (I,JDX) * Y / X 
          FAIL_AZFP(I,JDX)=FAIL_AZFP(I,JDX) * Y / X 
        END IF 
      END DO 
C 
C  End of JDX loop 
C 
      END DO 
C 
      RETURN 
      END 
