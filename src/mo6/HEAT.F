      INTEGER FUNCTION HEAT(TFAHR) 
C 
C  This function calculates the heat index/apparent 
C  temperature HEAT (F) given the temperature TFAHR (F) 
C  and absolute humidity ABSHUM (grain/lb). The heat 
C  index is what the human body tries to compensate 
C  for via perspiration during the summer (equivalent 
C  to winter's wind chill temperature). 
C  Algorithm source: Meisner and Graves, "Apparent 
C  Temperature", Weatherwise, August 1985 
C 
C 
C Called by ACCF. 
C 
C Calls RELHUM, SLTRYM, SLTRYS. 
C 
C Changes: 
C  
C  04 Aug 00 @EPA-bg  bug203 removed absolute humidity from 
C            parameter list of HEAT and RELHUM; used common, 
C            ALUHIN.I, in RELHUM.FOR to pass ABSHUM 
C  21 Jan 00 @DynTel-km  1-013, New Function 
C 
C Input on call: 
C 
C   argument list : TFAHR, ABSHUM 
C 
C Output on return: 
C 
C   function      : HEAT 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C T_CEL            R        Temperature in degrees Celsius (C) 
C E                R        Relative vapor pressure 
C ES               R        Saturation vapor pressure 
C DF               R        Flag - moderate/severe sultriness 
C RH               R        Relative humidity (%) 
C PINF             R        Temporary variable 
C HIDX             R        Heat index in degrees Celsius (C) 
C TCONV            R        Temperature conversion constant 
C FTOC             R        F to C conversion constant 
C CTOF             R        C to F conversion constant 
C 
C Notes : 
C 
C 
      IMPLICIT NONE 
C 
C     Declare parameter list 
C 
      REAL, INTENT(IN)   :: TFAHR 
C 
C     Declare external functions 
C 
      REAL RELHUM 
C 
C     Delcare local variables 
C 
      REAL RH 
      REAL T_CEL 
      REAL ES 
      REAL E 
      REAL PINF 
      REAL DF 
      REAL HIDX 
      REAL TCONV 
      REAL FTOC 
      REAL CTOF 
C 
C  Initialize variables 
C 
      TCONV=32. 
      FTOC=5./9. 
      CTOF=9./5. 
      DF=.0 
      HIDX=.0 
C 
C  Get relative humidity 
C 
      RH=RELHUM(TFAHR) 
C 
C  Compute saturation vapor pressure ES and relative 
C  vapor pressure E from relative humidity RH 
C 
      T_CEL=(TFAHR-TCONV)*FTOC 
      ES=6.11*10.**((7.567*T_CEL)/(239.7+T_CEL)) 
      E=.01*RH*ES 
C 
      PINF=.1*E 
C 
C  Compute "moderate" heat index (C) 
C 
      CALL SLTRYM(T_CEL,DF,PINF,HIDX) 
C 
C  Compute "severe" heat index (C) 
C 
      IF(DF.LT.0.) CALL SLTRYS(T_CEL,PINF,HIDX) 
C 
C  Convert to Fahrenheit (F) 
C 
      HEAT=TCONV+CTOF*HIDX 
C 
      RETURN 
      END 
C 
C  Subroutine to compute "moderate" heat index 
C 
      SUBROUTINE SLTRYM(TC,DF,PINF,HIDX) 
C 
C  Subroutine SLTRYM calculates "moderate" heat index 
C  HIDX (C) given the temperature TC (C) and relative 
C  vapor pressure PINF. If heat index is "severe" 
C  it returns with negative value of DF. 
C 
C 
C Called by HEAT. 
C 
C Calls 
C 
C Changes: 
C  
C  21 Jan 00 @DynTel-km  1-013, New Function 
C 
C Input on call: 
C 
C   argument list : TC, DF, PINF, HIDX 
C 
C Output on return: 
C 
C   function      : HIDX, DF 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C Q, RS, PB, TB    R        Constant variables 
C ZS, EHC, PHI2    R        Constant variables 
C R, CHC, C1, C2   R        Constant variables 
C C3, C4, C5, C6   R        Constant variables 
C C7, C8, C9, C10  R        Constant variables 
C W1, W4, W5       R        Constant variables 
C HER, ERA, QV     R        Temporary variables 
C EZA, HR, ARA     R        Temporary variables 
C AZA, Q2U, QJ     R        Temporary variables 
C K, L, F, RF      R        Temporary variables 
C W2, W3, W6, W7   R        Temporary variables 
C 
C Notes : 
C 
      REAL, INTENT(IN)      :: TC 
      REAL, INTENT(INOUT)   :: DF 
      REAL, INTENT(IN)      :: PINF 
      REAL, INTENT(INOUT)   :: HIDX 
C 
      REAL Q, RS, PB, TB, ZS, EHC, PHI2, R, CHC 
      REAL HER, ERA, QV, EZA, HR, ARA, AZA, Q2U, QJ 
      REAL K, L, F, RF, W2, W3, W6, W7 
      REAL C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 
      REAL W1, W4, W5 
C 
C  Initialize variables 
C 
      Q=180. 
      RS=.0387 
      PB=5.65 
      TB=37. 
      ZS=.0521 
      EHC=17.4 
      PHI2=.84 
      R=.124 
      CHC=11.6 
      C1=4.18 
      C2=.036 
      C3=.143 
      C4=.00112 
      C5=.0168 
      C6=.060606 
      C7=3.35 
      C8=.049 
      C9=.5 
      C10=4.05 
      W1=.2016 
      W4=159.0984 
      W5=37. 
C 
C  Compute "moderate" heat index 
C 
      HER=C1+C2*TC 
      ERA=1./(EHC+HER) 
      QV=Q*(C3-C4*TC-C5*PINF) 
      EZA=C6/EHC 
      HR=C7+C8*TC 
      ARA=1./(CHC+HR) 
      AZA=C6/CHC 
      Q2U=((TB-TC)+(PB-PINF)*ERA/(ZS-EZA))/(RS+ERA) 
      QJ=(Q-QV-(1.-PHI2)*Q2U)/PHI2 
      K=(RS+ARA)+(ZS+AZA)/R-((TB-TC)+(PB-PINF)/R)/QJ 
      L=(RS+ARA)*(ZS+AZA) 
      L=(L-((TB-TC)*(ZS+AZA)+(PB-PINF)*ARA)/QJ)/R 
      F=K*K-4.*L 
      IF(F .GE. 0.0) THEN 
         RF=C9*(-K+SQRT(F)) 
         DF=60.*RF 
         IF(DF .GE. 0.0) THEN 
            W2=(1.-PHI2)/(RS+ERA) 
            W3=PHI2/(RS+RF+ARA) 
            W6=C10*ERA/(ZS+EZA) 
            W7=C10*(RF+ARA)/(ZS+R*RF+AZA) 
            HIDX=(-W4+W2*(W5+W6)+W3*(W5+W7))/(W1+W2+W3) 
         END IF 
C 
      ELSE     ! "severe" heat index 
         DF = -1.0 
      END IF 
C 
      RETURN 
      END 
C 
C  Subroutine to compute "severe" heat index 
C 
      SUBROUTINE SLTRYS(TC,PINF,HIDX) 
C 
C  Subroutine SLTRYS calculates "severe" heat index 
C  HIDX (C) given the temperature TC (C) and relative 
C  vapor pressure PINF. 
C 
C 
C Called by HEAT. 
C 
C Calls 
C 
C Changes: 
C  
C  21 Jan 00 @DynTel-km  1-013, New Function 
C 
C Input on call: 
C 
C   argument list : TC, PINF, HIDX 
C 
C Output on return: 
C 
C   function      : HIDX 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C IT               I        Loop control variable 
C ITMAX            I        Max number of iterations 
C Q, PB, TB, HC    R        Constant variables 
C ZS, EHC, PHI2    R        Constant variables 
C C1, C2, C3, C4   R        Constant variables 
C C5, C6, C7, C8   R        Constant variables 
C C9, C10, N1, N2  R        Constant variables 
C N5               R        Constant variable 
C HR, RA, ZA, QV   R        Temporary variables 
C QU, ZS, R3, C    R        Temporary variables 
C N3, N4           R        Temporary variables 
C 
C Notes : 
C 
C 
      REAL, INTENT(IN)     :: TC 
      REAL, INTENT(IN)     :: PINF 
      REAL, INTENT(INOUT)  :: HIDX 
C 
      REAL Q, RS, PB, TB, HC 
      REAL HR, RA, ZA, QV, QU, ZS, R3, C, N3, N4 
      REAL C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 
      REAL N1, N2, N5 
C 
      INTEGER IT, ITMAX 
C 
C 
C  Initialize variables 
C 
      Q=180. 
      RS=.0387 
      PB=5.65 
      TB=37. 
      HC=12.3 
      C1=4.1 
      C2=.028 
      C3=.143 
      C4=.00112 
      C5=.0168 
      C6=600000. 
      C7=.2 
      C8=.060606 
      C9=.5 
      C10=4.05 
      N1=159.0984 
      N2=37. 
      N5=.2016 
      ITMAX=10 
C 
C  Compute "severe" heat index 
C 
      HR=C1+C2*TC 
      RA=1./(HC+HR) 
      ZA=C8/HC 
      QV=Q*(C3-C4*TC-C5*PINF) 
      QU=Q-QV 
      DO IT=1,ITMAX 
         ZS=((PB-PINF)*RA)/(QU*(RS+RA)-(TB-TC))-ZA 
         IF (ZS.LT.0.) ZS=0. 
         R3=(ZS/C6)**C7 
         C=ABS(RS-R3) 
         IF (C.LE..0001) EXIT 
         RS=C9*(RS+R3) 
      END DO 
      N3=C10*RA/(ZS+ZA) 
      N4=(RS+RA) 
      HIDX=(-N1+(N2+N3)/N4)/(N5+1./N4) 
C 
      RETURN 
      END 
