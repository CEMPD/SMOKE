      REAL FUNCTION PMSCF(IV,SPEED,SADJ,METHOD,SCFC) 
C 
C  This function computes the PM speed correction factors, based 
C  on the EMFACT2000v2.02 model. 
C 
C  Changes : 
C 
C  04 Feb 01 AIR Task 11: New Function 
C 
C     Called by:  PMEMFAC 
C 
C     Calls: none 
C 
C  Input on call: 
C    Parameter list: IV,SPEED,SADJ,METHOD,SCFC 
C    Common blocks: 
C    /VTNAMS/ VTHDG2B,VTHDG8B,VTHDD2B,VTHDD8B,VTGASBUS,VTSCHOOL 
C 
C  Output on return:  
C    Function PMSCF 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C  A             R    EMFAC SCF Coefficient A 
C  B             R    EMFAC SCF Coefficient B 
C  C             R    EMFAC SCF Coefficient C 
C  D             R    EMFAC SCF Coefficient D 
C  DSPD          R    Difference between speed and reference speed 
C  METHOD        i    Method used to compute SCF 
C  RSPD          R    Reference Speed 
C  SCFC          R    EMFAC SCF Coefficient A, B, C, D 
C  SF            R    SCF at input speed 
C  SFADJ         R    SCF at reference speed 
C  SPD           R    Speed, with lower/upper limits applied 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I' 
C 
C  Declare parameter list 
C 
      INTEGER IV 
      INTEGER METHOD 
      REAL    SADJ 
      REAL    SCFC(4) 
      REAL    SPEED 
C 
C  Declare local variables 
C 
      REAL    A 
      REAL    B 
      REAL    C 
      REAL    D 
      REAL    DSPD 
      REAL    RSPD 
      REAL    SF 
      REAL    SFADJ 
      REAL    SPD 
C 
C  Assign coefficients 
C 
      A=SCFC(1) 
      B=SCFC(2) 
      C=SCFC(3) 
      D=SCFC(4) 
C 
C  Assign speeds. Valid speed range is 5.0+ mph, unless using  
C  special reference speed of 27.4 mph, then valid range is  
C  2.5 to 65.0 mph. 
C 
      SPD =  SPEED 
      RSPD = SADJ 
C 
      IF(SADJ.GE.27.) THEN 
        SPD=MAX( 2.5,SPEED) 
        SPD=MIN(65.0,SPEED) 
      ELSE 
        SPD=MAX( 5.0,SPEED) 
      ENDIF 
C 
C  Select the method/formulae to use 
C 
      SELECT CASE (METHOD) 
C 
C  METHOD 1 
C  If vehicle class is any type of bus or truck 8501 pounds or 
C  heavier, then SCF at reference speed is 1.0.   
C 
      CASE (1) 
C 
        SF = EXP(A + SPD * (B + SPD * (C + SPD * D))) 
C 
        IF(IV.GE.VTHDGV2B.AND.IV.LE.VTHDGV8B .OR. 
     *     IV.GE.VTHDDV2B.AND.IV.LE.VTHDDV8B .OR. 
     *     IV.GE.VTGASBUS.AND.IV.LE.VTSCHOOL) THEN 
          SFADJ = 1.0 
        ELSE 
          SFADJ = EXP(A + RSPD * (B + RSPD * (C + RSPD * D))) 
        ENDIF 
C 
        IF(SFADJ.GT.0.) THEN 
          PMSCF= SF / SFADJ 
        ELSE 
          PMSCF = 1.0 
        ENDIF 
C 
C  METHOD 2 
C 
      CASE (2) 
C 
        DSPD = SPD - RSPD 
C   
        PMSCF = EXP(DSPD * (A + DSPD * (B + DSPD * (C + DSPD * D)))) 
C 
C  METHOD 3 
C 
      CASE (3) 
C 
        SF = A + SPD * (B + SPD * (C + SPD * D)) 
C  
        SFADJ = A + RSPD * (B + RSPD * (C + RSPD * D)) 
C 
        IF(SFADJ.GT.0.) THEN 
          PMSCF = SF / SFADJ 
        ELSE 
          PMSCF = 1.0 
        ENDIF 
C 
      END SELECT 
C 
      END 
