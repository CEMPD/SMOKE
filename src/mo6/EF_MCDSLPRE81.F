      REAL FUNCTION EF_MCDSLPRE81(MY,IP,IV,VMTAGE,ISR) 
C 
C  EF_MCDSLPRE81 returns the basic RUNNING emission factor. Uncorrected for 
C  temperature, methane offset, temperature, reformulated fuels, 
C  inspection/maintenance and tampering etc. 
C 
C  Called by EF_DSPTCH 
C 
C  Calls IERPTR 
C 
C  Changes: (Last change first) 
C 
C  08 Sep 99 @DynTel-bag 2-660 Removed this code from BEREF.FOR, 
C            added start emissions, and renamed to EF_MCDSLPRE81.FOR. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,IP,IV,VMTAGE 
C    common blocks: 
C    /BASEQ1/ REBZML,REBDR 
C    /BASEQ5/ ERUZML,ERUDR 
C    /REGION/ IREJN 
C 
C  Output on return: 
C 
C    function: EF_MCDSLPRE81 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C  ZPOINT   R    zero mileage point = intercept of ef curve 
C  SLOPE1   R    slope of ef curve  
C  VMTAGE   R    vehicle mileage in units of 1/10,000 miles 
C 
C  Notes: 
C      
      IMPLICIT NONE 
      INCLUDE 'BASEQ1.I' 
      INCLUDE 'BASEQ5.I' 
      INCLUDE 'REGION.I' 
C 
C     Declare external functions 
C 
      INTEGER IERPTR 
C 
C     Declare parameter list 
C 
      INTEGER MY 
      INTEGER IP 
      INTEGER IV 
      INTEGER ISR 
      REAL    VMTAGE 
C 
C     Declare local variables 
C 
      REAL    ZPOINT 
      REAL    SLOPE1 
      INTEGER LMY 
      INTEGER KEYER 
      INTEGER IGER 
      INTEGER, PARAMETER :: RUN = 1, START = 2 
C 
      LMY=MY 
C 
C  Emission factors can be calculated from data in BD04, BD40 or from  
C  user input data. The integer variable KEYER distinguishes between  
C  these cases. IGER is an integer index for the model year. 
C     
      IGER=IERPTR(RUN,LMY,IP,IV,KEYER) 
C 
      SELECT CASE(ISR) 
C 
      CASE(RUN) 
C 
C  block data input BD04 
C 
      IF(KEYER.EQ.1) THEN 
        ZPOINT=REBZML( IGER, IP, IV, IREJN) 
        SLOPE1=REBDR ( IGER, IP, IV, IREJN) 
        EF_MCDSLPRE81=ZPOINT+SLOPE1*VMTAGE  
      END IF 
C 
C  user input 
C 
      IF(KEYER.EQ.2) THEN 
        ZPOINT=ERUZML(IGER,IP,IV,IREJN) 
        SLOPE1=ERUDR (IGER,IP,IV,IREJN) 
        EF_MCDSLPRE81=ZPOINT+SLOPE1*VMTAGE 
      END IF 
C 
      CASE(START) 
C 
C  BD04 block data input 
C 
      IF(KEYER.EQ.1) THEN  
        ZPOINT=SEBZML( IGER, IP, IV, IREJN ) 
        SLOPE1=SEBDR ( IGER, IP, IV, IREJN ) 
        EF_MCDSLPRE81=ZPOINT+SLOPE1*VMTAGE 
      END IF 
C 
C  user input 
C 
      IF(KEYER.EQ.2) THEN 
        ZPOINT=ESUZML( IGER, IP, IV, IREJN ) 
        SLOPE1=ESUDR ( IGER, IP, IV, IREJN ) 
        EF_MCDSLPRE81=ZPOINT+SLOPE1*VMTAGE 
      END IF 
C 
      END SELECT 
C 
      RETURN 
      END 
