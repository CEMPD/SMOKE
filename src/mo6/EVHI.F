      REAL FUNCTION EVHI(KEYEQ,RIVAL) 
C 
C  EVHI calculates the high altitude, 1972-76 model year, hot soak 
C  for LDGVs and LDGT12s given the fuel RVP or Wade Index, respectvely. 
C 
C  Called by EVMAIN. 
C 
C  Changes : (Last change first) 
C 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  24 Sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: RIVAL 
C    common blocks: 
C    /CITRV2/ RVP090,RVP115 
C    /EVAHIA/ HIHS 
C 
C  Output on return: 
C 
C    function: EVHI 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  R1       R    low RVP or diurnal Index 
C  R2       R    high RVP or diurnal Index 
C 
C  Notes: 
C 
C  Modified for MOBILE4.1 to add support for pressure/purge fail. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'EVAHIA.I' 
C 
      REAl RIVAL, R1, R2       
      INTEGER KEYEQ 
C 
C  Hot soak rate call 
C 
      R1=RVP090 
      R2=RVP115 
      EVHI=(HIHS(KEYEQ,2)*(RIVAL-R1)-HIHS(KEYEQ,1)*(RIVAL-R2))/(R2-R1) 
      RETURN 
C 
      END 
