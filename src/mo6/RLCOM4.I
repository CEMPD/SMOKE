C 
C  Changes: (Last change first) 
C 
C 13 Feb 01 AIR Task 11: New include. 
 
C  Common block dictionary: 
C 
C   Name   Type             Description 
C  ------  ----  -------------------------------------- 
C  /RLCOM4/ 
C 
C  ROADFEMC   R    road fuel economy rates for motorcycles 
C 
C  The array is initialized in BD28P5. 
C 
C 
      REAL ROADFEMC 
C 
      COMMON /RLCOM4/ ROADFEMC(12) 
