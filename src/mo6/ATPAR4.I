C 
C  ATPAR4 contains several ATP variables. 
C 
C  Changes: (Last change first) 
C 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  CAN       R    Anti-tampering rate for evap cannister disablement 
C  CAP       R    Anti-tampering rate for gas cap removal 
C  EFF       R    Basic anti-tampering effectiveness rates. 
C  EFFSAV    R    Saves the default values for EFF. 
C  PCV       R    Anti-tampering rate for PCV system disablement. 
C 
C 
C  Notes: 
C 
C 24 Feb 96 @DynTel-MLA 2-611 Added EFFSAV to save default values for EFF.       
C 
C The initialization of the data can be found in BD25 
C 
      REAL CAN,PCV,CAP,EFF, EFFSAV 
      COMMON /ATPAR4/ CAN(2),PCV(2),CAP(2),EFF(2,3,2,15), 
     &                 EFFSAV(2,3,2,15) 
