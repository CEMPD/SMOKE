C 
C  Changes: (Last change first) 
C 
C  8 Apr 99 @DynTel-DDJ 2-696 Added hourly index to RLRATE array 
C 14 Apr 97 - DynTel@MLA  Expanded array RLRATE to accomodate the new GAS  
C           BUS vehicle type. 
C  3 Dec 96 DynTel@MLA Expanded array RLRATE to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 
      REAL RLRATE 
      COMMON /RLCOM3/ RLRATE(25,14,24) 
