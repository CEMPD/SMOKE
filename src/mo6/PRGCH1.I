C 
C  Changes: (Last change first) 
C 
C 16 Apr 97 - DynTel@MLA  Expanded array LVPRG to accomodate the new 
C             GAS BUS vehicle type. 
C 12 Dec 96 DynTel@MLA Expanded array LVPRG to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 
C The data is initialized in GETPRG.FOR 
C 
      INTEGER LPGSY, LPG1ST, LPGLST, LVPRG 
      COMMON /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG(14) 
