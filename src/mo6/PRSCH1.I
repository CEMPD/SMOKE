C 
C  Changes: (Last change first) 
C 
C 16 Apr 97 - DynTel@MLA  Expanded array LVPRS to accomodate the new 
C             GAS BUS vehicle type. 
C 12 Dec 96 DynTel@MLA Expanded array LVPRS to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 
C The data is initialized in GETPRS.FOR 
C 
      INTEGER LPRSY, LPR1ST, LPRLST, LVPRS 
C 
      COMMON /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS(14) 
