C 
C  Changes: (Last change first) 
C 
C 27 Nov 97 DynTel@MLA Expanded array LVTFLG to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded array LVTFLG to accomodate the new GAS  
C             BUS vehicle type. 
C 
C The initialization of the data can be found in BD25. 
C 
      INTEGER LAPSY,LAP1ST,LAPLST,LVTFLG 
      COMMON /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(14) 
