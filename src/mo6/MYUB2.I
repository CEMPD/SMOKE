C 
C  Changes: (Last change first) 
C 
C 26 Nov 97 DynTel@MLA Expanded array MYCCEI to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded array MYCCEI to accomodate the new GAS  
C           BUS vehicle type. 
C 
C 
C The data is initialized in BD25. 
C 
      INTEGER MYCCEI 
      COMMON /MYUB2/  MYCCEI(7,14) 
