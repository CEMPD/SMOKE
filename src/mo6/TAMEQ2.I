C 
C  Changes: (Last change first) 
C 
C 26 Nov 97 DynTel@MLA Expanded array MYGTAM to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded array MYGTAM to accomodate the new GAS  
C             BUS vehicle type. 
C 
C 
C The data is initialized in BD25. 
C 
      INTEGER MYGTAM 
      REAL F50K 
      COMMON /TAMEQ2/ MYGTAM(3,14),F50K                                                               
