C 
C  Changes: (Last change first) 
C 
C 26 Nov 97 DynTel@MLA Expanded array PBBTR to accomodate the new Heavy Duty 
C           Gasoline and LDGT vehicle categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded arrays PBBTR to accomodate the new GAS  
C           BUS vehicle type. 
C 
C 
C The data is initialized in BD25. 
C 
      INTEGER LPOD 
      REAL PBBTR 
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,14)                                                           
