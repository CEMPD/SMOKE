C 
C  Changes: (Last change first) 
C 
C 26 Nov 97 DynTel@MLA Expanded arrays MYGTS, MYEGR, and MYVTS to accomodate  
C           the new Heavy Duty Gasoline and LDGT vehicle categories. 
C 10 Apr 97 - DynTel@MLA  Expanded arrays MYEGR, MYTGS, and MYVTS to  
C             accomodate the new GAS BUS vehicle type. 
C 
C 
C 
C The data is initialized in BD25. 
C 
      INTEGER MYTGS, MYEGR, MYVTS 
      COMMON /MYUB1/  MYTGS(20,6,14),MYEGR(4,14),MYVTS(3,2,14)                                          
