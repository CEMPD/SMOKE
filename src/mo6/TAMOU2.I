C 
C  Changes: (Last change first) 
C 26 Nov 97 DynTel@MLA Expanded arrays GCONLY and TOB to accomodate the new  
C           Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded arrays GCONLY and TOB to accomodate the  
C            new GAS BUS vehicle types. 
C 
C 
C The data is initialized in BD25. 
C 
      REAL TOB, GCONLY 
      COMMON /TAMOU2/ TOB(25,14),GCONLY(25,14)                                                         
