C 
C  Changes: (Last change first) 
C 
C 19 Oct 00 - @EPA-bag    Removed obsolete IOBMY and IVOB 
C 14 Apr 97 - DynTel@MLA  Expanded array IVOB and S2EFF to accomodate the  
C           new GAS BUS vehicle type. 
C  3 Dec 96 DynTel@MLA Expanded arrays IVOB and S2EFF to accomodate the 
C           new Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD28. 
C 
      INTEGER IS2SY, NPHASE 
      REAL S2EFF 
C                       
      COMMON /RLCOM2/ IS2SY,NPHASE,S2EFF(14) 
                      
