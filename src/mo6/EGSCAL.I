C 
C  Changes: (Last change first) 
C 
C 26 Nov 97 DynTel@MLA Expanded array TGS to accomodate the new Heavy Duty 
C           Gasoline and LDGT vehicle categories. 
C 
C 16 Apr 97 - DynTel@MLA  Expanded array TGS to accomodate the new GAS BUS 
C             vehicle type. 
C 
C 
C The data is initilized in BD25. 
C 
      REAL AER,TGS,EGS 
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(20,6,14),EGS(7,2) 
